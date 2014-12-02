// rm matrix-cache.exe; mpic++ matrix-cache.cpp -o matrix-cache.exe -mavx -O3 -Wall -std=c++0x; ./matrix-cache.exe
// salloc -n 1 -p sandy; mpic++ matrix-cache.cpp -o matrix-cache -mavx -Wall -std=c++0x; srun ./matrix-cache
#include <iostream>
#include <cassert>
#include <cstdio>
#include <mpi.h>
#include <ctime>
#include <cstdint>  // uintxx_t
#include <cstdlib>
#include <fstream>  // for teeStream class
#include <cmath>
#include <string>
#include <sstream>  // stringstream manipulation to add numbers
#include <unistd.h> // sleep
#include <iomanip>  // setprecision

using namespace std;

const uint32_t N_MAX_DIM     = 2000; //300mb limit -> 100mb per matrix, double is 8 byte -> 13107200 elements -> 3620 Elements
const uint32_t N_POINTS      = 100;
const uint32_t N_REPETITIONS = 10;

/* if this is 1 files will be flushed for in-situ analysis and other feedback *
 * will be printed on the console                                             */

/* outstream similar to cout but outputs into files */
class teeStream{
    public:
        ofstream fileStream;
        teeStream( void ) { }
        teeStream( const char * filename ) {
            fileStream.open( filename, std::ofstream::out | std::ofstream::app );
        }
        void open( string filename ) {
            fileStream.open( filename, std::ofstream::out | std::ofstream::app );
        }
        ~teeStream(void) {
            fileStream.close();
        }
        template <class T> teeStream& operator<< (T val) {
            fileStream << val;
            cout << val;
            fileStream.flush();
            return *this;
        }
        teeStream& operator<< (ostream& (*pfun)(ostream&)) {
            pfun(fileStream);
            pfun(cout);
            return *this;
        }
};

teeStream tout;

uint64_t GetCPUCycle() {
    // used ressources: Gcc manual, http://www.strchr.com/performance_measurements_with_rdtsc
    // ToDo: Read gcc manual 6.4.2 more thoroughly
    // to be used for time measurements (a bit too exact, it also would need cpu clocks per second) or as a random seed
    uint32_t lo,hi;
    asm volatile (
        ".intel_syntax noprefix \n"
        "xor eax,eax            \n"
        "cpuid                  \n" // complete preceding instructions in pipeline
        "rdtsc                  \n"
        ".att_syntax            \n"
        //"mov [_hi],eax          \n"
        : "=d" (hi), "=a" (lo)  /* out */
        : /* no input */
        : "cc" /* cluttered registers */
    );
    return (uint64_t(hi) << 32) + uint64_t(lo);
}


template <typename MATRIX_DATATYPE>
void MatrixMult( const uint32_t matrixsize, const MATRIX_DATATYPE* A, const MATRIX_DATATYPE* B, MATRIX_DATATYPE* C ) {
    memset( C, 0, matrixsize * matrixsize * sizeof(MATRIX_DATATYPE) );
    for (uint32_t i=0; i<matrixsize; i++) {     //go down every line
        for (uint32_t k=0; k<matrixsize; k++) {     //go down every column in *this and every line in mat
            MATRIX_DATATYPE r = A[ i * matrixsize + k ];
            for (uint32_t j=0; j<matrixsize; j++) {     //go down every column
                C[ i * matrixsize + j ] += r * B[ k * matrixsize + j ];
            }
        }
    }
}

int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    
    /* Create Timestamp and rank strings for filenames */
    time_t t = time(0);   // get time now
    struct tm * now = localtime( &t );
    stringstream timestamp;
    timestamp << 1900 + now->tm_year << "-" << 1 + now->tm_mon << "-"
                << now->tm_mday << "_" << now->tm_hour << "-" << now->tm_min
                << "_";
    string timesr = timestamp.str();

    // Open File for times
    FILE * datafile = NULL;
    datafile = fopen( (timesr+string("data.txt")).c_str(), "w" );
    fprintf( datafile, "# Square Matrix Size\tdatatype size/Bytes\tTime for Matrixmult\tStdDev for Time measurement\n" );
    fflush(datafile);
    
    // Open Log File, Save some Stats
    tout.open( timesr+string(".txt") );
    tout << "Processor Name: " << processor_name << endl;
    double   t0    = MPI_Wtime();
    uint64_t nc0   = GetCPUCycle();
    sleep(10);     
    uint64_t nc1   = GetCPUCycle();
    double   t1    = MPI_Wtime();
    double cpufreq = (nc1-nc0)/(t1-t0);
    tout << "MPI_Wtime resolution: %e" << MPI_Wtick() << endl;
    tout << "CPU Frequency : (" << cpufreq/1e9 << " +- "
         << cpufreq/1e9 * MPI_Wtick() / (t1-t0) << ") GHz\n";
    
    /* Parameters and Levels:                                                *
     *     Datatype    : uint8_t, uint16_t, uint32_t, uint64_t, uint128_t    *
     *                   float, double                                       *
     *     Memory Order: row major, column major                             *
     *     Problem size: different kind of MB sized Matrices                 */
    #define MDTYPE double
    MDTYPE * A = (MDTYPE*) malloc( N_MAX_DIM * N_MAX_DIM * sizeof(MDTYPE) );
    MDTYPE * B = (MDTYPE*) malloc( N_MAX_DIM * N_MAX_DIM * sizeof(MDTYPE) );
    MDTYPE * C = (MDTYPE*) malloc( N_MAX_DIM * N_MAX_DIM * sizeof(MDTYPE) );

    // Fill Matrices with Random Data
    cout << "Filling Matrices with Random Data\n";
    srand(nc1);
    for (uint32_t i = 0; i < N_MAX_DIM; i++)
        for (uint32_t j = 0; j < N_MAX_DIM; j++) {
            A[ i * N_MAX_DIM + j ] = (MDTYPE) rand();
            B[ i * N_MAX_DIM + j ] = (MDTYPE) rand();
        }
    
    // Create Array with levels for Matrix-Size
    uint32_t * NLevels = (uint32_t*) malloc( N_POINTS * sizeof(uint32_t) );
    const float N_MULT = pow( N_MAX_DIM, (float) 1. / N_POINTS );
    uint32_t lastLevel = 0, k=0, NWritten=0;
    while (NWritten < N_POINTS) {
        uint32_t N = floor( pow( N_MULT, k ) );
        k++;
        if (N > N_MAX_DIM)
            break;
        if (N != lastLevel) { // dismiss double values
            NLevels[NWritten] = N;
            NWritten++;
            lastLevel = N;
        }
    }
    
    // Do actual Benchmark
    for (uint32_t m = 1; m <= NWritten-1 ; m++) {
        uint32_t matrixsize = NLevels[m];
        cout << "matrixsize: " << matrixsize << flush;
         
        uint64_t t0, t1, sum=0, sum2=0;
        for (uint32_t t = 0; t < N_REPETITIONS; t++) {
            // wait a bit for every time measurement trying to reduce variations 
            t0 = GetCPUCycle();
            while ( GetCPUCycle() - t0 < 1e4 ) {}

            t0       = GetCPUCycle();
            MatrixMult( matrixsize, A, B, C );
            t1       = GetCPUCycle();
                
            sum  += t1-t0;
            sum2 += (t1-t0)*(t1-t0);
        }

        const double mean   = double(sum) / N_REPETITIONS / cpufreq;
        const double stddev = sqrt( double(sum2-sum)/(N_REPETITIONS*(N_REPETITIONS-1)) ) /cpufreq;
        fprintf( datafile, "%u\t%lu\t%e\t%e\n", matrixsize, sizeof(MDTYPE), mean, stddev );
        cout << " => Needed: " << mean << " +- " << stddev << " sec\n";
        fflush(datafile);
    } // loop over matrix size*/

    free(A); free(B); free(C);// free(NLevels);
    if (datafile != NULL) fclose(datafile);
    MPI_Finalize();
    return 0;
}
