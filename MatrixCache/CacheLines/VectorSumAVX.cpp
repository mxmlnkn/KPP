/* ToDo:
    - Besseres Sampling überlegen (Gesamt-Matrix 10MB mitsamt Stride und dann
      halt nur so und so viel Stride möglich (kreisförmig)
    - 3D-Plot
    - Korrekte Umrechnung Stride in Byte, einzeichnen Caches
    - Caching genauer überlegen => welcher maximaler stride sinnvoll?
        . nur ein Element in Cache ladbar? -> wird sowieso nur CacheLine geladen
        . also nur bis Cacheline-Size sinnvoll? Aber warum dann Unterschied
          Stride 64 und 128 -> Mindestens 2048Byte lange CL?
    - GFlops so wenig, warum? -> doch mit Matrix-Mul probieren, oder gar mit
      Matrix-Potenz, um mal Peak Flop auszureizen?
    - Doch wieder auf logarithmisch umsteigen :, damit wenigstens bis 10^8 für
      Stride 0 möglich ist
*/

//KurzForm: rm VectorSumAVX.exe; mpic++ VectorSumAVX.cpp -o VectorSumAVX.exe -Wall -std=c++0x -mavx -fopenmp -O3 -DNDEBUG

// salloc -p sandy --nodes 1 --exclusive --ntasks-per-node 4 --time 120
// rm VectorSumAVX.exe; mpic++ VectorSumAVX.cpp -o VectorSumAVX.exe -Wall -std=c++0x -mavx -fopenmp -O3 -DNDEBUG; srun -n 1 ./VectorSumAVX.exe

#include <iostream>
#include <cassert>
#include <cstdio>
#include <mpi.h>
#include <omp.h>     // for omp_set_num_threads
#include <ctime>
#include <cstdint>   // uintxx_t */
#include <cstdlib>   // srand,...
#include <fstream>   // for teeStream class
#include <cmath>
#include <cstring>   // memset
#include <string>
#include <sstream>   // stringstream manipulation to add numbers
#include <unistd.h>  // sleep

using namespace std;


#define NUMBER_OF_THREADS 16
#define MDTYPE double
#define VECTORSIZE_BYTES 32
// avx is 256 bit = 4 double = 32 byte vector of four single floats
const uint16_t VECTORSIZE_ELEMENTS = VECTORSIZE_BYTES / sizeof(MDTYPE);
typedef double v4sf __attribute__((vector_size(VECTORSIZE_BYTES)));

union d4vector {
  v4sf v;
  MDTYPE d[VECTORSIZE_ELEMENTS] __attribute__((aligned(VECTORSIZE_BYTES)));
};

const uint32_t N_MIN_LENGTH            = 1024;    //1;     // in multiples of VECTORSIZE_BYTES
const uint32_t N_MAX_LENGTH            = 224*1024;  // in multiples of VECTORSIZE_BYTES
const double   REQUIRED_RELATIVE_ERROR = 0.05;
const uint32_t N_REPETITIONS_OFFSET    = 0;
const uint32_t NTH_SAVE_DETAILED       = 100;
const uint32_t N_MIN_STRIDE            = 0;    // Stride in Multiples of VECTORSIZE_BYTES
const uint32_t N_MAX_STRIDE            = 10;

uint64_t GetCPUCycle() {
    uint32_t lo,hi;
    asm volatile (
        ".intel_syntax noprefix \n"
        "rdtscp                 \n"
        ".att_syntax            \n"
        : "=d" (hi), "=a" (lo)  /* out */
        : /* no input */
        : "cc","rcx" /* clobbered registers */
    );
    return (uint64_t(hi) << 32) + uint64_t(lo);
}


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


template <typename BASE_DATATYPE>
void VectorAdd( const uint32_t length, const uint32_t stride,
                const BASE_DATATYPE* A, const BASE_DATATYPE* B, BASE_DATATYPE* C )
                // stride and length in VECTORSIZE_BYTES
{
    // N_MAX_DIM must be multiple of 4 because of avx! (just pad the rest with zeros
    d4vector * const Av = (d4vector*) A;
    d4vector * const Bv = (d4vector*) B;
    d4vector * const Cv = (d4vector*) C;
    #pragma omp parallel for num_threads(NUMBER_OF_THREADS)
    for (uint32_t i=0; i<length; i++)
        Cv[i*(1+stride)].v = Av[i*(1+stride)].v + Bv[i*(1+stride)].v;
}

/******************************************************************************/

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
    fprintf( datafile, "# Array Length\tdatatype size/Bytes\tTime for Matrixmult\tStdDev for Time measurement\tStride/Bytes\n" );
    fflush(datafile);

    // Open Log File, Save some Stats
    tout.open( timesr+string(".txt") );
    tout << "Processor Name: " << processor_name << endl;
    double   t0    = MPI_Wtime();
    uint64_t nc0   = GetCPUCycle();
    while( MPI_Wtime() - t0 < 10. ) {}
    uint64_t nc1   = GetCPUCycle();
    double   t1    = MPI_Wtime();
    double cpufreq = (nc1-nc0)/(t1-t0);
    tout << ", t0=" << t0 << ", t1=" << t1 << ", nc0=" << nc0 << ", nc1=" << nc1 << endl;
    tout << "MPI_Wtime resolution: " << MPI_Wtick() << endl;
    tout << "CPU Frequency : (" << cpufreq/1e9 << " +- "
         << cpufreq/1e9 * MPI_Wtick() / (t1-t0) << ") GHz\n" << flush;

    /* Allocate Matrices */
    void     * Am = malloc( N_MAX_LENGTH * ( 1 + N_MAX_STRIDE ) * sizeof(MDTYPE) + (VECTORSIZE_BYTES-1) );
    void     * Bm = malloc( N_MAX_LENGTH * ( 1 + N_MAX_STRIDE ) * sizeof(MDTYPE) + (VECTORSIZE_BYTES-1) );
    void     * Cm = malloc( N_MAX_LENGTH * ( 1 + N_MAX_STRIDE ) * sizeof(MDTYPE) + (VECTORSIZE_BYTES-1) );
    d4vector * Av = (d4vector*) ( ((uintptr_t)Am+(VECTORSIZE_BYTES-1)) & ~ (uintptr_t)(VECTORSIZE_BYTES-1) );
    d4vector * Bv = (d4vector*) ( ((uintptr_t)Bm+(VECTORSIZE_BYTES-1)) & ~ (uintptr_t)(VECTORSIZE_BYTES-1) );
    d4vector * Cv = (d4vector*) ( ((uintptr_t)Cm+(VECTORSIZE_BYTES-1)) & ~ (uintptr_t)(VECTORSIZE_BYTES-1) );
    // if above is not aligned to 32 byte, avx-access will segfault !!!
    MDTYPE   * A  = (MDTYPE*) Av;
    MDTYPE   * B  = (MDTYPE*) Bv;
    MDTYPE   * C  = (MDTYPE*) Cv;
    if (Av==NULL or Bv==NULL or Cv==NULL)
        return 1;

    // Fill Matrices with Random Data
    cout << "Filling Vectors\n" << flush;
    for (uint32_t i = 0; i < N_MAX_LENGTH*(1+N_MAX_STRIDE); i++)
    {
        A[i] = (MDTYPE) rand();
        B[i] = (MDTYPE) rand();
    }

    FILE * detailedlogfile = NULL;
    detailedlogfile = fopen( (timesr+string("Times.txt")).c_str(), "w" );
    fprintf( detailedlogfile, "#MatrixSize / Elements \tCPU clocks 1\tCPU clocks 2\t...\n" );

    // Do actual Benchmark
    uint32_t benchmarksDone = 0;
    for (uint32_t stride = N_MIN_STRIDE; stride <= N_MAX_STRIDE; stride+=1) {
        const double stride_t0 = MPI_Wtime();
        for (uint32_t length = N_MIN_LENGTH; length <= N_MAX_LENGTH; length+=1024)
        {
            uint64_t sum=0, sum2=0;
            uint32_t NRepetitions = 0;

            while(true)
            {
                uint64_t t0 = GetCPUCycle();
                VectorAdd( length, stride, A, B, C );
                uint64_t t1 = GetCPUCycle();

                if (t < N_REPETITIONS_OFFSET)
                    continue;
                if (benchmarksDone % NTH_SAVE_DETAILED == 0)
                    fprintf( detailedlogfile, "%lu\t", t1-t0 );

                sum  += t1-t0;
                sum2 += (t1-t0)*(t1-t0);
                NRepetitions++;

                if ( NRepetitions >= 3 ) {
                    const double mean   = double(sum) / NRepetitions / cpufreq;
                    const double stddev = sqrt( double(sum2-sum)/(NRepetitions*(NRepetitions-1)) ) /cpufreq;
                    if (stddev/mean <= REQUIRED_RELATIVE_ERROR )
                        break;
                }
            }
            if (benchmarksDone % NTH_SAVE_DETAILED == 0) {
                fprintf( detailedlogfile, "\n" );
                fflush ( detailedlogfile );
            }

            const double mean   = double(sum) / NRepetitions / cpufreq;
            const double stddev = sqrt( double(sum2-sum)/(NRepetitions*(NRepetitions-1)) ) /cpufreq;
            fprintf( datafile, "%u\t%lu\t%e\t%e\t%u\n", length, sizeof(MDTYPE), mean, stddev, stride );
            //tout << " => Needed: " << mean << " +- " << stddev << " sec";
            //float gflops = 2.*length/mean/1e9;
            //tout << " => " << gflops << " +- " << stddev/mean * gflops << " GFlops";
            //tout << " N_Repetitions: " << NRepetitions << endl;
            fflush(datafile);

            benchmarksDone++;
        } // loop over matrix size*/
        
        const double stride_t1 = MPI_Wtime();
        tout << "Stride:" << stride << " took " << stride_t1-stride_t0 << " seconds" << endl << flush;
    }

    if (detailedlogfile != NULL) fclose(detailedlogfile);
    free(Am); free(Bm); free(Cm);
    if (datafile != NULL) fclose(datafile);
    MPI_Finalize();
    return 0;
}
