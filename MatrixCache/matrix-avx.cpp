/* TODO:
    - Test also for Matrices with sizes not module 4 (double vector size avx)
      by filling the rest with 0.
        -> is fmul 0,b faster than fmul a,b ?
    - Compare 2n x n vs. n x 2n matrices -> different stride -> slower?
        -> implement arbitrary m x n also if not mod 4
    - parallelize memset (C,0), e.g. with double sum in 2nd innermost loop
    - automatize Matrix-Correctness-Check: 
           A = ( 1 2 ) , B = ( 1 1 )  => C = A*B = ( 1 1 ) * (1*1+2*2)
               ( 1 2 )       ( 2 2 )               ( 1 1 )
        => problematic if not reallocated and reinitialized every timethe size
           changes ...
*/


// rm matrix-avx.exe; rm matrix-avx.o; mpic++ -g -c matrix-avx.cpp -o matrix-avx.o -Wall -std=c++0x -mno-mmx -mno-sse2 -mno-sse3 -mno-sse4 -mno-sse4.1 -mno-sse4.2 -mavx -mno-avx2 -fopenmp -DNDEBUG; mpic++ matrix-avx.o -o matrix-avx.exe -fopenmp; ./matrix-avx.exe
// objdump -dS ./example1.o  | grep -22 c.v | tail -25
// watch out! no proper debugging information in objdump if compiled with -O3 !
// too many libraries (e.g. mpi.h) may also complicate the objdump disassembled reading

//KurzForm: rm matrix-avx.exe; rm matrix-avx.o; mpic++ matrix-avx.cpp -o matrix-avx.exe -Wall -std=c++0x -mavx -fopenmp -O3 -DNDEBUG

// salloc -p sandy --nodes 1 --exclusive --ntasks-per-node 4 --time 120
// mpic++ matrix-avx.cpp -o matrix-avx.exe -Wall -std=c++0x -mavx -fopenmp -O3 -DNDEBUG
// srun -n 1 ./matrix-avx.exe

/*
    c.v = a.v + b.v;
  b6:   48 8b 83 a0 00 00 00    mov    0xa0(%rbx),%rax
  bd:   48 89 43 40             mov    %rax,0x40(%rbx)
  c1:   48 8b 83 a8 00 00 00    mov    0xa8(%rbx),%rax
  c8:   48 89 43 48             mov    %rax,0x48(%rbx)
  cc:   48 8b 83 b0 00 00 00    mov    0xb0(%rbx),%rax
  d3:   48 89 43 50             mov    %rax,0x50(%rbx)
  d7:   48 8b 83 b8 00 00 00    mov    0xb8(%rbx),%rax
  de:   48 89 43 58             mov    %rax,0x58(%rbx)
  e2:   48 8b 83 80 00 00 00    mov    0x80(%rbx),%rax
  e9:   48 89 43 20             mov    %rax,0x20(%rbx)
  ed:   48 8b 83 88 00 00 00    mov    0x88(%rbx),%rax
  f4:   48 89 43 28             mov    %rax,0x28(%rbx)
  f8:   48 8b 83 90 00 00 00    mov    0x90(%rbx),%rax
  ff:   48 89 43 30             mov    %rax,0x30(%rbx)
 103:   48 8b 83 98 00 00 00    mov    0x98(%rbx),%rax
 10a:   48 89 43 38             mov    %rax,0x38(%rbx)
 10e:   dd 43 40                fldl   0x40(%rbx)
 111:   dd 43 20                fldl   0x20(%rbx)
 114:   de c1                   faddp  %st,%st(1)
 116:   dd 43 48                fldl   0x48(%rbx)
 119:   dd 43 28                fldl   0x28(%rbx)
 11c:   de c1                   faddp  %st,%st(1)
 11e:   dd 43 50                fldl   0x50(%rbx)
 121:   dd 43 30                fldl   0x30(%rbx)
 124:   de c1                   faddp  %st,%st(1)
 126:   dd 43 58                fldl   0x58(%rbx)
 129:   dd 43 38                fldl   0x38(%rbx)
 12c:   de c1                   faddp  %st,%st(1)
 12e:   d9 cb                   fxch   %st(3)
 130:   dd 5d d0                fstpl  -0x30(%rbp)
 133:   d9 c9                   fxch   %st(1)
 135:   dd 5d d8                fstpl  -0x28(%rbp)
 138:   dd 5d e0                fstpl  -0x20(%rbp)
 13b:   dd 5d e8                fstpl  -0x18(%rbp)
 13e:   48 8b 45 d0             mov    -0x30(%rbp),%rax
 142:   48 89 03                mov    %rax,(%rbx)
 145:   48 8b 45 d8             mov    -0x28(%rbp),%rax
 149:   48 89 43 08             mov    %rax,0x8(%rbx)
 14d:   48 8b 45 e0             mov    -0x20(%rbp),%rax
 151:   48 89 43 10             mov    %rax,0x10(%rbx)
 155:   48 8b 45 e8             mov    -0x18(%rbp),%rax
 159:   48 89 43 18             mov    %rax,0x18(%rbx)
 15d:   48 8b 03                mov    (%rbx),%rax
 160:   48 89 43 60             mov    %rax,0x60(%rbx)
 164:   48 8b 43 08             mov    0x8(%rbx),%rax
 168:   48 89 43 68             mov    %rax,0x68(%rbx)
 16c:   48 8b 43 10             mov    0x10(%rbx),%rax
 170:   48 89 43 70             mov    %rax,0x70(%rbx)
 174:   48 8b 43 18             mov    0x18(%rbx),%rax
 178:   48 89 43 78             mov    %rax,0x78(%rbx)

    printf("%e, %e, %e, %e\n", c.d[0], c.d[1], c.d[2], c.d[3]);
*/
// was machen die ganzen mov oben ???

// -mno-mmx -mno-sse -mno-sse2 -mno-sse3 -mno-sse4 -mno-sse4.1 -mno-sse4.2 -mno-avx2 -mavx
/*
    c.v = a.v + b.v;
  97:   c5 fd 28 4b 40          vmovapd 0x40(%rbx),%ymm1
  9c:   c5 fd 28 43 20          vmovapd 0x20(%rbx),%ymm0
  a1:   c5 f5 58 c0             vaddpd %ymm0,%ymm1,%ymm0
  a5:   c5 fd 29 03             vmovapd %ymm0,(%rbx)

    printf("%e, %e, %e, %e\n", c.d[0], c.d[1], c.d[2], c.d[3]);
*/

/* Notebook:
    -mno-avx:
        Cycles needed: 1010645126
        Cycles needed: 1361462411
        Cycles needed: 1717592451
        Cycles needed: 1010677637
        Cycles needed: 1181628134
    -mavx
        */

#include <iostream>
#include <cassert>
#include <cstdio>
#include <mpi.h>
#include <omp.h>   // for omp_set_num_threads
#include <ctime>
#include <cstdint>  // uintxx_t */
#include <cstdlib>   // srand,...
#include <fstream>  // for teeStream class
#include <cmath>
#include <cstring>    // memset
#include <string> 
#include <sstream>  // stringstream manipulation to add numbers
#include <unistd.h> // sleep
//#include <iomanip>  // setprecision

using namespace std;


#define MDTYPE double
// avx is 256 bit = 4 double = 32 byte vector of four single floats
#define VECTORSIZE_BYTES 32
const uint16_t VECTORSIZE_ELEMENTS = VECTORSIZE_BYTES / sizeof(MDTYPE);
#define NUMBER_OF_THREADS 16
typedef double v4sf __attribute__((vector_size(VECTORSIZE_BYTES)));  

union d4vector {
  v4sf v;
  MDTYPE d[VECTORSIZE_BYTES/sizeof(MDTYPE)] __attribute__((aligned(VECTORSIZE_BYTES)));
};

const uint32_t N_MIN_DIM     = 1;
const uint32_t N_MAX_DIM     = 3000; //300mb limit -> 100mb per matrix, double is 8 byte -> 13107200 elements -> 3620 Elements
const uint32_t N_POINTS      = 1000;
const bool     LOGARITHMIC   = false;
const uint32_t N_REPETITIONS_MAX = 1000;
const uint32_t N_REPETITIONS_MIN = 10;
const uint32_t MIN_TIME_SECONDS  = 5;
const uint32_t N_REPETITIONS_OFFSET = 0;
const uint32_t NTH_SAVE_DETAILED = 1;

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
        : "cc","ebx","ecx" /* clobbered registers */
    );
    return (uint64_t(hi) << 32) + uint64_t(lo);
}


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

template <typename MATRIX_DATATYPE>
void MatrixMemcpy( const uint32_t matrixsize, const MATRIX_DATATYPE* A, MATRIX_DATATYPE* C ) {
    memset( C, 0, matrixsize * matrixsize * sizeof(MATRIX_DATATYPE) );
    #pragma omp parallel for num_threads(NUMBER_OF_THREADS)
    for (int i = 0; i<NUMBER_OF_THREADS; i++)
        memcpy( &(C[matrixsize * matrixsize / NUMBER_OF_THREADS]),
                &(A[matrixsize * matrixsize / NUMBER_OF_THREADS * i ]),
                matrixsize * matrixsize / NUMBER_OF_THREADS * sizeof(MATRIX_DATATYPE) );
}

template <typename MATRIX_DATATYPE>
void MatrixMultOpenMPavx( const uint32_t matrixsize, const MATRIX_DATATYPE* A, const MATRIX_DATATYPE* B, MATRIX_DATATYPE* C ) {
    // N_MAX_DIM must be multiple of 4 because of avx! (just pad the rest with zeros
    d4vector * const Bv = (d4vector*) B;
    d4vector * const Cv = (d4vector*) C;
    uint32_t N_DIM_AVX = matrixsize / (VECTORSIZE_ELEMENTS); 
    //tout << "MatrixSize: " << matrixsize << " mod " << VECTORSIZE_ELEMENTS << endl;
    assert( matrixsize % VECTORSIZE_ELEMENTS == 0 );
    memset( C, 0, matrixsize * matrixsize * sizeof(MDTYPE) );
    #pragma omp parallel for num_threads(NUMBER_OF_THREADS)
    for (uint32_t i=0; i<matrixsize; i++) {     //go down every line
        for (uint32_t k=0; k<matrixsize; k++) {     //go down every column in *this and every line in mat
            MDTYPE r = A[ i * matrixsize + k ];
            for (uint32_t j=0; j<N_DIM_AVX; j++) {     //go down every column
                Cv[ i * N_DIM_AVX + j].v += r * Bv[ k * N_DIM_AVX + j ].v;
            }
        }
    }
}

void PrintMatrix( const uint32_t matrixsize, const MATRIX_DATATYPE* data )
    for (int i=0; i < matrixsize; i++) {
        cout << "[";
        for (int j=0; j < matrixsize; j++) {
            cout << setw(4) << data[i*matrixsize+j] << " ";
        }
        cout << "]" << endl;
    }
    return;

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
    fprintf( datafile, "# Square Matrix Size\tdatatype size/Bytes\tTime for Matrixmult\tStdDev for Time measurement\n" );
    fflush(datafile);
    
    // Open Log File, Save some Stats
    tout.open( timesr+string(".txt") );
    tout << "Processor Name: " << processor_name << endl;
    double   t0    = double(clock()) / (double)CLOCKS_PER_SEC; //MPI_Wtime();
    uint64_t nc0   = GetCPUCycle();
    while( GetCPUCycle()-nc0 < 10*3e9 ) {}
    uint64_t nc1   = GetCPUCycle();
    double   t1    = double(clock()) / (double)CLOCKS_PER_SEC;
    double cpufreq = (nc1-nc0)/(t1-t0);
    tout << "clock=" << clock() << ", CLOCKS_PER_SEC=" << CLOCKS_PER_SEC;
    tout << ", t0=" << t0 << ", t1=" << t1 << ", nc0=" << nc0 << ", nc1=" << nc1 << endl;
    tout << "MPI_Wtime resolution: " << MPI_Wtick() << endl;
    tout << "CPU Frequency : (" << cpufreq/1e9 << " +- "
         << cpufreq/1e9 * MPI_Wtick() / (t1-t0) << ") GHz\n";
    
    
    /* Parameters and Levels:                                                *
     *     Datatype    : uint8_t, uint16_t, uint32_t, uint64_t, uint128_t    *
     *                   float, double                                       *
     *     Memory Order: row major, column major                             *
     *     Problem size: different kind of MB sized Matrices                 */
    void     * Am = malloc( N_MAX_DIM * N_MAX_DIM * sizeof(MDTYPE) + (VECTORSIZE_BYTES-1) );
    void     * Bm = malloc( N_MAX_DIM * N_MAX_DIM * sizeof(MDTYPE) + (VECTORSIZE_BYTES-1) );
    void     * Cm = malloc( N_MAX_DIM * N_MAX_DIM * sizeof(MDTYPE) + (VECTORSIZE_BYTES-1) );
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
    cout << "Filling Matrices with Random Data\n";
    srand(GetCPUCycle());
    for (uint32_t i = 0; i < N_MAX_DIM; i++)
        for (uint32_t j = 0; j < N_MAX_DIM; j++) {
            A[ i * N_MAX_DIM + j ] = i; //(MDTYPE) rand();
            B[ i * N_MAX_DIM + j ] = j; //(MDTYPE) rand();
        }
    // first row of A is in memory, then row 2 follows
    //tout << "B:" << endl;
    //PrintMatrix( 
    
    //Cv[0].v = Bv[0].v; // this line didn't work after GetCPUCycle, because latter forgot to declare ebx and ecx as clobbered !!!
    

    // Create Array with levels for Matrix-Size
    uint32_t * NLevels = (uint32_t*) malloc( N_POINTS * sizeof(uint32_t) );
    const float N_MULT = pow( N_MAX_DIM-N_MIN_DIM , (float) 1. / N_POINTS );
    uint32_t lastLevel = 0, k=0, NWritten=0;
    while (NWritten < N_POINTS) {
        if (LOGARITHMIC) {
            uint32_t N = N_MIN_DIM + floor( pow( N_MULT, k ) );
            // round to next multiple of 4 for double (4 vector in avx)
            N = (N / VECTORSIZE_ELEMENTS) * VECTORSIZE_ELEMENTS;
            k++;
            if (N > N_MAX_DIM)
                break;
            if (N != lastLevel) { // dismiss double values
                NLevels[NWritten] = N;
                NWritten++;
                lastLevel = N;
            }
        } else {
            if ( 4*(k+1) > N_MAX_DIM )
                break;
            NLevels[NWritten++] = 4*(++k);
        }
    }
    
    /*
    nc0 = GetCPUCycle();
    uint32_t matrixsize = N_MAX_DIM;

    nc1 = GetCPUCycle(); 
    cout << "Cycles needed: " << nc1-nc0
         << " => " << 2*matrixsize*matrixsize*matrixsize/((nc1-nc0)/cpufreq)/1e9 << " Gflops\n";
    */
    
    FILE * detailedlogfile = NULL;
    detailedlogfile = fopen( (timesr+string("Times.txt")).c_str(), "w" );
    fprintf( detailedlogfile, "#MatrixSize / Elements \tCPU clocks 1\tCPU clocks 2\t...\n" );
    
    
    // Do actual Benchmark
    for (uint32_t m = 0; m <= NWritten-1 ; m++) 
    {
        uint32_t matrixsize = NLevels[m];
        tout << "matrixsize: " << matrixsize << flush;
        if (m % NTH_SAVE_DETAILED == 0)
            fprintf( detailedlogfile, "%u\t", matrixsize );
         
        uint64_t sum=0, sum2=0;
        uint64_t clocks_begin=GetCPUCycle();
        uint32_t NRepetitions = 0;
        
        for (uint32_t t=0; t < N_REPETITIONS_MAX+N_REPETITIONS_OFFSET; t++)
        {
            uint64_t t0 = GetCPUCycle();
            MatrixMultOpenMPavx( matrixsize, A, B, C );
            //MatrixMemcpy( matrixsize, A, C );
            uint64_t t1 = GetCPUCycle();
            
            // Check Matrix for Errors
            /*double multiplier = pow( matrixsize, 3)/3. + pow( matrixsize, 2)/2. + matrixsize/6.;
            for (uint32_t i = 0; i < N_MAX_DIM; i++)
                for (uint32_t j = 0; j < N_MAX_DIM; j++) 
                    if ( fabs( C[ i * N_MAX_DIM + j ] - multiplier ) > 1e-5) {
                        tout << "Result Matrix Wrong !!! C[" << i << "," << j << "]=" << C[ i * N_MAX_DIM + j ] << "!=" << multiplier << endl;
                        return 1;
                    }
            */
            
            if (t < N_REPETITIONS_OFFSET)
                continue;
            if (m % NTH_SAVE_DETAILED == 0) 
                fprintf( detailedlogfile, "%lu\t", t1-t0 );
                
            sum  += t1-t0;
            sum2 += (t1-t0)*(t1-t0);
            NRepetitions++;
            
            if ( NRepetitions >= N_REPETITIONS_MIN and
            (GetCPUCycle()-clocks_begin) > MIN_TIME_SECONDS*cpufreq )
                break;
        }
        if (k % NTH_SAVE_DETAILED == 0) {
            fprintf( detailedlogfile, "\n" );
            fflush ( detailedlogfile );
        }

        const double mean   = double(sum) / NRepetitions / cpufreq;
        const double stddev = sqrt( double(sum2-sum)/(NRepetitions*(NRepetitions-1)) ) /cpufreq;
        fprintf( datafile, "%u\t%lu\t%e\t%e\n", matrixsize, sizeof(MDTYPE), mean, stddev );
        tout << " => Needed: " << mean << " +- " << stddev << " sec";
        float gflops = 2.*pow(matrixsize,3)/mean/1e9;
        tout << " => " << gflops << " +- " << stddev/mean * gflops << " GFlops";
        tout << " N_Repetitions: " << NRepetitions << endl;
        fflush(datafile);
    } // loop over matrix size*/

    
    /*
    union d4vector a, b, c;
    a.d[0] = 1; a.d[1] = 2; a.d[2] = 3; a.d[3] = 4;
    b.d[0] = 5; b.d[1] = 6; b.d[2] = 7; b.d[3] = 8;
    c.v = a.v + b.v; 
    printf("%e, %e, %e, %e\n", c.d[0], c.d[1], c.d[2], c.d[3]);
    printf("%e, %e, %e, %e\n", a.d[0]+b.d[0], a.d[1]+b.d[1], a.d[2]+b.d[2], a.d[3]+b.d[3]);
    */
    
    if (detailedlogfile != NULL) fclose(detailedlogfile);
    free(NLevels);
    free(Am); free(Bm); free(Cm);
    if (datafile != NULL) fclose(datafile);
    MPI_Finalize();
    return 0;
}
