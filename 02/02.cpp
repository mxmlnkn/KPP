// rm 02.exe; mpic++ 02.cpp -o 02.exe -Wall -std=c++0x; mpirun --bind-to core -n 2 ./02.exe; python plot.py
// Learned: MPI_STATUS_IGNORE, MPI_ANY_TAG, Communicator, MPI_STATUS, MPI_Get_count
// --bind-to-core necessray because of rtdsc, but this option is also default since OpenMPI 1.7.4
// with 4 threads it takes around 1m, so 100 threads should be possible, because the time also decreases with more threads

//??? wie viel wissen darf ich bei der integration reinstecken?... Integrationsgrenzen, Auch noch finden => Broadcast ???

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

const uint32_t N_LOG_VALUES    = 200;
const uint32_t N_MAX_INTERVALS = 1e4;
const uint32_t LIN_UP_TO       = 500;
const uint32_t N_REPETITIONS   = 100;

#define INT_PRECISION double
#define VERBOSE 1

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
            #if VERBOSE == 1
                cout << val;
            #endif
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


/* if this is 1 files will be flushed for in-situ analysis and other feedback *
 * will be printed on the console                                             */

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


class Function1D {
public:
    Function1D(void) {}
    ~Function1D(void) {}
    virtual double operator()(double x) const = 0;
};
class func : public Function1D {
public:
    func(void) {}
    ~func(void) {}
    double operator()(double x) const {
        return 3*sqrt(x) - x*x + 4*x - 6;
    }
    double antiderivative(double x) const {
        return 2*x*sqrt(x) - x*x*x/3 + 2*x*x - 6*x;
    }
} f;

/* f must be a reference or else there are problems with abstract classes */
INT_PRECISION IntCenter ( const Function1D & f, const INT_PRECISION a, 
                          const INT_PRECISION b, const uint64_t N )
{
    /* N - Number of Intervals, meaning we need N+1 sampling points in 1D */
    if (N==0)
        return 0;
    const INT_PRECISION dx  = (b-a)/N;
    INT_PRECISION       sum = 0;
    for (uint64_t i=0; i<N; i++) {
        const INT_PRECISION xa = a + dx*i;
        const INT_PRECISION xb = a + dx*(i+1);
        sum += f( (xa+xb)/2 );
    }
    return sum * dx;
}

INT_PRECISION IntCenterNew ( const Function1D & f, const INT_PRECISION a, 
                          const INT_PRECISION b, const uint64_t N )
{
    /* N - Number of Intervals, meaning we need N+1 sampling points in 1D */
    if (N==0)
        return 0;
    const INT_PRECISION dx  = (b-a)/N;
    INT_PRECISION       sum = 0;
    for (uint64_t i=0; i<N; i++)
        sum += f( a + dx*(i+0.5) );
    return sum * dx;
}

INT_PRECISION IntTrapeze (const Function1D & f, const INT_PRECISION a, 
                          const INT_PRECISION b, const uint64_t N )
{
    /* N - Number of Intervals */
    if (N==0)
        return 0;
    const INT_PRECISION dx  = (b-a)/N;
    INT_PRECISION       sum = 0;
    for (uint64_t i=0; i<N; i++)
    {
        INT_PRECISION xa=a+dx*i;      //left sampling point
        INT_PRECISION xb=a+dx*(i+1);  //right sampling point

        sum += (f(xa)+f(xb))/2;
    }
    return sum * dx;
}

INT_PRECISION IntTrapezeNew (const Function1D & f, const INT_PRECISION a, 
                          const INT_PRECISION b, const uint64_t N )
{
    /* N - Number of Intervals */
    if (N==0)
        return 0;
    const INT_PRECISION dx  = (b-a)/N;
    INT_PRECISION       sum = 0;
    for (uint64_t i=0; i<=N; i++)
        sum += f( a + dx*i );
    return ( sum - 0.5*( f(a) + f(b) ) ) * dx;
}

int main(int argc, char* argv[])
{
    /* Initialize MPI */
    MPI_Init(NULL, NULL);
    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
    stringstream ranksr; ranksr << rank;

    tout.open( timesr+string("log_rank_")+(ranksr.str())+string(".txt") );
    tout << "Rank: " << rank << "/" << size << "Processor Name: " << processor_name << endl;

    /* Measure CPU Speed and save in file */
    //FILE* statfile = fopen( ( timesr + string("stats") + string("_rank_")
    //                          + ranksr.str() + string(".txt") ).c_str(), "w" );
    double   t0   = MPI_Wtime();
    uint64_t nc0  = GetCPUCycle();
    if (rank == 0) tout << "Measuring CPU Time ." << flush;
    for (uint16_t i=0; i<1; i++) {
        sleep( 1 ); // wait 1 second
        if (rank == 0) tout << "." << flush;
    }
    if (rank == 0) tout << "Done\n";
    uint64_t nc1     = GetCPUCycle();
    double   t1      = MPI_Wtime();
    double   cpufreq = (nc1-nc0)/(t1-t0);
    //if (rank == 0) tout << "MPI_Wtime resolution: " << MPI_Wtick() << endl;
    if (rank == 0) printf( "CPU Frequency : %lu Cycles / %e = (%e +- %e) GHz\n",
            nc1-nc0, t1-t0, cpufreq/1e9, cpufreq/1e9 * MPI_Wtick() / (t1-t0) );
    //fprintf( statfile, "# CPU Frequency / Hz\tsigma / Hz\n" );
    //fprintf( statfile, "%e\t%e\n", cpufreq, cpufreq * MPI_Wtick() / (t1-t0) );
    //fclose( statfile );

    /* Find Integration Interval with Newton's Method (in parallel?) -> N-Sektion */
    func f;
    double a = 1, b = 4;
    double correct_int = f.antiderivative(b) - f.antiderivative(a);
    if (rank==0) tout << "Correct Result: " << correct_int << endl << flush;

    /* Open File for relative Errors ( integration points ) */
    FILE * datafile = NULL;
    if (rank == 0) {
        datafile = fopen( (timesr+string("data.txt")).c_str(), "w" );
        fprintf( datafile, "# integration intervals\tprocessors used\ttime / "
                 "CPU clocks\tsigma CPU clocks\tresult\trelative Error\n" );
        fflush(datafile);
    }

for (uint32_t worldsize = 1; worldsize <= uint32_t(size); worldsize++) {
    tout << "Worldsize: " << worldsize << endl;
    MPI_Barrier( MPI_COMM_WORLD );

    /* Create new subset Communicator */
    MPI_Comm subcomm;

    if (worldsize > 1 && worldsize < uint32_t(size)) {
        MPI_Group worldgroup, subgroup;
        MPI_Comm_group(MPI_COMM_WORLD, &worldgroup); // get old group

        int * ranks = (int*)malloc( sizeof(int)*size );
        for (int i=0; i<size; i++)
            ranks[i] = i;

        if ( uint32_t(rank) < worldsize )
            MPI_Group_incl( worldgroup, worldsize, ranks, &subgroup);
        else
            MPI_Group_incl( worldgroup, size-worldsize, &(ranks[worldsize]), &subgroup);

        // all threads of MPI_COMM_WORLD need to call this routine !
        MPI_Comm_create(MPI_COMM_WORLD, subgroup, &subcomm);
        free(ranks);
    } else if ( worldsize == uint32_t(size) )
        subcomm = MPI_COMM_WORLD;
    if ( uint32_t(rank) >= worldsize ) continue;

    const double N_MULT = pow( N_MAX_INTERVALS - LIN_UP_TO, (double) 1. / N_LOG_VALUES );
    /**************** Loop through different interval numbers *****************/
    for (uint64_t kN = 1; kN <= LIN_UP_TO + N_LOG_VALUES; kN++) {
        uint64_t N;
        if ( kN > LIN_UP_TO )
            N = LIN_UP_TO + floor( pow( N_MULT, kN - LIN_UP_TO ) );
        else
            N = kN;
        //if (rank == 0) tout << "Integration Intervals: " << N << endl << flush;

        double   dx, loc_a, loc_b, loc_int, glob_int;
        uint64_t loc_N, modulo, N_before_me, DBG_glob_N, t0, t1, sum=0, sum2=0;

        for (uint32_t t = 0; t < N_REPETITIONS; t++) {
            /* wait a bit for every time measurement trying to reduce variations */
            t0 = GetCPUCycle();
            while ( GetCPUCycle() - t0 < 1e4 ) {}
            
            if (worldsize == 1) {
                t0       = GetCPUCycle();
                glob_int = IntTrapezeNew( f, a, b, N );
                t1       = GetCPUCycle();
            } else {
                MPI_Barrier( subcomm );
                t0       = GetCPUCycle();

                /**************************************************************
                 * Distribute integration intervals between threads:          *
                 *   e.g. modulo == 2 then rank 0 and 1 receive extra work    *
                 *        integration intervals before a certain rank:        *
                 *        rank 0: 0, rank 1: loc_N+1, rank 2: 2*(loc_N+1),    *
                 *        rank 3: 2*(loc_N+1) + loc_N                         *
                 **************************************************************/
                loc_N       = uint32_t(N) / uint32_t(worldsize);
                modulo      = N - loc_N * worldsize;
                assert( loc_N * worldsize <= N );
                loc_N       = loc_N + ( uint64_t(rank) < modulo );
                N_before_me = ( uint64_t(rank) < modulo ) ? rank * loc_N
                                                          : rank*loc_N + modulo;

                #if DEBUG == 1
                    MPI_Allreduce( &loc_N, &DBG_glob_N, 1, MPI_UINT64_T, MPI_SUM, subcomm );
                    if ( !(DBG_glob_N == N) )
                        tout << "DBG_glob_N: " << DBG_glob_N << " != " << N << " :N" << endl;
                    assert( DBG_glob_N == N );
                #endif

                dx      = (b-a)/N;
                loc_a   = a +  N_before_me       *dx;
                loc_b   = a + (N_before_me+loc_N)*dx;

                loc_int = IntTrapezeNew( f, loc_a, loc_b, int(ceil(N/worldsize)) );
                MPI_Reduce( &loc_int, &glob_int, 1, MPI_DOUBLE, MPI_SUM, 0, subcomm );

                t1 = GetCPUCycle();

                #if DEBUG == 1
                    /* Pass Beacon to let every thread output its value in order */
                    char beacon = 0;
                    tout << flush;
                    MPI_Barrier( subcomm );
                    if (rank == 0) {
                        tout << "[Rank " << rank << "] Local Integration Value: " << loc_int << endl;
                        MPI_Send( &beacon, 1, MPI_CHAR, rank+1, 0, subcomm );
                    } else {
                        MPI_Recv( &beacon, 1, MPI_CHAR, rank-1, 0, subcomm, MPI_STATUS_IGNORE );
                        tout << "[Rank " << rank << "] Local Integration Value: " << loc_int << endl;
                        if ( uint32_t(rank) < worldsize-1 )
                            MPI_Send( &beacon, 1, MPI_CHAR, rank+1, 0, subcomm );
                    }
                    tout << flush;
                    MPI_Barrier( subcomm );
                #endif
            }
            sum  += t1-t0;
            sum2 += (t1-t0)*(t1-t0);
        }
        //integration intervals\tprocessors used\ttime / CPU clocks\tsigma CPU clocks\tresult\trelative Error
        if (rank == 0) {
            fprintf( datafile, "%lu\t%u\t%e\t%e\t%e\t%e\n", N, uint32_t(worldsize),
            double(sum) / N_REPETITIONS / cpufreq,
            sqrt( double(sum2-sum)/(N_REPETITIONS*(N_REPETITIONS-1)) ) /cpufreq,
            glob_int, abs(glob_int - correct_int) / correct_int );
            fflush(datafile);
        }
    } // loop over intervall numbers
} // loop over worldsize

    MPI_Finalize();
    if (rank == 0) if (datafile != NULL) fclose(datafile);
    return 0;
}
