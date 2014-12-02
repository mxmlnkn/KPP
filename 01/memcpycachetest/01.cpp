// rm 01.exe; mpic++ 01.cpp -o 01.exe -Wall -std=c++0x; mpirun --bind-to-core -n 2 ./01.exe; python plot.py
// Learned: MPI_STATUS_IGNORE, MPI_ANY_TAG, Communicator, MPI_STATUS, MPI_Get_count
// --bind-to core necessesary because of rtdsc, but this option is also default since OpenMPI 1.7.4


/* todo: cronjob
zwischen 1200 und 12300 byte nachrichtenlänge genauer schauen!
    sowieso nochmal 2 cores auf einem durchführen ...
    
    nachrichtenlänge aus array nehmen, den man vorher anlegt
    voraussichtliche restzeit ausreichen und anzeigen (in abhängigkeit von 10gb/s)
    oder halt einfach prozent + aktuelle zeit anzeigen, damit ichs selber timen kann
    */
    
/* Datenpunkte:
bis 300mb hoch bzw. 1e8 statt 1e7 (100mb)

Heim-PC
		L1-Cache: 2 x 32 KB instruction caches
				  2 x 32 KB data caches
		L2-Cache: 2 x 256 KB
        L3-Cache 3MB
        
Node:
    L1: 8 x 32 KB data caches
    L2: 8 x 256 KB
    L3: 20 MB

*/

/*
Redo with:
    Openmpi 1.8.3
    bind-to core
    exclusive
    logging of node-name!
    more data points at certain areas, see above
    
ToDo:
    Memcpy muss nur ein Thread machen!
    MPI_Wtime ist nutzlos wenn nicht kumulativ, das weißn ich ja nun -> kleinere Dateien
    
*/

#include <iostream>
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


/******************************************************************************
 * This Benchmark sends messages of varying length and measures the time.     *
 * Parameters:                                                                *
 *    N_REPETITIONS - Number of times the time is being measured              *
 *    MSG_MAX_SIZE  - Maximal message length to be measured                   *
 *    BEGIN_WITH_LONGEST_MSG - this tests the times of messages with diff.    *
 *        lengths, but begins with the longest. The reasoning here is, that   *
 *        data sent already could be cached for it may be accessed later on   *
 *    LIN_UP_TO     - Up to this message length all lengths will be tested    *
 *    N_LOG_VALUES  - Number of Values after LIN_UP_TO to be measured.        *
 *                    These Values will be spaced logarithmically, e.g.:      *
 *                    100, 1000, 10000 Bytes                                  *
 *    USE_MEMCPY    - Use memcpy to benchmark the RAM itself and not MPI      *
 *    N_OFFSET      - The first N_OFFSET time measurements will be ignored    *
 *                    in order to get stationary measurements                 *
 *    PREVENT_CACHING - the buffer will be overwritten with random data.      *
 *         you can't cache something unknown, except if you cache it while    *
 *         it is being created                                                *
 *    NTH_SAVE_DETAILED - save every THIS the whole N_REPETITIONS measure-    *
 *         ments not being evaluated statistically yet                        *
 ******************************************************************************/
void Benchmark( int rank, volatile uint8_t * msg, volatile uint8_t * msg2, string filename,
    bool PREVENT_CACHING, bool USE_MEMCPY, bool BEGIN_WITH_LONGEST_MSG,
    uint16_t N_REPETITIONS = 100, uint16_t N_OFFSET = 1,
    uint64_t MSG_MAX_SIZE = 1e8, uint64_t LIN_UP_TO = 10,
    uint64_t N_LOG_VALUES = 1000, uint64_t NTH_SAVE_DETAILED = 2000 )
{
    const float N_MULT = pow( MSG_MAX_SIZE, (float) 1. / N_LOG_VALUES );
    srand( GetCPUCycle() );
    
    /* Measure CPU Frequency, we will need this later*/
    double   t0   = MPI_Wtime();
    uint64_t nc0  = GetCPUCycle();
    sleep( 1 ); // wait 1 second
    uint64_t nc1  = GetCPUCycle();
    double   t1   = MPI_Wtime();
    double   freq = (nc1-nc0)/(t1-t0);

    FILE * logfile         = NULL;
    FILE * detailedlogfile = NULL;
    FILE * logcumulative   = NULL;
    if (rank == 0) {
        logfile = fopen( (filename+string(".txt")).c_str(), "w" );
        fprintf( logfile, "# message length / bytes\ttime / CPU clocks\tsigma CPU clocks\t MPI_Wtime\tsigma MPI_Wtime\n" );
        detailedlogfile = fopen( (filename+string("_Times.txt")).c_str(), "w" );
        fprintf( detailedlogfile, "#Length / Byte\tCPU clocks 1\tCPU clocks 2\t...\n" );
        if (!PREVENT_CACHING && !USE_MEMCPY) {
            logcumulative = fopen( (filename+string("_Cumulative.txt")).c_str(), "w" );
            fprintf( logcumulative, "# message length / bytes\tMPI_Wtime / s\n" );
        }
    }

    uint16_t k = BEGIN_WITH_LONGEST_MSG ? N_LOG_VALUES + LIN_UP_TO : 0;
    while(true) {
        /* Get the next message length to test. Stop if longer than permitted */
        k = BEGIN_WITH_LONGEST_MSG ? k-1 : k+1;
        uint64_t length = 1;
        if ( k > LIN_UP_TO ) {
            length = floor( pow( N_MULT, k-LIN_UP_TO ) );
            if (length < LIN_UP_TO) // exclude values already calculated
                continue;
        } else
            length = k;
        if (length > MSG_MAX_SIZE || length == 0)
            break;

        if ( (rank == 0) && (k % NTH_SAVE_DETAILED == 0) ) {
            fprintf( detailedlogfile, "%lu\t", length );
        }
        /* some temporary variables for doing statistics later on             */
        uint64_t nc_sum = 0, nc_sum2 = 0;
        double   t_sum  = 0, t_sum2  = 0;

        /* Rank 0 sends random message to rank 1. As MPI_Send blocks the      *
         * program until it knows that the message was received, time         *
         * measurement is trivial. This should measure basically ping         *
         * time + transmission protocol (TCP-IP?) overhead.                   */
        for (uint16_t i = N_REPETITIONS + N_OFFSET; i > 0; i--) {
            uint64_t nc0, nc1; // number of cpu cycles
            double   t0, t1;   // time measured with MPI_Wtime

            /* before send fill buffer with random numbers to prevent caching */
            if (PREVENT_CACHING)
                for (uint32_t i = 0; i < length; i++)
                    msg[i] = uint8_t(rand() % 256);
            MPI_Barrier( MPI_COMM_WORLD );
            if ( rank == 1 && !USE_MEMCPY )
                MPI_Recv( const_cast<uint8_t *>(msg), length, MPI_BYTE, 
                          rank-1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
                          
            // All Threads to memcpy in parallel
            if (USE_MEMCPY) {
                    nc0 = GetCPUCycle();
                    memcpy ( const_cast<uint8_t *>(msg) ,
                             const_cast<uint8_t *>(msg2), length );
                    nc1 = GetCPUCycle();
            } else if (rank == 0) {
                /* sleep 1ms until receiving mechanism is set up => Benchmark *
                 * takes 0.001*N_REPETITIONS*(N_LIN_UP_TO+N_LOG_VALUES) longer*/
                    t0  = MPI_Wtime(); while( MPI_Wtime()-t0 < 0.001 ) {}
                    nc0 = GetCPUCycle();
                    MPI_Send( const_cast<uint8_t *>(msg), length, MPI_BYTE,
                              rank+1, i, MPI_COMM_WORLD);
                    nc1 = GetCPUCycle();
            }
            // only use results after some initial tests
            if (i <= N_REPETITIONS) {
                nc_sum2 += (nc1-nc0)*(nc1-nc0);
                nc_sum  += (nc1-nc0);
            }
            if ( (rank == 0) && (k % NTH_SAVE_DETAILED == 0) ) {
                fprintf( detailedlogfile, "%lu\t", nc1-nc0 );
            }
        }

        const uint16_t N = N_REPETITIONS;
        float nc_mean    = float(nc_sum) / float(N);
        float nc_stddev  = sqrt( double(nc_sum2-nc_sum)/(N*(N-1)) );
        
        if ( rank == 0 ) {
            if (k % NTH_SAVE_DETAILED == 0) {
                fprintf( detailedlogfile, "\n" );
                fflush ( detailedlogfile );
            }
            fprintf( logfile, "%lu\t%f\t%f\n", length, nc_mean, nc_stddev);
            fflush(logfile);
#if VERBOSE == 1
            if ( k % NTH_SAVE_DETAILED == 0 ) {
                tout << length << " bytes in " << nc_mean << " +- " << nc_stddev << " clocks\n";
            }
#endif
        }
        
        

        /* Do the time measurements accumulative */
        if (!PREVENT_CACHING && !USE_MEMCPY) {
            MPI_Barrier( MPI_COMM_WORLD ); // don't accidentally steal message from prior for loop
            volatile uint64_t repetitions;
            if (rank == 0) {
                const double min_time_to_measure = 0.001; //MPI_Wtick() * 1e3; // measure at least 1000 discrete timesteps
                repetitions = max( uint64_t(100), uint64_t(ceil( min_time_to_measure * freq / nc_mean)) );
                MPI_Send( const_cast<uint64_t *>(&repetitions), 1, MPI_UINT64_T, rank+1, 0, MPI_COMM_WORLD );
            } else if (rank == 1)
                MPI_Recv( const_cast<uint64_t *>(&repetitions), 1, MPI_UINT64_T, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            MPI_Barrier( MPI_COMM_WORLD ); // prevent one thread waiting for the other thread while the time is already being measured
            
            double t0,t1;
            t0 = MPI_Wtime();
            for (uint64_t i = 0; i < repetitions; i++) {
                if ( rank == 1 ) {
                    MPI_Recv( const_cast<uint8_t *>(msg), length, MPI_BYTE,
                              rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send( const_cast<uint8_t *>(msg), length, MPI_BYTE,
                              rank-1, 0, MPI_COMM_WORLD);
                } else if ( rank == 0 ) {
                    MPI_Send( const_cast<uint8_t *>(msg), length, MPI_BYTE,
                              rank+1, 0, MPI_COMM_WORLD);
                    MPI_Recv( const_cast<uint8_t *>(msg), length, MPI_BYTE,
                              rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                // MPI_Barrier( MPI_COMM_WORLD ); // without this, at least on windows more and more and more RAM is being allocated (1.5gb instead of 18mb theoretically ... it seems some messages get lost in the buffer of MPI, I am not sure why ...
            }
            t1  = MPI_Wtime();
            if (rank==0) {
                fprintf( logcumulative, "%lu\t%e\n", length, 
                         (t1-t0)/(2*repetitions) );
                if ( k % NTH_SAVE_DETAILED == 0 ) {
                    tout << setprecision(16) << endl << "t0: " << t0 << "  Time needed: " << t1-t0 << endl << setprecision(5);
                    fflush( logcumulative ); // this hangs sometimes .... not sure why file shouldn't be closed by other thread prematurely
                }
            }
            MPI_Barrier( MPI_COMM_WORLD ); //without this the program sometimes hangs between to couts !!!
        }

    }
    MPI_Barrier( MPI_COMM_WORLD );
    if (rank == 0) {
        /* if other rank closes file while rank 0 is still writing the        *
         * program will hang                                                  */
        if (logfile         != NULL) fclose(logfile);
        if (detailedlogfile != NULL) fclose(detailedlogfile);
        if (logcumulative   != NULL) fclose(logcumulative);
    }
}


int main(int argc, char* argv[]) {
    // comodo Firewall notices communication on 127.0.0.1:56629 (TCP) (port changes)
    MPI_Init(NULL, NULL);
    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // create time prefix string yyyy-mm-dd_
    time_t t = time(0);   // get time now
    struct tm * now = localtime( &t );
    stringstream timestamp;
    timestamp << 1900 + now->tm_year << "-" << 1 + now->tm_mon << "-"
                << now->tm_mday << "_" << now->tm_hour << "-" << now->tm_min
                << "_";
    string timesr = timestamp.str();
    
    stringstream ranksr; ranksr << rank;
    tout.open( timesr+string("log_rank_")+(ranksr.str())+string(".txt") );
    
    tout << "Processor Name: " << processor_name << endl;
    
    // print some infos relevant to timing
    if (rank == 0) {
        FILE* statfile = fopen( (timesr+string("stats.txt")).c_str(), "w" );
        
        double   t0   = MPI_Wtime();
        uint64_t nc0  = GetCPUCycle();
        tout << "Measuring CPU Time ." << flush;
        for (uint16_t i=0; i<10; i++) {
            sleep( 1 ); // wait 1 second
            tout << "." << flush;
        }
        tout << "Done\n";
        uint64_t nc1  = GetCPUCycle();
        double   t1   = MPI_Wtime();
        double   freq = (nc1-nc0)/(t1-t0);
        printf( "MPI_Wtime resolution: %e\n", MPI_Wtick() );
        printf( "CPU Frequency : %lu Cycles / %e = (%e +- %e) GHz\n",
                nc1-nc0, t1-t0, freq/1e9, freq/1e9 * MPI_Wtick() / (t1-t0) );
                
        fprintf( statfile, "# CPU Frequency / Hz\tsigma / Hz\n" );
        fprintf( statfile, "%e\t%e\n", freq, freq * MPI_Wtick() / (t1-t0) );
        fclose( statfile );
    }
    
    
    // create buffer for random data
    const uint32_t msg_size = 1e9;
    volatile uint8_t * msg  = (uint8_t*) malloc( sizeof(uint8_t) * msg_size );
    volatile uint8_t * msg2 = (uint8_t*) malloc( sizeof(uint8_t) * msg_size );

    //void Benchmark( uint8_t * msg, uint8_t * msg2, string filename,
    //    bool PREVENT_CACHING, bool USE_MEMCPY, bool BEGIN_WITH_LONGEST_MSG,
    
    //used RAM as per Windows: ~18mb
    // if (rank == 0) tout << "============== Caching ==============" << endl;
    // Benchmark( rank, msg, msg2, timesr+string("Caching"             ), false, false, false );
    // if (rank == 0) tout << "============== Caching_LongestFirst ==============" << endl;
    // Benchmark( rank, msg, msg2, timesr+string("Caching_LongestFirst"), false, false, true  );
    // if (rank == 0) tout << "============== NoCache_LongestFirst ==============" << endl;
    // Benchmark( rank, msg, msg2, timesr+string("NoCache_LongestFirst"), true , false, true  );
    
    if (rank == 0) tout << "============== Memcpy_Caching ==============" << endl;
    Benchmark( rank, msg, msg2, timesr+string("Memcpy_Caching"             ), false, true, false );
    // if (rank == 0) tout << "============== Memcpy_Caching_LongestFirst ==============" << endl;
    // Benchmark( rank, msg, msg2, timesr+string("Memcpy_Caching_LongestFirst"), false, true, true  );
    // if (rank == 0) tout << "============== Memcpy_NoCache_LongestFirst ==============" << endl;
    // Benchmark( rank, msg, msg2, timesr+string("Memcpy_NoCache_LongestFirst"), true , true, true  );

    free( const_cast<uint8_t*>(msg ) );
    free( const_cast<uint8_t*>(msg2) );
    
    MPI_Finalize();
    return 0;
}
