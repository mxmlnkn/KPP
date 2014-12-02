//rm 01.exe; mpic++ 01.cpp -o 01.exe -Wall; mpirun -n 2 ./01.exe

#include <iostream>
#include <cstdio>
#include <mpi.h>
#include <ctime>
#include <cctype>
#include <cstdlib>
#include <cmath>

using namespace std;


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



int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    // comodo Firewall notices communication on 127.0.0.1:56629 (TCP) (port changes)
    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Parameters
    const uint32_t N_VALUES      = 6;  // double message size every step 1,2,4,8,16,32 bytes
    const uint16_t N_MULT        = 16; // 1,256,256^2,...
    const uint16_t N_REPETITIONS = 50; // Repeat measurement to reduce statistical errors
    const uint16_t N_OFFSET      = 0;  // discard first N_OFFSET measurements, because the state is not stationary yet
    
    // Open Log File with measured Times
    char filename[256];
    sprintf( filename, "times_rank_%i.txt", rank );
    FILE* logfile = fopen( filename, "w" );
    fprintf( logfile, "# message length/Bytes\tTime Send/clocks\t Time Ping Pong/clocks (clocks per second: %u\n", CLOCKS_PER_SEC );
    
    // create buffer with random data
    const uint32_t msg_size = ceil(pow(N_MULT,N_VALUES));
    uint8_t* msg = (uint8_t*) malloc( sizeof(uint8_t) * ceil(pow(N_MULT,N_VALUES)) );
    printf( "msg pointer: %p\n", (void*) msg );
    srand( GetCPUCycle() );
    for (uint32_t i = 0; i < msg_size; i++)
        msg[i] = uint8_t(rand() % 256);
    
    clock_t  clocks_start, clocks_end;
    uint64_t cycles_start, cycles_end, cycles_diff;
    MPI_Status status;
    
    
    // display buffer in order of the rank number (first rank informs next rank that it has finished writing. Other threads wait for a "go!"-message of their respective prior thread. Without doing this the output could be completely jumbled up. 
    /*int beacon = 51471+rank; //51471 ~ 'start' in leetspeak
    while(true)
        if (beacon == 51471) {
            cout << "[Rank " << rank << "] ";
            for (uint32_t i = 0; i < msg_size; i++) {
                printf("%X",msg[i]);
            }
            printf("\n");
            if (rank < size-1)
                MPI_Send( &beacon, 1, MPI_UINT8_T, rank+1, 51471, MPI_COMM_WORLD);
            break;
        } else
            MPI_Recv( &beacon, 1, MPI_UINT8_T, rank-1, 51471, MPI_COMM_WORLD, &status);
    */
    for (uint32_t length = 1; length < msg_size; length *= N_MULT)
    {
        MPI_Barrier( MPI_COMM_WORLD ); // we don't one sending 256 bytes and the other one waiting for 512 Bytes
        uint64_t sum = 0, sum2 = 0;
        
        for (uint16_t i = 0; i < N_REPETITIONS; i++) {
            // Rank 0 sends random message to rank 1. As MPI_Send blocks the program until it know that the message was received, time measurement is trivial. This should measure basically ping time + transmission protocol (TCP-IP) overhead.
            MPI_Barrier( MPI_COMM_WORLD );
            if (rank == 0) { // rank 0 sends messages
                clocks_start = clock();
                cycles_start = GetCPUCycle();
                MPI_Send( &msg, length, MPI_UINT8_T, rank+1, i, MPI_COMM_WORLD); // standard blocking Send
                cycles_end = GetCPUCycle();
                clocks_end = clock();
            } else { // rank 1 receives messages
                clocks_start = clock();
                cycles_start = GetCPUCycle();
                MPI_Recv( &msg, length, MPI_UINT8_T, rank-1, i, MPI_COMM_WORLD, &status);
                cycles_end = GetCPUCycle();
                clocks_end = clock();
            }
            // only forward results to output after some initial tests
            if (i > N_OFFSET) {
                cycles_diff = cycles_end - cycles_start;
                sum2       += cycles_diff*cycles_diff;
                sum        += cycles_diff;
                fprintf( logfile, "%u\t%u\t%lu\n", length, uint32_t(clocks_end-clocks_start), cycles_diff );
            }
            // after every benchmark fill message buffer with new random numbers to prevent caching
            for (uint32_t i = 0; i < length; i++)
                msg[i] = uint8_t(rand() % 256);
        }
        
        const uint16_t N = N_REPETITIONS - N_OFFSET;
        float mean = float(sum) / float(N);
        float stddev = sqrt( double(sum2-sum)/(N*(N-1)) );
        
        MPI_Barrier( MPI_COMM_WORLD ); // we don't one sending 256 bytes and the other one waiting for 512 Bytes
        if (rank == 0) {
            printf( "[Rank %u] %u Bytes take %.1f +- %.1f cycles ", rank, length, mean, stddev );
            if (rank == 0) printf( "to send\n" );
            else           printf( "to receive\n" );
        }
    }
    
    cout << "[Rank " << rank << "] " << "Clocks/s: " << CLOCKS_PER_SEC << "Cycles consumed by the whole program: " << clock() << endl;
    
    fclose( logfile );
    MPI_Finalize();
    return 0;
}
}