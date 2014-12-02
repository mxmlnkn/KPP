// rm pingpong_memleak.exe; mpic++ pingpong_memleak.cpp -o pingpong_memleak.exe -Wall -std=c++0x; mpirun -n 2 ./pingpong_memleak.exe

#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    volatile char msg[256];

    for (uint64_t i = 0; i < 1e9; i++) {
        MPI_Send( const_cast<char*>(msg), 256, MPI_CHAR, rank, 0, MPI_COMM_WORLD);
        MPI_Recv( const_cast<char*>(msg), 256, MPI_CHAR, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("%i\n",i);
        MPI_Barrier( MPI_COMM_WORLD );
        for (uint32_t k = 0; k < 256; k++)
            msg[k]++;
    }
    
    MPI_Finalize();
    return 0;
}
