#include <iostream>
#include <ctime>
#include <mpi.h>
#include <unistd.h> // sleep

using namespace std;

int main(int argc, char* argv[]) 
{
    MPI_Init(NULL,NULL);
    cout << "Wtick: " << MPI_Wtick();
    for (int i=0; i<1e7; i++) {
        clock_t   t0  = MPI_Wtime();
        sleep(0.1);
        //sleep(1);     // SLEEP WON'T BE COUNTED FOR CLOCK!!! FUCK, SAME SEEMS TO BE TRUE FOR MPI_WTIME
        // ALSO MPI_WTIME might not be working with -mno-sse!!!
        clock_t   t1    = MPI_Wtime();
        cout << "clock=" << clock() << ", CLOCKS_PER_SEC=" << CLOCKS_PER_SEC;
        cout << ", t0=" << t0 << ", t1=" << t1 << ", dt=" << double(t1-t0)/CLOCKS_PER_SEC << endl;
    }
    return 0;
}

/*
Hypatia@Hypatia-PC /cygdrive/d/Studium/7TH SEMESTER/KPP/MatrixCache
$ mpic++ matrix-diverror.cpp -o matrix-diverror.exe -Wall -mno-mmx -std=c++0x; ./matrix-diverror.exe
clock=0, CLOCKS_PER_SEC=1000, t0=0, t1=0
MPI_Wtime resolution: 1e-06

Hypatia@Hypatia-PC /cygdrive/d/Studium/7TH SEMESTER/KPP/MatrixCache
$ mpic++ matrix-diverror.cpp -o matrix-diverror.exe -Wall -mno-sse -std=c++0x; ./matrix-diverror.exe
clock=0, CLOCKS_PER_SEC=1000, t0=8.37245e-314, t1=8.37245e-314
MPI_Wtime resolution: 0

Hypatia@Hypatia-PC /cygdrive/d/Studium/7TH SEMESTER/KPP/MatrixCache
$ mpic++ -g matrix-diverror.cpp -o matrix-diverror.exe -Wall -std=c++0x -mno-mmx -mno-sse -mno-sse2 -mno-sse3 -mno-sse4; ./matrix-diverror.exe
clock=15, CLOCKS_PER_SEC=1000, t0=8.37245e-314, t1=8.37245e-314
MPI_Wtime resolution: 0

*/