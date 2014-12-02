/*

rm 03-general.exe; mpic++ 03-general.cpp -o 03-general.exe -Wall -std=c++0x; mpirun -n 9 ./03-general.exe

salloc -n 16
mpic++ 03-general.cpp -o 03-general.exe -Wall -std=c++0x
srun -n 16 ./03-general.exe

Output for COMM_DIRECTION == 0 or 1
    [Rank 0] Sum of ranks on same y-axis: 24
    [Rank 1] Sum of ranks on same y-axis: 28
    [Rank 2] Sum of ranks on same y-axis: 32
    [Rank 3] Sum of ranks on same y-axis: 36

Output for COMM_DIRECTION == 2 or 3
    [Rank 0] Sum of ranks on same y-axis: 6
    [Rank 4] Sum of ranks on same y-axis: 22
    [Rank 8] Sum of ranks on same y-axis: 38
    [Rank 12] Sum of ranks on same y-axis: 54

          -----------------
         | 3 | 7 | 11 | 15 | 36
         |---+---+----+----|
     =>  | 2 | 6 | 10 | 14 | 32
         |---+---+----+----|
         | 1 | 5 |  9 | 13 | 28
         |---+---+----+----|
         | 0 | 4 |  8 | 12 | 24
          -----------------
           6  22   38   54
*/

#include <mpi.h>
#include <iostream>
#include <cstdlib>  // malloc

#define COMM_DIRECTION    0  // (LEFT,RIGHT,DOWN,UP) = (0,1,2,3)
#define BITMASK_AXIS      uint8_t(256+128+64+32+16+8+4+2)
#define BITMASK_DIRECTION 1

using namespace std;

#define DEBUG_FINDMAXAREA 0
#define DEBUG 1

void findMinCircumference( int worldsize, int NDimensions, int * dimSizes )
{
    if (worldsize < NDimensions)
    {
        cerr << "Please start at least " << NDimensions << " Threads!" << endl;
        return;
    }

    int * sizesTmp  = reinterpret_cast<int*>( malloc( sizeof(int)*NDimensions ) );
    int   circumTmp = 2*worldsize+2*1;
    for (int j=0; j<NDimensions; j++)
        sizesTmp[j] = 1;
    sizesTmp[0] = worldsize;
    memcpy( dimSizes, sizesTmp, sizeof(int)*NDimensions );

#if DEBUG_FINDMAXAREA == 1
    cout << endl << "Find Maximum area with " << worldsize << " elements in "
         << NDimensions << " dimensions :" << endl;
    cout << "Initial Circumference: " << circumTmp << endl;
#endif

    /* Brute Force ... */
    while(true) {
        /* Increment 'Key' */
        bool carrier = false;
        int  loc     = 0;
        do {
            sizesTmp[loc]++;

            /* Calc area */
            int circumference = 0;
            int area          = 1;
            for (int j=0; j<NDimensions; j++) {
                area *= sizesTmp[j];
                circumference += 2*sizesTmp[j];
                /* this is not correct for higher than 2D! */
            }

            /* Carrier Flag */
            if (area > worldsize) {
                carrier = true;
                sizesTmp[loc++] = 1;
            } else if ( area == worldsize ) {
                carrier = false;
                if ( circumference < circumTmp) {
                    circumTmp = circumference;
                    memcpy( dimSizes, sizesTmp, sizeof(int)*NDimensions );
#if DEBUG_FINDMAXAREA == 1
                    cout << "Found smaller Circumference of length: "
                         << circumference  << endl << flush;
#endif
                }
            } else {
                carrier = true;
                loc     = 0;
            }

            /* Debug Print out */
#if DEBUG_FINDMAXAREA == 1
            for (int j=0; j<NDimensions; j++)
                cout << sizesTmp[j] << " ";
            cout << " => Area: " << area << " Circumference: " << circumference
                 << endl << flush;
#endif
        } while(carrier and loc < NDimensions);
        /* after this loop the 'key' has been incremented by one >valid< step *
         * or it has reached the last configuration                           */
#if DEBUG_FINDMAXAREA == 1
        cout << "One valid increment done" << endl;
#endif

        if ( loc >= NDimensions )
            break;
    }

    free( sizesTmp );
}

int main( int argc, char **argv ) {

    /* Initialize MPI */
    MPI_Init(NULL, NULL);
    int rank,worldsize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    /* Create 2D-Torus */
    MPI_Comm  commTorus;
    const int NDimensions = 2;
    const int periodic[]  = {true,true};
          int dimSizes[NDimensions];
    findMinCircumference( worldsize, NDimensions, dimSizes ); // does the same as MPI_Dims_create
#if DEBUG == 1
    printf( "Torus Size: (%i,%i)\n", dimSizes[0], dimSizes[1] );
#endif
    MPI_Dims_create( worldsize, NDimensions, dimSizes );
#if DEBUG == 1
    cout << "Dimension calculated by MPI: ";
    for (int j=0; j<NDimensions; j++)
        cout << dimSizes[j] << " ";
    cout << endl;
#endif

    MPI_Cart_create(MPI_COMM_WORLD, NDimensions, dimSizes, periodic, true , &commTorus);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int coords[NDimensions];
    MPI_Cart_coords( commTorus, rank, NDimensions, coords);
    
    int sourcerank, destrank;
    int axis      = (COMM_DIRECTION & BITMASK_AXIS) >> 1;
    int direction = (COMM_DIRECTION & BITMASK_DIRECTION);
    if (direction == 0) direction = -1;
    MPI_Cart_shift ( commTorus, axis, direction, &sourcerank, &destrank);
#if DEBUG == 1
    printf( "[Rank %i] coord=(%i,%i) => src:%i dest:%i\n", rank, coords[0], coords[1], sourcerank, destrank );
#endif

    /* Cycle on y-axis, initiated by y=0 threads */
    int beacon;
    MPI_Barrier( commTorus );
    if ( coords[axis] == 0 )
    {
        beacon = rank;
        MPI_Send( &beacon, 1, MPI_INT, destrank  , 123, commTorus );
        MPI_Recv( &beacon, 1, MPI_INT, sourcerank, 123, commTorus, MPI_STATUS_IGNORE );
        printf( "[Rank %i] Sum of ranks on same y-axis: %i\n", rank, beacon );
    }
    else
    {
        MPI_Recv( &beacon, 1, MPI_INT, sourcerank, 123, commTorus, MPI_STATUS_IGNORE );
        beacon += rank;
        MPI_Send( &beacon, 1, MPI_INT, destrank  , 123, commTorus );
    }

    MPI_Finalize();
}