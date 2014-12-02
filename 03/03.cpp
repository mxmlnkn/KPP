/*

rm 03.exe; mpic++ 03.cpp -o 03.exe -Wall -std=c++0x; mpirun -n 9 ./03.exe

salloc -n 16
mpic++ 03.cpp -o 03.exe -Wall -std=c++0x
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
do
	srun -n $i ./03.exe
	echo ""
done


Output:

[Rank 0] Sum of ranks on same y-axis: 0
[Rank 1] Sum of ranks on same y-axis: 1

[Rank 2] Sum of ranks on same y-axis: 2
[Rank 1] Sum of ranks on same y-axis: 1
[Rank 0] Sum of ranks on same y-axis: 0

[Rank 0] Sum of ranks on same y-axis: 1
[Rank 2] Sum of ranks on same y-axis: 5

[Rank 1] Sum of ranks on same y-axis: 1
[Rank 4] Sum of ranks on same y-axis: 4
[Rank 3] Sum of ranks on same y-axis: 3
[Rank 2] Sum of ranks on same y-axis: 2
[Rank 0] Sum of ranks on same y-axis: 0

[Rank 4] Sum of ranks on same y-axis: 9
[Rank 0] Sum of ranks on same y-axis: 1
[Rank 2] Sum of ranks on same y-axis: 5

[Rank 5] Sum of ranks on same y-axis: 5
[Rank 4] Sum of ranks on same y-axis: 4
[Rank 1] Sum of ranks on same y-axis: 1
[Rank 6] Sum of ranks on same y-axis: 6
[Rank 2] Sum of ranks on same y-axis: 2
[Rank 3] Sum of ranks on same y-axis: 3
[Rank 0] Sum of ranks on same y-axis: 0

[Rank 2] Sum of ranks on same y-axis: 5
[Rank 0] Sum of ranks on same y-axis: 1
[Rank 6] Sum of ranks on same y-axis: 13
[Rank 4] Sum of ranks on same y-axis: 9

[Rank 6] Sum of ranks on same y-axis: 21
[Rank 3] Sum of ranks on same y-axis: 12
[Rank 0] Sum of ranks on same y-axis: 3

[Rank 2] Sum of ranks on same y-axis: 5
[Rank 6] Sum of ranks on same y-axis: 13
[Rank 8] Sum of ranks on same y-axis: 17
[Rank 4] Sum of ranks on same y-axis: 9
[Rank 0] Sum of ranks on same y-axis: 1

[Rank 2] Sum of ranks on same y-axis: 2
[Rank 1] Sum of ranks on same y-axis: 1
[Rank 7] Sum of ranks on same y-axis: 7
[Rank 8] Sum of ranks on same y-axis: 8
[Rank 9] Sum of ranks on same y-axis: 9
[Rank 10] Sum of ranks on same y-axis: 10
[Rank 6] Sum of ranks on same y-axis: 6
[Rank 4] Sum of ranks on same y-axis: 4
[Rank 3] Sum of ranks on same y-axis: 3
[Rank 5] Sum of ranks on same y-axis: 5
[Rank 0] Sum of ranks on same y-axis: 0

[Rank 9] Sum of ranks on same y-axis: 30
[Rank 3] Sum of ranks on same y-axis: 12
[Rank 6] Sum of ranks on same y-axis: 21
[Rank 0] Sum of ranks on same y-axis: 3

[Rank 2] Sum of ranks on same y-axis: 2
[Rank 1] Sum of ranks on same y-axis: 1
[Rank 9] Sum of ranks on same y-axis: 9
[Rank 10] Sum of ranks on same y-axis: 10
[Rank 12] Sum of ranks on same y-axis: 12
[Rank 11] Sum of ranks on same y-axis: 11
[Rank 8] Sum of ranks on same y-axis: 8
[Rank 4] Sum of ranks on same y-axis: 4
[Rank 5] Sum of ranks on same y-axis: 5
[Rank 6] Sum of ranks on same y-axis: 6
[Rank 3] Sum of ranks on same y-axis: 3
[Rank 0] Sum of ranks on same y-axis: 0
[Rank 7] Sum of ranks on same y-axis: 7

[Rank 2] Sum of ranks on same y-axis: 5
[Rank 10] Sum of ranks on same y-axis: 21
[Rank 4] Sum of ranks on same y-axis: 9
[Rank 6] Sum of ranks on same y-axis: 13
[Rank 12] Sum of ranks on same y-axis: 25
[Rank 8] Sum of ranks on same y-axis: 17
[Rank 0] Sum of ranks on same y-axis: 1

[Rank 12] Sum of ranks on same y-axis: 39
[Rank 3] Sum of ranks on same y-axis: 12
[Rank 6] Sum of ranks on same y-axis: 21
[Rank 0] Sum of ranks on same y-axis: 3
[Rank 9] Sum of ranks on same y-axis: 30

[Rank 12] Sum of ranks on same y-axis: 54
[Rank 4] Sum of ranks on same y-axis: 22
[Rank 0] Sum of ranks on same y-axis: 6
[Rank 8] Sum of ranks on same y-axis: 38

*/

#include <mpi.h>
#include <iostream>
#include <cstdlib>  // malloc



using namespace std;

#define DEBUG_FINDMAXAREA 0
#define DEBUG 0

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
    const int NDimensions = 2;           // 2D-Cartesian Topology
    const int periodic[]  = {true,true}; // Make Torus through periodic boundaries
          int dimSizes[NDimensions];
    findMinCircumference( worldsize, NDimensions, dimSizes ); // does the same as MPI_Dims_create

    MPI_Cart_create(MPI_COMM_WORLD, NDimensions, dimSizes, periodic, true , &commTorus);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* if reorder argument of MPI_Cart_create == false, then this should also *
     * be callable from prior that command and still get the correct rank     */

    int coords[2];
    /* Get Process coords from rank (I find it actually easier, to implement  *
     * such a linear->N-D mapping myself, therefore I really don't see the    *
     * pros in using MPI internal functions aside from code uniformity        */
    MPI_Cart_coords( commTorus, rank, NDimensions, coords);

#if DEBUG == 1
    /* Pass Beacon to let every thread output its value in order */
    {char beacon = 0;
    MPI_Barrier( MPI_COMM_WORLD );
    if (rank == 0) {
        /* Message Body Begin */
        cout << "[Rank:" << rank << "] <-> (";
        for (int j=0; j<NDimensions-1; j++)
            cout << coords[j] << ",";
        cout << coords[NDimensions-1] << ")" << endl;
        /*  Message Body End  */
        MPI_Send( &beacon, 1, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD );
    } else {
        MPI_Recv( &beacon, 1, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        /* Message Body Begin */
        cout << "[Rank:" << rank << "] <-> (";
        for (int j=0; j<NDimensions-1; j++)
            cout << coords[j] << ",";
        cout << coords[NDimensions-1] << ")" << endl;
        /*  Message Body End  */
        if ( rank < worldsize-1 )
            MPI_Send( &beacon, 1, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD );
    }
    cout << flush;
    MPI_Barrier( MPI_COMM_WORLD );}
#endif

    int neighbors[2*NDimensions];
    /* Shift axis 0 upward by 1, commTorus already has saved our rank! */
    MPI_Cart_shift ( commTorus, 0, 1, &neighbors[0], &neighbors[1]);
    /* Shift axis 1 upward by 1 */
    MPI_Cart_shift ( commTorus, 1, 1, &neighbors[2], &neighbors[3]);
#if DEBUG == 1
    printf( "rank= %d coords= %d %d  neighbors(u,d,l,r)= %d %d %d %d\n",
            rank,coords[0],coords[1],neighbors[0],neighbors[1],neighbors[2],
            neighbors[3]);
#endif

    /* Example Output let's us derive direction names                         *
     *   [Rank:0] <-> (0,0)                                                   *
     *   [Rank:1] <-> (0,1)          -----------                              *
     *   [Rank:2] <-> (0,2)         | 2 | 5 | 8 |                             *
     *   [Rank:3] <-> (1,0)         |-----------|                             *
     *   [Rank:4] <-> (1,1)      => | 1 | 4 | 7 |                             *
     *   [Rank:5] <-> (1,2)         |-----------|                             *
     *   [Rank:6] <-> (2,0)         | 0 | 3 | 6 |                             *
     *   [Rank:7] <-> (2,1)          -----------                              *
     *   [Rank:8] <-> (2,2)                                                   *
     *   rank= 7 coords= 2 1  neighbors(u,d,l,r)= 4 1 6 8                     *
     *   rank= 0 coords= 0 0  neighbors(u,d,l,r)= 6 3 2 1                     *
     *   rank= 3 coords= 1 0  neighbors(u,d,l,r)= 0 6 5 4                     *
     *   rank= 1 coords= 0 1  neighbors(u,d,l,r)= 7 4 0 2                     *
     *   rank= 2 coords= 0 2  neighbors(u,d,l,r)= 8 5 1 0                     *
     *   rank= 5 coords= 1 2  neighbors(u,d,l,r)= 2 8 4 3                     *
     *   rank= 6 coords= 2 0  neighbors(u,d,l,r)= 3 0 8 7                     *
     *   rank= 8 coords= 2 2  neighbors(u,d,l,r)= 5 2 7 6                     *
     *   rank= 4 coords= 1 1  neighbors(u,d,l,r)= 1 7 3 5                     *
     * => So source basically contains the inverse of the shift in a certain  *
     *    direction. And shift=-1 basically switches source and destination   */

    const int left   = neighbors[0];
    const int right  = neighbors[1];
    const int down   = neighbors[2];
    const int up     = neighbors[3];
    
    const int xcoord = coords[0];
    const int ycoord = coords[1];

    /* Cycle on y-axis, initiated by y=0 threads */
    int beacon;
    MPI_Barrier( commTorus );
    if ( ycoord == 0 )
    {
        beacon = rank;
        MPI_Send( &beacon, 1, MPI_INT, up  , xcoord, commTorus );
        MPI_Recv( &beacon, 1, MPI_INT, down, xcoord, commTorus, MPI_STATUS_IGNORE );
        printf( "[Rank %i] Sum of ranks on same y-axis: %i\n", rank, beacon );
    }
    else
    {
        MPI_Recv( &beacon, 1, MPI_INT, down, xcoord, commTorus, MPI_STATUS_IGNORE );
        beacon += rank;
        MPI_Send( &beacon, 1, MPI_INT, up  , xcoord, commTorus );
    }
    
    /* Example output following the one in the comment above *
     *   [Rank 0] Sum of ranks on same y-axis: 3             *
     *   [Rank 3] Sum of ranks on same y-axis: 12            *
     *   [Rank 6] Sum of ranks on same y-axis: 21            *
     *      => correct :)                                    */

    MPI_Finalize();
}