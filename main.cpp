#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include "mpi.h"
#include <sys/mman.h>
#include "umvm.h"
Matrix m;
int main(int argc, char* argv[])
{
    int rank, P;
    int MaxX, MaxY, X, Y;
    int opt;

    // parameters for MPI_Cart_*
    int Dimensions[2] = {-1, -1};
    int Periods[2] = {1, 1}; 
    int CartesianCoords[2];
    
    MaxX = MaxY =  -1;
 
    MPI_Comm Cartesian;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

 
    while((opt = getopt(argc, argv, "X:Y:")) != -1)
        switch(opt)
        {
        case 'X':
            MaxX = strtoul(optarg, NULL, 0);
            break;
        case 'Y':
            MaxY = strtoul(optarg, NULL, 0);
            break;
        default:
            if(rank == 0) printf( "Unknown option. \n");
            
        }
   if((MaxX == -1)||(MaxY == -1))
    {
        if(rank == 0) printf("Please, specify grid size by -X and -Y options.\n");
        MPI_Finalize();
        return 0;
    }
   
    if( MaxX * MaxY != P)
    {
        if(rank == 0) printf("X*Y should be equal to the number of processes.\n");
        MPI_Finalize();
        return 0;
    }
    
    Dimensions[0] = MaxX;
    Dimensions[1] = MaxY;
    
    if(MPI_Cart_create(MPI_COMM_WORLD, 2, Dimensions, Periods, 1, &Cartesian) != MPI_SUCCESS)
        if(rank == 0) printf("Failed to create Cartesian communicator\n"); 
    
    
    MPI_Cart_coords(Cartesian, rank, 2, CartesianCoords);
    printf("Hello Matrix! My coords are %d, %d \n", CartesianCoords[0], CartesianCoords[1]);
    
    MPI_Finalize();
    return 1;
}
