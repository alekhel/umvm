#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>
#include "umvm.h"
int main(int argc, char* argv[])
{
    int rank, P;
    int MaxX, MaxY, X, Y; //Grid parameters and coords
    unsigned int N, M; //Matrix height and width
    unsigned int Weight; // Average row weight    
    int opt;
    
    int H;
    double StartTime, EndTime; 
    // parameters for MPI_Cart_*
    int Dimensions[2] = {-1, -1};
    int Periods[2] = {1, 1}; 
    int CartesianCoords[2];
    Matrix Strip, Columns, Block;
    MaxX = MaxY = N = M = Weight = 0;
 
    MPI_Comm Cartesian;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

 
    while((opt = getopt(argc, argv, "X:Y:W:M:N:")) != -1)
        switch(opt)
        {
        case 'X':
            MaxX = strtoul(optarg, NULL, 0);
            break;
        case 'Y':
            MaxY = strtoul(optarg, NULL, 0);
            break;
        case 'W':
            Weight = strtoul(optarg, NULL, 0);
            break;
        case 'N':
            N = strtoul(optarg, NULL, 0);
            break;
        case 'M':
            M = strtoul(optarg, NULL, 0);
            break;
         default:
            if(rank == 0) printf( "Unknown option. \n");
            
        }
    if(ParameterSanityCheck(P, MaxX, MaxY, N, M, Weight))
        return -1;

    Dimensions[0] = MaxX;
    Dimensions[1] = MaxY;
    if(MPI_Cart_create(MPI_COMM_WORLD, 2, Dimensions, Periods, 0, &Cartesian) != MPI_SUCCESS)
        if(rank == 0) printf("Failed to create Cartesian communicator\n"); 
    MPI_Cart_coords(Cartesian, rank, 2, CartesianCoords);
   MPI_Finalize();
    return 1;
}
