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
    
   Matrix Block;
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
  
    GenerateMatrix(Block, Cartesian, P, MaxX, MaxY, Weight, Weight, N, M);    
    MPI_Finalize();
    return 1;
}
