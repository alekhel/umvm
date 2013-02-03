#include "umvm.h"
#include "umvm_internals.h"
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>
int main(int argc, char* argv[])
{
    int rank, P;
    int MaxX, MaxY, X, Y; //Grid parameters and coords
    Ind N, M; //Matrix height and width
    int Weight; // Average row weight    
    int opt;
    int ITER_H;
    int Type;
    int UsageType = -1;
    char *Folder, *Prefix;
    Matrix Block;
    MaxX = MaxY = N = M = Weight = Type = ITER_H = 0;
    Folder = Prefix = NULL;
    MPI_Comm Cartesian;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

 
    while((opt = getopt(argc, argv, "t:f:p:X:Y:W:M:N:I:")) != -1)
        switch(opt)
        {
        case 't':
            UsageType = strtoul(optarg, NULL, 0);
            break;
        case 'f':
            Folder = optarg;
            break;
        case 'p':
            Prefix = optarg;
            break;
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
         case 'I':
            ITER_H = strtoul(optarg, NULL, 0);
            break;
         default:
            if(rank == 0) printf( "Unknown option. \n");
            
        }
    printf("options parsed\n");
    if(UsageType == 0)
    {
        MaxX = (1<<MaxX);
        MaxY = (1<<MaxY);
        N = (1<<N);
        M = (1<<M);
        if(ParameterSanityCheck(P, MaxX, MaxY, N, M, Weight))
            return -1;
  
        GenerateMatrix(Block, Cartesian, P, MaxX, MaxY, Weight, Weight, N, M, ITER_H);    
        StoreMatrixToFolder(Folder, Prefix, Block, Type,  P, MaxX, MaxY, Weight, Weight, N, M, Cartesian);    
    }
    if(UsageType == 1)
    {
        LoadMatrixFromFolder(Folder, Prefix, Block, Type,  P, MaxX, MaxY, Weight, Weight, N, M, Cartesian);    
    }
    
    if((UsageType != 0)&&(UsageType!=1)&&(rank == 0))
        printf("Please specify -t option: 0 to generate and store new matrix, 1 to load existing matrix.\n");
    
    printf("[main] My rank is %d, have %d elements.\n", rank, CountElements(Block));
    MPI_Finalize();
    return 1;
}
