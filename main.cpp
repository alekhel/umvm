#include "umvm.h"
#include "umvm_internals.h"
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>
int main(int argc, char* argv[])
{
    int rank, P, MaxP;
    int MaxX, MaxY, X, Y; //Grid parameters and coords
    Ind N, M; //Matrix height and width
    int Weight; // Average row weight    
    int opt;
    int ITER_H;
    int Type;
    double StartTime, EndTime;
    int UsageType = -1;
    char *Folder, *Prefix;
    Matrix Block;
    MaxX = MaxY = N = M = Weight = Type = ITER_H = 0;
    Folder = Prefix = NULL;
    MPI_Comm Cartesian;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MaxP);

 
    while((opt = getopt(argc, argv, "t:f:p:X:Y:W:M:N:")) != -1)
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
    ITER_H = 1;//TODO
    MaxX = (1<<MaxX);
    MaxY = (1<<MaxY);
    N = (1<<N);
    M = (1<<M);
    if(UsageType == 0)
    {
               if(ParameterSanityCheck(UsageType, MaxP, MaxX, MaxY, N, M, Weight))
            return -1;
  
        GenerateMatrix(Block, Cartesian, MaxP, MaxX, MaxY, Weight, Weight, N, M, ITER_H);    
        
        StoreMatrixToFolder(Folder, Prefix, Block, Type,  MaxP, MaxX, MaxY, Weight, Weight, N, M, Cartesian);    
    }
    if(UsageType == 1)
    {
        
        LoadIterationForThreeProcessTypes(Folder, Prefix, Block, Type,  MaxP, MaxX, MaxY, Weight, Weight, N, M, Cartesian); 
//        PrintMatrixStructure(Block);
        MPI_Barrier(MPI_COMM_WORLD);
        MulBroadcast3(1, Block, Type, Cartesian, MaxP, MaxX, MaxY, N, M, Weight);
    }
    if((UsageType != 0)&&(UsageType!=1)&&(rank == 0))
        printf("Please specify -t option: 0 to generate and store new matrix, 1 to load existing matrix.\n");
  // printf("MyRank is %d, I have %d elements. \n", rank, CountElements(Block)); 
  //  PrintMatrixStructure(Block);
  //  if(rank == 7)
    //    PrintMatrix(Block);
    MPI_Finalize();
    return 1;
}
