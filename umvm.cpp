#include "umvm.h"
#include "umvm_internals.h"
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "math.h"
int ParameterSanityCheck(int t, int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if((MaxX == 0)||(MaxY == 0))
    {
        if(rank == 0) printf("Please, specify grid size by -X and -Y options.\n");
        MPI_Finalize();
        return -1;
    }
    if((N == 0)||(M == 0))
    {
        if(rank == 0) printf("Please, specify matrix height and width by -N and -M options.\n");
        MPI_Finalize();
        return -1;
    }
    if(Weight == 0)
    {
        if(rank == 0) printf("Please, specify average row weight by -W option.\n");
        MPI_Finalize();
        return -1;
    }
    if(Weight > M)
    {
        if(rank == 0) printf("Average row weight should be less than matrix width.\n");
        MPI_Finalize();
        return -1;
    }
    if( (t == 0)&&(MaxX * MaxY != P))
    {
        if(rank == 0) printf("X*Y should be equal to the number of processes.\n");
        MPI_Finalize();
        return -1;
    }
    if( (t == 0)&& (N%P || M%P) )
    {
        if(rank == 0) printf("N and M should be divisible by number of processes.\n");
        MPI_Finalize();
        return -1;
    }
    return 0;
}

int CountElements(Matrix &m)
{
    Matrix::iterator it;
    int res = 0;
    for(it = m.begin(); it != m.end(); it++)
        res+=it->second.size();
    return res;
}

int GenerateMatrix(Matrix &MyBlock, MPI_Comm &Cartesian, int P, int MaxX, int MaxY, int MinWeight, int MaxWeight,Ind N, Ind M, int ITER_H)
/*
Every process generates a strip of N/P rows with uniformly distributed non-zeros in each row.
Each row contains between MinWeight and MaxWeight nonzeroes.
Than the strips are broken in chunks, M/MaxY columns in each chunk, and sent to corresponding processes of Cartesian
communicator, wich size is MaxX*MaxY.
Finally each process of cartesian communicator has a contiguous matrix block with (N/MaxX) columns and (M/MaxY) rows.
MyBlock and Cartesian are output parameters.
Returns 0 on success.
*/
{
    int rank;
    Matrix Strip, Columns;
    int  H = N/P;
    int Weight = (MaxWeight+MinWeight)/2;
    double StartTime = MPI_Wtime();
    double EndTime;
    double GenerationTime, TranspositionTime, DistributionTime;
    GenerationTime = TranspositionTime = DistributionTime = 0;
    int Dimensions[2] = {MaxX, MaxY};
    int Periods[2] = {1, 1}; 
    int IterationsNum = 1;
    double Ratio = double(MaxX)/double(MaxY);
    
    if(MPI_Cart_create(MPI_COMM_WORLD, 2, Dimensions, Periods, 0, &Cartesian) != MPI_SUCCESS)
        if(rank == 0) printf("Failed to create Cartesian communicator\n"); 
  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Ind BigN = N;
    if(H%ITER_H)
        IterationsNum = ceil((H)/float(ITER_H));
    else
        IterationsNum = H/ITER_H;
    if (rank == 0)    
         printf ("[Generation] Started %d iterations: ", IterationsNum);
    for(int i = 0; i < IterationsNum; i ++ )
    {
        int h;
        
        h = ITER_H;
        if(i == IterationsNum - 1)
            if(H%ITER_H)
                h = H%ITER_H;
        Strip.clear();
        GenerateStripRowwise( h, M, Weight, Weight, Strip, i*ITER_H*P + rank*h, 0);
        
        MPI_Barrier(MPI_COMM_WORLD);
        EndTime = MPI_Wtime();
        GenerationTime += (EndTime - StartTime);
        StartTime = EndTime;

        RowwiseToColumnwise(Strip, Columns);
   
        MPI_Barrier(MPI_COMM_WORLD);
        EndTime = MPI_Wtime();
        TranspositionTime += (EndTime - StartTime);
        StartTime = EndTime;
   
        DistributeMatrixChunks( P, MaxX, MaxY, Weight, N, M, Columns, MyBlock, Cartesian);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        EndTime = MPI_Wtime();
        DistributionTime += (EndTime - StartTime);
        StartTime = EndTime;
        if(rank == 0)
            printf("%d ", i+1);
    }
    if(rank == 0)
        printf("\n");
    if(rank == 0)
        printf("[Generation]P = %d, N = %lu, W = %d, X/Y = %f. Initialization took %f seconds.\n", P, N, Weight, Ratio, GenerationTime);
     if(rank == 0)
        printf("[Generation]P = %d, N = %lu, W = %d, X/Y = %f.  Transposition took %f seconds.\n", P, N, Weight, Ratio, TranspositionTime);
     if(rank == 0)
        printf("[Generation]P = %d, N = %lu, W = %d, X/Y = %f.   Distribution took %f seconds.\n",  P, N, Weight, Ratio, DistributionTime);
    return 0; 
}

int StoreMatrixToFolder(char *DirName,char *FileNamePrefix, 
                        Matrix &Block, int Type, int P, int MaxX, int MaxY, 
                        int MinWeight, int MaxWeight,  Ind N, Ind M,  
                        MPI_Comm Cartesian)
{
    int rank, size; 
    int MyCoords[2];
    FILE *Out;
    struct stat st;
    char *FileName;
    float StartTime = 1;
    float EndTime = 2;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Cart_coords(Cartesian, rank, 2, MyCoords);
   
    if(rank == 0)
        if(stat(DirName, &st) == -1)
            if(mkdir(DirName, 0755))
            {
                perror("[WriteMatrixToFolder] Could not create directory.\n");
                return 2;
            }
    MPI_Barrier(Cartesian);
    StartTime = MPI_Wtime();
    size = strlen(DirName) + 1 + strlen(FileNamePrefix) + 22;
    FileName = new char[size];
    sprintf(FileName, "%s/%s_%d_%d", DirName, FileNamePrefix, MyCoords[0], MyCoords[1]);
    
    size = CountElements(Block) + 2*(N/MaxX)*MaxWeight+ 2;
    int *Buf = new  int[size];
    int size1 = size;
    size = SerializeChunk(Block.begin(), Block.end(), size, Type, Buf);
    Out = fopen(FileName, "w");
    fwrite(Buf, sizeof(int), size, Out);
    fclose(Out);
    if(rank == 0)
    {
        sprintf(FileName, "%s/%s_config", DirName, FileNamePrefix);    
        Out = fopen(FileName, "w");
        fprintf(Out, "%d %d %d %d %d %lu %lu\n", P, MaxX, MaxY, MinWeight, MaxWeight, N, M);
        fclose(Out);
    }
    MPI_Barrier(Cartesian);
    EndTime = MPI_Wtime();
    if(rank == 0)
        printf("[StoreMatrixToFolder] Storing took %f seconds.\n", EndTime - StartTime);
    delete [] Buf;
    delete [] FileName;
    return 0;
}

int LoadMatrixFromFolder(char *DirName,char *FileNamePrefix, 
                        Matrix &Block, int &Type, int &MaxP, int &MaxX, int &MaxY, 
                        int &MinWeight, int &MaxWeight,  Ind &N, Ind &M,  
                        MPI_Comm Cartesian)
/*Reads parameters from $DirName/config, creates Cartesian communicator, loads corresponding blocks en each process.
 AdditionalDimensions is needed, if you use matrix for multipliclication with three types of nodes: workers, vector distributors, result producers.
 If AdditionalDimensions are non-zero, P should be (X+AdditionalDimensions)*(Y+AdditionalDimensions)!!
 If you want just to test Load efficiency, or mesure matrix characteristics, AdditionalDimensions should be zero. 
 * */
{
    int rank, size; 
    int MyCoords[2];
    FILE *In;
    struct stat st;
    char *FileName;
    int *Buf;
    MPI_Group World, ReadersGroup;
    MPI_Comm ReadersComm;
    int P;
    size = strlen(DirName) + 1 + strlen(FileNamePrefix) + 22;
    FileName = new char[size];
  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    sprintf(FileName, "%s/%s_config", DirName, FileNamePrefix);    
    In = fopen(FileName, "r");
    if(In == NULL)
    {
        printf("[LoadMatrixFromFolder] My rank is %d, can not open %s!\n", rank, FileName);
        return 1;
    }
    if(!fscanf(In, "%d %d %d %d %d %lu %lu\n", &P, &MaxX, &MaxY, &MinWeight, &MaxWeight, &N, &M))
    {
        printf("[LoadMatrixFromFolder] My rank is %d, can not read from config file %s\n", rank, FileName);    
        return 1;
    }
    fclose(In);
    if((MaxP < P)&& (rank == 0))
    {
        printf ("[LoadMatrix] Need at least %d processes, only %d given.\n", P, MaxP);
        return 1;  
    }

    MPI_Cart_coords(Cartesian, rank, 2, MyCoords);
    if((MyCoords[0] < MaxX)&&(MyCoords[1] < MaxY)) 
    {
    
        size = strlen(DirName) + 1 + strlen(FileNamePrefix) + 22;
        FileName = new char[size];
    
        if(!sprintf(FileName, "%s/%s_%d_%d", DirName, FileNamePrefix, MyCoords[0], MyCoords[1]))
        {
            printf("[LoadMatrixFromFolder] My rank is %d, sprintf failed.\n", rank);    
            return 1;
        }
        In = fopen(FileName, "r");
        if(!fread(&size, sizeof(int), 1, In))
        {
         
            printf("[LoadMatrixFromFolder] My rank is %d, fread size failed.\n", rank);
            return 1;
        } 
        struct stat info;
        fstat(fileno(In), &info);
        Buf = (int*) mmap(NULL, info.st_size , PROT_READ, MAP_PRIVATE, fileno(In), 0);
        DeserializeChunk(Buf, Type, Block);
        munmap(Buf, size*sizeof(int));
        fclose(In);
        delete [] FileName;
    }
    return 0;   
}
int LoadIterationForThreeProcessTypes(char *DirName,char *FileNamePrefix, 
                                      Matrix &Block, int &ProcessType, int &MaxP, int &MaxX, int &MaxY, 
                                      int &MinWeight, int &MaxWeight,  Ind &N, Ind &M,  
                                      MPI_Comm &Cartesian)
/*
 *Create Cartesian Communicator (MaxX+1)*(MaxY+1). Load matrix from DirName/FileNamePrefix into processes that are in top left corner of this communicator and assign to this processes ProcessType = 0.
Processes with cart coords [0:MaxX+1, MaxY] are ResultColletors. They are assigned blocks of size [1, N/MaxX] of zeroes and ProcessType 1;
Processes with cart coords [MaxX,  0:MaxY+1] are RightDistributors. They are assigned blocks of size [1, M/MaxY] of ones and ProcessType 2.
The last one process [MaxX+1, MaxY+1] has ProcessType 3.
 * */
{
    int Dimensions[2] = {MaxX+1, MaxY+1};
    int MatrixMaxX, MatrixMaxY;
    int Periods[2] = {1, 1}; 
    int MyCoords[2] = {-1, -1};
    int ChunkType;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(MaxP < (MaxX+1)*(MaxY+1))
    {
        if(rank == 0)
            printf("[LoadIteration] Need at least %d processes, %d given.\n", (MaxX+1)*(MaxY+1), MaxP);
        return 1;
    }
    if(MPI_Cart_create(MPI_COMM_WORLD, 2, Dimensions, Periods, 0, &Cartesian) != MPI_SUCCESS)
        if(rank == 0) printf("Failed to create Cartesian communicator\n"); 
    
    MPI_Cart_coords(Cartesian,rank, 2,  MyCoords);
    LoadMatrixFromFolder(DirName, FileNamePrefix, 
                         Block, ChunkType, MaxP,  MatrixMaxX, MatrixMaxY, 
                         MinWeight,  MaxWeight, N, M,  
                         Cartesian);
    
    ProcessType = -1;
    if(MyCoords[1] < MatrixMaxY)
    {
        ProcessType = 0;
    }
    if(MyCoords[0] == MatrixMaxX && MyCoords[1] == MatrixMaxY)
    {
        ProcessType = 3;
        return 0;
    }
    
    if(MyCoords[1] == MatrixMaxY)
    {
        ProcessType = 1;
        Block.clear();
    }
    
    if(MyCoords[0] == MatrixMaxX)
    {
        ProcessType = 2;
        Block.clear();
        for(unsigned int i = 0; i < M/MatrixMaxY; i++)
           Block[0].insert(MyCoords[1]*M/MatrixMaxY + i);
    }
    MaxX = MatrixMaxX;
    MaxY = MatrixMaxY;
    return 0;
     
}

int MulBroadcast3(int IterNum, Matrix MyBlock, int Type,  MPI_Comm Cartesian, int P, int MaxX, int MaxY, Ind N, Ind M, int Weight)
/* Three types of processes: Workers, Right Distributors, ResHandlers. 
 * For this function, Matrix chunks on workers should be stored rowwise!*/
{
   
    int MyCoords[2] = {-1, -1};
    int rank = -1;
    int RightSendSize, ResSendSize;
    int *ResSendBuf, *ResRecieveBuf, *RightBuf;
    int ChunkType, ierr;
    Matrix  Right; //In this version only one line is stored in this matrix.
    MPI_Comm HorizontalComm; //For sending result.
    MPI_Comm VerticalComm; //For retrieveing right hand vector.

    RightSendSize = (M/MaxY)*2 + 20; //Upper estimate of maximum send size of right hand part.
    ResSendSize = N/MaxX + 20; //20 + Number of partial result elements, produced by each Worker.
    ResSendBuf = new int[ResSendSize];
    ResRecieveBuf = new int[ResSendSize];//MPI_Reduce needs two non aliased buffers. 
    RightBuf = new int[RightSendSize];
    if((ResRecieveBuf == NULL)||(ResSendBuf == NULL)||(RightBuf == NULL))
    {
        printf("[MulBroadcast3] My rank is %d, I've failed to allocate %d elements.\n", 
                rank, 2*ResSendSize + RightSendSize);
        return 1;
    }

    

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Cart_coords(Cartesian,rank, 2,  MyCoords);
    MPI_Comm_split(Cartesian, MyCoords[0], MyCoords[1], &HorizontalComm);// MaxX communicators, processes [0..MaxY] in each
    MPI_Comm_split(Cartesian, MyCoords[1], MyCoords[0], &VerticalComm);// MaxY communicators, processes [0..MaxX] in each
    
    memset(ResRecieveBuf,0, ResSendSize*sizeof(int));
    memset(ResSendBuf,0, ResSendSize*sizeof(int));
    MPI_Barrier(Cartesian);
    float Ratio = ((float)MaxX)/MaxY;
    double StartTime, EndTime;
    
    for(int Iter = 0; Iter < IterNum; Iter++)
    {
        /*Getting right hand vector, part of job for Workers and RightDistributors*/
        StartTime = MPI_Wtime();
        if((Type == 0)||(Type == 2))
        {

            if(MyCoords[0] == MaxX)
                SerializeChunk(MyBlock.begin(), MyBlock.end(), RightSendSize, 0, RightBuf );
       
             MPI_Bcast(RightBuf, RightSendSize, MPI_INT, MaxX, VerticalComm);
      
            if(MyCoords[0] != MaxX)
                DeserializeChunk(RightBuf, ChunkType, Right);     
        }
        EndTime = MPI_Wtime();
        if(rank == 0)
            printf("[MulBroadcast3] P = %d, N = %lu, W = %d, X/Y = %f. Right part distribution took %f\n", P, N, Weight, Ratio,EndTime - StartTime);
        StartTime = EndTime;
         /*Calculate partial result and put its image into ResSendBuf.*/
        if(Type == 0)
            for(Matrix::iterator it = MyBlock.begin(); it != MyBlock.end(); it++)
                {
                    int ind, r;
                    r =  MulLine(it->second, Right[0])%2; 
                    ind = it->first - MyCoords[0]*N/MaxX;
                    ResSendBuf[ind] = r;
                }
        MPI_Barrier(HorizontalComm);
        EndTime = MPI_Wtime();
        if(rank == 0)
            printf("[MulBroadcast3] P = %d, N = %lu, W = %d, X/Y = %f. Partial calculations took    %f\n",  P, N, Weight, Ratio,EndTime - StartTime);
        StartTime = EndTime;
   
        /*Now Workers send partial results to ResultHandlers*/
        if((Type == 0)||(Type == 1))
            MPI_Reduce(ResSendBuf, ResRecieveBuf, ResSendSize, MPI_INT, MPI_LXOR, MaxY, HorizontalComm);
          
        /*ResultHandlers put results from ResRecieveBuf in MyBlock[0] with needed offset.*/
        if(Type == 1)
        {
            MyBlock[0].clear();
            for(unsigned int i = 0; i < N/MaxX; i++)
            {
                if(ResRecieveBuf[i])
                    MyBlock[0].insert(i+MyCoords[0]*N/MaxX);
            }
         }
        EndTime = MPI_Wtime();
        if(rank == 0)
            printf("[MulBroadcast3] P = %d, N = %lu, W = %d, X/Y = %f. Final calculations took      %f\n", P, N, Weight, Ratio,EndTime - StartTime );
        StartTime = EndTime;
 
        /*ResultHandlers send the result to RightDistributors for the next iteration*/
           if((Type == 1)||(Type == 2))
            RedistributeResVector3(MyBlock, Type, MyCoords[0], MyCoords[1], N, MaxX, MaxY, ResSendBuf,ResSendSize, Cartesian);
    
        MPI_Barrier(VerticalComm);
        EndTime = MPI_Wtime();
        if(rank == 0)
            printf("[MulBroadcast3] P = %d, N = %lu, W = %d, X/Y = %f. Result redistribution took   %f\n",  P, N, Weight, Ratio,EndTime - StartTime);
        StartTime = EndTime;
 
    }
    delete [] ResRecieveBuf;
    delete [] ResSendBuf;
    delete [] RightBuf;
    return 0;
}

