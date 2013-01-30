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
int ParameterSanityCheck(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight)
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
    if( MaxX * MaxY != P)
    {
        if(rank == 0) printf("X*Y should be equal to the number of processes.\n");
        MPI_Finalize();
        return -1;
    }
    if( N%P || M%P )
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

int GenerateMatrix(Matrix &MyBlock, MPI_Comm &Cartesian, int P, int MaxX, int MaxY, int MinWeight, int MaxWeight,Ind N, Ind M)
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
    int Dimensions[2] = {MaxX, MaxY};
    int Periods[2] = {1, 1}; 
    double Ratio = double(MaxX)/double(MaxY);
    
    if(MPI_Cart_create(MPI_COMM_WORLD, 2, Dimensions, Periods, 0, &Cartesian) != MPI_SUCCESS)
        if(rank == 0) printf("Failed to create Cartesian communicator\n"); 
  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
    GenerateStripRowwise(rank*H, (rank+1)*H, 0, M, Weight, Weight, Strip);
    MPI_Barrier(MPI_COMM_WORLD);
    
    EndTime = MPI_Wtime();
    if(rank == 0)
        printf("[Generation]P = %d, N = %lu, W = %d, X/Y = %f. Initialization took %f seconds.\n", P, N, Weight, Ratio, EndTime - StartTime);
    StartTime = EndTime;

    RowwiseToColumnwise(Strip, Columns);
   
    MPI_Barrier(MPI_COMM_WORLD);
    EndTime = MPI_Wtime();
    if(rank == 0)
        printf("[Generation]P = %d, N = %lu, W = %d, X/Y = %f.  Transposition took %f seconds.\n", P, N, Weight, Ratio, EndTime - StartTime);
    StartTime = EndTime;
   
    DistributeMatrixChunks( P, MaxX, MaxY, Weight, N, M, Columns, MyBlock, Cartesian);
    MPI_Barrier(MPI_COMM_WORLD);
    EndTime = MPI_Wtime();

    if(rank == 0)
        printf("[Generation]P = %d, N = %lu, W = %d, X/Y = %f.   Distribution took %f seconds.\n",  P, N, Weight, Ratio, EndTime - StartTime);
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
    
    size = strlen(DirName) + 1 + strlen(FileNamePrefix) + 22;
    FileName = new char[size];
    sprintf(FileName, "%s/%s_%d_%d", DirName, FileNamePrefix, MyCoords[0], MyCoords[1]);
    
    size = CountElements(Block) + 2*(N/MaxX)*MaxWeight+ 2;
    int *Buf = new  int[size];
    size = SerializeChunk(Block.begin(), Block.end(), size, Type, Buf);
    printf("Serialized size %d, elements number %d\n",size, CountElements(Block) );
    Out = fopen(FileName, "w");
    fwrite(Buf, sizeof(int), size, Out);
    fclose(Out);
    if(rank == 0)
    {
        sprintf(FileName, "%s/config", DirName);    
        Out = fopen(FileName, "w");
        fprintf(Out, "%d %d %d %d %d %lu %lu\n", P, MaxX, MaxY, MinWeight, MaxWeight, N, M);
        fclose(Out);
    }
    delete [] Buf;
    delete [] FileName;
    return 0;
}

int LoadMatrixFromFolder(char *DirName,char *FileNamePrefix, 
                        Matrix &Block, int &Type, int &P, int &MaxX, int &MaxY, 
                        int &MinWeight, int &MaxWeight,  Ind &N, Ind &M,  
                        MPI_Comm &Cartesian)
/*Reads parameters from $DirName/config, creates Cartesian communicator, loads corresponding blocks en each process.*/
{
    int rank, size; 
    int MyCoords[2];
    FILE *In;
    struct stat st;
    char *FileName;
    int *Buf;
    size = strlen(DirName) + 1 + strlen(FileNamePrefix) + 22;
    FileName = new char[size];
  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    sprintf(FileName, "%s/config", DirName);    
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
    
    int Dimensions[2] = {MaxX, MaxY};
    int Periods[2] = {1, 1}; 
    if(MPI_Cart_create(MPI_COMM_WORLD, 2, Dimensions, Periods, 0, &Cartesian) != MPI_SUCCESS)
        if(rank == 0) printf("Failed to create Cartesian communicator\n"); 
    MPI_Cart_coords(Cartesian, rank, 2, MyCoords);
 
    size = strlen(DirName) + 1 + strlen(FileNamePrefix) + 22;
    FileName = new char[size];
    
    if(!sprintf(FileName, "%s/%s_%d_%d", DirName, FileNamePrefix, MyCoords[0], MyCoords[1]))
    {
        printf("[LoadMatrixFromFolder] My rank is %d, sprintf failed.\n", rank);    
        return 1;
    }
 
    
    
    In = fopen(FileName, "r");
    if(!fread(&size, sizeof(int), 1, In))
     
    
    Buf = (int*) mmap(NULL, size*sizeof(int), PROT_READ, MAP_PRIVATE, fileno(In), 0);
    DeserializeChunk(Buf, Type, Block);
    munmap(Buf, size*sizeof(int));
    fclose(In);
    delete [] FileName;
    return 0;   
}

