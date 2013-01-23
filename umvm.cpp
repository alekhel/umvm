#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <string.h>
#include "umvm.h"
#include "umvm_internals.h"
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

int GenerateMatrix(Matrix &MyBlock, MPI_Comm &Cartesian, int P, Ind MaxX, Ind MaxY, int MinWeight, int MaxWeight,Ind N, Ind M)
/*
Every process generates a strip of N/P rows with uniformly distributed non-zeros in each row.
Each row contains between MinWeight and MaxWeight nonzeroes.
Than the strips are broken in chunks, M/MaxY columns in each chunk, and sent to corresponding processes of Cartesian
communicator, wich size is MaxX*MaxY.
Finally each process of cartesian communicator has a contiguous matrix block with (N/MaxX) columns and (M/MaxY) rows.
MyBlock and Cartesian are output parameters.
*/
{
    int rank;
    Matrix Strip, Columns, Block;
    int  H = N/P;
    int Weight = (MaxWeight+MinWeight)/2;
    double StartTime = MPI_Wtime();
    double EndTime;
    int Dimensions[2] = {MaxX, MaxY};
    int Periods[2] = {1, 1}; 
 
    if(MPI_Cart_create(MPI_COMM_WORLD, 2, Dimensions, Periods, 0, &Cartesian) != MPI_SUCCESS)
        if(rank == 0) printf("Failed to create Cartesian communicator\n"); 
  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
    GenerateStripRowwise(rank*H, (rank+1)*H, 0, M, Weight, Weight, Strip);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    EndTime = MPI_Wtime();
    if(rank == 0)
        printf("[main] Generation of strip took %f seconds.\n", EndTime - StartTime);
    StartTime = EndTime;

    RowwiseToColumnwise(Strip, Columns);
   
    MPI_Barrier(MPI_COMM_WORLD);
    EndTime = MPI_Wtime();
    if(rank == 0)
        printf("[main] Transposition took %f seconds.\n",   EndTime - StartTime);
    StartTime = EndTime;
   
    DistributeMatrixChunks( P, MaxX, MaxY, Weight, N, M, Columns, Block, Cartesian);
    MPI_Barrier(MPI_COMM_WORLD);
    EndTime = MPI_Wtime();
    if(rank == 0)
        printf("[main] Distribution took %f seconds.\n",  EndTime - StartTime);
  
}
