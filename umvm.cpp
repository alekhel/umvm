#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <string.h>
#include "umvm.h"

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
};

void TryGenerateStripRowwise(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight)
{
    Matrix strip;
    int rank;
    int H;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    H = N/P;
    GenerateStripRowwise(rank*H, (rank+1)*H, 0, M, Weight, Weight, strip);
    for (int i = rank*H; i < (rank+1)*H; i++)
        for(unsigned int j = 0;  j < strip[i].size(); j ++ )
            printf ("My rank is %d, my  strip[%d][%d] = %7lu\n", rank, i, j,  strip[i][j]);
}

void RowwiseToColumnwise(Matrix Rows, Matrix &Columns)
{
    Columns.clear();
    for(Matrix::iterator it = Rows.begin(); it != Rows.end(); it++)
        for(unsigned int i = 0; i < it->second.size(); i++)
            Columns[it->second[i]].push_back(it->first);
}

void TryRowwiseToColumnwise(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight)
{
    Matrix Columns, Rows;
    int rank;
    int H;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    H = N/P;
    GenerateStripRowwise(rank*H, (rank+1)*H, 0, M, Weight, Weight, Rows);
    RowwiseToColumnwise(Rows, Columns);

    #ifdef DEBUG_PRINT
    for(Matrix::iterator it = Columns.begin(); it != Columns.end(); it++)
        for(unsigned int j = 0;  j < it->second.size(); j ++ )
              printf ("My rank is %d, my  Columns[%lu][%u] = %7lu\n", rank, it->first, j,  it->second[j]);
    #endif
}


int SerializeChunk(Matrix::iterator Start, Matrix::iterator End, unsigned int Size, int Type,  int Res[])
/*Res should be already allocated! 
 Size defines max buffer size. 
 Type is 0 for Columnwise blocks and 1 for Rowwise blocks
 Returns used space, 0 if there is not enough space in buffer and buffer checks are enabled.
 Res[0] is first dimension size.
 Res[1] is Type.
 Then there is array of first dimension indexes, size Res[0].
 Then there is array of Sum_{i = 0, i = j} Length[i].
 Then there are arrays of second dimension data.
*/
{
    memset(Res, 0, Size * sizeof(int));
    int ColumnsCount = 0;
    int Offset = 0;
    int HeaderSize = 2;
    unsigned int UsedSpace = HeaderSize;
 
    for(Matrix::iterator it = Start; it != End; it++)
    {
        #ifdef BUFFER_CHECKS
        if(UsedSpace == Size)
        {
            printf("Buffer overflow in SerializeChunk while writing first dimension.\n");
            return 0;
        }
        #endif
        Res[UsedSpace++] = it->first; 
        ColumnsCount++;
    }
    
    Res[0] = ColumnsCount;
    Res[1] = Type;
    #ifdef BUFFER_CHECKS
    if(UsedSpace+ColumnsCount >= Size)
    {
        printf("Buffer overflow in SerializeChunk while writing offsets.\n");
        return 0;
    }
    #endif
    for(Matrix::iterator it = Start; it != End; it++)
    {
        
        #ifdef BUFFER_CHECKS
        if( HeaderSize + 2*ColumnsCount + Offset + it->second.size() > Size)
        {
            printf("Buffer overflow in SerializeChunk while writing second dimension. Needed %lu, given %u\n",
                    HeaderSize + 2*ColumnsCount + Offset + it->second.size(), Size );
            return 0;
        }
        #endif
 
        for(unsigned int i = 0; i < it->second.size(); i++)
        {
            Res[HeaderSize + 2*ColumnsCount + Offset + i] = it->second[i];    
        }
        Offset += it->second.size();
        Res[UsedSpace++] = Offset;
    }
    return UsedSpace+Offset;
}


void TrySerializeChunk(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight)
{
    Matrix strip;
    int rank;
    int H;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    H = N/P;
    GenerateStripRowwise(rank*H, (rank+1)*H, 0, M, Weight, Weight, strip);
    int  Buf[Weight*H+2*H+2];
    if(!SerializeChunk(strip.begin(), strip.end(), Weight*H+2*H+2, 1, Buf))    
    {
        printf("Buffer overflow in SeriaizeChunk. SerializeChunk returned 0.\n");
        return;
    }
    if(rank == 0)
        for(unsigned int i = 0; i < Weight*H +2*H +2; i++)
            printf("%d ", Buf[i]);
    printf("\n");
}

void DistributeMatrixChunks(int CartX, int CartY, int P, int MaxX, int MaxY, 
                            Ind  N, Ind M, Matrix &Columns, MPI_Comm Cartesian)
{
    int rank;
    int StripH, BlockH, BlockW;
    int SendCount[MaxX][MaxY];
    int SendSize[MaxX][MaxY];
    int RecieveCount[MaxX][MaxY];
    int RecieveSize[MaxX][MaxY];
    int MyMaxSendSize = 0;
    int CurSendSize = 0;
    int HeaderSize = 2;
    memset(SendCount, 0, MaxX*MaxY*sizeof(int));
    memset(SendSize, HeaderSize*sizeof(int), MaxX*MaxY*sizeof(int)); 
    memset(RecieveCount, 0, MaxX*MaxY*sizeof(int));
    memset(RecieveSize, 0, MaxX*MaxY*sizeof(int));
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    StripH = N/P;
    BlockH = N/MaxX;
    BlockW = M/MaxY;


    for(Matrix::iterator it = Columns.begin(); it != Columns.end(); it++)
    {
        SendCount[StripH*rank/BlockH][it->first/BlockW] |= 1;
        CurSendSize =  (2 + it->second.size())*sizeof(Ind); 
        SendSize[StripH*rank/BlockH][it->first/BlockW] += CurSendSize;
        if(MyMaxSendSize < CurSendSize) MyMaxSendSize = CurSendSize;
    }
    MPI_Allreduce(SendCount, RecieveCount, MaxX*MaxY, MPI_INT, MPI_SUM, Cartesian);    
    #ifdef DEBUG_PRINT
    printf("My coords are [%d, %d], I'll  recieve %d messages.\n", CartX, CartY, RecieveCount[CartX][CartY]);
    for (int i = 0; i < MaxX; i++)
        for(int j = 0; j < MaxY; j++)
            printf("My coords are [%d, %d], I'll send %d messages to [%d, %d].\n",
                    CartX, CartY, SendCount[i][j], i, j);
    #endif
    MPI_Allreduce(SendSize, RecieveSize, MaxX*MaxY, MPI_INT, MPI_MAX, Cartesian);    
    

}


