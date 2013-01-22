#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <string.h>
#include "umvm_internals.h"
#include "umvm_debug.h"

void TryGenerateStripRowwise(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight)
{
    Matrix strip;
    int rank;
    int H;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    H = N/P;
    GenerateStripRowwise(rank*H, (rank+1)*H, 0, M, Weight, Weight, strip);
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

    for(Matrix::iterator it = Columns.begin(); it != Columns.end(); it++)
        for(std::set<Ind>::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++)
              printf ("My rank is %d, I have [%lu,%lu]\n", rank, it->first, *it1);
}



void TrySerializeChunk(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight)
{
    Matrix strip;
    int rank;
    int H;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    H = N/P;
    GenerateStripRowwise(rank*H, (rank+1)*H, 0, M, Weight, Weight, strip);
    int  *Buf = new int [Weight*H+2*H+2];
    if(!SerializeChunk(strip.begin(), strip.end(), Weight*H+2*H+2, 1, Buf))    
    {
        delete [] Buf;
        printf("Buffer overflow in SeriaizeChunk. SerializeChunk returned 0.\n");
        return;
    }
    if(rank == 0)
        for(unsigned int i = 0; i < Weight*H +2*H +2; i++)
            printf("%d ", Buf[i]);
    printf("\n");
    delete [] Buf;
}


void TryDeserializeChunk(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight)
{
    Matrix strip, Chunk;
    int rank;
    int H;
    int Type;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    H = N/P;
    GenerateStripRowwise(rank*H, (rank+1)*H, 0, M, Weight, Weight, strip);
    PrintMatrixStructure(strip);
    int  *Buf = new int [Weight*H+2*H+2];
    printf("Serializing...\n");
    if(!SerializeChunk(strip.begin(), strip.end(), Weight*H+2*H+2, 1, Buf))    
    {
        delete [] Buf;
        printf("Buffer overflow in SeriaizeChunk. SerializeChunk returned 0.\n");
        return;
    }
    DeserializeChunk(Buf, Type, Chunk);
    PrintMatrixStructure(Chunk);
    delete [] Buf;
}


int TryGetAndInsertChunk(int P, int MaxX, int MaxY, int MaxWeight, Ind  N, Ind M, MPI_Comm Cartesian)
/*Modelling chunk distribution in 1 process, debug purposes only.*/
{
   
   int size;
    int StripH, BlockH, BlockW;
    Matrix::iterator ChunkStart, ChunkEnd;
    int MaxSendSize = 0;
    int SendSize = 0;
    int DestCoords[2];
    int MyCoords[2];
    int Type;
    Matrix Chunk, Columns, Strip;
    MPI_Status status;
    Matrix Blocks[P];
  
    StripH = N/P;
    BlockH = N/MaxX;
    BlockW = M/MaxY;
    MaxSendSize = BlockH*BlockW + 2*BlockH*MaxWeight+ 2+1;
  
    int *SendBuf = new (std::nothrow) int[MaxSendSize];
    int *RecieveBuf = new (std::nothrow) int[MaxSendSize];
    
    memset(RecieveBuf, 0, MaxSendSize*sizeof(int));
    int ElementsBefore = 0;
    int ElementsAfter = 0;
    for(int rank = 0; rank < P; rank++)
    {
    Columns.clear();
    GenerateStripRowwise(rank*StripH, (rank+1)*StripH, 0, M, MaxWeight, MaxWeight, Strip);
    RowwiseToColumnwise(Strip, Columns);
    ElementsBefore += CountElements(Columns);
    printf("Rank %d generated %d elements\n", rank, CountElements(Columns));
    {
            for(int j = 0; j < P; j++)
                if(rank != j)
                {
                    size = 1;
                    if(GetChunkDestinatedToXY(j/MaxX, j%MaxY, rank,  P, MaxX, MaxY, N, M, 
                                              Columns, ChunkStart, ChunkEnd))
                    {
                       
                       Matrix mm(ChunkStart, ChunkEnd);
                       printf("Rank %d sent %d elements to %d\n", rank, CountElements(mm), j);
                       if(! (size = SerializeChunk(ChunkStart, ChunkEnd, MaxSendSize, 0, SendBuf)))
                           return 2;
                       printf("Sent size was %d, max is %d\n", size, MaxSendSize);
                    }
                    else                            
                    {
                        SendBuf[0]= 0;
                    }
                    memcpy( RecieveBuf, SendBuf, sizeof(int)*size); 
                    if(RecieveBuf[0])
                    {
                        if(!DeserializeChunk(RecieveBuf, Type, Chunk))
                            return 3;
                        
                        AddChunkToBlock(Blocks[j], Chunk.begin(), Chunk.end());
                        printf("Rank %d recieved %d elements from %d, and have now %d\n", 
                                j, CountElements(Chunk), rank, CountElements(Blocks[j]));
                       printf("RecievedBuf ");
                       for (int l = 0; l < size; l++)
                           printf("%d ",RecieveBuf[l]);
                       printf("\n");
                       Chunk.clear();
                    }
           
                
                }
        
            
        if(GetChunkDestinatedToXY(rank/MaxX, rank%MaxY, rank, P, MaxX, MaxY, N, M, Columns, ChunkStart, ChunkEnd))
            AddChunkToBlock(Blocks[rank], ChunkStart, ChunkEnd);
        Matrix mm(ChunkStart, ChunkEnd);
        printf("Rank %d recieved %d elements from %d, and have now %d\n", 
                                rank, CountElements(mm), rank, CountElements(Blocks[rank]));
          
        }
    }
    for(int k = 0; k < P; k++)
        ElementsAfter += CountElements(Blocks[k]);
    printf("[TryGetAndInsert] Before %d, After %d\n", ElementsBefore, ElementsAfter);
    delete [] RecieveBuf;
    delete [] SendBuf;
}


