#include "umvm.h"
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <string.h>
#include "math.h"
#include <stdio.h>
#include <time.h>
#ifdef C11RANDOM
#include <random>
#else
#include <cstdlib>
#endif



#ifdef C11RANDOM
int GenerateStripRowwise( int Height, int Width, int MinWeight, int MaxWeight, 
                          Matrix &m, Ind RowOffset = 0, Ind ColumnOffset = 0)
/*Attention, columns indexes are not sorted by default*/
/*The first version does not use MaxWeight parameter*/
{
    
    m.clear();
    std::default_random_engine generator(4);
    std::uniform_int_distribution<unsigned int> distribution(0, Width);
    for (Ind i = 0; i < Height; i++)
    {
        for(int j = 0; j < MinWeight; j++)
        {
            m[i + RowOffset].insert(ColumnOffset + distribution(generator));
        }     
    }
    return 0;
}
#else
int GenerateStripRowwise( int Height, int Width, int MinWeight, int MaxWeight, 
                          Matrix &m, Ind RowOffset = 0, Ind ColumnOffset = 0)
/*Attention, columns indexes are not sorted by default*/
/*The first version does not use MaxWeight parameter*/
{
    
    m.clear();
    srand(time(NULL)^RowOffset); 
    for (int i = 0; i < Height; i++)
    {
        for(int j = 0; j < MinWeight; j++)
        {
            m[i + RowOffset].insert( ColumnOffset +  floor(Width*(double(rand())/RAND_MAX)));
        }     
    }
    return 0;
}
#endif

void RowwiseToColumnwise(Matrix Rows, Matrix &Columns)
{
    Columns.clear();
    for(Matrix::iterator it = Rows.begin(); it != Rows.end(); it++)
        {
        for(Line::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++)
            Columns[*it1].insert(it->first);
        }
}
int GetMaxBufferSize(Matrix::iterator Start, Matrix::iterator End)
/*Estimates maximum buffer size (Inds number) needed to serialize a chunk*/
{
Matrix m(Start, End);
return  CountElements(m) + 2*m.size() + 20;    
}

int GetMaxBufferSize(Matrix &Block)
/*Estimates maximum buffer size (Inds number) needed to serialize a chunk*/
{
return  CountElements(Block) + 2*Block.size() + 20;    
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
 Variables are named as if first dimension was Columns, second dimension was Rows.
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
 
        int i = 0;
        for(Line::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++, i++)
        {
            Res[HeaderSize + 2*ColumnsCount + Offset + i] = *it1;    
        }
        Offset += it->second.size();
        Res[UsedSpace++] = Offset;
    }
    return UsedSpace+Offset;
}

int DeserializeChunk(int Buf[], int &Type, Matrix &Chunk)
/*Res should be already allocated! 
 Size defines max buffer size. 
 Type is 0 for Columnwise blocks and 1 for Rowwise blocks
 Returns used space, 0 if there is not enough space in buffer and buffer checks are enabled.
 Res[0] is first dimension size.
 Res[1] is Type.
 Then there is array of first dimension indexes, size Res[0].
 Then there is array of Sum_{i = 0, i = j} Length[i].
 Then there are arrays of second dimension data.
 Variables are named as if first dimension was Columns, second dimension was Rows.
Returns 0 on failure, Size otherwise.
*/
{
    unsigned int ColumnsCount = 0;
    unsigned int HeaderSize = 2;
    unsigned int Size;
    Chunk.clear();
   
    ColumnsCount = Buf[0];
    if(!ColumnsCount)
        return 0;
    Type = Buf[1];
    Size = Buf[HeaderSize + 2*ColumnsCount - 1]  + HeaderSize + 2*ColumnsCount;
    #ifdef BUFFER_CHECKS
    if(HeaderSize+ColumnsCount*2 >= Size)
    {
        printf("[Deserialize] Columns %d, LastOffset Index %d, Last Offset %d, Size %d\n",ColumnsCount, HeaderSize + 2*ColumnsCount - 1,  Buf[HeaderSize + 2*ColumnsCount - 1], Size );
        return 0;
    }
    #endif
    int j = 0;
    for(unsigned int i = 0; i < ColumnsCount; i++)
    {
        
        while(j++ < Buf[i + HeaderSize + ColumnsCount])
            Chunk[Buf[i+ HeaderSize]].insert(Buf[j + HeaderSize + 2*ColumnsCount - 1]);
        j--;
        #ifdef BUFFER_CHECKS
        if((unsigned int)j >= Size)
        {
            printf("Buffer overflow in DeserializeChunk while reading elements, %d elements read.\n", 
                    j - HeaderSize - 2*ColumnsCount);
            return 0;
        }
        #endif
 
    }
    return Size;
}


void PrintMatrixStructure(Matrix &m)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for(Matrix::iterator it = m.begin(); it != m.end(); it++)
    {
        printf("[PrintMatrix] My rank is %d, I have %lu elements from  line %lu\n", rank, it->second.size(), it->first);
    }
}
void PrintMatrix(Matrix &m)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for(Matrix::iterator it = m.begin(); it != m.end(); it++)
    {
        printf("[PrintMatrix] My rank is %d, I have %lu elements from  line %lu:", rank, it->second.size(), it->first);
        for(Line::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++ )
            printf(" %lu ", *it1);
         printf("\n");
    }
}

void PrintSerialized(Matrix &m, char* prefix)
{
    int s = 3*CountElements(m) + 2;
    int Buf[s];
    s = SerializeChunk(m.begin(), m.end(), s, 0, Buf);
    printf("%s", prefix);
    for(int i = 0; i < s; i++)
        printf("%d ", Buf[i]);
    printf("\n");
}

int GetChunkDestinatedToXY(int X, int Y, int rank, int P, int MaxX, int MaxY, 
                            Ind N, Ind M, Matrix &Columns, 
                            Matrix::iterator &ChunkStart, Matrix::iterator &ChunkEnd)
/*
Chunk is an output parameter.
This version distributes equal number of columns, but only across processes that should recieve columns from this one.
iReturns 0 if there is no elements destinated to XY, 1 otherwise.
Example for P = 4, MaxX = MaxY = 2 (Numbers are ranks of origin generator processes)
1 1 1 1        1 1 2 2
2 2 2 2  ---\  1 1 2 2
3 3 3 3  ---/  3 3 4 4
4 4 4 4        3 3 4 4
*/
{
    int StripH, BlockH, BlockW;
    Ind Lower, Upper;
    StripH = N/P;
    BlockH = N/MaxX;
    BlockW = M/MaxY;
    Matrix::iterator Start, End;
    if(X != StripH*rank/BlockH)
        return 0;
 
    if(Columns.end() != (Start = Columns.lower_bound(Y*BlockW)))
      {
          End = Columns.lower_bound((Y+1)*BlockW);
          Upper = (--End)->first;
          Lower = Start->first;
          if(( Upper >= (unsigned int)(Y*BlockW)) && (Lower <= Upper)) 
          {
                    ChunkStart = Start;
                    ChunkEnd = ++End; 
                    return 1;
          }
      }
    return 0;    
}

void AddChunkToBlock(Matrix &Block, Matrix::iterator Start, Matrix::iterator End)
{
    for(Matrix::iterator it = Start; it != End; it++)
       Block[it->first].insert(it->second.begin(), it->second.end());
}

int DistributeMatrixChunks(int P, int MaxX, int MaxY, int MaxWeight, Ind  N, Ind M, Matrix &Columns, Matrix &Block, MPI_Comm Cartesian)
/*
Each process sends a message to all others, containing 0 if it has no elements destinated to reciever, packed chunk.
otherwise. Packages with elements start with number of first dimension entries, could not be 0.
Elements to send are determined by GetChunkDestinatedToXY function.
The main problem of this function is memory cost.
Columns are erased as soon, as they are sent.
return values:
0 success
1 not enough memory
2 SerializeChunk failed
3 DeserializeChunk failed
*/
{
    int rank, size;
    int StripH, BlockH, BlockW;
    Matrix::iterator ChunkStart, ChunkEnd;
    int MaxSendSize = 0;
    int SendSize = 0;
    int DestCoords[2];
    int MyCoords[2];
    int Type;
    Matrix Chunk;
    MPI_Status status;

    StripH = N/P;
    BlockH = N/MaxX;
    BlockW = M/MaxY;
    MaxSendSize = CountElements(Columns) + 2*Columns.size() + 20;
    int *SendBuf = new (std::nothrow) int[MaxSendSize];
    int *RecieveBuf = new (std::nothrow) int[MaxSendSize];
    if(RecieveBuf == NULL)
    {
        printf("[DistributeMatrixChunks] My rank is %d, I've failed to allocate %d elements to RecieveBuffer\n", 
                rank, MaxSendSize);
        return 1;
    }
    memset(RecieveBuf, 0, MaxSendSize*sizeof(int));
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for(int i = 0; i < P; i++)
        if(rank == i)   
        {
            for(int j = 0; j < P; j++)
                if(i != j)
                {
                    MPI_Cart_coords(Cartesian, j, 2, DestCoords);
                    if(GetChunkDestinatedToXY(DestCoords[0], DestCoords[1], rank,  P, MaxX, MaxY, N, M, 
                                              Columns, ChunkStart, ChunkEnd))
                    {
                       SendSize = MaxSendSize;
                       if(SendBuf == NULL)
                           printf("Not enough memory, need %d.\n", SendSize);
                       if(! (size = SerializeChunk(ChunkStart, ChunkEnd, SendSize, 0, SendBuf)))
                           return 2;
                       Columns.erase(ChunkStart, ChunkEnd);
                       MPI_Send(SendBuf, size, MPI_INT, j, 0 , MPI_COMM_WORLD);
                    }
                    else                            
                    {
                        SendBuf[0]= 0;
                        MPI_Send(SendBuf, 1, MPI_INT, j, 0 , MPI_COMM_WORLD);
                    }
                }
        }
        else
        {
            MPI_Recv(RecieveBuf, MaxSendSize, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            if(RecieveBuf[0])
            {
               if(!DeserializeChunk(RecieveBuf, Type, Chunk))
                   return 3;
               AddChunkToBlock(Block, Chunk.begin(), Chunk.end());
               Chunk.clear();
            }
        }
            
    MPI_Cart_coords(Cartesian, rank, 2, MyCoords);
    if(GetChunkDestinatedToXY(MyCoords[0], MyCoords[1], rank, P, MaxX, MaxY, N, M, Columns, ChunkStart, ChunkEnd))
        AddChunkToBlock(Block, ChunkStart, ChunkEnd);
    delete [] RecieveBuf;
    delete [] SendBuf;
}

int MulLine(Line &v1, Line &v2)
{
    Line res;
    std::set_intersection( v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter( res, res.begin() ) );
    return res.size();

}

int LineFromIntBuf(int Buf[], Ind Offset, int BufSize, Line &v1)
{
   v1.clear();
   for(int i = 0; i < BufSize; i++)
   {
       v1.insert(Offset+Buf[i]);
   }

    return 0;
}

