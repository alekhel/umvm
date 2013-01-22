#include <iostream>
#include <algorithm>
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
}

void TryGenerateStripRowwise(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight)
{
    Matrix strip;
    int rank;
    int H;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    H = N/P;
    GenerateStripRowwise(rank*H, (rank+1)*H, 0, M, Weight, Weight, strip);
   // for (int i = rank*H; i < (rank+1)*H; i++)
     //   for(unsigned int j = 0;  j < strip[i].size(); j ++ )
       //     printf ("My rank is %d, my  strip[%d][%d] = %7lu\n", rank, i, j,  strip[i][j]);
}

void RowwiseToColumnwise(Matrix Rows, Matrix &Columns)
{
    Columns.clear();
    for(Matrix::iterator it = Rows.begin(); it != Rows.end(); it++)
        {
        for(std::set<Ind>::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++)
            Columns[*it1].insert(it->first);
        }
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
        for(std::set<Ind>::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++, i++)
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
    int ColumnsCount = 0;
    int HeaderSize = 2;
    int Size;
    Chunk.clear();
   
    ColumnsCount = Buf[0];
    if(!ColumnsCount)
        return 0;
    Type = Buf[1];
    Size = Buf[HeaderSize + 2*ColumnsCount - 1]  + HeaderSize + 2*ColumnsCount;
    #ifdef BUFFER_CHECKS
    if(HeaderSize+ColumnsCount*2 >= Size)
    {
        printf("Buffer overflow in DeserializeChunk while reading offsets. Size is %d\n", Size);
        return 0;
    }
    #endif
    int j = 0;
    for(int i = 0; i < ColumnsCount; i++)
    {
        
        while(j++ < Buf[i + HeaderSize + ColumnsCount])
            Chunk[Buf[i+ HeaderSize]].insert(Buf[j + HeaderSize + 2*ColumnsCount - 1]);
        j--;
        #ifdef BUFFER_CHECKS
        if(j >= Size)
        {
            printf("Buffer overflow in DeserializeChunk while reading elements, %d elements read.\n", 
                    j - HeaderSize - 2*ColumnsCount);
            return 0;
        }
        #endif
 
    }
    return Size;
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


int CountElements(Matrix &m)
{
    Matrix::iterator it;
    int res = 0;
    for(it = m.begin(); it != m.end(); it++)
        res+=it->second.size();
    return res;
}

void PrintMatrixStructure(Matrix &m)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for(Matrix::iterator it = m.begin(); it != m.end(); it++)
    {
        printf("[PrintMatrix] My rank is %d, I have elements %lu from  row %lu\n", rank, it->second.size(), it->first);
    }
}
void PrintMatrix(Matrix &m)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for(Matrix::iterator it = m.begin(); it != m.end(); it++)
    {
        printf("[PrintMatrix] My rank is %d, I have elements %lu from  row %lu:", rank, it->second.size(), it->first);
        for(std::set<Ind>::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++ )
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


int GetChunkDestinatedToXY(int X, int Y, int rank, int P, int MaxX, int MaxY, Ind N, Ind M, Matrix &Columns, 
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
    Block.clear();

    StripH = N/P;
    BlockH = N/MaxX;
    BlockW = M/MaxY;
    MaxSendSize = CountElements(Columns) + 2*BlockH*MaxWeight+ 2;
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
                       MPI_Send(SendBuf, size, MPI_INT, j, 0 , MPI_COMM_WORLD);
                    }
                    else                            
                    {
                       // int SendBuf[1];
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


