#include "umvm.h"
int GenerateStripRowwise( int Height, int Width, int MinWeight, int MaxWeight, 
                          Matrix &m, Ind RowOffset = 0, Ind ColumnOffset = 0);

void RowwiseToColumnwise(Matrix Rows, Matrix &Columns);

int DistributeMatrixChunks(int P, int MaxX, int MaxY, int MaxWeight, 
                            Ind  N, Ind M, Matrix &Columns, Matrix &Block, MPI_Comm Cartesian);

int SerializeChunk(Matrix::iterator Start, Matrix::iterator End, unsigned int Size, int Type,  int Res[]);

int DeserializeChunk(int Buf[], int &Type, Matrix &Chunk);

void AddChunkToBlock(Matrix &Block, Matrix::iterator Start, Matrix::iterator End);

void PrintMatrixStructure(Matrix &m);

void PrintMatrix(Matrix &m);

int GetChunkDestinatedToXY(int X, int Y, int rank, int P, int MaxX, int MaxY, Ind N, Ind M, 
                           Matrix &Columns, Matrix::iterator &ChunkStart, Matrix::iterator &ChunkEnd);
int GenerateStripRowwise(Ind StartRow, Ind EndRow, Ind StartColumn, Ind EndColumn, 
                          int MinWeight, int MaxWeight, Matrix &m);

int GetNeededBufferSize(Matrix::iterator Start, Matrix::iterator End);


int GetMaxBufferSize(Matrix::iterator Start, Matrix::iterator End);
