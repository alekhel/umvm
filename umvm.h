#include "mpi.h"
#include <map>
#include <set>
//#define C11RANDOM // uncomment if your compiler supports -std=c++11
#define DEBUG_PRINT
#define BUFFER_CHECKS // comment only if the speedup is really necessary
typedef unsigned long int Ind;
typedef std::map < Ind, std::set <Ind> > Matrix; 

int GenerateStripRowwise(Ind StartRow, Ind EndRow, Ind StartColumn, Ind EndColumn, 
                          int MinWeight, int MaxWeight, Matrix &m);

void TryGenerateStripRowwise(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);

int ParameterSanityCheck(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);

void RowwiseToColumnwise(Matrix Rows, Matrix &Columns);

void TryRowwiseToColumnwise(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);

int DistributeMatrixChunks(int P, int MaxX, int MaxY, int MaxWeight, 
                            Ind  N, Ind M, Matrix &Columns, Matrix &Block, MPI_Comm Cartesian);

int SerializeChunk(Matrix::iterator Start, Matrix::iterator End, unsigned int Size, int Type,  int Res[]);

int DeserializeChunk(int Buf[], int &Type, Matrix &Chunk);

void TrySerializeChunk(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);

void TryDeserializeChunk(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);

void PrintMatrixStructure(Matrix &m);

void PrintMatrix(Matrix &m);

int CountElements(Matrix &m);

int WriteMatrixToFolder(char name[], Matrix &m);

int ReadMatrixFromFolder(char name[], Matrix &m);

void AddChunkToBlock(Matrix &Block, Matrix::iterator Start, Matrix::iterator End);
int TryGetAndInsertChunk(int P, int MaxX, int MaxY, int MaxWeight, Ind  N, Ind M, MPI_Comm Cartesian);
