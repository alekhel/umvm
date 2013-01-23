#include "mpi.h"
#include <map>
#include <set>
//#define C11RANDOM // uncomment if your compiler supports -std=c++11
#define DEBUG_PRINT
#define BUFFER_CHECKS // comment only if the speedup is really necessary
typedef unsigned long int Ind;
typedef std::map < Ind, std::set <Ind> > Matrix; 

int ParameterSanityCheck(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);

int CountElements(Matrix &m);

int WriteMatrixToFolder(char name[], Matrix &m);

int ReadMatrixFromFolder(char name[], Matrix &m);

int GenerateMatrix(Matrix &MyBlock, MPI_Comm &Cartesian, int P, Ind MaxX, Ind MaxY, 
                   int MinWeight, int MaxWeight, Ind N, Ind M);

