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

int StoreMatrixToFolder(char *DirName,char *FileNamePrefix, 
                        Matrix &Block, int Type, int P, int MaxX, int MaxY, 
                        int MinWeight, int MaxWeight,  Ind N, Ind M,  
                        MPI_Comm Cartesian);

int LoadMatrixFromFolder(char *DirName,char *FileNamePrefix, 
                        Matrix &Block, int &Type, int &P, int &MaxX, int &MaxY, 
                        int &MinWeight, int &MaxWeight,  Ind &N, Ind &M,  
                        MPI_Comm &Cartesian);

int GenerateMatrix(Matrix &MyBlock, MPI_Comm &Cartesian, int P, int MaxX, int MaxY, 
                    int MinWeight, int MaxWeight, Ind N, Ind M);
