#include <map>
#include <vector>
// #define C11RANDOM // uncomment if your compiler supports -std=c++11
typedef unsigned long int Ind;
typedef std::map < Ind, std::vector <Ind> > Matrix; 

int GenerateStripRowwise(Ind StartRow, Ind EndRow, Ind StartColumn, Ind EndColumn, 
                          int MinWeight, int MaxWeight, Matrix &m);

void TryGenerateStripRowwise(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);

int ParameterSanityCheck(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);
