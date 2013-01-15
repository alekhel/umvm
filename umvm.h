#include <map>
#include <vector>
typedef unsigned long int Ind;
typedef std::map < Ind, std::vector <Ind> > Matrix; 
int GenerateStripRowwise(Ind StartRow, Ind EndRow, Ind StartColumn, Ind EndColumn, 
                          int MinWeight, int MaxWeight, Matrix &m);
