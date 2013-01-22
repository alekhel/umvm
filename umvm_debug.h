#include "umvm.h"
void TryGenerateStripRowwise(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);
void TryRowwiseToColumnwise(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);
void TrySerializeChunk(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);
void TryDeserializeChunk(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight);
int TryGetAndInsertChunk(int P, int MaxX, int MaxY, int MaxWeight, Ind  N, Ind M, MPI_Comm Cartesian);
