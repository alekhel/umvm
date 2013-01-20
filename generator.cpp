#include "umvm.h"
#include "math.h"
#include <stdio.h>
#include <time.h>
#ifdef C11RANDOM
#include <random>
#else
#include <cstdlib>
#endif



#ifdef C11RANDOM
int GenerateStripRowwise(Ind StartRow, Ind EndRow, Ind StartColumn, Ind EndColumn, 
                          int MinWeight, int MaxWeight, Matrix &m)
/*Attention, columns indexes are not sorted by default*/
/*The first version does not use MaxWeight parameter*/
{
    
    m.clear();
    int H = EndRow - StartRow;
    std::default_random_engine generator(4);
    std::uniform_int_distribution<unsigned int> distribution(0, EndColumn - StartColumn);
    for (Ind i = StartRow; i < EndRow; i++)
    {
        for(int j = 0; j < MinWeight; j++)
        {
            m[i].insert( StartColumn + distribution(generator));
        }     
    }
    return 0;
}
#else
int GenerateStripRowwise(Ind StartRow, Ind EndRow, Ind StartColumn, Ind EndColumn, 
                          int MinWeight, int MaxWeight, Matrix &m)
/*Attention, columns indexes are not sorted by default*/
/*The first version does not use MaxWeight parameter*/
{
    
    m.clear();
    int H = EndRow - StartRow;
    srand(time(NULL)^StartRow); 
    for (Ind i = StartRow; i < EndRow; i++)
    {
        for(int j = 0; j < MinWeight; j++)
        {
            m[i].insert( StartColumn +  floor((EndColumn - StartColumn)*(double(rand())/RAND_MAX)));
        }     
    }
    return 0;
}
#endif
