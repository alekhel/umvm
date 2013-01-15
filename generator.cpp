#include "umvm.h"
#include <random>

int GenerateStripRowwise(Ind StartRow, Ind EndRow, Ind StartColumn, Ind EndColumn, 
                          int MinWeight, int MaxWeight, Matrix &m)
/*Attention, columns indexes are not sorted by default*/
/*The first version does not use MaxWeight parameter*/
{
    
    m.clear();
    std::default_random_engine generator(4);
    std::uniform_int_distribution<unsigned int> distribution(0, EndColumn - StartColumn);
    for (Ind i = StartRow; i < EndRow; i++)
    {
        for(int j = 0; j < MinWeight; j++)
        {
            m[StartRow + j].push_back( StartColumn + distribution(generator));
        }     
    }
    return 0;
}
