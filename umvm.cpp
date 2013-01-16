#include "mpi.h"
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <sys/mman.h>
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
};

void TryGenerateStripRowwise(int P, int MaxX, int MaxY, Ind  N, Ind M, unsigned int Weight)
{
    Matrix strip;
    int rank;
    int H;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    H = N/P;
    GenerateStripRowwise(rank*H, (rank+1)*H, 0, M, Weight, Weight, strip);
    for (int i = rank*H; i < (rank+1)*H; i++)
        for(unsigned int j = 0;  j < strip[i].size(); j ++ )
            printf ("My rank is %d, my  strip[%d][%d] = %7lu\n", rank, i, j,  strip[i][j]);
}
