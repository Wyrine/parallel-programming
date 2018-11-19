#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <mpi.h>

#define N 128
//#define DEBUG

	double 
genA(int i, int j)
{
	const double x = 0.0001;
	return x*(1+i+j); 
}

	double
genB(int i, int j)
{
	const double x = -0.0005;
	return x*(-1-i-j);
}

	void
matmul(double **A, double **B, double **C, int dims)
{
	int i, j, k;
	for(i = 0; i < dims; i++)
		for(j = 0; j < dims; j++)
			for(k = 0; k < dims; k++)
				C[i][j] += (A[i][k] * B[k][j]);
}

	int
main(int argc, char **argv)
{
	int i, j;
	int subsects;
	int world, rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	subsects = sqrt(world);
	int dims = N / subsects;
	double A[dims][dims], B[dims][dims], C[dims][dims];

	if(rank == 0)
	{
		for(i = 0; i < N; i++)
		{
			for(j = 0; j < N; j++)
			{
				A[i][j] = genA(i,j); 
				B[i][j] = genB(i,j);
#ifdef DEBUG
				printf("rank: %d, A[%d][%d] = %lf, B[%d][%d] = %lf\n", rank, i, j, A[i][j],
						j, i, B[i][j]);
#endif
			}
		}
	}
	MPI_Finalize();
	return 0;
}
