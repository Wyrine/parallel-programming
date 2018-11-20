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
build_matrices(double ***A, double ***B, double ***C, int rank, int world)
{
	int subsects = sqrt(world);
	int dims = N / subsects;
	int i, j;
	if(rank == 0)
	{
		A[0] = calloc(sizeof(double), N);
		B[0] = calloc(sizeof(double), N);
		C[0] = calloc(sizeof(double), N);
		for(i = 0; i < N; i++)
		{
			A[0][i] = calloc(sizeof(double), N); 
			B[0][i] = calloc(sizeof(double), N); 
			C[0][i] = calloc(sizeof(double), N); 
			for(j = 0; j < N; j++)
			{
				A[0][i][j] = genA(i,j); 
				B[0][i][j] = genB(i,j);
#ifdef DEBUG
				printf("rank: %d, A[%d][%d] = %lf, B[%d][%d] = %lf\n", rank, i, j, A[0][i][j],
						j, i, B[0][i][j]);
#endif
			}
		}
	}
	else
	{
		A[0] = calloc(sizeof(double), dims);
		B[0] = calloc(sizeof(double), dims);
		C[0] = calloc(sizeof(double), dims);
		for(i = 0; i < dims; i++)
		{
			A[0][i] =  calloc(sizeof(double), dims);
			B[0][i] =  calloc(sizeof(double), dims);
			C[0][i] =  calloc(sizeof(double), dims);
		}
	}
	MPI_Scatter(
	return dims;
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

	double **A, **B, **C;
	int dims = build_matrices(&A, &B, &C, rank, world);
	printf("rank %d, val: %lf\n", rank, A[0][0]);

	MPI_Finalize();
	return 0;
}
