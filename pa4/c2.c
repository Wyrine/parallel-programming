#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

#define N (256)
#define low (0.0)
#define high (1.0)

double getElt()
{
	double scale = rand() / (double) RAND_MAX; /* [0, 1.0] */
	return low + scale * ( high - low );      /* [low, high] */
}

void initMat(double A[], int negate)
{
	int i,j;
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			A[i*N + j] = (negate) ? -1.0 * getElt() : getElt();
}

double **create_mat(int sz)
{
	//do one malloc call so that we can use scatter
	double **A = malloc( sz* (sizeof(double*) + sz * sizeof(double)));
	int i;
	int * const data = A + sz;
	for(i = 0; i < sz; i++)
		A[i] = data + i * sz;
	return A;
}

void printMat(double A[], int r, int c)
{
	int i, j;
}

	int
main(int argc, char **argv)
{
	int p,i,j,k, rank;
	int coord[2];
	int local;
	MPI_Comm mesh2D;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	int dims[] = {0,0};
	int qperiodic[] = {0,0};

	//create 2D mesh of processors
	MPI_Dims_create(p, 2, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, qperiodic, 1, &mesh2D); 
	MPI_Comm_rank(mesh2D, &rank);
	MPI_Cart_coords(mesh2D, rank, 2, coord);

	local = N / sqrt(p);
	//double A[N][N], B[N][N], C[N][N];
	double *A, *B, *C;
	if(rank == 0) 
	{
		srand(0);
		A = malloc(N*N*sizeof(double));
		initMat(A, 0);
		B = malloc(N*N*sizeof(double));
		initMat(B, 1);
		C = calloc(sizeof(double), N*N);
	}
	else
	{
		A = malloc(local*local*sizeof(double));
		B = malloc(local*local*sizeof(double));
		C = calloc(sizeof(double), local*local);
	}
	MPI_Scatter(A, local*local, MPI_DOUBLE, A, local*local, MPI_DOUBLE, 0, mesh2D);
	MPI_Scatter(B, local*local, MPI_DOUBLE, B, local*local, MPI_DOUBLE, 0, mesh2D);
	if(rank == 0)
	{
		printf("rank: %d, (%d, %d)\n", rank, coord[0], coord[1]);
/*
		for(i = 0; i < local; i++)
			for(j = local; j < N; j++)
				fprintf(stderr, "A[%d][%d] = %lf\n", i, j%local, A[i*N + j]);
*/
	}
	if(rank == 1)
	{
		printf("rank: %d, (%d, %d)\n", rank, coord[0], coord[1]);
/*
		for(i = 0; i < local; i++)
			for(j = 0; j < local; j++)
				fprintf(stdout, "A[%d][%d] = %lf\n", i, j, A[i*local + j]);
*/
	}

	free(A);free(B);free(C);
	MPI_Finalize();
	return 0;
}
