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

//if negate is set, then we are doing matrix B and so set the indices to be negative
void fill(double *A, int local, int negate)
{
	int i, j;
	for(i = 0; i < local; i++)
		for(j = 0; j < local; j++)
			A[i*local + j] = (negate) ? -1.0 * getElt() : getElt();
}

//pseudo-scatter function
double* initMat(int negate, int rank, int p, MPI_Comm mesh2D)
{
	int root_p = sqrt(p);
	int local = N / root_p;
	int send_rank;
	double *A = malloc(sizeof(double) * local * local); // 2D array local x local

	//if non-zero rank, then just receive the block matrix and return it
	if(rank > 0)
	{
		MPI_Status status;
		MPI_Recv(A, local*local, MPI_DOUBLE, 0, 0, mesh2D, &status);
		return A;
	}
	int i,j;
	//otherwise, create root_p x root_p submatrices
	for(i = root_p -1; i >= 0; i--)
	{
		for(j = root_p -1; j >= 0; j--)
		{
			int coords[] = {i, j};
			fill(A, local, negate);
			//get the rank of the processor at (i, j)
			MPI_Cart_rank(mesh2D, coords, &send_rank);
			//and send to that processor unless it's me
			if(send_rank > 0)
				MPI_Send(A, local*local, MPI_DOUBLE, send_rank, 0, mesh2D);
		}
	}
	//fill(A, local, negate);
	return A;
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
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, qperiodic, 0, &mesh2D); 
	MPI_Comm_rank(mesh2D, &rank);
	MPI_Cart_coords(mesh2D, rank, 2, coord);

	srand(0);
	local = N / sqrt(p);
	double *A, *B, *C;

	A = initMat(0, rank, p, mesh2D);
	B = initMat(1, rank, p, mesh2D);
	C = calloc(sizeof(double), local*local);

	free(A);free(B);free(C);
	MPI_Finalize();
	return 0;
}
