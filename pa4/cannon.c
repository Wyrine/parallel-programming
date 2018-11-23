#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

#define N (2)
#define low (0.0)
#define high (1.0)

	double 
getElt()
{
	double scale = rand() / (double) RAND_MAX; /* [0, 1.0] */
	return low + scale * ( high - low );      /* [low, high] */
}

//if negate is set, then we are doing matrix B and so set the indices to be negative
	void 
fill(double *A, int local, int negate)
{
	int i, j;
	for(i = 0; i < local; i++)
		for(j = 0; j < local; j++)
			A[i*local + j] = (negate) ? -1.0 * getElt() : getElt();
}

//pseudo-scatter function -- create the block matrices and then send them to
//the proper processes
	double* 
initMat(int negate, int rank, int p, MPI_Comm mesh2D)
{
	int root_p = sqrt(p);
	int local = N / root_p;
	int send_rank;
	int i,j;
	// 2D array local x local
	double *A = malloc(sizeof(double) * local * local); 

	//if non-zero rank, then just receive the block matrix and return it
	if(rank > 0)
	{
		MPI_Status status;
		MPI_Recv(A, local*local, MPI_DOUBLE, 0, 0, mesh2D, &status);
		return A;
	}
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
			printf("Sending to %d, coords: (%d, %d)\n", send_rank, coords[0], coords[1]);
			if(send_rank > 0)
				MPI_Send(A, local*local, MPI_DOUBLE, send_rank, 0, mesh2D);
		}
	}
	return A;
}

//performs matrix multiplication and stores results into C
	void 
matmul(double *A, double *B, double *C, int local)
{
	int i, j, k;
	for(i=0; i < local; i++)
		for(j = 0; j < local; j++)
			for(k = 0; k < local; k++)
				C[i*local + j] += (A[i*local + k]*B[k*local + j]);
}

	void 
initial_align(double *A, double *B, int local, int coords[], MPI_Comm mesh2D)
{
	MPI_Status s;
	int left, right, up, down;
	/* shift block matrices of A left by row index in mesh */
	MPI_Cart_shift(mesh2D, 1, coords[0], &left, &right); 
	MPI_Sendrecv_replace(A, local*local, MPI_DOUBLE, left, 0, right, 0, mesh2D, &s);

	/* shift block matrices of B up by column index in mesh */
	MPI_Cart_shift(mesh2D, 0, coords[1], &up, &down);
	MPI_Sendrecv_replace(B, local*local, MPI_DOUBLE, up, 0, down, 0, mesh2D, &s);
}
	void
cannon(double *A, double *B, double *C, int local, int p, int root_p, MPI_Comm mesh2D)
{
	int rank;
	MPI_Comm_rank(mesh2D, &rank);
	int i, left, right, up, down;
	for(i = 1; i < root_p; i++)
	{
		MPI_Cart_shift(mesh2D, 1, 1, &left, &right);
		MPI_Sendrecv_replace(A, local*local, MPI_DOUBLE, left, 0, right, 0, mesh2D, MPI_STATUS_IGNORE);

		MPI_Cart_shift(mesh2D, 0, 1, &up, &down);
		MPI_Sendrecv_replace(B, local*local, MPI_DOUBLE, up, 0, down, 0, mesh2D, MPI_STATUS_IGNORE);

		matmul(A, B, C, local);
	}
}

	void
printMat(double *C, int rank)
{
	int i, j;
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
		{
			printf("rank %d, C[%d][%d] = %lf\n", rank, i, j, C[i*N + j]);
		}
}

	int
main(int argc, char **argv)
{
	//one dimensional array representation of a 2D array
	double *A, *B, *C;
	int p, rank, root_p;
	int coord[2];
	int local;
	MPI_Comm mesh2D;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	int dims[] = {0,0};
	//set up cyclic connections between edge processors
	int qperiodic[] = {1,1};

	//create 2D mesh of processors
	MPI_Dims_create(p, 2, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, qperiodic, 1, &mesh2D); 
	MPI_Comm_rank(mesh2D, &rank);
	MPI_Cart_coords(mesh2D, rank, 2, coord);

	srand(0);
	root_p = sqrt(p);
	local = N / root_p;

	A = initMat(0, rank, p, mesh2D);
	B = initMat(1, rank, p, mesh2D);
	if(rank != 0)
		C = calloc(sizeof(double), local*local);
	else //for gather at the end
		C = calloc(sizeof(double), N*N);

	/* do the initial setup and do first matmul */
	initial_align(A, B, local, coord, mesh2D);
	matmul(A, B, C, local);

	cannon(A, B, C, local, p, root_p, mesh2D);
	//gather to proc 0
	double tmp[N*N];
	MPI_Gather(C, local*local, MPI_DOUBLE, tmp, local*local, MPI_DOUBLE, 0, mesh2D);
	if(rank == 0)
		printMat(tmp, rank);
	free(A);free(B);free(C);
	MPI_Finalize();
	return 0;
}
