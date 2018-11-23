#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

#define N (256)
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
	/* shift block matrices of A left */
	MPI_Cart_shift(mesh2D, 0, coords[0], &left, &right); 
	MPI_Sendrecv_replace(A, local*local, MPI_DOUBLE, left, 0, right, 0, mesh2D, &s);

	/* shift block matrices of B up */
	MPI_Cart_shift(mesh2D, 1, coords[1], &up, &down);
	MPI_Sendrecv_replace(B, local*local, MPI_DOUBLE, up, 0, down, 0, mesh2D, &s);
#if 0
	for(i = 0; i < root_p; i++)
	{
		/*
		   int MPI_Cart_shift(MPI_Comm comm, int direction, int disp,
		   int *rank_source, int *rank_dest)
		 */
		MPI_Cart_shift(mesh2D, 0, i, &left, &right);
		/*       int MPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype,
				 int dest, int sendtag, int source, int recvtag, MPI_Comm comm,
				 MPI_Status *status)
		 */
	}
#endif
}

	int
main(int argc, char **argv)
{
	//one dimensional array representation of a 2D array
	double *A, *B, *C;
	int p, rank;
	int coord[2];
	int local;
	MPI_Comm mesh2D;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	int dims[] = {0,0};
	//want cyclic connections for communication
	int qperiodic[] = {1,1};

	//create 2D mesh of processors
	MPI_Dims_create(p, 2, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, qperiodic, 0, &mesh2D); 
	MPI_Comm_rank(mesh2D, &rank);
	MPI_Cart_coords(mesh2D, rank, 2, coord);

	srand(0);
	local = N / sqrt(p);

	A = initMat(0, rank, p, mesh2D);
	B = initMat(1, rank, p, mesh2D);
	C = calloc(sizeof(double), local*local);
	/* do the initial setup and do first matmul */
	initial_align(A, B, local, coord, mesh2D);
	matmul(A, B, C, local);



	free(A);free(B);free(C);
	MPI_Finalize();
	return 0;
}
