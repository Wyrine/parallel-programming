/*  Cannon's Algorithm
Authors: Kirolos Shahat and Grace Zhao
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <string.h>
#include "mpi.h"

#define N (256)
#define low (0.0)
#define high (1.0)

//Generate a random double between 0.0 and 1.0
	double 
getElt()
{
	double scale = rand() / (double) RAND_MAX; /* [0, 1.0] */
	return low + scale * ( high - low );      /* [low, high] */
}

//pseudo-scatter function -- create the block matrices and then send them to
//the proper processes
	double* 
initMat(int negate, int rank, int p, MPI_Comm mesh2D)
{
	int root_p = sqrt(p);
	int local = N / root_p;
	double *A[root_p][root_p];
	int i,j;

	//if non-zero rank, then just receive the block matrix and return it
	if(rank > 0)
	{
		A[0][0] = malloc(sizeof(double) * local*local);
		MPI_Recv(A[0][0], local*local, MPI_DOUBLE, 0, 0, mesh2D, MPI_STATUS_IGNORE);
		return A[0][0];
	}

	/* the rest only occurs on root processor -- processor 0 */	

	/* allocation step */
	for(i = 0; i < root_p; i++)
		for(j = 0; j < root_p; j++)
			A[i][j] = malloc(sizeof(double) * local * local);

	/* fill matrices */
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			A[i / local][j / local][(i%local)*local + j % local] = (negate) ? -1.0 * getElt(): getElt();

	/* scatter */
	int coords[2];
	int send_rank;
	for(i = 0; i < root_p; i++)
	{
		coords[0] = i;
		for(j = 0; j < root_p; j++)
		{
			//if the processor is not the root 
			if(i !=  0 || j != 0)
			{
				coords[1] = j;
				MPI_Cart_rank(mesh2D, coords, &send_rank);
				MPI_Send(A[i][j], local*local, MPI_DOUBLE, send_rank, 0, mesh2D);
				free(A[i][j]);
			}
		}
	}
	//return the submatrix belonging to this processor
	return A[0][0];
}

//performs serial matrix multiplication and stores result into C 
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
	int left, right, up, down;
	/* shift block matrices of A left by row index in mesh */
	MPI_Cart_shift(mesh2D, 1, coords[0], &left, &right); 
	MPI_Sendrecv_replace(A, local*local, MPI_DOUBLE, left, 0, right, 0, mesh2D, MPI_STATUS_IGNORE);

	/* shift block matrices of B up by column index in mesh */
	MPI_Cart_shift(mesh2D, 0, coords[1], &up, &down);
	MPI_Sendrecv_replace(B, local*local, MPI_DOUBLE, up, 0, down, 0, mesh2D, MPI_STATUS_IGNORE);
}

	void
cannon(double *A, double *B, double *C, int local, int p, int root_p, MPI_Comm mesh2D)
{
	int rank;
	MPI_Comm_rank(mesh2D, &rank);
	int i, left, right, up, down;
	for(i = 1; i < root_p; i++)
	{
		//shift A values left by one
		MPI_Cart_shift(mesh2D, 1, 1, &left, &right);
		MPI_Sendrecv_replace(A, local*local, MPI_DOUBLE, left, 0, right, 0, mesh2D, MPI_STATUS_IGNORE);

		//shift B values up by one
		MPI_Cart_shift(mesh2D, 0, 1, &up, &down);
		MPI_Sendrecv_replace(B, local*local, MPI_DOUBLE, up, 0, down, 0, mesh2D, MPI_STATUS_IGNORE);

		//matmul the new results
		matmul(A, B, C, local);
	}
}

	void
gather(double ***C, int local, int p, int root_p, MPI_Comm mesh2D)
{
	int coords[2];
	int i, j, recv_rank;
	/* Gather to root processor */
	for(i = 0; i < root_p; i++)
	{
		coords[0] = i;
		for(j = 0; j < root_p; j++)
		{
			//Do the following work if the rank is not root
			if(i != 0 || j != 0)
			{
				coords[1] = j;
				MPI_Cart_rank(mesh2D, coords, &recv_rank);
				MPI_Recv(C[i][j], local*local, MPI_DOUBLE, recv_rank, 0, mesh2D, MPI_STATUS_IGNORE);
			}
		}
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
	int i, j;
	double t;
	MPI_Comm mesh2D;
	MPI_Init(&argc, &argv);

	t = MPI_Wtime();
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
	//for processor coords and then the block matrix 
	double ***C_Total;
	//Get my local A matrix
	A = initMat(0, rank, p, mesh2D);
	//And local B matrix
	B = initMat(1, rank, p, mesh2D);
	if(rank == 0) //for gathering all at the end
	{
		//set up block matrices of C
		C_Total = malloc(sizeof(double**) * root_p);
		for(i = 0; i < root_p; i++)
		{
			C_Total[i] = malloc(sizeof(double*) * root_p);
			for(j = 0; j < root_p; j++)
				C_Total[i][j] = calloc(sizeof(double), local*local);
		}
		C = C_Total[0][0];
	}
	//otherwise I just need a local block of the matrix
	else C = calloc(sizeof(double), local*local);

	/* Do cannon's algorithm */ 
	initial_align(A, B, local, coord, mesh2D);
	matmul(A, B, C, local);
	cannon(A, B, C, local, p, root_p, mesh2D);

	/* The following does a gather onto processor 0 and then prints C from there */
	if(rank == 0)
	{
		//gather to proc 0
		gather(C_Total, local, p, root_p, mesh2D);
		t = MPI_Wtime() - t;
		//print C matrix 
		for(i = 0; i < N; i++)
			for(j = 0; j < N; j++)
				printf("C[%d][%d] = %lf\n", i, j, C_Total[i / local][j/local][(i%local)*local + (j%local)]);
		printf("Total time taken (excluding printing) = %lf\n", t);
		//Then free the submatrices
		for(i = 0; i < root_p; i++)
		{
			for(j = 0; j < root_p; j++)
				free(C_Total[i][j]);
			free(C_Total[i]);
		}
		free(C_Total);
	}
	else
	{
		//send to rank 0
		MPI_Send(C, local*local, MPI_DOUBLE, 0, 0, mesh2D);
		//and free my submatrix
		free(C);
	}
	//Free the A and B blocks
	free(A); free(B);
	MPI_Finalize();
	return 0;
}
