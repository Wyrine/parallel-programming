#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define X 0.5
#define N 64

//takes the highest power x of left and multiplies it throughout right
	void 
myProd(double *left, double *right, int *len, MPI_Datatype *dtype)
{
	double xpow = left[*len-1];
	int i;
	for(i = 0; i < *len; i++) right[i] *= xpow;
}

// initializes the values of x pre-scan
	void 
initial_work(double *x, int rank, int p)
{
	int i;

	x[0] = X;
	if(rank == 0) x[0] = 1;

	for(i = 1; i < N / p; i++) x[i] = X * x[i-1];
}

	int
main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int rank, p, i, els_per_rank;
	double res;
	MPI_Op op; 

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Op_create(&myProd, 0, &op);	

	els_per_rank = N / p;
	double a[N];
	double x[els_per_rank];
	if(rank == 0)
	{
		srand(10);
		for(i = 0; i < N; i++) a[i] = (double)rand() / (double)RAND_MAX;
	}
	//send the values of a to the other ranks
	MPI_Scatter(a, els_per_rank, MPI_DOUBLE, a, els_per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//get all of the powers of x needed
	initial_work(x, rank, p); 	
	MPI_Scan(x, x, els_per_rank, MPI_DOUBLE, op, MPI_COMM_WORLD);

	//perform dot product of this rank's x's and a's
	for(i = 0, res = 0; i < els_per_rank; i++) res += x[i] * a[i];
	//get the sum of the dot product results on all of the ranks 
	MPI_Allreduce(&res, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

	printf("Rank %d, final result: %.17g\n", rank, res);

	MPI_Op_free(&op);
	MPI_Finalize();	
	return 0;
}
