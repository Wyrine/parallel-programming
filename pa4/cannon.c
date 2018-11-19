#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <mpi.h>

#define N 128

double 
genElt(int i, int j, int rank, int A)
{
	const double x = 0.0001;
	double rv = x + rank*i*x + rank*j*x;
	return (A) ? rv : -1.0 * rv; 
}

int
main(int argc, char **argv)
{
	double A[N][N], B[N][N], C[N][N];
	int world, rank;
	int subsects;
	MPI_Init(NULL, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &world);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	subsects = sqrt(world);
	
	MPI_Finalize();
	return 0;
}
