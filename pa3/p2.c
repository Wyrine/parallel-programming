#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

	double
double_rand( double min, double max )
{
	double scale = rand() / (double) RAND_MAX; /* [0, 1.0] */
	return min + scale * ( max - min );      /* [min, max] */
}

	int
main(int argc, char **argv)
{
	const int N = 128;
	double x[N], b[N], A[N][N], tmp[N];
	int i, j, work, rank, world;
	MPI_Init(&argc, &argv);
	srand(0);

	MPI_Comm_size(MPI_COMM_WORLD, &world);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	work = N / world;
	for(i = 0; i < work; i++)
	{
		x[i+rank*work] = double_rand(-1.0, 1.0);
		for(j = 0; j < work; j++)
		{
			A[i+rank*work][j] = double_rand(-1.0, 1.0);
		}
	} 
	MPI_Allgather(x+work*rank, work, MPI_DOUBLE, tmp, work, MPI_DOUBLE,MPI_COMM_WORLD);
	for(i = 0; i < work; i++)
	{	
		b[i+rank*work] = 0.0;
		for(j = 0; j < N; j++)
		{
			b[i+rank*work] += (A[i+rank*work][j] * tmp[j]);
		}
		printf("\trank %d: b[%d] = %lf\n", rank, i+rank*work, b[i+rank*work]);
	}
	MPI_Finalize();
	return 0;
}
