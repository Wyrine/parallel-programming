#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

	double
double_rand( double min, double max , int i)
{
#if 0
	double scale = rand() / (double) RAND_MAX; /* [0, 1.0] */
	return min + scale * ( max - min );      /* [min, max] */
#endif
	const double x = 0.0001;
	return x + x*i;
}

	int
main(int argc, char **argv)
{
	const int N = 128;
	double x[N], b[N], A[N][N];
	double start, stop;
	int i, j, work, rank, world;
	MPI_Init(&argc, &argv);
	srand(0);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world);
	start = MPI_Wtime();
	work = N / world; 
	double tmp[work];
	for(i = 0; i < work; i++)
	{
		int rel_idx = i+rank*work;
		tmp[i] = double_rand(-1.0, 1.0, rel_idx);
		for(j = 0; j < N; j++)
			A[rel_idx][j] = double_rand(-1.0, 1.0, rel_idx);
	} 
	MPI_Allgather(tmp, work, MPI_DOUBLE, x, work, MPI_DOUBLE,MPI_COMM_WORLD);
	for(i = 0; i < work; i++)
	{	
		int rel_idx = i+rank*work;
		b[rel_idx] = 0.0;
		for(j = 0; j < N; j++)
			b[rel_idx] += (A[rel_idx][j] * x[j]);
		printf("\trank %d: b[%d] = %lf\n", rank, rel_idx, b[rel_idx]);
	}
	printf("Rank: %d, time: %lf\n", rank, MPI_Wtime() - start);
	MPI_Finalize();
	return 0;
}
