#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

	double
double_rand( double min, double max )
{
	double scale = rand() / (double) RAND_MAX; /* [0, 1.0] */
	return min + scale * ( max - min );      /* [min, max] */
}

	int
main(int argc, char **argv)
{
	srand(0);
	double x[128];
	double b[128];
	double A[128][128];
	int i, j;
	clock_t start = clock();
	for(i = 0; i < 128; i++)
	{
		x[i] = double_rand(-1.0, 1.0);
		for(j = 0; j < 128; j++)
		{
			A[i][j] = double_rand(-1.0, 1.0);
		}
		fprintf(stderr, "x[%d] = %lf\n", i, x[i]);
	}
	printf("solution vector:\n");
	for(i = 0; i < 128; i++)
	{	
		b[i] = 0.0;
		for(j = 0; j < 128; j++)
		{
			b[i] += (A[i][j] * x[j]);
		}
		printf("\tb[%d] = %lf\n", i, b[i]);
	}
	printf("Total time: %lf\n", (clock() - start) / (double) CLOCKS_PER_SEC);
	return 0;
}
