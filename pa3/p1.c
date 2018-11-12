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
	srand(0);
	double x[128];
	double A[128][128];
	int i, j;
	for(i = 0; i < 128; i++)
	{
		x[i] = double_rand(-1.0, 1.0);
		for(j = 0; j < 128; j++)
		{
			A[i][j] = double_rand(-1.0, 1.0);
		}
	}
	return 0;
}
