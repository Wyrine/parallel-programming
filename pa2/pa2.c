#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

int
main(int argc, char **argv)
{
	int rank;
	int worldsz;
	int i;
	double a[64];
	const double x = 0.5;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldsz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0)
	{
		srand(time(NULL));
		for(i=0; i < 64; i++)
			a[i] = (double)rand() / (double)RAND_MAX;
	}
	MPI_Bcast(a, 64, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	
	
	MPI_Finalize();	
	return 0;
}
