#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

	int
main(int argc, char **argv)
{
	int rank, worldsz, i, j;
	long long x;
	clock_t start, stop;
	double total;

	char *arr = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldsz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank == 0)
		printf("Problem 1.2 timings\n");

	for(i = 0; i < 30; i++)
	{
		x = 1;
		for(j = 0; j < i; j++) x *= 2;
		arr = calloc(x, 1);

		start = clock();
		if(rank == 0)
		{
			MPI_Send(arr, x, MPI_SIGNED_CHAR, 1, 0, MPI_COMM_WORLD);
			stop = clock();
			total = (double) (stop - start) / CLOCKS_PER_SEC;
			fprintf(stderr, "\tSend message size: %lld, time(secs): %lf\n", x, total);
		}
		else if(rank == 1)
		{
			MPI_Recv(arr, x, MPI_SIGNED_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			stop = clock();
			total = (double) (stop - start) / CLOCKS_PER_SEC;
			fprintf(stderr, "\tRecv message size: %lld, time(secs): %lf\n", x, total);
		}
		free(arr);
	}
	MPI_Finalize();
}
