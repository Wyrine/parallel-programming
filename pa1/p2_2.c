#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define S_VAL (50)

	int
main(int argc, char **argv)
{
	int worldsz, rank, r_val;
	int s_val;
	int *arr = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldsz);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	arr = malloc(sizeof(int) * worldsz);
	s_val = S_VAL + rank;

	MPI_Allgather(&s_val, 1, MPI_INT, arr, 1, MPI_INT,  MPI_COMM_WORLD);
	//gather all values to processor 0 and print them from p0
	for(int i = 0; i < worldsz; i++)
		printf("Problem 2.2, rank is %d, arr[%d] = %d\n", rank, i, arr[i]);

	free(arr);
	arr = NULL;
	MPI_Finalize();
}
