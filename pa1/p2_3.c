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
	
	s_val = S_VAL + rank;

	MPI_Reduce(&s_val, &r_val, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if(rank == 0)
		printf("Problem 2.3 - sum, rank = %d, sum = %d\n", rank, r_val);

	MPI_Reduce(&s_val, &r_val, 1, MPI_INT, MPI_PROD, 0, MPI_COMM_WORLD);
	if(rank == 0)
		printf("Problem 2.3 - product, rank = %d, prod = %d\n", rank, r_val);

	MPI_Finalize();
}
