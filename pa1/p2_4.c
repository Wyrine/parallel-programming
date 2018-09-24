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

	MPI_Allreduce(&s_val, &r_val, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	printf("Problem 2.4 - sum, rank = %d, sum = %d\n", rank, r_val);

	MPI_Allreduce(&s_val, &r_val, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
	printf("Problem 2.4 - product, rank = %d, prod = %d\n", rank, r_val);


	MPI_Finalize();
}
