//Problem 1 - send to right
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define S_VAL (50)

int
main(int argc, char **argv)
{
	int worldsz, rank, r_val;
	int s_val;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldsz);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 	
	
	s_val = S_VAL + rank;

	if(rank == 0)
	{
		MPI_Send(&s_val, 1, MPI_INT, (rank+1) % worldsz, 0, MPI_COMM_WORLD);
		MPI_Recv(&r_val, 1, MPI_INT, worldsz-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("Problem 1.1 - right, rank %d received from %d value %d\n", rank, worldsz-1, r_val);
	}
	else
	{
		MPI_Recv(&r_val, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("Problem 1.1 - right, rank %d received from %d value %d\n", rank, rank-1, r_val);
		MPI_Send(&s_val, 1, MPI_INT, (rank+1) % worldsz, 0, MPI_COMM_WORLD);
	}
	MPI_Finalize();
}
