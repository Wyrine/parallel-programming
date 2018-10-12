#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define X 0.5
void myProd(double *left, double *right, int *len, MPI_Datatype *dtype)
{
	right[0] = right[0] *X + left[1];
	right[1] = right[1] *X + left[1];
	left[1] += right[1];
}

double initial_work(double *a, int rank, int p)
{
	double prev_x, total;
	int i;
	int n = 64;
	if(rank == 0)
	{
		prev_x = 1;
		total = a[0];
	}
	else
	{
		prev_x = X;
		total = a[n/p * rank];
	}
	for(i = 1; i < n / p; i++)
	{
		prev_x *= X;
		total = total + a[n/p * rank + i] * prev_x;	
	}
	return total;
}

int
main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int j, rank, p, i, els_per_rank;
	double res[2], tmp[2];
	MPI_Op op; 
	
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Op_create(&myProd, 0, &op);	

	els_per_rank = 64 / p;
	double a[els_per_rank];
	double a_to_send[els_per_rank];

	if(rank == 0)
	{
		srand(1);
		for(j = 0; j < p; j++)
		{
			for(i = 0; i < els_per_rank; i++)
			{
				if(j == 0)
					a[i] = (double)rand() / (double)RAND_MAX;
				else
					a_to_send[i] = (double)rand() / (double)RAND_MAX;
			}
			if(j != 0)
				MPI_Send(a_to_send, els_per_rank, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
		}
	}
	else
		MPI_Recv(a, els_per_rank, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	res[0] = initial_work(a, rank, p); 	
	res[1] = res[0];
	MPI_Scan(res, tmp, 2, MPI_DOUBLE, op, MPI_COMM_WORLD);
	printf("Rank %d, Final result: %lf\n", rank, tmp[1]);
	/*if(rank == p-1)
	{
		printf("Final result: %lf\n", res);
	}*/
	
	MPI_Op_free(&op);
	MPI_Finalize();	
	return 0;
}
