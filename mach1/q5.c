#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv) {
	if (argc != 2) return 1;

	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;
	int tag = 100;

	int root = size-1; // Make the last rank the root 

	double log_size = log(size)/log(2);
	if (floor(log_size) != ceil(log_size)) {
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
	}

	int n_tot = atoi(argv[1]);
	//if (n_tot < 2*size) {
	//	MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
	//}
	int n = ceil((double) n_tot/size); // length of each vector
	
	// lower and upper limit to sum
	int l = rank*n + 1;
	int u = (rank+1)*n;
	if (u > n_tot) {
		if (l > n_tot) {
			u = l;
		} else {
			u = n_tot;
		}
	}

	printf("I am %d, with (%d,%d)\n", rank, l, u);
	double sum = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	/*if (rank == root) {
		printf("%d, %d\n", n, n_tot);
	}*/

	MPI_Finalize();
	return 0;
}