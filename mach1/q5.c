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

	double log_size = log(size)/log(2);
	if (floor(log_size) != ceil(log_size)) {
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
	}

	int n_tot = atoi(argv[1]);
	int n = ceil((double) n_tot/size); // length of each vector
	int l = rank*n;
	int u = (rank+1)*n - 1;
	if (l >= n_tot) {
		l = u = 0;
	} else if (u >= n_tot) {
		u = n_tot - 1;
	}
	printf("I am %d, with (%d,%d)\n", rank, l, u);

	MPI_Finalize();
	return 0;
}