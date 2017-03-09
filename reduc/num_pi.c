#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "zeta.h"

int main(int argc, char **argv) {
	if (argc != 2) return 1;

	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		double log_size = log(size)/log(2);
		if (floor(log_size) != ceil(log_size)) {
			MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
		}
	}

	unsigned long int n_tot, n, l, u;
	n_tot = atoi(argv[1]);
	n = ceil((double) n_tot/size); // length of each vector
	l = rank*n + 1;
	u = (rank+1)*n;
	double sum;
	if (l > n_tot) {
		l = u = 0;
		sum = 0;
	} else {
		if (u > n_tot) {
			u = n_tot;
		}
		sum = zeta(l, u);
	}

	double num_pi; 
	MPI_Allreduce(&sum, &num_pi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	num_pi = sqrt(6*num_pi);
	
	if (rank == 0) {
		printf("Accuracy: %.16g.\n", fabs(M_PI - num_pi));
	}

	MPI_Finalize();
	return 0;
}