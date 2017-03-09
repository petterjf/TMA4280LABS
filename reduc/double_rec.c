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
	int tag = 100;
	MPI_Status status;

	double log_size = log(size)/log(2);
	if (rank == 0) {
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

	int q;
	double sum_q;
	for (int d = 0; d < log_size; d++) {
		q = rank ^ (int) pow(2, d);
		MPI_Send(&sum, 1, MPI_DOUBLE, q, tag, MPI_COMM_WORLD);
		MPI_Recv(&sum_q, 1, MPI_DOUBLE, q, tag, MPI_COMM_WORLD, &status);
		sum += sum_q;
	}
	
	double num_pi = sqrt(6*sum);
	if (rank == 0) {
		printf("Accuracy: %.16g.\n", fabs(M_PI - num_pi));
	}

	MPI_Finalize();
	return 0;
}