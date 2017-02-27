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
		clock_t time = clock();

		double log_size = log(size)/log(2);
		if (floor(log_size) != ceil(log_size)) {
			MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
		}
	}

	int n_tot = atoi(argv[1]);

	int n, l, u;
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

	double * sums;
	if (rank == 0) {
		sums = (double *) malloc(sizeof(double)*size);
	}
	MPI_Gather(&sum, 1, MPI_DOUBLE, sums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double num_pi = 0; 
		for (int i = 0; i < size; i++) {
			num_pi += sums[i];
		}

		num_pi = sqrt(6*num_pi);
		double acc = fabs(M_PI - num_pi);
		double time = (clock() - time)/CLOCKS_PER_SEC; 
		printf("Accuracy: %.17g. Time: %f ms.\n", acc, time);
	}
	MPI_Finalize();
	return 0;
}