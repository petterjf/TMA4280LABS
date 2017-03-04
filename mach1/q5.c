#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "mach.h"

int main(int argc, char **argv) {
	if (argc != 2) return 1;

	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double time;
	if (rank == 0) {
		time = MPI_Wtime();
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
		sum = 4*mach(1.0/5.0, l, u) - mach(1.0/239.0, l, u);
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
		
		num_pi = 4*num_pi;
		double acc = fabs(M_PI - num_pi);

		time = 1e3*(MPI_Wtime() - time);
		printf("Accuracy: %.15f. Time: %f ms.\n", acc, time);
	}

	MPI_Finalize();
	return 0;
}