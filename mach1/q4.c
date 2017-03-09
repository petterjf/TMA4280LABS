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

	if (size < 2) return 1;

	unsigned long n_tot, n;
	n_tot = atoi(argv[1]);
	n = ceil((double) n_tot/(size-1)); // Max length of a vector
	double * vec;
	double sum = 0;

	if (rank == 0) {
		vec = (double *) malloc(sizeof(double)*n_tot);

		for (unsigned long int i = 1; i <= n_tot; i++) {
			vec[i-1] = 4*pow(-1, i-1)*pow(1.0/5.0, 2*i-1)/(2*i-1) - pow(-1, i-1)*pow(1.0/239.0, 2*i-1)/(2*i-1);
		}

		for (int i = 1; i < size; i++) {
			MPI_Send(&vec[(i-1)*n], n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
		}

		double c_sum; // collect sum
		for (int i = 1; i < size; i++) {
			MPI_Recv(&c_sum, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
			sum += c_sum;
		}

		double num_pi = 4*sum;
		printf("%.16g\n", num_pi);
	} else {
		vec = (double *) malloc(sizeof(double)*n);

		MPI_Recv(vec, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		for (unsigned long int i = 0; i < n; i++) {
			sum += vec[i];
		}
		MPI_Send(&sum, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}