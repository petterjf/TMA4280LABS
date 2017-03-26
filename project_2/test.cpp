#include <mpi.h>

int main(int argc, char **argv) {
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int *fu = new int[size];
	int *fu2 = new int[size];
	for (int i = 0; i < size; i++) {
		fu[i] = rank;
	}

	MPI_Alltoall(fu, 1, MPI_INT, fu2, 1, MPI_INT, MPI_COMM_WORLD);

	for (int i = 0; i < size; i++) {
		printf("%d ", fu2[i]);
	}
	printf("fra %d\n", rank);

	MPI_Finalize();
	return 0;
}