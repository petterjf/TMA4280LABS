#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

int main(int argc, char **argv) {
	if (argc != 2) return 1;

	int size, rank, tag;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	int n, nn, n_tot, m;
	n_tot = atoi(argv[1]]);
	n = n_tot/size;
	m = n - 1;
	nn = 4*n;
	real h = 1.0/n_tot;

	// Check that the number of internal nodes in a direction is a power of two
	double log_n = log(n)/log(2);
	if (floor(log_n) != ceil(log_n)) {
		return 1;
	}

	// Check that the number of processes is a power of two
	double log_size = log(size)/log(2);
	if (floor(log_size) != ceil(log_size)) {
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
	}

	// Grid points
	real *grid = mk_1D_array(n+1, false);
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = (rank*n + i) * h;
    }

    // The diagonal of the eigenvalue matrix of T
    real *diag = mk_1D_array(m, false);
    for (size_t i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((rank*n+i+1) * PI / n_tot));
    }

    // Initialize the right hand side data
    real **b = mk_2D_array(m*size, m, false);
    real **bt = mk_2D_array(m, m*size, false);
    real *z = mk_1D_array(nn, false);
    for (size_t i = 0; i < m*size; i++) {
        for (size_t j = 0; j < m; j++) {
            b[i][j] = h * h * rhs(grid[i], grid[j]);
        }
    }

	MPI_Finalize();
	return 0;
}

real rhs(real x, real y) {
    return 2 * (y - y*y + x - x*x);
}

void transpose(real **bt, real **b, size_t m) {
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = b[j][i];
        }
    }
}

real *mk_1D_array(size_t n, bool zero) {
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

real **mk_2D_array(size_t n1, size_t n2, bool zero) {
    real **ret = (real **)malloc(n1 * sizeof(real *));

    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }

    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}