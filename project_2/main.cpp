#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mpi.h>
#include "mk2DArray.h"

double rhs(double x, double y);
extern "C" void fst_(double *v, size_t *n, double *w, size_t *nn);
extern "C" void fstinv_(double *v, size_t *n, double *w, size_t *nn);

void print_matrix(auto **, size_t, size_t);
void print_vector(auto *, size_t);


int main(int argc, char **argv) {
   if (argc != 2) return 1;

   int size, rank;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   size_t m, n, nn;
   n = atoi(argv[1]); // number of internal vertical/horizontal nodes
   m = n/size; // number of horizontal nodes for each process
   nn = 4*(n+1); // Fourier coefficients
   double h = 1.0/(n+1); // constant grid size

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

   // vectors and matrices
   auto *grid_x = new double[m]; 
   auto *grid_y = new double[n]; 
   auto *diag = new double[n];
   auto **b1 = mk2DArray<double>(m,n); // matrix to transform
   auto **b2 = mk2DArray<double>(n,m); // matrix to send
   auto **b3 = mk2DArray<double>(n,m);
   auto *z = new double[nn];

   // x values for the internal nodes
   for (size_t i = 0; i < m; i++) {
      grid_x[i] = (i+1+rank*m)*h;
   }
   // y values for the internal nodes
   for (size_t i = 0; i < n; i++) {
      grid_y[i] = (i+1)*h;
   }

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         b1[i][j] = h*h*rhs(grid_x[i], grid_y[j]);
      }
   }
  
   for (size_t i = 0; i < m; i++) {
      fst_(b1[i], &n, z, &nn);
   }

   // sub-matrices 
   for (size_t i = 0; i < size; i++) {
      for (size_t j = 0; j < m; j++) {
         for (size_t k = 0; k < m; k++) {   
            b2[i*m+j][k] = b1[j][i*m+k];
         }
      }
   }

   // send the parts to the right processes
   MPI_Alltoall(*b2, m*m, MPI_DOUBLE, *b3, m*m, MPI_DOUBLE, MPI_COMM_WORLD);

   // transpose
   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         b1[i][j] = b3[j][i];
      }
   }

   // inverse transform
   for (size_t i = 0; i < m; i++) {
      fstinv_(b1[i], &n, z, &nn);
   }

   
   // The diagonal of the eigenvalue matrix of T
   for (size_t i = 0; i < n; i++) {
      diag[i] = 2.0 * (1.0 - cos((i+1)*M_PI/(n+1)));
   }
   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         b1[i][j] = b1[i][j]/(diag[rank*m+i] + diag[j]);
      }
   }

   for (size_t i = 0; i < m; i++) {
      fst_(b1[i], &n, z, &nn);
   }

   // sub-matrices 
   for (size_t i = 0; i < size; i++) {
      for (size_t j = 0; j < m; j++) {
         for (size_t k = 0; k < m; k++) {   
            b2[i*m+j][k] = b1[j][i*m+k];
         }
      }
   }

   // send the parts to the right processes
   MPI_Alltoall(*b2, m*m, MPI_DOUBLE, *b3, m*m, MPI_DOUBLE, MPI_COMM_WORLD);

   // transpose
   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         b1[i][j] = b3[j][i];
      }
   }

   // inverse transform
   for (size_t i = 0; i < m; i++) {
      fstinv_(b1[i], &n, z, &nn);
   }

   // Calculate maximal value of solution
   double u_max = 0.0;
   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         u_max = u_max > b1[i][j] ? u_max : b1[i][j];
      }
   }

   // Get the maximum from all processes
   MPI_Allreduce(&u_max, &u_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   
   // only one process needs to print the result
   if (rank == 0) {
      printf("u_max = %e\n", u_max);
   }

   del2DArray(b1);
   del2DArray(b2);
   delete [] z;
   delete [] diag;
   delete [] grid_x;
   delete [] grid_y;

   MPI_Finalize();
   return 0;
}

double rhs(double x, double y) {
   return 2 * (y - y*y + x - x*x);
}

void print_matrix(auto **b, size_t m, size_t n) {
   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         std::cout << b[i][j] << " ";
      }
      std::cout << std::endl;
   }
}


void print_vector(auto *b, size_t n) {
   for (size_t i = 0; i < n; i++) {
      std::cout << b[i] << " ";
   }
   std::cout << std::endl;
}
