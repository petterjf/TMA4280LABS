#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mpi.h>
#include <omp.h>

double rhs(double x, double y);
double anal_soln(double x, double y);
bool is_pow2(size_t n);
void transpose(double **bt, double **b, size_t m, size_t n);

template <typename T>
T** mk_2D_array(size_t m, size_t n);

template <typename T>
void del_2D_array(T **arr);

extern "C" {
   void fst_(double *v, size_t *n, double *w, size_t *nn);
   void fstinv_(double *v, size_t *n, double *w, size_t *nn);
}

int main(int argc, char **argv) {
   if (argc != 2) return 1;

   int size, rank, thrd;
   MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thrd);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   size_t m, n, nn;
   n = atoi(argv[1]); // number of internal vertical/horizontal nodes
   m = n/size; // number of vertical nodes for each process
   nn = 4*n; // Fourier coefficients
   double h = 1.0/n; // constant grid size

   // check that the number of processes and internal nodes is a power of two
   if (!(is_pow2(size) && is_pow2(n))) {
      std::cout << "n and P both need to be powers of two\n";
      MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
   }

   double time;
   if (rank == 0) {
      time = MPI_Wtime();
   }

   // vectors and matrices
   double *grid_x = new double[m]; 
   double *grid_y = new double[n]; 
   double *diag = new double[n];
   double **soln = mk_2D_array<double>(m,n); // the analytical solution
   double **b1 = mk_2D_array<double>(m,n); // matrix to transform
   double **b2 = mk_2D_array<double>(n,m); // matrix to send
   double **b3 = mk_2D_array<double>(n,m); // receive matrix

   double max_error = 0;
   int num_threads;
   #pragma omp parallel
   {
      num_threads = omp_get_num_threads();
      double * z = new double[nn];

      // x values for the internal nodes
      #pragma omp for
      for (size_t i = 0; i < m; i++) {
         grid_x[i] = (i+1+rank*m)*h;
      }
      // y values for the internal nodes
      #pragma omp for
      for (size_t i = 0; i < n; i++) {
         grid_y[i] = (i+1)*h;
      }

      #pragma omp for collapse(2)
      for (size_t i = 0; i < m; i++) {
         for (size_t j = 0; j < n; j++) {
            b1[i][j] = h*h*rhs(grid_x[i], grid_y[j]);
            soln[i][j] = anal_soln(grid_x[i], grid_y[j]);
         }
      }

      #pragma omp for
      for (size_t i = 0; i < m; i++) {
         fst_(b1[i], &n, z, &nn);
      }

      // sub-matrices
      #pragma omp for collapse(3)
      for (size_t i = 0; i < (size_t)size; i++) {
         for (size_t j = 0; j < m; j++) {
            for (size_t k = 0; k < m; k++) {   
               b2[i*m+j][k] = b1[j][i*m+k];
            }
         }
      }

      // send the parts to the right processes
      #pragma omp master
      MPI_Alltoall(*b2, m*m, MPI_DOUBLE, *b3, m*m, MPI_DOUBLE, MPI_COMM_WORLD);
      
      #pragma omp barrier
      transpose(b1, b3, n, m);

      // inverse transform
      #pragma omp for
      for (size_t i = 0; i < m; i++) {
         fstinv_(b1[i], &n, z, &nn);
      }

      // The diagonal of the eigenvalue matrix of T
      #pragma omp for
      for (size_t i = 0; i < n; i++) {
         diag[i] = 2.0 * (1.0 - cos((i+1)*M_PI/n));
      }

      #pragma omp for collapse(2)
      for (size_t i = 0; i < m; i++) {
         for (size_t j = 0; j < n; j++) {
            b1[i][j] = b1[i][j]/(diag[rank*m+i] + diag[j]);
         }
      }

      #pragma omp for
      for (size_t i = 0; i < m; i++) {
         fst_(b1[i], &n, z, &nn);
      }

      // sub-matrices 
      #pragma omp for collapse(3)
      for (size_t i = 0; i < (size_t)size; i++) {
         for (size_t j = 0; j < m; j++) {
            for (size_t k = 0; k < m; k++) {   
               b2[i*m+j][k] = b1[j][i*m+k];
            }
         }
      }

      // send the parts to the right processes
      #pragma omp master
      MPI_Alltoall(*b2, m*m, MPI_DOUBLE, *b3, m*m, MPI_DOUBLE, MPI_COMM_WORLD);

      #pragma omp barrier
      transpose(b1, b3, n, m);

      // inverse transform
      #pragma omp for
      for (size_t i = 0; i < m; i++) {
         fstinv_(b1[i], &n , z, &nn);
      }

      // exclude the last row if we are at in last process, because that row only contains gibberish
      #pragma omp master
      if (rank == size-1) {
         m--;
      }

      #pragma omp for collapse(2) reduction(max:max_error)
      for (size_t i = 0; i < m; i++) {
         for (size_t j = 0; j < n-1; j++) {
            double error = fabs(b1[i][j] - soln[i][j]);
            max_error = max_error > error ? max_error : error;
         }
      }
      delete [] z;
   }

   // get the maximum from all processes
   MPI_Allreduce(&max_error, &max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

   del_2D_array(b1);
   del_2D_array(b2);
   del_2D_array(b3);
   del_2D_array(soln);
   delete [] diag;
   delete [] grid_x;
   delete [] grid_y;

   // we have to wait on everyone
   MPI_Barrier(MPI_COMM_WORLD);
   // only one process needs to print the result
   if (rank == 0) {
      std::cout << "t = " << num_threads << ", P = " << size << ", n = " << n << std::endl;
      std::cout << "max_error = " << std::scientific << max_error << std::endl;
      time = MPI_Wtime() - time;
      std::cout << "time = " << std::defaultfloat << time << " sec\n";
   }

   MPI_Finalize();
   return 0;
}

double rhs(double x, double y) {
   return 5*M_PI*M_PI*anal_soln(x, y);
}

double anal_soln(double x, double y) {
   return sin(M_PI*x)*sin(2*M_PI*y);
}

bool is_pow2(size_t n) {
   return (n != 0) && ((n & (n - 1)) == 0);
}

void transpose(double **bt, double **b, size_t m, size_t n) {
   // m and n are the number of rows and columns, respectively, for b
   #pragma omp for collapse(2)
   for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < m; j++) {
         bt[i][j] = b[j][i];
      }
   }
}

// Modified version of the code found at http://stackoverflow.com/a/21944048
template <typename T>
T** mk_2D_array(size_t m, size_t n) {
   T** ptr = new T*[m];  // allocate pointers
   T* pool = new T[m*n];  // allocate pool
   for (size_t i = 0; i < m; ++i, pool += n)
       ptr[i] = pool;
   return ptr;
}

template <typename T>
void del_2D_array(T **arr) {
   delete [] arr[0];  // remove the pool
   delete [] arr;     // remove the pointers
}