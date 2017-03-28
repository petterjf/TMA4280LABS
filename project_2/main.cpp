#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mpi.h>

double rhs(double x, double y);
double anal_soln(double x, double y);
bool is_pow2(double x);
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

   int size, rank;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   size_t m, n, nn;
   n = atoi(argv[1]); // number of internal vertical/horizontal nodes
   m = n/size; // number of vertical nodes for each process
   nn = 4*n; // Fourier coefficients
   double h = 1.0/n; // constant grid size

   // check that the number of processes and internal nodes is a power of two
   if (!(is_pow2(size) && is_pow2(n))) {
      MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
   }

   // vectors and matrices
   auto *grid_x = new double[m]; 
   auto *grid_y = new double[n]; 
   auto *diag = new double[n];
   auto *z = new double[nn];
   auto **soln = mk_2D_array<double>(m,n); // the analytical solution
   auto **b1 = mk_2D_array<double>(m,n); // matrix to transform
   auto **b2 = mk_2D_array<double>(n,m); // matrix to send
   auto **b3 = mk_2D_array<double>(n,m); // receive matrix

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
         soln[i][j] = anal_soln(grid_x[i], grid_y[j]);
      }
   }
  
   for (size_t i = 0; i < m; i++) {
      fst_(b1[i], &n, z, &nn);
   }

   // sub-matrices 
   for (size_t i = 0; i < (size_t)size; i++) {
      for (size_t j = 0; j < m; j++) {
         for (size_t k = 0; k < m; k++) {   
            b2[i*m+j][k] = b1[j][i*m+k];
         }
      }
   }

   // send the parts to the right processes
   MPI_Alltoall(*b2, m*m, MPI_DOUBLE, *b3, m*m, MPI_DOUBLE, MPI_COMM_WORLD);

   transpose(b1, b3, n, m);

   // inverse transform
   for (size_t i = 0; i < m; i++) {
      fstinv_(b1[i], &n, z, &nn);
   }

   // The diagonal of the eigenvalue matrix of T
   for (size_t i = 0; i < n; i++) {
      diag[i] = 2.0 * (1.0 - cos((i+1)*M_PI/n));
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
   for (size_t i = 0; i < (size_t)size; i++) {
      for (size_t j = 0; j < m; j++) {
         for (size_t k = 0; k < m; k++) {   
            b2[i*m+j][k] = b1[j][i*m+k];
         }
      }
   }

   // send the parts to the right processes
   MPI_Alltoall(*b2, m*m, MPI_DOUBLE, *b3, m*m, MPI_DOUBLE, MPI_COMM_WORLD);

   transpose(b1, b3, n, m);

   // inverse transform
   for (size_t i = 0; i < m; i++) {
      fstinv_(b1[i], &n , z, &nn);
   }

   // error for the process
   double loc_max_error = 0;

   // exclude the last row if we are at the last element, because that row only contains gibberish
   if (rank == size-1) {
      m--;
   }

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n-1; j++) {
         double error = fabs(b1[i][j] - soln[i][j]);
         loc_max_error = loc_max_error > error ? loc_max_error : error;
      }
   }

   // get the maximum from all processes
   double max_error;
   MPI_Allreduce(&loc_max_error, &max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   
   // only one process needs to print the result
   if (rank == 0) {
      std::cout << "max_error = " << std::scientific << max_error << std::endl;
   }

   del_2D_array(b1);
   del_2D_array(b2);
   del_2D_array(b3);
   del_2D_array(soln);
   delete [] z;
   delete [] diag;
   delete [] grid_x;
   delete [] grid_y;

   MPI_Finalize();
   return 0;
}

double rhs(double x, double y) {
   return 5*M_PI*M_PI*anal_soln(x, y);
}

double anal_soln(double x, double y) {
   return sin(M_PI*x)*sin(2*M_PI*y);
}

bool is_pow2(double x) {
   x = log(x)/log(2);
   return (floor(x) == ceil(x));
}

void transpose(double **bt, double **b, size_t m, size_t n) {
   // m and n are the number of rows and columns, respectively, for b
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