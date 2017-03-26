// Modified version of the code found at http://stackoverflow.com/a/21944048
#include <cstdlib>

template <typename T>
T** mk2DArray(size_t m, size_t n) {
   T** ptr = new T*[m];  // allocate pointers
   T* pool = new T[m*n];  // allocate pool
   for (size_t i = 0; i < m; ++i, pool += n)
       ptr[i] = pool;
   return ptr;
}

template <typename T>
void del2DArray(T **arr) {
   delete [] arr[0];  // remove the pool
   delete [] arr;     // remove the pointers
}