#ifndef MK2DARRAY_H
#define MK2DARRAY_H

#include <cstdlib>

template <typename T>
T** mkDArray(size_t, size_t);

template <typename T>
void del2DArray(T**);

#include "mk2DArray.cpp"

#endif // MK2DARRAY_H