#ifndef __TESTING_ADIOS2_CUDA_ROUTINES_H__
#define __TESTING_ADIOS2_CUDA_ROUTINES_H__

#include <cuda_runtime.h>

void cuda_increment(int M, int N, int offset, float *vec, float val);

#endif
