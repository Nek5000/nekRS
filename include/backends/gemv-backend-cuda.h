#if !defined(__GEMV_BACKEND_CUDA_H__)
#define __GEMV_BACKEND_CUDA_H__

#include "gemv-backend.h"

#include <cuda_runtime.h>

static inline void check_cuda_runtime_(cudaError_t err, const char *file,
                                       const unsigned line) {
  if (err == cudaSuccess) return;
  fprintf(stderr, "CUDA runtime error: %s in file: %s line: %u\n",
          cudaGetErrorString(err), file, line);
  exit(EXIT_FAILURE);
}

#define check_cuda_runtime(call) check_cuda_runtime_((call), __FILE__, __LINE__)

static inline void cuda_copy(void *dest, const void *src, size_t count,
                             gemv_direction_t direction) {
  enum cudaMemcpyKind kind = cudaMemcpyDefault;
  switch (direction) {
  case GEMV_D2H: kind = cudaMemcpyDeviceToHost; break;
  case GEMV_H2D: kind = cudaMemcpyHostToDevice; break;
  default:
    gemv_log(GEMV_ERROR, "cuda_copy: Invalid direction = %d", direction);
    break;
  }

  check_cuda_runtime(cudaMemcpy(dest, src, count, kind));
}

#endif // __GEMV_BACKEND_CUDA_H__
