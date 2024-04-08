#include "backends/gemv-backend-cuda.h"

static void *d_A = NULL;
static unsigned n = 0, m = 0;
static int initialized = 0;

static void cuda_init(const struct gemv_t *gemv) {
  check_cuda_runtime(cudaSetDevice(gemv->device));

  n = gemv->n, m = gemv->m;
  check_cuda_runtime(cudaMalloc((void **)&d_A, n * m * sizeof(double)));

#if 0
  check_cuda_runtime(
      cudaMemcpy(d_A, A, n * n * sizeof(float), cudaMemcpyHostToDevice));
#endif

  initialized = 1;
}

static void cuda_gemv(void *d_y, const void *d_x, const struct gemv_t *gemv) {}

static void cuda_finalize(void) {
  check_cuda_runtime(cudaFree(d_A)), d_A = NULL;
}

void gemv_register_cuda(void) {
  gemv_backend_register("cuda", cuda_init, cuda_copy, cuda_gemv, cuda_finalize);
}
