#include <assert.h>
#include <stddef.h>
#include <stdio.h>

#include "gemv.h"

#define N 8

int main(int argc, char *argv[]) {
  struct gemv_t *handle = gemv_init(&argc, &argv);

  double A[N * N];
  for (unsigned i = 0; i < N * N; i++) A[i] = 1;
  gemv_set_matrix(handle, N, N, A);

  double h_x[N], h_y[N];
  for (unsigned i = 0; i < N; i++) h_x[i] = i;

  gemv_init_session(handle);

  double *x, *y;
  gemv_device_malloc(&x, N), gemv_device_malloc(&y, N);

  gemv_copy(x, h_x, N, GEMV_H2D);

  gemv_run(y, x, handle);

  gemv_copy(h_y, y, N, GEMV_D2H);

  gemv_device_free(&x), gemv_device_free(&y);

  gemv_finalize_session();

  for (unsigned i = 0; i < N; i++) printf("y[%d] = %lf\n", i, h_y[i]);

  gemv_finalize(&handle);

  return 0;
}
