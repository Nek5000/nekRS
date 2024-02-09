#include <cstdio>
#include <cstdlib>

#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <lapacke.h>
#include <platform.hpp>

#define check_hip_runtime(call)                                                \
  {                                                                            \
    hipError_t err = (call);                                                   \
    if (err != hipSuccess) {                                                   \
      fprintf(stderr, "HIP runtime error: %s\n", hipGetErrorString(err));      \
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);                                 \
    }                                                                          \
  }

#include "crs_box_impl.hpp"

static int initialized = 0;
static hipblasHandle_t handle = NULL;
static int nr = 0;
static double *d_A_inv = NULL;
static float *d_A_inv_f32 = NULL;
static void *d_r = NULL, *d_x = NULL;
static void *h_r = NULL, *h_x = NULL;

static double one = 1.0, zero = 0.0;
static float one_f32 = 1.0f, zero_f32 = 0.0f;

static uint gs_n;
static occa::memory o_gs_off, o_gs_idx;
static occa::memory o_cx;
static occa::kernel gatherRHS;

static void asm1_gs_setup(int un, int *u2c) {
  struct map_t {
    uint u, c;
  };

  struct array map;
  array_init(struct map_t, &map, un);

  struct map_t m;
  for (uint i = 0; i < un; i++) {
    if (u2c[i] < 0)
      continue;
    m.u = i, m.c = u2c[i];
    array_cat(struct map_t, &map, &m, 1);
  }

  if (map.n == 0) {
    gs_n = 0;
    return;
  }

  buffer bfr;
  buffer_init(&bfr, 1024);
  sarray_sort_2(struct map_t, map.ptr, map.n, c, 0, u, 0, &bfr);
  buffer_free(&bfr);

  gs_n = 1;
  struct map_t *pm = (struct map_t *)map.ptr;
  uint c = pm[0].c;
  for (uint i = 1; i < map.n; i++) {
    if (pm[i].c != c) {
      gs_n++;
      c = pm[i].c;
    }
  }

  dlong *gs_off = tcalloc(dlong, gs_n + 1);
  dlong *gs_idx = tcalloc(dlong, map.n);
  gs_off[0] = 0, gs_idx[0] = pm[0].u;
  uint count = 0;
  for (uint i = 1; i < map.n; i++) {
    if (pm[i].c != pm[i - 1].c) {
      count++;
      gs_off[count] = i;
    }
    gs_idx[i] = pm[i].u;
  }
  count++;

  o_gs_off = platform->device.malloc((gs_n + 1) * sizeof(dlong));
  o_gs_off.copyFrom(gs_off, sizeof(dlong) * (gs_n + 1), 0);

  o_gs_idx = platform->device.malloc(map.n * sizeof(dlong));
  o_gs_idx.copyFrom(gs_idx, sizeof(dlong) * map.n, 0);

  free(gs_off), free(gs_idx);

  array_free(&map);
}

void asm1_setup(struct csr *A, unsigned null_space, struct box *box,
                occa::kernel &gatherRHS_) {
  nr = A->nr;
  double *B = tcalloc(double, A->nr * A->nr);
  for (uint i = 0; i < A->nr; i++) {
    for (uint j = A->offs[i]; j < A->offs[i + 1]; j++)
      B[i * A->nr + A->cols[j] - A->base] = A->vals[j];
  }

  int *ipiv = tcalloc(int, A->nr);
  int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, A->nr, A->nr, B, A->nr, ipiv);
  if (info != 0) {
    fprintf(stderr, "dgetrf failed !\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, A->nr, B, A->nr, ipiv);
  if (info != 0) {
    fprintf(stderr, "dgetri failed !\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  free(ipiv);

  double *A_inv = tcalloc(double, A->nr * A->nr);
  for (uint i = 0; i < A->nr; i++) {
    for (uint j = 0; j < A->nr; j++)
      A_inv[i * A->nr + j] = B[i * A->nr + j];
  }

  float *A_inv_f32 = tcalloc(float, A->nr * A->nr);
  for (uint i = 0; i < A->nr; i++)
    for (uint j = 0; j < A->nr; j++)
      A_inv_f32[i * A->nr + j] = (float)B[i * A->nr + j];

  check_hip_runtime(hipMalloc(&d_A_inv, A->nr * A->nr * sizeof(double)));
  check_hip_runtime(hipMemcpy(d_A_inv, A_inv, A->nr * A->nr * sizeof(double),
                              hipMemcpyHostToDevice));

  check_hip_runtime(hipMalloc(&d_A_inv_f32, A->nr * A->nr * sizeof(float)));
  check_hip_runtime(hipMemcpy(d_A_inv_f32, A_inv_f32,
                              A->nr * A->nr * sizeof(float),
                              hipMemcpyHostToDevice));

  free(B), free(A_inv), free(A_inv_f32);

  h_r = calloc(A->nr, sizeof(double));
  h_x = calloc(A->nr, sizeof(double));
  check_hip_runtime(hipMalloc(&d_r, A->nr * sizeof(double)));
  check_hip_runtime(hipMalloc(&d_x, A->nr * sizeof(double)));

  o_cx = platform->device.malloc(A->nr * sizeof(float));
  gatherRHS = gatherRHS_;
  asm1_gs_setup(box->un, box->u2c);

  hipblasCreate(&handle);

  initialized = 1;
}

void asm1_solve(float *x, struct box *box, occa::memory &o_r) {
  if (!initialized) {
    fprintf(stderr, "asm1 is not initialized !\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  gatherRHS(gs_n, o_gs_off, o_gs_idx, o_r, o_cx);
  hipblasSgemv(handle, HIPBLAS_OP_T, nr, nr, &one_f32, d_A_inv_f32, nr,
               (float *)o_cx.ptr(), 1, &zero_f32, (float *)d_x, 1);
  check_hip_runtime(
      hipMemcpy(h_x, d_x, nr * sizeof(float), hipMemcpyDeviceToHost));

  float *h_x_f32 = (float *)h_x;
  for (uint i = 0; i < box->un; i++) {
    if (box->u2c[i] >= 0)
      x[i] = h_x_f32[box->u2c[i]];
    else
      x[i] = 0;
  }

  for (uint i = box->un; i < box->sn; i++)
    x[i] = 0;
}

void asm1_solve_float(float *x, struct box *box, const float *r) {
  float *h_r_ = (float *)h_r;
  for (uint i = 0; i < nr; i++)
    h_r_[i] = 0;
  for (uint i = 0; i < box->sn; i++) {
    if (box->u2c[i] >= 0)
      h_r_[box->u2c[i]] += r[i];
  }

  check_hip_runtime(
      hipMemcpy(d_r, h_r_, nr * sizeof(float), hipMemcpyHostToDevice));

  hipblasSgemv(handle, HIPBLAS_OP_T, nr, nr, &one_f32, d_A_inv_f32, nr,
               (float *)d_r, 1, &zero_f32, (float *)d_x, 1);

  check_hip_runtime(
      hipMemcpy(h_x, d_x, nr * sizeof(float), hipMemcpyDeviceToHost));

  float *h_x_ = (float *)h_x;
  for (uint i = 0; i < box->sn; i++) {
    if (box->u2c[i] >= 0)
      x[i] = h_x_[box->u2c[i]];
    else
      x[i] = 0;
  }
}

void asm1_solve_double(double *x, struct box *box, const double *r) {
  double *h_r_ = (double *)h_r;
  for (uint i = 0; i < nr; i++)
    h_r_[i] = 0;
  for (uint i = 0; i < box->sn; i++) {
    if (box->u2c[i] >= 0)
      h_r_[box->u2c[i]] += r[i];
  }

  check_hip_runtime(
      hipMemcpy(d_r, h_r_, nr * sizeof(double), hipMemcpyHostToDevice));

  hipblasDgemv(handle, HIPBLAS_OP_T, nr, nr, &one, d_A_inv, nr, (double *)d_r,
               1, &zero, (double *)d_x, 1);

  check_hip_runtime(
      hipMemcpy(h_x, d_x, nr * sizeof(double), hipMemcpyDeviceToHost));

  double *h_x_ = (double *)h_x;
  for (uint i = 0; i < box->sn; i++) {
    if (box->u2c[i] >= 0)
      x[i] = h_x_[box->u2c[i]];
    else
      x[i] = 0;
  }
}

void asm1_solve(void *x, struct box *box, const void *r) {
  if (!initialized) {
    fprintf(stderr, "GPU BLAS not initialized.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if (box->dom == gs_float) {
    asm1_solve_float((float *)x, box, (const float *)r);
    return;
  }

  if (box->dom == gs_double) {
    asm1_solve_double((double *)x, box, (const double *)r);
    return;
  }
}

void asm1_free(struct box *box) {
  check_hip_runtime(hipFree(d_A_inv));
  check_hip_runtime(hipFree(d_A_inv_f32));
  check_hip_runtime(hipFree(d_r));
  check_hip_runtime(hipFree(d_x));
  free(h_r), h_r = NULL;
  free(h_x), h_x = NULL;

  nr = 0;
  initialized = 0;

  hipblasDestroy(handle);
}
