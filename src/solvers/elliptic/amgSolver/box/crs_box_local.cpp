#include <cstdio>
#include <cstdlib>

#include <lapacke.h>
#include <platform.hpp>

#include "crs_box_impl.hpp"
#include "gemv.h"

static uint gs_n;
static occa::memory o_gs_off, o_gs_idx;

static void setup_gather_rhs(int un, int *u2c) {
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

static int nr = 0;
static struct gemv_t *gemv = NULL;

static void setup_gemv(const struct csr *A) {
  nr = A->nr;

  double *B = tcalloc(double, A->nr * A->nr);
  for (uint i = 0; i < A->nr; i++) {
    for (uint j = A->offs[i]; j < A->offs[i + 1]; j++)
      B[i * A->nr + A->cols[j] - A->base] = A->vals[j];
  }

  int *ipiv = tcalloc(int, A->nr);
  int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, A->nr, A->nr, B, A->nr, ipiv);
  if (info != 0) {
    fprintf(stderr, "dgetrf failed !\n"), fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, A->nr, B, A->nr, ipiv);
  if (info != 0) {
    fprintf(stderr, "dgetri failed !\n"), fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  free(ipiv);

  double *A_inv = tcalloc(double, A->nr * A->nr);
  for (uint i = 0; i < A->nr; i++) {
    for (uint j = 0; j < A->nr; j++)
      A_inv[i * A->nr + j] = B[i * A->nr + j];
  }

  gemv = gemv_init(NULL, NULL);
  gemv_set_matrix(gemv, A_inv);
  gemv_set_backend(gemv, platform->device.mode().c_str());

  free(B), free(A_inv);
}

static int initialized = 0;
static float *h_r = NULL, *h_x = NULL;
static occa::memory o_r, o_x, o_cx;

void asm1_setup(struct csr *A, unsigned null_space, struct box *box) {
  // Setup local assembly of the rhs.
  setup_gather_rhs(box->un, box->u2c);

  // Setup the gemv.
  setup_gemv(A);

  // Allocate work arrays.
  h_r = (float *)calloc(A->nr, sizeof(float));
  h_x = (float *)calloc(A->nr, sizeof(float));
  o_cx = platform->device.malloc(A->nr * sizeof(float));
  o_x = platform->device.malloc(A->nr * sizeof(float));
  o_r = platform->device.malloc(A->nr * sizeof(float));

  initialized = 1;
}

void asm1_solve(float *x, struct box *box, occa::memory &o_r) {
  if (!initialized) {
    fprintf(stderr, "asm1 is not initialized !\n"), fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if (gs_n == 0)
    return;

  platform->gatherRHSKernel(gs_n, o_gs_off, o_gs_idx, o_r, o_cx);

  // hipblasSgemv(handle, HIPBLAS_OP_T, nr, nr, &one_f32, d_A_inv_f32, nr,
  //              (float *)o_cx.ptr(), 1, &zero_f32, (float *)d_x, 1);
  // check_hip_runtime(
  //     hipMemcpy(h_x, d_x, nr * sizeof(float), hipMemcpyDeviceToHost));

  gemv_run(o_x.ptr(), gemv, o_cx.ptr());
  gemv_copy(h_x, o_x.ptr(), sizeof(float) * nr, GEMV_D2H);

  for (uint i = 0; i < box->un; i++) {
    if (box->u2c[i] >= 0)
      x[i] = h_x[box->u2c[i]];
    else
      x[i] = 0;
  }

  for (uint i = box->un; i < box->sn; i++)
    x[i] = 0;
}

void asm1_solve(float *x, struct box *box, const float *r) {
  if (!initialized) {
    fprintf(stderr, "asm1 is not initialized !\n"), fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if (gs_n == 0)
    return;

  for (uint i = 0; i < nr; i++)
    h_r[i] = 0;
  for (uint i = 0; i < box->sn; i++) {
    if (box->u2c[i] >= 0)
      h_r[box->u2c[i]] += r[i];
  }

  // check_hip_runtime(
  //     hipMemcpy(d_r, h_r, nr * sizeof(float), hipMemcpyHostToDevice));

  // hipblasSgemv(handle, HIPBLAS_OP_T, nr, nr, &one_f32, d_A_inv_f32, nr,
  //              (float *)d_r, 1, &zero_f32, (float *)d_x, 1);

  // check_hip_runtime(
  //     hipMemcpy(h_x, d_x, nr * sizeof(float), hipMemcpyDeviceToHost));

  gemv_copy(o_r.ptr(), h_r, sizeof(float) * nr, GEMV_H2D);
  gemv_run(o_x.ptr(), gemv, o_r.ptr());
  gemv_copy(h_x, o_x.ptr(), sizeof(float) * nr, GEMV_D2H);

  for (uint i = 0; i < box->sn; i++) {
    if (box->u2c[i] >= 0)
      x[i] = h_x[box->u2c[i]];
    else
      x[i] = 0;
  }
}

void asm1_free(struct box *box) {
  free(h_r), h_r = NULL;
  free(h_x), h_x = NULL;
  gs_n = nr = initialized = 0;
  gemv_finalize(&gemv);
}
