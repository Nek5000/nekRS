#include <cassert>
#include <float.h>
#include <limits.h>
#include <malloc.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "crs_box_impl.hpp"

#if defined(ENABLE_BLAS)
#include <cblas.h>
#include <lapacke.h>
#endif

struct cholmod_csr *fac_A_ll = NULL;
double *A_ll_inv = NULL, *y_inv = NULL;
float *A_ll_inv_f32 = NULL, *y_inv_f32 = NULL;

/*
  portable log base 2

  does a binary search to find leading order bit

  UINT_BITS = number of bits in a uint
  BITS(0) = UINT_BITS
  BITS(i) = half of BITS(i-1), rounded up
  MASK(i) = bitmask with BITS(i) 1's followed by BITS(i) 0's
*/

static unsigned lg(uint v) {
  unsigned r = 0;
#define UINT_BITS (sizeof(uint) * CHAR_BIT)
#define BITS(i) ((UINT_BITS + (1 << (i)) - 1) >> (i))
#define MASK(i) ((((uint)1 << BITS(i)) - 1) << BITS(i))
#define CHECK(i)                                                               \
  if ((BITS(i) != 1) && (v & MASK(i)))                                         \
  v >>= BITS(i), r += BITS(i)
  CHECK(1);
  CHECK(2);
  CHECK(3);
  CHECK(4);
  CHECK(5);
  CHECK(6);
  CHECK(7);
  CHECK(8);
  CHECK(9); /* this covers up to 1024-bit uints */
  if (v & 2)
    ++r;
  return r;
#undef UINT_BITS
#undef BITS
#undef MASK
#undef CHECK
}

struct yale_mat {
  uint i, j;
  double v;
};

/* the tuple list describing the condensed dofs:
   [(separator level, share count, global id)] */
struct dof {
  ulong id;
  uint level, count;
};

static double cholesky_time = 0;
static double xxt_time = 0;
static double qqt_time = 0;
static double local_time = 0;

#define T double
#define SUFFIX _double
#define gs_domain gs_double
#include "crs_xxt_impl.hpp"
#undef T
#undef SUFFIX
#undef gs_domain

#define T float
#define SUFFIX _float
#define gs_domain gs_double
#include "crs_xxt_impl.hpp"
#undef T
#undef SUFFIX
#undef gs_domain

struct xxt {
  gs_dom dom;
  void *solver;
};

struct xxt *crs_xxt_setup(uint n, const ulong *id, uint nz, const uint *Ai,
                          const uint *Aj, const double *A, uint null_space,
                          const struct comm *comm, gs_dom dom) {
  struct xxt *xxt = tcalloc(struct xxt, 1);
  xxt->dom = dom;
  switch (dom) {
  case gs_double:
    xxt->solver =
        (void *)crs_xxt_setup_double(n, id, nz, Ai, Aj, A, null_space, comm);
    break;
  case gs_float:
    xxt->solver =
        (void *)crs_xxt_setup_float(n, id, nz, Ai, Aj, A, null_space, comm);
    break;
  default:
    fprintf(stderr, "Domain %u is not supported.\n", dom);
    exit(EXIT_FAILURE);
    break;
  }

  return xxt;
}

void crs_xxt_solve(void *x, struct xxt *xxt, const void *b) {
  switch (xxt->dom) {
  case gs_double:
    crs_xxt_solve_double((double *)x, (struct xxt_double *)xxt->solver,
                         (const double *)b);
    break;
  case gs_float:
    crs_xxt_solve_float((float *)x, (struct xxt_float *)xxt->solver,
                        (const float *)b);
    break;
  default:
    fprintf(stderr, "Domain %u is not supported.\n", xxt->dom);
    exit(EXIT_FAILURE);
    break;
  }
}

void crs_xxt_stats(struct xxt *xxt) {
  switch (xxt->dom) {
  case gs_double:
    crs_xxt_stats_double((struct xxt_double *)xxt->solver);
    break;
  case gs_float:
    crs_xxt_stats_float((struct xxt_float *)xxt->solver);
    break;
  default:
    fprintf(stderr, "Domain %u is not supported.\n", xxt->dom);
    exit(EXIT_FAILURE);
    break;
  }
}

void crs_xxt_times(double *cholesky_time_, double *local_time_,
                   double *xxt_time_, double *qqt_time_) {
  cholesky_time_[0] = cholesky_time;
  local_time_[0] = local_time;
  xxt_time_[0] = xxt_time;
  qqt_time_[0] = qqt_time;
}

void crs_xxt_free(struct xxt *xxt) {
  switch (xxt->dom) {
  case gs_double:
    crs_xxt_free_double((struct xxt_double *)xxt->solver);
    break;
  case gs_float:
    crs_xxt_free_float((struct xxt_float *)xxt->solver);
    break;
  default:
    fprintf(stderr, "Domain %u is not supported.\n", xxt->dom);
    exit(EXIT_FAILURE);
    break;
  }
  free(xxt);
}
