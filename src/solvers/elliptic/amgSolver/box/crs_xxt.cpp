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

static uint initialized = 0;
static uint num_vertices = 0, num_elements = 0;
static int remote_size = 0;
static const MPI_Comm *inter_comm = NULL;
static int *send_counts = NULL;
static int *send_disps = NULL;
static int *element_to_offset = NULL;
static struct xxt *xxt = NULL;
static double *wrk = NULL;

void crs_xxt_setup_inter_comm(uint n, const ulong *id, uint nnz, const uint *Ai,
                              const uint *Aj, const double *A, uint null_space,
                              const MPI_Comm *inter_comm_, gs_dom dom,
                              struct comm *local_comm) {
  // Assertions and check if the communicator is an inter-communicator.
  {
    assert(!initialized);

    int test_inter_comm;
    MPI_Comm_test_inter(*inter_comm_, &test_inter_comm);
    assert(test_inter_comm);

    // Get the local size.
    int size;
    MPI_Comm_size(*inter_comm_, &size);
    assert(size == 1);

    // Get the remote size. if remote_size == 1, fall back to XXT or CHOLMOD.
    MPI_Comm_remote_size(*inter_comm_, &remote_size);
    if (remote_size == 1) {
      xxt = crs_xxt_setup(n, id, nnz, Ai, Aj, A, null_space, local_comm, dom);
      initialized = 1;
      return;
    }

    inter_comm = inter_comm_;
    send_counts = tcalloc(int, remote_size);
    send_disps = tcalloc(int, remote_size + 1);
  }

  {
    // 0. Let's broadcast to remote processes that we are doing a setup.
    int operation = 0;
    MPI_Bcast(&operation, 1, MPI_INT, MPI_ROOT, *inter_comm);

    // 1. Let's broadcast the number of dofs.
    MPI_Bcast(&n, 1, MPI_UNSIGNED, MPI_ROOT, *inter_comm);

    // 2. Let's broadcast the number of non-zeros.
    MPI_Bcast(&nnz, 1, MPI_UNSIGNED, MPI_ROOT, *inter_comm);

    // 3. Let's broadcast the null space.
    MPI_Bcast(&null_space, 1, MPI_UNSIGNED, MPI_ROOT, *inter_comm);
  }

  // 4. Let's scatter global ids to remote processes and run parRSB to partition
  // the system in order to reduce the non-zeros in XXT.
  num_vertices = nnz / n;
  num_elements = n / num_vertices;
  {
    const uint element_size = num_elements / remote_size;
    const uint element_remainder = num_elements % remote_size;

    for (uint i = 0; i < remote_size; i++) {
      const uint num_elements_i = element_size + (i < element_remainder);
      send_counts[i] = num_elements_i * num_vertices;
      send_disps[i + 1] = send_disps[i] + send_counts[i];
    }

    long long *send_buffer = tcalloc(long long, n);
    for (uint i = 0; i < n; i++)
      send_buffer[i] = id[i];

    MPI_Scatterv(send_buffer, send_counts, send_disps, MPI_LONG_LONG_INT, NULL,
                 0, MPI_LONG_LONG_INT, MPI_ROOT, *inter_comm);
    free(send_buffer);
  }

  // 5. Let's get the partition info from workers. We will have to adjust
  // send_counts and send_disps.
  int *partition = tcalloc(int, num_elements);
  {
    for (uint i = 0; i < remote_size; i++) {
      send_counts[i] /= num_vertices;
      send_disps[i] /= num_vertices;
    }

    MPI_Gatherv(NULL, 0, MPI_INT, partition, send_counts, send_disps, MPI_INT,
                MPI_ROOT, *inter_comm);
  }

  // 6. Now let's calculate the send_counts and send_disps based on the new
  // partition. Then scatter the send_counts (coarse_size) to the workers.
  {
    for (uint i = 0; i < remote_size; i++) {
      send_counts[i] = 0;
      send_disps[i] = 0;
    }

    for (uint i = 0; i < num_elements; i++) {
      assert(partition[i] < remote_size);
      assert(partition[i] >= 0);
      send_counts[partition[i]] += num_vertices;
    }

    MPI_Scatter(send_counts, 1, MPI_INT, NULL, 0, MPI_INT, MPI_ROOT,
                *inter_comm);

    for (uint i = 0; i < remote_size; i++)
      send_disps[i + 1] = send_disps[i] + send_counts[i];
  }

  // 7. Create element to process map based on the partition.
  {
    struct map_t {
      uint idx, p;
    };

    struct array map;
    array_init(struct map_t, &map, num_elements);

    for (uint i = 0; i < num_elements; i++) {
      struct map_t m = {.idx = i, .p = partition[i]};
      array_cat(struct map_t, &map, &m, 1);
    }

    buffer bfr;
    buffer_init(&bfr, num_elements * sizeof(struct map_t));
    sarray_sort_2(struct map_t, map.ptr, map.n, p, 0, idx, 0, &bfr);
    buffer_free(&bfr);

    element_to_offset = tcalloc(int, num_elements);
    struct map_t *pm = (struct map_t *)map.ptr;
    for (uint i = 0; i < num_elements; i++)
      element_to_offset[pm[i].idx] = i;

    free(partition);
    array_free(&map);
  }

  // 8. Now let's send the global ids and non-zeros to the remote workers.
  {
    size_t max_unit_size = sizeof(long long);
    if (max_unit_size < sizeof(double))
      max_unit_size = sizeof(double);

    void *send_buffer = (void *)calloc(max_unit_size, nnz);

    long long *send_buffer_ll = (long long *)send_buffer;
    for (uint i = 0; i < num_elements; i++) {
      const uint idx = element_to_offset[i] * num_vertices;
      for (uint j = 0; j < num_vertices; j++)
        send_buffer_ll[idx + j] = id[i * num_vertices + j];
    }

    MPI_Scatterv(send_buffer_ll, send_counts, send_disps, MPI_LONG_LONG_INT,
                 NULL, 0, MPI_LONG_LONG_INT, MPI_ROOT, *inter_comm);

    double *send_buffer_dbl = (double *)send_buffer;
    const uint num_vertices_sq = num_vertices * num_vertices;
    for (uint i = 0; i < num_elements; i++) {
      const uint idx = element_to_offset[i] * num_vertices_sq;
      for (uint j = 0; j < num_vertices_sq; j++)
        send_buffer_dbl[idx + j] = A[i * num_vertices_sq + j];
    }

    for (uint i = 0; i < remote_size; i++) {
      send_counts[i] *= num_vertices;
      send_disps[i] *= num_vertices;
    }

    MPI_Scatterv(send_buffer_dbl, send_counts, send_disps, MPI_DOUBLE, NULL, 0,
                 MPI_DOUBLE, MPI_ROOT, *inter_comm);

    for (uint i = 0; i < remote_size; i++) {
      send_counts[i] /= num_vertices;
      send_disps[i] /= num_vertices;
    }

    free(send_buffer);
  }

  // 9. Now let's allocate the work array.
  { wrk = tcalloc(double, n); }

  initialized = 1;
}

void crs_xxt_solve_inter_comm(double *x, const double *rhs) {
  assert(initialized);

  if (remote_size == 1) {
    crs_xxt_solve(x, xxt, rhs);
    return;
  }

  // Let remove processes know that we are doing a solve.
  int operation = 1;
  MPI_Bcast(&operation, 1, MPI_INT, MPI_ROOT, *inter_comm);

  for (uint i = 0; i < num_elements; i++) {
    for (uint j = 0; j < num_vertices; j++)
      wrk[element_to_offset[i] * num_vertices + j] = rhs[i * num_vertices + j];
  }

  // Send the RHS to the remote processes.
  MPI_Scatterv(wrk, send_counts, send_disps, MPI_DOUBLE, NULL, 0, MPI_DOUBLE,
               MPI_ROOT, *inter_comm);

  // Get the solution from the remote processes.
  MPI_Gatherv(NULL, 0, MPI_DOUBLE, wrk, send_counts, send_disps, MPI_DOUBLE,
              MPI_ROOT, *inter_comm);

  for (uint i = 0; i < num_elements; i++) {
    for (uint j = 0; j < num_vertices; j++)
      x[i * num_vertices + j] = wrk[element_to_offset[i] * num_vertices + j];
  }
}

void crs_xxt_finalize_inter_comm() {
  assert(initialized);

  if (remote_size == 1) {
    crs_xxt_free(xxt);
    initialized = 0;
    return;
  }

  // Let's broadcast to remote processes that we are done.
  int operation = 2;
  MPI_Bcast(&operation, 1, MPI_INT, MPI_ROOT, *inter_comm);

  inter_comm = NULL;
  free(send_counts), send_counts = NULL;
  free(send_disps), send_disps = NULL;
  free(element_to_offset), element_to_offset = NULL;
  free(wrk), wrk = NULL;

  initialized = 0;
}
