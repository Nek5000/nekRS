#include "nekrs_crs.hpp"
#include "crs_box.hpp"

struct crs {
  uint un, type;
  struct comm c;
  gs_dom dom;
  float *wrk;
  void *x, *rhs;
  void *solver;
};

static struct crs *crs = NULL;

void jl_setup(uint type, uint n, const ulong *id, uint nnz, const uint *Ai,
              const uint *Aj, const double *A, uint null, gs_dom dom,
              MPI_Comm comm, const MPI_Comm *inter_comm) {
  if (crs != NULL) {
    if (crs->c.id == 0) {
      fprintf(stderr, "jl_setup: coarse solver is already initialized.\n");
      fflush(stderr);
    }
    return;
  }

  crs = tcalloc(struct crs, 1);

  comm_init(&crs->c, comm);
  crs->type = type;
  crs->un = n;

  crs->dom = dom;
  const char *tmp = getenv("NEKRS_CRS_DOM");
  if (tmp && strncmp(tmp, "gs_double", 32) == 0)
    crs->dom = gs_double;
  if (tmp && strncmp(tmp, "gs_float", 32) == 0)
    crs->dom = gs_float;

  size_t usize;
  switch (dom) {
  case gs_double:
    usize = sizeof(double);
    break;
  case gs_float:
    usize = sizeof(float);
    break;
  default:
    fprintf(stderr, "jl_setup: unknown gs_dom = %d.\n", dom);
    MPI_Abort(comm, EXIT_FAILURE);
    break;
  }

  crs->x = calloc(usize, 2 * n);
  crs->rhs = (void *)((char *)crs->x + n * usize);
  crs->wrk = tcalloc(float, crs->un);

  struct comm *c = &crs->c;
  switch (type) {
  case JL_XXT:
    crs->solver = (void *)crs_xxt_setup(n, id, nnz, Ai, Aj, A, null, c, dom);
    break;
  case JL_BOX:
    crs->solver =
        (void *)crs_box_setup(n, id, nnz, Ai, Aj, A, null, c, inter_comm, dom);
    break;
  default:
    break;
  }
}

#define DOMAIN_SWITCH(dom, macro)                                              \
  {                                                                            \
    switch (dom) {                                                             \
    case gs_double:                                                            \
      macro(double);                                                           \
      break;                                                                   \
    case gs_float:                                                             \
      macro(float);                                                            \
      break;                                                                   \
    }                                                                          \
  }

void jl_solve(occa::memory o_x, occa::memory o_rhs) {
  o_rhs.copyTo(crs->wrk, crs->un * sizeof(float), 0);
#define copy_from_buf(T)                                                       \
  {                                                                            \
    T *rhs = (T *)crs->rhs;                                                    \
    for (uint i = 0; i < crs->un; i++)                                         \
      rhs[i] = crs->wrk[i];                                                    \
  }
  DOMAIN_SWITCH(crs->dom, copy_from_buf);
#undef copy_from_buf

  switch (crs->type) {
  case JL_XXT:
    crs_xxt_solve(crs->x, (struct xxt *)crs->solver, crs->rhs);
    break;
  case JL_BOX:
    crs_box_solve(crs->x, (struct box *)crs->solver, crs->rhs);
    break;
  default:
    break;
  }

#define copy_to_buf(T)                                                         \
  {                                                                            \
    T *x = (T *)crs->x;                                                        \
    for (uint i = 0; i < crs->un; i++)                                         \
      crs->wrk[i] = x[i];                                                      \
  }
  DOMAIN_SWITCH(crs->dom, copy_to_buf);
#undef copy_to_buf
  o_x.copyFrom(crs->wrk, crs->un * sizeof(float), 0);
}

void jl_solve2(occa::memory o_x, occa::memory o_rhs) {
  crs_box_solve2(o_x, (struct box *)crs->solver, o_rhs);
}

void jl_free() {
  if (crs == NULL)
    return;

  switch (crs->type) {
  case JL_XXT:
    crs_xxt_free((struct xxt *)crs->solver);
    break;
  case JL_BOX:
    crs_box_free((struct box *)crs->solver);
    break;
  default:
    break;
  }

  comm_free(&crs->c);
  free(crs->x), free(crs->wrk), free(crs), crs = NULL;
}

#undef DOMAIN_SWITCH
