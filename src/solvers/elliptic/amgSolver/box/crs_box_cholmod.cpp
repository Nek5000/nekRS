#include "crs_box_impl.hpp"
#include <assert.h>

#if defined(ENABLE_CHOLMOD)
#include <cholmod.h>

// Callback function to be called if there is an error.
static void error_handler(int status, const char *file, int line,
                          const char *message) {
  fprintf(stderr, "CHOLMOD error: file: %s line: %d status: %d: %s\n", file,
          line, status, message);
  fflush(stderr);
}

struct cholmod_csr {
  unsigned nr;
  cholmod_common cm;
  cholmod_sparse *A;
  cholmod_factor *L;
  cholmod_dense *r, *X, *Y, *Z;
  gs_dom dom;
};

static int asm1_initialized = 0;

#define TT double
#define SUFFIX _double
#define PRECISION CHOLMOD_DOUBLE
#include "crs_box_cholmod_impl.hpp"
#undef PRECISION
#undef SUFFIX
#undef TT

#define TT float
#define SUFFIX _float
#define PRECISION CHOLMOD_SINGLE
#include "crs_box_cholmod_impl.hpp"
#undef PRECISION
#undef SUFFIX
#undef TT

static void csr_finalize(struct cholmod_csr *B) {
  if (!B)
    return;
  cholmod_free_sparse(&B->A, &B->cm);
  cholmod_free_factor(&B->L, &B->cm);
  cholmod_free_dense(&B->r, &B->cm);
  cholmod_finish(&B->cm);
  free(B);
}

void asm1_cholmod_setup(struct csr *A, unsigned null_space, struct box *box) {
  if (asm1_initialized)
    return;

  switch (box->dom) {
  case gs_double:
    box->ss = (void *)csr_setup_double(A, null_space, &box->global);
    break;
  case gs_float:
    box->ss = (void *)csr_setup_float(A, null_space, &box->global);
    break;
  default:
    break;
  }

  asm1_initialized = 1;
}

void asm1_cholmod_solve(void *x, struct box *box, const void *r) {
  if (!asm1_initialized) {
    fprintf(stderr, "CHOLMOD solve is not initialized.\n");
    exit(EXIT_FAILURE);
  }

  switch (box->dom) {
  case gs_double:
    solve_double((double *)x, box, (double *)r);
    break;
  case gs_float:
    solve_float((float *)x, box, (float *)r);
    break;
  default:
    break;
  }
}

void asm1_cholmod_free(struct box *box) {
  if (asm1_initialized) {
    csr_finalize((struct cholmod_csr *)box->ss), box->ss = NULL;
    asm1_initialized = 0;
  }
}

struct cholmod_csr *sparse_cholmod_factor(uint n, const uint *Arp,
                                          const uint *Aj, const void *A,
                                          gs_dom dom, buffer *bfr) {
  switch (dom) {
  case gs_double:
    return sparse_cholmod_factor_double(n, Arp, Aj, (double *)A, bfr);
    break;
  case gs_float:
    return sparse_cholmod_factor_float(n, Arp, Aj, (float *)A, bfr);
    break;
  default:
    break;
  }
  return NULL;
}

void sparse_cholmod_solve(void *x, struct cholmod_csr *factor, const void *r) {
  switch (factor->dom) {
  case gs_double:
    return sparse_cholmod_solve_double((double *)x, factor, (double *)r);
    break;
  case gs_float:
    return sparse_cholmod_solve_float((float *)x, factor, (float *)r);
    break;
  default:
    break;
  }
}

void sparse_cholmod_free(struct cholmod_csr *factor) {
  switch (factor->dom) {
  case gs_double:
    return sparse_cholmod_free_double(factor);
    break;
  case gs_float:
    return sparse_cholmod_free_float(factor);
    break;
  default:
    break;
  }
}

#else  // ENABLE_CHOLMOD
void asm1_cholmod_setup(struct csr *A, unsigned null_space, struct box *box) {
  fprintf(stderr, "CHOLMOD not enabled !\n");
  exit(EXIT_FAILURE);
}

void asm1_cholmod_solve(void *x, struct box *box, const void *r) {
  fprintf(stderr, "CHOLMOD not enabled !\n");
  exit(EXIT_FAILURE);
}

void asm1_cholmod_free(struct box *box) {
  fprintf(stderr, "CHOLMOD not enabled !\n");
  exit(EXIT_FAILURE);
}

struct cholmod_csr *sparse_cholmod_factor(uint n, const uint *Arp,
                                          const uint *Aj, const void *A,
                                          gs_dom dom, buffer *bfr) {
  fprintf(stderr, "CHOLMOD not enabled !\n");
  exit(EXIT_FAILURE);
  return NULL;
}

void sparse_cholmod_solve(void *x, struct cholmod_csr *factor, const void *r) {
  fprintf(stderr, "CHOLMOD not enabled !\n");
  exit(EXIT_FAILURE);
}

static void sparse_cholmod_free(struct cholmod_csr *factor) {
  fprintf(stderr, "CHOLMOD not enabled !\n");
  exit(EXIT_FAILURE);
}
#endif // ENABLE_CHOLMOD
