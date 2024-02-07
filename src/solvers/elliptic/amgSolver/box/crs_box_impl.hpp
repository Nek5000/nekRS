#if !defined(_CRS_BOX_IMPL_HPP_)
#define _CRS_BOX_IMPL_HPP_

#include "crs_box.hpp"
#include "crs_box_csr.hpp"
#include "crs_box_timer.hpp"

#define BOX_DOMAIN_SWITCH(dom, macro)                                          \
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

#define BOX_XXT 1
#define BOX_CHOLMOD 2
#define BOX_GPU_BLAS 4

struct box {
  // User input.
  uint un;
  uint ncr;
  gs_dom dom;
  unsigned algo, mult;
  // Computed data structures for Schwarz.
  struct comm global, local;
  struct gs_data *gsh;
  double *inv_mul;
  uint cn;
  sint *u2c;
  // Work arrays.
  uint sn;
  void *sx, *srhs;
  buffer bfr;
  // Pointer to the asm1 solver.
  void *ss;
};

void box_debug(const int verbose, const char *fmt, ...);

// ASM1: CHOLMOD, redundant API interface
struct cholmod_csr;
struct cholmod_csr *sparse_cholmod_factor(uint n, const uint *Arp,
                                          const uint *Aj, const void *A,
                                          gs_dom dom, buffer *bfr);
void sparse_cholmod_solve(void *x, struct cholmod_csr *factor, const void *r);
void sparse_cholmod_free(struct cholmod_csr *factor);

// ASM1: CHOLMOD API interface.
void asm1_cholmod_setup(struct csr *A, unsigned null_space, struct box *box);
void asm1_cholmod_solve(void *x, struct box *box, const void *r);
void asm1_cholmod_free(struct box *box);

// ASM1: GPU BLAS interface.
void asm1_gpu_blas_setup(struct csr *A, unsigned null_space, struct box *box,
                         occa::kernel &gatherRHSKernel);
void asm1_gpu_blas_solve(void *x, struct box *box, const void *r);
void asm1_gpu_blas_solve(float *x, struct box *box, occa::memory &o_r);
void asm1_gpu_blas_free(struct box *box);

#endif
