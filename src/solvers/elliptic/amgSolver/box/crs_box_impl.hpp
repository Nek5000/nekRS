#if !defined(_CRS_BOX_IMPL_HPP_)
#define _CRS_BOX_IMPL_HPP_

#include "crs_box.hpp"

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

#define BOX_GPU_BLAS 4

// Structure to keep track of data structures in the box solver.
struct box {
  // User input.
  uint un, ncr;
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

// CSR matrix used to setup the local solver.
struct csr {
  uint base, nr, *offs, *cols;
  double *vals;
};

struct csr *csr_setup(const unsigned nz, const unsigned *const ia,
                      const unsigned *const ja, const double *const va,
                      const sint *const u2c, const double tol, buffer *bfr);

void csr_free(struct csr *A);

// Local solver.
void asm1_setup(struct csr *A, unsigned null_space, struct box *box);
void asm1_solve(void *x, struct box *box, occa::memory &o_r);
void asm1_solve(void *x, struct box *box, const void *r);
void asm1_free(struct box *box);

// Timer routines.
enum BOX_METRIC {
  COPY_RHS = 0,
  CRS_DSAVG1 = 1,
  ASM1 = 2,
  CRS_DSAVG2 = 3,
  MULT_RHS_UPDATE = 4,
  COPY_TO_NEK5000 = 5,
  MAP_VTX_TO_BOX = 6,
  ASM2 = 7,
  MAP_BOX_TO_VTX = 8,
  COPY_FROM_NEK5000 = 9,
  CRS_DSAVG3 = 10,
  COPY_SOLUTION = 11,
  NONE = 100
};

void timer_init();
void timer_tic(const struct comm *c);
void timer_toc(BOX_METRIC m);
void timer_dump(struct comm *c, unsigned interval);
void timer_print(struct comm *c, unsigned interval);

#endif
