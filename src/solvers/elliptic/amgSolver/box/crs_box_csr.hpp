#if !defined(_CRS_BOX_CSR_HPP_)
#define _CRS_BOX_CSR_HPP_

#include "crs_box.hpp"

struct csr {
  uint base, nr, *offs, *cols;
  double *vals;
};

struct csr *csr_setup(const unsigned nz, const unsigned *const ia,
                      const unsigned *const ja, const double *const va,
                      const sint *const u2c, const double tol, buffer *bfr);

void csr_free(struct csr *A);

#endif
