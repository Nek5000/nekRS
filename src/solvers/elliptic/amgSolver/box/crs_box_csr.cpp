#include "crs_box_csr.hpp"
#include <stdlib.h>

struct csr *csr_setup(const unsigned nz, const unsigned *const ia,
                      const unsigned *const ja, const double *const va,
                      const sint *const u2c, const double tol, buffer *bfr) {
  // Assemble the local matrix now.
  struct mij_t {
    uint r, c;
    double v;
  };

  struct array mijs;
  array_init(struct mij_t, &mijs, nz);

  struct mij_t mij;
  for (uint z = 0; z < nz; z++) {
    sint i = u2c[ia[z]], j = u2c[ja[z]];
    if (i < 0 || j < 0 || fabs(va[z]) < tol)
      continue;
    mij.r = i, mij.c = j, mij.v = va[z];
    array_cat(struct mij_t, &mijs, &mij, 1);
  }

  struct array aijs;
  array_init(struct mij_t, &aijs, nz);

  sarray_sort_2(struct mij_t, mijs.ptr, mijs.n, r, 0, c, 0, bfr);
  if (mijs.n > 0) {
    struct mij_t *pm = (struct mij_t *)mijs.ptr;
    uint s = 0;
    for (uint i = 1; i < mijs.n; i++) {
      if ((pm[i].r != pm[s].r) || (pm[i].c != pm[s].c)) {
        array_cat(struct mij_t, &aijs, &pm[s], 1);
        s = i;
      } else {
        pm[s].v += pm[i].v;
      }
    }
    array_cat(struct mij_t, &aijs, &pm[s], 1);
  }
  array_free(&mijs);

  // Convert the assembled matrix to a CSR matrix.
  struct csr *A = tcalloc(struct csr, 1);
  A->base = 0;
  A->nr = 0;
  struct mij_t *pa = (struct mij_t *)aijs.ptr;
  if (aijs.n > 0)
    A->nr = pa[aijs.n - 1].r + 1;

  // Allocate the csr data structures
  A->offs = tcalloc(uint, A->nr + 1);
  A->cols = tcalloc(uint, aijs.n);
  A->vals = tcalloc(double, aijs.n);
  A->offs[0] = 0;
  for (uint i = 0, r = 0; r < A->nr; r++) {
    uint j = i;
    while (j < aijs.n && pa[j].r == r) {
      A->cols[j] = pa[j].c, A->vals[j] = pa[j].v;
      j++;
    }
    A->offs[r + 1] = A->offs[r] + j - i;
    i = j;
  }
  array_free(&aijs);

  return A;
}

void csr_free(struct csr *A) {
  if (!A)
    return;
  free(A->offs), free(A->cols), free(A->vals);
  free(A);
}
