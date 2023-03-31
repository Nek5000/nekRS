#ifndef _PARRSB_COARSE_IMPL_H_
#define _PARRSB_COARSE_IMPL_H_

#include "coarse.h"

uint unique_ids(sint *perm, ulong *uid, uint n, const ulong *ids, buffer *bfr);

struct coarse {
  unsigned type;       // type = schur-2-lvl, schur-3-lvl
  unsigned null_space; // Is there a null space or not
  uint un;             // User vector size
  uint cn;   // Compressed (ignoring duplicates and zero global ids) vector size
  uint an;   // Assembled size -- this is the local size of the assmebled coarse
             // matrix
  sint *u2c; // Mapping from user vector to compress vector
  struct gs_data *c2a; // Mapping from compressed vector to assmbled vector
  buffer bfr;

  ulong s[3], ng[3];
  uint n[3];
  struct comm c;
  void *solver;
};

int schur_setup(struct coarse *crs, struct array *eij, struct crystal *cr,
                buffer *bfr);
int schur_solve(scalar *x, struct coarse *crs, scalar *b, scalar tol,
                buffer *bfr);
int schur_free(struct coarse *crs);

#endif
