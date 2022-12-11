#ifndef _PARRSB_COARSE_H_
#define _PARRSB_COARSE_H_

#include "gslib.h"
#include "mat.h"

struct coarse;

// API for the Laplacian (which involves solving for the dual graph)
struct coarse *coarse_setup(unsigned nelt, unsigned nv, const long long *vtx,
                            const scalar *coord, unsigned null_space,
                            unsigned type, struct comm *c);
void coarse_solve(scalar *x, struct coarse *crs, scalar *b, scalar tol);
void coarse_free(struct coarse *crs);

// Alternative API for a general matrix
#define crs_parrsb_setup PREFIXED_NAME(crs_parrsb_setup)
#define crs_parrsb_solve PREFIXED_NAME(crs_parrsb_solve)
#define crs_parrsb_free PREFIXED_NAME(crs_parrsb_free)

struct coarse *crs_parrsb_setup(uint n, const ulong *id, uint nz,
                                const uint *Ai, const uint *Aj, const scalar *A,
                                unsigned null_space, unsigned type,
                                const struct comm *comm);
void crs_parrsb_solve(scalar *x, struct coarse *crs, scalar *b, scalar tol);
void crs_parrsb_free(struct coarse *crs);

#endif
