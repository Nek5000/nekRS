#ifndef _PARRSB_IMPL_H_
#define _PARRSB_IMPL_H_

#include "parRSB.h"
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef scalar
#undef scalar
#endif
#define scalar double

#ifdef SCALAR_MAX
#undef SCALAR_MAX
#endif
#define SCALAR_MAX DBL_MAX

#ifdef gs_scalar
#undef gs_scalar
#endif
#define gs_scalar gs_double

#define SCALAR_TOL 1e-12

#define MAXDIM 3 // Maximum dimension of the mesh
#define MAXNV 8  // Maximum number of vertices per element

//------------------------------------------------------------------------------
// RCB / RIB
// `struct rcb_element` is used for RCB and RIB partitioning.
// `struct rsb_element` should be a superset of `struct rcb_element`
struct rcb_element {
  uint proc, origin, seq;
  ulong globalId;
  scalar coord[MAXDIM], fiedler;
};

int rcb(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);
int rib(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);

//------------------------------------------------------------------------------
// RSB
//
struct rsb_element {
  uint proc, origin, seq;
  ulong globalId;
  scalar coord[MAXDIM], fiedler;
  slong vertices[MAXNV];
};

//------------------------------------------------------------------------------
// Find number of components
//
uint get_components(sint *component, struct array *elems, unsigned nv,
                    struct comm *c, buffer *buf, int verbose);
uint get_components_v2(sint *component, struct array *elems, unsigned nv,
                       const struct comm *ci, buffer *bfr, int verbose);

//------------------------------------------------------------------------------
// Laplacian
//
#define GS 1
#define CSR 2
#define CSC 4

struct laplacian;
struct laplacian *laplacian_init(struct rsb_element *elems, uint nel, int nv,
                                 int type, struct comm *c, buffer *buf);
int laplacian(scalar *v, struct laplacian *l, scalar *u, buffer *buf);
void laplacian_free(struct laplacian *l);

//------------------------------------------------------------------------------
// Misc
//
int log2ll(long long n);
void parrsb_barrier(struct comm *c);

#endif
