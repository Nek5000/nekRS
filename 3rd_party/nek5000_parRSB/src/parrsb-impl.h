#ifndef _PARRSB_IMPL_H_
#define _PARRSB_IMPL_H_

#define _POSIX_C_SOURCE 200809L

#include <float.h>
#include <stdlib.h>

#include "parRSB.h"

#ifdef scalar
#undef scalar
#endif
#define scalar double

#ifdef SCALAR_MAX
#undef SCALAR_MAX
#endif
#define SCALAR_MAX DBL_MAX
#define SCALAR_TOL 1e-12

#define MAXDIM 3 // Maximum dimension of the mesh.
#define MAXNV 8  // Maximum number of vertices per element.

//------------------------------------------------------------------------------
// RCB / RIB.
// `struct rcb_element` is used for RCB and RIB partitioning.
struct rcb_element {
  uint proc, origin;
  ulong globalId;
  scalar coord[MAXDIM], fiedler;
};

int rcb(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);
int rib(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);

//------------------------------------------------------------------------------
// RSB.
// `struct rsb_element` = `struct rcb_element` + vertices. Order is important.
struct rsb_element {
  uint proc, origin;
  ulong globalId;
  scalar coord[MAXDIM], fiedler;
  slong vertices[MAXNV];
};

void rsb(struct array *elements, int nv, const parrsb_options *const options,
         const struct comm *comms, buffer *bfr);

//------------------------------------------------------------------------------
// Find number of components.
//
uint get_components(sint *component, struct array *elems, unsigned nv,
                    struct comm *c, buffer *buf, int verbose);
uint get_components_v2(sint *component, struct array *elems, unsigned nv,
                       const struct comm *ci, buffer *bfr, int verbose);

//------------------------------------------------------------------------------
// Dump partition statistics.
//
void parrsb_dump_stats_start(const uint nv_);

void parrsb_dump_stats(const struct comm *const gc, const struct comm *const lc,
                       const struct array *const elems, buffer *bfr);

void parrsb_dump_stats_end(const struct comm *const gc, const char *prefix);

uint parrsb_get_neighbors(const struct array *const elems, const unsigned nv,
                          const struct comm *const gc,
                          const struct comm *const lc, buffer *bfr);

//------------------------------------------------------------------------------
// Laplacian.
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
// Misc.
//
int log2ll(long long n);

void parrsb_barrier(struct comm *c);

void parrsb_print(const struct comm *c, int verbose, const char *fmt, ...);
#endif
