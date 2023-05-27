#ifndef _PARRSB_MULTIGRID_H_
#define _PARRSB_MULTIGRID_H_

#include "mat.h"

struct mg;
struct mg *mg_setup(const struct par_mat *M, const int factor, const int sagg,
                    struct crystal *cr, buffer *bfr);
void mg_vcycle(scalar *u, scalar *rhs, struct mg *d, struct comm *c,
               buffer *bfr);
void mg_free(struct mg *d);

#endif
