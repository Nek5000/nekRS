#if !defined(_CRS_BOX_HPP_)
#define _CRS_BOX_HPP_

// Inclusion of occa.hpp here is a workaround for a weird issue. Fix it before
// releasing rhe code.
#include "occa.hpp"
#define OMPI_SKIP_MPICXX 1
#include "gslib.h"

struct xxt;
struct xxt *crs_xxt_setup(uint n, const ulong *id, uint nz, const uint *Ai,
                          const uint *Aj, const double *A, uint null_space,
                          const struct comm *comm, gs_dom dom);
void crs_xxt_solve(void *x, struct xxt *data, const void *b);
void crs_xxt_stats(struct xxt *data);
void crs_xxt_times(double *cholesky, double *local, double *xxt, double *qqt);
void crs_xxt_free(struct xxt *data);

void crs_xxt_setup_inter_comm(uint n, const ulong *id, uint nz, const uint *Ai,
                              const uint *Aj, const double *A, uint null_space,
                              const MPI_Comm *inter_comm, gs_dom dom,
                              struct comm *local_comm);
void crs_xxt_solve_inter_comm(double *x, const double *rhs);
void crs_xxt_finalize_inter_comm(void);

struct box;
struct box *crs_box_setup(uint n, const ulong *id, uint nnz, const uint *Ai,
                          const uint *Aj, const double *A, uint null,
                          const struct comm *comm, const MPI_Comm *inter_comm,
                          gs_dom dom);
void crs_box_solve(void *x, struct box *data, const void *b);
void crs_box_solve_go_gs(occa::memory &o_x, struct box *data,
                         occa::memory &o_rhs);
void crs_box_solve2(occa::memory &o_x, struct box *data, occa::memory &o_rhs);
void crs_box_free(struct box *data);

#endif
