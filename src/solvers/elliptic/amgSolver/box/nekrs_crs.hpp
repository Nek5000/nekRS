#if !defined(_NEKRS_CRS_HPP_)
#define _NEKRS_CRS_HPP_

#include "elliptic.h"
#include "gslib.h"
#include "platform.hpp"
#include <string>

#define JL_XXT 1
#define JL_BOX 2

void jl_setup_aux(uint *ntot, ulong **gids, uint *nnz, uint **ia, uint **ja,
                  double **a, elliptic_t *elliptic, elliptic_t *ellipticf);

void jl_setup(uint type, uint n, const ulong *id, uint nnz, const uint *Ai,
              const uint *Aj, const double *A, uint null, gs_dom dom,
              MPI_Comm comm, const MPI_Comm *inter_comm);

void jl_solve(occa::memory o_x, occa::memory o_rhs);

void jl_solve2(occa::memory o_x, occa::memory o_rhs);

void jl_free();

#endif
