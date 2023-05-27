#ifndef _PARRSB_ILU_H_
#define _PARRSB_ILU_H_

#include "mat.h"

typedef struct {
  // ILU type: ILU(0), ILUC, etc.
  int type;
  // Verbose level: 0, 1, etc.
  int verbose;
  // Use pivoting or not: 0 or 1
  int pivot;
  // 1st dropping rule: An entry a_ij is dropped abs(a_ij) < tol
  scalar tol;
  // 2nd dropping rule: Entries are dropped so that total nnz per row/col < p
  unsigned int nnz_per_row;
} ilu_options;

struct ilu;
struct ilu *ilu_setup(const uint n, const int nv, const long long *vtx,
                      const ilu_options *options, MPI_Comm comm);
void ilu_free(struct ilu *ilu);

#endif
