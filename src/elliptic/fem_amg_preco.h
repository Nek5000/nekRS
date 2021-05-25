#ifndef fem_amg_preco_hpp_
#define fem_amg_preco_hpp_

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SEMFEMData{
  long long* Ai;
  long long* Aj;
  double* Av;
  int nnz;
  long long rowStart;
  long long rowEnd;
  long long* dofMap;
} SEMFEMData;

SEMFEMData* fem_amg_setup(const int N_, const int n_elem_, 
                   double *x_m_, double *y_m_, double *z_m_, 
                   double *pmask_,
                   MPI_Comm comm,
                   long long int* gatherGlobalNodes
                   );

#ifdef __cplusplus
}
#endif

#endif
