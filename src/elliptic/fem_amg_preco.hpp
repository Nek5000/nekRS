#ifndef fem_amg_preco_hpp_
#define fem_amg_preco_hpp_
#include <utility>
#include <map>
#include "gslib.h"

struct SEMFEMData{
  long long* Ai;
  long long* Aj;
  double* Av;
  long long nnz;
  long long rowStart;
  long long rowEnd;
  long long* dofMap;
  ~SEMFEMData(){
    free(Ai);
    free(Aj);
    free(Av);
    free(dofMap);
  }
};

SEMFEMData* fem_amg_setup(const sint *n_x_, const sint *n_y_, const sint *n_z_, 
                   const sint *n_elem_,
                   double *x_m_, double *y_m_, double *z_m_, 
                   double *pmask_,
                   MPI_Comm comm,
                   long long int* gatherGlobalNodes
                   );

#endif