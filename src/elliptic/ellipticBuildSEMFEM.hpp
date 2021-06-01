#ifndef fem_amg_preco_hpp_
#define fem_amg_preco_hpp_

#include <platform.hpp>
#include <mpi.h>

struct SEMFEMData{
  long long* Ai;
  long long* Aj;
  double* Av;
  int nnz;
  long long rowStart;
  long long rowEnd;
  long long* dofMap;
};

SEMFEMData* ellipticBuildSEMFEM(const int N_, const int n_elem_, 
                   occa::memory _o_x, occa::memory _o_y, occa::memory _o_z,
                   double *pmask_,
                   MPI_Comm comm,
                   long long int* gatherGlobalNodes);


#endif
