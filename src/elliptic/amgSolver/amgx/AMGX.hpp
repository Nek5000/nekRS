#ifndef AMGX_H
#define AMGX_H

#include <mpi.h>

#ifdef ENABLE_AMGX
#include <amgx_c.h>
#else
typedef void* AMGX_vector_handle;
typedef void* AMGX_matrix_handle;
typedef void* AMGX_solver_handle;
typedef void* AMGX_config_handle;
typedef void* AMGX_resources_handle;
#endif

int AMGXenabled();
void AMGXfinalize();

class AMGX_t
{
public:

  ~AMGX_t();

  AMGX_t(const int nLocalRows, const int nnz,
         const long long *rows, const long long *cols, const double *values, /* COO */ 
         const int null_space, const MPI_Comm comm, int deviceID,
         int useFP32, int MPIDIRECT, const char* cfgFile);
  int solve(void *rhs, void *x);

private:
  MPI_Comm comm;
  int nLocalRows;
  AMGX_vector_handle AmgXP;
  AMGX_vector_handle AmgXRHS;
  AMGX_matrix_handle AmgXA;
  AMGX_solver_handle solver;
  AMGX_config_handle cfg;
  AMGX_resources_handle rsrc;
};

#endif
