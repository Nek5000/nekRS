#ifndef AMGX_H
#define AMGX_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

int AMGXsetup(const int nLocalRows, const int nnz,
              const long long *rows, const long long *cols, const double *values, /* COO */ 
              const int null_space, const MPI_Comm comm, int deviceID,
              int useFP32, int MPIDIRECT, const char* cfgFile);
int AMGXsolve(void *rhs, void *x);
void AMGXfree();
int AMGXenabled();

#ifdef __cplusplus
}
#endif

#endif
