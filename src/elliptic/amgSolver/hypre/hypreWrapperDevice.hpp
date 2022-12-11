#ifndef HYPRE_WRAPPER_DEVICE_H
#define HYPRE_WRAPPER_DEVICE_H

#include <mpi.h>
#include "occa.hpp"

namespace hypreWrapperDevice {

int enabled();
void finalize();

constexpr int NPARAM = 11;

#ifndef HYPRE_HEADER
// has to match HYPRE config
typedef int HYPRE_Int;
typedef long long int HYPRE_BigInt;
typedef float HYPRE_Real;
typedef void HYPRE_IJMatrix;
typedef void HYPRE_IJVector;
typedef void HYPRE_Solver;
#endif

class boomerAMG_t {

public:
  ~boomerAMG_t();

  boomerAMG_t(int nrows,
             int nz,
             const long long int *Ai,
             const long long int *Aj,
             const double *Av,
             int null_space,
             MPI_Comm ce,
             occa::device device,
             int useFP32,
             const double *param,
             int verbose);

  void solve(const occa::memory &o_b, const occa::memory &o_x);

private:
  MPI_Comm comm;
  int rank;
  occa::device device;
  int nRows;
  double params[NPARAM];
  HYPRE_Solver *solver;
  HYPRE_IJMatrix *A;
  HYPRE_IJVector *b;
  HYPRE_IJVector *x;
};

} // namespace hypreWrapperDevice

#endif
