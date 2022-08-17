#ifndef HYPRE_WRAPPER_H
#define HYPRE_WRAPPER_H

#include <mpi.h>
#include "occa.hpp"

#define BOOMERAMG_NPARAM 10

namespace hypreWrapper
{

// has to match HYPRE config
typedef int Int;
typedef long long int BigInt;
typedef float Real;
typedef void IJMatrix;

typedef Int HYPRE_Int;
typedef BigInt HYPRE_BigInt;
typedef IJMatrix HYPRE_IJMatrix;
typedef Real HYPRE_Real;

int BoomerAMGSetup(int nrows,
                   int nz, const long long int *Ai, const long long int *Aj, const double *Av,
                   int null_space, MPI_Comm ce, int Nthreads,
                   int useFP32, const double *param, int verbose);

int BoomerAMGSolve(void *b, void *x);

void Free();

void IJMatrixGetRowCounts
(
  HYPRE_IJMatrix  *matrix,
  HYPRE_Int       nrows,
  HYPRE_BigInt   *rows,
  HYPRE_Int      *ncols
);

void IJMatrixGetValues
(
  HYPRE_IJMatrix  *matrix,
  HYPRE_Int       nrows,
  HYPRE_Int      *ncols,
  HYPRE_BigInt   *rows,
  HYPRE_BigInt   *cols,
  HYPRE_Real  *values
);

int IJMatrixDestroy
(
  HYPRE_IJMatrix *matrix
);

int IJMatrixAddToValues
(
  HYPRE_IJMatrix      *matrix,
  HYPRE_Int            nrows,
  HYPRE_Int           *ncols,
  const HYPRE_BigInt  *rows,
  const HYPRE_BigInt  *cols,
  const HYPRE_Real *values
);

int IJMatrixCreate
(
  MPI_Comm        comm,
  HYPRE_BigInt    ilower,
  HYPRE_BigInt    iupper,
  HYPRE_BigInt    jlower,
  HYPRE_BigInt    jupper,
  HYPRE_IJMatrix *matrix
);

int IJMatrixSetObjectType
(
  HYPRE_IJMatrix *matrix
);

int IJMatrixInitialize
(
  HYPRE_IJMatrix *matrix
);

int IJMatrixAssemble
(
  HYPRE_IJMatrix *matrix
);

} // namespace

namespace hypreWrapperDevice
{

// has to match HYPRE config
typedef int HYPRE_Int;
typedef long long int HYPRE_BigInt;
typedef float HYPRE_Real;
typedef void HYPRE_IJMatrix;

int BoomerAMGSetup(int nrows, int nz,
                   const long long int * Ai, const long long int * Aj, const double * Av,
                   int null_space, MPI_Comm ce, occa::device device,
                   int useFP32, const double *param, int verbose);

int BoomerAMGSolve(const occa::memory& o_b, const occa::memory& o_x);

void Free();

} // namespace

#endif
