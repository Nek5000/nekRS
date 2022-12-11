#ifndef HYPRE_WRAPPER_H
#define HYPRE_WRAPPER_H

#include <mpi.h>

namespace hypreWrapper {

constexpr int NPARAM = 11;

#ifndef HYPRE_HEADER
// has to match HYPRE config
typedef int Int;
typedef long long int BigInt;
typedef float Real;

typedef void IJMatrix;
typedef void IJVector;
typedef void Solver;

typedef Int HYPRE_Int;
typedef BigInt HYPRE_BigInt;
typedef IJMatrix HYPRE_IJMatrix;
typedef IJVector HYPRE_IJVector;
typedef Solver HYPRE_Solver;
typedef Real HYPRE_Real;
#endif

void finalize();

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
             int Nthreads,
             int useFP32,
             const double *param,
             int verbose);

  void solve(void *b, void *x);

private:
  MPI_Comm comm;
  int nRows;
  int Nthreads;
  int rank;
  double params[NPARAM];
  HYPRE_Solver *solver;
  HYPRE_IJMatrix *A;
  HYPRE_IJVector *b;
  HYPRE_IJVector *x;
};

class IJ_t {

public:
  IJ_t();

  ~IJ_t();

  int MatrixGetRowCounts(HYPRE_Int nrows, HYPRE_BigInt *rows, HYPRE_Int *ncols);

  int MatrixGetValues(HYPRE_Int nrows,
                        HYPRE_Int *ncols,
                        HYPRE_BigInt *rows,
                        HYPRE_BigInt *cols,
                        HYPRE_Real *values);

  int MatrixAddToValues(HYPRE_Int nrows,
                          HYPRE_Int *ncols,
                          const HYPRE_BigInt *rows,
                          const HYPRE_BigInt *cols,
                          const HYPRE_Real *values);

  int MatrixCreate(MPI_Comm comm,
                     HYPRE_BigInt ilower,
                     HYPRE_BigInt iupper,
                     HYPRE_BigInt jlower,
                     HYPRE_BigInt jupper);

  int MatrixSetObjectType();

  int MatrixInitialize();

  int MatrixAssemble();

private:
  HYPRE_IJMatrix *A;
};



} // namespace hypreWrapper

#endif
