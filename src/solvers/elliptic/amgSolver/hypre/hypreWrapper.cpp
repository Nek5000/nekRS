#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <vector>
#include <map>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

#include "hypreWrapper.hpp"

static int HYPREinit = 0;

#define DEBUG 0

namespace hypreWrapper {

__attribute__((visibility("default")))
boomerAMG_t::boomerAMG_t(int _nRows,
                   int nz,
                   const long long int *Ai,
                   const long long int *Aj,
                   const double *Av,
                   int null_space,
                   MPI_Comm ce,
                   int _Nthreads,
                   int useFP32,
                   const double *param,
                   int verbose)
{
  nRows = _nRows;
  Nthreads = _Nthreads;

  MPI_Comm_dup(ce, &comm);
  MPI_Comm_rank(comm, &rank);

  if (sizeof(HYPRE_Real) != ((useFP32) ? sizeof(float) : sizeof(double))) {
    if (rank == 0) {
      printf("hypreWrapperDevice: HYPRE floating point precision does not match!\n");
      MPI_Abort(comm, 1);
    }
  }

  if ((int)param[0]) {
    for (int i = 0; i < NPARAM; i++)
      params[i] = param[i + 1];
  }
  else {
    params[0] = 8;    /* coarsening */
    params[1] = 6;    /* interpolation */
    params[2] = 1;    /* number of cycles */
    params[3] = 16;   /* smoother for crs level */
    params[4] = 3;    /* sweeps */
    params[5] = 16;   /* smoother */
    params[6] = 1;    /* sweeps   */
    params[7] = 0.25; /* threshold */
    params[8] = 0.0;  /* non galerkin tolerance */
    params[9] = 0.0;  /* agressive coarsening */
    params[10] = 2;  /* chebyRelaxOrder */
  }

  // Setup matrix
  long long rowStart = nRows;
  MPI_Scan(MPI_IN_PLACE, &rowStart, 1, MPI_LONG_LONG, MPI_SUM, comm);
  rowStart -= nRows;

  const HYPRE_BigInt ilower = (HYPRE_BigInt)rowStart;
  const HYPRE_BigInt iupper = ilower + (HYPRE_BigInt)nRows - 1;

  if(!HYPREinit) HYPRE_Init();

  A = new HYPRE_IJMatrix();
  HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, A);
  HYPRE_IJMatrixSetObjectType(*A, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(*A);

  std::map<HYPRE_BigInt, std::vector<std::pair<HYPRE_BigInt, HYPRE_Real>>> rowToColAndVal;
  for (int i = 0; i < nz; i++) {
    HYPRE_BigInt mati = (HYPRE_BigInt)(Ai[i]);
    HYPRE_BigInt matj = (HYPRE_BigInt)(Aj[i]);
    HYPRE_Real matv = (HYPRE_Real)Av[i];
    rowToColAndVal[mati].emplace_back(std::make_pair(matj, matv));
  }

  const HYPRE_Int rowsToSet = rowToColAndVal.size();
  std::vector<HYPRE_Int> ncols(rowsToSet);
  std::vector<HYPRE_BigInt> rows(rowsToSet);
  std::vector<HYPRE_BigInt> cols(nz);
  std::vector<HYPRE_Real> vals(nz);

  unsigned rowCtr = 0;
  unsigned colCtr = 0;
  for (auto &&rowAndColValPair : rowToColAndVal) {
    const auto &row = rowAndColValPair.first;
    const auto &colAndValues = rowAndColValPair.second;

    rows[rowCtr] = row;
    ncols[rowCtr] = colAndValues.size();

    for (auto &&colAndValue : colAndValues) {
      const auto &col = colAndValue.first;
      const auto &val = colAndValue.second;
      cols[colCtr] = col;
      vals[colCtr] = val;
      ++colCtr;
    }
    ++rowCtr;
  }

  HYPRE_IJMatrixSetValues(*A, rowsToSet, ncols.data(), rows.data(), cols.data(), vals.data());

  HYPRE_IJMatrixAssemble(*A);

  if(DEBUG)
    HYPRE_IJMatrixPrint(*A, "matrix.dat");

  // Setup solver
  solver = new HYPRE_Solver(); 
  HYPRE_BoomerAMGCreate(solver);

  HYPRE_BoomerAMGSetCoarsenType(*solver, params[0]);
  HYPRE_BoomerAMGSetInterpType(*solver, params[1]);

  HYPRE_BoomerAMGSetChebyOrder(*solver, params[10]);
  // HYPRE_BoomerAMGSetChebyFraction(*solver, 0.2);

  if (params[5] > 0) {
    HYPRE_BoomerAMGSetCycleRelaxType(*solver, params[5], 1);
    HYPRE_BoomerAMGSetCycleRelaxType(*solver, params[5], 2);
  }
  HYPRE_BoomerAMGSetCycleRelaxType(*solver, 9, 3);

  HYPRE_BoomerAMGSetCycleNumSweeps(*solver, params[6], 1);
  HYPRE_BoomerAMGSetCycleNumSweeps(*solver, params[6], 2);
  HYPRE_BoomerAMGSetCycleNumSweeps(*solver, 1, 3);

  if (null_space) {
    HYPRE_BoomerAMGSetMinCoarseSize(*solver, 2);
    HYPRE_BoomerAMGSetCycleRelaxType(*solver, params[3], 3);
    HYPRE_BoomerAMGSetCycleNumSweeps(*solver, params[4], 3);
  }

  HYPRE_BoomerAMGSetStrongThreshold(*solver, params[7]);

  if (params[8] > 1e-3) {
    HYPRE_BoomerAMGSetNonGalerkinTol(*solver, params[8]);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(*solver, 0.0, 0);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(*solver, 0.01, 1);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(*solver, 0.05, 2);
  }

  HYPRE_BoomerAMGSetAggNumLevels(*solver, params[9]);

  HYPRE_BoomerAMGSetMaxIter(*solver, params[2]); // number of V-cycles
  HYPRE_BoomerAMGSetTol(*solver, 0);

  if (DEBUG)
    HYPRE_BoomerAMGSetPrintLevel(*solver, 3);
  else
    HYPRE_BoomerAMGSetPrintLevel(*solver, 1);

  // Create and initialize rhs and solution vectors
  b = new HYPRE_IJVector();
  HYPRE_IJVectorCreate(comm, ilower, iupper, b);
  HYPRE_IJVectorSetObjectType(*b, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(*b);
  HYPRE_IJVectorAssemble(*b);

  x = new HYPRE_IJVector();
  HYPRE_IJVectorCreate(comm, ilower, iupper, x);
  HYPRE_IJVectorSetObjectType(*x, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(*x);
  HYPRE_IJVectorAssemble(*x);

  // Perform AMG setup
  HYPRE_ParVector par_b;
  HYPRE_ParVector par_x;
  HYPRE_IJVectorGetObject(*b, (void **)&par_b);
  HYPRE_IJVectorGetObject(*x, (void **)&par_x);
  HYPRE_ParCSRMatrix par_A;
  HYPRE_IJMatrixGetObject(*A, (void **)&par_A);

#ifdef _OPENMP
  const int NthreadsSave = omp_get_max_threads();
  omp_set_num_threads(Nthreads);
  omp_set_num_threads(Nthreads);
#endif

  const int err = HYPRE_BoomerAMGSetup(*solver, par_A, par_b, par_x);
  if (err > 0) {
    if (rank == 0)
      printf("HYPRE_BoomerAMGSetup failed with %d!\n", err);
    MPI_Abort(comm, err);
  }

#ifdef _OPENMP
  omp_set_num_threads(NthreadsSave);
#endif
  HYPREinit = 1;
}

void __attribute__((visibility("default"))) boomerAMG_t::solve(void *bin, void *xin)
{
  HYPRE_IJVectorUpdateValues(*x, nRows, NULL, (HYPRE_Real *)xin, 1);
  HYPRE_IJVectorUpdateValues(*b, nRows, NULL, (HYPRE_Real *)bin, 1);

  HYPRE_ParVector par_x;
  HYPRE_ParVector par_b;
  HYPRE_ParCSRMatrix par_A;

  HYPRE_IJVectorGetObject(*x, (void **)&par_x);
  HYPRE_IJVectorGetObject(*b, (void **)&par_b);
  HYPRE_IJMatrixGetObject(*A, (void **)&par_A);

#ifdef _OPENMP
  const int _Nthreads = omp_get_max_threads();
  omp_set_num_threads(Nthreads);
#endif

#if 0
  HYPRE_IJVectorPrint(*b, "b.dat");
  HYPRE_IJVectorPrint(*x, "x.dat");
#endif

  const int err = HYPRE_BoomerAMGSolve(*solver, par_A, par_b, par_x);
  if (err > 0) {
    if (rank == 0)
      printf("HYPRE_BoomerAMGSolve failed with %d!\n", err);
    MPI_Abort(comm, err);
  }
#ifdef _OPENMP
  omp_set_num_threads(_Nthreads);
#endif

  HYPRE_IJVectorGetValues(*x, nRows, NULL, (HYPRE_Real *)xin);
}

__attribute__((visibility("default"))) boomerAMG_t::~boomerAMG_t()
{
  if(solver) HYPRE_BoomerAMGDestroy(*solver);
  if(A) HYPRE_IJMatrixDestroy(*A);
  if(x) HYPRE_IJVectorDestroy(*x);
  if(b) HYPRE_IJVectorDestroy(*b);
}

void __attribute__((visibility("default"))) finalize()
{
  if(!HYPREinit) return;

  HYPRE_Finalize();
}

__attribute__((visibility("default"))) IJ_t::IJ_t() {}

__attribute__((visibility("default"))) IJ_t::~IJ_t()
{
  if(A) HYPRE_IJMatrixDestroy(*A);
}

int __attribute__((visibility("default")))
IJ_t::MatrixGetRowCounts(HYPRE_Int nrows, HYPRE_BigInt *rows, HYPRE_Int *ncols)
{
  return HYPRE_IJMatrixGetRowCounts(*(HYPRE_IJMatrix *)A, nrows, rows, ncols);
}

int __attribute__((visibility("default")))
IJ_t::MatrixGetValues(HYPRE_Int nrows,
                               HYPRE_Int *ncols,
                               HYPRE_BigInt *rows,
                               HYPRE_BigInt *cols,
                               HYPRE_Complex *values)
{
  return HYPRE_IJMatrixGetValues(*(HYPRE_IJMatrix *)A, nrows, ncols, rows, cols, values);
}

int __attribute__((visibility("default"))) 
IJ_t::MatrixAddToValues(HYPRE_Int nrows,
                                 HYPRE_Int *ncols,
                                 const HYPRE_BigInt *rows,
                                 const HYPRE_BigInt *cols,
                                 const HYPRE_Complex *values)
{
  return HYPRE_IJMatrixAddToValues(*(HYPRE_IJMatrix *)A, nrows, ncols, rows, cols, values);
}

int __attribute__((visibility("default"))) 
IJ_t::MatrixCreate(MPI_Comm comm,
                            HYPRE_BigInt ilower,
                            HYPRE_BigInt iupper,
                            HYPRE_BigInt jlower,
                            HYPRE_BigInt jupper)
{
  A = new HYPRE_IJMatrix(); 
  return HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, A);
}

int __attribute__((visibility("default"))) IJ_t::MatrixSetObjectType()
{
  return HYPRE_IJMatrixSetObjectType(*(HYPRE_IJMatrix *)A, HYPRE_PARCSR);
}

int __attribute__((visibility("default"))) IJ_t::MatrixInitialize()
{
  return HYPRE_IJMatrixInitialize(*(HYPRE_IJMatrix *)A);
}

int __attribute__((visibility("default"))) IJ_t::MatrixAssemble()
{
  return HYPRE_IJMatrixAssemble(*(HYPRE_IJMatrix *)A);
}

} // namespace hypreWrapper
