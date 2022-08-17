#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <vector>
#include <map>

#include "omp.h"

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

#define BOOMERAMG_NPARAM 10
static double boomerAMGParam[BOOMERAMG_NPARAM];

struct hypre_data {
  MPI_Comm comm;
  HYPRE_Solver solver;
  HYPRE_IJMatrix A;
  HYPRE_IJVector b;
  HYPRE_IJVector x;
  HYPRE_BigInt ilower;
  HYPRE_Real *bb;
  int nRows;
  int Nthreads;
};
static hypre_data *data;

namespace hypreWrapper
{

int __attribute__((visibility("default"))) BoomerAMGSetup(int nrows,
                   int nz, const long long int *Ai, const long long int *Aj, const double *Av,
                   int null_space, MPI_Comm ce, int Nthreads,
                   int useFP32, const double *param, int verbose)
{

  data = (hypre_data*) malloc(sizeof(struct hypre_data));

  data->Nthreads = Nthreads;   

  MPI_Comm comm;
  MPI_Comm_dup(ce, &comm);
  data->comm = comm;
  int rank;
  MPI_Comm_rank(comm,&rank);

  if(sizeof(HYPRE_Real) != ((useFP32) ? sizeof(float) : sizeof(double))) {
    if(rank == 0) printf("hypreWrapperDevice: HYPRE floating point precision does not match!\n");
    MPI_Abort(ce, 1);
  } 

  if ((int) param[0]) {
    for (int i = 0; i < BOOMERAMG_NPARAM; i++)
      boomerAMGParam[i] = param[i+1]; 
  } else {
    boomerAMGParam[0]  = 8;    /* coarsening */
    boomerAMGParam[1]  = 6;    /* interpolation */
    boomerAMGParam[2]  = 1;    /* number of cycles */
    boomerAMGParam[3]  = 16;   /* smoother for crs level */
    boomerAMGParam[4]  = 3;    /* sweeps */
    boomerAMGParam[5]  = 16;   /* smoother */
    boomerAMGParam[6]  = 1;    /* sweeps   */
    boomerAMGParam[7]  = 0.25; /* threshold */
    boomerAMGParam[8]  = 0.0;  /* non galerkin tolerance */
  }

  // Setup matrix
  long long rowStart = nrows;
  MPI_Scan(MPI_IN_PLACE, &rowStart, 1, MPI_LONG_LONG, MPI_SUM, ce);
  rowStart -= nrows;

  data->nRows = nrows;
  HYPRE_BigInt ilower = (HYPRE_BigInt) rowStart;
  data->ilower = ilower;
  HYPRE_BigInt iupper = ilower + (HYPRE_BigInt) nrows - 1; 

  HYPRE_IJMatrixCreate(comm,ilower,iupper,ilower,iupper,&data->A);
  HYPRE_IJMatrix A_ij = data->A;
  HYPRE_IJMatrixSetObjectType(A_ij,HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(A_ij);

  std::map<HYPRE_BigInt, std::vector<std::pair<HYPRE_BigInt, HYPRE_Real>>> rowToColAndVal;
  for(int i=0; i<nz; i++) 
  {
    HYPRE_BigInt mati = (HYPRE_BigInt)(Ai[i]);
    HYPRE_BigInt matj = (HYPRE_BigInt)(Aj[i]);
    HYPRE_Real matv = (HYPRE_Real) Av[i]; 
    rowToColAndVal[mati].emplace_back(std::make_pair(matj, matv));
  }

  const HYPRE_Int rowsToSet = rowToColAndVal.size();
  std::vector<HYPRE_Int> ncols(rowsToSet);
  std::vector<HYPRE_BigInt> rows(rowsToSet);
  std::vector<HYPRE_BigInt> cols(nz);
  std::vector<HYPRE_Real> vals(nz);

  unsigned rowCtr = 0;
  unsigned colCtr = 0;
  for(auto&& rowAndColValPair : rowToColAndVal){
    const auto & row = rowAndColValPair.first;
    const auto & colAndValues = rowAndColValPair.second;

    rows[rowCtr] = row;
    ncols[rowCtr] = colAndValues.size();

    for(auto&& colAndValue : colAndValues){
      const auto & col = colAndValue.first;
      const auto & val = colAndValue.second;
      cols[colCtr] = col;
      vals[colCtr] = val;
      ++colCtr;
    }
    ++rowCtr;
  }

  HYPRE_IJMatrixSetValues(A_ij, rowsToSet, ncols.data(), rows.data(), cols.data(), vals.data());

  HYPRE_IJMatrixAssemble(A_ij);

#if 0
  HYPRE_IJMatrixPrint(A_ij, "matrix.dat");
#endif

  // Setup solver
  HYPRE_BoomerAMGCreate(&data->solver);
  HYPRE_Solver solver = data->solver;

  HYPRE_BoomerAMGSetCoarsenType(solver,boomerAMGParam[0]);
  HYPRE_BoomerAMGSetInterpType(solver,boomerAMGParam[1]);

  //HYPRE_BoomerAMGSetChebyFraction(solver, 0.2); 

  if (boomerAMGParam[5] > 0) {
    HYPRE_BoomerAMGSetCycleRelaxType(solver, boomerAMGParam[5], 1);
    HYPRE_BoomerAMGSetCycleRelaxType(solver, boomerAMGParam[5], 2);
  } 
  HYPRE_BoomerAMGSetCycleRelaxType(solver, 9, 3);

  HYPRE_BoomerAMGSetCycleNumSweeps(solver, boomerAMGParam[6], 1);
  HYPRE_BoomerAMGSetCycleNumSweeps(solver, boomerAMGParam[6], 2);
  HYPRE_BoomerAMGSetCycleNumSweeps(solver, 1, 3);

  if (null_space) {
    HYPRE_BoomerAMGSetMinCoarseSize(solver, 2);
    HYPRE_BoomerAMGSetCycleRelaxType(solver, boomerAMGParam[3], 3);
    HYPRE_BoomerAMGSetCycleNumSweeps(solver, boomerAMGParam[4], 3);
  }

  HYPRE_BoomerAMGSetStrongThreshold(solver,boomerAMGParam[7]);

  if (boomerAMGParam[8] > 1e-3) {
    HYPRE_BoomerAMGSetNonGalerkinTol(solver,boomerAMGParam[8]);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(solver,0.0 , 0);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(solver,0.01, 1);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(solver,0.05, 2);
  }

  HYPRE_BoomerAMGSetAggNumLevels(solver, boomerAMGParam[9]); 

  HYPRE_BoomerAMGSetMaxIter(solver,boomerAMGParam[2]); // number of V-cycles
  HYPRE_BoomerAMGSetTol(solver,0);

  if(verbose)
    HYPRE_BoomerAMGSetPrintLevel(solver,3);
  else
    HYPRE_BoomerAMGSetPrintLevel(solver,1);

  // Create and initialize rhs and solution vectors
  HYPRE_IJVectorCreate(comm,ilower,iupper,&data->b);
  HYPRE_IJVector b = data->b;
  HYPRE_IJVectorSetObjectType(b,HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(b);
  HYPRE_IJVectorAssemble(b);

  HYPRE_IJVectorCreate(comm,ilower,iupper,&data->x);
  HYPRE_IJVector x = data->x;
  HYPRE_IJVectorSetObjectType(x,HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(x);
  HYPRE_IJVectorAssemble(x);

  // Perform AMG setup
  HYPRE_ParVector par_b;
  HYPRE_ParVector par_x;
  HYPRE_IJVectorGetObject(b,(void**) &par_b);
  HYPRE_IJVectorGetObject(x,(void**) &par_x);
  HYPRE_ParCSRMatrix par_A;
  HYPRE_IJMatrixGetObject(data->A,(void**) &par_A);

  const int _Nthreads = omp_get_max_threads();
  omp_set_num_threads(data->Nthreads);
  HYPRE_BoomerAMGSetup(solver,par_A,par_b,par_x);
  omp_set_num_threads(_Nthreads);

  return 0;
}

int __attribute__((visibility("default"))) BoomerAMGSolve(void *b, void *x)
{
  int err;

#if 1
  HYPRE_IJVectorUpdateValues(data->x,data->nRows,NULL,(HYPRE_Real*) x, 1);
#else
  HYPRE_IJVectorSetValues(data->x,data->nRows,NULL,(HYPRE_Real*) o_x.ptr());
  HYPRE_IJVectorAssemble(data->x);
#endif

#if 1
  HYPRE_IJVectorUpdateValues(data->b,data->nRows,NULL,(HYPRE_Real*) b, 1);
#else
  HYPRE_IJVectorSetValues(data->b,data->nRows,NULL,(HYPRE_Real*) o_b.ptr());
  HYPRE_IJVectorAssemble(data->b);
#endif

  HYPRE_ParVector par_x;
  HYPRE_ParVector par_b;
  HYPRE_ParCSRMatrix par_A;

  HYPRE_IJVectorGetObject(data->x,(void **) &par_x);
  HYPRE_IJVectorGetObject(data->b,(void **) &par_b);
  HYPRE_IJMatrixGetObject(data->A,(void **) &par_A);

  const int _Nthreads = omp_get_max_threads();
  omp_set_num_threads(data->Nthreads);

#if 0
  HYPRE_IJVectorPrint(data->b, "b.dat");
  HYPRE_IJVectorPrint(data->x, "x.dat");
#endif

  err = HYPRE_BoomerAMGSolve(data->solver,par_A,par_b,par_x);
  if(err > 0) { 
    int rank;
    MPI_Comm_rank(data->comm,&rank);
    if(rank == 0) printf("HYPRE_BoomerAMGSolve failed!\n");
    return 1;
  }
  omp_set_num_threads(_Nthreads);

  HYPRE_IJVectorGetValues(data->x,data->nRows,NULL,(HYPRE_Real*) x);

  return 0; 
}

void __attribute__((visibility("default"))) Free()
{
  HYPRE_BoomerAMGDestroy(data->solver);
  HYPRE_IJMatrixDestroy(data->A);
  HYPRE_IJVectorDestroy(data->x);
  HYPRE_IJVectorDestroy(data->b);
  HYPRE_Finalize();
  free(data);
}

int __attribute__((visibility("default"))) IJMatrixGetRowCounts
(
  void  *matrix,
  HYPRE_Int       nrows,
  HYPRE_BigInt   *rows,
  HYPRE_Int      *ncols
)
{
  return HYPRE_IJMatrixGetRowCounts(*(HYPRE_IJMatrix*) matrix, nrows, rows, ncols);
}

int __attribute__((visibility("default"))) IJMatrixGetValues
(
  void           *matrix,
  HYPRE_Int       nrows,
  HYPRE_Int      *ncols,
  HYPRE_BigInt   *rows,
  HYPRE_BigInt   *cols,
  HYPRE_Complex  *values
)
{
  return HYPRE_IJMatrixGetValues(*(HYPRE_IJMatrix*) matrix, nrows, ncols, rows, cols, values);
}

int __attribute__((visibility("default"))) IJMatrixDestroy
(
  void *matrix
)
{
  return HYPRE_IJMatrixDestroy(*(HYPRE_IJMatrix*) matrix);
}

int __attribute__((visibility("default"))) IJMatrixAddToValues
(
  void                *matrix,
  HYPRE_Int            nrows,
  HYPRE_Int           *ncols,
  const HYPRE_BigInt  *rows,
  const HYPRE_BigInt  *cols,
  const HYPRE_Complex *values
)
{
  return HYPRE_IJMatrixAddToValues(*(HYPRE_IJMatrix*) matrix, nrows, ncols, rows, cols, values);
}

int __attribute__((visibility("default"))) IJMatrixCreate
(
  MPI_Comm        comm,
  HYPRE_BigInt    ilower,
  HYPRE_BigInt    iupper,
  HYPRE_BigInt    jlower,
  HYPRE_BigInt    jupper,
  void *matrix
)
{
  return HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, (HYPRE_IJMatrix*) matrix);
}

int __attribute__((visibility("default"))) IJMatrixSetObjectType
(
  void *matrix
)
{
  return HYPRE_IJMatrixSetObjectType(*(HYPRE_IJMatrix*) matrix, HYPRE_PARCSR);
}

int __attribute__((visibility("default"))) IJMatrixInitialize
(
  void *matrix
)
{
  return HYPRE_IJMatrixInitialize(*(HYPRE_IJMatrix*) matrix);
}

int __attribute__((visibility("default"))) IJMatrixAssemble
(
  void *matrix
)
{
  return HYPRE_IJMatrixAssemble(*(HYPRE_IJMatrix*) matrix);
}

}
