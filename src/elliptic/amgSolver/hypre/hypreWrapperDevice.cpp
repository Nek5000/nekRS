#include <mpi.h>
#include <cstdio>
#include <map>
#include <vector>
#include "occa.hpp"

#ifdef ENABLE_HYPRE_GPU 

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_utilities.h"

#define BOOMERAMG_NPARAM 10
static double boomerAMGParam[BOOMERAMG_NPARAM];

class hypre_data_t;

struct hypre_data_t {
  MPI_Comm comm;
  HYPRE_Solver solver;
  HYPRE_IJMatrix A;
  HYPRE_IJVector b;
  HYPRE_IJVector x;
  HYPRE_BigInt iupper;
  HYPRE_BigInt ilower;
  occa::device device;
  int nRows;
};
static hypre_data_t *data;

namespace hypreWrapperDevice
{

int __attribute__((visibility("default")))
BoomerAMGSetup(int nrows, int nz,
               const long long int * Ai, const long long int * Aj, const double * Av,
               int null_space, MPI_Comm ce, occa::device device,
               int useFP32, const double *param, int verbose)
{
  MPI_Comm comm;
  MPI_Comm_dup(ce, &comm);

  int rank;
  MPI_Comm_rank(comm,&rank);

  if(sizeof(HYPRE_Real) != ((useFP32) ? sizeof(float) : sizeof(double))) {
    if(rank == 0) printf("hypreWrapperDevice: HYPRE floating point precision does not match!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  } 

  if ((int) param[0]) {
    for (int i = 0; i < BOOMERAMG_NPARAM; i++)
      boomerAMGParam[i] = param[i+1]; 
  } else {
    boomerAMGParam[0]  = 8;    /* coarsening */
    boomerAMGParam[1]  = 6;    /* interpolation */
    boomerAMGParam[2]  = 1;    /* number of cycles */
    boomerAMGParam[3]  = 8;    /* smoother for crs level */
    boomerAMGParam[4]  = 3;    /* sweeps */
    boomerAMGParam[5]  = 8;    /* smoother */
    boomerAMGParam[6]  = 1;    /* sweeps   */
    boomerAMGParam[7]  = 0.25; /* threshold */
    boomerAMGParam[8]  = 0.0;  /* non galerkin tolerance */
    boomerAMGParam[9]  = 0;    /* agressive coarsening */
  }

  data = new hypre_data_t();
  data->comm = comm;
  data->device = device;
  data->nRows = nrows;
  data->ilower = data->nRows;
  MPI_Scan(MPI_IN_PLACE, &data->ilower, 1, MPI_LONG_LONG, MPI_SUM, ce);
  data->ilower -= data->nRows;
  data->iupper = (data->ilower + data->nRows) - 1; 

  HYPRE_Init();

  HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
  HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);

  HYPRE_Int spgemm_use_vendor = 0; 
#if defined(HYPRE_USING_HIP) || defined(HYPRE_USING_SYCL)
  use_vendor = 1;

  HYPRE_Int  spgemm_alg = 1;
  HYPRE_Int  spgemm_rowest_mtd = 3;
  HYPRE_Int  spgemm_rowest_nsamples = 32;
  HYPRE_Real spgemm_rowest_mult = 1.5;
  char       spgemm_hash_type = 'L';

  hypre_SetSpGemmAlgorithm(spgemm_alg);
  hypre_SetSpGemmRownnzEstimateMethod(spgemm_rowest_mtd);
  hypre_SetSpGemmRownnzEstimateNSamples(spgemm_rowest_nsamples);
  hypre_SetSpGemmRownnzEstimateMultFactor(spgemm_rowest_mult);
  hypre_SetSpGemmHashType(spgemm_hash_type);
#endif
  HYPRE_SetSpGemmUseVendor(spgemm_use_vendor);
  HYPRE_SetSpMVUseVendor(1);

  HYPRE_SetUseGpuRand(1);

  // Setup matrix
  {
    HYPRE_IJMatrixCreate(comm,data->ilower,data->iupper,data->ilower,data->iupper,&data->A);
    HYPRE_IJMatrixSetObjectType(data->A,HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(data->A);

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

    auto o_ncols = data->device.malloc(ncols.size() * sizeof(HYPRE_Int),ncols.data());
    auto o_rows = data->device.malloc(rows.size() * sizeof(HYPRE_BigInt),rows.data());
    auto o_cols = data->device.malloc(cols.size() * sizeof(HYPRE_BigInt),cols.data());
    auto o_vals = data->device.malloc(vals.size() * sizeof(HYPRE_Real),vals.data());

    HYPRE_IJMatrixSetValues(data->A, 
                              rowsToSet /* values for nrows */, 
                              (HYPRE_Int*) o_ncols.ptr() /* cols for each row */, 
                              (HYPRE_BigInt*) o_rows.ptr(),
                              (HYPRE_BigInt*) o_cols.ptr(),
                              (HYPRE_Real*) o_vals.ptr());

    o_ncols.free();
    o_rows.free();
    o_cols.free();
    o_vals.free();

    HYPRE_IJMatrixAssemble(data->A);
#if 0
    HYPRE_IJMatrixPrint(data->A, "matrix.dat");
#endif
  }

  // Setup solver
  HYPRE_BoomerAMGCreate(&data->solver);

  HYPRE_BoomerAMGSetCoarsenType(data->solver,boomerAMGParam[0]);
  HYPRE_BoomerAMGSetInterpType(data->solver,boomerAMGParam[1]);

  HYPRE_BoomerAMGSetModuleRAP2(data->solver, 1);
  HYPRE_BoomerAMGSetKeepTranspose(data->solver, 1);

  //HYPRE_BoomerAMGSetChebyFraction(data->solver, 0.2); 

  if (boomerAMGParam[5] > 0) {
    HYPRE_BoomerAMGSetCycleRelaxType(data->solver, boomerAMGParam[5], 1);
    HYPRE_BoomerAMGSetCycleRelaxType(data->solver, boomerAMGParam[5], 2);
  } 
  HYPRE_BoomerAMGSetCycleRelaxType(data->solver, 9, 3);

  HYPRE_BoomerAMGSetCycleNumSweeps(data->solver, boomerAMGParam[6], 1);
  HYPRE_BoomerAMGSetCycleNumSweeps(data->solver, boomerAMGParam[6], 2);
  HYPRE_BoomerAMGSetCycleNumSweeps(data->solver, 1, 3);

  if (null_space) {
    HYPRE_BoomerAMGSetMinCoarseSize(data->solver, 2);
    HYPRE_BoomerAMGSetCycleRelaxType(data->solver, boomerAMGParam[3], 3);
    HYPRE_BoomerAMGSetCycleNumSweeps(data->solver, boomerAMGParam[4], 3);
  }

  HYPRE_BoomerAMGSetStrongThreshold(data->solver,boomerAMGParam[7]);

// NonGalerkin not supported yet
#if 0
  if (boomerAMGParam[8] > 1e-3) {
    HYPRE_BoomerAMGSetNonGalerkinTol(data->solver,boomerAMGParam[8]);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(data->solver,0.0 , 0);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(data->solver,0.01, 1);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(data->solver,0.05, 2);
  }
#endif
  
  if(boomerAMGParam[9] > 0) {
    HYPRE_BoomerAMGSetAggNumLevels(data->solver, boomerAMGParam[9]); 
    HYPRE_BoomerAMGSetAggInterpType(data->solver, 5);
    //HYPRE_BoomerAMGSetNumPaths(data->solver, 10);
  }

  HYPRE_BoomerAMGSetMaxIter(data->solver, boomerAMGParam[2]);
  HYPRE_BoomerAMGSetTol(data->solver,0);

  if(verbose)
    HYPRE_BoomerAMGSetPrintLevel(data->solver,3);
  else
    HYPRE_BoomerAMGSetPrintLevel(data->solver,1);

  // Create and initialize rhs and solution vectors
  HYPRE_IJVectorCreate(comm,data->ilower,data->iupper,&data->b);
  HYPRE_IJVectorSetObjectType(data->b,HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(data->b);
  HYPRE_IJVectorAssemble(data->b);

  HYPRE_IJVectorCreate(comm,data->ilower,data->iupper,&data->x);
  HYPRE_IJVectorSetObjectType(data->x,HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(data->x);
  HYPRE_IJVectorAssemble(data->x);

  // Perform AMG setup
  HYPRE_ParVector par_b;
  HYPRE_ParVector par_x;
  HYPRE_ParCSRMatrix par_A;
  HYPRE_IJVectorGetObject(data->b,(void**) &par_b);
  HYPRE_IJVectorGetObject(data->x,(void**) &par_x);
  HYPRE_IJMatrixGetObject(data->A,(void**) &par_A);
  HYPRE_BoomerAMGSetup(data->solver,par_A,par_b,par_x);

  return 0;
}

int __attribute__((visibility("default")))
BoomerAMGSolve(const occa::memory& o_b, const occa::memory& o_x)
{
  data->device.finish(); // input buffers ready

#if 1
  HYPRE_IJVectorUpdateValues(data->x,data->nRows,NULL,(HYPRE_Real*) o_x.ptr(), 1);
#else
  HYPRE_IJVectorSetValues(data->x,data->nRows,NULL,(HYPRE_Real*) o_x.ptr());
  HYPRE_IJVectorAssemble(data->x);
#endif

#if 1
  HYPRE_IJVectorUpdateValues(data->b,data->nRows,NULL,(HYPRE_Real*) o_b.ptr(), 1);
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

#if 0
  HYPRE_IJVectorPrint(data->b, "b.dat");
  HYPRE_IJVectorPrint(data->x, "x.dat");
#endif

  const int retVal = HYPRE_BoomerAMGSolve(data->solver,par_A,par_b,par_x);
  if(retVal > 0) { 
    int rank;
    MPI_Comm_rank(data->comm,&rank);
    if(rank == 0) printf("HYPRE_BoomerAMGSolve failed with retVal=%d!\n", retVal);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // sync copy (blocks host until buffer is ready) 
  HYPRE_IJVectorGetValues(data->x,data->nRows,NULL,(HYPRE_Real*) o_x.ptr());

  return 0; 
}

void __attribute__((visibility("default"))) 
Free()
{
  HYPRE_BoomerAMGDestroy(data->solver);
  HYPRE_IJMatrixDestroy(data->A);
  HYPRE_IJVectorDestroy(data->x);
  HYPRE_IJVectorDestroy(data->b);
  HYPRE_Finalize();
  free(data);
}

} // namespace

#else

namespace hypreWrapperDevice
{

int __attribute__((visibility("default"))) 
BoomerAMGSetup(int nrows, int nz,
               const long long int * Ai, const long long int * Aj, const double * Av,
               int null_space, MPI_Comm ce, occa::device device,
               int useFP32, const double *param, int verbose)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);  
  if(rank == 0) printf("ERROR: Recompile with HYPRE GPU support!\n");
  return 1;
}

int __attribute__((visibility("default"))) 
BoomerAMGSolve(const occa::memory& o_b, const occa::memory& o_x)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);  
  if(rank == 0) printf("ERROR: Recompile with HYPRE GPU support!\n");
  return 1;
}

void __attribute__((visibility("default"))) 
Free()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);  
  if(rank == 0) printf("ERROR: Recompile with HYPRE GPU support!\n");
}

} // namespace

#endif
