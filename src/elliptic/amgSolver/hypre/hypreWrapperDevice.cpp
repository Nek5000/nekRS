#include <mpi.h>
#include <cstdio>
#include <map>
#include <vector>

#ifdef ENABLE_HYPRE_GPU
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_utilities.h"
#endif

#include "hypreWrapperDevice.hpp"

#define DEBUG 0
#ifdef ENABLE_HYPRE_GPU

static int HYPREinit = 0;

namespace hypreWrapperDevice {

__attribute__((visibility("default"))) 
boomerAMG_t::boomerAMG_t(int _nRows,
                         int nz,
                         const long long int *Ai,
                         const long long int *Aj,
                         const double *Av,
                         int null_space,
                         MPI_Comm ce,
                         occa::device _device,
                         int useFP32,
                         const double *param,
                         int verbose)
{
  MPI_Comm_dup(ce, &comm);

  device = _device;
  nRows = _nRows;

  MPI_Comm_rank(comm, &rank);

  if (sizeof(HYPRE_Real) != ((useFP32) ? sizeof(float) : sizeof(double))) {
    if (rank == 0)
      printf("hypreWrapperDevice: HYPRE floating point precision does not match!\n");
    MPI_Abort(comm, 1);
  }

  if ((int)param[0]) {
    for (int i = 0; i < NPARAM; i++)
      params[i] = param[i + 1];
  }
  else {
    params[0] = 8;    /* coarsening */
    params[1] = 6;    /* interpolation */
    params[2] = 1;    /* number of cycles */
    params[3] = 8;    /* smoother for crs level */
    params[4] = 3;    /* sweeps */
    params[5] = 8;    /* smoother */
    params[6] = 1;    /* sweeps   */
    params[7] = 0.25; /* threshold */
    params[8] = 0.0;  /* non galerkin tolerance */
    params[9] = 0;    /* agressive coarsening */
    params[10] = 2;    /* chebyRelaxOrder */
  }

  long long rowStart = nRows;
  MPI_Scan(MPI_IN_PLACE, &rowStart, 1, MPI_LONG_LONG, MPI_SUM, comm);
  rowStart -= nRows;

  const HYPRE_BigInt ilower = (HYPRE_BigInt)rowStart;
  const HYPRE_BigInt iupper = ilower + (HYPRE_BigInt)nRows - 1;

  if(!HYPREinit) {
    HYPRE_Init();

    HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
    HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);

    HYPRE_Int spgemm_use_vendor = 0;
#if defined(HYPRE_USING_HIP) || defined(HYPRE_USING_SYCL)
    use_vendor = 1;
 
    HYPRE_Int spgemm_alg = 1;
    HYPRE_Int spgemm_rowest_mtd = 3;
    HYPRE_Int spgemm_rowest_nsamples = 32;
    HYPRE_Real spgemm_rowest_mult = 1.5;
    char spgemm_hash_type = 'L';
 
    hypre_SetSpGemmAlgorithm(spgemm_alg);
    hypre_SetSpGemmRownnzEstimateMethod(spgemm_rowest_mtd);
    hypre_SetSpGemmRownnzEstimateNSamples(spgemm_rowest_nsamples);
    hypre_SetSpGemmRownnzEstimateMultFactor(spgemm_rowest_mult);
    hypre_SetSpGemmHashType(spgemm_hash_type);
#endif
    HYPRE_SetSpGemmUseVendor(spgemm_use_vendor);
    HYPRE_SetSpMVUseVendor(1);

    HYPRE_SetUseGpuRand(1);
  }


  // Setup matrix
  {
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

    auto o_ncols = device.malloc(ncols.size() * sizeof(HYPRE_Int), ncols.data());
    auto o_rows = device.malloc(rows.size() * sizeof(HYPRE_BigInt), rows.data());
    auto o_cols = device.malloc(cols.size() * sizeof(HYPRE_BigInt), cols.data());
    auto o_vals = device.malloc(vals.size() * sizeof(HYPRE_Real), vals.data());

    HYPRE_IJMatrixSetValues(*A,
                            rowsToSet /* values for nrows */,
                            (HYPRE_Int *)o_ncols.ptr() /* cols for each row */,
                            (HYPRE_BigInt *)o_rows.ptr(),
                            (HYPRE_BigInt *)o_cols.ptr(),
                            (HYPRE_Real *)o_vals.ptr());

    o_ncols.free();
    o_rows.free();
    o_cols.free();
    o_vals.free();

    HYPRE_IJMatrixAssemble(*A);

    if(DEBUG)
      HYPRE_IJMatrixPrint(*A, "matrix.dat");
  }

  solver = new HYPRE_Solver();
  HYPRE_BoomerAMGCreate(solver);

  HYPRE_BoomerAMGSetCoarsenType(*solver, params[0]);
  HYPRE_BoomerAMGSetInterpType(*solver, params[1]);

  HYPRE_BoomerAMGSetModuleRAP2(*solver, 1);
  HYPRE_BoomerAMGSetKeepTranspose(*solver, 1);

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

// NonGalerkin not supported yet
#if 0
  if (params[8] > 1e-3) {
    HYPRE_BoomerAMGSetNonGalerkinTol(*solver,params[8]);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(*solver,0.0 , 0);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(*solver,0.01, 1);
    HYPRE_BoomerAMGSetLevelNonGalerkinTol(*solver,0.05, 2);
  }
#endif

  if (params[9] > 0) {
    HYPRE_BoomerAMGSetAggNumLevels(*solver, params[9]);
    HYPRE_BoomerAMGSetAggInterpType(*solver, 5);
    // HYPRE_BoomerAMGSetNumPaths(*solver, 10);
  }

  HYPRE_BoomerAMGSetMaxIter(*solver, params[2]);
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
  HYPRE_ParCSRMatrix par_A;
  HYPRE_IJVectorGetObject(*b, (void **)&par_b);
  HYPRE_IJVectorGetObject(*x, (void **)&par_x);
  HYPRE_IJMatrixGetObject(*A, (void **)&par_A);
  const int err = HYPRE_BoomerAMGSetup(*solver, par_A, par_b, par_x);
  if (err > 0) {
    if (rank == 0)
      printf("HYPRE_BoomerAMGDeviceSetup failed with %d!\n", err);
    MPI_Abort(comm, err);
  }

  HYPREinit = 1;
}

void __attribute__((visibility("default")))
boomerAMG_t::solve(const occa::memory &o_b, const occa::memory &o_x)
{
  device.finish(); // input buffers ready

  HYPRE_IJVectorUpdateValues(*x, nRows, NULL, (HYPRE_Real *)o_x.ptr(), 1);
  HYPRE_IJVectorUpdateValues(*b, nRows, NULL, (HYPRE_Real *)o_b.ptr(), 1);

  HYPRE_ParVector par_x;
  HYPRE_ParVector par_b;
  HYPRE_ParCSRMatrix par_A;

  HYPRE_IJVectorGetObject(*x, (void **)&par_x);
  HYPRE_IJVectorGetObject(*b, (void **)&par_b);
  HYPRE_IJMatrixGetObject(*A, (void **)&par_A);

#if 0
  HYPRE_IJVectorPrint(b, "b.dat");
  HYPRE_IJVectorPrint(x, "x.dat");
#endif

  const int err = HYPRE_BoomerAMGSolve(*solver, par_A, par_b, par_x);
  if (err > 0) {
    if (rank == 0)
      printf("HYPRE_BoomerAMGDeviceSolve failed with %d!\n", err);
    MPI_Abort(comm, err);
  }

  // sync copy (blocks host until buffer is ready)
  HYPRE_IJVectorGetValues(*x, nRows, NULL, (HYPRE_Real *)o_x.ptr());
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

int __attribute__((visibility("default"))) enabled() { return 1; }

} // namespace hypreWrapperDevice

#else

namespace hypreWrapperDevice {

__attribute__((visibility("default"))) 
boomerAMG_t::boomerAMG_t(int nrows,
                   int nz,
                   const long long int *Ai,
                   const long long int *Aj,
                   const double *Av,
                   int null_space,
                   MPI_Comm ce,
                   occa::device device,
                   int useFP32,
                   const double *param,
                   int verbose)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    printf("ERROR: Recompile with HYPRE GPU support!\n");
  MPI_Abort(MPI_COMM_WORLD, 1);
}

void __attribute__((visibility("default")))
boomerAMG_t::solve(const occa::memory &o_b, const occa::memory &o_x)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    printf("ERROR: Recompile with HYPRE GPU support!\n");
  MPI_Abort(MPI_COMM_WORLD, 1);
}

__attribute__((visibility("default"))) boomerAMG_t::~boomerAMG_t()
{
}

void __attribute__((visibility("default"))) finalize()
{
}

int __attribute__((visibility("default"))) enabled() { return 0; }

} // namespace hypreWrapperDevice

#endif
