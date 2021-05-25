#include "elliptic.h"
#include "platform.hpp"
#include <string>
#include "gslib.h"
#include "fem_amg_preco.hpp"
#include "boomerAMG.h"
#include "HYPRE_utilities.h"
#if 0
#include "amgx.h"
#endif

void ellipticSEMFEMSetup(elliptic_t* elliptic)
{
  mesh_t* mesh = elliptic->mesh;
  const int N = mesh->Nq;
  const int Nelements = mesh->Nelements;
  const int nullspace = (int) elliptic->allNeumann;
  double nullParams[] = {0}; // use default parameters
  double* mask = (double*) malloc(N*N*N*Nelements*sizeof(double));
  for(int i = 0; i < N*N*N*Nelements; ++i) mask[i] = 1.0;
  for(dlong n = 0; n < elliptic->Nmasked; n++){
    mask[elliptic->maskIds[n]] = 0.0;
  }
  
  SEMFEMData* data = fem_amg_setup(
    &N,
    &N,
    &N,
    &Nelements,
    mesh->x,
    mesh->y,
    mesh->z,
    mask,
    platform->comm.mpiComm,
    mesh->globalIds
  );

  const int useFP32 = elliptic->options.compareArgs("SEMFEM SOLVER PRECISION", "FP32");
  const int sizeType = useFP32 ? sizeof(float) : sizeof(dfloat);
  const long long numRows = data->rowEnd - data->rowStart + 1;
  elliptic->o_dofMap = platform->device.malloc(numRows * sizeof(long long), data->dofMap);
  elliptic->o_SEMFEMBuffer1 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeType);
  elliptic->o_SEMFEMBuffer2 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeType);

  elliptic->numRowsSEMFEM = numRows;


  occa::properties SEMFEMKernelProps = platform->kernelInfo;
  if(useFP32){
    SEMFEMKernelProps["defines/" "pfloat"] = "float";
  } else {
    SEMFEMKernelProps["defines/" "pfloat"] = "double";
  }
  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = install_dir + "/okl/elliptic/";
  std::string filename = oklpath + "ellipticGather.okl";
  elliptic->gatherKernel = platform->device.buildKernel(
    filename,
    "gather",
    SEMFEMKernelProps
  );
  filename = oklpath + "ellipticScatter.okl";
  elliptic->scatterKernel = platform->device.buildKernel(
    filename,
    "scatter",
    SEMFEMKernelProps
  );

  if(elliptic->options.compareArgs("SEMFEM SOLVER", "BOOMERAMG")){
    if(useFP32){
      if(sizeof(HYPRE_Real) != sizeof(float)){
        if(platform->comm.mpiRank == 0)
          printf("HYPRE has not been built to support FP32.\n");
        ABORT(EXIT_FAILURE);
      }
    } else {
      if(sizeof(HYPRE_Real) != sizeof(double)){
        if(platform->comm.mpiRank == 0)
          printf("HYPRE has not been built to support FP64.\n");
        ABORT(EXIT_FAILURE);
      }
    }
    boomerAMGSetup(
      numRows,
      data->nnz,
      data->Ai,
      data->Aj,
      data->Av,
      nullspace,
      platform->comm.mpiComm,
      1,
      platform->device.id(),
      &nullParams[0]
    );
  }
  else if(elliptic->options.compareArgs("SEMFEM SOLVER", "AMGX")){
#if 1
    if(platform->comm.mpiRank == 0){
      printf("AMGX solver is not yet supported!\n");
    }
    ABORT(EXIT_FAILURE);
#else
    AMGXsetup(
      numRows,
      data->nnz,
      data->Ai,
      data->Aj,
      data->Av,
      nullspace,
      platform->comm.mpiComm,
      platform->device.id(),
      useFP32,
      nullptr // use default parameters
    );
#endif
  } else {
    if(platform->comm.mpiRank == 0){
      std::string amgSolver;
      elliptic->options.getArgs("SEMFEM SOLVER", amgSolver);
      printf("SEMFEM SOLVER %s is not supported!\n", amgSolver.c_str());
    }
    ABORT(EXIT_FAILURE);
  }

  delete data;
}