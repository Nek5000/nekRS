#include "elliptic.h"
#include "platform.hpp"
#include <string>
#include "gslib.h"
#include "fem_amg_preco.hpp"
#include "boomerAMG.h"
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

  const long long numRows = data->rowEnd - data->rowStart + 1;
  elliptic->o_dofMap = platform->device.malloc(numRows * sizeof(long long), data->dofMap);
  elliptic->o_SEMFEMBuffer1 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeof(dfloat));
  elliptic->o_SEMFEMBuffer2 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeof(dfloat));

  elliptic->numRowsSEMFEM = numRows;


  std::string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const std::string oklpath = install_dir + "/okl/elliptic/";
  std::string filename = oklpath + "ellipticPreSEMFEM.okl";
  elliptic->preSEMFEMKernel = platform->device.buildKernel(
    filename,
    "preSEMFEM",
    platform->kernelInfo
  );
  filename = oklpath + "ellipticPostSEMFEM.okl";
  elliptic->postSEMFEMKernel = platform->device.buildKernel(
    filename,
    "postSEMFEM",
    platform->kernelInfo
  );

  if(elliptic->options.compareArgs("SEMFEM SOLVER", "BOOMERAMG")){
    boomerAMGSetup(
      numRows,
      data->nnz,
      data->Ai,
      data->Aj,
      data->Av,
      nullspace,
      platform->comm.mpiComm,
      1,
      platform->device.device_id(),
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
      platform->device.device_id(),
      0, // do not use FP32
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