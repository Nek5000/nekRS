#include <string>
#include "platform.hpp"
#include "gslib.h"
#include "boomerAMG.h"
#include "elliptic.h"
#include "fem_amg_preco.h"
#include "amgx.h"

void ellipticSEMFEMSetup(elliptic_t* elliptic)
{
  mesh_t* mesh = elliptic->mesh;
  double nullParams[] = {0}; // use default parameters
  double* mask = (double*) malloc(mesh->Np*mesh->Nelements*sizeof(double));
  for(int i = 0; i < mesh->Np*mesh->Nelements; ++i) mask[i] = 1.0;
  for(dlong n = 0; n < elliptic->Nmasked; n++){
    mask[elliptic->maskIds[n]] = 0.0;
  }
  
  SEMFEMData* data = fem_amg_setup(
    mesh->Nq,
    mesh->Nelements,
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
    boomerAMGSetup(
      numRows,
      data->nnz,
      data->Ai,
      data->Aj,
      data->Av,
      (int) elliptic->allNeumann,
      platform->comm.mpiComm,
      1, /* Nthreads */
      platform->device.id(),
      0, /* use FP32 */
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
      (int) elliptic->allNeumann,
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

void ellipticSEMFEMSolve(elliptic_t* elliptic, occa::memory& o_r, occa::memory& o_z)
{
  mesh_t* mesh = elliptic->mesh;

  occa::memory& o_buffer = elliptic->o_SEMFEMBuffer1;
  occa::memory& o_buffer2 = elliptic->o_SEMFEMBuffer2;

  elliptic->gatherKernel(
    elliptic->numRowsSEMFEM,
    elliptic->o_dofMap,
    o_r,
    o_buffer
  );

  platform->linAlg->fill(elliptic->Nfields * elliptic->Ntotal, 0.0, o_z);

  if(elliptic->options.compareArgs("SEMFEM SOLVER", "BOOMERAMG")){
    // TODO: *NOT* device compatible
    boomerAMGSolve(o_buffer2.ptr(), o_buffer.ptr());
  } else {
    //AMGXsolve(o_buffer2.ptr(), o_buffer.ptr());
  }

  elliptic->scatterKernel(
    elliptic->numRowsSEMFEM,
    elliptic->o_dofMap,
    o_buffer2,
    o_z
  );

  oogs::startFinish(o_z, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
}
