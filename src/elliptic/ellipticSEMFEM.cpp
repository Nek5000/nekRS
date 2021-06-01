#include <string>
#include "platform.hpp"
#include "gslib.h"
#include "boomerAMG.h"
#include "elliptic.h"
#include "ellipticBuildSEMFEM.hpp"
#include "amgx.h"

void ellipticSEMFEMSetup(elliptic_t* elliptic)
{

  const int useFP32 = elliptic->options.compareArgs("SEMFEM SOLVER PRECISION", "FP32");
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

  MPI_Barrier(platform->comm.mpiComm);
  double tStart = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("setup SEMFEM preconditioner ... \n"); fflush(stdout);

  mesh_t* mesh = elliptic->mesh;
  double* mask = (double*) malloc(mesh->Np*mesh->Nelements*sizeof(double));
  for(int i = 0; i < mesh->Np*mesh->Nelements; ++i) mask[i] = 1.0;
  for(dlong n = 0; n < elliptic->Nmasked; n++){
    mask[elliptic->maskIds[n]] = 0.0;
  }
  
  SEMFEMData* data = ellipticBuildSEMFEM(
    mesh->Nq,
    mesh->Nelements,
    mesh->o_x,
    mesh->o_y,
    mesh->o_z,
    mask,
    platform->comm.mpiComm,
    mesh->globalIds
  );

  if(platform->options.compareArgs("BUILD ONLY", "TRUE")) return;

  const int sizeType = useFP32 ? sizeof(float) : sizeof(dfloat);
  const long long numRows = data->rowEnd - data->rowStart + 1;
  elliptic->o_dofMap = platform->device.malloc(numRows * sizeof(long long), data->dofMap);
  elliptic->o_SEMFEMBuffer1 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeType);
  elliptic->o_SEMFEMBuffer2 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeType);

  elliptic->numRowsSEMFEM = numRows;

  int setupRetVal;
  if(elliptic->options.compareArgs("SEMFEM SOLVER", "BOOMERAMG")){
      double settings[BOOMERAMG_NPARAM+1];
      settings[0]  = 1;    /* custom settings              */
      settings[1]  = 8;    /* coarsening                   */
      settings[2]  = 6;    /* interpolation                */
      settings[3]  = 1;    /* number of cycles             */
      settings[4]  = 18;   /* smoother for crs level       */
      settings[5]  = 3;    /* number of coarse sweeps      */
      settings[6]  = 18;   /* smoother                     */
      settings[7]  = 1;    /* number of sweeps             */
      settings[8]  = 0.25; /* strong threshold             */
      settings[9]  = 0.1;  /* non galerkin tol             */
      settings[10] = 0;    /* aggressive coarsening levels */

      if(platform->device.mode() != "Serial") {
        if(platform->comm.mpiRank == 0) printf("HYPRE has not been built with GPU support!\n");
        MPI_Barrier(platform->comm.mpiComm);
        ABORT(1);
      } 

      setupRetVal = boomerAMGSetup(
        numRows,
        data->nnz,
        data->Ai,
        data->Aj,
        data->Av,
        (int) elliptic->allNeumann,
        platform->comm.mpiComm,
        1, /* Nthreads */
        platform->device.id(),
        0, /* use FP32 - hardwired as no runtime switch is available */
        settings 
      );
  }
  else if(elliptic->options.compareArgs("SEMFEM SOLVER", "AMGX")){
    if(platform->device.mode() != "CUDA") {
      if(platform->comm.mpiRank == 0) printf("AmgX only supports CUDA!\n");
      MPI_Barrier(platform->comm.mpiComm);
      ABORT(1);
    } 
      
    string configFile;
    elliptic->options.getArgs("AMGX CONFIG FILE", configFile);
    char *cfg = NULL;
    if(configFile.size()) cfg = (char*) configFile.c_str();
    setupRetVal = AMGXsetup(
      numRows,
      data->nnz,
      data->Ai,
      data->Aj,
      data->Av,
      (int) elliptic->allNeumann,
      platform->comm.mpiComm,
      platform->device.id(),
      useFP32,
      std::stoi(getenv("NEKRS_GPU_MPI")),
      cfg);
  }
  else {
    if(platform->comm.mpiRank == 0){
      std::string amgSolver;
      elliptic->options.getArgs("SEMFEM SOLVER", amgSolver);
      printf("SEMFEM SOLVER %s is not supported!\n", amgSolver.c_str());
    }
    ABORT(EXIT_FAILURE);
  }
  if(setupRetVal) MPI_Abort(platform->comm.mpiComm, 1);

  delete data;
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStart); fflush(stdout);
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
    boomerAMGSolve(o_buffer2.ptr(), o_buffer.ptr());
  } else {
    AMGXsolve(o_buffer2.ptr(), o_buffer.ptr());
  }

  elliptic->scatterKernel(
    elliptic->numRowsSEMFEM,
    elliptic->o_dofMap,
    o_buffer2,
    o_z
  );

  oogs::startFinish(o_z, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
}
