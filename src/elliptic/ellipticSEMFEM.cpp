#include <string>
#include "platform.hpp"
#include "gslib.h"
#include "boomerAMG.h"
#include "elliptic.h"
#include "ellipticBuildSEMFEM.hpp"
#include "amgx.h"

namespace{
occa::kernel gatherKernel;
occa::kernel scatterKernel;
occa::memory o_dofMap;
occa::memory o_SEMFEMBuffer1;
occa::memory o_SEMFEMBuffer2;
double* SEMFEMBuffer1_h_d;
double* SEMFEMBuffer2_h_d;
dlong numRowsSEMFEM;
}

void ellipticSEMFEMSetup(elliptic_t* elliptic)
{
  const int verbose = (platform->options.compareArgs("VERBOSE","TRUE")) ? 1: 0;
  const int useFP32 = elliptic->options.compareArgs("SEMFEM SOLVER PRECISION", "FP32");
  gatherKernel = platform->kernels.get("gather");
  scatterKernel = platform->kernels.get("scatter");

  MPI_Barrier(platform->comm.mpiComm);
  double tStart = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("setup SEMFEM preconditioner ... \n"); fflush(stdout);

  mesh_t* mesh = elliptic->mesh;
  double* mask = (double*) malloc(mesh->Np*mesh->Nelements*sizeof(double));
  for(int i = 0; i < mesh->Np*mesh->Nelements; ++i) mask[i] = 1.0;
  if(elliptic->Nmasked > 0){
    dlong* maskIds = (dlong*) calloc(elliptic->Nmasked, sizeof(dlong));
    elliptic->o_maskIds.copyTo(maskIds, elliptic->Nmasked * sizeof(dlong));
    for (dlong i = 0; i < elliptic->Nmasked; i++) mask[maskIds[i]] = 0.;
    free(maskIds);
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

  const int sizeType = useFP32 ? sizeof(float) : sizeof(dfloat);
  const long long numRows = data->rowEnd - data->rowStart + 1;
  o_dofMap = platform->device.malloc(numRows * sizeof(long long), data->dofMap);
  o_SEMFEMBuffer1 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeType);
  o_SEMFEMBuffer2 = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal,sizeType);

  const bool useDevice = platform->options.compareArgs("AMG SOLVER LOCATION", "GPU");
  if(!useDevice){
    SEMFEMBuffer1_h_d = (dfloat*) calloc(elliptic->Nfields * elliptic->Ntotal, sizeof(dfloat));
    SEMFEMBuffer2_h_d = (dfloat*) calloc(elliptic->Nfields * elliptic->Ntotal, sizeof(dfloat));
  }

  numRowsSEMFEM = numRows;

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
      settings[9]  = 0.05; /* non galerkin tol             */
      settings[10] = 0;    /* aggressive coarsening levels */

      platform->options.getArgs("BOOMERAMG COARSEN TYPE", settings[1]);
      platform->options.getArgs("BOOMERAMG INTERPOLATION TYPE", settings[2]);
      platform->options.getArgs("BOOMERAMG SMOOTHER TYPE", settings[6]);
      platform->options.getArgs("BOOMERAMG SMOOTHER SWEEPS", settings[7]);
      platform->options.getArgs("BOOMERAMG ITERATIONS", settings[3]);
      platform->options.getArgs("BOOMERAMG STRONG THRESHOLD", settings[8]);
      platform->options.getArgs("BOOMERAMG NONGALERKIN TOLERANCE" , settings[9]);
      platform->options.getArgs("BOOMERAMG AGGRESSIVE COARSENING LEVELS" , settings[10]);


      if(useFP32)
      {
        if(platform->comm.mpiRank == 0) printf("HYPRE does not support FP32!\n");
        MPI_Barrier(platform->comm.mpiComm);
        ABORT(1);
      }

      if(platform->device.mode() != "Serial" && useDevice) {
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
        useDevice ? platform->device.id() : -1,
        0, /* do not use FP32 - hardwired as no runtime switch is available */
        settings,
        verbose 
      );
  }
  else if(elliptic->options.compareArgs("SEMFEM SOLVER", "AMGX")){
    if(platform->device.mode() != "CUDA") {
      if(platform->comm.mpiRank == 0) printf("AmgX only supports CUDA!\n");
      MPI_Barrier(platform->comm.mpiComm);
      ABORT(1);
    } 
      
    std::string configFile;
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

  free(data);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStart); fflush(stdout);
}

void ellipticSEMFEMSolve(elliptic_t* elliptic, occa::memory& o_r, occa::memory& o_z)
{
  mesh_t* mesh = elliptic->mesh;

  occa::memory& o_buffer = o_SEMFEMBuffer1;
  occa::memory& o_buffer2 = o_SEMFEMBuffer2;

  gatherKernel(
    numRowsSEMFEM,
    o_dofMap,
    o_r,
    o_buffer
  );

  platform->linAlg->fill(mesh->Np * mesh->Nelements, 0.0, o_z);

  if(elliptic->options.compareArgs("SEMFEM SOLVER", "BOOMERAMG")){

    const bool useDevice = platform->options.compareArgs("AMG SOLVER LOCATION", "GPU");
    if(platform->device.mode() != "Serial" && !useDevice)
    {
      o_buffer.copyTo(SEMFEMBuffer1_h_d, elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));
      boomerAMGSolve(SEMFEMBuffer2_h_d, SEMFEMBuffer1_h_d);
      o_buffer2.copyFrom(SEMFEMBuffer2_h_d, elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));

    } else {
      boomerAMGSolve(o_buffer2.ptr(), o_buffer.ptr());
    }
  } else {
    AMGXsolve(o_buffer2.ptr(), o_buffer.ptr());
  }

  scatterKernel(
    numRowsSEMFEM,
    o_dofMap,
    o_buffer2,
    o_z
  );

  oogs::startFinish(o_z, 1, mesh->Np * mesh->Nelements, ogsDfloat, ogsAdd, elliptic->oogs);
}
