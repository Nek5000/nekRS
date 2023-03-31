#include "nrssys.hpp"
#include "platform.hpp"
#include "elliptic.h"
#include "SEMFEMSolver.hpp"

static occa::kernel gatherKernel;
static occa::kernel scatterKernel;

SEMFEMSolver_t::SEMFEMSolver_t(elliptic_t* elliptic_)
{
  MPI_Barrier(platform->comm.mpiComm);
  double tStart = MPI_Wtime();
  if(platform->comm.mpiRank == 0)
    printf("setup SEMFEM solver ... \n"); fflush(stdout);

  elliptic = elliptic_;

  const int verbose = (platform->options.compareArgs("VERBOSE","TRUE")) ? 1: 0;
  const int useFP32 = elliptic->options.compareArgs("COARSE SOLVER PRECISION", "FP32");
  const bool useDevice = elliptic->options.compareArgs("COARSE SOLVER LOCATION", "DEVICE");

  if(!gatherKernel.isInitialized()) gatherKernel = platform->kernels.get("gather");
  if(!scatterKernel.isInitialized()) scatterKernel = platform->kernels.get("scatter");

  mesh_t* mesh = elliptic->mesh;
  double* mask = (double*) malloc(mesh->Np*mesh->Nelements*sizeof(double));
  for(int i = 0; i < mesh->Np*mesh->Nelements; ++i) mask[i] = 1.0;
  if(elliptic->Nmasked > 0){
    dlong* maskIds = (dlong*) calloc(elliptic->Nmasked, sizeof(dlong));
    elliptic->o_maskIds.copyTo(maskIds, elliptic->Nmasked * sizeof(dlong));
    for (dlong i = 0; i < elliptic->Nmasked; i++) mask[maskIds[i]] = 0.;
    free(maskIds);
  }

  // here we assume lambda0 is constant (in space and time)
  // use first entry of o_lambda as representative value
  pfloat lambda0;
  elliptic->o_lambda0.copyTo(&lambda0, sizeof(pfloat));

  auto hypreIJ = new hypreWrapper::IJ_t();
  matrix_t* matrix = build(
    mesh->Nq,
    mesh->Nelements,
    mesh->o_x,
    mesh->o_y,
    mesh->o_z,
    mask,
    lambda0,
    *hypreIJ,
    platform->comm.mpiComm,
    mesh->globalIds
  );
  free(mask);


  const dlong numRows = matrix->rowEnd - matrix->rowStart + 1;
  numRowsSEMFEM = numRows;

  o_dofMap = platform->device.malloc(numRows * sizeof(long long), matrix->dofMap);

  o_SEMFEMBuffer1 = platform->device.malloc(numRows * sizeof(pfloat));
  o_SEMFEMBuffer2 = platform->device.malloc(numRows * sizeof(pfloat));
  if(!useDevice){
    SEMFEMBuffer1_h_d = (pfloat*) calloc(numRows, sizeof(pfloat));
    SEMFEMBuffer2_h_d = (pfloat*) calloc(numRows, sizeof(pfloat));
  }

  if(elliptic->options.compareArgs("COARSE SOLVER", "BOOMERAMG")){
      double settings[hypreWrapper::NPARAM+1];
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
      settings[11] = 2;    /* chebyRelaxOrder */

      if(elliptic->options.compareArgs("MULTIGRID SEMFEM", "TRUE")) {
        settings[4]  = 16;
        settings[6]  = 16;
      }  

      platform->options.getArgs("BOOMERAMG COARSEN TYPE", settings[1]);
      platform->options.getArgs("BOOMERAMG INTERPOLATION TYPE", settings[2]);
      platform->options.getArgs("BOOMERAMG COARSE SMOOTHER TYPE", settings[4]);
      platform->options.getArgs("BOOMERAMG SMOOTHER TYPE", settings[6]);
      platform->options.getArgs("BOOMERAMG SMOOTHER SWEEPS", settings[7]);
      platform->options.getArgs("BOOMERAMG ITERATIONS", settings[3]);
      platform->options.getArgs("BOOMERAMG STRONG THRESHOLD", settings[8]);
      platform->options.getArgs("BOOMERAMG NONGALERKIN TOLERANCE" , settings[9]);
      platform->options.getArgs("BOOMERAMG AGGRESSIVE COARSENING LEVELS" , settings[10]);
      platform->options.getArgs("BOOMERAMG CHEBYSHEV RELAX ORDER" , settings[11]);

      if(platform->device.mode() != "Serial" && useDevice) {
        boomerAMG = new hypreWrapperDevice::boomerAMG_t(
                        numRows,
                        matrix->nnz,
                        matrix->Ai,
                        matrix->Aj,
                        matrix->Av,
                        (int) elliptic->allNeumann,
                        platform->comm.mpiComm,
                        platform->device.occaDevice(),
                        useFP32,
                        settings,
                        verbose);
      } else {
        boomerAMG = new hypreWrapper::boomerAMG_t(
          numRows,
          matrix->nnz,
          matrix->Ai,
          matrix->Aj,
          matrix->Av,
          (int) elliptic->allNeumann,
          platform->comm.mpiComm,
          1, /* Nthreads */
          useFP32,
          settings,
          verbose 
        );
      }
  }
  else if(elliptic->options.compareArgs("COARSE SOLVER", "AMGX")){
    nrsCheck(platform->device.mode() != "CUDA", platform->comm.mpiComm, EXIT_FAILURE,
             "AmgX only supports CUDA!\n", "");
      
    std::string configFile;
    elliptic->options.getArgs("AMGX CONFIG FILE", configFile);
    char *cfg = NULL;
    if(configFile.size()) cfg = (char*) configFile.c_str();
    AMGX = new AMGX_t(
      numRows,
      matrix->nnz,
      matrix->Ai,
      matrix->Aj,
      matrix->Av,
      (int) elliptic->allNeumann,
      platform->comm.mpiComm,
      platform->device.id(),
      useFP32,
      std::stoi(getenv("NEKRS_GPU_MPI")),
      cfg);
  }
  else {
    std::string amgSolver;
    elliptic->options.getArgs("COARSE SOLVER", amgSolver);
    nrsAbort(platform->comm.mpiComm, EXIT_FAILURE,
             "COARSE SOLVER %s is not supported!\n", amgSolver.c_str(), "");
  }

  free(matrix);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStart); fflush(stdout);
}

SEMFEMSolver_t::~SEMFEMSolver_t()
{
  const auto useDevice = elliptic->options.compareArgs("COARSE SOLVER LOCATION", "DEVICE");
  if(boomerAMG) {
    if(useDevice)
      delete (hypreWrapperDevice::boomerAMG_t*) this->boomerAMG;
    else
      delete (hypreWrapper::boomerAMG_t*) this->boomerAMG;
  }
  if(AMGX) delete AMGX;

  o_dofMap.free();
  o_SEMFEMBuffer1.free();
  o_SEMFEMBuffer2.free();
}

void SEMFEMSolver_t::run(occa::memory& o_r, occa::memory& o_z)
{
  mesh_t* mesh = elliptic->mesh;

  const bool useDevice = elliptic->options.compareArgs("COARSE SOLVER LOCATION", "DEVICE");

  occa::memory& o_bufr = o_SEMFEMBuffer1;
  occa::memory& o_bufz = o_SEMFEMBuffer2;

  // E->T
  gatherKernel(
    numRowsSEMFEM,
    o_dofMap,
    o_r,
    o_bufr
  );

  platform->linAlg->pfill(elliptic->Nfields * elliptic->fieldOffset, 0.0, o_z);

  if(elliptic->options.compareArgs("COARSE SOLVER", "BOOMERAMG")){

    if(!useDevice)
    {
      o_bufr.copyTo(SEMFEMBuffer1_h_d, numRowsSEMFEM * sizeof(pfloat));
      auto boomerAMG = (hypreWrapper::boomerAMG_t*) this->boomerAMG;
      boomerAMG->solve(SEMFEMBuffer1_h_d, SEMFEMBuffer2_h_d);
      o_bufz.copyFrom(SEMFEMBuffer2_h_d, numRowsSEMFEM * sizeof(pfloat));

    } else {
      auto boomerAMG = (hypreWrapperDevice::boomerAMG_t*) this->boomerAMG;
      boomerAMG->solve(o_bufr, o_bufz);
    }

  } else if(elliptic->options.compareArgs("COARSE SOLVER", "AMGX") && useDevice){

    AMGX->solve(o_bufr.ptr(), o_bufz.ptr());

  } else {

    nrsAbort(platform->comm.mpiComm, EXIT_FAILURE,
             "Trying to call an unknown SEMFEM solver!\n", "");

  }

  // T->E
  scatterKernel(
    numRowsSEMFEM,
    o_dofMap,
    o_bufz,
    o_z
  );

  oogs::startFinish(o_z, 1, 0, ogsPfloat, ogsAdd, elliptic->oogs);
}
