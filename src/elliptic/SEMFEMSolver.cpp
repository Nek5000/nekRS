#include "elliptic.h"
#include "SEMFEMSolver.hpp"

SEMFEMSolver_t::SEMFEMSolver_t(elliptic_t *elliptic_)
{
  elliptic = elliptic_;
  const dfloat lambda0 = elliptic->lambda0Avg;

  auto mesh = elliptic->mesh;

  const int verbose = (platform->options.compareArgs("VERBOSE", "TRUE")) ? 1 : 0;
  const int useFP32 = elliptic->options.compareArgs("COARSE SOLVER PRECISION", "FP32");
  const bool useDevice = elliptic->options.compareArgs("COARSE SOLVER LOCATION", "DEVICE");

  MPI_Barrier(platform->comm.mpiComm);
  double tStart = MPI_Wtime();
  if (platform->comm.mpiRank == 0) {
    printf("setup SEMFEM solver (lambdaAvg=%g) ... \n", lambda0);
  }
  fflush(stdout);

  const auto mask = [&]() {
    std::vector<int> mask(mesh->Nlocal, 1);
    if (elliptic->Nmasked > 0) {
      std::vector<int> maskIds(elliptic->o_maskIds.size());
      elliptic->o_maskIds.copyTo(maskIds.data());
      for (dlong i = 0; i < maskIds.size(); i++) {
        mask[maskIds[i]] = 0;
      }
    }
    return mask;
  }();

  const auto matrix = buildMatrix(mesh->Nq,
                                  mesh->Nelements,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mask,
                                  lambda0,
                                  mesh->ogs->hostGsh);

  const auto numRows = static_cast<int>(matrix.rowEnd - matrix.rowStart + 1);

  o_dofMap = platform->device.malloc<dlong>(matrix.dofMap.size());
  o_dofMap.copyFrom(matrix.dofMap.data());

  if (elliptic->options.compareArgs("COARSE SOLVER", "BOOMERAMG")) {
    double settings[hypreWrapper::NPARAM + 1];
    settings[0] = 1;    /* custom settings              */
    settings[1] = 8;    /* coarsening                   */
    settings[2] = 6;    /* interpolation                */
    settings[3] = 1;    /* number of cycles             */
    settings[4] = 18;   /* smoother for crs level       */
    settings[5] = 3;    /* number of coarse sweeps      */
    settings[6] = 18;   /* smoother                     */
    settings[7] = 1;    /* number of sweeps             */
    settings[8] = 0.25; /* strong threshold             */
    settings[9] = 0.05; /* non galerkin tol             */
    settings[10] = 0;   /* aggressive coarsening levels */
    settings[11] = 1;   /* chebyRelaxOrder */
    settings[12] = 0.3; /* chebyRelaxFraction */

    if (elliptic->options.compareArgs("MULTIGRID SEMFEM", "TRUE") && elliptic->mesh->N == 1) {
      settings[6] = 16;
      settings[4] = settings[6];
    }

    platform->options.getArgs("BOOMERAMG COARSEN TYPE", settings[1]);
    platform->options.getArgs("BOOMERAMG INTERPOLATION TYPE", settings[2]);
    platform->options.getArgs("BOOMERAMG COARSE SMOOTHER TYPE", settings[4]);
    platform->options.getArgs("BOOMERAMG SMOOTHER TYPE", settings[6]);
    platform->options.getArgs("BOOMERAMG SMOOTHER SWEEPS", settings[7]);
    platform->options.getArgs("BOOMERAMG ITERATIONS", settings[3]);
    platform->options.getArgs("BOOMERAMG STRONG THRESHOLD", settings[8]);
    platform->options.getArgs("BOOMERAMG NONGALERKIN TOLERANCE", settings[9]);
    platform->options.getArgs("BOOMERAMG AGGRESSIVE COARSENING LEVELS", settings[10]);
    platform->options.getArgs("BOOMERAMG CHEBYSHEV RELAX ORDER", settings[11]);
    platform->options.getArgs("BOOMERAMG CHEBYSHEV FRACTION", settings[12]);

    const auto numRows = o_dofMap.size();

    if (platform->device.mode() != "Serial" && useDevice) {
      boomerAMG = new hypreWrapperDevice::boomerAMG_t(numRows,
                                                      matrix.nnz,
                                                      matrix.Ai.data(),
                                                      matrix.Aj.data(),
                                                      matrix.Av.data(),
                                                      (int)elliptic->nullspace,
                                                      platform->comm.mpiComm,
                                                      platform->device.occaDevice(),
                                                      useFP32,
                                                      settings,
                                                      verbose);
    } else {
      boomerAMG = new hypreWrapper::boomerAMG_t(numRows,
                                                matrix.nnz,
                                                matrix.Ai.data(),
                                                matrix.Aj.data(),
                                                matrix.Av.data(),
                                                (int)elliptic->nullspace,
                                                platform->comm.mpiComm,
                                                1, /* Nthreads */
                                                useFP32,
                                                settings,
                                                verbose);
    }
  } else if (elliptic->options.compareArgs("COARSE SOLVER", "AMGX")) {
    nekrsCheck(platform->device.mode() != "CUDA",
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s\n",
               "AmgX only supports CUDA!");

    std::string configFile;
    elliptic->options.getArgs("AMGX CONFIG FILE", configFile);
    char *cfg = NULL;
    if (configFile.size()) {
      cfg = (char *)configFile.c_str();
    }
    AMGX = new AMGX_t(numRows,
                      matrix.nnz,
                      matrix.Ai.data(),
                      matrix.Aj.data(),
                      matrix.Av.data(),
                      (int)elliptic->nullspace,
                      platform->comm.mpiComm,
                      platform->device.id(),
                      useFP32,
                      std::stoi(getenv("NEKRS_GPU_MPI")),
                      cfg);
  } else {
    std::string amgSolver;
    elliptic->options.getArgs("COARSE SOLVER", amgSolver);
    nekrsAbort(platform->comm.mpiComm,
               EXIT_FAILURE,
               "COARSE SOLVER %s is not supported!\n",
               amgSolver.c_str());
  }

  if (platform->comm.mpiRank == 0) {
    printf("done (%gs)\n", MPI_Wtime() - tStart);
  }
  fflush(stdout);
}

SEMFEMSolver_t::~SEMFEMSolver_t()
{
  const auto useDevice = elliptic->options.compareArgs("COARSE SOLVER LOCATION", "DEVICE");
  if (boomerAMG) {
    if (useDevice) {
      delete (hypreWrapperDevice::boomerAMG_t *)this->boomerAMG;
    } else {
      delete (hypreWrapper::boomerAMG_t *)this->boomerAMG;
    }
  }
  if (AMGX) {
    delete AMGX;
  }

  o_dofMap.free();
}

void SEMFEMSolver_t::run(const occa::memory &o_r, occa::memory &o_z)
{
  const auto useDevice = elliptic->options.compareArgs("COARSE SOLVER LOCATION", "DEVICE");
  const dlong numRows = o_dofMap.size();

  auto o_rT = platform->deviceMemoryPool.reserve<pfloat>(numRows);
  auto o_zT = platform->deviceMemoryPool.reserve<pfloat>(numRows);

  static occa::kernel gatherKernel;
  if (!gatherKernel.isInitialized()) {
    gatherKernel = platform->kernelRequests.load("gather");
  }

  gatherKernel(numRows, o_dofMap, o_r, o_rT); // E->T

  if (elliptic->options.compareArgs("COARSE SOLVER", "BOOMERAMG")) {

    if (!useDevice) {
      static std::vector<pfloat> rT;
      if (rT.size() < numRows) {
        rT.resize(o_rT.size());
      }
      static std::vector<pfloat> zT;
      if (zT.size() < numRows) {
        zT.resize(o_zT.size());
      }

      o_rT.copyTo(rT.data());
      auto boomerAMG = (hypreWrapper::boomerAMG_t *)this->boomerAMG;
      boomerAMG->solve(rT.data(), zT.data());
      o_zT.copyFrom(zT.data());
    } else {
      auto boomerAMG = (hypreWrapperDevice::boomerAMG_t *)this->boomerAMG;
      boomerAMG->solve(o_rT, o_zT);
    }

  } else if (elliptic->options.compareArgs("COARSE SOLVER", "AMGX") && useDevice) {

    AMGX->solve(o_rT.ptr(), o_zT.ptr());

  } else {

    nekrsAbort(platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "Unknown solver!");
  }

  static occa::kernel scatterKernel;
  if (!scatterKernel.isInitialized()) {
    scatterKernel = platform->kernelRequests.load("scatter");
  }
  scatterKernel(numRows, o_dofMap, o_zT, o_z); // T->E

  oogs::startFinish(o_z, 1, 0, ogsPfloat, ogsAdd, elliptic->oogs);
}
