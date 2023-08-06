#include <math.h>
#include <limits>
#include <vector>
#include <numeric>

static int ciMode = 0;

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", std::string("0.05"));
  options.setArgs("DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("3e-3"));
  options.setArgs("DT", std::string("1e-4"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("PRESSURE INITIAL GUESS", "PREVIOUS");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");

  options.setArgs("BDF ORDER", "3");
  if (ciMode == 1 || ciMode == 3) {
    options.setArgs("NEKNEK CORRECTOR STEPS", "1");
    options.setArgs("NEKNEK BOUNDARY EXT ORDER", "3");
  } else if (ciMode == 2) {
    options.setArgs("NEKNEK CORRECTOR STEPS", "0");
    options.setArgs("NEKNEK BOUNDARY EXT ORDER", "1");
  }
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-08"));

  if (ciMode == 2) {
    int fail = 0;

    MPI_Allreduce(MPI_IN_PLACE, &fail, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);
    if (fail) {
      CIFAIL;
    }
  }
}

void ciTestPointInterp(nrs_t *nrs, dfloat time, int tstep)
{
  // query velocity solution from inner session
  std::vector<dfloat> x = {M_PI};
  std::vector<dfloat> y = {M_PI};
  std::vector<dfloat> z = {0.05};
  int n = x.size();

  auto o_xinterp = platform->device.malloc(n * sizeof(dfloat), x.data());
  auto o_yinterp = platform->device.malloc(n * sizeof(dfloat), y.data());
  auto o_zinterp = platform->device.malloc(n * sizeof(dfloat), z.data());
  auto o_UouterExact = platform->device.malloc(nrs->NVfields * n * sizeof(dfloat));
  auto o_Uouter = platform->device.malloc(nrs->NVfields * n * sizeof(dfloat));
  
  exactUVW(n, n, time, o_xinterp, o_yinterp, o_zinterp, o_UouterExact);

  auto interp = pointInterpolation_t(nrs);
  interp.setPoints(n, x.data(), y.data(), z.data());
  interp.find(pointInterpolation_t::VerbosityLevel::Detailed);
  interp.eval(nrs->NVfields, nrs->fieldOffset, nrs->o_U, n, o_Uouter);

  std::vector<dfloat> UouterExact(nrs->NVfields * n);
  o_UouterExact.copyTo(UouterExact.data(), nrs->NVfields * n * sizeof(dfloat));

  std::vector<dfloat> Uouter(nrs->NVfields * n);
  o_Uouter.copyTo(Uouter.data(), nrs->NVfields * n * sizeof(dfloat));

  // compare difference
  dfloat maxErrU = 0.0;
  dfloat maxErrV = 0.0;
  dfloat maxErrW = 0.0;
  for (int i = 0; i < n; i++) {
    maxErrU = std::max(maxErrU, fabs(UouterExact[0 * n + i] - Uouter[0 * n + i]));
    maxErrV = std::max(maxErrV, fabs(UouterExact[1 * n + i] - Uouter[1 * n + i]));
    maxErrW = std::max(maxErrW, fabs(UouterExact[2 * n + i] - Uouter[2 * n + i]));
  }

  auto tol = 1e-6;
  if(ciMode == 2) tol = 1e-3;
  bool passed = true;
  passed &= (maxErrU < tol);
  passed &= (maxErrV < tol);
  passed &= (maxErrW < tol);

  int allPassed = passed;
  MPI_Allreduce(MPI_IN_PLACE, &allPassed, 1, MPI_INT, MPI_MIN, platform->comm.mpiCommParent);

  if(!passed){
    if(platform->comm.mpiRank == 0){
      std::cout << "Point interpolation errors too large:\n";
      std::cout << "  maxErrU = " << maxErrU << "\n";
      std::cout << "  maxErrV = " << maxErrV << "\n";
      std::cout << "  maxErrW = " << maxErrW << "\n";
    }
    CIFAIL;
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  const int rank = platform->comm.mpiRank;
  
  if (!nrs->lastStep) {
    return;
  }
  
  ciTestPointInterp(nrs, time, tstep);

  auto *mesh = nrs->meshV;

  auto o_Uexact = platform->o_mempool.slice0;
  exactUVW(mesh->Nlocal, nrs->fieldOffset, time, mesh->o_x, mesh->o_y, mesh->o_z, o_Uexact);

  platform->linAlg->axpbyMany(mesh->Nlocal, nrs->NVfields, nrs->fieldOffset, 1.0, nrs->o_U, -1.0, o_Uexact);

  auto o_UerrX = o_Uexact + 0 * nrs->fieldOffset * sizeof(dfloat);
  auto o_UerrY = o_Uexact + 1 * nrs->fieldOffset * sizeof(dfloat);
  auto o_UerrZ = o_Uexact + 2 * nrs->fieldOffset * sizeof(dfloat);

  platform->linAlg->abs(mesh->Nlocal, o_UerrX);
  platform->linAlg->abs(mesh->Nlocal, o_UerrY);
  platform->linAlg->abs(mesh->Nlocal, o_UerrZ);

  const auto errVxGlobal = platform->linAlg->max(mesh->Nlocal, o_UerrX, platform->comm.mpiCommParent);
  const auto errVyGlobal = platform->linAlg->max(mesh->Nlocal, o_UerrY, platform->comm.mpiCommParent);
  const auto errVzGlobal = platform->linAlg->max(mesh->Nlocal, o_UerrZ, platform->comm.mpiCommParent);

  // scheme is first-order accurate
  dfloat expectedErrVx = 0.0;
  dfloat expectedErrVy = 0.0;
  dfloat expectedErrVz = 0.0;

  if (ciMode == 1 || ciMode == 3) {
    expectedErrVx = 5E-07;
    expectedErrVy = 5E-07;
    expectedErrVz = 2E-08;
  } else if (ciMode == 2) {
    expectedErrVx = 0.008;
    expectedErrVy = 0.007;
    expectedErrVz = 8E-05;
  }

  bool passed = true;
  passed &= errVxGlobal < 1.1 * expectedErrVx;
  passed &= errVyGlobal < 1.1 * expectedErrVy;
  passed &= errVzGlobal < 1.1 * expectedErrVz;

  int flag = passed;
  MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);

  if (flag) {
    (CIPASS);
  } else {
    if (platform->comm.mpiRank == 0) {
      printf("errVxGlobal = %g, expectedErrVx = %g, errVyGlobal = %g, expectedErrVy = %g, errVzGlobal = %g, "
             "expectedErrVz = %g\n",
             errVxGlobal,
             expectedErrVx,
             errVyGlobal,
             expectedErrVy,
             errVzGlobal,
             expectedErrVz);
    }
    (CIFAIL);
  }
}
