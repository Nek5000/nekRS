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
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-08"));

  if (ciMode == 2) {
    std::string casename = platform->options.getArgs("CASENAME");
    int fail = 0;
    if (casename == "inside") {
      fail = 1;
      if (platform->comm.mpiRank == 0) {
        std::cout << "Unexpected input file for CI mode" << std::endl;
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &fail, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);
    if (fail) {
      CIFAIL;
    }
  }
}

void ciTestPartitionOfUnity(nrs_t *nrs, double time, int tstep)
{
  auto mesh = nrs->meshV;
  auto minMaxFld = [&](occa::memory o_a) {
    const auto minFld = platform->linAlg->min(mesh->Nlocal, o_a, platform->comm.mpiCommParent);
    const auto maxFld = platform->linAlg->max(mesh->Nlocal, o_a, platform->comm.mpiCommParent);
    return std::make_tuple(minFld, maxFld);
  };

  const auto [minX, maxX] = minMaxFld(mesh->o_x);
  const auto [minY, maxY] = minMaxFld(mesh->o_y);
  const auto [minZ, maxZ] = minMaxFld(mesh->o_z);
  const auto expectedGlobalVolume = (maxX - minX) * (maxY - minY) * (maxZ - minZ);

  const auto o_partition = nrs->neknek->partitionOfUnity();
  auto o_volume = platform->o_memPool.reserve<dfloat>(mesh->Nlocal);

  platform->linAlg->axmyz(mesh->Nlocal, 1.0, o_partition, mesh->o_LMM, o_volume);

  const auto globalVolume = platform->linAlg->sum(mesh->Nlocal, o_volume, platform->comm.mpiCommParent);

  o_volume.free();

  const auto relErr = abs(globalVolume - expectedGlobalVolume) / expectedGlobalVolume;

  const auto tol = 5e-6;
  if (relErr > tol) {
    if (platform->comm.mpiRank == 0) {
      std::cout << "Partition of unity test failed: " << relErr << std::endl;
      std::cout << "Expected value : " << expectedGlobalVolume << std::endl;
      std::cout << "Actual   value : " << globalVolume << std::endl;
    }
    CIFAIL;
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }
}

void ciTestPointInterp(nrs_t *nrs, double time, int tstep)
{
  // query velocity solution from inner session
  std::string casename = platform->options.getArgs("CASENAME");
  std::vector<dfloat> x;
  std::vector<dfloat> y;
  std::vector<dfloat> z;
  if (casename == "outside") {
    x = {M_PI};
    y = {M_PI};
    z = {0.05};
  } else {
    x = {0.25 * M_PI, M_PI, M_PI, 1.75 * M_PI};
    y = {0.25 * M_PI, M_PI, M_PI, 1.75 * M_PI};
    z = {0.05, 0.05, 0.05, 0.05};
  }
  int n = x.size();

  auto o_xinterp = platform->device.malloc<dfloat>(n, x.data());
  auto o_yinterp = platform->device.malloc<dfloat>(n, y.data());
  auto o_zinterp = platform->device.malloc<dfloat>(n, z.data());
  auto o_UouterExact = platform->device.malloc<dfloat>(nrs->NVfields * n);
  auto o_Uouter = platform->device.malloc<dfloat>(nrs->NVfields * n);

  exactUVW(n, n, time, o_xinterp, o_yinterp, o_zinterp, o_UouterExact);

  auto interp = pointInterpolation_t(nrs);
  interp.setPoints(n, o_xinterp, o_yinterp, o_zinterp);
  interp.find(pointInterpolation_t::VerbosityLevel::Detailed);
  interp.eval(nrs->NVfields, nrs->fieldOffset, nrs->o_U, n, o_Uouter);

  std::vector<dfloat> UouterExact(nrs->NVfields * n);
  o_UouterExact.copyTo(UouterExact.data(), nrs->NVfields * n);

  std::vector<dfloat> Uouter(nrs->NVfields * n);
  o_Uouter.copyTo(Uouter.data(), nrs->NVfields * n);

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
  if (ciMode == 2) {
    tol = 1e-3;
  }
  bool passed = true;
  passed &= (maxErrU < tol);
  passed &= (maxErrV < tol);
  passed &= (maxErrW < tol);

  int allPassed = passed;
  MPI_Allreduce(MPI_IN_PLACE, &allPassed, 1, MPI_INT, MPI_MIN, platform->comm.mpiCommParent);

  if (!passed) {
    if (platform->comm.mpiRank == 0) {
      std::cout << "Point interpolation errors too large:\n";
      std::cout << "  maxErrU = " << maxErrU << "\n";
      std::cout << "  maxErrV = " << maxErrV << "\n";
      std::cout << "  maxErrW = " << maxErrW << "\n";
    }
    CIFAIL;
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }

  o_xinterp.free();
  o_yinterp.free();
  o_zinterp.free();
  o_UouterExact.free();
  o_Uouter.free();
}

void vizNekNekField(nrs_t *nrs)
{
  auto mesh = nrs->meshV;
  auto neknek = nrs->neknek;

  occa::memory o_null;

  const auto neknekFieldOffset = neknek->fieldOffset();

  std::vector<dfloat> Uviz(nrs->NVfields * nrs->fieldOffset, 0.0);
  std::vector<dfloat> neknekU(nrs->NVfields * neknekFieldOffset, 0.0);
  std::fill(Uviz.begin(), Uviz.end(), 0.0);
  std::fill(neknekU.begin(), neknekU.end(), 0.0);
  neknek->o_U().copyTo(neknekU.data(), nrs->NVfields * neknekFieldOffset);

  std::vector<dlong> neknekPointMap(mesh->Nlocal, -1);
  neknek->o_pointMap().copyTo(neknekPointMap.data(), mesh->Nlocal);

  for (int n = 0; n < mesh->Nlocal; ++n) {
    auto ip = neknekPointMap[n];
    if (ip != -1) {
      for (int fld = 0; fld < nrs->NVfields; ++fld) {
        Uviz[fld * nrs->fieldOffset + n] = neknekU[fld * neknekFieldOffset + ip];
      }
    } else {
      // not a neknek point, set to 0
      for (int fld = 0; fld < nrs->NVfields; ++fld) {
        Uviz[fld * nrs->fieldOffset + n] = 0.0;
      }
    }
  }

  auto o_Uviz = platform->device.malloc<dfloat>(nrs->NVfields * nrs->fieldOffset, Uviz.data());
  writeFld("viz", 0.0, 0, 1, 1, o_Uviz, o_null, o_null, 0);
}

void ciTestInitialStep(nrs_t *nrs, double time, int tstep)
{
  auto mesh = nrs->meshV;
  auto neknek = nrs->neknek;
  const auto neknekFieldOffset = neknek->fieldOffset();
  exactUVW(mesh->Nlocal, nrs->fieldOffset, time, mesh->o_x, mesh->o_y, mesh->o_z, nrs->o_U);
  
  neknek->updateBoundary(tstep, 2); // pass in stage == 2 to prevent lagging

  auto o_Uexp = platform->o_memPool.reserve<dfloat>(nrs->NVfields * neknekFieldOffset);
  auto npt = neknek->npt();
  exactUVW(npt, neknekFieldOffset, time, neknek->o_x(), neknek->o_y(), neknek->o_z(), o_Uexp);

  // compare errors
  platform->linAlg->axpbyMany(npt, nrs->NVfields, neknekFieldOffset, 1.0, neknek->o_U(), -1.0, o_Uexp);
  platform->linAlg->aydxMany(npt, nrs->NVfields, neknekFieldOffset, 1, 1.0, neknek->o_U(), o_Uexp);
  platform->linAlg->abs(nrs->NVfields * neknekFieldOffset, o_Uexp);
  auto maxRelErr = platform->linAlg->max(nrs->NVfields * neknekFieldOffset, o_Uexp, platform->comm.mpiComm);
  auto maxRelErrGlob =
      platform->linAlg->max(nrs->NVfields * neknekFieldOffset, o_Uexp, platform->comm.mpiCommParent);
  if (platform->comm.mpiRank == 0) {
    std::cout << "NekNek interpolation errors:\n";
    std::cout << "  maxRelErr = " << maxRelErr << "\n";
    std::cout << "  maxRelErrGlob = " << maxRelErrGlob << "\n";
  }

  const auto tol = 1e-5;
  bool passed = (maxRelErrGlob < tol);
  o_Uexp.free();

  if (!passed) {
    vizNekNekField(nrs);
    CIFAIL;
    nrsFinalize(nrs);
    exit(platform->exitValue);
  }
}

void ciTestErrors(nrs_t *nrs, double time, int tstep)
{
  if (tstep == 0) {
    ciTestInitialStep(nrs, time, tstep);
    ciTestPointInterp(nrs, time, tstep);
    ciTestPartitionOfUnity(nrs, time, tstep);
  }

  const int rank = platform->comm.mpiRank;

  if (!nrs->lastStep) {
    return;
  }

  auto *mesh = nrs->meshV;

  auto o_Uexact = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
  auto o_Sexact = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
  exactUVW(mesh->Nlocal, nrs->fieldOffset, time, mesh->o_x, mesh->o_y, mesh->o_z, o_Uexact);
  o_Sexact.copyFrom(o_Uexact, nrs->fieldOffset);

  platform->linAlg->axpbyMany(mesh->Nlocal, nrs->NVfields, nrs->fieldOffset, 1.0, nrs->o_U, -1.0, o_Uexact);
  platform->linAlg->axpby(mesh->Nlocal, 1.0, nrs->cds->o_S, -1.0, o_Sexact);

  auto o_UerrX = o_Uexact + 0 * nrs->fieldOffset;
  auto o_UerrY = o_Uexact + 1 * nrs->fieldOffset;
  auto o_UerrZ = o_Uexact + 2 * nrs->fieldOffset;

  platform->linAlg->abs(mesh->Nlocal, o_UerrX);
  platform->linAlg->abs(mesh->Nlocal, o_UerrY);
  platform->linAlg->abs(mesh->Nlocal, o_UerrZ);
  platform->linAlg->abs(mesh->Nlocal, o_Sexact);

  const auto errVxGlobal = platform->linAlg->max(mesh->Nlocal, o_UerrX, platform->comm.mpiCommParent);
  const auto errVyGlobal = platform->linAlg->max(mesh->Nlocal, o_UerrY, platform->comm.mpiCommParent);
  const auto errVzGlobal = platform->linAlg->max(mesh->Nlocal, o_UerrZ, platform->comm.mpiCommParent);
  const auto errSGlobal = platform->linAlg->max(mesh->Nlocal, o_Sexact, platform->comm.mpiCommParent);

  o_Uexact.free();
  o_Sexact.free();

  // scheme is first-order accurate
  dfloat expectedErrVx = 0.0;
  dfloat expectedErrVy = 0.0;
  dfloat expectedErrVz = 0.0;
  dfloat expectedErrS = 0.0;

  if (ciMode == 1 || ciMode == 3) {
    expectedErrVx = 5E-07;
    expectedErrVy = 5E-07;
    expectedErrVz = 3E-08;
    expectedErrS = 5E-07;
  } else if (ciMode == 2) {
    expectedErrVx = 0.008;
    expectedErrVy = 0.007;
    expectedErrVz = 8E-05;
    expectedErrS = 0.008;
  }

  bool passed = true;
  passed &= errVxGlobal < 1.1 * expectedErrVx;
  passed &= errVyGlobal < 1.1 * expectedErrVy;
  passed &= errVzGlobal < 1.1 * expectedErrVz;
  passed &= errSGlobal < 1.1 * expectedErrS;

  int flag = passed;
  MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MAX, platform->comm.mpiCommParent);

  if (flag) {
    (CIPASS);
  } else {
    if (platform->comm.mpiRank == 0) {
      printf("errVxGlobal = %g, expectedErrVx = %g, errVyGlobal = %g, expectedErrVy = %g, errVzGlobal = %g, "
             "expectedErrVz = %g, errSGlobal = %g, expectedErrS = %g\n",
             errVxGlobal,
             expectedErrVx,
             errVyGlobal,
             expectedErrVy,
             errVzGlobal,
             expectedErrVz,
             errSGlobal,
             expectedErrS);
    }
    (CIFAIL);
  }
}
