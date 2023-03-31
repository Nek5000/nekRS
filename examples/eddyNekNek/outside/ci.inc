#include <math.h>
#include <limits>
#include <vector>
#include <numeric>

static int ciMode = 0;

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("RESTART FROM FILE", std::string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", std::string("0.05"));
  options.setArgs("DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("3e-3"));
  options.setArgs("DT", std::string("1e-4"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("PRESSURE INITIAL GUESS", "PREVIOUS");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");

  options.setArgs("BDF ORDER", "3");
  if (ciMode == 1) {
    options.setArgs("NEKNEK CORRECTOR STEPS", "1");
    options.setArgs("NEKNEK BOUNDARY EXT ORDER", "3");
  }
  else if (ciMode == 2) {
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

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  if (!nrs->lastStep)
    return;
  const int rank = platform->comm.mpiRank;

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

  if (ciMode == 1) {
    expectedErrVx = 5E-07;
    expectedErrVy = 5E-07;
    expectedErrVz = 7E-09;
  }
  else if (ciMode == 2) {
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
  }
  else {
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
