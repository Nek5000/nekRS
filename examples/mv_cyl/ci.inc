#include <math.h>
#include <iostream>

static int ciMode = 0;

#define EPS 5e-2

double numCVODERHS = 0;

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("NUMBER OF SCALARS", std::string("1"));
  options.setArgs("END TIME", std::string("0.1"));
  options.setArgs("DT", std::string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("PRESSURE INITIAL GUESS", "PREVIOUS");

  options.setArgs("BDF ORDER", "3");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-09"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-10"));
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));

  options.setArgs("MESH SOLVER", "NONE");
  options.setArgs("MOVING MESH", "TRUE");

  platform->par->set("casedata", "p_rot", 0.0);

  if (ciMode == 2) {
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
  }
  if (ciMode == 3) {
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("MESH INITIAL GUESS", "PROJECTION");
    options.setArgs("MESH RESIDUAL PROJECTION VECTORS", std::string("1"));
    options.setArgs("MESH SOLVER", "PCG");
    options.setArgs("MESH BLOCK SOLVER", "TRUE");
    options.setArgs("MESH SOLVER TOLERANCE", "1e-6");
  }
  if (ciMode == 4) {
    nekrsAbort(comm, EXIT_FAILURE, "%s\n", "ciMode=4 not defined");
  }
  if (ciMode == 5) {
    platform->par->set("casedata", "p_rot", 0.5);
  }
  if (ciMode == 6 || ciMode == 8) {
    options.setArgs("MESH SOLVER", "PCG");
    options.setArgs("MESH BLOCK SOLVER", "TRUE");
    options.setArgs("MESH SOLVER TOLERANCE", "1e-6");
    options.setArgs("MESH INITIAL GUESS", "PROJECTION");
    options.setArgs("MESH RESIDUAL PROJECTION VECTORS", std::string("1"));
    platform->par->set("casedata", "p_rot", 0.5);
  }
  if (ciMode == 7 || ciMode == 8) {
    // cvode
    options.setArgs("CVODE", "TRUE");

    double dt0;
    options.getArgs("DT", dt0);

    options.setArgs("SCALAR00 CVODE ABSOLUTE TOLERANCE", "1e-08");
    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-06");
    options.setArgs("SCALAR00 SOLVER", "CVODE");
  }
}

void ciTestErrors(nrs_t *nrs, double time, int tstep)
{
  if (nrs->cds->cvode) {
    numCVODERHS += nrs->cds->cvode->numRHSEvals();
    numCVODERHS += nrs->cds->cvode->numLinIters();
  }

  if (!nrs->lastStep) {
    return;
  }

  numCVODERHS /= tstep;

  const int rank = platform->comm.mpiRank;
  mesh_t *mesh = nrs->mesh;

  nrs->copyToNek(time, tstep);
  nek::userchk();

  auto err = nek::ptr<double>("errors");
  const int Nfields = 6;
  dfloat expectedErr[Nfields];

  expectedErr[0] = 0.2465620E-06;
  expectedErr[1] = 0.6938894E-16;
  expectedErr[2] = 0.2973762E-05;
  expectedErr[3] = 0.5820622E-06;
  expectedErr[4] = 0.2465620E-05;
  expectedErr[5] = 0.4929804E-06;

  int expectedSIter = 9;
  int expectedPIter = 9;
  int expectedVIter = 8;

  if (ciMode == 1) {
  }
  if (ciMode == 2) {
    expectedErr[5] = 0.3E-05;
  }
  if (ciMode == 3) {
    expectedErr[5] = 0.3E-05;
  }
  if (ciMode == 5) {
    // expectedVIter = 9; // no extrap
  }
  if (ciMode == 6) {
    // expectedVIter = 8; // no extrap
  }
  if (ciMode == 7 || ciMode == 8) {
    expectedErr[5] = 0.35E-05;
  }

  auto errVal = abs(err[0] - expectedErr[0]) / expectedErr[0];
  CiEvalTest("vol (LinfErr: " + to_string_f(errVal) + ")", errVal < EPS);

  errVal = abs(err[2] - expectedErr[2]) / expectedErr[2];
  CiEvalTest("p0th (LinfErr: " + to_string_f(errVal) + ")", errVal < EPS);

  errVal = abs(err[3] - expectedErr[3]) / expectedErr[3];
  CiEvalTest("dpdt (LinfErr: " + to_string_f(errVal) + ")", errVal < EPS);

  errVal = abs(err[4] - expectedErr[4]) / expectedErr[4];
  CiEvalTest("yPison (LinfErr: " + to_string_f(errVal) + ")", errVal < EPS);

  errVal = abs(err[5] - expectedErr[5]) / expectedErr[5];
  CiEvalTest("T (LinfErr: " + to_string_f(errVal) + ")", errVal < EPS);

  errVal = abs(expectedVIter - nrs->uvwSolver->Niter());
  CiEvalTest("velIter (LinfErr: " + to_string_f(errVal) + ")", errVal < 2);

  errVal = abs(expectedPIter - nrs->pSolver->Niter());
  CiEvalTest("pIter (LinfErr: " + to_string_f(errVal) + ")", errVal < 2);
}
