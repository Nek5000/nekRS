#include <math.h>

static int ciMode = 0;

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-3

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("7"));
  options.setArgs("RESTART FROM FILE", string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("END TIME", string("10"));
  options.setArgs("DT", string("2e-2"));
  options.setArgs("SUBCYCLING STEPS", string("0"));
  if (ciMode == 2) options.setArgs("SUBCYCLING STEPS", string("1"));
  options.setArgs("TIME INTEGRATOR", "TOMBO2");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-06"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-04"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", string("1e-06"));
  options.setArgs("VARIABLEPROPERTIES", "TRUE");
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  if (!nrs->lastStep) return;
 
  const int rank = platform->comm.mpiRank;

  nek::ocopyToNek(time, tstep);
  nek::userchk();

  double *err = (double *) nek::scPtr(1);
  double vxErr, sErr;
  vxErr = abs((err[0] - 2.06559)/err[0]);
  sErr  = abs((err[1] - 28.3833)/err[1]);

  if (rank == 0)
    printf("relative error to target: vx=%g s=%g\n", vxErr, sErr);

  (vxErr < EPS && sErr < EPS) ? (PASS) : (FAIL); 
}
