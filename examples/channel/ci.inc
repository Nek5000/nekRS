#include <math.h>

static int ciMode = 0;

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("7"));
  options.setArgs("RESTART FROM FILE", string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", string("0.1"));
  options.setArgs("DENSITY", string("1"));
  options.setArgs("END TIME", string("0.1"));
  options.setArgs("DT", string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", string("0"));
  options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
  options.setArgs("TIME INTEGRATOR", "TOMBO3");
  options.setArgs("ADVECTION", "FALSE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-08"));
  options.setArgs("STRESSFORMULATION", "TRUE");
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  if (!nrs->lastStep) return;
 
  const int rank = platform->comm.mpiRank;
 
  nek::ocopyToNek(time, tstep);
  nek::userchk();

  double *err = (double *) nek::scPtr(5);

  double vxErr = abs((err[0] - 9.5942E-09)/err[0]);
  double vyErr = abs((err[1] - 2.3012E-07)/err[1]);

  if (rank == 0)
    printf("relative error to target: vx=%g vy=%g\n", 
           vxErr, vyErr);

  (vxErr < EPS && vyErr < EPS) ? (PASS) : (FAIL); 
}
