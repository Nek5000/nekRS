#include <math.h>

static int ciMode = 0;

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("RESTART FROM FILE", std::string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", std::string("10.0"));
  options.setArgs("DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("0.3"));
  options.setArgs("DT", std::string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
  options.setArgs("TIME INTEGRATOR", "TOMBO3");
  options.setArgs("ADVECTION", "FALSE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-08"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-06"));
  options.setArgs("STRESSFORMULATION", "TRUE");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE");

  if(ciMode == 1){
    platform->par->set("casedata", "p_rot", 0.0);
  }
  if(ciMode == 2){
    platform->par->set("casedata", "p_rot", 0.7853981633974483);
  }
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  if (!nrs->lastStep) return;
 
  const int rank = platform->comm.mpiRank;
 
  nek::ocopyToNek(time, tstep);
  nek::userchk();

  double *err = (double *) nek::scPtr(1);

  double vxErr = 0.0;
  double vyErr = 0.0;
  vxErr = abs((err[0] - 1.1207E-08)/err[0]);
  vyErr = abs((err[1] - 1.8415E-11));

  const double absTol = 1e-8;

  if (rank == 0)
    printf("relative error to target: vx=%g vy=%g\n", 
           vxErr, vyErr);

  (vxErr < EPS && vyErr < absTol) ? (PASS) : (FAIL); 
}
