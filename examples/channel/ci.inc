#include <math.h>

static int ciMode = 0;

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("RESTART FROM FILE", std::string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", std::string("0.1"));
  options.setArgs("DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("0.1"));
  options.setArgs("DT", std::string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
  options.setArgs("BDF ORDER", "3");
  options.setArgs("ADVECTION", "FALSE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-08"));
  options.setArgs("VELOCITY STRESSFORMULATION", "TRUE");
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

  double *err = (double *) nek::scPtr(5);

  const double vxErr = abs((err[0] - 9.5942E-09)/err[0]);
  const double vyErr = abs((err[1] - 2.3012E-07)/err[1]);

  if (rank == 0)
    printf("relative error to target: vx=%g vy=%g\n", 
           vxErr, vyErr);

  (vxErr < EPS && vyErr < EPS) ? (CIPASS) : (CIFAIL); 
}
