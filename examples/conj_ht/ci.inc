#include <math.h>

static int ciMode = 0;

#define EPS 1e-3

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("RESTART FROM FILE", std::string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("END TIME", std::string("10"));
  options.setArgs("DT", std::string("2e-2"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  if (ciMode == 2) options.setArgs("SUBCYCLING STEPS", std::string("1"));
  options.setArgs("BDF ORDER", "2");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-06"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-04"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-06"));
  options.setArgs("VARIABLEPROPERTIES", "TRUE");

  // set up CVODE solver details
  if(ciMode == 3){
    options.setArgs("CVODE", "TRUE");
    options.setArgs("CVODE INTEGRATOR", "BDF");
    options.setArgs("CVODE MAX TIMESTEPPER ORDER", "3");
    options.setArgs("CVODE EPS LIN", "0.05");
    options.setArgs("CVODE STOP TIME", "TRUE");

    const auto hmax = 1;
    options.setArgs("CVODE HMAX", to_string_f(hmax));

    const auto maxSteps = 10000;
    options.setArgs("CVODE MAX STEPS", std::to_string(maxSteps));

    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-06");
    options.setArgs("CVODE ABSOLUTE TOLERANCE", "1e-08");
    options.setArgs("CVODE GMR VECTORS", "20");

    options.setArgs("SCALAR00 SOLVER", "CVODE");

  }
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

  (vxErr < EPS && sErr < EPS) ? (CIPASS) : (CIFAIL); 
}
