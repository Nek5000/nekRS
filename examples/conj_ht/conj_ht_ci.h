#include <math.h>

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-3

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("7"));
  options.setArgs("RESTART FROM FILE", string("0"));
  options.setArgs("TSTEPS FOR SOLUTION OUTPUT", "0");
  options.setArgs("FINAL TIME", string("10"));
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

void ciTestErrors(ins_t *ins, dfloat time, int tstep)
{
  if (tstep != ins->NtimeSteps) return;
 
  const int rank = ins->mesh->rank;

  nek_ocopyFrom(time, tstep);
  nek_userchk();

  double *norm = nekData.cbscnrs;
  double vxErr, sErr;

  switch (ciMode) {
    // cross compare solution to nek5000
    case 1: vxErr = abs((norm[0] - 2.06559)/norm[0]);
            sErr  = abs((norm[1] - 28.3833)/norm[1]);
            break;
  }

  if (rank == 0)
    printf("relative error to target: vx=%g s=%g\n", vxErr, sErr);

  (vxErr < EPS && sErr < EPS) ? (PASS) : (FAIL); 
}
