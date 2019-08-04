#include <math.h>

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-3

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("7"));
  options.setArgs("VISCOSITY", string("0.05"));
  options.setArgs("FINAL TIME", string("0.1"));
  options.setArgs("DT", string("1e-4"));
  options.setArgs("SUBCYCLING STEPS", string("0"));
  if (ciMode == 2) options.setArgs("SUBCYCLING STEPS", string("1"));
  options.setArgs("TIME INTEGRATOR", "TOMBO2");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-09"));
}

void ciTestErrors(ins_t *ins, dfloat time, int tstep)
{
  if (tstep != ins->NtimeSteps) return;
 
  const int rank = ins->mesh->rank;
 
  nek_ocopyFrom(ins, time, tstep);
  nek_userchk();

  double *err = nekData.cbscnrs;
  double vxErr, vyErr, prErr;

  switch (ciMode) {
    case 1: vxErr = abs((err[0] - 6.079791E-07)/err[0]);
            vyErr = abs((err[1] - 7.650646E-07)/err[1]);
            prErr = abs((err[2] - 1.935742E-05)/err[2]);
            break;
    case 2: vxErr = abs((err[0] - 6.582217E-07)/err[0]);
            vyErr = abs((err[1] - 7.860122E-07)/err[1]);
            prErr = abs((err[2] - 2.012561E-05)/err[2]);
            break;
  }

  if (rank == 0)
    printf("relative error to target: vx=%g vy=%g pr=%g\n", vxErr, vyErr, prErr);

  (vxErr < EPS && vyErr < EPS && prErr < EPS) ? (PASS) : (FAIL); 
}
