#include <math.h>

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-3

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("7"));
  options.setArgs("RESTART FROM FILE", string("0"));
  options.setArgs("VISCOSITY", string("0.01"));
  options.setArgs("FINAL TIME", string("1.0"));
  options.setArgs("DT", string("2e-4"));
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
  double vxErr, prErr;

  switch (ciMode) {
    case 1: vxErr = abs((err[0] - 1.018684E-04)/err[0]);
            prErr = abs((err[1] - 5.148630E-04)/err[1]);
            break;
    case 2: vxErr = abs((err[0] - 1.014353E-04)/err[0]);
            prErr = abs((err[1] - 5.127882E-04)/err[1]);
            break;
  }

  if (rank == 0)
    printf("relative error to target: vx=%g pr=%g\n", vxErr, prErr);

  (vxErr < EPS && prErr < EPS) ? (PASS) : (FAIL); 
}
