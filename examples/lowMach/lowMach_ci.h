#include <math.h>

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-2

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("7"));
  options.setArgs("RESTART FROM FILE", string("0"));
  options.setArgs("TSTEPS FOR SOLUTION OUTPUT", "0");
  options.setArgs("FINAL TIME", string("0.3"));
  options.setArgs("DT", string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", string("0"));
  if (ciMode == 2) options.setArgs("SUBCYCLING STEPS", string("1"));
  options.setArgs("TIME INTEGRATOR", "TOMBO2");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-08"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("VARIABLEPROPERTIES", "TRUE");
}

void ciTestErrors(ins_t *ins, dfloat time, int tstep)
{
  if (tstep != ins->NtimeSteps) return;
 
  const int rank = ins->mesh->rank;

  ins->o_div.copyTo(ins->div);
  dlong Nlocal = ins->mesh->Nelements * ins->mesh->Np;
  memcpy(nekData.qtl, ins->div, sizeof(dfloat)*Nlocal);
 
  nek_ocopyFrom(time, tstep);
  nek_userchk();

  double *err = nekData.cbscnrs;
  double vxErr, prErr, sErr;

  switch (ciMode) {
    case 1: vxErr = abs((err[0] - 1.59E-05)/err[0]);
            prErr = abs((err[1] - 3.19E-03)/err[1]);
            sErr  = abs((err[2] - 3.07E-07)/err[2]);
            break;
    case 2: vxErr = abs((err[0] - 2.19E-05)/err[0]);
            prErr = abs((err[1] - 3.22E-03)/err[1]);
            sErr  = abs((err[2] - 1.05E-06)/err[2]);
            break;
  }

  if (rank == 0)
    printf("relative error to target: vx=%g pr=%g s=%g\n", vxErr, prErr, sErr);

  (vxErr < EPS && prErr < EPS && sErr < EPS) ? (PASS) : (FAIL); 
}
