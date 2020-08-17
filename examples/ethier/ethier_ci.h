#include <math.h>

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("9"));
  options.setArgs("RESTART FROM FILE", string("0"));
  options.setArgs("TSTEPS FOR SOLUTION OUTPUT", "0");
  options.setArgs("VISCOSITY", string("0.01"));
  options.setArgs("DENSITY", string("1"));
  options.setArgs("NUMBER OF SCALARS", string("2"));
  options.setArgs("SCALAR00 DIFFUSIVITY", string("0.01"));
  options.setArgs("SCALAR00 DENSITY", string("1"));
  options.setArgs("SCALAR01 DIFFUSIVITY", string("0.01"));
  options.setArgs("SCALAR01 DENSITY", string("1"));
  options.setArgs("FINAL TIME", string("0.2"));
  options.setArgs("DT", string("2e-3"));
  options.setArgs("SUBCYCLING STEPS", string("0"));
  if (ciMode == 2) {
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", string("1"));
  }
  options.setArgs("TIME INTEGRATOR", "TOMBO3");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-08"));
  options.setArgs("PRESSURE RESIDUAL PROJECTION", "FALSE");
  options.setArgs("SCALAR00 SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("VARIABLEPROPERTIES", "FALSE");
}

void ciTestErrors(ins_t *ins, dfloat time, int tstep)
{
  if (tstep != ins->NtimeSteps) return;
 
  const int rank = ins->mesh->rank;
 
  nek_ocopyFrom(time, tstep);
  nek_userchk();

  double *err = (double *) nek_scPtr(1);

  const double vxErr = abs((err[0] - 2.46E-10)/err[0]);
  const double prErr = abs((err[1] - 7.64E-10)/err[1]);
  double s1Err, s2Err;
  
  const int pIterErr = abs(ins->NiterP - 4);
  int velIterErr;

  switch (ciMode) {
    case 1 : velIterErr = abs(ins->NiterU - 10);
             s1Err = abs((err[2] - 9.75E-12)/err[2]);
             s2Err = abs((err[3] - 9.75E-12)/err[3]);
             break;
    case 2 : velIterErr = abs(ins->NiterU - 10);
             s1Err = abs((err[2] - 1.69E-11)/err[2]);
             s2Err = abs((err[3] - 1.69E-11)/err[3]);
             break;
  }

  if (rank == 0)
    printf("relative error to target: vx=%g pr=%g s1=%g s2=%g velIter=%d pIter=%d\n", 
           vxErr, prErr, s1Err, s2Err, velIterErr, pIterErr);


  (vxErr < EPS && prErr < EPS && s1Err < EPS && s2Err < EPS && 
  velIterErr <= 1 && pIterErr <= 2) ? (PASS) : (FAIL); 
}
