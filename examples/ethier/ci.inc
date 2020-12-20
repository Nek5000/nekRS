#include <math.h>

static int ciMode = 0;

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
  options.setArgs("END TIME", string("0.06"));
  options.setArgs("DT", string("2e-3"));
  options.setArgs("SUBCYCLING STEPS", string("0"));
  options.setArgs("PRESSURE RESIDUAL PROJECTION", "FALSE");

  options.setArgs("VELOCITY BLOCK SOLVER", "FALSE");
  options.setArgs("SCALAR INITIAL GUESS DEFAULT","PREVIOUS");
  options.setArgs("VELOCITY INITIAL GUESS DEFAULT","PREVIOUS");

  if (ciMode == 2) {
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", string("1"));
    options.setArgs("PRESSURE RESIDUAL PROJECTION", "TRUE");
  }
  if (ciMode == 3) {
    options.setArgs("PRESSURE RESIDUAL PROJECTION", "TRUE");
    options.setArgs("VELOCITY RESIDUAL PROJECTION", "TRUE");
  }
  if (ciMode == 4) {
    options.setArgs("END TIME", string("0.2"));
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", string("1"));
    options.setArgs("PRESSURE RESIDUAL PROJECTION", "TRUE");
    options.setArgs("VELOCITY RESIDUAL PROJECTION", "TRUE");
  }
  options.setArgs("TIME INTEGRATOR", "TOMBO3");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-08"));
  options.setArgs("SCALAR00 SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("VARIABLEPROPERTIES", "FALSE");
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  if (!nrs->lastStep) return;
 
  const int rank = nrs->mesh->rank;
 
  nek_ocopyFrom(time, tstep);
  nek_userchk();

  double *err = (double *) nek_scPtr(1);

  double vxErr, prErr;
  double s1Err, s2Err;
  
  int pIterErr;
  int velIterErr;

  switch (ciMode) {
    case 1 : velIterErr = abs(nrs->NiterU - 10);
             s1Err = abs((err[2] - 5.25E-12)/err[2]);
             s2Err = abs((err[3] - 6.09E-12)/err[3]);
             pIterErr = abs(nrs->NiterP - 4);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 7.17E-10)/err[1]);
             break;
    case 2 : velIterErr = abs(nrs->NiterU - 10);
             s1Err = abs((err[2] - 6.49E-12)/err[2]);
             s2Err = abs((err[3] - 7.28E-12)/err[3]);
             pIterErr = abs(nrs->NiterP - 1);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 8.38E-10)/err[1]);
             break;
    case 3 : velIterErr = abs(nrs->NiterU - 4);
             s1Err = abs((err[2] - 5.25E-12)/err[2]);
             s2Err = abs((err[3] - 6.09E-12)/err[3]);
             pIterErr = abs(nrs->NiterP - 1);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 8.34E-10)/err[1]);
             break;
    case 4 : velIterErr = abs(nrs->NiterU - 1);
             s1Err = abs((err[2] - 1.71E-11)/err[2]);
             s2Err = abs((err[3] - 2.00E-11)/err[3]);
             pIterErr = abs(nrs->NiterP - 1);
             vxErr = abs((err[0] - 3.65E-10)/err[0]);
             prErr = abs((err[1] - 6.71E-10)/err[1]);
             break;
  }

  if (rank == 0)
    printf("relative error to target: vx=%g pr=%g s1=%g s2=%g velIter=%d pIter=%d\n", 
           vxErr, prErr, s1Err, s2Err, velIterErr, pIterErr);


  (vxErr < EPS && prErr < EPS && s1Err < EPS && s2Err < EPS && 
  velIterErr <= 1 && pIterErr <= 2) ? (PASS) : (FAIL); 
}
