#include <math.h>
#include <iostream>

static int ciMode = 0;

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("7"));
  options.setArgs("RESTART FROM FILE", string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("NUMBER OF SCALARS", string("0"));
  options.setArgs("NUMBER TIMESTEPS", "30");
  options.setArgs("DT", string("16e-4"));
  options.setArgs("SUBCYCLING STEPS", string("0"));
  options.setArgs("PRESSURE RESIDUAL PROJECTION", "FALSE");

  options.setArgs("TIME INTEGRATOR", "TOMBO2");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-09"));
  
  if(ciMode == 2){
    options.setArgs("SUBCYCLING STEPS", string("1"));
  }
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  if (!nrs->lastStep) return;
 
  const int rank = platform->comm.mpiRank;
 
  nek::ocopyToNek(time, tstep);
  nek::userchk();

  double *err = (double *) nek::scPtr(1);
  const int Nfields = 2;
  dfloat expectedErr[Nfields];
  if(ciMode == 1)
  {
    expectedErr[0] = 1.585353E-03;
    expectedErr[1] = 1.806126E-03;
  }
  if(ciMode == 2)
  {
    expectedErr[0] = 2.834266E-03;
    expectedErr[1] = 4.305412E-03;
  }

  const dfloat floor = 1e-12; // values below this threshold are ignored
  const int expectedUIter = 10;
  const int expectedVIter = 10;
  const int expectedWIter = 2;
  const int expectedPIter = 14;

  bool passed = true;

  if(rank == 0) printf("relative error to target: ");
  for(int i = 0 ; i < Nfields; ++i){
    if(abs(err[i]) > floor){
      const dfloat relErr = abs(err[i]-expectedErr[i])/expectedErr[i];
      if(rank == 0) printf("err[%d]=%g ", i, relErr);
      if(relErr > EPS) passed = false;
    }
  }
  int uIterErr = abs(expectedUIter - nrs->uSolver->Niter);
  int vIterErr = abs(expectedVIter - nrs->vSolver->Niter);
  int wIterErr = abs(expectedWIter - nrs->wSolver->Niter);
  int pIterErr = abs(expectedPIter - nrs->pSolver->Niter);
  if(rank == 0) printf("uIter=%d, vIter=%d, wIter=%d, pIter=%d\n", uIterErr, vIterErr, wIterErr, pIterErr);
  if(uIterErr >= 2) passed = false;
  if(vIterErr >= 2) passed = false;
  if(wIterErr >= 2) passed = false;
  if(pIterErr >= 2) passed = false;

  passed ? (PASS) : (FAIL);
}
