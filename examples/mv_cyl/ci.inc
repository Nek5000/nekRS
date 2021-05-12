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
  options.setArgs("NUMBER OF SCALARS", string("1"));
  options.setArgs("END TIME", string("0.1"));
  options.setArgs("DT", string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", string("0"));
  options.setArgs("PRESSURE RESIDUAL PROJECTION", "FALSE");

  options.setArgs("TIME INTEGRATOR", "TOMBO3");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-09"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-10"));
  options.setArgs("SCALAR00 SOLVER TOLERANCE", string("1e-12"));
  
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
  const int Nfields = 6;
  dfloat expectedErr[Nfields];
  int expectedPIter;
  if(ciMode == 1)
  {
    expectedErr[0] = 0.2465620E-06;
    expectedErr[1] = 0.6938894E-16;
    expectedErr[2] = 0.3044039E-05;
    expectedErr[3] = 0.2144404E-05;
    expectedErr[4] = 0.2465620E-05;
    expectedErr[5] = 0.1080412E-05;
    expectedPIter = 6;
  }
  if(ciMode == 2)
  {
    expectedErr[0] = 0.2465620E-06;
    expectedErr[1] = 0.6938894E-16;
    expectedErr[2] = 0.3044039E-05;
    expectedErr[3] = 0.2144404E-05;
    expectedErr[4] = 0.2465620E-05;
    expectedErr[5] = 0.2443839E-05;
    expectedPIter = 6;
  }

  const dfloat floor = 1e-12; // values below this threshold are ignored
  const int expectedVIter = 3;

  double vxErr, prErr;
  double s1Err;
  
  int pIterErr;
  int velIterErr;

  bool passed = true;

  if(rank == 0) printf("relative error to target: ");
  for(int i = 0 ; i < Nfields; ++i){
    if(abs(err[i]) > floor){
      const dfloat relErr = abs(err[i]-expectedErr[i])/expectedErr[i];
      if(rank == 0) printf("err[%d]=%g ", i, relErr);
      if(relErr > EPS) passed = false;
    }
  }
  velIterErr = abs(expectedVIter - nrs->uvwSolver->Niter);
  pIterErr = abs(expectedPIter - nrs->pSolver->Niter);
  if(rank == 0) printf("velIter=%d pIter=%d\n", velIterErr, pIterErr);
  if(velIterErr >= 2) passed = false;
  if(pIterErr >= 2) passed = false;

  passed ? (PASS) : (FAIL);
}
