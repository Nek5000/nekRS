#include <math.h>
#include <iostream>

static int ciMode = 0;

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("RESTART FROM FILE", std::string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("NUMBER OF SCALARS", std::string("1"));
  options.setArgs("END TIME", std::string("0.1"));
  options.setArgs("DT", std::string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("PRESSURE INITIAL GUESS", "PREVIOUS");

  options.setArgs("TIME INTEGRATOR", "TOMBO3");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-09"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-10"));
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("MESH SOLVER", "USER");

  platform->par->set("casedata", "p_rot", 0.0);
  
  if(ciMode == 2){
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
  }
  if(ciMode == 3){
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("MESH SOLVER", "ELASTICITY");
    options.setArgs("MESH SOLVER TOLERANCE", "1e-6");
  }
  if(ciMode == 4){
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("MESH INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("MESH SOLVER", "ELASTICITY");
    options.setArgs("MESH SOLVER TOLERANCE", "1e-6");
  }
  if(ciMode == 5){
    platform->par->set("casedata", "p_rot", 0.5);
  }
  if(ciMode == 6){
    options.setArgs("MESH SOLVER", "ELASTICITY");
    options.setArgs("MESH SOLVER TOLERANCE", "1e-12");
    platform->par->set("casedata", "p_rot", 0.5);
  }
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  if (!nrs->lastStep) return;
 
  const int rank = platform->comm.mpiRank;
  mesh_t* mesh = nrs->meshV;
 
  nek::ocopyToNek(time, tstep);
  nek::userchk();

  double *err = (double *) nek::scPtr(1);
  const int Nfields = 6;
  dfloat expectedErr[Nfields];

  expectedErr[0] = 0.2465620E-06;
  expectedErr[1] = 0.6938894E-16;
  expectedErr[2] = 0.2973762E-05;
  expectedErr[3] = 0.5820565E-06;
  expectedErr[4] = 0.2465620E-05;
  expectedErr[5] = 0.2558128E-05;

  int expectedPIter = 5;
  if(ciMode == 1)
  {
  }
  if(ciMode == 2)
  {
    expectedErr[5] = 0.1042903E-05;
  }
  if(ciMode == 3)
  {
    expectedErr[5] = 0.1042903E-05;
  }
  if(ciMode == 4)
  {
    expectedErr[5] = 0.1042903E-05;
  }
  if(ciMode == 5)
  {
  }
  if(ciMode == 6)
  {
  }

  const dfloat absTol = 1e-8; // values below this threshold are ignored
  const int expectedVIter = 3;

  double vxErr, prErr;
  double s1Err;
  
  int pIterErr;
  int velIterErr;

  bool passed = true;

  if(rank == 0) printf("relative error to target: ");
  for(int i = 0 ; i < Nfields; ++i){
    if(abs(err[i]) > absTol){
      const dfloat relErr = abs(err[i]-expectedErr[i])/expectedErr[i];
      if(rank == 0) printf("err[%d]=%g ", i, relErr);
      if(relErr > EPS && abs(expectedErr[i]) < abs(err[i])) passed = false;
    }
  }
  velIterErr = abs(expectedVIter - nrs->uvwSolver->Niter);
  pIterErr = abs(expectedPIter - nrs->pSolver->Niter);
  if(rank == 0) printf("velIter=%d pIter=%d\n", velIterErr, pIterErr);
  if(velIterErr >= 3) passed = false;
  if(pIterErr >= 2) passed = false;

  if(ciMode == 3)
  {
    const int expectedMeshIter =  47;
    const int meshIterErr = abs(nrs->meshSolver->Niter - expectedMeshIter);
    if(meshIterErr >= 3) passed = false;
  }

  passed ? (PASS) : (FAIL);
}
