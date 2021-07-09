#include <math.h>
#include <limits>

static int ciMode = 0;

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("9"));
  options.setArgs("RESTART FROM FILE", string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
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
  if (ciMode == 5) {
    options.setArgs("END TIME", string("0.2"));
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", string("0"));
    options.setArgs("MOVING MESH", string("TRUE"));
  }
  if (ciMode == 6) {
    options.setArgs("END TIME", string("0.2"));
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", string("1"));
    options.setArgs("MOVING MESH", string("TRUE"));
  }
  if (ciMode == 7) {
    options.setArgs("PRESSURE RESIDUAL PROJECTION", "TRUE");
    options.setArgs("PRESSURE MAXIMUM ITERATIONS", "1000");
    options.setArgs("PRESSURE PRECONDITIONER", "JACOBI");
    options.setArgs("VELOCITY RESIDUAL PROJECTION", "TRUE");
    options.setArgs("END TIME", string("0.012"));
  }

  options.setArgs("TIME INTEGRATOR", "TOMBO3");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-10"));
  options.setArgs("SCALAR00 SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("VARIABLEPROPERTIES", "FALSE");
}

dfloat sum(dfloat const * const array, const int size, MPI_Comm comm)
{
  dfloat sumr = 0.0;
  for(int i = 0 ; i < size; ++i)
    sumr += array[i];
  MPI_Allreduce(MPI_IN_PLACE, &sumr, 1, MPI_DFLOAT, MPI_SUM, comm);
  return sumr;
}
dfloat max(dfloat const * const array, const int size, MPI_Comm comm)
{
  dfloat maxr = -9e30;
  for(int i = 0 ; i < size; ++i)
    maxr = (array[i] > maxr) ? array[i] : maxr;
  MPI_Allreduce(MPI_IN_PLACE, &maxr, 1, MPI_DFLOAT, MPI_MAX, comm);
  return maxr;
}
dfloat min(dfloat const * const array, const int size, MPI_Comm comm)
{
  dfloat minr = 9e30;
  for(int i = 0 ; i < size; ++i)
    minr = (array[i] < minr) ? array[i] : minr;
  MPI_Allreduce(MPI_IN_PLACE, &minr, 1, MPI_DFLOAT, MPI_MIN, comm);
  return minr;
}
void ciTestLinAlg(const int N)
{
  linAlg_t* linAlg = platform->linAlg;
  MPI_Comm comm = platform->comm.mpiComm;
  int rank = platform->comm.mpiRank;
  dfloat* x = (dfloat*) calloc(N, sizeof(dfloat));
  for(int i = 0 ; i < N; ++i)
    x[i] = drand48();
  occa::memory o_x = platform->device.malloc(N*sizeof(dfloat), x);
  const dfloat referenceSum = sum(x,N,comm);
  const dfloat referenceMax = max(x,N,comm);
  const dfloat referenceMin = min(x,N,comm);
  const dfloat err_sum = abs(
    (referenceSum - linAlg->sum(N, o_x, comm))/referenceSum
  );
  const dfloat err_max = abs(
    (referenceMax - linAlg->max(N, o_x, comm))/referenceMax
  );
  const dfloat err_min = abs(
    (referenceMin - linAlg->min(N, o_x, comm))/referenceMin
  );
  if(rank == 0) printf("linAlg errs (N=%d): max=%g, min=%g, sum=%g\n", N, err_max, err_min, err_sum);
  free(x);
  o_x.free();
  const dfloat testTol = 50. * std::numeric_limits<dfloat>::epsilon();
  if(err_max > testTol || err_min > testTol || err_sum > testTol)
    FAIL;
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  const int rank = platform->comm.mpiRank;
  if(tstep == 1 && ciMode != 7){
    int NiterP = nrs->pSolver->Niter;
    const int expectedNiterP = 6;
    const int pIterErr = abs(NiterP - expectedNiterP);
    if(pIterErr >= 2) {
      if(rank==0){
        printf("Expected %d P iterations in first step, but got %d!\n", expectedNiterP, NiterP);
      }
      FAIL;
    }
  }
  if (!nrs->lastStep) return;
 

  ciTestLinAlg(BLOCKSIZE / 16);
  ciTestLinAlg(BLOCKSIZE / 8);
  ciTestLinAlg(BLOCKSIZE / 4);
  ciTestLinAlg(BLOCKSIZE / 2);
  ciTestLinAlg(BLOCKSIZE);
  ciTestLinAlg(2 * BLOCKSIZE);
  ciTestLinAlg(4 * BLOCKSIZE);
  ciTestLinAlg(8 * BLOCKSIZE);
  ciTestLinAlg(16 * BLOCKSIZE);

  nek::ocopyToNek(time, tstep);
  nek::userchk();

  double *err = (double *) nek::scPtr(1);

  double vxErr, prErr;
  double s1Err, s2Err;
  
  int pIterErr;
  int velIterErr;

  int NiterU;

  if(platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE"))
    NiterU = nrs->uvwSolver->Niter;
  else
    NiterU = nrs->uSolver->Niter;

  int NiterP = nrs->pSolver->Niter;

  int NiterS01 = nrs->cds->solver[0]->Niter;
  int NiterS02 = nrs->cds->solver[1]->Niter;
  int s01IterErr, s02IterErr;

  switch (ciMode) {
    case 1 : velIterErr = abs(NiterU - 10);
             s1Err = abs((err[2] - 5.25E-12)/err[2]);
             s2Err = abs((err[3] - 6.09E-12)/err[3]);
             pIterErr = abs(NiterP - 2);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 7.02E-010)/err[1]);
             s01IterErr = abs(NiterS01 - 2);
             s02IterErr = abs(NiterS02 - 2);
             break;
    case 2 : velIterErr = abs(NiterU - 10);
             s1Err = abs((err[2] - 6.49E-12)/err[2]);
             s2Err = abs((err[3] - 7.28E-12)/err[3]);
             pIterErr = abs(NiterP - 2);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 7.11E-10)/err[1]);
             s01IterErr = abs(NiterS01 - 2);
             s02IterErr = abs(NiterS02 - 2);
             break;
    case 3 : velIterErr = abs(NiterU - 4);
             s1Err = abs((err[2] - 5.25E-12)/err[2]);
             s2Err = abs((err[3] - 6.09E-12)/err[3]);
             pIterErr = abs(NiterP - 2);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 7.11E-10)/err[1]);
             s01IterErr = abs(NiterS01 - 2);
             s02IterErr = abs(NiterS02 - 2);
             break;
    case 4 : velIterErr = abs(NiterU - 1);
             s1Err = abs((err[2] - 1.71E-11)/err[2]);
             s2Err = abs((err[3] - 2.00E-11)/err[3]);
             pIterErr = abs(NiterP - 1);
             vxErr = abs((err[0] - 3.65E-10)/err[0]);
             prErr = abs((err[1] - 6.71E-10)/err[1]);
             s01IterErr = abs(NiterS01 - 2);
             s02IterErr = abs(NiterS02 - 2);
             break;
    case 5 : velIterErr = abs(NiterU - 12);
             s1Err = abs((err[2] - 1.53E-06)/err[2]);
             s2Err = abs((err[3] - 1.73E-06)/err[3]);
             pIterErr = abs(NiterP - 12);
             vxErr = abs((err[0] - 2.57E-07)/err[0]);
             prErr = abs((err[1] - 2.16E-06)/err[1]);
             s01IterErr = abs(NiterS01 - 5);
             s02IterErr = abs(NiterS02 - 5);
             break;
    case 6 : velIterErr = abs(NiterU - 12);
             s1Err = abs((err[2] - 1.27E-06)/err[2]);
             s2Err = abs((err[3] - 1.41E-06)/err[3]);
             pIterErr = abs(NiterP - 12);
             vxErr = abs((err[0] - 2.43E-07)/err[0]);
             prErr = abs((err[1] - 2.89E-06)/err[1]);
             s01IterErr = abs(NiterS01 - 5);
             s02IterErr = abs(NiterS02 - 5);
             break;
    case 7 : velIterErr = abs(NiterU - 5);
             s1Err = abs((err[2] - 2E-13)/err[2]);
             s2Err = abs((err[3] - 2E-13)/err[3]);
             pIterErr = abs(NiterP - 409);
             vxErr = abs((err[0] - 1.4E-10)/err[0]);
             prErr = abs((err[1] - 8.7E-9)/err[1]);
             s01IterErr = abs(NiterS01 - 2);
             s02IterErr = abs(NiterS02 - 2);
             break;

     }

  if (rank == 0)
    printf("relative error to target: vx=%g pr=%g s1=%g s2=%g velIter=%d pIter=%d s01Iter=%d s02Iter=%d\n", 
           vxErr, prErr, s1Err, s2Err, velIterErr, pIterErr, s01IterErr, s02IterErr);

  (vxErr < EPS && prErr < EPS && s1Err < EPS && s2Err < EPS && 
   velIterErr <= 1 && pIterErr <= 2 && s01IterErr <= 1 && s02IterErr <= 1) ? (PASS) : (FAIL); 
}
