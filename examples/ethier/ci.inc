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
void ciTestLinAlg(nrs_t *nrs, const int N)
{
  linAlg_t* linAlg = nrs->linAlg;
  MPI_Comm comm = nrs->mesh->comm;
  int rank = nrs->mesh->rank;
  dfloat* x = (dfloat*) calloc(N, sizeof(dfloat));
  for(int i = 0 ; i < N; ++i)
    x[i] = drand48();
  occa::memory o_x = nrs->mesh->device.malloc(N*sizeof(dfloat), x);
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
  if (!nrs->lastStep) return;
 
  const int rank = nrs->mesh->rank;
 
  nek_ocopyFrom(time, tstep);
  nek_userchk();

  double *err = (double *) nek_scPtr(1);

  double vxErr, prErr;
  double s1Err, s2Err;
  
  int pIterErr;
  int velIterErr;

  ciTestLinAlg(nrs, BLOCKSIZE / 16);
  ciTestLinAlg(nrs, BLOCKSIZE / 8);
  ciTestLinAlg(nrs, BLOCKSIZE / 4);
  ciTestLinAlg(nrs, BLOCKSIZE / 2);
  ciTestLinAlg(nrs, BLOCKSIZE);
  ciTestLinAlg(nrs, 2 * BLOCKSIZE);
  ciTestLinAlg(nrs, 4 * BLOCKSIZE);
  ciTestLinAlg(nrs, 8 * BLOCKSIZE);
  ciTestLinAlg(nrs, 16 * BLOCKSIZE);

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
