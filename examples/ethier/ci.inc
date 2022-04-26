#include <math.h>
#include <limits>
#include <vector>
#include <numeric>
#include <array>

static int ciMode = 0;

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("9"));
  options.setArgs("RESTART FROM FILE", std::string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", std::string("0.01"));
  options.setArgs("DENSITY", std::string("1"));
  options.setArgs("NUMBER OF SCALARS", std::string("2"));
  options.setArgs("SCALAR00 DIFFUSIVITY", std::string("0.01"));
  options.setArgs("SCALAR00 DENSITY", std::string("1"));
  options.setArgs("SCALAR01 DIFFUSIVITY", std::string("0.01"));
  options.setArgs("SCALAR01 DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("0.06"));
  options.setArgs("DT", std::string("2e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("PRESSURE INITIAL GUESS", "PREVIOUS");

  options.setArgs("VELOCITY BLOCK SOLVER", "FALSE");
  options.setArgs("SCALAR INITIAL GUESS","PREVIOUS");
  options.setArgs("VELOCITY INITIAL GUESS","PREVIOUS");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");

  if (ciMode == 2) {
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 3) {
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 4) {
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 5) {
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("0"));
  }
  if (ciMode == 6) {
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
  }
  if (ciMode == 7) {
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("PRESSURE MAXIMUM ITERATIONS", "1000");
    options.setArgs("PRESSURE PRECONDITIONER", "JACOBI");
    options.setArgs("END TIME", std::string("0.012"));
  }
  if (ciMode == 8) {
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("VARIABLE DT", "TRUE");
    options.setArgs("TARGET CFL", "0.5");
  }
  if (ciMode == 9) {
    options.setArgs("ADVECTION TYPE", "CONVECTIVE");
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 10) {
    options.setArgs("ADVECTION TYPE", "CONVECTIVE");
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 11) {
    options.setArgs("PRESSURE MULTIGRID SMOOTHER", "DAMPEDJACOBI,CHEBYSHEV");
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }

  options.setArgs("TIME INTEGRATOR", "TOMBO3");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", std::string("1e-10"));
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("VARIABLEPROPERTIES", "FALSE");
}

void surfaceArea(nrs_t* nrs, int BID, dfloat& SA)
{
  auto* mesh = nrs->meshV;

  SA = 0.0;

  for(int elem = 0; elem < mesh->Nelements; ++elem){
    for(int face = 0; face < mesh->Nfaces; ++face){
      for(int fp = 0; fp < mesh->Nfp; ++fp){
        const auto base = mesh->Nsgeo * (mesh->Nfaces * mesh->Nfp * elem + mesh->Nfp * face + fp);
        /* volume index of face node */
        const auto n = mesh->faceNodes[face * mesh->Nfp + fp];
        const auto WsJ = mesh->sgeo[base + WSJID];
        if(mesh->EToB[elem * mesh->Nfaces + face] == BID){
            SA += WsJ;
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &SA, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
}

void ciTestAvgBoundary(nrs_t* nrs)
{
  const int rank = platform->comm.mpiRank;
  auto * mesh = nrs->meshV;
  const int BID = 1;
  const auto testTol = 50. * std::numeric_limits<dfloat>::epsilon();

  constexpr int N = 5;

  dfloat SA = 0.0;

  surfaceArea(nrs, BID, SA);

  // surface area of unit sphere
  const auto SAref = 4. * M_PI;
  if(std::abs((SA - SAref)/SAref) > 0.02){
    if(rank == 0){
      std::cout << "ciTestAvgBoundary: surface area of unit sphere is " << SAref << " instead of " << SA << std::endl;
    }
    FAIL;
  }

  occa::memory o_fld = platform->device.malloc(N * nrs->fieldOffset * sizeof(dfloat));
  platform->linAlg->fill(N*nrs->fieldOffset, SA, o_fld);

  std::array<dfloat, N> results;
  mesh->avgBoundaryValue(BID, N, nrs->fieldOffset, o_fld, results.data());

  for(auto&& v : results){
    if(std::abs(v - SA) > testTol){
      if(rank == 0){
        std::cout << "ciTestAvgBoundary: FAILED" << std::endl;
      }
      FAIL;
    }
  }
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
dfloat norm2(dfloat const * const array, const int size, MPI_Comm comm)
{
  dfloat norm = 0.0;
  for(int i = 0 ; i < size; ++i)
    norm += array[i] * array[i];
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, comm);
  return sqrt(norm);
}
dfloat norm1(dfloat const * const array, const int size, MPI_Comm comm)
{
  dfloat norm = 0.0;
  for(int i = 0 ; i < size; ++i)
    norm += fabs(array[i]);
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, comm);
  return norm;
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

  const dfloat testTol = 50. * std::numeric_limits<dfloat>::epsilon();
  auto evaluate_method = [&](auto linAlgMethod, auto referenceMethod, std::string id)
  {
    const dfloat reference = referenceMethod(x,N,comm);
    const dfloat value = linAlgMethod(N, o_x, comm);
    const dfloat err = abs(
      (value - reference)/reference
    );
    if(rank == 0) printf("linAlg errs (N=%d): %s=%g\n", N, id.c_str(), err);
    if(err > testTol) FAIL;
  };

  evaluate_method(
    [&](auto N, auto o_x, auto comm) { 
      return linAlg->sum  (N, o_x, comm); 
    },
    sum, "sum");
  evaluate_method(
    [&](auto N, auto o_x, auto comm) { 
      return linAlg->max  (N, o_x, comm); 
    },
    max, "max");
  evaluate_method(
    [&](auto N, auto o_x, auto comm) { 
      return linAlg->min  (N, o_x, comm); 
    },
    min, "min");
  evaluate_method(
    [&](auto N, auto o_x, auto comm) { 
      return linAlg->norm2(N, o_x, comm); 
    },
    norm2, "norm2");
  evaluate_method(
    [&](auto N, auto o_x, auto comm) { 
      return linAlg->norm1(N, o_x, comm); 
    },
    norm1, "norm1");
  free(x);
  o_x.free();
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  const int rank = platform->comm.mpiRank;
  if(tstep == 1 && ciMode != 7){
    int NiterP = nrs->pSolver->Niter;

    // mode 11 uses Chebyshev+Jacobi smoothing
    const int expectedNiterP = (ciMode == 11) ? 13 : 6;
    const int pIterErr = abs(NiterP - expectedNiterP);
    if(pIterErr >= 2) {
      if(rank==0){
        printf("Expected %d P iterations in first step, but got %d!\n", expectedNiterP, NiterP);
      }
      FAIL;
    }
  }
  if (!nrs->lastStep) return;

  ciTestAvgBoundary(nrs);

  ciTestLinAlg(1);
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
             s1Err = abs((err[2] - 4.76E-12)/err[2]);
             s2Err = abs((err[3] - 5.54E-12)/err[3]);
             pIterErr = abs(NiterP - 4);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 7.02E-010)/err[1]);
             s01IterErr = abs(NiterS01 - 10);
             s02IterErr = abs(NiterS02 - 10);
             break;
    case 2 : 
    case 11:
             velIterErr = abs(NiterU - 10);
             s1Err = abs((err[2] - 6.11E-12)/err[2]);
             s2Err = abs((err[3] - 6.84E-12)/err[3]);
             pIterErr = abs(NiterP - 1);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 7.11E-10)/err[1]);
             s01IterErr = abs(NiterS01 - 10);
             s02IterErr = abs(NiterS02 - 10);
             break;
    case 3 : velIterErr = abs(NiterU - 4);
             s1Err = abs((err[2] - 5.25E-12)/err[2]);
             s2Err = abs((err[3] - 6.09E-12)/err[3]);
             pIterErr = abs(NiterP - 2);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 7.11E-10)/err[1]);
             s01IterErr = abs(NiterS01 - 10);
             s02IterErr = abs(NiterS02 - 10);
             break;
    case 4 : velIterErr = abs(NiterU - 1);
             s1Err = abs((err[2] - 1.71E-11)/err[2]);
             s2Err = abs((err[3] - 2.00E-11)/err[3]);
             pIterErr = abs(NiterP - 1);
             vxErr = abs((err[0] - 3.65E-10)/err[0]);
             prErr = abs((err[1] - 6.71E-10)/err[1]);
             s01IterErr = abs(NiterS01 - 10);
             s02IterErr = abs(NiterS02 - 10);
             break;
    case 5 : velIterErr = abs(NiterU - 12);
             s1Err = abs((err[2] - 1.53E-06)/err[2]);
             s2Err = abs((err[3] - 1.73E-06)/err[3]);
             pIterErr = abs(NiterP - 12);
             vxErr = abs((err[0] - 2.57E-07)/err[0]);
             prErr = abs((err[1] - 2.16E-06)/err[1]);
             s01IterErr = abs(NiterS01 - 10);
             s02IterErr = abs(NiterS02 - 10);
             break;
    case 6 : velIterErr = abs(NiterU - 12);
             s1Err = abs((err[2] - 1.27E-06)/err[2]);
             s2Err = abs((err[3] - 1.41E-06)/err[3]);
             pIterErr = abs(NiterP - 12);
             vxErr = abs((err[0] - 2.43E-07)/err[0]);
             prErr = abs((err[1] - 2.89E-06)/err[1]);
             s01IterErr = abs(NiterS01 - 10);
             s02IterErr = abs(NiterS02 - 10);
             break;
    case 7 : velIterErr = abs(NiterU - 5);
             s1Err = abs((err[2] - 2E-13)/err[2]);
             s2Err = abs((err[3] - 2E-13)/err[3]);
             
             // just check for convergence
             pIterErr = (NiterP < 999) ? 0 : 1000;

             vxErr = abs((err[0] - 1.4E-10)/err[0]);
             prErr = abs((err[1] - 8.7E-9)/err[1]);
             s01IterErr = abs(NiterS01 - 10);
             s02IterErr = abs(NiterS02 - 10);
             break;
    case 8 : velIterErr = abs(NiterU - 10);
             s1Err = abs((err[2] - 1.03E-11)/err[2]);
             s2Err = abs((err[3] - 1.34E-11)/err[3]);
             pIterErr = abs(NiterP - 2);
             vxErr = abs((err[0] - 3.38E-10)/err[0]);
             prErr = abs((err[1] - 7.8E-10)/err[1]);
             s01IterErr = abs(NiterS01 - 10);
             s02IterErr = abs(NiterS02 - 10);
             break;
    case 9 : velIterErr = abs(NiterU - 10);
             s1Err = abs((err[2] - 1.03E-11)/err[2]);
             s2Err = abs((err[3] - 1.34E-11)/err[3]);
             pIterErr = abs(NiterP - 1);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 6.98E-10)/err[1]);
             s01IterErr = abs(NiterS01 - 10);
             s02IterErr = abs(NiterS02 - 10);
             break;
    case 10 : velIterErr = abs(NiterU - 10);
             s1Err = abs((err[2] - 1.03E-11)/err[2]);
             s2Err = abs((err[3] - 1.34E-11)/err[3]);
             pIterErr = abs(NiterP - 1);
             vxErr = abs((err[0] - 2.78E-10)/err[0]);
             prErr = abs((err[1] - 6.98E-10)/err[1]);
             s01IterErr = abs(NiterS01 - 10);
             s02IterErr = abs(NiterS02 - 10);
             break;
     }

  if (ciMode == 4){
    dfloat memoryUsage = platform->device.occaDevice().memoryAllocated()/1e9;
    const dfloat expectedMemoryUsage = 0.0555274;
    if(memoryUsage > 1.01 * expectedMemoryUsage && platform->comm.mpiCommSize == 2)
    {
      if(platform->comm.mpiRank == 0){
        printf("Expected memory usage exceeded!\n"
          "Expected : %g GB\n"
          "Actual   : %g GB\n",
          expectedMemoryUsage, memoryUsage);
        fflush(stdout);
      }
      FAIL;
    }
  }
  if (ciMode == 8){
    const int expectedSteps = 68;
    const int errSteps = abs(expectedSteps - tstep);
    if (rank == 0)
      printf("nSteps err variable dt:  %d\n", errSteps);
    (errSteps <= 2) ? ({/** no-op **/}) : (FAIL);
  }

  if(ciMode == 5 || ciMode == 6)
  {
    if(platform->options.compareArgs("MOVING MESH", "FALSE"))
    {
      if(rank == 0)
        printf("Unable to test moving mesh in ethier. Please re-run the case with mv_ethier!\n");
      (FAIL);
    }
  }

  if (rank == 0)
    printf("relative error to target: vx=%g pr=%g s1=%g s2=%g velIter=%d pIter=%d s01Iter=%d s02Iter=%d\n", 
           vxErr, prErr, s1Err, s2Err, velIterErr, pIterErr, s01IterErr, s02IterErr);
  
  std::vector<double> relErrors = {vxErr, prErr, s1Err, s2Err};
  std::vector<double> absErrors = {abs(err[0]), abs(err[1]), abs(err[2]), abs(err[3])};
  bool passTest = true;
  const double absTol = 5e-11; // values below this threshold are ignored
  for(size_t i = 0; i < relErrors.size(); ++i){

    const double absError = absErrors[i];
    const double relError = relErrors[i];
    if(absError > absTol){
      passTest &= relError < EPS;
    }
  }

  passTest &= velIterErr <= 1;
  passTest &= pIterErr <= 2;
  passTest &= s01IterErr <= 1;
  passTest &= s02IterErr <= 1;

  passTest ? (PASS) : (FAIL);

}
