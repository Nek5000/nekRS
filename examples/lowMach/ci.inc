#include <math.h>

static int ciMode = 0;

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-2

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("7"));
  options.setArgs("RESTART FROM FILE", string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("END TIME", string("0.3"));
  options.setArgs("DT", string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", string("0"));
  if (ciMode == 2) options.setArgs("SUBCYCLING STEPS", string("1"));
  options.setArgs("TIME INTEGRATOR", "TOMBO2");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-12"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-08"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", string("1e-12"));
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  if (!nrs->lastStep) return;
 
  const int rank = platform->comm.mpiRank;

  nrs->o_div.copyTo(nrs->div);
  dlong Nlocal = nrs->meshV->Nelements * nrs->meshV->Np;
  memcpy(nekData.qtl, nrs->div, sizeof(dfloat)*Nlocal);
 
  nek::ocopyToNek(time, tstep);
  nek::userchk();

  double *err = (double *) nek::scPtr(1);
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
