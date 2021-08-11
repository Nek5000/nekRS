#include <math.h>
#include <limits>

static int ciMode = 0;

#define PASS { if (rank == 0) printf("TESTS passed \n"); MPI_Finalize(); exit(0); }
#define FAIL { if (rank == 0) printf("TESTS failed!\n"); MPI_Finalize(); exit(2); }

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", string("3"));
  options.setArgs("RESTART FROM FILE", string("0"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", std::to_string(1.0 / 19000.0));
  options.setArgs("DENSITY", string("1"));
  options.setArgs("END TIME", string("6e-1"));
  options.setArgs("DT", string("6e-3"));
  options.setArgs("SUBCYCLING STEPS", string("1"));
  options.setArgs("PRESSURE INITIAL GUESS", "PREVIOUS");

  options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
  options.setArgs("VELOCITY INITIAL GUESS","PREVIOUS");

  options.setArgs("TIME INTEGRATOR", "TOMBO2");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("VELOCITY SOLVER TOLERANCE", string("1e-6"));
  options.setArgs("PRESSURE SOLVER TOLERANCE", string("1e-4"));
  options.setArgs("VARIABLEPROPERTIES", "FALSE");
}

void ciTestErrors(nrs_t *nrs, dfloat time, int tstep)
{
  const int rank = platform->comm.mpiRank;
  mesh_t* mesh = nrs->meshV;
  occa::memory o_Uz = nrs->o_U + 2*nrs->fieldOffset*sizeof(dfloat);
  const dfloat ubar = platform->linAlg->innerProd(mesh->Nlocal, o_Uz, mesh->o_LMM, platform->comm.mpiComm)/mesh->volume;

  const dfloat expected_ubar = 0.9;
  const dfloat TOL = 5e-5;

  const dfloat err = abs(ubar - expected_ubar);

  if(rank == 0)
    printf("abs err in bulk velocity: %g\n", err);

  if(abs(ubar-expected_ubar) > TOL)
    (FAIL);
  else
    (PASS);
}
