#if !defined(nekrs_runtime_hpp_)
#define nekrs_runtime_hpp_

#include "nrs.hpp"

namespace timeStepper {

void step(nrs_t* nrs, dfloat time, dfloat dt, int tstep);
void coeffs(nrs_t *nrs, dfloat dt, int tstep);
void makef(nrs_t* nrs, dfloat time, int tstep, occa::memory o_FU, occa::memory o_BF);
occa::memory velocityStrongSubCycle(nrs_t* nrs, int nEXT, dfloat time, occa::memory o_U);
occa::memory velocityStrongSubCycleMovingMesh(nrs_t* nrs, int nEXT, dfloat time, occa::memory o_U);
void fluidSolve(nrs_t* nrs, dfloat time, occa::memory o_U, int stage);
void meshSolve(nrs_t* nrs, dfloat time, occa::memory o_U, int stage);

void makeq(nrs_t* nrs, dfloat time, int tstep, occa::memory o_FS, occa::memory o_BF);
occa::memory scalarStrongSubCycleMovingMesh(cds_t* cds, int nEXT, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S);
occa::memory scalarStrongSubCycle(cds_t* cds, int nEXT, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S);
void scalarSolve(nrs_t* nrs, dfloat time, occa::memory o_S, int stage);
void printInfo(nrs_t* nrs, dfloat time, int tstep, double tElapsedStep, double tElapsed);
}

#endif
