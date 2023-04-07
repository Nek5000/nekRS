#if !defined(nekrs_runtime_hpp_)
#define nekrs_runtime_hpp_

#include "nrs.hpp"

namespace timeStepper {

void adjustDt(nrs_t *nrs, int tstep);

void initStep(nrs_t *nrs, dfloat time, dfloat dt, int tstep);
bool runStep(nrs_t *nrs, std::function<bool(int)> convergenceCheck, int stage);
void finishStep(nrs_t *nrs);
void advectionFlops(mesh_t *mesh, int Nfields);

void setDt(nrs_t *nrs, dfloat dt, int tstep);
void makef(nrs_t* nrs, dfloat time, int tstep, occa::memory o_FU, occa::memory o_BF);
occa::memory velocityStrongSubCycle(nrs_t* nrs, int nEXT, dfloat time, occa::memory o_U);
occa::memory velocityStrongSubCycleMovingMesh(nrs_t* nrs, int nEXT, dfloat time, occa::memory o_U);
void fluidSolve(nrs_t* nrs, dfloat time, occa::memory o_P, occa::memory o_U, int stage, int tstep);

void makeq(nrs_t *nrs, dfloat time, int tstep, occa::memory o_FS,
           occa::memory o_BF);
occa::memory scalarStrongSubCycleMovingMesh(cds_t *cds, int nEXT, dfloat time,
                                            int is, occa::memory o_U,
                                            occa::memory o_S);
occa::memory scalarStrongSubCycle(cds_t *cds, int nEXT, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S);
void scalarSolveCvode(nrs_t *nrs, dfloat tn, dfloat time, occa::memory o_S, int stage, int tstep);
void scalarSolve(nrs_t *nrs, dfloat time, occa::memory o_S, int stage);
void printInfo(nrs_t *nrs, dfloat time, int tstep, bool printStepInfo, bool printVerboseInfo);
}

#endif
