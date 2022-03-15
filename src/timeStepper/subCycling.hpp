#if !defined(nekrs_subcycle_hpp_)
#define nekrs_subcycle_hpp_

#include "nrs.hpp"

occa::memory velocitySubCycle(nrs_t* nrs, int nEXT, dfloat time, occa::memory o_U);
occa::memory velocitySubCycleMovingMesh(nrs_t* nrs, int nEXT, dfloat time, occa::memory o_U);
occa::memory scalarSubCycleMovingMesh(cds_t *cds, int nEXT, dfloat time,
                                      int is, occa::memory o_U,
                                      occa::memory o_S);
occa::memory scalarSubCycle(cds_t *cds, int nEXT, dfloat time, int is,
                            occa::memory o_U, occa::memory o_S);

#endif
