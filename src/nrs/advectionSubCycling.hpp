#if !defined(nekrs_subcycle_hpp_)
#define nekrs_subcycle_hpp_

#include "nekrsSys.hpp"
#include "mesh.h"

occa::memory advectionSubcyclingRK(mesh_t *_mesh, mesh_t *meshV,
                                   double time, dfloat *dt, int Nsubsteps, dfloat *coeffBDF, int nEXT,
                                   int nFields, occa::kernel _opKernel, oogs_t *_gsh,
                                   dlong _meshOffset, dlong _fieldOffset, dlong cubatureOffset, dlong fieldOffsetSum,
                                   occa::memory o_Urst, occa::memory o_U);

#endif
