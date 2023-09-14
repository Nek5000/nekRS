#if !defined(nekrs_applyDirichlet_hpp_)
#define nekrs_applyDirichlet_hpp_

#include "nrs.hpp"
void applyDirichlet(nrs_t *nrs, double time);
void applyDirichletVelocity(nrs_t *nrs, double time, occa::memory& o_U,occa::memory& o_Ue,occa::memory& o_P);
void applyDirichletMesh(nrs_t *nrs, double time, occa::memory& o_UM, occa::memory& o_UMe, occa::memory& o_U);
void applyDirichletScalars(nrs_t *nrs, double time, occa::memory& o_S, occa::memory& o_Se);

#endif
