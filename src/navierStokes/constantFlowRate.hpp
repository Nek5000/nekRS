#if !defined(nekrs_constant_flow_rate_hpp_)
#define nekrs_constant_flow_rate_hpp_

#include "nrs.hpp"

namespace ConstantFlowRate{
bool apply(nrs_t *nrs, int tstep, dfloat time);
void compute(nrs_t *nrs, dfloat lengthScale, dfloat time);
bool checkIfRecompute(nrs_t* nrs, int tstep);
void printInfo(mesh_t* mesh, bool verboseInfo);
dfloat scaleFactor();
}

#endif
