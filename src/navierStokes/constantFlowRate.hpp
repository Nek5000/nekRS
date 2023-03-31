#if !defined(nekrs_constant_flow_rate_hpp_)
#define nekrs_constant_flow_rate_hpp_

#include "nrs.hpp"

namespace ConstantFlowRate{
bool adjust(nrs_t *nrs, int tstep, dfloat time);
void printInfo(mesh_t* mesh, bool verboseInfo);
dfloat scaleFactor();
}

#endif
