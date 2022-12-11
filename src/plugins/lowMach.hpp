#if !defined(nekrs_lowMach_hpp_)
#define nekrs_lowMach_hpp_

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace lowMach
{
void setup(nrs_t* nrs, dfloat gamma0);
void buildKernel(occa::properties kernelInfo);
void qThermalIdealGasSingleComponent(dfloat time, occa::memory o_div);
void dpdt(occa::memory o_FU);
}

#endif
