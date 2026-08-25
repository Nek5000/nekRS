#if !defined(nekrs_Radiation_hpp_)
#define nekrs_Radiation_hpp_

#include "nrs.hpp"

// Monte Carlo radiative view-factor computation between specified boundary-ID
// groups, enabled via the [RADIATION] .par section. setup() runs once (pure
// geometry, independent of the flow field). step() couples the resulting
// view-factor matrix to the temperature field via a gray-diffuse radiosity
// solve, recomputing the radiative flux BC every RADIATION updateFrequency
// steps; call it from UDF_ExecuteStep. See src/core/plugins/Radiation.cpp.
namespace Radiation
{
void buildKernel(occa::properties kernelInfo);
void setup();
void step(double time, int tstep);
}

#endif
