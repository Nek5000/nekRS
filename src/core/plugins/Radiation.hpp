#if !defined(nekrs_Radiation_hpp_)
#define nekrs_Radiation_hpp_

#include "nrs.hpp"

// Monte Carlo radiative view-factor computation between specified boundary-ID
// groups, enabled via the [RADIATION] .par section. Runs once at setup (pure
// geometry, independent of the flow field); does not participate in the time
// loop. See src/core/plugins/Radiation.cpp for the implementation.
namespace Radiation
{
void buildKernel(occa::properties kernelInfo);
void setup();
}

#endif
