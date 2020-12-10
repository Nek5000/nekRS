#if !defined(nekrs_tombo_hpp_)
#define nekrs_tombo_hpp_

#include "nrs.hpp"

namespace tombo
{
occa::memory pressureSolve(nrs_t* nrs, dfloat time);
occa::memory velocitySolve(nrs_t* nrs, dfloat time);
}

#endif
