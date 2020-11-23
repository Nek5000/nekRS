#if !defined(nekrs_tombo_hpp_)
#define nekrs_tombo_hpp_

#include "nrs.hpp"

namespace tombo
{
occa::memory pressureSolve(ins_t* ins, dfloat time);
occa::memory velocitySolve(ins_t* ins, dfloat time);
}

#endif
