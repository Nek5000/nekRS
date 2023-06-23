#if !defined(nekrs_subcycle_hpp_)
#define nekrs_subcycle_hpp_

#include "nrs.hpp"

occa::memory velocitySubCycle(nrs_t* nrs, int nEXT, double time);
occa::memory scalarSubCycle(nrs_t* nrs, int nEXT, double time, int scalarIdx);

#endif
