#if !defined(nekrs_lowMach_hpp_)
#define nekrs_lowMach_hpp_ 

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace lowMach
{
// alphaRef := p0thRef/(rhoRef * cpRef * TRef)
// use alphaRef = 1 when solving for a dimensional formulation 
void setup(nrs_t* nrs, dfloat alphaRef, occa::memory& o_beta, occa::memory& o_kappa);
void buildKernel(occa::properties kernelInfo);
void qThermalSingleComponent(dfloat time, occa::memory& o_div);
void dpdt(occa::memory& o_FU);
}

#endif
