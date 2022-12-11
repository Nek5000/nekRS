#if !defined(nekrs_RANSktau_hpp_)
#define nekrs_RANSktau_hpp_

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace RANSktau
{
void buildKernel(occa::properties kernelInfo);
void updateSourceTerms();
void setup(nrs_t* nrsIn, dfloat mue, dfloat rho, int startIndex);
void setup(nrs_t* nrsIn, dfloat mue, dfloat rho, int startIndex, const dfloat* coeffIn);
void updateProperties();
occa::memory o_mue_t();
}

#endif
