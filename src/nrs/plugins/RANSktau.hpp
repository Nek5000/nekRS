#if !defined(nekrs_RANSktau_hpp_)
#define nekrs_RANSktau_hpp_

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace RANSktau
{
void buildKernel(occa::properties kernelInfo);
void updateSourceTerms();
void setup(int startFieldIndex);
void updateProperties();
const deviceMemory<dfloat> o_mue_t();
}

#endif
