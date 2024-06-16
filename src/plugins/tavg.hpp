#if !defined(nekrs_tavg_hpp_)
#define nekrs_tavg_hpp_

/*
     Statistics can be obtained from runtime averages:

     <X>    := AVG(X)
     <X'Y'> := AVG(X*Y) - AVG(X)*AVG(Y)
*/

#include "nekInterfaceAdapter.hpp"

#include <vector>

namespace tavg
{
void buildKernel(occa::properties kernelInfo);
void run(double time);
void setup(dlong fieldOffset, const std::vector< std::vector<deviceMemory<dfloat>> >& fields);
void outfld();
void outfld(int outXYZ, int FP64);
void reset();
deviceMemory<dfloat> o_avg();
}

#endif
