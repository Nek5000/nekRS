#if !defined(nekrs_tavg_hpp_)
#define nekrs_tavg_hpp_

/*
     Statistics can be obtained from runtime averages:

     <X>    := AVG(X)
     <X'Y') := AVG(X*Y) - AVG(X)*AVG(Y)
*/

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

#include <vector>

namespace tavg
{
typedef std::vector< std::vector<occa::memory> > fields;

void buildKernel(occa::properties kernelInfo);
void run(dfloat time);
void setup(nrs_t *nrs_, const fields& fields);
void setup(nrs_t* nrs_);
void outfld();
void outfld(int outXYZ, int FP64);
void reset();
occa::memory o_avg();
}

#endif
