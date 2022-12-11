/*
     Compute expected value E(X) 
     Statistics can be obtained by:

     avg(X)   := E(X)
     var(X)   := E(X*X) - E(X)*E(X)
     cov(X,Y) := E(X*Y) - E(X)*E(Y)

     Note: E-operator is linear, in the sense that the expected
           value is given by E(X) = 1/N * sum[ E(X)_i ], where E(X)_i
           is the expected value of the sub-ensemble i (i=1...M).
 */

#if !defined(nekrs_tavg_hpp_)
#define nekrs_tavg_hpp_

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
occa::memory userFieldAvg();
}

#endif
