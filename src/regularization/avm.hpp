#if !defined(nekrs_avm_hpp_)
#define nekrs_avm_hpp_

#include "cds.hpp"
#include "nrs.hpp"
namespace avm{
void setup(cds_t* nrs);
void apply(cds_t* nrs, const dfloat time, const dlong scalarIndex, occa::memory o_S);
}

#endif
