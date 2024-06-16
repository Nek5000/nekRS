#if !defined(nekrs_avm_hpp_)
#define nekrs_avm_hpp_

#include "cds.hpp"
namespace avm{
void setup(cds_t* cds);
void apply(cds_t* cds, const double time, const dlong scalarIndex, occa::memory o_S);
}

#endif
