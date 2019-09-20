#if !defined(nekrs_velrecycling_hpp_)
#define nekrs_velrecycling_hpp_

#include "nekrs.hpp"
void buildKernel(ins_t *ins);
void copy();
void setup(ins_t *ins_, occa::memory o_wrk_, const hlong eOffset, const int bID_,
           const dfloat wbar_);

#endif
