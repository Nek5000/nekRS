#if !defined(nekrs_tombo_hpp_)
#define nekrs_tombo_hpp_

#include "nekrs.hpp"

namespace tombo {

void pressureSolve(ins_t *ins, dfloat time, occa::memory o_wrk6, occa::memory o_Pnew);
void velocitySolve(ins_t *ins, dfloat time, occa::memory o_wrk6, occa::memory o_Unew);

}

#endif
