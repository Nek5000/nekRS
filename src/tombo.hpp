#if !defined(nekrs_tombo_hpp_)
#define nekrs_tombo_hpp_

#include "nekrs.hpp"

namespace tombo {

void pressureRhs(ins_t *ins, dfloat time, occa::memory o_rhsP);
void pressureSolve(ins_t *ins, dfloat time, occa::memory o_rhsP, occa::memory o_rkP);
void velocityRhs(ins_t *ins, dfloat time, occa::memory o_rhsU);
void velocitySolve(ins_t *ins, dfloat time, occa::memory o_rhsU, occa::memory o_Uhat);
void qthermal(ins_t *ins, dfloat time, occa::memory o_qtl);

}

#endif
