#if !defined(nekrs_tombo_hpp_)
#define nekrs_tombo_hpp_

#include "nekrs.hpp"
void pressureRhs(ins_t *ins, dfloat time);
void pressureSolve(ins_t *ins, dfloat time, occa::memory o_rkP);
void insGradient(ins_t *ins, dfloat time, occa::memory o_P, occa::memory o_GP);
void velocityRhs(ins_t *ins, dfloat time);
void velocitySolve(ins_t *ins, dfloat time, occa::memory o_Uhat);
void subCycle(ins_t *ins, dfloat time, int Nstages, occa::memory o_U, occa::memory o_Ud);
void advection(ins_t *ins, dfloat time);
void curlCurl(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_NC);
void insDivergence(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_DU);

#endif
