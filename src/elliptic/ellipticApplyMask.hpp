#ifndef ellipticApplyMask_hpp
#define ellipticApplyMask_hpp

#include <string>
#include "occa.hpp"

class elliptic_t;

void ellipticApplyMask(elliptic_t *solver, occa::memory &o_x, std::string precision);
void ellipticApplyMask(elliptic_t *solver,
                       dlong Nelements,
                       dlong Nmasked,
                       const occa::memory &o_elementList,
                       const occa::memory &o_maskIds,
                       occa::memory &o_x,
                       std::string precision);
#endif
