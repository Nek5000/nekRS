#ifndef ellipticApplyMask_hpp
#define ellipticApplyMask_hpp

#include <string>
#include "occa.hpp"

class elliptic_t;

void applyMask(elliptic_t *solver, occa::memory &o_x, std::string precision);
void applyMaskInterior(elliptic_t *solver, occa::memory &o_x, std::string precision);
void applyMaskExterior(elliptic_t *solver, occa::memory &o_x, std::string precision);
#endif