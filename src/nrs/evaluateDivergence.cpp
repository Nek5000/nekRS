#include "nrs.hpp"

void nrs_t::evaluateDivergence(const double time) 
{
  if (this->userDivergence) {
    platform->timer.tic("udfDiv", 1);
    platform->linAlg->fill(this->mesh->Nlocal, 0.0, this->o_div);
    this->userDivergence(time);
    platform->timer.toc("udfDiv");
  }
}
