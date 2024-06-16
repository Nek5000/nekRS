#include "nrs.hpp"
#include "platform.hpp"
#include "linAlg.hpp"

void nrs_t::Qcriterion(const occa::memory &o_U, occa::memory &o_Q)
{
  auto o_SijOij = this->strainRotationRate(o_U);

  static occa::kernel kernel;
  if (!kernel.isInitialized()) {
    kernel = platform->kernelRequests.load("nrs-Qcriterion");
  }
  kernel(this->meshV->Nlocal, this->fieldOffset, this->o_div, o_SijOij, o_Q);
}

occa::memory nrs_t::Qcriterion(const occa::memory &o_U)
{
  auto o_Q = platform->o_memPool.reserve<dfloat>(this->meshV->Nlocal);
  Qcriterion(o_U, o_Q);
  return o_Q;
}

occa::memory nrs_t::Qcriterion()
{
  auto o_Q = platform->o_memPool.reserve<dfloat>(this->meshV->Nlocal);
  Qcriterion(this->o_U, o_Q);
  return o_Q;
}
