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
  kernel(mesh->Nlocal, this->fieldOffset, this->o_div, o_SijOij, o_Q);
}

void nrs_t::Qcriterion(occa::memory &o_Q)
{
  Qcriterion(o_U, o_Q);
}

occa::memory nrs_t::Qcriterion(const occa::memory &o_U)
{
  auto o_Q = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
  Qcriterion(o_U, o_Q);
  return o_Q;
}

occa::memory nrs_t::Qcriterion()
{
  auto o_Q = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
  Qcriterion(this->o_U, o_Q);
  return o_Q;
}
