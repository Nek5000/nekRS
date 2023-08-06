#include "nrs.hpp"
#include "platform.hpp"
#include "linAlg.hpp"
#include "postProcessing.hpp"

void postProcessing::Qcriterion(nrs_t *nrs, occa::memory& o_Q)
{
  occa::memory o_SijOij = platform->o_mempool.slice2;
  strainRotationRate(nrs, true, true, o_SijOij); 

  auto kernel = platform->kernels.get("Qcriterion");
  kernel(nrs->meshV->Nlocal, nrs->fieldOffset, nrs->o_div, o_SijOij, o_Q);
}
