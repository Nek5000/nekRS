#include "nrs.hpp"
#include "platform.hpp"
#include "linAlg.hpp"
#include "postProcessing.hpp"

void postProcessing::Qcriterion(nrs_t *nrs, const occa::memory& o_U, occa::memory& o_Q)
{
  occa::memory o_SijOij = 
      platform->o_memPool.reserve<dfloat>(3 * nrs->NVfields * static_cast<size_t>(nrs->fieldOffset));
  strainRotationRate(nrs, true, true, o_U, o_SijOij); 

  auto kernel = platform->kernels.get("Qcriterion");
  kernel(nrs->meshV->Nlocal, nrs->fieldOffset, nrs->o_div, o_SijOij, o_Q);
}
