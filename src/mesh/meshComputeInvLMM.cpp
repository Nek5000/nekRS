#include "mesh.h"
#include "linAlg.hpp"
#include "platform.hpp"
void mesh_t::computeInvLMM()
{
  o_invLMM.copyFrom(o_LMM, Nlocal);
  oogs::startFinish(o_invLMM, 1, 0, ogsDfloat, ogsAdd, oogs);
  platform->linAlg->ady(Nelements * Np, 1.0, o_invLMM);

  if (!o_invAJwTimesInvDegree.isInitialized()) {
    o_invAJwTimesInvDegree = platform->device.malloc<dfloat>(Nlocal);
  }
  platform->linAlg->axmyz(Nlocal, 1 / volume, o_invLMM, ogs->o_invDegree, o_invAJwTimesInvDegree);
}
