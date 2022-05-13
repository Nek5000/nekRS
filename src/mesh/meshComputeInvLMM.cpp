#include "mesh.h"
#include "linAlg.hpp"
#include "platform.hpp"
void mesh_t::computeInvLMM()
{
  o_invLMM.copyFrom(o_LMM, Nelements * Np * sizeof(dfloat));
  oogs::startFinish(o_invLMM, 1, 0, ogsDfloat, ogsAdd, oogs);
  platform->linAlg->ady(Nelements * Np, 1.0, o_invLMM);
}