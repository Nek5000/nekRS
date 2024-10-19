#include "nrs.hpp"

static occa::memory _strainRotationRate(nrs_t *nrs, bool rotationRate, const occa::memory &o_U, bool smooth)
{
  auto mesh = nrs->mesh;

  const int nFields = (rotationRate) ? 2 * nrs->NVfields + nrs->NVfields : 2 * nrs->NVfields;

  auto o_SO = platform->deviceMemoryPool.reserve<dfloat>(nFields * nrs->fieldOffset);

  nrs->SijOijKernel(mesh->Nelements,
                    nrs->fieldOffset,
                    static_cast<int>(rotationRate),
                    static_cast<int>(smooth),
                    mesh->o_vgeo,
                    mesh->o_D,
                    o_U,
                    o_SO);

  if (smooth) {
    oogs::startFinish(o_SO, nFields, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);

    platform->linAlg->axmyMany(mesh->Nlocal, nFields, nrs->fieldOffset, 0, 1.0, mesh->o_invLMM, o_SO);
  }

  return o_SO;
}

occa::memory nrs_t::strainRotationRate(bool smooth)
{
  return _strainRotationRate(this, true, this->o_U, smooth);
}

occa::memory nrs_t::strainRotationRate(const occa::memory &o_U, bool smooth)
{
  return _strainRotationRate(this, true, o_U, smooth);
}

occa::memory nrs_t::strainRate(bool smooth)
{
  return _strainRotationRate(this, false, this->o_U, smooth);
}

occa::memory nrs_t::strainRate(const occa::memory &o_U, bool smooth)
{
  return _strainRotationRate(this, false, o_U, smooth);
}
