#include "nrs.hpp"
#include "mesh.h"

occa::memory meshSolve(nrs_t* nrs, double time, int stage)
{
  mesh_t *mesh = nrs->_mesh;
  linAlg_t* linAlg = platform->linAlg;

  auto o_rhs = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
  platform->linAlg->fill(nrs->NVfields * nrs->fieldOffset, 0, o_rhs);

  platform->timer.tic("meshSolve", 1);

  nrs->setEllipticCoeffKernel(
    mesh->Nlocal,
    1.0,
    0 * nrs->fieldOffset,
    nrs->fieldOffset,
    0,
    nrs->o_meshMue,
    nrs->o_meshRho,
    o_NULL,
    nrs->o_ellipticCoeff);

  auto o_U = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
  if (platform->options.compareArgs("MESH INITIAL GUESS", "EXTRAPOLATION") && stage == 1)
    o_U.copyFrom(mesh->o_Ue);
  else
    o_U.copyFrom(mesh->o_U);

  ellipticSolve(nrs->meshSolver, o_rhs, o_U);

  platform->timer.toc("meshSolve");

  return o_U;
}
