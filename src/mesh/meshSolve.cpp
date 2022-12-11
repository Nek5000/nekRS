#include "nrs.hpp"
#include "mesh.h"
void meshSolve(nrs_t* nrs, dfloat time, occa::memory o_U, int stage)
{
  mesh_t* mesh = nrs->meshV;
  linAlg_t* linAlg = platform->linAlg;

  platform->timer.tic("meshSolve", 1);
  nrs->setEllipticCoeffKernel(
    mesh->Nlocal,
    1.0,
    0 * nrs->fieldOffset,
    nrs->fieldOffset,
    nrs->o_meshMue,
    nrs->o_meshRho,
    nrs->o_ellipticCoeff);

  occa::memory o_Unew = [&](nrs_t* nrs, dfloat time, int stage){
    mesh_t* mesh = nrs->meshV;
    oogs_t* gsh = nrs->gsh;

    platform->linAlg->fill(nrs->NVfields*nrs->fieldOffset, 0, platform->o_mempool.slice3);

    const occa::memory& o_U0 = nrs->meshSolver->options.compareArgs("INITIAL GUESS", "EXTRAPOLATION") && stage == 1 ?   
                               mesh->o_Ue : mesh->o_U;
    platform->o_mempool.slice0.copyFrom(o_U0, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
    ellipticSolve(nrs->meshSolver, platform->o_mempool.slice3, platform->o_mempool.slice0);

    // enforce C0
    oogs::startFinish(platform->o_mempool.slice0, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);
    platform->linAlg->axmyMany(
      mesh->Nlocal,
      nrs->NVfields,
      nrs->fieldOffset,
      0,
      1.0,
      nrs->meshSolver->o_invDegree,
      platform->o_mempool.slice0
    );

    return platform->o_mempool.slice0;
  }(nrs, time, stage);
  o_U.copyFrom(o_Unew, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
  platform->timer.toc("meshSolve");
}
