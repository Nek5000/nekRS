/*

   The MIT License (MIT)

   Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

 */

#include "elliptic.h"
#include "ellipticPrecon.h"
#include "ellipticMultiGrid.h"
#include "linAlg.hpp"

void ellipticUpdateJacobi(elliptic_t *elliptic, occa::memory &o_invDiagA)
{
  auto mesh = elliptic->mesh;
  dfloat flopCount = 0.0;
  if (elliptic->mgLevel) {
    elliptic->ellipticBlockBuildDiagonalPfloatKernel(mesh->Nelements,
                                                     elliptic->Nfields,
                                                     elliptic->fieldOffset,
                                                     elliptic->loffset,
                                                     mesh->o_ggeoPfloat,
                                                     mesh->o_DPfloat,
                                                     mesh->o_DTPfloat,
                                                     elliptic->o_lambda0,
                                                     elliptic->o_lambda1,
                                                     o_invDiagA);
    flopCount += 12 * mesh->Nq + 12;
    flopCount += (elliptic->poisson) ? 0.0 : 2.0;
    flopCount *= static_cast<double>(mesh->Nlocal) * elliptic->Nfields;
    flopCount *= 0.5;

    oogs::startFinish(o_invDiagA,
                      elliptic->Nfields,
                      elliptic->fieldOffset,
                      ogsPfloat,
                      ogsAdd,
                      elliptic->oogs);

    const pfloat one = 1.0;
    platform->linAlg->padyMany(mesh->Nlocal, elliptic->Nfields, elliptic->fieldOffset, one, o_invDiagA);
  } else {
    elliptic->ellipticBlockBuildDiagonalKernel(mesh->Nelements,
                                               elliptic->Nfields,
                                               elliptic->fieldOffset,
                                               elliptic->loffset,
                                               mesh->o_ggeo,
                                               mesh->o_D,
                                               mesh->o_DT,
                                               elliptic->o_lambda0,
                                               elliptic->o_lambda1,
                                               o_invDiagA);

    flopCount += 12 * mesh->Nq + 12;
    flopCount += (elliptic->poisson) ? 0.0 : 2.0;
    flopCount *= static_cast<double>(mesh->Nlocal) * elliptic->Nfields;

    oogs::startFinish(o_invDiagA,
                      elliptic->Nfields,
                      elliptic->fieldOffset,
                      ogsDfloat,
                      ogsAdd,
                      elliptic->oogs);

    const dfloat one = 1.0;
    platform->linAlg->adyMany(mesh->Nlocal, elliptic->Nfields, elliptic->fieldOffset, one, o_invDiagA);
  }

  platform->flopCounter->add(elliptic->name + " ellipticUpdateJacobi", flopCount);
}

void ellipticUpdateJacobi(elliptic_t *ellipticBase)
{
  setupAide& options = ellipticBase->options;

  if(options.compareArgs("PRECONDITIONER", "MULTIGRID") &&
     options.compareArgs("MULTIGRID SMOOTHER", "DAMPEDJACOBI")) {
    precon_t *precon = ellipticBase->precon;
    MGSolver_t::multigridLevel** levels = precon->MGSolver->levels;

    for(int levelIndex = 0; levelIndex < ellipticBase->nLevels; levelIndex++) {
      auto mgLevel = dynamic_cast<pMGLevel*>(levels[levelIndex]);
      auto elliptic = mgLevel->elliptic;
      auto mesh = elliptic->mesh;
  
      const bool coarsestLevel = (levelIndex == ellipticBase->nLevels - 1);
      if(coarsestLevel && elliptic->options.compareArgs("MULTIGRID COARSE SOLVE", "TRUE"))
        continue;

      ellipticUpdateJacobi(elliptic, mgLevel->o_invDiagA); 
    }
  } else if(options.compareArgs("PRECONDITIONER", "JACOBI")) {
    precon_t *precon = ellipticBase->precon;
    ellipticUpdateJacobi(ellipticBase, precon->o_invDiagA); 
  }
}
