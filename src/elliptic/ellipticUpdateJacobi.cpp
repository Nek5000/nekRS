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
#include "linAlg.hpp"

void ellipticUpdateJacobi(elliptic_t *elliptic, occa::memory &o_invDiagA)
{
  dfloat flopCount = 0.0;
  mesh_t *mesh = elliptic->mesh;
  setupAide& options = elliptic->options;

  const dlong Nlocal = mesh->Np * mesh->Nelements;

  if(elliptic->mgLevel)
    elliptic->ellipticBlockBuildDiagonalPfloatKernel(mesh->Nelements,
                                                     elliptic->Nfields,
                                                     elliptic->fieldOffset,
                                                     elliptic->loffset,
                                                     mesh->o_ggeoPfloat,
                                                     mesh->o_DPfloat,
                                                     mesh->o_DTPfloat,
                                                     elliptic->o_lambda,
                                                     o_invDiagA);
  else
    elliptic->ellipticBlockBuildDiagonalKernel(mesh->Nelements,
                                               elliptic->Nfields,
                                               elliptic->fieldOffset,
                                               elliptic->loffset,
                                               mesh->o_ggeo,
                                               mesh->o_D,
                                               mesh->o_DT,
                                               elliptic->o_lambda,
                                               o_invDiagA /* pfloat */);

  flopCount += 12 * mesh->Nq + 12;
  flopCount += (elliptic->poisson) ? 0.0 : 2.0;
  flopCount *= static_cast<double>(mesh->Nlocal) * elliptic->Nfields;
  if(elliptic->mgLevel) flopCount *= 0.5;

  oogs::startFinish(o_invDiagA, elliptic->Nfields, elliptic->fieldOffset, ogsPfloat, ogsAdd, elliptic->oogs);

  const pfloat one = 1.0;
  platform->linAlg->padyMany(Nlocal, elliptic->Nfields, elliptic->fieldOffset, one, o_invDiagA);
  platform->flopCounter->add(elliptic->name + " ellipticUpdateJacobi", flopCount);
}
