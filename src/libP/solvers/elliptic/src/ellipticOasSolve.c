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

void ellipticOasSolve(elliptic_t* elliptic, dfloat lambda,
                      occa::memory &o_r, occa::memory &o_z)
{
  precon_t* precon = elliptic->precon;
  mesh_t* mesh = elliptic->mesh;

  elliptic_t* elliptic1 = (elliptic_t*) precon->ellipticOneRing; // should rename
  mesh_t* mesh1 = elliptic1->mesh;

  elliptic_t* ellipticOasCoarse = (elliptic_t*) (precon->ellipticOasCoarse);
  mesh_t* meshCoarse   = ellipticOasCoarse->mesh;

  // TW: possibility these device have difference queues
  mesh1->device.finish();
  mesh->device.finish();

  // 1. collect patch residual
  ellipticOneRingExchange(elliptic, elliptic1, mesh1->Np * sizeof(dfloat), o_r, elliptic1->o_r);

  // mask local overlapping patch residual
  if (elliptic1->Nmasked)
    mesh1->maskKernel(elliptic1->Nmasked, elliptic1->o_maskIds, elliptic1->o_r);

  // TW: possibility these device have difference queues
  mesh1->device.finish();
  mesh->device.finish();

  ellipticPreconditioner(elliptic1, lambda, elliptic1->o_r, elliptic1->o_x); // may need to zero o_x

  // gathering over all patches - so have to remove local multiplicity
  elliptic1->dotMultiplyKernel(mesh1->Nelements * mesh1->Np,
                               elliptic1->ogs->o_invDegree,
                               elliptic1->o_x,
                               elliptic1->o_z);

  // sum up overlapping patches
  ogsGatherScatter(elliptic1->o_z, ogsDfloat, ogsAdd, elliptic->precon->oasOgs);

  o_z.copyFrom(elliptic1->o_z, mesh->Nelements * mesh->Np * sizeof(dfloat), 0);

  // 2. solve coarse problem
  //   a. call solver
  //   b. prolongate (watch out for +=)

  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();

  precon->oasRestrictionKernel(meshCoarse->Nelements,
                               precon->o_oasRestrictionMatrix,
                               o_r, ellipticOasCoarse->o_r);

  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();

  // why do I Have to do (1/deg)*S*G*o_rCoarse here ? ---------->
  ogsGatherScatter(ellipticOasCoarse->o_r, ogsDfloat, ogsAdd, ellipticOasCoarse->ogs);

  ellipticOasCoarse->dotMultiplyKernel(meshCoarse->Nelements * meshCoarse->Np,
                                       meshCoarse->ogs->o_invDegree,
                                       ellipticOasCoarse->o_r,
                                       ellipticOasCoarse->o_r);

  if (ellipticOasCoarse->Nmasked)
    meshCoarse->maskKernel(ellipticOasCoarse->Nmasked,
                           ellipticOasCoarse->o_maskIds,
                           ellipticOasCoarse->o_r);

  // <----------

  ellipticPreconditioner(ellipticOasCoarse, lambda, ellipticOasCoarse->o_r, ellipticOasCoarse->o_x);

  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();

  // prolongate to QN (note kernel expects restriction matrix)
  // do we need to weight the sum against patches?
  precon->oasProlongationKernel(mesh->Nelements, precon->o_oasRestrictionMatrix,
                                ellipticOasCoarse->o_x, o_z);

  mesh1->device.finish();
  mesh->device.finish();
  meshCoarse->device.finish();

  if (elliptic->Nmasked)
    mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_z);
}
