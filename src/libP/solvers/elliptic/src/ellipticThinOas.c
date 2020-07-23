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

// z = P*r
void ellipticThinOas(elliptic_t* elliptic, dfloat lambda, occa::memory &o_r, occa::memory &o_z)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  dlong Nelements = mesh->Nelements;
  dlong totalHaloPairs = mesh->totalHaloPairs;
  dlong Ndata = mesh->Nfp * totalHaloPairs;

  // extract halo for thin kernel
  elliptic->oasHaloGetKernel(totalHaloPairs,
                             1, // number of fields
                             Np * (Nelements + totalHaloPairs), // field size including halo eelments
                             elliptic->o_oasHaloElementList,
                             elliptic->o_oasHaloGetNodeIds,
                             o_r,
                             elliptic->o_oasHaloBuffer);

  // copy extracted halo to HOST
  elliptic->o_oasHaloBuffer.copyTo(elliptic->oasSendBuffer, Ndata * sizeof(dfloat), 0);// zero offset

  // start halo exchange
  meshHaloExchangeStart(mesh,
                        mesh->Nfp * sizeof(dfloat),
                        elliptic->oasSendBuffer,
                        elliptic->oasRecvBuffer);

  // finish halo exchange
  meshHaloExchangeFinish(mesh);

  // copy halo data to on device halo arary
  elliptic->o_oasHaloBuffer.copyFrom(elliptic->oasRecvBuffer, Ndata * sizeof(dfloat), 0);  // zero offset

  // populate halo (with offset
  elliptic->oasHaloPutKernel(totalHaloPairs,
                             1, // number of fields
                             Np * (Nelements + totalHaloPairs), // field size including halo eelments
                             elliptic->o_oasHaloElementList,
                             elliptic->o_oasHaloPutNodeIds,
                             elliptic->o_oasHaloBuffer,
                             o_r); // place incoming halo data into end buffer

  // element-wise fast-directional-approximate-inverse
  elliptic->oasPreconditionerKernel(Nelements,
                                    elliptic->o_oasMapP,
                                    elliptic->o_oasForward,
                                    elliptic->o_oasBack,
                                    elliptic->o_oasDiagInvOp,
                                    o_r,
                                    elliptic->o_oasTmp); // need oasNq^3*Nelements

  // gather scatter to o_z need to stack
  ogsGatherScatter(elliptic->o_oasTmp, ogsDfloat, ogsAdd, elliptic->oasOgs);

  // copy back to o_z
  o_z.copyFrom(elliptic->o_oasTmp, mesh->Np * mesh->Nelements * sizeof(dfloat),0); // 0 offset
}
