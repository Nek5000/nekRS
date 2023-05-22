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

#include "mesh.h"
//#include "ogsInterface.h"

void meshFree(mesh_t* mesh)
{
  if(mesh->EX) free(mesh->EX);   // coordinates of vertices for each element
  if(mesh->EY) free(mesh->EY);
  if(mesh->EZ) free(mesh->EZ);

  if(mesh->EToV) free(mesh->EToV);   // element-to-vertex connectivity
  if(mesh->EToE) free(mesh->EToE);   // element-to-element connectivity
  if(mesh->EToF) free(mesh->EToF);   // element-to-(local)face connectivity
  if(mesh->EToP) free(mesh->EToP);   // element-to-partition/process connectivity
  if(mesh->EToB) free(mesh->EToB);   // element-to-boundary condition type

  if(mesh->elementInfo) free(mesh->elementInfo);   //type of element

  // MPI halo exchange info
  if(mesh->haloElementList) free(mesh->haloElementList);   // sorted list of elements to be sent in halo exchange
  if(mesh->NhaloPairs) free(mesh->NhaloPairs);        // number of elements worth of data to send/recv

  if(mesh->haloGetNodeIds) free(mesh->haloGetNodeIds);   // volume node ids of outgoing halo nodes
  if(mesh->haloPutNodeIds) free(mesh->haloPutNodeIds);   // volume node ids of incoming halo nodes

  if(mesh->haloSendRequests) free(mesh->haloSendRequests);
  if(mesh->haloRecvRequests) free(mesh->haloRecvRequests);

  // CG gather-scatter info
  if(mesh->globalIds) free(mesh->globalIds);
  //if(mesh->gsh) ogsHostFree(mesh->gsh);
  //if(  mesh->hostGsh) ogsHostFree(  mesh->hostGsh);// gslib struct pointer
  if(  mesh->ogs) ogsFree(  mesh->ogs); //occa gs pointer

  // list of elements that are needed for global gather-scatter
  if(mesh->globalGatherElementList) free(mesh->globalGatherElementList);

  // list of elements that are not needed for global gather-scatter
  if(mesh->localGatherElementList) free(mesh->localGatherElementList);

  // volume node info
  if(mesh->r) free(mesh->r);
  if(mesh->s) free(mesh->s);
  if(mesh->t) free(mesh->t);      // coordinates of local nodes
  if(mesh->x) free(mesh->x);
  if(mesh->y) free(mesh->y);
  if(mesh->z) free(mesh->z);      // coordinates of physical nodes
  if(mesh->vertexNodes) free(mesh->vertexNodes);

  if(mesh->D) free(mesh->D);   // 1D differentiation matrix (for tensor-product)
  if(mesh->gllz) free(mesh->gllz);   // 1D GLL quadrature nodes
  if(mesh->gllw) free(mesh->gllw);   // 1D GLL quadrature weights

  // face node info
  if(mesh->faceNodes) free(mesh->faceNodes);   // list of element reference interpolation nodes on element faces
  if (mesh->vmapM)
    free(mesh->vmapM);                               // list of volume nodes that are face nodes
  if(mesh->faceVertices) free(mesh->faceVertices);   // list of mesh vertices on each face

  if(mesh->cubr) free(mesh->cubr);
  if(mesh->cubs) free(mesh->cubs);
  if(mesh->cubt) free(mesh->cubt);
  if(mesh->cubw) free(mesh->cubw);   // coordinates and weights of local cubature nodes
  if(mesh->cubx) free(mesh->cubx);
  if(mesh->cuby) free(mesh->cuby);
  if(mesh->cubz) free(mesh->cubz);      // coordinates of physical nodes
  if(mesh->cubInterp) free(mesh->cubInterp);   // interpolate from W&B to cubature nodes
  if(mesh->cubProject) free(mesh->cubProject);   // projection matrix from cubature nodes to W&B nodes
  if(mesh->cubD) free(mesh->cubD);         // 1D differentiation matrix
  if(mesh->cubDiffInterp) free(mesh->cubDiffInterp);       // 1D weak differentiation matrix
  if(mesh->cubDW) free(mesh->cubDW);       // 1D weak differentiation matrix
  if(mesh->cubDWmatrices) free(mesh->cubDWmatrices);

  // surface integration node info
  if(mesh->intInterp) free(mesh->intInterp);   // interp from surface node to integration nodes
  if(mesh->intx) free(mesh->intx);
  if(mesh->inty) free(mesh->inty);
  if(mesh->intz) free(mesh->intz);   // coordinates of suface integration nodes

  //degree raising and lowering interpolation matrices
  if(mesh->interpRaise) free(mesh->interpRaise);
  if(mesh->interpLower) free(mesh->interpLower);

}
