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
#include "ogsInterface.h"

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

  // boundary faces
  if(mesh->boundaryInfo) free(mesh->boundaryInfo);   // list of boundary faces (type, vertex-1, vertex-2, vertex-3)

  // MPI halo exchange info
  if(mesh->haloElementList) free(mesh->haloElementList);   // sorted list of elements to be sent in halo exchange
  if(mesh->NhaloPairs) free(mesh->NhaloPairs);        // number of elements worth of data to send/recv

  if(mesh->haloGetNodeIds) free(mesh->haloGetNodeIds);   // volume node ids of outgoing halo nodes
  if(mesh->haloPutNodeIds) free(mesh->haloPutNodeIds);   // volume node ids of incoming halo nodes

  if(mesh->haloSendRequests) free(mesh->haloSendRequests);
  if(mesh->haloRecvRequests) free(mesh->haloRecvRequests);

  // CG gather-scatter info
  if(mesh->globalIds) free(mesh->globalIds);
  if(mesh->maskedGlobalIds) free(mesh->maskedGlobalIds);
  if(  mesh->gsh) ogsHostFree(  mesh->gsh);
  if(  mesh->hostGsh) ogsHostFree(  mesh->hostGsh);// gslib struct pointer
  if(  mesh->ogs) ogsFree(  mesh->ogs);//occa gs pointer

  // list of elements that are needed for global gather-scatter
  if(mesh->globalGatherElementList) free(mesh->globalGatherElementList);

  // list of elements that are not needed for global gather-scatter
  if(mesh->localGatherElementList) free(mesh->localGatherElementList);

  //list of fair pairs
  if(mesh->EToFPairs) free(mesh->EToFPairs);
  if(mesh->FPairsToE) free(mesh->FPairsToE);
  if(mesh->FPairsToF) free(mesh->FPairsToF);

  if(mesh->vgeo) free(mesh->vgeo);

  if(mesh->ggeo) free(mesh->ggeo);

  // volume node info
  if(mesh->r) free(mesh->r);
  if(mesh->s) free(mesh->s);
  if(mesh->t) free(mesh->t);      // coordinates of local nodes
  if(mesh->Dr) free(mesh->Dr);
  if(mesh->Ds) free(mesh->Ds);
  if(mesh->Dt) free(mesh->Dt);   // collocation differentiation matrices
  if(mesh->Dmatrices) free(mesh->Dmatrices);
  if(mesh->MM) free(mesh->MM);
  if(mesh->invMM) free(mesh->invMM);             // reference mass matrix
  if(mesh->invLMM) free(mesh->invLMM);
  if(mesh->Srr) free(mesh->Srr);
  if(mesh->Srs) free(mesh->Srs);
  if(mesh->Srt) free(mesh->Srt);   //element stiffness matrices
  if(mesh->Ssr) free(mesh->Ssr);
  if(mesh->Sss) free(mesh->Sss);
  if(mesh->Sst) free(mesh->Sst);
  if(mesh->Str) free(mesh->Str);
  if(mesh->Sts) free(mesh->Sts);
  if(mesh->Stt) free(mesh->Stt);
  if(mesh->Smatrices) free(mesh->Smatrices);
  if(mesh->x) free(mesh->x);
  if(mesh->y) free(mesh->y);
  if(mesh->z) free(mesh->z);      // coordinates of physical nodes
  if(mesh->vertexNodes) free(mesh->vertexNodes);

  if(mesh->D) free(mesh->D);   // 1D differentiation matrix (for tensor-product)
  if(mesh->gllz) free(mesh->gllz);   // 1D GLL quadrature nodes
  if(mesh->gllw) free(mesh->gllw);   // 1D GLL quadrature weights

  if(mesh->gjr) free(mesh->gjr);
  if(mesh->gjw) free(mesh->gjw);   // 1D nodes and weights for Gauss Jacobi quadature
  if(mesh->gjI) free(mesh->gjI);
  if(mesh->gjD) free(mesh->gjD);   // 1D GLL to Gauss node interpolation and differentiation matrices
  if(mesh->gjD2) free(mesh->gjD2);       // 1D GJ to GJ node differentiation

  if(mesh->oasForward) free(mesh->oasForward);
  if(mesh->oasBack) free(mesh->oasBack);
  if(mesh->oasDiagOp) free(mesh->oasDiagOp);

  if(mesh->oasForwardDg) free(mesh->oasForwardDg);
  if(mesh->oasBackDg) free(mesh->oasBackDg);
  if(mesh->oasDiagOpDg) free(mesh->oasDiagOpDg);

  if(mesh->rmapP) free(mesh->rmapP);

  //reference patch inverse (for OAS precon)
  if(mesh->invAP) free(mesh->invAP);

  // face node info
  if(mesh->faceNodes) free(mesh->faceNodes);   // list of element reference interpolation nodes on element faces
  if(mesh->vmapM) free(mesh->vmapM);       // list of volume nodes that are face nodes
  if(mesh->vmapP) free(mesh->vmapP);       // list of volume nodes that are paired with face nodes
  if(mesh->mapP) free(mesh->mapP);       // list of surface nodes that are paired with -ve surface  nodes
  if(mesh->faceVertices) free(mesh->faceVertices);   // list of mesh vertices on each face

  if(mesh->LIFT) free(mesh->LIFT);   // lift matrix
  if(mesh->FMM) free(mesh->FMM);    // Face Mass Matrix
  if(mesh->sMT) free(mesh->sMT);   // surface mass (MM*LIFT)^T

  if(mesh->sgeo) free(mesh->sgeo);

  // field info for PDE solver
  if(mesh->q) free(mesh->q);      // solution data array
  if(mesh->fQM) free(mesh->fQM);
  if(mesh->fQP) free(mesh->fQP);   //solution trace arrays
  if(mesh->rhsq) free(mesh->rhsq);
  if(mesh->rhsq2) free(mesh->rhsq2);
  if(mesh->rhsq3) free(mesh->rhsq3);   // right hand side data array
  if(mesh->resq) free(mesh->resq);   // residual data array (for LSERK time-stepping)

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
  if(mesh->cubDrW) free(mesh->cubDrW);      // 'r' weak differentiation matrix
  if(mesh->cubDsW) free(mesh->cubDsW);      // 's' weak differentiation matrix
  if(mesh->cubDtW) free(mesh->cubDtW);      // 't' weak differentiation matrix
  if(mesh->cubDWmatrices) free(mesh->cubDWmatrices);

  if(mesh->cubvgeo) free(mesh->cubvgeo);    //volume geometric data at cubature points
  if(mesh->cubsgeo) free(mesh->cubsgeo);    //surface geometric data at cubature points
  if(mesh->cubggeo) free(mesh->cubggeo);    //second type volume geometric data at cubature points

  // c2 at cubature points (for wadg)
  if(mesh->c2) free(mesh->c2);

  //source injection
  if(mesh->sourceq) free(mesh->sourceq);
  if(mesh->MRABsourceNelements) free(mesh->MRABsourceNelements);
  if(mesh->sourceElements) free(mesh->sourceElements);

  // surface integration node info
  if(mesh->intInterp) free(mesh->intInterp);   // interp from surface node to integration nodes
  if(mesh->intLIFT) free(mesh->intLIFT);     // lift from surface integration nodes to W&B volume nodes
  if(mesh->intx) free(mesh->intx);
  if(mesh->inty) free(mesh->inty);
  if(mesh->intz) free(mesh->intz);   // coordinates of suface integration nodes

  // Bernstein-Bezier info
  if(mesh->VB) free(mesh->VB);
  if(mesh->invVB) free(mesh->invVB);   // Bernstein Vandermonde matrices
  if(mesh->BBMM) free(mesh->BBMM);
  if(mesh->invVB1D) free(mesh->invVB1D);
  if(mesh->invVB2D) free(mesh->invVB2D);
  if(mesh->D0ids) free(mesh->D0ids);
  if(mesh->D1ids) free(mesh->D1ids);
  if(mesh->D2ids) free(mesh->D2ids);
  if(mesh->D3ids) free(mesh->D3ids);   // Bernstein deriv matrix indices
  if(mesh->Dvals) free(mesh->Dvals);   // Bernstein deriv matrix values
  if(mesh->D0Tids) free(mesh->D0Tids);
  if(mesh->D1Tids) free(mesh->D1Tids);
  if(mesh->D2Tids) free(mesh->D2Tids);
  if(mesh->D3Tids) free(mesh->D3Tids);   // Bernstein transpose deriv matrix indices
  if(mesh->DTvals) free(mesh->DTvals);   // Bernstein transpose deriv matrix values
  if(mesh->VBq) free(mesh->VBq);
  if(mesh->PBq) free(mesh->PBq);   // cubature interpolation/projection matrices
  if(mesh->L0ids) free(mesh->L0ids);   // L0 matrix ids
  if(mesh->L0vals) free(mesh->L0vals);   // L0 values (L0 tridiagonal in 2D)
  if(mesh->ELids) free(mesh->ELids);   // lift reduction matrix indices
  if(mesh->ELvals) free(mesh->ELvals);   // lift reduction matrix values
  if(mesh->BBRaiseids) free(mesh->BBRaiseids);   //Bernstein elevate matrix indices
  if(mesh->BBRaiseVals) free(mesh->BBRaiseVals);   //Bernstein elevate matrix values
  if(mesh->BBLower) free(mesh->BBLower);   //Berstein projection matrix.

  //degree raising and lowering interpolation matrices
  if(mesh->interpRaise) free(mesh->interpRaise);
  if(mesh->interpLower) free(mesh->interpLower);

  //sparse basis info
  if(mesh->sparseV) free(mesh->sparseV);
  if(mesh->invSparseV) free(mesh->invSparseV);
  if(mesh->invSparseV) free(mesh->invSparseV);
  if(mesh->sparseMM) free(mesh->sparseMM);
  if(mesh->FaceModes) free(mesh->FaceModes);
  if(mesh->sparseStackedNZ) free(mesh->sparseStackedNZ);
  if(mesh->sparseSrrT) free(mesh->sparseSrrT);
  if(mesh->sparseSrsT) free(mesh->sparseSrsT);
  if(mesh->sparseSssT) free(mesh->sparseSssT);
  if(mesh->Ind) free(mesh->Ind);

  if(mesh->mmapM) free(mesh->mmapM);
  if(mesh->mmapP) free(mesh->mmapP);
  if(mesh->mmapS) free(mesh->mmapS);
  if(mesh->mapSgn) free(mesh->mapSgn);

  if(mesh->MRABlevel) free(mesh->MRABlevel);
  if(mesh->MRABNelements) free(mesh->MRABNelements);
  if(mesh->MRABNhaloElements) free(mesh->MRABNhaloElements);
  if(mesh->MRABelementIds) free(mesh->MRABelementIds);
  if(mesh->MRABhaloIds) free(mesh->MRABhaloIds);
  if(mesh->MRABshiftIndex) free(mesh->MRABshiftIndex);

  if(mesh->MRABpmlNelements) free(mesh->MRABpmlNelements);
  if(mesh->MRABpmlNhaloElements) free(mesh->MRABpmlNhaloElements);
  if(mesh->MRABpmlElementIds) free(mesh->MRABpmlElementIds);
  if(mesh->MRABpmlIds) free(mesh->MRABpmlIds);
  if(mesh->MRABpmlHaloElementIds) free(mesh->MRABpmlHaloElementIds);
  if(mesh->MRABpmlHaloIds) free(mesh->MRABpmlHaloIds);

  if(mesh->errtmp) free(mesh->errtmp);

  if(mesh->plotEToV) free(mesh->plotEToV);        // triangulation of plot nodes
  if(mesh->plotR) free(mesh->plotR);
  if(mesh->plotS) free(mesh->plotS);
  if(mesh->plotT) free(mesh->plotT);   // coordinates of plot nodes in reference element
  if(mesh->plotInterp) free(mesh->plotInterp);      // warp & blend to plot node interpolation matrix

  if(mesh->contourEToV) free(mesh->contourEToV);
  if(mesh->contourVX) free(mesh->contourVX);
  if(mesh->contourVY) free(mesh->contourVY);
  if(mesh->contourVZ) free(mesh->contourVZ);
  if(mesh->contourInterp) free(mesh->contourInterp);
  if(mesh->contourInterp1) free(mesh->contourInterp1);
  if(mesh->contourFilter) free(mesh->contourFilter);

  //SEMFEM data
  if(mesh->FEMEToV) free(mesh->FEMEToV);
  if(mesh->rFEM) free(mesh->rFEM);
  if(mesh->sFEM) free(mesh->sFEM);
  if(mesh->tFEM) free(mesh->tFEM);
  if(mesh->SEMFEMInterp) free(mesh->SEMFEMInterp);

  if(mesh->pmlElementList) free(mesh->pmlElementList);   // deprecated

  if(mesh->invTau) free(mesh->invTau);   // deprecated in Boltzmann
  if(mesh->probeR) free(mesh->probeR);
  if(mesh->probeS) free(mesh->probeS);
  if(mesh->probeT) free(mesh->probeT);
  if(mesh->probeElementIds) free(mesh->probeElementIds);
  if(mesh->probeIds) free(mesh->probeIds);
  if(mesh->probeI) free(mesh->probeI);
}
