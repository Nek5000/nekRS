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

//returns the ipdg patch A matrix for element eM
void BuildLocalIpdgPatchAxTri2D(elliptic_t* elliptic,
                                mesh_t* mesh,
                                int basisNp,
                                dfloat* basis,
                                dfloat lambda,
                                dfloat* MS,
                                dlong eM,
                                dfloat* A);
void BuildLocalIpdgPatchAxQuad2D(elliptic_t* elliptic,
                                 mesh_t* mesh,
                                 dfloat lambda,
                                 dfloat* B,
                                 dfloat* Br,
                                 dfloat* Bs,
                                 dlong eM,
                                 dfloat* A);
void BuildLocalIpdgPatchAxTet3D(elliptic_t* elliptic,
                                mesh_t* mesh,
                                dfloat lambda,
                                dfloat* MS,
                                dlong eM,
                                dfloat* A);
void BuildLocalIpdgPatchAxHex3D(elliptic_t* elliptic,
                                mesh_t* mesh,
                                dfloat lambda,
                                dfloat* B,
                                dfloat* Br,
                                dfloat* Bs,
                                dfloat* Bt,
                                dlong eM,
                                dfloat* A);

//returns the C0FEM patch A matrix for element eM
void BuildLocalContinuousPatchAxTri2D(elliptic_t* elliptic,
                                      mesh_t* mesh,
                                      dfloat lambda,
                                      dlong eM,
                                      dfloat* A);
void BuildLocalContinuousPatchAxQuad2D(elliptic_t* elliptic,
                                       mesh_t* mesh,
                                       dfloat lambda,
                                       dlong eM,
                                       dfloat* B,
                                       dfloat* Br,
                                       dfloat* Bs,
                                       dfloat* A);
void BuildLocalContinuousPatchAxTet3D(elliptic_t* elliptic,
                                      mesh_t* mesh,
                                      dfloat lambda,
                                      dlong eM,
                                      dfloat* A);
void BuildLocalContinuousPatchAxHex3D(elliptic_t* elliptic,
                                      mesh_t* mesh,
                                      dfloat lambda,
                                      dlong eM,
                                      dfloat* B,
                                      dfloat* Br,
                                      dfloat* Bs,
                                      dfloat* Bt,
                                      dfloat* A);

void ellipticBuildLocalPatchesTri2D(elliptic_t* elliptic, dfloat lambda, dfloat rateTolerance,
                                    dlong* Npatches, dlong** patchesIndex, dfloat** patchesInvA);
void ellipticBuildLocalPatchesQuad2D(elliptic_t* elliptic, dfloat lambda, dfloat rateTolerance,
                                     dlong* Npatches, dlong** patchesIndex, dfloat** patchesInvA);
void ellipticBuildLocalPatchesTet3D(elliptic_t* elliptic, dfloat lambda, dfloat rateTolerance,
                                    dlong* Npatches, dlong** patchesIndex, dfloat** patchesInvA);
void ellipticBuildLocalPatchesHex3D(elliptic_t* elliptic, dfloat lambda, dfloat rateTolerance,
                                    dlong* Npatches, dlong** patchesIndex, dfloat** patchesInvA);

void ellipticBuildLocalPatches(elliptic_t* elliptic, dfloat lambda, dfloat rateTolerance,
                               dlong* Npatches, dlong** patchesIndex, dfloat** patchesInvA)
{
  switch(elliptic->elementType) {
  case TRIANGLES:
    ellipticBuildLocalPatchesTri2D(elliptic,
                                   lambda,
                                   rateTolerance,
                                   Npatches,
                                   patchesIndex,
                                   patchesInvA);
    break;
  case QUADRILATERALS:
    ellipticBuildLocalPatchesQuad2D(elliptic,
                                    lambda,
                                    rateTolerance,
                                    Npatches,
                                    patchesIndex,
                                    patchesInvA);
    break;
  case TETRAHEDRA:
    ellipticBuildLocalPatchesTet3D(elliptic,
                                   lambda,
                                   rateTolerance,
                                   Npatches,
                                   patchesIndex,
                                   patchesInvA);
    break;
  case HEXAHEDRA:
    ellipticBuildLocalPatchesHex3D(elliptic,
                                   lambda,
                                   rateTolerance,
                                   Npatches,
                                   patchesIndex,
                                   patchesInvA);
    break;
  }
}

void ellipticBuildLocalPatchesTri2D(elliptic_t* elliptic, dfloat lambda, dfloat rateTolerance,
                                    dlong* Npatches, dlong** patchesIndex, dfloat** patchesInvA)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  // surface mass matrices MS = MM*LIFT
  dfloat* MS = (dfloat*) calloc(mesh->Nfaces * mesh->Nfp * mesh->Nfp,sizeof(dfloat));
  for (int f = 0; f < mesh->Nfaces; f++)
    for (int n = 0; n < mesh->Nfp; n++) {
      int fn = mesh->faceNodes[f * mesh->Nfp + n];

      for (int m = 0; m < mesh->Nfp; m++) {
        dfloat MSnm = 0;

        for (int i = 0; i < mesh->Np; i++)
          MSnm += mesh->MM[fn + i * mesh->Np] *
                  mesh->LIFT[i * mesh->Nfp * mesh->Nfaces + f * mesh->Nfp + m];

        MS[m + n * mesh->Nfp + f * mesh->Nfp * mesh->Nfp]  = MSnm;
      }
    }

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  *patchesIndex = (dlong*) calloc(mesh->Nelements, sizeof(dlong));

  //temp patch storage
  dfloat* patchA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* invRefAA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));

  (*Npatches) = 1;
  dlong refPatches = 0;

  //build a mini mesh struct for the reference patch
  mesh_t* refMesh = (mesh_t*) calloc(1,sizeof(mesh_t));
  memcpy(refMesh,mesh,sizeof(mesh_t));

  //vertices of reference patch
  dfloat V1x = -1., V2x = 1., V3x =        0.;
  dfloat V1y =  0., V2y = 0., V3y =  sqrt(3.);

  refMesh->Nelements = 1;

  refMesh->EX = (dfloat*) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EY = (dfloat*) calloc(mesh->Nverts,sizeof(dfloat));

  refMesh->EX[0] = V1x;
  refMesh->EY[0] = V1y;
  refMesh->EX[1] = V2x;
  refMesh->EY[1] = V2y;
  refMesh->EX[2] = V3x;
  refMesh->EY[2] = V3y;

  refMesh->EToV = (hlong*) calloc(mesh->Nverts, sizeof(hlong));

  refMesh->EToV[0] = 0;
  refMesh->EToV[1] = 1;
  refMesh->EToV[2] = 2;

  refMesh->EToB = (int*) calloc(mesh->Nfaces,sizeof(int));
  for (int n = 0; n < mesh->Nfaces; n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesTri2D(refMesh, mesh->N);
  meshPhysicalNodesTri2D(refMesh);
  meshGeometricFactorsTri2D(refMesh);
  meshConnectFaceNodes2D(refMesh);
  meshSurfaceGeometricFactorsTri2D(refMesh);

  int basisNp = mesh->Np;
  dfloat* basis;
  if(options.compareArgs("BASIS","BERN")) {
    basis = mesh->VB;
  } else {// default to degree N Lagrange basis
    basis = (dfloat*) calloc(basisNp * basisNp, sizeof(dfloat));
    for(int n = 0; n < basisNp; ++n)
      basis[n + n * basisNp] = 1;
  }

  //start with reference patch
  dfloat* refPatchInvA = *patchesInvA;
  if (options.compareArgs("DISCRETIZATION","IPDG"))
    BuildLocalIpdgPatchAxTri2D(elliptic, refMesh, basisNp, basis, lambda, MS, 0, refPatchInvA);
  else if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
    BuildLocalContinuousPatchAxTri2D(elliptic, refMesh, lambda, 0, refPatchInvA);

  matrixInverse(mesh->Np, refPatchInvA);

  // loop over all elements
  for(dlong eM = 0; eM < mesh->Nelements; ++eM) {
    //build the patch A matrix for this element
    if (options.compareArgs("DISCRETIZATION","IPDG"))
      BuildLocalIpdgPatchAxTri2D(elliptic, mesh, basisNp, basis, lambda, MS, eM, patchA);
    else if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
      BuildLocalContinuousPatchAxTri2D(elliptic, mesh, lambda, eM, refPatchInvA);

    dlong eP0 = mesh->EToE[eM * mesh->Nfaces + 0];
    dlong eP1 = mesh->EToE[eM * mesh->Nfaces + 1];
    dlong eP2 = mesh->EToE[eM * mesh->Nfaces + 2];

    if(eP0 >= 0 && eP1 >= 0 && eP2 >= 0) { //check if this is an interior patch
      refPatchInvA = *patchesInvA;

      //hit the patch with the reference inverse
      for(int n = 0; n < mesh->Np; ++n)
        for(int m = 0; m < mesh->Np; ++m) {
          invRefAA[n * mesh->Np + m] = 0.;
          for (int k = 0; k < mesh->Np; k++)
            invRefAA[n * mesh->Np + m] += refPatchInvA[n * mesh->Np + k] * patchA[k * mesh->Np + m];
        }

      dfloat cond = matrixConditionNumber(mesh->Np,invRefAA);
      dfloat rate = (sqrt(cond) - 1.) / (sqrt(cond) + 1.);

      // printf("Element %d's conditioned patch reports cond = %g and rate = %g \n", eM, cond, rate);

      if (rate < rateTolerance) {
        (*patchesIndex)[eM] = 0;
        refPatches++;
        continue;
      }
    }
    ++(*Npatches);
    *patchesInvA = (dfloat*) realloc(*patchesInvA,
                                     (*Npatches) * mesh->Np * mesh->Np * sizeof(dfloat));

    matrixInverse(mesh->Np, patchA);

    //copy inverse into patchesInvA
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        dlong id = ((*Npatches) - 1) * mesh->Np * mesh->Np + n * mesh->Np + m;
        (*patchesInvA)[id] = patchA[n * mesh->Np + m];
      }

    (*patchesIndex)[eM] = (*Npatches) - 1;
  }

  printf("saving " dlongFormat " full patches\n",*Npatches);
  printf("using " dlongFormat " reference patches\n", refPatches);

  free(refMesh);
  free(patchA);
  free(invRefAA);
  free(MS);
}

void ellipticBuildLocalPatchesQuad2D(elliptic_t* elliptic, dfloat lambda, dfloat rateTolerance,
                                     dlong* Npatches, dlong** patchesIndex, dfloat** patchesInvA)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  // build some monolithic basis arrays
  dfloat* B  = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* Br = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* Bs = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nj = 0; nj < mesh->N + 1; ++nj)
    for(int ni = 0; ni < mesh->N + 1; ++ni) {
      int node = 0;

      for(int j = 0; j < mesh->N + 1; ++j)
        for(int i = 0; i < mesh->N + 1; ++i) {
          if(nj == j && ni == i)
            B[mode * mesh->Np + node] = 1;
          if(nj == j)
            Br[mode * mesh->Np + node] = mesh->D[ni + mesh->Nq * i];
          if(ni == i)
            Bs[mode * mesh->Np + node] = mesh->D[nj + mesh->Nq * j];

          ++node;
        }
      ++mode;
    }

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  *patchesIndex = (dlong*) calloc(mesh->Nelements, sizeof(dlong));

  //temp patch storage
  dfloat* patchA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* invRefAA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));

  (*Npatches) = 1;
  dlong refPatches = 0;

  //build a mini mesh struct for the reference patch
  mesh_t* refMesh = (mesh_t*) calloc(1,sizeof(mesh_t));
  memcpy(refMesh,mesh,sizeof(mesh_t));

  //vertices of reference patch
  dfloat V1x = -1., V2x =  1., V3x =  1., V4x = -1.;
  dfloat V1y = -1., V2y = -1., V3y =  1., V4y =  1.;

  refMesh->Nelements = 1;

  refMesh->EX = (dfloat*) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EY = (dfloat*) calloc(mesh->Nverts,sizeof(dfloat));

  refMesh->EX[0] = V1x;
  refMesh->EY[0] = V1y;
  refMesh->EX[1] = V2x;
  refMesh->EY[1] = V2y;
  refMesh->EX[2] = V3x;
  refMesh->EY[2] = V3y;
  refMesh->EX[3] = V4x;
  refMesh->EY[3] = V4y;

  refMesh->EToV = (hlong*) calloc(mesh->Nverts, sizeof(hlong));

  refMesh->EToV[0] = 0;
  refMesh->EToV[1] = 1;
  refMesh->EToV[2] = 2;
  refMesh->EToV[3] = 3;

  refMesh->EToB = (int*) calloc(mesh->Nfaces,sizeof(int));
  for (int n = 0; n < mesh->Nfaces; n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesQuad2D(refMesh, mesh->N);
  meshPhysicalNodesQuad2D(refMesh);
  meshGeometricFactorsQuad2D(refMesh);
  meshConnectFaceNodes2D(refMesh);
  meshSurfaceGeometricFactorsQuad2D(refMesh);

  //start with reference patch
  dfloat* refPatchInvA = *patchesInvA;
  if (options.compareArgs("DISCRETIZATION","IPDG"))
    BuildLocalIpdgPatchAxQuad2D(elliptic, refMesh, lambda, B,Br,Bs, 0, refPatchInvA);
  else if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
    BuildLocalContinuousPatchAxQuad2D(elliptic, refMesh, lambda, 0, B,Br,Bs, refPatchInvA);

  matrixInverse(mesh->Np, refPatchInvA);

  // loop over all elements
  for(dlong eM = 0; eM < mesh->Nelements; ++eM) {
    //build the patch A matrix for this element
    if (options.compareArgs("DISCRETIZATION","IPDG"))
      BuildLocalIpdgPatchAxQuad2D(elliptic, mesh, lambda, B,Br,Bs, eM, patchA);
    else if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
      BuildLocalContinuousPatchAxQuad2D(elliptic, mesh, lambda, eM, B,Br,Bs, refPatchInvA);

    dlong eP0 = mesh->EToE[eM * mesh->Nfaces + 0];
    dlong eP1 = mesh->EToE[eM * mesh->Nfaces + 1];
    dlong eP2 = mesh->EToE[eM * mesh->Nfaces + 2];
    dlong eP3 = mesh->EToE[eM * mesh->Nfaces + 3];

    int fP0 = mesh->EToF[eM * mesh->Nfaces + 0];
    int fP1 = mesh->EToF[eM * mesh->Nfaces + 1];
    int fP2 = mesh->EToF[eM * mesh->Nfaces + 2];
    int fP3 = mesh->EToF[eM * mesh->Nfaces + 3];

    if(eP0 >= 0 && eP1 >= 0 && eP2 >= 0 && eP3 >= 0) { //check if this is an interior patch
      refPatchInvA = *patchesInvA;

      //hit the patch with the reference inverse
      for(int n = 0; n < mesh->Np; ++n)
        for(int m = 0; m < mesh->Np; ++m) {
          invRefAA[n * mesh->Np + m] = 0.;
          for (int k = 0; k < mesh->Np; k++)
            invRefAA[n * mesh->Np + m] += refPatchInvA[n * mesh->Np + k] * patchA[k * mesh->Np + m];
        }

      dfloat cond = matrixConditionNumber(mesh->Np,invRefAA);
      dfloat rate = (sqrt(cond) - 1.) / (sqrt(cond) + 1.);

      // printf("Element %d's conditioned patch reports cond = %g and rate = %g \n", eM, cond, rate);

      if (rate < rateTolerance) {
        (*patchesIndex)[eM] = 0;
        refPatches++;
        continue;
      }
    }
    ++(*Npatches);
    *patchesInvA = (dfloat*) realloc(*patchesInvA,
                                     (*Npatches) * mesh->Np * mesh->Np * sizeof(dfloat));

    matrixInverse(mesh->Np, patchA);

    //copy inverse into patchesInvA
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        dlong id = ((*Npatches) - 1) * mesh->Np * mesh->Np + n * mesh->Np + m;
        (*patchesInvA)[id] = patchA[n * mesh->Np + m];
      }

    (*patchesIndex)[eM] = (*Npatches) - 1;
  }

  printf("saving " dlongFormat " full patches\n",*Npatches);
  printf("using " dlongFormat " reference patches\n", refPatches);

  free(refMesh);
  free(patchA);
  free(invRefAA);
  free(B);
  free(Br);
  free(Bs);
}

void ellipticBuildLocalPatchesTet3D(elliptic_t* elliptic, dfloat lambda, dfloat rateTolerance,
                                    dlong* Npatches, dlong** patchesIndex, dfloat** patchesInvA)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  // surface mass matrices MS = MM*LIFT
  dfloat* MS = (dfloat*) calloc(mesh->Nfaces * mesh->Nfp * mesh->Nfp,sizeof(dfloat));
  for (int f = 0; f < mesh->Nfaces; f++)
    for (int n = 0; n < mesh->Nfp; n++) {
      int fn = mesh->faceNodes[f * mesh->Nfp + n];

      for (int m = 0; m < mesh->Nfp; m++) {
        dfloat MSnm = 0;

        for (int i = 0; i < mesh->Np; i++)
          MSnm += mesh->MM[fn + i * mesh->Np] *
                  mesh->LIFT[i * mesh->Nfp * mesh->Nfaces + f * mesh->Nfp + m];

        MS[m + n * mesh->Nfp + f * mesh->Nfp * mesh->Nfp]  = MSnm;
      }
    }

  (*Npatches) = 1;
  dlong refPatches = 0;

  //build a mini mesh struct for the reference patch
  mesh_t* refMesh = (mesh_t*) calloc(1,sizeof(mesh_t));
  memcpy(refMesh,mesh,sizeof(mesh_t));

  //vertices of reference patch
  dfloat V1x = -1., V2x = 1., V3x =        0., V4x = 0;
  dfloat V1y =  0., V2y = 0., V3y =  sqrt(3.), V4y = 1. / sqrt(3.);
  dfloat V1z =  0., V2z = 0., V3z =        0., V4z = 2 * sqrt(6.) / 3.;

  refMesh->Nelements = 1;

  refMesh->EX = (dfloat*) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EY = (dfloat*) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EZ = (dfloat*) calloc(mesh->Nverts,sizeof(dfloat));

  refMesh->EX[0] = V1x;
  refMesh->EY[0] = V1y;
  refMesh->EZ[0] = V1z;
  refMesh->EX[1] = V2x;
  refMesh->EY[1] = V2y;
  refMesh->EZ[1] = V2z;
  refMesh->EX[2] = V3x;
  refMesh->EY[2] = V3y;
  refMesh->EZ[2] = V3z;
  refMesh->EX[3] = V4x;
  refMesh->EY[3] = V4y;
  refMesh->EZ[3] = V4z;

  refMesh->EToV = (hlong*) calloc(mesh->Nverts, sizeof(hlong));

  refMesh->EToV[0] = 0;
  refMesh->EToV[1] = 1;
  refMesh->EToV[2] = 2;
  refMesh->EToV[3] = 3;

  refMesh->EToB = (int*) calloc(mesh->Nfaces,sizeof(int));
  for (int n = 0; n < mesh->Nfaces; n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesTet3D(refMesh, mesh->N);
  meshPhysicalNodesTet3D(refMesh);
  meshGeometricFactorsTet3D(refMesh);
  meshConnectFaceNodes3D(refMesh);
  meshSurfaceGeometricFactorsTet3D(refMesh);

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  *patchesIndex = (dlong*) calloc(mesh->Nelements, sizeof(dlong));

  //temp patch storage
  dfloat* patchA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* invRefAA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));

  //start with reference patch
  dfloat* refPatchInvA = *patchesInvA;
  if (options.compareArgs("DISCRETIZATION","IPDG"))
    BuildLocalIpdgPatchAxTet3D(elliptic, refMesh, lambda, MS, 0, refPatchInvA);
  else if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
    BuildLocalContinuousPatchAxTet3D(elliptic, refMesh, lambda, 0, refPatchInvA);

  matrixInverse(mesh->Np, refPatchInvA);

  dfloat maxRate = 0.;
  dfloat maxCond = 0.;

  // loop over all elements
  for(dlong eM = 0; eM < mesh->Nelements; ++eM) {
    //build the patch A matrix for this element
    if (options.compareArgs("DISCRETIZATION","IPDG"))
      BuildLocalIpdgPatchAxTet3D(elliptic, mesh, lambda, MS, eM, refPatchInvA);
    else if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
      BuildLocalContinuousPatchAxTet3D(elliptic, refMesh, lambda, eM, refPatchInvA);

    dlong eP0 = mesh->EToE[eM * mesh->Nfaces + 0];
    dlong eP1 = mesh->EToE[eM * mesh->Nfaces + 1];
    dlong eP2 = mesh->EToE[eM * mesh->Nfaces + 2];
    dlong eP3 = mesh->EToE[eM * mesh->Nfaces + 3];

    int fP0 = mesh->EToF[eM * mesh->Nfaces + 0];
    int fP1 = mesh->EToF[eM * mesh->Nfaces + 1];
    int fP2 = mesh->EToF[eM * mesh->Nfaces + 2];
    int fP3 = mesh->EToF[eM * mesh->Nfaces + 3];

    if(eP0 >= 0 && eP1 >= 0 && eP2 >= 0 && eP3 >= 0) { //check if this is an interior patch
      refPatchInvA = *patchesInvA;

      //hit the patch with the reference inverse
      for(int n = 0; n < mesh->Np; ++n)
        for(int m = 0; m < mesh->Np; ++m) {
          invRefAA[n * mesh->Np + m] = 0.;
          for (int k = 0; k < mesh->Np; k++)
            invRefAA[n * mesh->Np + m] += refPatchInvA[n * mesh->Np + k] * patchA[k * mesh->Np + m];
        }

      dfloat cond = matrixConditionNumber(mesh->Np,invRefAA);
      dfloat rate = (sqrt(cond) - 1.) / (sqrt(cond) + 1.);

      //printf("Element %d's conditioned patch reports cond = %g and rate = %g \n", eM, cond, rate);
      maxRate = mymax(rate,maxRate);
      maxCond = mymax(cond,maxCond);

      if (rate < rateTolerance) {
        (*patchesIndex)[eM] = 0;
        refPatches++;
        continue;
      }
    }
    ++(*Npatches);
    *patchesInvA = (dfloat*) realloc(*patchesInvA,
                                     (*Npatches) * mesh->Np * mesh->Np * sizeof(dfloat));

    matrixInverse(mesh->Np, patchA);

    //copy inverse into patchesInvA
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        int id = ((*Npatches) - 1) * mesh->Np * mesh->Np + n * mesh->Np + m;
        (*patchesInvA)[id] = patchA[n * mesh->Np + m];
      }

    (*patchesIndex)[eM] = (*Npatches) - 1;
  }

  printf("saving " dlongFormat " full patches\n",*Npatches);
  printf("using " dlongFormat " reference patches\n", refPatches);
  printf("Max condition number = %g, and slowest CG convergence rate = %g\n", maxCond, maxRate);

  free(refMesh);
  free(patchA);
  free(invRefAA);
  free(MS);
}

void ellipticBuildLocalPatchesHex3D(elliptic_t* elliptic, dfloat lambda, dfloat rateTolerance,
                                    dlong* Npatches, dlong** patchesIndex, dfloat** patchesInvA)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  // build some monolithic basis arrays
  dfloat* B  = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* Br = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* Bs = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* Bt = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nk = 0; nk < mesh->N + 1; ++nk)
    for(int nj = 0; nj < mesh->N + 1; ++nj)
      for(int ni = 0; ni < mesh->N + 1; ++ni) {
        int node = 0;

        for(int k = 0; k < mesh->N + 1; ++k)
          for(int j = 0; j < mesh->N + 1; ++j)
            for(int i = 0; i < mesh->N + 1; ++i) {
              if(nk == k && nj == j && ni == i)
                B[mode * mesh->Np + node] = 1;
              if(nj == j && nk == k)
                Br[mode * mesh->Np + node] = mesh->D[ni + mesh->Nq * i];
              if(ni == i && nk == k)
                Bs[mode * mesh->Np + node] = mesh->D[nj + mesh->Nq * j];
              if(ni == i && nj == j)
                Bt[mode * mesh->Np + node] = mesh->D[nk + mesh->Nq * k];

              ++node;
            }

        ++mode;
      }

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  *patchesIndex = (dlong*) calloc(mesh->Nelements, sizeof(dlong));

  //temp patch storage
  dfloat* patchA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* invRefAA = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));

  (*Npatches) = 1;
  dlong refPatches = 0;

  //build a mini mesh struct for the reference patch
  mesh_t* refMesh = (mesh_t*) calloc(1,sizeof(mesh_t));
  memcpy(refMesh,mesh,sizeof(mesh_t));

  //vertices of reference patch
  dfloat V1x = -1., V2x =  1., V3x =  1., V4x = -1., V5x = -1., V6x =  1., V7x =  1., V8x = -1.;
  dfloat V1y = -1., V2y = -1., V3y =  1., V4y =  1., V5y = -1., V6y = -1., V7y =  1., V8y =  1.;
  dfloat V1z = -1., V2z = -1., V3z = -1., V4z = -1., V5z =  1., V6z =  1., V7z =  1., V8z =  1.;

  refMesh->Nelements = 1;

  refMesh->EX = (dfloat*) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EY = (dfloat*) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EZ = (dfloat*) calloc(mesh->Nverts,sizeof(dfloat));

  refMesh->EX[0] = V1x;
  refMesh->EY[0] = V1y;
  refMesh->EZ[0] = V1z;
  refMesh->EX[1] = V2x;
  refMesh->EY[1] = V2y;
  refMesh->EZ[1] = V2z;
  refMesh->EX[2] = V3x;
  refMesh->EY[2] = V3y;
  refMesh->EZ[2] = V3z;
  refMesh->EX[3] = V4x;
  refMesh->EY[3] = V4y;
  refMesh->EZ[3] = V4z;
  refMesh->EX[4] = V5x;
  refMesh->EY[4] = V5y;
  refMesh->EZ[4] = V5z;
  refMesh->EX[5] = V6x;
  refMesh->EY[5] = V6y;
  refMesh->EZ[5] = V6z;
  refMesh->EX[6] = V7x;
  refMesh->EY[6] = V7y;
  refMesh->EZ[6] = V7z;
  refMesh->EX[7] = V8x;
  refMesh->EY[7] = V8y;
  refMesh->EZ[7] = V8z;

  refMesh->EToV = (hlong*) calloc(mesh->Nverts, sizeof(hlong));

  refMesh->EToV[0] = 0;
  refMesh->EToV[1] = 1;
  refMesh->EToV[2] = 2;
  refMesh->EToV[3] = 3;
  refMesh->EToV[4] = 4;
  refMesh->EToV[5] = 5;
  refMesh->EToV[6] = 6;
  refMesh->EToV[7] = 7;

  refMesh->EToB = (int*) calloc(mesh->Nfaces,sizeof(int));
  for (int n = 0; n < mesh->Nfaces; n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesHex3D(refMesh, mesh->N);
  meshPhysicalNodesHex3D(refMesh);
  meshGeometricFactorsHex3D(refMesh);
  meshConnectFaceNodes3D(refMesh);
  meshSurfaceGeometricFactorsHex3D(refMesh);

  //start with reference patch
  dfloat* refPatchInvA = *patchesInvA;
  if (options.compareArgs("DISCRETIZATION","IPDG"))
    BuildLocalIpdgPatchAxHex3D(elliptic, refMesh, lambda, B,Br,Bs,Bt, 0, refPatchInvA);
  else if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
    BuildLocalContinuousPatchAxHex3D(elliptic, refMesh, lambda,  0, B,Br,Bs,Bt, refPatchInvA);
  matrixInverse(mesh->Np, refPatchInvA);

  // loop over all elements
  for(dlong eM = 0; eM < mesh->Nelements; ++eM) {
    //build the patch A matrix for this element
    if (options.compareArgs("DISCRETIZATION","IPDG"))
      BuildLocalIpdgPatchAxHex3D(elliptic, mesh, lambda, B,Br,Bs,Bt, eM,  refPatchInvA);
    else if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
      BuildLocalContinuousPatchAxHex3D(elliptic, mesh, lambda,  eM, B,Br,Bs,Bt, refPatchInvA);

    dlong eP0 = mesh->EToE[eM * mesh->Nfaces + 0];
    dlong eP1 = mesh->EToE[eM * mesh->Nfaces + 1];
    dlong eP2 = mesh->EToE[eM * mesh->Nfaces + 2];
    dlong eP3 = mesh->EToE[eM * mesh->Nfaces + 3];
    dlong eP4 = mesh->EToE[eM * mesh->Nfaces + 4];
    dlong eP5 = mesh->EToE[eM * mesh->Nfaces + 5];

    int fP0 = mesh->EToF[eM * mesh->Nfaces + 0];
    int fP1 = mesh->EToF[eM * mesh->Nfaces + 1];
    int fP2 = mesh->EToF[eM * mesh->Nfaces + 2];
    int fP3 = mesh->EToF[eM * mesh->Nfaces + 3];
    int fP4 = mesh->EToF[eM * mesh->Nfaces + 4];
    int fP5 = mesh->EToF[eM * mesh->Nfaces + 5];

    if(eP0 >= 0 && eP1 >= 0 && eP2 >= 0 && eP3 >= 0 && eP4 >= 0 && eP5 >= 0) { //check if this is an interior patch
      refPatchInvA = *patchesInvA;

      //hit the patch with the reference inverse
      for(int n = 0; n < mesh->Np; ++n)
        for(int m = 0; m < mesh->Np; ++m) {
          invRefAA[n * mesh->Np + m] = 0.;
          for (int k = 0; k < mesh->Np; k++)
            invRefAA[n * mesh->Np + m] += refPatchInvA[n * mesh->Np + k] * patchA[k * mesh->Np + m];
        }

      dfloat cond = matrixConditionNumber(mesh->Np,invRefAA);
      dfloat rate = (sqrt(cond) - 1.) / (sqrt(cond) + 1.);

      // printf("Element %d's conditioned patch reports cond = %g and rate = %g \n", eM, cond, rate);

      if (rate < rateTolerance) {
        (*patchesIndex)[eM] = 0;
        refPatches++;
        continue;
      }
    }
    ++(*Npatches);
    *patchesInvA = (dfloat*) realloc(*patchesInvA,
                                     (*Npatches) * mesh->Np * mesh->Np * sizeof(dfloat));

    matrixInverse(mesh->Np, patchA);

    //copy inverse into patchesInvA
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        dlong id = ((*Npatches) - 1) * mesh->Np * mesh->Np + n * mesh->Np + m;
        (*patchesInvA)[id] = patchA[n * mesh->Np + m];
      }

    (*patchesIndex)[eM] = (*Npatches) - 1;
  }

  printf("saving " dlongFormat " full patches\n",*Npatches);
  printf("using " dlongFormat " reference patches\n", refPatches);

  free(refMesh);
  free(patchA);
  free(invRefAA);
  free(B);
  free(Br);
  free(Bs);
}

//returns the ipdg patch A matrix for element eM
void BuildLocalIpdgPatchAxTri2D(elliptic_t* elliptic,
                                mesh_t* mesh,
                                int basisNp,
                                dfloat* basis,
                                dfloat lambda,
                                dfloat* MS,
                                dlong eM,
                                dfloat* A)
{
  dlong vbase = eM * mesh->Nvgeo;
  dfloat drdx = mesh->vgeo[vbase + RXID];
  dfloat drdy = mesh->vgeo[vbase + RYID];
  dfloat dsdx = mesh->vgeo[vbase + SXID];
  dfloat dsdy = mesh->vgeo[vbase + SYID];
  dfloat J = mesh->vgeo[vbase + JID];

  dfloat* Ae = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));

  /* start with stiffness matrix  */
  for(int n = 0; n < mesh->Np; ++n)
    for(int m = 0; m < mesh->Np; ++m) {
      Ae[n * mesh->Np + m]  = J * lambda * mesh->MM[n * mesh->Np + m];
      Ae[n * mesh->Np + m] += J * drdx * drdx * mesh->Srr[n * mesh->Np + m];
      Ae[n * mesh->Np + m] += J * drdx * dsdx * mesh->Srs[n * mesh->Np + m];
      Ae[n * mesh->Np + m] += J * dsdx * drdx * mesh->Ssr[n * mesh->Np + m];
      Ae[n * mesh->Np + m] += J * dsdx * dsdx * mesh->Sss[n * mesh->Np + m];

      Ae[n * mesh->Np + m] += J * drdy * drdy * mesh->Srr[n * mesh->Np + m];
      Ae[n * mesh->Np + m] += J * drdy * dsdy * mesh->Srs[n * mesh->Np + m];
      Ae[n * mesh->Np + m] += J * dsdy * drdy * mesh->Ssr[n * mesh->Np + m];
      Ae[n * mesh->Np + m] += J * dsdy * dsdy * mesh->Sss[n * mesh->Np + m];
    }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m)
        Ae[n * mesh->Np + m] += elliptic->allNeumannPenalty * elliptic->allNeumannScale *
                                elliptic->allNeumannScale;
  }

  for (int fM = 0; fM < mesh->Nfaces; fM++) {
    // load surface geofactors for this face
    dlong sid = mesh->Nsgeo * (eM * mesh->Nfaces + fM);
    dfloat nx = mesh->sgeo[sid + NXID];
    dfloat ny = mesh->sgeo[sid + NYID];
    dfloat sJ = mesh->sgeo[sid + SJID];
    dfloat hinv = mesh->sgeo[sid + IHID];

    int bc = mesh->EToB[fM + mesh->Nfaces * eM]; //raw boundary flag

    dfloat penalty = elliptic->tau * hinv;

    int bcD = 0, bcN = 0;
    int bcType = 0;

    if(bc > 0) bcType = elliptic->BCType[bc];        //find its type (Dirichlet/Neumann)

    // this needs to be double checked (and the code where these are used)
    if(bcType == 1) { // Dirichlet
      bcD = 1;
      bcN = 0;
    } else if(bcType == 2) { // Neumann
      bcD = 0;
      bcN = 1;
    }

    // mass matrix for this face
    dfloat* MSf = MS + fM * mesh->Nfp * mesh->Nfp;

    // penalty term just involves face nodes
    for(int n = 0; n < mesh->Nfp; ++n)
      for(int m = 0; m < mesh->Nfp; ++m) {
        int nM = mesh->faceNodes[fM * mesh->Nfp + n];
        int mM = mesh->faceNodes[fM * mesh->Nfp + m];

        // OP11 = OP11 + 0.5*( gtau*mmE )
        dfloat MSfnm = sJ * MSf[n * mesh->Nfp + m];
        Ae[nM * mesh->Np + mM] += 0.5 * (1. - bcN) * (1. + bcD) * penalty * MSfnm;
      }

    // now add differential surface terms
    for(int n = 0; n < mesh->Nfp; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        int nM = mesh->faceNodes[fM * mesh->Nfp + n];

        for(int i = 0; i < mesh->Nfp; ++i) {
          int iM = mesh->faceNodes[fM * mesh->Nfp + i];

          dfloat MSfni = sJ * MSf[n * mesh->Nfp + i]; // surface Jacobian built in

          dfloat DxMim = drdx * mesh->Dr[iM * mesh->Np + m] + dsdx * mesh->Ds[iM * mesh->Np + m];
          dfloat DyMim = drdy * mesh->Dr[iM * mesh->Np + m] + dsdy * mesh->Ds[iM * mesh->Np + m];

          // OP11 = OP11 + 0.5*( - mmE*Dn1)
          Ae[nM * mesh->Np + m] += -0.5 * nx * (1 + bcD) * (1 - bcN) * MSfni * DxMim;
          Ae[nM * mesh->Np + m] += -0.5 * ny * (1 + bcD) * (1 - bcN) * MSfni * DyMim;
        }
      }

    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Nfp; ++m) {
        int mM = mesh->faceNodes[fM * mesh->Nfp + m];

        for(int i = 0; i < mesh->Nfp; ++i) {
          int iM = mesh->faceNodes[fM * mesh->Nfp + i];

          dfloat MSfim = sJ * MSf[i * mesh->Nfp + m];

          dfloat DxMin = drdx * mesh->Dr[iM * mesh->Np + n] + dsdx * mesh->Ds[iM * mesh->Np + n];
          dfloat DyMin = drdy * mesh->Dr[iM * mesh->Np + n] + dsdy * mesh->Ds[iM * mesh->Np + n];

          // OP11 = OP11 + (- Dn1'*mmE );
          Ae[n * mesh->Np + mM] +=  -0.5 * nx * (1 + bcD) * (1 - bcN) * DxMin * MSfim;
          Ae[n * mesh->Np + mM] +=  -0.5 * ny * (1 + bcD) * (1 - bcN) * DyMin * MSfim;
        }
      }
  }

  for(int j = 0; j < basisNp; ++j)
    for(int i = 0; i < basisNp; ++i) {
      dfloat val = 0;
      for (int n = 0; n < mesh->Np; n++)
        for (int m = 0; m < mesh->Np; m++)
          val += basis[n * basisNp + j] * Ae[n * mesh->Np + m] * basis[m * basisNp + i];

      A[i + j * basisNp] = val;
    }

  free(Ae);
}

//returns the continuous C0 patch A matrix for element eM
void BuildLocalContinuousPatchAxTri2D(elliptic_t* elliptic,
                                      mesh_t* mesh,
                                      dfloat lambda,
                                      dlong eM,
                                      dfloat* A)
{
  dlong gbase = eM * mesh->Nggeo;
  dfloat Grr = mesh->ggeo[gbase + G00ID];
  dfloat Grs = mesh->ggeo[gbase + G01ID];
  dfloat Gss = mesh->ggeo[gbase + G11ID];
  dfloat J   = mesh->ggeo[gbase + GWJID];

  /* start with stiffness matrix  */
  for(int n = 0; n < mesh->Np; ++n) {
    if (elliptic->mapB[n + eM * mesh->Np] != 1) { //dont fill rows for masked nodes
      for(int m = 0; m < mesh->Np; ++m) {
        if (elliptic->mapB[m + eM * mesh->Np] != 1) {//dont fill rows for masked nodes
          A[n * mesh->Np + m] = J * lambda * mesh->MM[m + n * mesh->Np];
          A[n * mesh->Np + m] += Grr * mesh->Srr[m + n * mesh->Np];
          A[n * mesh->Np + m] += Grs * mesh->Srs[m + n * mesh->Np];
          A[n * mesh->Np + m] += Grs * mesh->Ssr[m + n * mesh->Np];
          A[n * mesh->Np + m] += Gss * mesh->Sss[m + n * mesh->Np];
        } else {
          A[n * mesh->Np + m] = 0;
        }
      }
    } else {
      A[n + n * mesh->Np] = 1; //just put a 1 so A is invertable
    }
  }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n = 0; n < mesh->Np; ++n)
      if (elliptic->mapB[n + eM * mesh->Np] != 1) { //dont fill rows for masked nodes
        for(int m = 0; m < mesh->Np; ++m) {
          if (elliptic->mapB[m + eM * mesh->Np] == 1) continue; //skip masked nodes
          A[n * mesh->Np + m] += elliptic->allNeumannPenalty * elliptic->allNeumannScale *
                                 elliptic->allNeumannScale;
        }
      }
  }
}

//returns the patch A matrix for element eM
void BuildLocalIpdgPatchAxQuad2D(elliptic_t* elliptic, mesh_t* mesh, dfloat lambda,
                                 dfloat* B, dfloat* Br, dfloat* Bs, dlong eM, dfloat* A)
{
  /* start with stiffness matrix  */
  for(int n = 0; n < mesh->Np; ++n)
    for(int m = 0; m < mesh->Np; ++m) {
      A[n * mesh->Np + m] = 0;

      // (grad phi_n, grad phi_m)_{D^e}
      for(int i = 0; i < mesh->Np; ++i) {
        dlong base = eM * mesh->Np * mesh->Nvgeo + i;
        dfloat drdx = mesh->vgeo[base + mesh->Np * RXID];
        dfloat drdy = mesh->vgeo[base + mesh->Np * RYID];
        dfloat dsdx = mesh->vgeo[base + mesh->Np * SXID];
        dfloat dsdy = mesh->vgeo[base + mesh->Np * SYID];
        dfloat JW   = mesh->vgeo[base + mesh->Np * JWID];

        int idn = n * mesh->Np + i;
        int idm = m * mesh->Np + i;
        dfloat dlndx = drdx * Br[idn] + dsdx * Bs[idn];
        dfloat dlndy = drdy * Br[idn] + dsdy * Bs[idn];
        dfloat dlmdx = drdx * Br[idm] + dsdx * Bs[idm];
        dfloat dlmdy = drdy * Br[idm] + dsdy * Bs[idm];
        A[n * mesh->Np + m] += JW * (dlndx * dlmdx + dlndy * dlmdy);
        A[n * mesh->Np + m] += lambda * JW * B[idn] * B[idm];
      }

      for (int fM = 0; fM < mesh->Nfaces; fM++)
        // accumulate flux terms for negative and positive traces
        for(int i = 0; i < mesh->Nfp; ++i) {
          int vidM = mesh->faceNodes[i + fM * mesh->Nfp];

          // grab vol geofacs at surface nodes
          dlong baseM = eM * mesh->Np * mesh->Nvgeo + vidM;
          dfloat drdxM = mesh->vgeo[baseM + mesh->Np * RXID];
          dfloat drdyM = mesh->vgeo[baseM + mesh->Np * RYID];
          dfloat dsdxM = mesh->vgeo[baseM + mesh->Np * SXID];
          dfloat dsdyM = mesh->vgeo[baseM + mesh->Np * SYID];

          // grab surface geometric factors
          dlong base = mesh->Nsgeo * (eM * mesh->Nfp * mesh->Nfaces + fM * mesh->Nfp + i);
          dfloat nx = mesh->sgeo[base + NXID];
          dfloat ny = mesh->sgeo[base + NYID];
          dfloat wsJ = mesh->sgeo[base + WSJID];
          dfloat hinv = mesh->sgeo[base + IHID];

          // form negative trace terms in IPDG
          int idnM = n * mesh->Np + vidM;
          int idmM = m * mesh->Np + vidM;

          dfloat dlndxM = drdxM * Br[idnM] + dsdxM * Bs[idnM];
          dfloat dlndyM = drdyM * Br[idnM] + dsdyM * Bs[idnM];
          dfloat ndotgradlnM = nx * dlndxM + ny * dlndyM;
          dfloat lnM = B[idnM];

          dfloat dlmdxM = drdxM * Br[idmM] + dsdxM * Bs[idmM];
          dfloat dlmdyM = drdyM * Br[idmM] + dsdyM * Bs[idmM];
          dfloat ndotgradlmM = nx * dlmdxM + ny * dlmdyM;
          dfloat lmM = B[idmM];

          dfloat penalty = elliptic->tau * hinv;
          int bc = mesh->EToB[fM + mesh->Nfaces * eM]; //raw boundary flag

          int bcD = 0, bcN = 0;
          int bcType = 0;

          if(bc > 0) bcType = elliptic->BCType[bc];        //find its type (Dirichlet/Neumann)

          // this needs to be double checked (and the code where these are used)
          if(bcType == 1) { // Dirichlet
            bcD = 1;
            bcN = 0;
          } else if(bcType == 2) { // Neumann
            bcD = 0;
            bcN = 1;
          }

          A[n * mesh->Np + m] += -0.5 * (1 + bcD) * (1 - bcN) * wsJ * lnM * ndotgradlmM;  // -(ln^-, N.grad lm^-)
          A[n * mesh->Np + m] += -0.5 * (1 + bcD) * (1 - bcN) * wsJ * ndotgradlnM * lmM;  // -(N.grad ln^-, lm^-)
          A[n * mesh->Np + m] += +0.5 * (1 + bcD) * (1 - bcN) * wsJ * penalty * lnM * lmM; // +((tau/h)*ln^-,lm^-)
        }
    }
}

void BuildLocalContinuousPatchAxQuad2D(elliptic_t* elliptic, mesh_t* mesh, dfloat lambda,
                                       dlong eM, dfloat* B, dfloat* Br, dfloat* Bs, dfloat* A)
{
  for (int ny = 0; ny < mesh->Nq; ny++)
    for (int nx = 0; nx < mesh->Nq; nx++) {
      if (elliptic->mapB[nx + ny * mesh->Nq + eM * mesh->Np] != 1) {
        for (int my = 0; my < mesh->Nq; my++)
          for (int mx = 0; mx < mesh->Nq; mx++) {
            if (elliptic->mapB[mx + my * mesh->Nq + eM * mesh->Np] == 1) continue;

            int id;
            int iid = (nx + ny * mesh->Nq) * mesh->Np + mx + my * mesh->Nq;
            A[iid] = 0;

            if (ny == my) {
              for (int k = 0; k < mesh->Nq; k++) {
                id = k + ny * mesh->Nq;
                dfloat Grr = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G00ID * mesh->Np];

                A[iid] += Grr * mesh->D[nx + k * mesh->Nq] * mesh->D[mx + k * mesh->Nq];
              }
            }

            id = mx + ny * mesh->Nq;
            dfloat Grs = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G01ID * mesh->Np];
            A[iid] += Grs * mesh->D[nx + mx * mesh->Nq] * mesh->D[my + ny * mesh->Nq];

            id = nx + my * mesh->Nq;
            dfloat Gsr = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G01ID * mesh->Np];
            A[iid] += Gsr * mesh->D[mx + nx * mesh->Nq] * mesh->D[ny + my * mesh->Nq];

            if (nx == mx) {
              for (int k = 0; k < mesh->Nq; k++) {
                id = nx + k * mesh->Nq;
                dfloat Gss = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G11ID * mesh->Np];

                A[iid] += Gss * mesh->D[ny + k * mesh->Nq] * mesh->D[my + k * mesh->Nq];
              }
            }

            if ((nx == mx) && (ny == my)) {
              id = nx + ny * mesh->Nq;
              dfloat JW = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + GWJID * mesh->Np];
              A[iid] += JW * lambda;
            }
          }
      } else {
        int iid = (nx + ny * mesh->Nq) * mesh->Np + nx + ny * mesh->Nq;
        A[iid] = 1; //just put a 1 so A is invertable
      }
    }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n = 0; n < mesh->Np; ++n)
      if (elliptic->mapB[n + eM * mesh->Np] != 1) { //dont fill rows for masked nodes
        for(int m = 0; m < mesh->Np; ++m) {
          if (elliptic->mapB[m + eM * mesh->Np] == 1) continue; //skip masked nodes
          A[n * mesh->Np + m] += elliptic->allNeumannPenalty * elliptic->allNeumannScale *
                                 elliptic->allNeumannScale;
        }
      }
  }
}

//returns the patch A matrix for element eM
void BuildLocalIpdgPatchAxTet3D(elliptic_t* elliptic,
                                mesh_t* mesh,
                                dfloat lambda,
                                dfloat* MS,
                                dlong eM,
                                dfloat* A)
{
  dlong vbase = eM * mesh->Nvgeo;
  dfloat drdx = mesh->vgeo[vbase + RXID];
  dfloat drdy = mesh->vgeo[vbase + RYID];
  dfloat drdz = mesh->vgeo[vbase + RZID];
  dfloat dsdx = mesh->vgeo[vbase + SXID];
  dfloat dsdy = mesh->vgeo[vbase + SYID];
  dfloat dsdz = mesh->vgeo[vbase + SZID];
  dfloat dtdx = mesh->vgeo[vbase + TXID];
  dfloat dtdy = mesh->vgeo[vbase + TYID];
  dfloat dtdz = mesh->vgeo[vbase + TZID];
  dfloat J = mesh->vgeo[vbase + JID];

  dfloat G00 = drdx * drdx + drdy * drdy + drdz * drdz;
  dfloat G01 = drdx * dsdx + drdy * dsdy + drdz * dsdz;
  dfloat G02 = drdx * dtdx + drdy * dtdy + drdz * dtdz;

  dfloat G10 = dsdx * drdx + dsdy * drdy + dsdz * drdz;
  dfloat G11 = dsdx * dsdx + dsdy * dsdy + dsdz * dsdz;
  dfloat G12 = dsdx * dtdx + dsdy * dtdy + dsdz * dtdz;

  dfloat G20 = dtdx * drdx + dtdy * drdy + dtdz * drdz;
  dfloat G21 = dtdx * dsdx + dtdy * dsdy + dtdz * dsdz;
  dfloat G22 = dtdx * dtdx + dtdy * dtdy + dtdz * dtdz;

  /* start with stiffness matrix  */
  for(int n = 0; n < mesh->Np; ++n)
    for(int m = 0; m < mesh->Np; ++m) {
      A[n * mesh->Np + m]  = J * lambda * mesh->MM[n * mesh->Np + m];
      A[n * mesh->Np + m] += J * G00 * mesh->Srr[n * mesh->Np + m];
      A[n * mesh->Np + m] += J * G01 * mesh->Srs[n * mesh->Np + m];
      A[n * mesh->Np + m] += J * G02 * mesh->Srt[n * mesh->Np + m];
      A[n * mesh->Np + m] += J * G10 * mesh->Ssr[n * mesh->Np + m];
      A[n * mesh->Np + m] += J * G11 * mesh->Sss[n * mesh->Np + m];
      A[n * mesh->Np + m] += J * G12 * mesh->Sst[n * mesh->Np + m];
      A[n * mesh->Np + m] += J * G20 * mesh->Str[n * mesh->Np + m];
      A[n * mesh->Np + m] += J * G21 * mesh->Sts[n * mesh->Np + m];
      A[n * mesh->Np + m] += J * G22 * mesh->Stt[n * mesh->Np + m];
    }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m)
        A[n * mesh->Np + m] += elliptic->allNeumannPenalty * elliptic->allNeumannScale *
                               elliptic->allNeumannScale;
  }

  for (int fM = 0; fM < mesh->Nfaces; fM++) {
    // load surface geofactors for this face
    dlong sid = mesh->Nsgeo * (eM * mesh->Nfaces + fM);
    dfloat nx = mesh->sgeo[sid + NXID];
    dfloat ny = mesh->sgeo[sid + NYID];
    dfloat nz = mesh->sgeo[sid + NZID];
    dfloat sJ = mesh->sgeo[sid + SJID];
    dfloat hinv = mesh->sgeo[sid + IHID];

    int bc = mesh->EToB[fM + mesh->Nfaces * eM]; //raw boundary flag

    dfloat penalty = elliptic->tau * hinv;

    int bcD = 0, bcN = 0;
    int bcType = 0;

    if(bc > 0) bcType = elliptic->BCType[bc];        //find its type (Dirichlet/Neumann)

    // this needs to be double checked (and the code where these are used)
    if(bcType == 1) { // Dirichlet
      bcD = 1;
      bcN = 0;
    } else if(bcType == 2) { // Neumann
      bcD = 0;
      bcN = 1;
    }

    // mass matrix for this face
    dfloat* MSf = MS + fM * mesh->Nfp * mesh->Nfp;

    // penalty term just involves face nodes
    for(int n = 0; n < mesh->Nfp; ++n)
      for(int m = 0; m < mesh->Nfp; ++m) {
        int nM = mesh->faceNodes[fM * mesh->Nfp + n];
        int mM = mesh->faceNodes[fM * mesh->Nfp + m];

        // OP11 = OP11 + 0.5*( gtau*mmE )
        dfloat MSfnm = sJ * MSf[n * mesh->Nfp + m];
        A[nM * mesh->Np + mM] += 0.5 * (1. - bcN) * (1. + bcD) * penalty * MSfnm;
      }

    // now add differential surface terms
    for(int n = 0; n < mesh->Nfp; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        int nM = mesh->faceNodes[fM * mesh->Nfp + n];

        for(int i = 0; i < mesh->Nfp; ++i) {
          int iM = mesh->faceNodes[fM * mesh->Nfp + i];

          dfloat MSfni = sJ * MSf[n * mesh->Nfp + i]; // surface Jacobian built in

          dfloat DxMim = drdx * mesh->Dr[iM * mesh->Np + m] + dsdx * mesh->Ds[iM * mesh->Np + m] +
                         dtdx * mesh->Dt[iM * mesh->Np + m];
          dfloat DyMim = drdy * mesh->Dr[iM * mesh->Np + m] + dsdy * mesh->Ds[iM * mesh->Np + m] +
                         dtdy * mesh->Dt[iM * mesh->Np + m];
          dfloat DzMim = drdz * mesh->Dr[iM * mesh->Np + m] + dsdz * mesh->Ds[iM * mesh->Np + m] +
                         dtdz * mesh->Dt[iM * mesh->Np + m];

          // OP11 = OP11 + 0.5*( - mmE*Dn1)
          A[nM * mesh->Np + m] += -0.5 * nx * (1 + bcD) * (1 - bcN) * MSfni * DxMim;
          A[nM * mesh->Np + m] += -0.5 * ny * (1 + bcD) * (1 - bcN) * MSfni * DyMim;
          A[nM * mesh->Np + m] += -0.5 * nz * (1 + bcD) * (1 - bcN) * MSfni * DzMim;
        }
      }

    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Nfp; ++m) {
        int mM = mesh->faceNodes[fM * mesh->Nfp + m];

        for(int i = 0; i < mesh->Nfp; ++i) {
          int iM = mesh->faceNodes[fM * mesh->Nfp + i];

          dfloat MSfim = sJ * MSf[i * mesh->Nfp + m];

          dfloat DxMin = drdx * mesh->Dr[iM * mesh->Np + n] + dsdx * mesh->Ds[iM * mesh->Np + n] +
                         dtdx * mesh->Dt[iM * mesh->Np + n];
          dfloat DyMin = drdy * mesh->Dr[iM * mesh->Np + n] + dsdy * mesh->Ds[iM * mesh->Np + n] +
                         dtdy * mesh->Dt[iM * mesh->Np + n];
          dfloat DzMin = drdz * mesh->Dr[iM * mesh->Np + n] + dsdz * mesh->Ds[iM * mesh->Np + n] +
                         dtdz * mesh->Dt[iM * mesh->Np + n];

          // OP11 = OP11 + (- Dn1'*mmE );
          A[n * mesh->Np + mM] +=  -0.5 * nx * (1 + bcD) * (1 - bcN) * DxMin * MSfim;
          A[n * mesh->Np + mM] +=  -0.5 * ny * (1 + bcD) * (1 - bcN) * DyMin * MSfim;
          A[n * mesh->Np + mM] +=  -0.5 * nz * (1 + bcD) * (1 - bcN) * DzMin * MSfim;
        }
      }
  }
}

//returns the continuous C0 patch A matrix for element eM
void BuildLocalContinuousPatchAxTet3D(elliptic_t* elliptic,
                                      mesh_t* mesh,
                                      dfloat lambda,
                                      dlong eM,
                                      dfloat* A)
{
  dlong gbase = eM * mesh->Nggeo;
  dfloat Grr = mesh->ggeo[gbase + G00ID];
  dfloat Grs = mesh->ggeo[gbase + G01ID];
  dfloat Grt = mesh->ggeo[gbase + G02ID];
  dfloat Gss = mesh->ggeo[gbase + G11ID];
  dfloat Gst = mesh->ggeo[gbase + G12ID];
  dfloat Gtt = mesh->ggeo[gbase + G22ID];
  dfloat J   = mesh->ggeo[gbase + GWJID];

  /* start with stiffness matrix  */
  for(int n = 0; n < mesh->Np; ++n) {
    if (elliptic->mapB[n + eM * mesh->Np] != 1) { //dont fill rows for masked nodes
      for(int m = 0; m < mesh->Np; ++m) {
        if (elliptic->mapB[m + eM * mesh->Np] != 1) {//dont fill rows for masked nodes
          A[n * mesh->Np + m] = J * lambda * mesh->MM[m + n * mesh->Np];
          A[n * mesh->Np + m] += Grr * mesh->Srr[m + n * mesh->Np];
          A[n * mesh->Np + m] += Grs * mesh->Srs[m + n * mesh->Np];
          A[n * mesh->Np + m] += Grt * mesh->Srt[m + n * mesh->Np];
          A[n * mesh->Np + m] += Grs * mesh->Ssr[m + n * mesh->Np];
          A[n * mesh->Np + m] += Gss * mesh->Sss[m + n * mesh->Np];
          A[n * mesh->Np + m] += Gst * mesh->Sst[m + n * mesh->Np];
          A[n * mesh->Np + m] += Grt * mesh->Str[m + n * mesh->Np];
          A[n * mesh->Np + m] += Gst * mesh->Sts[m + n * mesh->Np];
          A[n * mesh->Np + m] += Gtt * mesh->Stt[m + n * mesh->Np];
        } else {
          A[n * mesh->Np + m] = 0;
        }
      }
    } else {
      A[n + n * mesh->Np] = 1; //just put a 1 so A is invertable
    }
  }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n = 0; n < mesh->Np; ++n)
      if (elliptic->mapB[n + eM * mesh->Np] != 1) { //dont fill rows for masked nodes
        for(int m = 0; m < mesh->Np; ++m) {
          if (elliptic->mapB[m + eM * mesh->Np] == 1) continue; //skip masked nodes
          A[n * mesh->Np + m] += elliptic->allNeumannPenalty * elliptic->allNeumannScale *
                                 elliptic->allNeumannScale;
        }
      }
  }
}

//returns the patch A matrix for element eM
void BuildLocalIpdgPatchAxHex3D(elliptic_t* elliptic, mesh_t* mesh, dfloat lambda,
                                dfloat* B, dfloat* Br, dfloat* Bs, dfloat* Bt, dlong eM, dfloat* A)
{
  /* start with stiffness matrix  */
  for(int n = 0; n < mesh->Np; ++n)
    for(int m = 0; m < mesh->Np; ++m) {
      A[n * mesh->Np + m] = 0;

      // (grad phi_n, grad phi_m)_{D^e}
      for(int i = 0; i < mesh->Np; ++i) {
        dlong base = eM * mesh->Np * mesh->Nvgeo + i;
        dfloat drdx = mesh->vgeo[base + mesh->Np * RXID];
        dfloat drdy = mesh->vgeo[base + mesh->Np * RYID];
        dfloat drdz = mesh->vgeo[base + mesh->Np * RZID];
        dfloat dsdx = mesh->vgeo[base + mesh->Np * SXID];
        dfloat dsdy = mesh->vgeo[base + mesh->Np * SYID];
        dfloat dsdz = mesh->vgeo[base + mesh->Np * SZID];
        dfloat dtdx = mesh->vgeo[base + mesh->Np * TXID];
        dfloat dtdy = mesh->vgeo[base + mesh->Np * TYID];
        dfloat dtdz = mesh->vgeo[base + mesh->Np * TZID];
        dfloat JW   = mesh->vgeo[base + mesh->Np * JWID];

        int idn = n * mesh->Np + i;
        int idm = m * mesh->Np + i;
        dfloat dlndx = drdx * Br[idn] + dsdx * Bs[idn] + dtdx * Bt[idn];
        dfloat dlndy = drdy * Br[idn] + dsdy * Bs[idn] + dtdy * Bt[idn];
        dfloat dlndz = drdz * Br[idn] + dsdz * Bs[idn] + dtdz * Bt[idn];
        dfloat dlmdx = drdx * Br[idm] + dsdx * Bs[idm] + dtdx * Bt[idm];
        dfloat dlmdy = drdy * Br[idm] + dsdy * Bs[idm] + dtdy * Bt[idm];
        dfloat dlmdz = drdz * Br[idm] + dsdz * Bs[idm] + dtdz * Bt[idm];
        A[n * mesh->Np + m] += JW * (dlndx * dlmdx + dlndy * dlmdy + dlndz * dlmdz);
        A[n * mesh->Np + m] += lambda * JW * B[idn] * B[idm];
      }

      for (int fM = 0; fM < mesh->Nfaces; fM++)
        // accumulate flux terms for negative and positive traces
        for(int i = 0; i < mesh->Nfp; ++i) {
          int vidM = mesh->faceNodes[i + fM * mesh->Nfp];

          // grab vol geofacs at surface nodes
          dlong baseM = eM * mesh->Np * mesh->Nvgeo + vidM;
          dfloat drdxM = mesh->vgeo[baseM + mesh->Np * RXID];
          dfloat drdyM = mesh->vgeo[baseM + mesh->Np * RYID];
          dfloat drdzM = mesh->vgeo[baseM + mesh->Np * RZID];
          dfloat dsdxM = mesh->vgeo[baseM + mesh->Np * SXID];
          dfloat dsdyM = mesh->vgeo[baseM + mesh->Np * SYID];
          dfloat dsdzM = mesh->vgeo[baseM + mesh->Np * SZID];
          dfloat dtdxM = mesh->vgeo[baseM + mesh->Np * TXID];
          dfloat dtdyM = mesh->vgeo[baseM + mesh->Np * TYID];
          dfloat dtdzM = mesh->vgeo[baseM + mesh->Np * TZID];

          // grab surface geometric factors
          dlong base = mesh->Nsgeo * (eM * mesh->Nfp * mesh->Nfaces + fM * mesh->Nfp + i);
          dfloat nx = mesh->sgeo[base + NXID];
          dfloat ny = mesh->sgeo[base + NYID];
          dfloat nz = mesh->sgeo[base + NZID];
          dfloat wsJ = mesh->sgeo[base + WSJID];
          dfloat hinv = mesh->sgeo[base + IHID];

          // form negative trace terms in IPDG
          int idnM = n * mesh->Np + vidM;
          int idmM = m * mesh->Np + vidM;

          dfloat dlndxM = drdxM * Br[idnM] + dsdxM * Bs[idnM] + dtdxM * Bt[idnM];
          dfloat dlndyM = drdyM * Br[idnM] + dsdyM * Bs[idnM] + dtdyM * Bt[idnM];
          dfloat dlndzM = drdzM * Br[idnM] + dsdzM * Bs[idnM] + dtdzM * Bt[idnM];
          dfloat ndotgradlnM = nx * dlndxM + ny * dlndyM + nz * dlndzM;
          dfloat lnM = B[idnM];

          dfloat dlmdxM = drdxM * Br[idmM] + dsdxM * Bs[idmM] + dtdxM * Bt[idmM];
          dfloat dlmdyM = drdyM * Br[idmM] + dsdyM * Bs[idmM] + dtdyM * Bt[idmM];
          dfloat dlmdzM = drdzM * Br[idmM] + dsdzM * Bs[idmM] + dtdzM * Bt[idmM];
          dfloat ndotgradlmM = nx * dlmdxM + ny * dlmdyM + nz * dlmdzM;
          dfloat lmM = B[idmM];

          dfloat penalty = elliptic->tau * hinv;
          int bc = mesh->EToB[fM + mesh->Nfaces * eM]; //raw boundary flag

          int bcD = 0, bcN = 0;
          int bcType = 0;

          if(bc > 0) bcType = elliptic->BCType[bc];        //find its type (Dirichlet/Neumann)

          // this needs to be double checked (and the code where these are used)
          if(bcType == 1) { // Dirichlet
            bcD = 1;
            bcN = 0;
          } else if(bcType == 2) { // Neumann
            bcD = 0;
            bcN = 1;
          }

          A[n * mesh->Np + m] += -0.5 * (1 + bcD) * (1 - bcN) * wsJ * lnM * ndotgradlmM;  // -(ln^-, N.grad lm^-)
          A[n * mesh->Np + m] += -0.5 * (1 + bcD) * (1 - bcN) * wsJ * ndotgradlnM * lmM;  // -(N.grad ln^-, lm^-)
          A[n * mesh->Np + m] += +0.5 * (1 + bcD) * (1 - bcN) * wsJ * penalty * lnM * lmM; // +((tau/h)*ln^-,lm^-)
        }
    }
}

void BuildLocalContinuousPatchAxHex3D(elliptic_t* elliptic,
                                      mesh_t* mesh,
                                      dfloat lambda,
                                      dlong eM,
                                      dfloat* B,
                                      dfloat* Br,
                                      dfloat* Bs,
                                      dfloat* Bt,
                                      dfloat* A)
{
  for (int nz = 0; nz < mesh->Nq; nz++)
    for (int ny = 0; ny < mesh->Nq; ny++)
      for (int nx = 0; nx < mesh->Nq; nx++) {
        int idn = nx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
        if (elliptic->mapB[idn + eM * mesh->Np] != 1) {
          for (int mz = 0; mz < mesh->Nq; mz++)
            for (int my = 0; my < mesh->Nq; my++)
              for (int mx = 0; mx < mesh->Nq; mx++) {
                int idm = mx + my * mesh->Nq + mz * mesh->Nq * mesh->Nq;
                int iid = idn * mesh->Np + idm;
                if (elliptic->mapB[idm + eM * mesh->Np] == 1) continue;

                int id;
                A[iid] = 0;

                if ((ny == my) && (nz == mz)) {
                  for (int k = 0; k < mesh->Nq; k++) {
                    id = k + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                    dfloat Grr = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G00ID * mesh->Np];

                    A[iid] += Grr * mesh->D[nx + k * mesh->Nq] * mesh->D[mx + k * mesh->Nq];
                  }
                }

                if (nz == mz) {
                  id = mx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                  dfloat Grs = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G01ID * mesh->Np];
                  A[iid] += Grs * mesh->D[nx + mx * mesh->Nq] * mesh->D[my + ny * mesh->Nq];

                  id = nx + my * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                  dfloat Gsr = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G01ID * mesh->Np];
                  A[iid] += Gsr * mesh->D[mx + nx * mesh->Nq] * mesh->D[ny + my * mesh->Nq];
                }

                if (ny == my) {
                  id = mx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                  dfloat Grt = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G02ID * mesh->Np];
                  A[iid] += Grt * mesh->D[nx + mx * mesh->Nq] * mesh->D[mz + nz * mesh->Nq];

                  id = nx + ny * mesh->Nq + mz * mesh->Nq * mesh->Nq;
                  dfloat Gst = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G02ID * mesh->Np];
                  A[iid] += Gst * mesh->D[mx + nx * mesh->Nq] * mesh->D[nz + mz * mesh->Nq];
                }

                if ((nx == mx) && (nz == mz)) {
                  for (int k = 0; k < mesh->Nq; k++) {
                    id = nx + k * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                    dfloat Gss = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G11ID * mesh->Np];

                    A[iid] += Gss * mesh->D[ny + k * mesh->Nq] * mesh->D[my + k * mesh->Nq];
                  }
                }

                if (nx == mx) {
                  id = nx + my * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                  dfloat Gst = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G12ID * mesh->Np];
                  A[iid] += Gst * mesh->D[ny + my * mesh->Nq] * mesh->D[mz + nz * mesh->Nq];

                  id = nx + ny * mesh->Nq + mz * mesh->Nq * mesh->Nq;
                  dfloat Gts = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G12ID * mesh->Np];
                  A[iid] += Gts * mesh->D[my + ny * mesh->Nq] * mesh->D[nz + mz * mesh->Nq];
                }

                if ((nx == mx) && (ny == my)) {
                  for (int k = 0; k < mesh->Nq; k++) {
                    id = nx + ny * mesh->Nq + k * mesh->Nq * mesh->Nq;
                    dfloat Gtt = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + G22ID * mesh->Np];

                    A[iid] += Gtt * mesh->D[nz + k * mesh->Nq] * mesh->D[mz + k * mesh->Nq];
                  }
                }

                if ((nx == mx) && (ny == my) && (nz == mz)) {
                  id = nx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                  dfloat JW = mesh->ggeo[eM * mesh->Np * mesh->Nggeo + id + GWJID * mesh->Np];
                  A[iid] += JW * lambda;
                }
              }
        } else {
          int iid = idn * mesh->Np + idn;
          A[iid] = 1; //just put a 1 so A is invertable
        }
      }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n = 0; n < mesh->Np; ++n)
      if (elliptic->mapB[n + eM * mesh->Np] != 1) { //dont fill rows for masked nodes
        for(int m = 0; m < mesh->Np; ++m) {
          if (elliptic->mapB[m + eM * mesh->Np] == 1) continue; //skip masked nodes
          A[n * mesh->Np + m] += elliptic->allNeumannPenalty * elliptic->allNeumannScale *
                                 elliptic->allNeumannScale;
        }
      }
  }
}