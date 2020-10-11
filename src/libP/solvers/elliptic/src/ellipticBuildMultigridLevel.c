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

// create elliptic and mesh structs for multigrid levels
elliptic_t* ellipticBuildMultigridLevel(elliptic_t* baseElliptic, int Nc, int Nf)
{
  const int serial = baseElliptic->options.compareArgs("THREAD MODEL", "SERIAL");

  elliptic_t* elliptic = new elliptic_t();

  memcpy(elliptic,baseElliptic,sizeof(elliptic_t));

  int buildOnly  = 0;
  if(elliptic->options.compareArgs("BUILD ONLY", "TRUE")) buildOnly = 1;

  //populate the mini-mesh using the mesh struct
  mesh_t* mesh = new mesh_t();

#ifndef OCCA_VERSION_1_0
  memcpy(mesh,baseElliptic->mesh,sizeof(mesh_t));
  fflush(stdout);

#else

  mesh->rank = baseElliptic->mesh->rank;
  mesh->size = baseElliptic->mesh->size;

  MPI_Comm_dup(baseElliptic->mesh->comm, &(mesh->comm));

  mesh->dim = baseElliptic->mesh->dim;
  mesh->Nverts        = baseElliptic->mesh->Nverts;
  mesh->Nfaces        = baseElliptic->mesh->Nfaces;

  mesh->NfaceVertices = baseElliptic->mesh->NfaceVertices;

  mesh->Nfields = baseElliptic->mesh->Nfields;

  mesh->Nnodes = baseElliptic->mesh->Nnodes;
  mesh->EX = baseElliptic->mesh->EX; // coordinates of vertices for each element
  mesh->EY = baseElliptic->mesh->EY;
  mesh->EZ = baseElliptic->mesh->EZ;

  mesh->Nelements = baseElliptic->mesh->Nelements;
  mesh->EToV = baseElliptic->mesh->EToV; // element-to-vertex connectivity
  mesh->EToE = baseElliptic->mesh->EToE; // element-to-element connectivity
  mesh->EToF = baseElliptic->mesh->EToF; // element-to-(local)face connectivity
  mesh->EToP = baseElliptic->mesh->EToP; // element-to-partition/process connectivity
  mesh->EToB = baseElliptic->mesh->EToB; // element-to-boundary condition type

  mesh->elementInfo = baseElliptic->mesh->elementInfo;

  // boundary faces
  mesh->NboundaryFaces = baseElliptic->mesh->NboundaryFaces;
  mesh->boundaryInfo = baseElliptic->mesh->boundaryInfo;

  // MPI halo exchange info
  mesh->totalHaloPairs = baseElliptic->mesh->totalHaloPairs;
  mesh->haloElementList = baseElliptic->mesh->haloElementList;
  mesh->NhaloPairs = baseElliptic->mesh->NhaloPairs;
  mesh->NhaloMessages = baseElliptic->mesh->NhaloMessages;

  mesh->haloSendRequests = baseElliptic->mesh->haloSendRequests;
  mesh->haloRecvRequests = baseElliptic->mesh->haloRecvRequests;

  mesh->NinternalElements = baseElliptic->mesh->NinternalElements;
  mesh->NnotInternalElements = baseElliptic->mesh->NnotInternalElements;

  mesh->o_haloElementList = baseElliptic->mesh->o_haloElementList;
  mesh->o_haloBuffer      = baseElliptic->mesh->o_haloBuffer;
  mesh->o_internalElementIds    = baseElliptic->mesh->o_internalElementIds;
  mesh->o_notInternalElementIds = baseElliptic->mesh->o_notInternalElementIds;

  // volumeGeometricFactors;
  mesh->Nvgeo = baseElliptic->mesh->Nvgeo;
  mesh->vgeo = baseElliptic->mesh->vgeo;
  mesh->o_vgeo = baseElliptic->mesh->o_vgeo;

  // second order volume geometric factors
  mesh->Nggeo = baseElliptic->mesh->Nggeo;
  mesh->ggeo = baseElliptic->mesh->ggeo;
  mesh->o_ggeo = baseElliptic->mesh->o_ggeo;

  mesh->Nsgeo = baseElliptic->mesh->Nsgeo;
  mesh->sgeo = baseElliptic->mesh->sgeo;
  mesh->o_sgeo = baseElliptic->mesh->o_sgeo;

  // occa stuff
  mesh->device = baseElliptic->mesh->device;

  mesh->defaultStream = baseElliptic->mesh->defaultStream;
  mesh->dataStream = baseElliptic->mesh->dataStream;

  mesh->haloExtractKernel = baseElliptic->mesh->haloExtractKernel;
  mesh->addScalarKernel = baseElliptic->mesh->addScalarKernel;
  mesh->maskKernel = baseElliptic->mesh->maskKernel;
  mesh->sumKernel = baseElliptic->mesh->sumKernel;
#endif

  elliptic->mesh = mesh;

  setupAide options = elliptic->options;

  switch(elliptic->elementType) {
  case TRIANGLES:
    meshLoadReferenceNodesTri2D(mesh, Nc);
    if(elliptic->dim == 2) {
      meshPhysicalNodesTri2D(mesh);
      meshGeometricFactorsTri2D(mesh);
    }else{
      meshPhysicalNodesTri3D(mesh);
      meshGeometricFactorsTri3D(mesh);
    }
    break;
  case QUADRILATERALS: {
    meshLoadReferenceNodesQuad2D(mesh, Nc);
    if(elliptic->dim == 2) {
      meshPhysicalNodesQuad2D(mesh);
      meshGeometricFactorsQuad2D(mesh);
    }else{
      meshPhysicalNodesQuad3D(mesh);
      meshGeometricFactorsQuad3D(mesh);
    }
  } break;
  case TETRAHEDRA:
    meshLoadReferenceNodesTet3D(mesh, Nc);
    meshPhysicalNodesTet3D(mesh);
    break;
  case HEXAHEDRA:
    meshLoadReferenceNodesHex3D(mesh, Nc, 1);
    meshPhysicalNodesHex3D(mesh, buildOnly);
    meshGeometricFactorsHex3D(mesh);
    break;
  }

  // create halo extension for x,y arrays
  dlong totalHaloNodes = mesh->totalHaloPairs * mesh->Np;
  dlong localNodes     = mesh->Nelements * mesh->Np;
  // temporary send buffer
  dfloat* sendBuffer = (dfloat*) calloc(totalHaloNodes, sizeof(dfloat));

  // extend x,y arrays to hold coordinates of node coordinates of elements in halo
  mesh->x = (dfloat*) realloc(mesh->x, (localNodes + totalHaloNodes) * sizeof(dfloat));
  mesh->y = (dfloat*) realloc(mesh->y, (localNodes + totalHaloNodes) * sizeof(dfloat));
  meshHaloExchange(mesh, mesh->Np * sizeof(dfloat), mesh->x, sendBuffer, mesh->x + localNodes);
  meshHaloExchange(mesh, mesh->Np * sizeof(dfloat), mesh->y, sendBuffer, mesh->y + localNodes);

  if (elliptic->dim == 3) {
    mesh->z = (dfloat*) realloc(mesh->z, (localNodes + totalHaloNodes) * sizeof(dfloat));
    meshHaloExchange(mesh, mesh->Np * sizeof(dfloat), mesh->z, sendBuffer, mesh->z + localNodes);
  }

  switch(elliptic->elementType) {
  case TRIANGLES:
    if(elliptic->dim == 2)
      meshConnectFaceNodes2D(mesh);
    else
      meshConnectFaceNodes3D(mesh);
    break;
  case QUADRILATERALS: {
    if(elliptic->dim == 2) {
      if(!options.compareArgs("BOX DOMAIN", "TRUE")) {
        meshConnectFaceNodes2D(mesh);
        meshSurfaceGeometricFactorsQuad2D(mesh);
      }else {
        if(mesh->rank == 0) printf("WARNING: connecting periodic box\n");
        dfloat XMIN = -1, XMAX = +1; // default bi-unit cube
        dfloat YMIN = -1, YMAX = +1;
        options.getArgs("BOX XMIN", XMIN);
        options.getArgs("BOX YMIN", YMIN);
        options.getArgs("BOX XMAX", XMAX);
        options.getArgs("BOX YMAX", YMAX);
        meshConnectPeriodicFaceNodes2D(mesh, XMAX - XMIN, YMAX - YMIN);
        meshSurfaceGeometricFactorsQuad2D(mesh);
      }
    }else{
      meshConnectFaceNodes3D(mesh);
      meshSurfaceGeometricFactorsQuad3D(mesh);
    }
  } break;
  case TETRAHEDRA:
    meshConnectFaceNodes3D(mesh);
    break;
  case HEXAHEDRA:

    if(!options.compareArgs("BOX DOMAIN", "TRUE")) {
      meshConnectFaceNodes3D(mesh);
    }else {
      if(mesh->rank == 0)
        printf("WARNING: connecting periodic box\n");

      dfloat XMIN = -1, XMAX = +1; // default bi-unit cube
      dfloat YMIN = -1, YMAX = +1;
      dfloat ZMIN = -1, ZMAX = +1;

      options.getArgs("BOX XMIN", XMIN);
      options.getArgs("BOX YMIN", YMIN);
      options.getArgs("BOX ZMIN", ZMIN);

      options.getArgs("BOX XMAX", XMAX);
      options.getArgs("BOX YMAX", YMAX);
      options.getArgs("BOX ZMAX", ZMAX);

      meshConnectPeriodicFaceNodes3D(mesh, XMAX - XMIN, YMAX - YMIN, ZMAX - ZMIN);
    }
    meshSurfaceGeometricFactorsHex3D(mesh);
    break;
  }

  // global nodes
  meshParallelConnectNodes(mesh, 0, buildOnly);

  //dont need these once vmap is made
  free(mesh->x);
  free(mesh->y);
  if (elliptic->dim == 3)
    free(mesh->z);
  free(sendBuffer);

  dlong Ntotal = mesh->Np * mesh->Nelements;
  dlong Nblock = mymax(1,(Ntotal + blockSize - 1) / blockSize);
  dlong Nblock2 = mymax(1,(Nblock + blockSize - 1) / blockSize);

  elliptic->Nblock = Nblock;
  elliptic->Nblock2 = Nblock2;

  if (elliptic->elementType == TRIANGLES) {
    // build Dr, Ds, LIFT transposes
    dfloat* DrT = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    dfloat* DsT = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        DrT[n + m * mesh->Np] = mesh->Dr[n * mesh->Np + m];
        DsT[n + m * mesh->Np] = mesh->Ds[n * mesh->Np + m];
      }

    // build Dr, Ds transposes
    dfloat* DrsT = (dfloat*) calloc(2 * mesh->Np * mesh->Np, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        DrsT[n + m * mesh->Np] = mesh->Dr[n * mesh->Np + m];
        DrsT[n + m * mesh->Np + mesh->Np * mesh->Np] = mesh->Ds[n * mesh->Np + m];
      }

    dfloat* LIFTT = (dfloat*) calloc(mesh->Np * mesh->Nfaces * mesh->Nfp, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Nfaces * mesh->Nfp; ++m)
        LIFTT[n + m * mesh->Np] = mesh->LIFT[n * mesh->Nfp * mesh->Nfaces + m];

    //build element stiffness matrices
    dfloat* SrrT, * SrsT, * SsrT, * SssT;
    mesh->Srr = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Srs = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Ssr = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Sss = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    for (int n = 0; n < mesh->Np; n++)
      for (int m = 0; m < mesh->Np; m++)
        for (int k = 0; k < mesh->Np; k++)
          for (int l = 0; l < mesh->Np; l++) {
            mesh->Srr[m + n * mesh->Np] += mesh->Dr[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dr[m + k * mesh->Np];
            mesh->Srs[m + n * mesh->Np] += mesh->Dr[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Ds[m + k * mesh->Np];
            mesh->Ssr[m + n * mesh->Np] += mesh->Ds[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dr[m + k * mesh->Np];
            mesh->Sss[m + n * mesh->Np] += mesh->Ds[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Ds[m + k * mesh->Np];
          }
    SrrT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    SrsT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    SsrT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    SssT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    for (int n = 0; n < mesh->Np; n++)
      for (int m = 0; m < mesh->Np; m++) {
        SrrT[m + n * mesh->Np] = mesh->Srr[n + m * mesh->Np];
        SrsT[m + n * mesh->Np] = mesh->Srs[n + m * mesh->Np];
        SsrT[m + n * mesh->Np] = mesh->Ssr[n + m * mesh->Np];
        SssT[m + n * mesh->Np] = mesh->Sss[n + m * mesh->Np];
      }

    dfloat* ST = (dfloat*) calloc(3 * mesh->Np * mesh->Np, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        ST[n + m * mesh->Np + 0 * mesh->Np * mesh->Np] = mesh->Srr[n * mesh->Np + m];
        ST[n + m * mesh->Np + 1 * mesh->Np * mesh->Np] = mesh->Srs[n * mesh->Np + m] +
                                                         mesh->Ssr[n * mesh->Np + m];
        ST[n + m * mesh->Np + 2 * mesh->Np * mesh->Np] = mesh->Sss[n * mesh->Np + m];
      }

    // deriv operators: transpose from row major to column major
    int* D1ids = (int*) calloc(mesh->Np * 3,sizeof(int));
    int* D2ids = (int*) calloc(mesh->Np * 3,sizeof(int));
    int* D3ids = (int*) calloc(mesh->Np * 3,sizeof(int));
    dfloat* Dvals = (dfloat*) calloc(mesh->Np * 3,sizeof(dfloat));

    dfloat* VBq = (dfloat*) calloc(mesh->Np * mesh->cubNp,sizeof(dfloat));
    dfloat* PBq = (dfloat*) calloc(mesh->Np * mesh->cubNp,sizeof(dfloat));

    dfloat* L0vals = (dfloat*) calloc(mesh->Nfp * 3,sizeof(dfloat)); // tridiag
    int* ELids = (int*) calloc(1 + mesh->Np * mesh->max_EL_nnz,sizeof(int));
    dfloat* ELvals = (dfloat*) calloc(1 + mesh->Np * mesh->max_EL_nnz,sizeof(dfloat));

    for (int i = 0; i < mesh->Np; ++i)
      for (int j = 0; j < 3; ++j) {
        D1ids[i + j * mesh->Np] = mesh->D1ids[j + i * 3];
        D2ids[i + j * mesh->Np] = mesh->D2ids[j + i * 3];
        D3ids[i + j * mesh->Np] = mesh->D3ids[j + i * 3];
        Dvals[i + j * mesh->Np] = mesh->Dvals[j + i * 3];
      }

    for (int i = 0; i < mesh->cubNp; ++i)
      for (int j = 0; j < mesh->Np; ++j) {
        VBq[i + j * mesh->cubNp] = mesh->VBq[j + i * mesh->Np];
        PBq[j + i * mesh->Np] = mesh->PBq[i + j * mesh->cubNp];
      }


    for (int i = 0; i < mesh->Nfp; ++i)
      for (int j = 0; j < 3; ++j)
        L0vals[i + j * mesh->Nfp] = mesh->L0vals[j + i * 3];

    for (int i = 0; i < mesh->Np; ++i)
      for (int j = 0; j < mesh->max_EL_nnz; ++j) {
        ELids[i + j * mesh->Np] = mesh->ELids[j + i * mesh->max_EL_nnz];
        ELvals[i + j * mesh->Np] = mesh->ELvals[j + i * mesh->max_EL_nnz];
      }

    //BB mass matrix
    mesh->BBMM = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    for (int n = 0; n < mesh->Np; ++n)
      for (int m = 0; m < mesh->Np; ++m)
        for (int i = 0; i < mesh->Np; ++i)
          for (int j = 0; j < mesh->Np; ++j)
            mesh->BBMM[n + m * mesh->Np] += mesh->VB[m + j * mesh->Np] *
                                            mesh->MM[i + j * mesh->Np] * mesh->VB[n + i * mesh->Np];

    mesh->o_Dr = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                                     mesh->Dr);

    mesh->o_Ds = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                                     mesh->Ds);

    mesh->o_DrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                                      DrT);

    mesh->o_DsT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                                      DsT);

    mesh->o_Dmatrices = mesh->device.malloc(2 * mesh->Np * mesh->Np * sizeof(dfloat), DrsT);

    mesh->o_LIFT =
      mesh->device.malloc(mesh->Np * mesh->Nfaces * mesh->Nfp * sizeof(dfloat),
                          mesh->LIFT);

    mesh->o_LIFTT =
      mesh->device.malloc(mesh->Np * mesh->Nfaces * mesh->Nfp * sizeof(dfloat),
                          LIFTT);

    mesh->o_SrrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SrrT);
    mesh->o_SrsT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SrsT);
    mesh->o_SsrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SsrT);
    mesh->o_SssT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SssT);

    mesh->o_Smatrices = mesh->device.malloc(3 * mesh->Np * mesh->Np * sizeof(dfloat), ST);

    mesh->o_D1ids = mesh->device.malloc(mesh->Np * 3 * sizeof(int),D1ids);
    mesh->o_D2ids = mesh->device.malloc(mesh->Np * 3 * sizeof(int),D2ids);
    mesh->o_D3ids = mesh->device.malloc(mesh->Np * 3 * sizeof(int),D3ids);
    mesh->o_Dvals = mesh->device.malloc(mesh->Np * 3 * sizeof(dfloat),Dvals);

    mesh->o_BBMM = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),mesh->BBMM);

    mesh->o_VBq = mesh->device.malloc(mesh->Np * mesh->cubNp * sizeof(dfloat),VBq);
    mesh->o_PBq = mesh->device.malloc(mesh->Np * mesh->cubNp * sizeof(dfloat),PBq);

    mesh->o_L0vals = mesh->device.malloc(mesh->Nfp * 3 * sizeof(dfloat),L0vals);
    mesh->o_ELids =
      mesh->device.malloc(mesh->Np * mesh->max_EL_nnz * sizeof(int),ELids);
    mesh->o_ELvals =
      mesh->device.malloc(mesh->Np * mesh->max_EL_nnz * sizeof(dfloat),ELvals);

    free(DrT);
    free(DsT);
    free(LIFTT);
    free(SrrT);
    free(SrsT);
    free(SsrT);
    free(SssT);
    free(D1ids);
    free(D2ids);
    free(D3ids);
    free(Dvals);

    free(VBq);
    free(PBq);
    free(L0vals);
    free(ELids );
    free(ELvals);
  } else if (elliptic->elementType == QUADRILATERALS) { // This should be ok for Quad3D
    //lumped mass matrix
    mesh->MM = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    for (int j = 0; j < mesh->Nq; j++)
      for (int i = 0; i < mesh->Nq; i++) {
        int n = i + j * mesh->Nq;
        mesh->MM[n + n * mesh->Np] = mesh->gllw[i] * mesh->gllw[j];
      }

    mesh->o_D = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);
    mesh->o_Dmatrices = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);
    mesh->o_Smatrices = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D); //dummy

    mesh->o_vgeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nvgeo * mesh->Np * sizeof(dfloat),
                          mesh->vgeo);
    mesh->o_sgeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->Nfp * mesh->Nsgeo * sizeof(dfloat),
                          mesh->sgeo);
    mesh->o_ggeo =
      mesh->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo * sizeof(dfloat),
                          mesh->ggeo);

    mesh->o_LIFTT = baseElliptic->mesh->o_LIFTT; //dummy buffer
  } else if (elliptic->elementType == TETRAHEDRA) {
    // build Dr, Ds, LIFT transposes
    dfloat* DrT = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    dfloat* DsT = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    dfloat* DtT = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        DrT[n + m * mesh->Np] = mesh->Dr[n * mesh->Np + m];
        DsT[n + m * mesh->Np] = mesh->Ds[n * mesh->Np + m];
        DtT[n + m * mesh->Np] = mesh->Dt[n * mesh->Np + m];
      }

    // build Dr, Ds transposes
    dfloat* DrstT = (dfloat*) calloc(3 * mesh->Np * mesh->Np, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        DrstT[n + m * mesh->Np] = mesh->Dr[n * mesh->Np + m];
        DrstT[n + m * mesh->Np + mesh->Np * mesh->Np] = mesh->Ds[n * mesh->Np + m];
        DrstT[n + m * mesh->Np + 2 * mesh->Np * mesh->Np] = mesh->Dt[n * mesh->Np + m];
      }

    dfloat* LIFTT = (dfloat*) calloc(mesh->Np * mesh->Nfaces * mesh->Nfp, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Nfaces * mesh->Nfp; ++m)
        LIFTT[n + m * mesh->Np] = mesh->LIFT[n * mesh->Nfp * mesh->Nfaces + m];

    //build element stiffness matrices
    mesh->Srr = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Srs = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Srt = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Ssr = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Sss = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Sst = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Str = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Sts = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    mesh->Stt = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    for (int n = 0; n < mesh->Np; n++)
      for (int m = 0; m < mesh->Np; m++)
        for (int k = 0; k < mesh->Np; k++)
          for (int l = 0; l < mesh->Np; l++) {
            mesh->Srr[m + n * mesh->Np] += mesh->Dr[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dr[m + k * mesh->Np];
            mesh->Srs[m + n * mesh->Np] += mesh->Dr[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Ds[m + k * mesh->Np];
            mesh->Srt[m + n * mesh->Np] += mesh->Dr[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dt[m + k * mesh->Np];
            mesh->Ssr[m + n * mesh->Np] += mesh->Ds[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dr[m + k * mesh->Np];
            mesh->Sss[m + n * mesh->Np] += mesh->Ds[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Ds[m + k * mesh->Np];
            mesh->Sst[m + n * mesh->Np] += mesh->Ds[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dt[m + k * mesh->Np];
            mesh->Str[m + n * mesh->Np] += mesh->Dt[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dr[m + k * mesh->Np];
            mesh->Sts[m + n * mesh->Np] += mesh->Dt[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Ds[m + k * mesh->Np];
            mesh->Stt[m + n * mesh->Np] += mesh->Dt[n + l * mesh->Np] * mesh->MM[k + l * mesh->Np] *
                                           mesh->Dt[m + k * mesh->Np];
          }
    dfloat* SrrT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    dfloat* SrsT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    dfloat* SrtT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    dfloat* SsrT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    dfloat* SssT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    dfloat* SstT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    dfloat* StrT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    dfloat* StsT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    dfloat* SttT = (dfloat*) calloc(mesh->Np * mesh->Np,sizeof(dfloat));
    for (int n = 0; n < mesh->Np; n++)
      for (int m = 0; m < mesh->Np; m++) {
        SrrT[m + n * mesh->Np] = mesh->Srr[n + m * mesh->Np];
        SrsT[m + n * mesh->Np] = mesh->Srs[n + m * mesh->Np] + mesh->Ssr[n + m * mesh->Np];
        SrtT[m + n * mesh->Np] = mesh->Srt[n + m * mesh->Np] + mesh->Str[n + m * mesh->Np];
        SssT[m + n * mesh->Np] = mesh->Sss[n + m * mesh->Np];
        SstT[m + n * mesh->Np] = mesh->Sst[n + m * mesh->Np] + mesh->Sts[n + m * mesh->Np];
        SttT[m + n * mesh->Np] = mesh->Stt[n + m * mesh->Np];
      }
    dfloat* ST = (dfloat*) calloc(6 * mesh->Np * mesh->Np, sizeof(dfloat));
    for(int n = 0; n < mesh->Np; ++n)
      for(int m = 0; m < mesh->Np; ++m) {
        ST[n + m * mesh->Np + 0 * mesh->Np * mesh->Np] = mesh->Srr[n * mesh->Np + m];
        ST[n + m * mesh->Np + 1 * mesh->Np * mesh->Np] = mesh->Srs[n * mesh->Np + m] +
                                                         mesh->Ssr[n * mesh->Np + m];
        ST[n + m * mesh->Np + 2 * mesh->Np * mesh->Np] = mesh->Srt[n * mesh->Np + m] +
                                                         mesh->Str[n * mesh->Np + m];
        ST[n + m * mesh->Np + 3 * mesh->Np * mesh->Np] = mesh->Sss[n * mesh->Np + m];
        ST[n + m * mesh->Np + 4 * mesh->Np * mesh->Np] = mesh->Sst[n * mesh->Np + m] +
                                                         mesh->Sts[n * mesh->Np + m];
        ST[n + m * mesh->Np + 5 * mesh->Np * mesh->Np] = mesh->Stt[n * mesh->Np + m];
      }

    mesh->o_Dr = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                                     mesh->Dr);

    mesh->o_Ds = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                                     mesh->Ds);

    mesh->o_Dt = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                                     mesh->Dt);

    mesh->o_DrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                                      DrT);

    mesh->o_DsT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                                      DsT);

    mesh->o_DtT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat),
                                      DtT);

    mesh->o_Dmatrices = mesh->device.malloc(3 * mesh->Np * mesh->Np * sizeof(dfloat), DrstT);

    mesh->o_LIFT =
      mesh->device.malloc(mesh->Np * mesh->Nfaces * mesh->Nfp * sizeof(dfloat),
                          mesh->LIFT);

    mesh->o_LIFTT =
      mesh->device.malloc(mesh->Np * mesh->Nfaces * mesh->Nfp * sizeof(dfloat),
                          LIFTT);

    mesh->o_SrrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SrrT);
    mesh->o_SrsT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SrsT);
    mesh->o_SrtT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SrtT);
    mesh->o_SsrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SsrT);
    mesh->o_SssT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SssT);
    mesh->o_SstT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SstT);
    mesh->o_StrT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), StrT);
    mesh->o_StsT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), StsT);
    mesh->o_SttT = mesh->device.malloc(mesh->Np * mesh->Np * sizeof(dfloat), SttT);

    mesh->o_Smatrices = mesh->device.malloc(6 * mesh->Np * mesh->Np * sizeof(dfloat), ST);

    free(DrT);
    free(DsT);
    free(DtT);
    free(LIFTT);
    free(SrrT);
    free(SrsT);
    free(SrtT);
    free(SsrT);
    free(SssT);
    free(SstT);
    free(StrT);
    free(StsT);
    free(SttT);
  } else if (elliptic->elementType == HEXAHEDRA) {
    //lumped mass matrix
    mesh->MM = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
    dfloat* DT = (dfloat*) calloc(mesh->Nq * mesh->Nq, sizeof(dfloat));

    for (int j = 0; j < mesh->Nq; j++)
      for (int i = 0; i < mesh->Nq; i++)
        DT[j * mesh->Nq + i] = mesh->D[i * mesh->Nq + j];

    for (int k = 0; k < mesh->Nq; k++)
      for (int j = 0; j < mesh->Nq; j++)
        for (int i = 0; i < mesh->Nq; i++) {
          int n = i + j * mesh->Nq + k * mesh->Nq * mesh->Nq;
          mesh->MM[n + n * mesh->Np] = mesh->gllw[i] * mesh->gllw[j] * mesh->gllw[k];
        }

    mesh->o_D = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);
    mesh->o_Dmatrices = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);
    mesh->o_Smatrices = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), DT); // transpose(D)

    mesh->o_cubD = mesh->device.malloc(mesh->cubNq * mesh->cubNq * sizeof(dfloat), mesh->cubD);

    dfloat* cubInterpT = (dfloat*) calloc(mesh->cubNq * mesh->Nq, sizeof(dfloat));
    for(int n = 0; n < mesh->Nq; ++n)
      for(int m = 0; m < mesh->cubNq; ++m)
        cubInterpT[m + n * mesh->cubNq] = mesh->cubInterp[m * mesh->Nq + n];

    mesh->o_cubInterpT = mesh->device.malloc(mesh->cubNq * mesh->Nq * sizeof(dfloat), cubInterpT);

    free(cubInterpT);

    mesh->o_vgeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nvgeo * mesh->Np * sizeof(dfloat),
                          mesh->vgeo);

    mesh->o_sgeo =
      mesh->device.malloc(mesh->Nelements * mesh->Nfaces * mesh->Nfp * mesh->Nsgeo * sizeof(dfloat),
                          mesh->sgeo);

    mesh->o_ggeo =
      mesh->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo * sizeof(dfloat),
                          mesh->ggeo);

    mesh->o_cubggeo = mesh->o_ggeo; // dummy
    if(options.compareArgs("ELLIPTIC INTEGRATION", "CUBATURE"))
      mesh->o_cubggeo =
        mesh->device.malloc(mesh->Nelements * mesh->cubNp * mesh->Nggeo * sizeof(dfloat),
                            mesh->cubggeo);

    mesh->o_vmapM =
      mesh->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(dlong),
                          mesh->vmapM);

    mesh->o_vmapP =
      mesh->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(dlong),
                          mesh->vmapP);

    mesh->LIFT = baseElliptic->mesh->LIFT; //dummy buffer
    mesh->o_LIFTT = baseElliptic->mesh->o_LIFTT; //dummy buffer
  }

  //fill geometric factors in halo
  if(mesh->totalHaloPairs &&
     (elliptic->elementType == QUADRILATERALS || elliptic->elementType == HEXAHEDRA)) {
    dlong Nlocal = mesh->Np * mesh->Nelements;
    dlong Nhalo = mesh->totalHaloPairs * mesh->Np;
    dfloat* vgeoSendBuffer = (dfloat*) calloc(Nhalo * mesh->Nvgeo, sizeof(dfloat));

    // import geometric factors from halo elements
    mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (Nlocal + Nhalo) * mesh->Nvgeo * sizeof(dfloat));

    meshHaloExchange(mesh,
                     mesh->Nvgeo * mesh->Np * sizeof(dfloat),
                     mesh->vgeo,
                     vgeoSendBuffer,
                     mesh->vgeo + Nlocal * mesh->Nvgeo);

    mesh->o_vgeo =
      mesh->device.malloc((Nlocal + Nhalo) * mesh->Nvgeo * sizeof(dfloat), mesh->vgeo);
    free(vgeoSendBuffer);
  }

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(int),
                        mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(int),
                        mesh->vmapP);

  //set the normalization constant for the allNeumann Poisson problem on this coarse mesh
  hlong localElements = (hlong) mesh->Nelements;
  hlong totalElements = 0;
  MPI_Allreduce(&localElements, &totalElements, 1, MPI_HLONG, MPI_SUM, mesh->comm);
  elliptic->allNeumannScale = 1.0 / sqrt(mesh->Np * totalElements);

  elliptic->allNeumannPenalty = 0;
  elliptic->allNeumannScale = 0;

  elliptic->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));
  elliptic->o_tmp = mesh->device.malloc(Nblock * sizeof(dfloat), elliptic->tmp);
  elliptic->o_tmp2 = mesh->device.malloc(Nblock2 * sizeof(dfloat), elliptic->tmp);

  //tau
  if (elliptic->elementType == TRIANGLES ||
      elliptic->elementType == QUADRILATERALS) {
    elliptic->tau = 2.0 * (mesh->N + 1) * (mesh->N + 2) / 2.0;
    if(elliptic->dim == 3)
      elliptic->tau *= 1.5;
  }else {
    elliptic->tau = 2.0 * (mesh->N + 1) * (mesh->N + 3);
  }

  //setup an unmasked gs handle
  int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;
  meshParallelGatherScatterSetup(mesh, Ntotal, mesh->globalIds, mesh->comm, verbose);

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  elliptic->mapB = (int*) calloc(mesh->Nelements * mesh->Np,sizeof(int));
  for (dlong e = 0; e < mesh->Nelements; e++) {
    for (int n = 0; n < mesh->Np; n++) elliptic->mapB[n + e * mesh->Np] = 1E9;
    for (int f = 0; f < mesh->Nfaces; f++) {
      int bc = mesh->EToB[f + e * mesh->Nfaces];
      if (bc > 0) {
        for (int n = 0; n < mesh->Nfp; n++) {
          int BCFlag = elliptic->BCType[bc];
          int fid = mesh->faceNodes[n + f * mesh->Nfp];
          elliptic->mapB[fid + e * mesh->Np] = mymin(BCFlag,elliptic->mapB[fid + e * mesh->Np]);
        }
      }
    }
  }
  ogsGatherScatter(elliptic->mapB, ogsInt, ogsMin, mesh->ogs);

  //use the bc flags to find masked ids
  elliptic->Nmasked = 0;
  for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++) {
    if (elliptic->mapB[n] == 1E9)
      elliptic->mapB[n] = 0.;
    else if (elliptic->mapB[n] == 1)     //Dirichlet boundary
      elliptic->Nmasked++;
  }
  elliptic->o_mapB = mesh->device.malloc(mesh->Nelements * mesh->Np * sizeof(int), elliptic->mapB);

  elliptic->maskIds = (dlong*) calloc(elliptic->Nmasked, sizeof(dlong));
  elliptic->Nmasked = 0; //reset
  for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++)
    if (elliptic->mapB[n] == 1) elliptic->maskIds[elliptic->Nmasked++] = n;
  if (elliptic->Nmasked) elliptic->o_maskIds = mesh->device.malloc(
      elliptic->Nmasked * sizeof(dlong),
      elliptic->maskIds);

  //make a masked version of the global id numbering
  mesh->maskedGlobalIds = (hlong*) calloc(Ntotal,sizeof(hlong));
  memcpy(mesh->maskedGlobalIds, mesh->globalIds, Ntotal * sizeof(hlong));
  for (dlong n = 0; n < elliptic->Nmasked; n++)
    mesh->maskedGlobalIds[elliptic->maskIds[n]] = 0;

  //use the masked ids to make another gs handle
  elliptic->ogs = ogsSetup(Ntotal, mesh->maskedGlobalIds, mesh->comm, verbose, mesh->device);
  elliptic->o_invDegree = elliptic->ogs->o_invDegree;

  // HERE
  occa::properties kernelInfo = ellipticKernelInfo(mesh);

  //  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  // set kernel name suffix
  char* suffix;

  if(elliptic->elementType == TRIANGLES) {
    if(elliptic->dim == 2)
      suffix = strdup("Tri2D");
    else
      suffix = strdup("Tri3D");
  }
  if(elliptic->elementType == QUADRILATERALS) {
    if(elliptic->dim == 2)
      suffix = strdup("Quad2D");
    if(elliptic->dim == 3)
      suffix = strdup("Quad3D");
  }
  if(elliptic->elementType == TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(elliptic->elementType == HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  MPI_Barrier(mesh->comm);
  double tStartLoadKernel = MPI_Wtime();
  if(mesh->rank == 0)  printf("loading elliptic MG kernels ... "); fflush(stdout); 

  for (int r = 0; r < 2; r++) {
    MPI_Barrier(mesh->comm);

    if ((r == 0 && mesh->rank == 0) || (r == 1 && mesh->rank > 0)) {
      kernelInfo["defines/" "p_blockSize"] = blockSize;

      // add custom defines
      kernelInfo["defines/" "p_NpP"] = (mesh->Np + mesh->Nfp * mesh->Nfaces);
      kernelInfo["defines/" "p_Nverts"] = mesh->Nverts;

      int Nmax = mymax(mesh->Np, mesh->Nfaces * mesh->Nfp);
      kernelInfo["defines/" "p_Nmax"] = Nmax;

      int maxNodes = mymax(mesh->Np, (mesh->Nfp * mesh->Nfaces));
      kernelInfo["defines/" "p_maxNodes"] = maxNodes;

      int NblockV = mymax(1,maxNthreads / mesh->Np); // works for CUDA
      kernelInfo["defines/" "p_NblockV"] = NblockV;

      int one = 1; //set to one for now. TODO: try optimizing over these
      kernelInfo["defines/" "p_NnodesV"] = one;

      int NblockS = mymax(1,maxNthreads / maxNodes); // works for CUDA
      kernelInfo["defines/" "p_NblockS"] = NblockS;

      int NblockP = mymax(1,maxNthreads / (4 * mesh->Np)); // get close to maxNthreads threads
      kernelInfo["defines/" "p_NblockP"] = NblockP;

      int NblockG;
      if(mesh->Np <= 32) NblockG = ( 32 / mesh->Np );
      else NblockG = mymax(1,maxNthreads / mesh->Np);
      kernelInfo["defines/" "p_NblockG"] = NblockG;

      kernelInfo["defines/p_Nalign"] = USE_OCCA_MEM_BYTE_ALIGN;

/*
      //add standard boundary functions
      char* boundaryHeaderFileName;
      if (elliptic->dim == 2)
        boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary2D.h");
      else if (elliptic->dim == 3)
        boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary3D.h");
      kernelInfo["includes"] += boundaryHeaderFileName;
*/

      occa::properties AxKernelInfo = kernelInfo;
      sprintf(fileName, DELLIPTIC "/okl/ellipticAx%s.okl", suffix);
      sprintf(kernelName, "ellipticAx%s", suffix);
      if(serial) {
        AxKernelInfo["okl/enabled"] = false;
        sprintf(fileName,  DELLIPTIC "/okl/ellipticSerialAx%s.c", suffix);
      }
      elliptic->AxKernel = mesh->device.buildKernel(fileName,kernelName,AxKernelInfo);
      if(!strstr(pfloatString,dfloatString)){
        AxKernelInfo["defines/" "dfloat"] = pfloatString;
        sprintf(kernelName, "ellipticAx%s", suffix);
        elliptic->AxPfloatKernel = mesh->device.buildKernel(fileName,kernelName,AxKernelInfo);
        AxKernelInfo["defines/" "dfloat"] = dfloatString;
      }

      // check for trilinear
      if(elliptic->elementType != HEXAHEDRA) {
        sprintf(kernelName, "ellipticPartialAx%s", suffix);
      }else {
        if(elliptic->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
          sprintf(kernelName, "ellipticPartialAxTrilinear%s", suffix);
        else
          sprintf(kernelName, "ellipticPartialAx%s", suffix);
      }

      if(!serial) {
        elliptic->partialAxKernel = mesh->device.buildKernel(fileName,kernelName,AxKernelInfo);
        if(!strstr(pfloatString,dfloatString)) {
          AxKernelInfo["defines/" "dfloat"] = pfloatString;
          elliptic->partialAxPfloatKernel = mesh->device.buildKernel(fileName, kernelName, AxKernelInfo);
          AxKernelInfo["defines/" "dfloat"] = dfloatString;
        }
      }

/*
      // only for Hex3D - cubature Ax
      if(elliptic->elementType == HEXAHEDRA) {
        sprintf(fileName,  DELLIPTIC "/okl/ellipticCubatureAx%s.okl", suffix);

        sprintf(kernelName, "ellipticCubaturePartialAx%s", suffix);
        elliptic->partialCubatureAxKernel = mesh->device.buildKernel(fileName,
                                                                     kernelName,
                                                                     AxKernelInfo);
      }
*/

      if (options.compareArgs("BASIS", "BERN")) {
        sprintf(fileName, DELLIPTIC "/okl/ellipticGradientBB%s.okl", suffix);
        sprintf(kernelName, "ellipticGradientBB%s", suffix);

        elliptic->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialGradientBB%s", suffix);
        elliptic->partialGradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

/*
        sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdgBB%s.okl", suffix);
        sprintf(kernelName, "ellipticAxIpdgBB%s", suffix);
        elliptic->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialAxIpdgBB%s", suffix);
        elliptic->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
*/
      } else if (options.compareArgs("BASIS", "NODAL")) {
        sprintf(fileName, DELLIPTIC "/okl/ellipticGradient%s.okl", suffix);
        sprintf(kernelName, "ellipticGradient%s", suffix);

        elliptic->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialGradient%s", suffix);
        elliptic->partialGradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
/*
        sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdg%s.okl", suffix);
        sprintf(kernelName, "ellipticAxIpdg%s", suffix);
        elliptic->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialAxIpdg%s", suffix);
        elliptic->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
*/
      }
    }

    MPI_Barrier(mesh->comm);
  }

  MPI_Barrier(mesh->comm);
  if(mesh->rank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel); fflush(stdout);

  //new precon struct
  elliptic->precon = new precon_t();

  for (int r = 0; r < 2; r++) {
    MPI_Barrier(mesh->comm);

    if ((r == 0 && mesh->rank == 0) || (r == 1 && mesh->rank > 0)) {
      sprintf(fileName, DELLIPTIC "/okl/ellipticBlockJacobiPrecon.okl");
      sprintf(kernelName, "ellipticBlockJacobiPrecon");
      elliptic->precon->blockJacobiKernel =
        mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(kernelName, "ellipticPartialBlockJacobiPrecon");
      elliptic->precon->partialblockJacobiKernel = mesh->device.buildKernel(fileName,
                                                                            kernelName,
                                                                            kernelInfo);

      sprintf(fileName, DELLIPTIC "/okl/ellipticPatchSolver.okl");
      sprintf(kernelName, "ellipticApproxBlockJacobiSolver");
      elliptic->precon->approxBlockJacobiSolverKernel = mesh->device.buildKernel(fileName,
                                                                                 kernelName,
                                                                                 kernelInfo);

      //sizes for the coarsen and prolongation kernels. degree NFine to degree N
      int NqFine   = (Nf + 1);
      int NqCoarse = (Nc + 1);
      kernelInfo["defines/" "p_NqFine"] = Nf + 1;
      kernelInfo["defines/" "p_NqCoarse"] = Nc + 1;

      int NpFine, NpCoarse;
      switch(elliptic->elementType) {
      case TRIANGLES:
        NpFine   = (Nf + 1) * (Nf + 2) / 2;
        NpCoarse = (Nc + 1) * (Nc + 2) / 2;
        break;
      case QUADRILATERALS:
        NpFine   = (Nf + 1) * (Nf + 1);
        NpCoarse = (Nc + 1) * (Nc + 1);
        break;
      case TETRAHEDRA:
        NpFine   = (Nf + 1) * (Nf + 2) * (Nf + 3) / 6;
        NpCoarse = (Nc + 1) * (Nc + 2) * (Nc + 3) / 6;
        break;
      case HEXAHEDRA:
        NpFine   = (Nf + 1) * (Nf + 1) * (Nf + 1);
        NpCoarse = (Nc + 1) * (Nc + 1) * (Nc + 1);
        break;
      }
      kernelInfo["defines/" "p_NpFine"] = NpFine;
      kernelInfo["defines/" "p_NpCoarse"] = NpCoarse;

      int NblockVFine = maxNthreads / NpFine;
      int NblockVCoarse = maxNthreads / NpCoarse;
      kernelInfo["defines/" "p_NblockVFine"] = NblockVFine;
      kernelInfo["defines/" "p_NblockVCoarse"] = NblockVCoarse;

      // Use the same kernel with quads for the following kenels
      if(elliptic->dim == 3) {
        if(elliptic->elementType == QUADRILATERALS)
          suffix = strdup("Quad2D");
        if(elliptic->elementType == TRIANGLES)
          suffix = strdup("Tri2D");
      }

      sprintf(fileName, DELLIPTIC "/okl/ellipticPreconCoarsen%s.okl", suffix);
      sprintf(kernelName, "ellipticPreconCoarsen%s", suffix);
      elliptic->precon->coarsenKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(fileName, DELLIPTIC "/okl/ellipticPreconProlongate%s.okl", suffix);
      sprintf(kernelName, "ellipticPreconProlongate%s", suffix);
      elliptic->precon->prolongateKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }

  if(elliptic->elementType == HEXAHEDRA) {
    if(options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
      if(options.compareArgs("ELEMENT MAP", "TRILINEAR")) {
        // pack gllz, gllw, and elementwise EXYZ
        dfloat* gllzw = (dfloat*) calloc(2 * mesh->Nq, sizeof(dfloat));

        int sk = 0;
        for(int n = 0; n < mesh->Nq; ++n)
          gllzw[sk++] = mesh->gllz[n];
        for(int n = 0; n < mesh->Nq; ++n)
          gllzw[sk++] = mesh->gllw[n];

        elliptic->o_gllzw = mesh->device.malloc(2 * mesh->Nq * sizeof(dfloat), gllzw);
        free(gllzw);
      }
    }
  }

  if(!strstr(pfloatString,dfloatString)) {
    mesh->o_ggeoPfloat = mesh->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo * sizeof(pfloat));
    mesh->o_DmatricesPfloat = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(pfloat));
    mesh->o_SmatricesPfloat = mesh->device.malloc(mesh->Nq * mesh->Nq * sizeof(pfloat));
    elliptic->copyDfloatToPfloatKernel(mesh->Nelements * mesh->Np * mesh->Nggeo,
      elliptic->mesh->o_ggeoPfloat,
      mesh->o_ggeo);
    elliptic->copyDfloatToPfloatKernel(mesh->Nq * mesh->Nq,
      elliptic->mesh->o_DmatricesPfloat,
      mesh->o_Dmatrices);
    elliptic->copyDfloatToPfloatKernel(mesh->Nq * mesh->Nq,
      elliptic->mesh->o_SmatricesPfloat,
      mesh->o_Smatrices);
  }

  return elliptic;
}
