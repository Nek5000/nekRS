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
  memcpy(mesh,baseElliptic->mesh,sizeof(mesh_t));

  elliptic->mesh = mesh;

  setupAide options = elliptic->options;

  switch(elliptic->elementType) {
  case HEXAHEDRA:
    meshLoadReferenceNodesHex3D(mesh, Nc, 1);
    meshHaloSetup(mesh);
    meshPhysicalNodesHex3D(mesh, buildOnly);
    meshHaloPhysicalNodes(mesh);
    meshGeometricFactorsHex3D(mesh);

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
  meshParallelConnectNodes(mesh, buildOnly);

  //dont need these once vmap is made
  free(mesh->x);
  free(mesh->y);
  if (elliptic->dim == 3)
    free(mesh->z);

  dlong Ntotal = mesh->Np * mesh->Nelements;
  dlong Nblock = mymax(1,(Ntotal + BLOCKSIZE - 1) / BLOCKSIZE);
  dlong Nblock2 = mymax(1,(Nblock + BLOCKSIZE - 1) / BLOCKSIZE);

  elliptic->Nblock = Nblock;
  elliptic->Nblock2 = Nblock2;

  if (elliptic->elementType == HEXAHEDRA) {
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

    mesh->o_vmapM =
      mesh->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(dlong),
                          mesh->vmapM);

    mesh->o_vmapP =
      mesh->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(dlong),
                          mesh->vmapP);

    mesh->LIFT = baseElliptic->mesh->LIFT; //dummy buffer
    mesh->o_LIFTT = baseElliptic->mesh->o_LIFTT; //dummy buffer
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
  string suffix;
  if(elliptic->elementType == HEXAHEDRA)
    suffix = "Hex3D";

  string filename, kernelName;

  MPI_Barrier(mesh->comm);
  double tStartLoadKernel = MPI_Wtime();
  if(mesh->rank == 0) printf("loading elliptic MG kernels ... ");
  fflush(stdout);

  string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const string oklpath = install_dir + "/okl/elliptic/";

  for (int r = 0; r < 2; r++) {
    MPI_Barrier(mesh->comm);

    if ((r == 0 && mesh->rank == 0) || (r == 1 && mesh->rank > 0)) {
      kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;

      // add custom defines
      kernelInfo["defines/" "p_NpP"] = (mesh->Np + mesh->Nfp * mesh->Nfaces);
      kernelInfo["defines/" "p_Nverts"] = mesh->Nverts;

      int Nmax = mymax(mesh->Np, mesh->Nfaces * mesh->Nfp);
      kernelInfo["defines/" "p_Nmax"] = Nmax;

      int maxNodes = mymax(mesh->Np, (mesh->Nfp * mesh->Nfaces));
      kernelInfo["defines/" "p_maxNodes"] = maxNodes;

      int NblockV = mymax(1,BLOCKSIZE / mesh->Np); // works for CUDA
      kernelInfo["defines/" "p_NblockV"] = NblockV;

      int one = 1; //set to one for now. TODO: try optimizing over these
      kernelInfo["defines/" "p_NnodesV"] = one;

      int NblockS = mymax(1,BLOCKSIZE / maxNodes); // works for CUDA
      kernelInfo["defines/" "p_NblockS"] = NblockS;

      int NblockP = mymax(1,BLOCKSIZE / (4 * mesh->Np)); // get close to BLOCKSIZE threads
      kernelInfo["defines/" "p_NblockP"] = NblockP;

      int NblockG;
      if(mesh->Np <= 32) NblockG = ( 32 / mesh->Np );
      else NblockG = mymax(1,BLOCKSIZE / mesh->Np);
      kernelInfo["defines/" "p_NblockG"] = NblockG;

      kernelInfo["defines/p_Nalign"] = USE_OCCA_MEM_BYTE_ALIGN;

/*
      //add standard boundary functions
      char* boundaryHeaderFileName;
      if (elliptic->dim == 2)
        boundaryHeaderFileName = strdup(oklpath + "/data/ellipticBoundary2D.h");
      else if (elliptic->dim == 3)
        boundaryHeaderFileName = strdup(oklpath + "/data/ellipticBoundary3D.h");
      kernelInfo["includes"] += boundaryHeaderFileName;
 */

      occa::properties AxKernelInfo = kernelInfo;
      filename = oklpath + "ellipticAx" + suffix + ".okl";
      kernelName = "ellipticAx" + suffix;
      if(serial) {
        AxKernelInfo["okl/enabled"] = false;
        filename = oklpath + "ellipticSerialAx" + suffix + ".c";
      }
      elliptic->AxKernel = mesh->device.buildKernel(filename.c_str(),kernelName.c_str(),AxKernelInfo);
      if(!strstr(pfloatString,dfloatString)) {
        AxKernelInfo["defines/" "dfloat"] = pfloatString;
        kernelName = "ellipticAx" + suffix;
        elliptic->AxPfloatKernel = mesh->device.buildKernel(filename.c_str(),kernelName.c_str(),AxKernelInfo);
        AxKernelInfo["defines/" "dfloat"] = dfloatString;
      }

      // check for trilinear
      if(elliptic->elementType != HEXAHEDRA) {
        kernelName = "ellipticPartialAx" + suffix;
      }else {
        if(elliptic->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
          kernelName = "ellipticPartialAxTrilinear" + suffix;
        else
          kernelName = "ellipticPartialAx" + suffix;
      }

      if(!serial) {
        elliptic->partialAxKernel = mesh->device.buildKernel(filename.c_str(),kernelName.c_str(),AxKernelInfo);
        if(!strstr(pfloatString,dfloatString)) {
          AxKernelInfo["defines/" "dfloat"] = pfloatString;
          elliptic->partialAxPfloatKernel =
            mesh->device.buildKernel(filename.c_str(), kernelName.c_str(), AxKernelInfo);
          AxKernelInfo["defines/" "dfloat"] = dfloatString;
        }
      }

      if (options.compareArgs("BASIS", "NODAL")) {
        filename = oklpath + "ellipticGradient" + suffix + ".okl";
        kernelName = "ellipticGradient" + suffix;

        elliptic->gradientKernel = mesh->device.buildKernel(filename.c_str(),kernelName.c_str(),kernelInfo);

        kernelName = "ellipticPartialGradient" + suffix;
        elliptic->partialGradientKernel = mesh->device.buildKernel(filename.c_str(),kernelName.c_str(),kernelInfo);
/*
        sprintf(fileName, oklpath + "ellipticAxIpdg%s.okl", suffix);
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
  if(mesh->rank == 0) printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel);
  fflush(stdout);

  //new precon struct
  elliptic->precon = new precon_t();

  for (int r = 0; r < 2; r++) {
    MPI_Barrier(mesh->comm);

    if ((r == 0 && mesh->rank == 0) || (r == 1 && mesh->rank > 0)) {
      filename = oklpath + "ellipticBlockJacobiPrecon.okl";
      kernelName = "ellipticBlockJacobiPrecon";
      elliptic->precon->blockJacobiKernel =
        mesh->device.buildKernel(filename.c_str(),kernelName.c_str(),kernelInfo);

      kernelName = "ellipticPartialBlockJacobiPrecon";
      elliptic->precon->partialblockJacobiKernel = mesh->device.buildKernel(filename.c_str(),
                                                                            kernelName.c_str(),
                                                                            kernelInfo);

      filename = oklpath + "ellipticPatchSolver.okl";
      kernelName = "ellipticApproxBlockJacobiSolver";
      elliptic->precon->approxBlockJacobiSolverKernel = mesh->device.buildKernel(filename.c_str(),
                                                                                 kernelName.c_str(),
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

      int NblockVFine = BLOCKSIZE / NpFine;
      int NblockVCoarse = BLOCKSIZE / NpCoarse;
      kernelInfo["defines/" "p_NblockVFine"] = NblockVFine;
      kernelInfo["defines/" "p_NblockVCoarse"] = NblockVCoarse;

      // Use the same kernel with quads for the following kenels
      if(elliptic->dim == 3) {
        if(elliptic->elementType == QUADRILATERALS)
          suffix = strdup("Quad2D");
        if(elliptic->elementType == TRIANGLES)
          suffix = strdup("Tri2D");
      }

      filename = oklpath + "ellipticPreconCoarsen" + suffix + ".okl";
      kernelName = "ellipticPreconCoarsen" + suffix;
      elliptic->precon->coarsenKernel = mesh->device.buildKernel(filename.c_str(),kernelName.c_str(),kernelInfo);

      filename = oklpath + "ellipticPreconProlongate" + suffix + ".okl";
      kernelName = "ellipticPreconProlongate" + suffix;
      elliptic->precon->prolongateKernel = mesh->device.buildKernel(filename.c_str(),kernelName.c_str(),kernelInfo);
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
