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

void ellipticOasSetup(elliptic_t* elliptic, dfloat lambda,
                      occa::properties &kernelInfo)
{
  mesh_t* mesh = elliptic->mesh;

  setupAide options = elliptic->options;

  /* STAGE 1: build overlapping extended partition problem */

  /* build one ring patch extension using a single process MPI sub-communicator
     and store in elliptic->precon->ellipticOneRing */
  ellipticBuildOneRing(elliptic, lambda, kernelInfo);

  /* STAGE 2: build coarse problem */
  //set up the base level
  int Nc = 1;
  int Nf = mesh->N;

  // build coarsener
  int NqFine   = Nf + 1;
  int NqCoarse = Nc + 1;

  int NpFine   = (Nf + 1) * (Nf + 1) * (Nf + 1);
  int NpCoarse = (Nc + 1) * (Nc + 1) * (Nc + 1);

  int NblockVFine = maxNthreads / NpFine;
  int NblockVCoarse = maxNthreads / NpCoarse;

  elliptic_t* ellipticOasCoarse;

  mesh_t* meshN1 = new mesh_t[1];

  if (mesh->N > 1) { // assume
    // TW: should create   meshCloneHex3D(mesh_t *mesh)
    printf("=============BUILDING OAS COARSE LEVEL OF DEGREE %d==================\n", Nc);
    //    ellipticOasCoarse = ellipticBuildMultigridLevel(elliptic,Nc,%m    mesh_t *mesh1 = (mesh_t*) calloc(1, sizeof(mesh_t)); // check
    meshN1->N   = 1;

    meshN1->comm = mesh->comm;
    meshN1->rank = mesh->rank;
    meshN1->size = mesh->size;

    meshN1->dim = mesh->dim;
    meshN1->Nverts = mesh->Nverts;
    meshN1->Nfaces = mesh->Nfaces;
    meshN1->NfaceVertices = mesh->NfaceVertices;
    meshN1->Nnodes = mesh->Nnodes;
    meshN1->Nfields = mesh->Nfields;

    meshN1->faceVertices = mesh->faceVertices;
    meshN1->Nelements = mesh->Nelements;
    meshN1->EX = mesh->EX;
    meshN1->EY = mesh->EY;
    meshN1->EZ = mesh->EZ;

    meshN1->NboundaryFaces = mesh->NboundaryFaces;
    meshN1->boundaryInfo = mesh->boundaryInfo;
    meshN1->EToV = mesh->EToV;

    meshParallelConnect(meshN1);

    meshConnectBoundary(meshN1);

    meshLoadReferenceNodesHex3D(meshN1, meshN1->N); // degree 1

    meshPhysicalNodesHex3D(meshN1); // rely on trilinear map for hexes

    meshGeometricFactorsHex3D(meshN1);

    meshHaloSetup(meshN1); // nada

    meshConnectFaceNodes3D(meshN1);

    meshSurfaceGeometricFactorsHex3D(meshN1);

    meshParallelConnectNodes(meshN1); // data

    setupAide optionsN1 = elliptic->options; // check this
    // optionsN1.setArgs(string("PRECONDITIONER"),    string("MULTIGRID"));
    // optionsN1.setArgs(string("PRECONDITIONER"),    string("FULLALLMOND"));
    optionsN1.setArgs(string("KRYLOV SOLVER"),        string("PCG+NONBLOCKING+FLEXIBLE"));
    optionsN1.setArgs(string("PRECONDITIONER"),       string("MULTIGRID"));
    optionsN1.setArgs(string("MAXIMUM ITERATIONS"),   string("2"));
    optionsN1.setArgs(string("FIXED ITERATION COUNT"),string("TRUE"));
    optionsN1.setArgs(string("PARALMOND CYCLE"),      string("KCYCLE"));
    optionsN1.setArgs(string("POLYNOMIAL DEGREE"),    string("1"));

    occa::properties kernelInfoN1;

    kernelInfoN1["defines"].asObject();
    kernelInfoN1["includes"].asArray();
    kernelInfoN1["header"].asArray();
    kernelInfoN1["flags"].asObject();

#if 0
    optionsN1.setArgs(string("THREAD MODEL"),    string("SERIAL"));
    meshOccaSetup3D(meshN1, optionsN1, kernelInfoN1);
#else
    meshN1->device = mesh->device; // check this
    meshN1->defaultStream = mesh->defaultStream;
    meshN1->dataStream = mesh->dataStream;
    meshN1->computeStream = mesh->computeStream;
    meshN1->device.setStream(mesh->defaultStream);
    meshOccaPopulateDevice3D(meshN1, optionsN1, kernelInfoN1);
#endif

    //    std::cout << "KINFO LOOK AT THIS: " << kernelInfoN1 << std::endl;
    //    std::cout << "OPTIO LOOK AT THIS: " << optionsN1 << std::endl;

    // set up
    ellipticOasCoarse = ellipticSetup(meshN1, lambda, kernelInfoN1, optionsN1);
  }else{
    ellipticOasCoarse = elliptic;
  }

  dfloat* P    = (dfloat*) calloc(NqFine * NqCoarse,sizeof(dfloat));
  dfloat* R    = (dfloat*) calloc(NqFine * NqCoarse,sizeof(dfloat));

  // hard wire for linears
  for(int n = 0; n < NqFine; ++n) {
    P[n * NqCoarse + 0] = 0.5 * (1 - mesh->gllz[n]);
    P[n * NqCoarse + 1] = 0.5 * (1 + mesh->gllz[n]);
    R[0 * NqFine + n] = 0.5 * (1 - mesh->gllz[n]);
    R[1 * NqFine + n] = 0.5 * (1 + mesh->gllz[n]);
  }

  occa::memory o_R = elliptic->mesh->device.malloc(NqFine * NqCoarse * sizeof(dfloat), R);
  occa::memory o_P = elliptic->mesh->device.malloc(NqFine * NqCoarse * sizeof(dfloat), P);

  free(P);
  free(R);

  elliptic->precon->ellipticOasCoarse = ellipticOasCoarse;
  elliptic->precon->o_oasRestrictionMatrix = o_R;
  elliptic->precon->o_oasProlongationMatrix = o_P;

  // build degree 1 coarsening and prolongation matrices and kernels

  kernelInfo["defines/" "p_NqFine"] = Nf + 1;
  kernelInfo["defines/" "p_NqCoarse"] = Nc + 1;

  kernelInfo["defines/" "p_NpFine"] = NpFine;
  kernelInfo["defines/" "p_NpCoarse"] = NpCoarse;

  kernelInfo["defines/" "p_NblockVFine"] = NblockVFine;
  kernelInfo["defines/" "p_NblockVCoarse"] = NblockVCoarse;

  char* suffix;

  if(elliptic->elementType == HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  for (int r = 0; r < 2; r++) {
    if ((r == 0 && mesh->rank == 0) || (r == 1 && mesh->rank > 0)) {
      sprintf(fileName, DELLIPTIC "/okl/ellipticPreconCoarsen%s.okl", suffix);
      sprintf(kernelName, "ellipticPreconCoarsen%s", suffix);
      elliptic->precon->oasRestrictionKernel =
        mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(fileName, DELLIPTIC "/okl/ellipticPreconProlongate%s.okl", suffix);
      sprintf(kernelName, "ellipticPreconProlongate%s", suffix);
      elliptic->precon->oasProlongationKernel = mesh->device.buildKernel(fileName,
                                                                         kernelName,
                                                                         kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }
}
