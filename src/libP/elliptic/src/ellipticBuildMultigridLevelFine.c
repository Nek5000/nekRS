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

elliptic_t* ellipticBuildMultigridLevelFine(elliptic_t* baseElliptic)
{
  elliptic_t* elliptic = new elliptic_t();
  memcpy(elliptic, baseElliptic, sizeof(*baseElliptic));

  const int serial = baseElliptic->options.compareArgs("THREAD MODEL", "SERIAL");

  elliptic->var_coeff = 0;
  elliptic->lambda = (dfloat*) calloc(elliptic->Nfields, sizeof(dfloat)); // enforce lambda = 0

  mesh_t* mesh = elliptic->mesh;

  if(!strstr(pfloatString,dfloatString)) {
    mesh->o_ggeoPfloat = mesh->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo * sizeof(pfloat));
    mesh->o_DmatricesPfloat = mesh->device.malloc(mesh->Nq * mesh->Nq*sizeof(pfloat));
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

  char* suffix;
  occa::properties kernelInfo = ellipticKernelInfo(mesh);

  if(elliptic->elementType == TRIANGLES) {
    if(elliptic->dim == 2)
      suffix = strdup("Tri2D");
    else
      suffix = strdup("Tri3D");
  }
  if(elliptic->elementType == QUADRILATERALS) {
    if(elliptic->dim == 2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D");
  }
  if(elliptic->elementType == TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(elliptic->elementType == HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  for (int r = 0; r < 2; r++) {
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

      kernelInfo["defines/" "p_eNfields"] = elliptic->Nfields;
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

      if(elliptic->elementType != HEXAHEDRA) {
        sprintf(kernelName, "ellipticPartialAx%s", suffix);
      }else{
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
    }

    MPI_Barrier(mesh->comm);
  }

  return elliptic;
}
