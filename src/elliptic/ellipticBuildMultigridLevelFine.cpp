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
#include "platform.hpp"

elliptic_t* ellipticBuildMultigridLevelFine(elliptic_t* baseElliptic)
{
  platform_t* platform = platform_t::getInstance();
  elliptic_t* elliptic = new elliptic_t();
  memcpy(elliptic, baseElliptic, sizeof(*baseElliptic));

  const int serial = baseElliptic->options.compareArgs("THREAD MODEL", "SERIAL");

  elliptic->var_coeff = 0;
  elliptic->lambda = (dfloat*) calloc(elliptic->Nfields, sizeof(dfloat)); // enforce lambda = 0

  mesh_t* mesh = elliptic->mesh;

  if(!strstr(pfloatString,dfloatString)) {
    mesh->o_ggeoPfloat = platform->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo * sizeof(pfloat));
    mesh->o_DmatricesPfloat = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(pfloat));
    mesh->o_SmatricesPfloat = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(pfloat));

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

  string suffix;
  occa::properties kernelInfo = ellipticKernelInfo(mesh);
  if(elliptic->elementType == HEXAHEDRA)
    suffix = "Hex3D";

  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;

  // add custom defines
  kernelInfo["defines/" "p_NpP"] = (mesh->Np + mesh->Nfp * mesh->Nfaces);
  kernelInfo["defines/" "p_Nverts"] = mesh->Nverts;

  int Nmax = mymax(mesh->Np, mesh->Nfaces * mesh->Nfp);
  kernelInfo["defines/" "p_Nmax"] = Nmax;

  int maxNodes = mymax(mesh->Np, (mesh->Nfp * mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"] = maxNodes;

  int NblockV = mymax(1,BLOCKSIZE / mesh->Np);
  kernelInfo["defines/" "p_NblockV"] = NblockV;

  int one = 1; //set to one for now. TODO: try optimizing over these
  kernelInfo["defines/" "p_NnodesV"] = one;

  int NblockS = mymax(1,BLOCKSIZE / maxNodes);
  kernelInfo["defines/" "p_NblockS"] = NblockS;

  int NblockP = mymax(1,BLOCKSIZE / (4 * mesh->Np)); // get close to BLOCKSIZE threads
  kernelInfo["defines/" "p_NblockP"] = NblockP;

  int NblockG;
  if(mesh->Np <= 32) NblockG = ( 32 / mesh->Np );
  else NblockG = mymax(1,BLOCKSIZE / mesh->Np);
  kernelInfo["defines/" "p_NblockG"] = NblockG;

  kernelInfo["defines/" "p_eNfields"] = elliptic->Nfields;
  kernelInfo["defines/p_Nalign"] = USE_OCCA_MEM_BYTE_ALIGN;

  string filename, kernelName;

  string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const string oklpath = install_dir + "/okl/elliptic/";

  {
      occa::properties AxKernelInfo = kernelInfo;

      filename = oklpath + "ellipticAx" + suffix + ".okl";
      kernelName = "ellipticAx" + suffix;
      if(serial) {
        AxKernelInfo["okl/enabled"] = false;
        filename = oklpath + "ellipticSerialAx" + suffix + ".c";
      }
      elliptic->AxKernel = platform->device.buildKernel(filename.c_str(),kernelName.c_str(),AxKernelInfo);

      if(!strstr(pfloatString,dfloatString)) {
        AxKernelInfo["defines/" "dfloat"] = pfloatString;
        kernelName = "ellipticAx" + suffix;
        elliptic->AxPfloatKernel = platform->device.buildKernel(filename.c_str(),kernelName.c_str(),AxKernelInfo);
        AxKernelInfo["defines/" "dfloat"] = dfloatString;
      }

      if(elliptic->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
        kernelName = "ellipticPartialAxTrilinear" + suffix;
      else
        kernelName = "ellipticPartialAx" + suffix;

      if(!serial) {
        elliptic->partialAxKernel = platform->device.buildKernel(filename.c_str(),kernelName.c_str(),AxKernelInfo);
        if(!strstr(pfloatString,dfloatString)) {
          AxKernelInfo["defines/" "dfloat"] = pfloatString;
          elliptic->partialAxPfloatKernel =
            platform->device.buildKernel(filename.c_str(), kernelName.c_str(), AxKernelInfo);
          AxKernelInfo["defines/" "dfloat"] = dfloatString;
        }
      }
  }

  return elliptic;
}
