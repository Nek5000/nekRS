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

occa::properties ellipticKernelInfo(mesh_t* mesh)
{

  
  // info for kernel construction
  occa::properties kernelInfo = platform->kernelInfo;

  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  kernelInfo["defines/" "p_dim"] = mesh->dim;
  kernelInfo["defines/" "p_Nfields"] = mesh->Nfields;
  kernelInfo["defines/" "p_N"] = mesh->N;
  kernelInfo["defines/" "p_Nq"] = mesh->N + 1;
  kernelInfo["defines/" "p_Np"] = mesh->Np;
  kernelInfo["defines/" "p_Nfp"] = mesh->Nfp;
  kernelInfo["defines/" "p_Nfaces"] = mesh->Nfaces;
  kernelInfo["defines/" "p_NfacesNfp"] = mesh->Nfp * mesh->Nfaces;
  kernelInfo["defines/" "p_Nvgeo"] = mesh->Nvgeo;
  kernelInfo["defines/" "p_Nsgeo"] = mesh->Nsgeo;
  kernelInfo["defines/" "p_Nggeo"] = mesh->Nggeo;

  kernelInfo["defines/" "p_halfN"] = (int)((mesh->Nq + 1) / 2);

  kernelInfo["defines/" "p_NXID"] = NXID;
  kernelInfo["defines/" "p_NYID"] = NYID;
  kernelInfo["defines/" "p_NZID"] = NZID;
  kernelInfo["defines/" "p_SJID"] = SJID;
  kernelInfo["defines/" "p_IJID"] = IJID;
  kernelInfo["defines/" "p_WSJID"] = WSJID;
  kernelInfo["defines/" "p_IHID"] = IHID;

  kernelInfo["defines/" "p_G00ID"] = G00ID;
  kernelInfo["defines/" "p_G01ID"] = G01ID;
  kernelInfo["defines/" "p_G02ID"] = G02ID;
  kernelInfo["defines/" "p_G11ID"] = G11ID;
  kernelInfo["defines/" "p_G12ID"] = G12ID;
  kernelInfo["defines/" "p_G22ID"] = G22ID;
  kernelInfo["defines/" "p_GWJID"] = GWJID;

  kernelInfo["defines/" "p_RXID"] = RXID;
  kernelInfo["defines/" "p_SXID"] = SXID;
  kernelInfo["defines/" "p_TXID"] = TXID;

  kernelInfo["defines/" "p_RYID"] = RYID;
  kernelInfo["defines/" "p_SYID"] = SYID;
  kernelInfo["defines/" "p_TYID"] = TYID;

  kernelInfo["defines/" "p_RZID"] = RZID;
  kernelInfo["defines/" "p_SZID"] = SZID;
  kernelInfo["defines/" "p_TZID"] = TZID;

  kernelInfo["defines/" "p_JID"] = JID;
  kernelInfo["defines/" "p_JWID"] = JWID;

  return kernelInfo;
}
