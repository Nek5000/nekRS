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

namespace{
std::string gen_suffix(const elliptic_t * elliptic, const char * floatString)
{
  const std::string precision = std::string(floatString);
  if(precision.find(pfloatString) != std::string::npos){
    return std::string("_") + std::to_string(elliptic->mesh->N) + std::string("pfloat");
  }
  else{
    return std::string("_") + std::to_string(elliptic->mesh->N);
  }
  
}
}

elliptic_t* ellipticBuildMultigridLevelFine(elliptic_t* baseElliptic)
{
  
  elliptic_t* elliptic = new elliptic_t();
  memcpy(elliptic, baseElliptic, sizeof(*baseElliptic));

  elliptic->mgLevel = true;

  mesh_t* mesh = elliptic->mesh;
  ellipticBuildPreconditionerKernels(elliptic);

  elliptic->coeffField = baseElliptic->coeffField;
  elliptic->coeffFieldPreco = baseElliptic->coeffFieldPreco;

  elliptic->o_lambda = platform->device.malloc(mesh->Nlocal, sizeof(pfloat));
  platform->copyDfloatToPfloatKernel(mesh->Nlocal, baseElliptic->o_lambda, elliptic->o_lambda);

  pfloat *tmp = (pfloat*) calloc(mesh->Nlocal, sizeof(pfloat));
  for(int i = 0; i < mesh->Nlocal; i++) {
     tmp[i] = (pfloat) baseElliptic->ogs->invDegree[i];
  }
  elliptic->o_invDegree = platform->device.malloc(mesh->Nlocal * sizeof(pfloat), tmp);
  free(tmp);

  if(!strstr(pfloatString,dfloatString)) {
    mesh->o_ggeoPfloat = platform->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo ,  sizeof(pfloat));
    mesh->o_DPfloat = platform->device.malloc(mesh->Nq * mesh->Nq ,  sizeof(pfloat));
    mesh->o_DTPfloat = platform->device.malloc(mesh->Nq * mesh->Nq ,  sizeof(pfloat));

    platform->copyDfloatToPfloatKernel(mesh->Nelements * mesh->Np * mesh->Nggeo,
                                       mesh->o_ggeo,
                                       elliptic->mesh->o_ggeoPfloat);
    platform->copyDfloatToPfloatKernel(mesh->Nq * mesh->Nq,
                                       mesh->o_D,
                                       elliptic->mesh->o_DPfloat);
    platform->copyDfloatToPfloatKernel(mesh->Nq * mesh->Nq,
                                       mesh->o_DT,
                                       elliptic->mesh->o_DTPfloat);
  }

  std::string suffix = elliptic->coeffFieldPreco ? "CoeffHex3D" : "Hex3D";

  std::string kernelName;

  const std::string poissonPrefix = elliptic->poisson ? "poisson-" : "";

  {
      if(elliptic->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
        kernelName = "ellipticPartialAxTrilinear" + suffix;
      else
        kernelName = "ellipticPartialAx" + suffix;

      {
        const std::string kernelSuffix = gen_suffix(elliptic, dfloatString);
        elliptic->AxKernel = platform->kernels.get(poissonPrefix + kernelName + kernelSuffix);
      }
      {
        const std::string kernelSuffix = gen_suffix(elliptic, pfloatString);
        elliptic->AxPfloatKernel =
          platform->kernels.get(poissonPrefix + kernelName + kernelSuffix);
      }
  }

  return elliptic;
}
