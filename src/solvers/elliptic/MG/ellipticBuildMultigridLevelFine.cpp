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
  
  auto elliptic = new elliptic_t();
  memcpy(elliptic, baseElliptic, sizeof(*baseElliptic));

  auto mesh = new mesh_t();
  memcpy(mesh, baseElliptic->mesh, sizeof(*baseElliptic->mesh));
  elliptic->mesh = mesh;

  elliptic->mgLevel = true;

  ellipticBuildPreconditionerKernels(elliptic);

  elliptic->o_lambda0 = platform->device.malloc<pfloat>(mesh->Nlocal);
  if(elliptic->options.compareArgs("ELLIPTIC PRECO COEFF FIELD", "TRUE")) {
    platform->copyDfloatToPfloatKernel(mesh->Nlocal, baseElliptic->o_lambda0, elliptic->o_lambda0);
  } else {
    platform->linAlg->pfill(mesh->Nlocal, elliptic->lambda0Avg, elliptic->o_lambda0);
  }

  if(baseElliptic->poisson) { 
    elliptic->o_lambda1 = nullptr; 
  } else {
    elliptic->o_lambda1 = platform->device.malloc<pfloat>(mesh->Nlocal);
    platform->copyDfloatToPfloatKernel(mesh->Nlocal, baseElliptic->o_lambda1, elliptic->o_lambda1);
  }


  pfloat *tmp = (pfloat*) calloc(mesh->Nlocal, sizeof(pfloat));
  for(int i = 0; i < mesh->Nlocal; i++) {
     tmp[i] = (pfloat) baseElliptic->ogs->invDegree[i];
  }
  elliptic->o_invDegree = platform->device.malloc<pfloat>(mesh->Nlocal, tmp);
  free(tmp);

  if(!std::is_same<pfloat, dfloat>::value) {
    mesh->o_ggeo = platform->device.malloc<pfloat>(mesh->Nelements * mesh->Np * mesh->Nggeo);
    mesh->o_D = platform->device.malloc<pfloat>(mesh->Nq * mesh->Nq);
    mesh->o_DT = platform->device.malloc<pfloat>(mesh->Nq * mesh->Nq);

    platform->copyDfloatToPfloatKernel(mesh->Nelements * mesh->Np * mesh->Nggeo,
                                       baseElliptic->mesh->o_ggeo,
                                       mesh->o_ggeo);
    platform->copyDfloatToPfloatKernel(mesh->Nq * mesh->Nq,
                                       baseElliptic->mesh->o_D,
                                       mesh->o_D);
    platform->copyDfloatToPfloatKernel(mesh->Nq * mesh->Nq,
                                       baseElliptic->mesh->o_DT,
                                       mesh->o_DT);
  }

  std::string suffix = "CoeffHex3D";

  std::string kernelName;

  const std::string poissonPrefix = elliptic->poisson ? "poisson-" : "";

  if(elliptic->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
    kernelName = "ellipticPartialAxTrilinear" + suffix;
  else
    kernelName = "ellipticPartialAx" + suffix;

  const std::string kernelSuffix = gen_suffix(elliptic, pfloatString);
  elliptic->AxKernel = platform->kernels.get(poissonPrefix + kernelName + kernelSuffix);

  return elliptic;
}
