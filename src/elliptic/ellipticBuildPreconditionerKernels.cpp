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
#include <string>
#include "platform.hpp"
#include "linAlg.hpp"

void ellipticBuildPreconditionerKernels(elliptic_t* elliptic)
{
  
  mesh_t* mesh      = elliptic->mesh;

  std::string prefix = "Hex3D";
  std::string kernelName;

  const std::string orderSuffix = std::string("_") + std::to_string(mesh->N);

  {
    kernelName = "mask";
    mesh->maskKernel =
      platform->kernels.get(kernelName + orderSuffix);

    mesh->maskPfloatKernel =
      platform->kernels.get(kernelName + orderSuffix + "pfloat");

    kernelName = "updateSmoothedSolutionVec";
    elliptic->updateSmoothedSolutionVecKernel =
      platform->kernels.get(kernelName + orderSuffix);

    kernelName = "updateChebyshevSolutionVec";
    elliptic->updateChebyshevSolutionVecKernel =
      platform->kernels.get(kernelName + orderSuffix);

    kernelName = "updateIntermediateSolutionVec";
    elliptic->updateIntermediateSolutionVecKernel =
      platform->kernels.get(kernelName + orderSuffix);

    kernelName = "ellipticBlockBuildDiagonalHex3D";
    const std::string poissonPrefix = elliptic->poisson ? "poisson-" : "";
    elliptic->ellipticBlockBuildDiagonalKernel =
        platform->kernels.get(poissonPrefix + kernelName + orderSuffix);

    kernelName = "ellipticBlockBuildDiagonalPfloatHex3D";
    elliptic->ellipticBlockBuildDiagonalPfloatKernel =
        platform->kernels.get(poissonPrefix + kernelName + orderSuffix);

    elliptic->axmyzManyPfloatKernel = platform->kernels.get("axmyzManyPfloat");
  }
}
