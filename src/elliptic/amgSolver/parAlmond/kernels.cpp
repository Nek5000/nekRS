/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"
#include "platform.hpp"

namespace parAlmond {

int Nrefs = 0;

occa::kernel haloExtractKernel;

occa::kernel SpMVcsrKernel1;
occa::kernel SpMVcsrKernel2;
occa::kernel SpMVellKernel1;
occa::kernel SpMVellKernel2;
occa::kernel SpMVmcsrKernel1;
occa::kernel SpMVmcsrKernel2;

occa::kernel vectorSetKernel;
occa::kernel vectorScaleKernel;
occa::kernel vectorAddScalarKernel;
occa::kernel vectorAddKernel1;
occa::kernel vectorAddKernel2;
occa::kernel vectorDotStarKernel1;
occa::kernel vectorDotStarKernel2;
occa::kernel vectorInnerProdKernel;
occa::kernel kcycleCombinedOp1Kernel;
occa::kernel kcycleCombinedOp2Kernel;
occa::kernel kcycleWeightedCombinedOp1Kernel;
occa::kernel kcycleWeightedCombinedOp2Kernel;
occa::kernel vectorAddInnerProdKernel;
occa::kernel vectorAddWeightedInnerProdKernel;

void buildParAlmondKernels(MPI_Comm comm, occa::device device){
  // TODO: cleanup
}

void freeParAlmondKernels() {
  vectorDotStarKernel1.free();
  vectorDotStarKernel2.free();
}


} //namespace parAlmond
