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

#ifndef PARALMOND_KERNELS_HPP
#define PARALMOND_KERNELS_HPP

namespace parAlmond {

  void buildParAlmondKernels(MPI_Comm comm, occa::device device);

  void freeParAlmondKernels();

  extern int Nrefs;

  extern occa::kernel haloExtractKernel;

  extern occa::kernel SpMVcsrKernel1;
  extern occa::kernel SpMVcsrKernel2;
  extern occa::kernel SpMVellKernel1;
  extern occa::kernel SpMVellKernel2;
  extern occa::kernel SpMVmcsrKernel1;
  extern occa::kernel SpMVmcsrKernel2;

  extern occa::kernel vectorSetKernel;
  extern occa::kernel vectorScaleKernel;
  extern occa::kernel vectorAddScalarKernel;
  extern occa::kernel vectorAddKernel1;
  extern occa::kernel vectorAddKernel2;
  extern occa::kernel vectorDotStarKernel1;
  extern occa::kernel vectorDotStarKernel2;
  extern occa::kernel vectorInnerProdKernel;
  extern occa::kernel vectorAddInnerProdKernel;
  extern occa::kernel vectorAddWeightedInnerProdKernel;
  extern occa::kernel kcycleCombinedOp1Kernel;
  extern occa::kernel kcycleCombinedOp2Kernel;
  extern occa::kernel kcycleWeightedCombinedOp1Kernel;
  extern occa::kernel kcycleWeightedCombinedOp2Kernel;

} //namespace parAlmond

#endif