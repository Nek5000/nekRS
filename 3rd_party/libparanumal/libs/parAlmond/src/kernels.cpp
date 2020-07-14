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

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  double seed = (double) rank;
  srand48(seed);

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  if(sizeof(dlong)==4){
    kernelInfo["defines/" "dlong"]="int";
  }
  if(sizeof(dlong)==8){
    kernelInfo["defines/" "dlong"]="long long int";
  }

  if(sizeof(dfloat) == sizeof(double)){
    kernelInfo["defines/" "dfloat"]= "double";
    kernelInfo["defines/" "dfloat4"]= "double4";
  }
  else if(sizeof(dfloat) == sizeof(float)){
    kernelInfo["defines/" "dfloat"]= "float";
    kernelInfo["defines/" "dfloat4"]= "float4";
  }

  kernelInfo["defines/" "p_BLOCKSIZE"]= BLOCKSIZE;

  if(device.mode()=="OpenCL"){
    //kernelInfo["compiler_flags"] += "-cl-opt-disable";
  }

  if(device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo["compiler_flags"] += " --ftz=true ";
    kernelInfo["compiler_flags"] += " --prec-div=false ";
    kernelInfo["compiler_flags"] += " --prec-sqrt=false ";
    kernelInfo["compiler_flags"] += " --use_fast_math ";
    kernelInfo["compiler_flags"] += " --fmad=true "; // compiler option for cuda
  }

  if (rank==0) printf("Compiling parALMOND Kernels...");fflush(stdout);

  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {      

      SpMVcsrKernel1  = device.buildKernel(DPARALMOND"/okl/SpMVcsr.okl",  "SpMVcsr1",  kernelInfo);
      SpMVcsrKernel2  = device.buildKernel(DPARALMOND"/okl/SpMVcsr.okl",  "SpMVcsr2",  kernelInfo);
      SpMVellKernel1  = device.buildKernel(DPARALMOND"/okl/SpMVell.okl",  "SpMVell1",  kernelInfo);
      SpMVellKernel2  = device.buildKernel(DPARALMOND"/okl/SpMVell.okl",  "SpMVell2",  kernelInfo);
      SpMVmcsrKernel1 = device.buildKernel(DPARALMOND"/okl/SpMVmcsr.okl", "SpMVmcsr1", kernelInfo);
      SpMVmcsrKernel2 = device.buildKernel(DPARALMOND"/okl/SpMVmcsr.okl", "SpMVmcsr2", kernelInfo);

      vectorSetKernel = device.buildKernel(DPARALMOND"/okl/vectorSet.okl", "vectorSet", kernelInfo);
      vectorScaleKernel = device.buildKernel(DPARALMOND"/okl/vectorScale.okl", "vectorScale", kernelInfo);
      vectorAddScalarKernel = device.buildKernel(DPARALMOND"/okl/vectorAddScalar.okl", "vectorAddScalar", kernelInfo);
      vectorAddKernel1 = device.buildKernel(DPARALMOND"/okl/vectorAdd.okl", "vectorAdd1", kernelInfo);
      vectorAddKernel2 = device.buildKernel(DPARALMOND"/okl/vectorAdd.okl", "vectorAdd2", kernelInfo);
      vectorDotStarKernel1 = device.buildKernel(DPARALMOND"/okl/vectorDotStar.okl", "vectorDotStar1", kernelInfo);
      vectorDotStarKernel2 = device.buildKernel(DPARALMOND"/okl/vectorDotStar.okl", "vectorDotStar2", kernelInfo);
      vectorInnerProdKernel = device.buildKernel(DPARALMOND"/okl/vectorInnerProd.okl", "vectorInnerProd", kernelInfo);

      vectorAddInnerProdKernel = device.buildKernel(DPARALMOND"/okl/vectorAddInnerProd.okl", "vectorAddInnerProd", kernelInfo);
      vectorAddWeightedInnerProdKernel = device.buildKernel(DPARALMOND"/okl/vectorAddInnerProd.okl", "vectorAddWeightedInnerProd", kernelInfo);

      kcycleCombinedOp1Kernel = device.buildKernel(DPARALMOND"/okl/kcycleCombinedOp.okl", "kcycleCombinedOp1", kernelInfo);
      kcycleCombinedOp2Kernel = device.buildKernel(DPARALMOND"/okl/kcycleCombinedOp.okl", "kcycleCombinedOp2", kernelInfo);
      kcycleWeightedCombinedOp1Kernel = device.buildKernel(DPARALMOND"/okl/kcycleCombinedOp.okl", "kcycleWeightedCombinedOp1", kernelInfo);
      kcycleWeightedCombinedOp2Kernel = device.buildKernel(DPARALMOND"/okl/kcycleCombinedOp.okl", "kcycleWeightedCombinedOp2", kernelInfo);

      haloExtractKernel = device.buildKernel(DPARALMOND"/okl/haloExtract.okl", "haloExtract", kernelInfo);
    }
    MPI_Barrier(comm);
  }
  if(rank==0) printf("done.\n");
}

void freeParAlmondKernels() {

  haloExtractKernel.free();

  SpMVcsrKernel1.free();
  SpMVcsrKernel2.free();
  SpMVellKernel1.free();
  SpMVellKernel2.free();
  SpMVmcsrKernel1.free();
  SpMVmcsrKernel2.free();

  vectorSetKernel.free();
  vectorScaleKernel.free();
  vectorAddScalarKernel.free();
  vectorAddKernel1.free();
  vectorAddKernel2.free();
  vectorDotStarKernel1.free();
  vectorDotStarKernel2.free();
  vectorInnerProdKernel.free();
  kcycleCombinedOp1Kernel.free();
  kcycleCombinedOp2Kernel.free();
  kcycleWeightedCombinedOp1Kernel.free();
  kcycleWeightedCombinedOp2Kernel.free();
  vectorAddInnerProdKernel.free();
  vectorAddWeightedInnerProdKernel.free();

}


} //namespace parAlmond
