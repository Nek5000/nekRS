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

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  double seed = (double) rank;
  srand48(seed);

  
  occa::properties kernelInfo = platform->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));

  MPI_Barrier(comm);
  const double tStart = MPI_Wtime();
  if (rank==0) printf("loading parALMOND kernels ... ");fflush(stdout);

  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {      
      const string oklpath = install_dir + "/okl/parAlmond/";
      string filename;

      filename = oklpath + "vectorDotStar.okl";
      vectorDotStarKernel1 = device.buildKernel(filename.c_str(), "vectorDotStar1", kernelInfo);
      vectorDotStarKernel2 = device.buildKernel(filename.c_str(), "vectorDotStar2", kernelInfo);

      /*
      filename = oklpath + "SpMVcsr.okl";
      SpMVcsrKernel1  = device.buildKernel(filename.c_str(),  "SpMVcsr1",  kernelInfo);
      SpMVcsrKernel2  = device.buildKernel(filename.c_str(),  "SpMVcsr2",  kernelInfo);

      filename = oklpath + "SpMVell.okl";
      SpMVellKernel1  = device.buildKernel(filename.c_str(),  "SpMVell1",  kernelInfo);
      SpMVellKernel2  = device.buildKernel(filename.c_str(),  "SpMVell2",  kernelInfo);

      filename = oklpath + "SpMVmcsr.okl";
      SpMVmcsrKernel1 = device.buildKernel(filename.c_str(), "SpMVmcsr1", kernelInfo);
      SpMVmcsrKernel2 = device.buildKernel(filename.c_str(), "SpMVmcsr2", kernelInfo);

      filename = oklpath + "vectorSet.okl";
      vectorSetKernel = device.buildKernel(filename.c_str(), "vectorSet", kernelInfo);

      filename = oklpath + "vectorScale.okl";
      vectorScaleKernel = device.buildKernel(filename.c_str(), "vectorScale", kernelInfo);

      filename = oklpath + "vectorAddScalar.okl";
      vectorAddScalarKernel = device.buildKernel(filename.c_str(), "vectorAddScalar", kernelInfo);

      filename = oklpath + "vectorAdd.okl";
      vectorAddKernel1 = device.buildKernel(filename.c_str(), "vectorAdd1", kernelInfo);
      vectorAddKernel2 = device.buildKernel(filename.c_str(), "vectorAdd2", kernelInfo);

      filename = oklpath + "vectorInnerProd.okl";
      vectorInnerProdKernel = device.buildKernel(filename.c_str(), "vectorInnerProd", kernelInfo);

      filename = oklpath + "vectorAddInnerProd.okl";
      vectorAddInnerProdKernel = device.buildKernel(filename.c_str(), "vectorAddInnerProd", kernelInfo);
      vectorAddWeightedInnerProdKernel = device.buildKernel(filename.c_str(), "vectorAddWeightedInnerProd", kernelInfo);

      filename = oklpath + "kcycleCombinedOp.okl";
      kcycleCombinedOp1Kernel = device.buildKernel(filename.c_str(), "kcycleCombinedOp1", kernelInfo);
      kcycleCombinedOp2Kernel = device.buildKernel(filename.c_str(), "kcycleCombinedOp2", kernelInfo);
      kcycleWeightedCombinedOp1Kernel = device.buildKernel(filename.c_str(), "kcycleWeightedCombinedOp1", kernelInfo);
      kcycleWeightedCombinedOp2Kernel = device.buildKernel(filename.c_str(), "kcycleWeightedCombinedOp2", kernelInfo);

      filename = oklpath + "haloExtract.okl";
      haloExtractKernel = device.buildKernel(filename.c_str(), "haloExtract", kernelInfo);
      */
    }
    MPI_Barrier(comm);
  }
  MPI_Barrier(comm);
  if(rank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStart); fflush(stdout);
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
