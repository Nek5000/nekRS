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

#include "ogs.hpp"
#include "ogsKernels.hpp"

namespace ogs {

  int Nrefs = 0;

  void* hostBuf;
  size_t hostBufSize=0;

  void* haloBuf;
  occa::memory o_haloBuf;
  occa::memory h_haloBuf;

  occa::stream defaultStream;
  occa::stream dataStream;

  occa::kernel gatherScatterKernel_floatAdd;
  occa::kernel gatherScatterKernel_floatMul;
  occa::kernel gatherScatterKernel_floatMin;
  occa::kernel gatherScatterKernel_floatMax;
  occa::kernel gatherScatterKernel_doubleAdd;
  occa::kernel gatherScatterKernel_doubleMul;
  occa::kernel gatherScatterKernel_doubleMin;
  occa::kernel gatherScatterKernel_doubleMax;
  occa::kernel gatherScatterKernel_intAdd;
  occa::kernel gatherScatterKernel_intMul;
  occa::kernel gatherScatterKernel_intMin;
  occa::kernel gatherScatterKernel_intMax;
  occa::kernel gatherScatterKernel_longAdd;
  occa::kernel gatherScatterKernel_longMul;
  occa::kernel gatherScatterKernel_longMin;
  occa::kernel gatherScatterKernel_longMax;
  occa::kernel gatherScatterVecKernel_floatAdd;
  occa::kernel gatherScatterVecKernel_floatMul;
  occa::kernel gatherScatterVecKernel_floatMin;
  occa::kernel gatherScatterVecKernel_floatMax;
  occa::kernel gatherScatterVecKernel_doubleAdd;
  occa::kernel gatherScatterVecKernel_doubleMul;
  occa::kernel gatherScatterVecKernel_doubleMin;
  occa::kernel gatherScatterVecKernel_doubleMax;
  occa::kernel gatherScatterVecKernel_intAdd;
  occa::kernel gatherScatterVecKernel_intMul;
  occa::kernel gatherScatterVecKernel_intMin;
  occa::kernel gatherScatterVecKernel_intMax;
  occa::kernel gatherScatterVecKernel_longAdd;
  occa::kernel gatherScatterVecKernel_longMul;
  occa::kernel gatherScatterVecKernel_longMin;
  occa::kernel gatherScatterVecKernel_longMax;
  occa::kernel gatherScatterManyKernel_floatAdd;
  occa::kernel gatherScatterManyKernel_floatMul;
  occa::kernel gatherScatterManyKernel_floatMin;
  occa::kernel gatherScatterManyKernel_floatMax;
  occa::kernel gatherScatterManyKernel_doubleAdd;
  occa::kernel gatherScatterManyKernel_doubleMul;
  occa::kernel gatherScatterManyKernel_doubleMin;
  occa::kernel gatherScatterManyKernel_doubleMax;
  occa::kernel gatherScatterManyKernel_intAdd;
  occa::kernel gatherScatterManyKernel_intMul;
  occa::kernel gatherScatterManyKernel_intMin;
  occa::kernel gatherScatterManyKernel_intMax;
  occa::kernel gatherScatterManyKernel_longAdd;
  occa::kernel gatherScatterManyKernel_longMul;
  occa::kernel gatherScatterManyKernel_longMin;
  occa::kernel gatherScatterManyKernel_longMax;

  occa::kernel gatherKernel_floatAdd;
  occa::kernel gatherKernel_floatAddSelf;
  occa::kernel gatherKernel_floatMul;
  occa::kernel gatherKernel_floatMin;
  occa::kernel gatherKernel_floatMax;
  occa::kernel gatherKernel_doubleAdd;
  occa::kernel gatherKernel_doubleAddSelf;
  occa::kernel gatherKernel_doubleMul;
  occa::kernel gatherKernel_doubleMin;
  occa::kernel gatherKernel_doubleMax;
  occa::kernel gatherKernel_intAdd;
  occa::kernel gatherKernel_intMul;
  occa::kernel gatherKernel_intMin;
  occa::kernel gatherKernel_intMax;
  occa::kernel gatherKernel_longAdd;
  occa::kernel gatherKernel_longMul;
  occa::kernel gatherKernel_longMin;
  occa::kernel gatherKernel_longMax;
  occa::kernel gatherVecKernel_floatAdd;
  occa::kernel gatherVecKernel_floatMul;
  occa::kernel gatherVecKernel_floatMin;
  occa::kernel gatherVecKernel_floatMax;
  occa::kernel gatherVecKernel_doubleAdd;
  occa::kernel gatherVecKernel_doubleMul;
  occa::kernel gatherVecKernel_doubleMin;
  occa::kernel gatherVecKernel_doubleMax;
  occa::kernel gatherVecKernel_intAdd;
  occa::kernel gatherVecKernel_intMul;
  occa::kernel gatherVecKernel_intMin;
  occa::kernel gatherVecKernel_intMax;
  occa::kernel gatherVecKernel_longAdd;
  occa::kernel gatherVecKernel_longMul;
  occa::kernel gatherVecKernel_longMin;
  occa::kernel gatherVecKernel_longMax;
  occa::kernel gatherManyKernel_floatAdd;
  occa::kernel gatherManyKernel_floatMul;
  occa::kernel gatherManyKernel_floatMin;
  occa::kernel gatherManyKernel_floatMax;
  occa::kernel gatherManyKernel_doubleAdd;
  occa::kernel gatherManyKernel_doubleMul;
  occa::kernel gatherManyKernel_doubleMin;
  occa::kernel gatherManyKernel_doubleMax;
  occa::kernel gatherManyKernel_intAdd;
  occa::kernel gatherManyKernel_intMul;
  occa::kernel gatherManyKernel_intMin;
  occa::kernel gatherManyKernel_intMax;
  occa::kernel gatherManyKernel_longAdd;
  occa::kernel gatherManyKernel_longMul;
  occa::kernel gatherManyKernel_longMin;
  occa::kernel gatherManyKernel_longMax;

  occa::kernel scatterKernel_float;
  occa::kernel scatterKernel_double;
  occa::kernel scatterKernel_int;
  occa::kernel scatterKernel_long;
  occa::kernel scatterVecKernel_float;
  occa::kernel scatterVecKernel_double;
  occa::kernel scatterVecKernel_int;
  occa::kernel scatterVecKernel_long;
  occa::kernel scatterManyKernel_float;
  occa::kernel scatterManyKernel_double;
  occa::kernel scatterManyKernel_int;
  occa::kernel scatterManyKernel_long;
}


void ogs::initKernels(MPI_Comm comm, occa::device device) {

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  ogs::defaultStream = device.getStream();
  ogs::dataStream    = device.createStream();

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

  if (rank==0) printf("Compiling GatherScatter Kernels...");fflush(stdout);

  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {      

      ogs::gatherScatterKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatAdd", kernelInfo);
      ogs::gatherScatterKernel_floatMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatMul", kernelInfo);
      ogs::gatherScatterKernel_floatMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatMin", kernelInfo);
      ogs::gatherScatterKernel_floatMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatMax", kernelInfo);

      ogs::gatherScatterKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleAdd", kernelInfo);
      ogs::gatherScatterKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleMul", kernelInfo);
      ogs::gatherScatterKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleMin", kernelInfo);
      ogs::gatherScatterKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleMax", kernelInfo);

      ogs::gatherScatterKernel_intAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intAdd", kernelInfo);
      ogs::gatherScatterKernel_intMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intMul", kernelInfo);
      ogs::gatherScatterKernel_intMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intMin", kernelInfo);
      ogs::gatherScatterKernel_intMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intMax", kernelInfo);

      ogs::gatherScatterKernel_longAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longAdd", kernelInfo);
      ogs::gatherScatterKernel_longMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longMul", kernelInfo);
      ogs::gatherScatterKernel_longMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longMin", kernelInfo);
      ogs::gatherScatterKernel_longMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longMax", kernelInfo);

      ogs::gatherScatterVecKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatAdd", kernelInfo);
      ogs::gatherScatterVecKernel_floatMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatMul", kernelInfo);
      ogs::gatherScatterVecKernel_floatMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatMin", kernelInfo);
      ogs::gatherScatterVecKernel_floatMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatMax", kernelInfo);

      ogs::gatherScatterVecKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleAdd", kernelInfo);
      ogs::gatherScatterVecKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleMul", kernelInfo);
      ogs::gatherScatterVecKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleMin", kernelInfo);
      ogs::gatherScatterVecKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleMax", kernelInfo);

      ogs::gatherScatterVecKernel_intAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intAdd", kernelInfo);
      ogs::gatherScatterVecKernel_intMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intMul", kernelInfo);
      ogs::gatherScatterVecKernel_intMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intMin", kernelInfo);
      ogs::gatherScatterVecKernel_intMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intMax", kernelInfo);

      ogs::gatherScatterVecKernel_longAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longAdd", kernelInfo);
      ogs::gatherScatterVecKernel_longMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longMul", kernelInfo);
      ogs::gatherScatterVecKernel_longMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longMin", kernelInfo);
      ogs::gatherScatterVecKernel_longMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longMax", kernelInfo);

      ogs::gatherScatterManyKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatAdd", kernelInfo);
      ogs::gatherScatterManyKernel_floatMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatMul", kernelInfo);
      ogs::gatherScatterManyKernel_floatMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatMin", kernelInfo);
      ogs::gatherScatterManyKernel_floatMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatMax", kernelInfo);

      ogs::gatherScatterManyKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleAdd", kernelInfo);
      ogs::gatherScatterManyKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleMul", kernelInfo);
      ogs::gatherScatterManyKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleMin", kernelInfo);
      ogs::gatherScatterManyKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleMax", kernelInfo);

      ogs::gatherScatterManyKernel_intAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intAdd", kernelInfo);
      ogs::gatherScatterManyKernel_intMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intMul", kernelInfo);
      ogs::gatherScatterManyKernel_intMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intMin", kernelInfo);
      ogs::gatherScatterManyKernel_intMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intMax", kernelInfo);

      ogs::gatherScatterManyKernel_longAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longAdd", kernelInfo);
      ogs::gatherScatterManyKernel_longMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longMul", kernelInfo);
      ogs::gatherScatterManyKernel_longMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longMin", kernelInfo);
      ogs::gatherScatterManyKernel_longMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longMax", kernelInfo);



      ogs::gatherKernel_floatAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatAdd", kernelInfo);
      ogs::gatherKernel_floatAddSelf = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatAddSelf", kernelInfo);
      ogs::gatherKernel_floatMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatMul", kernelInfo);
      ogs::gatherKernel_floatMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatMin", kernelInfo);
      ogs::gatherKernel_floatMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatMax", kernelInfo);

      ogs::gatherKernel_doubleAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleAdd", kernelInfo);
      ogs::gatherKernel_doubleAddSelf = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleAddSelf", kernelInfo);
      ogs::gatherKernel_doubleMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleMul", kernelInfo);
      ogs::gatherKernel_doubleMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleMin", kernelInfo);
      ogs::gatherKernel_doubleMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleMax", kernelInfo);

      ogs::gatherKernel_intAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_intAdd", kernelInfo);
      ogs::gatherKernel_intMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_intMul", kernelInfo);
      ogs::gatherKernel_intMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_intMin", kernelInfo);
      ogs::gatherKernel_intMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_intMax", kernelInfo);

      ogs::gatherKernel_longAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_longAdd", kernelInfo);
      ogs::gatherKernel_longMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_longMul", kernelInfo);
      ogs::gatherKernel_longMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_longMin", kernelInfo);
      ogs::gatherKernel_longMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_longMax", kernelInfo);

      ogs::gatherVecKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatAdd", kernelInfo);
      ogs::gatherVecKernel_floatMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatMul", kernelInfo);
      ogs::gatherVecKernel_floatMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatMin", kernelInfo);
      ogs::gatherVecKernel_floatMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatMax", kernelInfo);

      ogs::gatherVecKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleAdd", kernelInfo);
      ogs::gatherVecKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleMul", kernelInfo);
      ogs::gatherVecKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleMin", kernelInfo);
      ogs::gatherVecKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleMax", kernelInfo);

      ogs::gatherVecKernel_intAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intAdd", kernelInfo);
      ogs::gatherVecKernel_intMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intMul", kernelInfo);
      ogs::gatherVecKernel_intMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intMin", kernelInfo);
      ogs::gatherVecKernel_intMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intMax", kernelInfo);

      ogs::gatherVecKernel_longAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longAdd", kernelInfo);
      ogs::gatherVecKernel_longMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longMul", kernelInfo);
      ogs::gatherVecKernel_longMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longMin", kernelInfo);
      ogs::gatherVecKernel_longMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longMax", kernelInfo);

      ogs::gatherManyKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatAdd", kernelInfo);
      ogs::gatherManyKernel_floatMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatMul", kernelInfo);
      ogs::gatherManyKernel_floatMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatMin", kernelInfo);
      ogs::gatherManyKernel_floatMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatMax", kernelInfo);

      ogs::gatherManyKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleAdd", kernelInfo);
      ogs::gatherManyKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleMul", kernelInfo);
      ogs::gatherManyKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleMin", kernelInfo);
      ogs::gatherManyKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleMax", kernelInfo);

      ogs::gatherManyKernel_intAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intAdd", kernelInfo);
      ogs::gatherManyKernel_intMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intMul", kernelInfo);
      ogs::gatherManyKernel_intMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intMin", kernelInfo);
      ogs::gatherManyKernel_intMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intMax", kernelInfo);

      ogs::gatherManyKernel_longAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longAdd", kernelInfo);
      ogs::gatherManyKernel_longMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longMul", kernelInfo);
      ogs::gatherManyKernel_longMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longMin", kernelInfo);
      ogs::gatherManyKernel_longMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longMax", kernelInfo);



      ogs::scatterKernel_float = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_float", kernelInfo);
      ogs::scatterKernel_double = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_double", kernelInfo);
      ogs::scatterKernel_int = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_int", kernelInfo);
      ogs::scatterKernel_long = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_long", kernelInfo);

      ogs::scatterVecKernel_float = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_float", kernelInfo);
      ogs::scatterVecKernel_double = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_double", kernelInfo);
      ogs::scatterVecKernel_int = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_int", kernelInfo);
      ogs::scatterVecKernel_long = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_long", kernelInfo);

      ogs::scatterManyKernel_float = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_float", kernelInfo);
      ogs::scatterManyKernel_double = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_double", kernelInfo);
      ogs::scatterManyKernel_int = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_int", kernelInfo);
      ogs::scatterManyKernel_long = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_long", kernelInfo);
    }
    MPI_Barrier(comm);
  }
  if(rank==0) printf("done.\n");
}

void ogs::freeKernels() {

  ogs::gatherScatterKernel_floatAdd.free();
  ogs::gatherScatterKernel_floatMul.free();
  ogs::gatherScatterKernel_floatMin.free();
  ogs::gatherScatterKernel_floatMax.free();
  ogs::gatherScatterKernel_doubleAdd.free();
  ogs::gatherScatterKernel_doubleMul.free();
  ogs::gatherScatterKernel_doubleMin.free();
  ogs::gatherScatterKernel_doubleMax.free();
  ogs::gatherScatterKernel_intAdd.free();
  ogs::gatherScatterKernel_intMul.free();
  ogs::gatherScatterKernel_intMin.free();
  ogs::gatherScatterKernel_intMax.free();
  ogs::gatherScatterKernel_longAdd.free();
  ogs::gatherScatterKernel_longMul.free();
  ogs::gatherScatterKernel_longMin.free();
  ogs::gatherScatterKernel_longMax.free();
  ogs::gatherScatterVecKernel_floatAdd.free();
  ogs::gatherScatterVecKernel_floatMul.free();
  ogs::gatherScatterVecKernel_floatMin.free();
  ogs::gatherScatterVecKernel_floatMax.free();
  ogs::gatherScatterVecKernel_doubleAdd.free();
  ogs::gatherScatterVecKernel_doubleMul.free();
  ogs::gatherScatterVecKernel_doubleMin.free();
  ogs::gatherScatterVecKernel_doubleMax.free();
  ogs::gatherScatterVecKernel_intAdd.free();
  ogs::gatherScatterVecKernel_intMul.free();
  ogs::gatherScatterVecKernel_intMin.free();
  ogs::gatherScatterVecKernel_intMax.free();
  ogs::gatherScatterVecKernel_longAdd.free();
  ogs::gatherScatterVecKernel_longMul.free();
  ogs::gatherScatterVecKernel_longMin.free();
  ogs::gatherScatterVecKernel_longMax.free();
  ogs::gatherScatterManyKernel_floatAdd.free();
  ogs::gatherScatterManyKernel_floatMul.free();
  ogs::gatherScatterManyKernel_floatMin.free();
  ogs::gatherScatterManyKernel_floatMax.free();
  ogs::gatherScatterManyKernel_doubleAdd.free();
  ogs::gatherScatterManyKernel_doubleMul.free();
  ogs::gatherScatterManyKernel_doubleMin.free();
  ogs::gatherScatterManyKernel_doubleMax.free();
  ogs::gatherScatterManyKernel_intAdd.free();
  ogs::gatherScatterManyKernel_intMul.free();
  ogs::gatherScatterManyKernel_intMin.free();
  ogs::gatherScatterManyKernel_intMax.free();
  ogs::gatherScatterManyKernel_longAdd.free();
  ogs::gatherScatterManyKernel_longMul.free();
  ogs::gatherScatterManyKernel_longMin.free();
  ogs::gatherScatterManyKernel_longMax.free();

  ogs::gatherKernel_floatAdd.free();
  ogs::gatherKernel_floatMul.free();
  ogs::gatherKernel_floatMin.free();
  ogs::gatherKernel_floatMax.free();
  ogs::gatherKernel_doubleAdd.free();
  ogs::gatherKernel_doubleMul.free();
  ogs::gatherKernel_doubleMin.free();
  ogs::gatherKernel_doubleMax.free();
  ogs::gatherKernel_intAdd.free();
  ogs::gatherKernel_intMul.free();
  ogs::gatherKernel_intMin.free();
  ogs::gatherKernel_intMax.free();
  ogs::gatherKernel_longAdd.free();
  ogs::gatherKernel_longMul.free();
  ogs::gatherKernel_longMin.free();
  ogs::gatherKernel_longMax.free();
  ogs::gatherVecKernel_floatAdd.free();
  ogs::gatherVecKernel_floatMul.free();
  ogs::gatherVecKernel_floatMin.free();
  ogs::gatherVecKernel_floatMax.free();
  ogs::gatherVecKernel_doubleAdd.free();
  ogs::gatherVecKernel_doubleMul.free();
  ogs::gatherVecKernel_doubleMin.free();
  ogs::gatherVecKernel_doubleMax.free();
  ogs::gatherVecKernel_intAdd.free();
  ogs::gatherVecKernel_intMul.free();
  ogs::gatherVecKernel_intMin.free();
  ogs::gatherVecKernel_intMax.free();
  ogs::gatherVecKernel_longAdd.free();
  ogs::gatherVecKernel_longMul.free();
  ogs::gatherVecKernel_longMin.free();
  ogs::gatherVecKernel_longMax.free();
  ogs::gatherManyKernel_floatAdd.free();
  ogs::gatherManyKernel_floatMul.free();
  ogs::gatherManyKernel_floatMin.free();
  ogs::gatherManyKernel_floatMax.free();
  ogs::gatherManyKernel_doubleAdd.free();
  ogs::gatherManyKernel_doubleMul.free();
  ogs::gatherManyKernel_doubleMin.free();
  ogs::gatherManyKernel_doubleMax.free();
  ogs::gatherManyKernel_intAdd.free();
  ogs::gatherManyKernel_intMul.free();
  ogs::gatherManyKernel_intMin.free();
  ogs::gatherManyKernel_intMax.free();
  ogs::gatherManyKernel_longAdd.free();
  ogs::gatherManyKernel_longMul.free();
  ogs::gatherManyKernel_longMin.free();
  ogs::gatherManyKernel_longMax.free();

  ogs::scatterKernel_float.free();
  ogs::scatterKernel_double.free();
  ogs::scatterKernel_int.free();
  ogs::scatterKernel_long.free();
  ogs::scatterVecKernel_float.free();
  ogs::scatterVecKernel_double.free();
  ogs::scatterVecKernel_int.free();
  ogs::scatterVecKernel_long.free();
  ogs::scatterManyKernel_float.free();
  ogs::scatterManyKernel_double.free();
  ogs::scatterManyKernel_int.free();
  ogs::scatterManyKernel_long.free();

  ogs::o_haloBuf.free();
  ogs::haloBuf = NULL;
}

