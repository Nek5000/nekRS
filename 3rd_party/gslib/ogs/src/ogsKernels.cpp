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

  occa::properties kernelInfo;

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
  occa::kernel gatherKernel_floatMul;
  occa::kernel gatherKernel_floatMin;
  occa::kernel gatherKernel_floatMax;
  occa::kernel gatherKernel_doubleAdd;
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

  ogs::kernelInfo["defines"].asObject();
  ogs::kernelInfo["includes"].asArray();
  ogs::kernelInfo["header"].asArray();
  ogs::kernelInfo["flags"].asObject();

  if(sizeof(dlong)==4){
   ogs::kernelInfo["defines/" "dlong"]="int";
  }
  if(sizeof(dlong)==8){
   ogs::kernelInfo["defines/" "dlong"]="long long int";
  }

  if(sizeof(dfloat) == sizeof(double)){
   ogs::kernelInfo["defines/" "dfloat"]= "double";
   ogs::kernelInfo["defines/" "dfloat4"]= "double4";
  }
  else if(sizeof(dfloat) == sizeof(float)){
   ogs::kernelInfo["defines/" "dfloat"]= "float";
   ogs::kernelInfo["defines/" "dfloat4"]= "float4";
  }


  if(device.mode()=="OpenCL"){
   //ogs::kernelInfo["compiler_flags"] += "-cl-opt-disable";
  }

  if(device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
   ogs::kernelInfo["compiler_flags"] += " --ftz=true ";
   ogs::kernelInfo["compiler_flags"] += " --prec-div=false ";
   ogs::kernelInfo["compiler_flags"] += " --prec-sqrt=false ";
   ogs::kernelInfo["compiler_flags"] += " --use_fast_math ";
   ogs::kernelInfo["compiler_flags"] += " --fmad=true "; // compiler option for cuda
  }

  if (rank==0) printf("Compiling GatherScatter Kernels...");fflush(stdout);

  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {      

      ogs::gatherScatterKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatAdd", ogs::kernelInfo);
      ogs::gatherScatterKernel_floatMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatMul", ogs::kernelInfo);
      ogs::gatherScatterKernel_floatMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatMin", ogs::kernelInfo);
      ogs::gatherScatterKernel_floatMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatMax", ogs::kernelInfo);

      ogs::gatherScatterKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleAdd", ogs::kernelInfo);
      ogs::gatherScatterKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleMul", ogs::kernelInfo);
      ogs::gatherScatterKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleMin", ogs::kernelInfo);
      ogs::gatherScatterKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleMax", ogs::kernelInfo);

      ogs::gatherScatterKernel_intAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intAdd", ogs::kernelInfo);
      ogs::gatherScatterKernel_intMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intMul", ogs::kernelInfo);
      ogs::gatherScatterKernel_intMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intMin", ogs::kernelInfo);
      ogs::gatherScatterKernel_intMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intMax", ogs::kernelInfo);

      ogs::gatherScatterKernel_longAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longAdd", ogs::kernelInfo);
      ogs::gatherScatterKernel_longMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longMul", ogs::kernelInfo);
      ogs::gatherScatterKernel_longMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longMin", ogs::kernelInfo);
      ogs::gatherScatterKernel_longMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longMax", ogs::kernelInfo);

      ogs::gatherScatterVecKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatAdd", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_floatMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatMul", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_floatMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatMin", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_floatMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatMax", ogs::kernelInfo);

      ogs::gatherScatterVecKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleAdd", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleMul", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleMin", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleMax", ogs::kernelInfo);

      ogs::gatherScatterVecKernel_intAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intAdd", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_intMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intMul", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_intMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intMin", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_intMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intMax", ogs::kernelInfo);

      ogs::gatherScatterVecKernel_longAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longAdd", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_longMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longMul", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_longMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longMin", ogs::kernelInfo);
      ogs::gatherScatterVecKernel_longMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longMax", ogs::kernelInfo);

      ogs::gatherScatterManyKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatAdd", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_floatMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatMul", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_floatMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatMin", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_floatMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatMax", ogs::kernelInfo);

      ogs::gatherScatterManyKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleAdd", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleMul", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleMin", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleMax", ogs::kernelInfo);

      ogs::gatherScatterManyKernel_intAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intAdd", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_intMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intMul", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_intMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intMin", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_intMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intMax", ogs::kernelInfo);

      ogs::gatherScatterManyKernel_longAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longAdd", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_longMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longMul", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_longMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longMin", ogs::kernelInfo);
      ogs::gatherScatterManyKernel_longMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longMax", ogs::kernelInfo);



      ogs::gatherKernel_floatAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatAdd", ogs::kernelInfo);
      ogs::gatherKernel_floatMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatMul", ogs::kernelInfo);
      ogs::gatherKernel_floatMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatMin", ogs::kernelInfo);
      ogs::gatherKernel_floatMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatMax", ogs::kernelInfo);

      ogs::gatherKernel_doubleAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleAdd", ogs::kernelInfo);
      ogs::gatherKernel_doubleMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleMul", ogs::kernelInfo);
      ogs::gatherKernel_doubleMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleMin", ogs::kernelInfo);
      ogs::gatherKernel_doubleMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleMax", ogs::kernelInfo);

      ogs::gatherKernel_intAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_intAdd", ogs::kernelInfo);
      ogs::gatherKernel_intMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_intMul", ogs::kernelInfo);
      ogs::gatherKernel_intMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_intMin", ogs::kernelInfo);
      ogs::gatherKernel_intMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_intMax", ogs::kernelInfo);

      ogs::gatherKernel_longAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_longAdd", ogs::kernelInfo);
      ogs::gatherKernel_longMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_longMul", ogs::kernelInfo);
      ogs::gatherKernel_longMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_longMin", ogs::kernelInfo);
      ogs::gatherKernel_longMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_longMax", ogs::kernelInfo);

      ogs::gatherVecKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatAdd", ogs::kernelInfo);
      ogs::gatherVecKernel_floatMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatMul", ogs::kernelInfo);
      ogs::gatherVecKernel_floatMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatMin", ogs::kernelInfo);
      ogs::gatherVecKernel_floatMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatMax", ogs::kernelInfo);

      ogs::gatherVecKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleAdd", ogs::kernelInfo);
      ogs::gatherVecKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleMul", ogs::kernelInfo);
      ogs::gatherVecKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleMin", ogs::kernelInfo);
      ogs::gatherVecKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleMax", ogs::kernelInfo);

      ogs::gatherVecKernel_intAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intAdd", ogs::kernelInfo);
      ogs::gatherVecKernel_intMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intMul", ogs::kernelInfo);
      ogs::gatherVecKernel_intMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intMin", ogs::kernelInfo);
      ogs::gatherVecKernel_intMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intMax", ogs::kernelInfo);

      ogs::gatherVecKernel_longAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longAdd", ogs::kernelInfo);
      ogs::gatherVecKernel_longMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longMul", ogs::kernelInfo);
      ogs::gatherVecKernel_longMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longMin", ogs::kernelInfo);
      ogs::gatherVecKernel_longMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longMax", ogs::kernelInfo);

      ogs::gatherManyKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatAdd", ogs::kernelInfo);
      ogs::gatherManyKernel_floatMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatMul", ogs::kernelInfo);
      ogs::gatherManyKernel_floatMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatMin", ogs::kernelInfo);
      ogs::gatherManyKernel_floatMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatMax", ogs::kernelInfo);

      ogs::gatherManyKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleAdd", ogs::kernelInfo);
      ogs::gatherManyKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleMul", ogs::kernelInfo);
      ogs::gatherManyKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleMin", ogs::kernelInfo);
      ogs::gatherManyKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleMax", ogs::kernelInfo);

      ogs::gatherManyKernel_intAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intAdd", ogs::kernelInfo);
      ogs::gatherManyKernel_intMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intMul", ogs::kernelInfo);
      ogs::gatherManyKernel_intMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intMin", ogs::kernelInfo);
      ogs::gatherManyKernel_intMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intMax", ogs::kernelInfo);

      ogs::gatherManyKernel_longAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longAdd", ogs::kernelInfo);
      ogs::gatherManyKernel_longMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longMul", ogs::kernelInfo);
      ogs::gatherManyKernel_longMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longMin", ogs::kernelInfo);
      ogs::gatherManyKernel_longMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longMax", ogs::kernelInfo);



      ogs::scatterKernel_float = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_float", ogs::kernelInfo);
      ogs::scatterKernel_double = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_double", ogs::kernelInfo);
      ogs::scatterKernel_int = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_int", ogs::kernelInfo);
      ogs::scatterKernel_long = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_long", ogs::kernelInfo);

      ogs::scatterVecKernel_float = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_float", ogs::kernelInfo);
      ogs::scatterVecKernel_double = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_double", ogs::kernelInfo);
      ogs::scatterVecKernel_int = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_int", ogs::kernelInfo);
      ogs::scatterVecKernel_long = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_long", ogs::kernelInfo);

      ogs::scatterManyKernel_float = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_float", ogs::kernelInfo);
      ogs::scatterManyKernel_double = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_double", ogs::kernelInfo);
      ogs::scatterManyKernel_int = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_int", ogs::kernelInfo);
      ogs::scatterManyKernel_long = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_long", ogs::kernelInfo);
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

