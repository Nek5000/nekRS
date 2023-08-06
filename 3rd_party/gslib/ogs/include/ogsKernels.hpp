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

#ifndef OGS_KERNELS_HPP
#define OGS_KERNELS_HPP 1

#include <limits>
#include "ogs.hpp"

namespace ogs {

  extern const int gatherNodesPerBlock;

  extern int Nrefs;

  extern void* hostBuf;
  extern size_t hostBufSize;

  extern void* haloBuf;
  extern occa::memory o_haloBuf;
  extern occa::memory h_haloBuf;

  extern occa::kernel gatherScatterNewKernel_floatAdd;
  extern occa::kernel gatherScatterNewKernel_doubleAdd;
  extern occa::kernel gatherScatterNewKernel_doubleMin;
  extern occa::kernel gatherScatterNewKernel_doubleMax;

  extern occa::kernel gatherScatterKernel_floatAdd;
  extern occa::kernel gatherScatterKernel_floatMul;
  extern occa::kernel gatherScatterKernel_floatMin;
  extern occa::kernel gatherScatterKernel_floatMax;

  extern occa::kernel gatherScatterKernel_doubleAdd;
  extern occa::kernel gatherScatterKernel_doubleMul;
  extern occa::kernel gatherScatterKernel_doubleMin;
  extern occa::kernel gatherScatterKernel_doubleMax;

  extern occa::kernel gatherScatterKernel_intAdd;
  extern occa::kernel gatherScatterKernel_intMul;
  extern occa::kernel gatherScatterKernel_intMin;
  extern occa::kernel gatherScatterKernel_intMax;

  extern occa::kernel gatherScatterKernel_longAdd;
  extern occa::kernel gatherScatterKernel_longMul;
  extern occa::kernel gatherScatterKernel_longMin;
  extern occa::kernel gatherScatterKernel_longMax;



  extern occa::kernel gatherScatterVecKernel_floatAdd;
  extern occa::kernel gatherScatterVecKernel_floatMul;
  extern occa::kernel gatherScatterVecKernel_floatMin;
  extern occa::kernel gatherScatterVecKernel_floatMax;

  extern occa::kernel gatherScatterVecKernel_doubleAdd;
  extern occa::kernel gatherScatterVecKernel_doubleMul;
  extern occa::kernel gatherScatterVecKernel_doubleMin;
  extern occa::kernel gatherScatterVecKernel_doubleMax;

  extern occa::kernel gatherScatterVecKernel_intAdd;
  extern occa::kernel gatherScatterVecKernel_intMul;
  extern occa::kernel gatherScatterVecKernel_intMin;
  extern occa::kernel gatherScatterVecKernel_intMax;

  extern occa::kernel gatherScatterVecKernel_longAdd;
  extern occa::kernel gatherScatterVecKernel_longMul;
  extern occa::kernel gatherScatterVecKernel_longMin;
  extern occa::kernel gatherScatterVecKernel_longMax;



  extern occa::kernel gatherScatterManyKernel_floatAdd;
  extern occa::kernel gatherScatterManyKernel_floatMul;
  extern occa::kernel gatherScatterManyKernel_floatMin;
  extern occa::kernel gatherScatterManyKernel_floatMax;

  extern occa::kernel gatherScatterManyKernel_doubleAdd;
  extern occa::kernel gatherScatterManyKernel_doubleMul;
  extern occa::kernel gatherScatterManyKernel_doubleMin;
  extern occa::kernel gatherScatterManyKernel_doubleMax;

  extern occa::kernel gatherScatterManyKernel_intAdd;
  extern occa::kernel gatherScatterManyKernel_intMul;
  extern occa::kernel gatherScatterManyKernel_intMin;
  extern occa::kernel gatherScatterManyKernel_intMax;

  extern occa::kernel gatherScatterManyKernel_longAdd;
  extern occa::kernel gatherScatterManyKernel_longMul;
  extern occa::kernel gatherScatterManyKernel_longMin;
  extern occa::kernel gatherScatterManyKernel_longMax;



  extern occa::kernel gatherKernel_floatAdd;
  extern occa::kernel gatherKernel_floatMul;
  extern occa::kernel gatherKernel_floatMin;
  extern occa::kernel gatherKernel_floatMax;

  extern occa::kernel gatherKernel_doubleAdd;
  extern occa::kernel gatherKernel_doubleMul;
  extern occa::kernel gatherKernel_doubleMin;
  extern occa::kernel gatherKernel_doubleMax;

  extern occa::kernel gatherKernel_intAdd;
  extern occa::kernel gatherKernel_intMul;
  extern occa::kernel gatherKernel_intMin;
  extern occa::kernel gatherKernel_intMax;

  extern occa::kernel gatherKernel_longAdd;
  extern occa::kernel gatherKernel_longMul;
  extern occa::kernel gatherKernel_longMin;
  extern occa::kernel gatherKernel_longMax;


  extern occa::kernel gatherVecKernel_floatAdd;
  extern occa::kernel gatherVecKernel_floatMul;
  extern occa::kernel gatherVecKernel_floatMin;
  extern occa::kernel gatherVecKernel_floatMax;

  extern occa::kernel gatherVecKernel_doubleAdd;
  extern occa::kernel gatherVecKernel_doubleMul;
  extern occa::kernel gatherVecKernel_doubleMin;
  extern occa::kernel gatherVecKernel_doubleMax;

  extern occa::kernel gatherVecKernel_intAdd;
  extern occa::kernel gatherVecKernel_intMul;
  extern occa::kernel gatherVecKernel_intMin;
  extern occa::kernel gatherVecKernel_intMax;

  extern occa::kernel gatherVecKernel_longAdd;
  extern occa::kernel gatherVecKernel_longMul;
  extern occa::kernel gatherVecKernel_longMin;
  extern occa::kernel gatherVecKernel_longMax;


  extern occa::kernel gatherManyKernel_floatAdd;
  extern occa::kernel gatherManyKernel_floatMul;
  extern occa::kernel gatherManyKernel_floatMin;
  extern occa::kernel gatherManyKernel_floatMax;

  extern occa::kernel gatherManyKernel_doubleAdd;
  extern occa::kernel gatherManyKernel_doubleMul;
  extern occa::kernel gatherManyKernel_doubleMin;
  extern occa::kernel gatherManyKernel_doubleMax;

  extern occa::kernel gatherManyKernel_intAdd;
  extern occa::kernel gatherManyKernel_intMul;
  extern occa::kernel gatherManyKernel_intMin;
  extern occa::kernel gatherManyKernel_intMax;

  extern occa::kernel gatherManyKernel_longAdd;
  extern occa::kernel gatherManyKernel_longMul;
  extern occa::kernel gatherManyKernel_longMin;
  extern occa::kernel gatherManyKernel_longMax;


  extern occa::kernel scatterKernel_float;
  extern occa::kernel scatterKernel_double;
  extern occa::kernel scatterKernel_int;
  extern occa::kernel scatterKernel_long;

  extern occa::kernel scatterVecKernel_float;
  extern occa::kernel scatterVecKernel_double;
  extern occa::kernel scatterVecKernel_int;
  extern occa::kernel scatterVecKernel_long;

  extern occa::kernel scatterManyKernel_float;
  extern occa::kernel scatterManyKernel_double;
  extern occa::kernel scatterManyKernel_int;
  extern occa::kernel scatterManyKernel_long;

  extern occa::stream defaultStream;
  extern occa::stream dataStream;

  void initKernels(MPI_Comm comm, occa::device device, ogsBuildKernel_t buildKernel, bool verbose = false);

  extern occa::properties kernelInfo;

  void freeKernels();
}

void occaGatherScatter(const  dlong Ngather,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v);

void occaGatherScatterVec(const  dlong Ngather,
                const int Nentries,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v);

void occaGatherScatterMany(const  dlong Ngather,
                const int Nentries,
                const dlong stride,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v);

void occaGather(const  dlong Ngather,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_gv);

void occaGatherVec(const  dlong Ngather,
                const int Nentries,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_gv);

void occaGatherMany(const  dlong Ngather,
                const int Nentries,
                const dlong stride,
                const dlong gstride,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_gv);

void occaScatter(const  dlong Nscatter,
                occa::memory o_scatterStarts,
                occa::memory o_scatterIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_sv);

void occaScatterVec(const  dlong Nscatter,
                const int Nentries,
                occa::memory o_scatterStarts,
                occa::memory o_scatterIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_sv);

void occaScatterMany(const  dlong Nscatter,
                const int Nentries,
                const dlong stride,
                const dlong sstride,
                occa::memory o_scatterStarts,
                occa::memory o_scatterIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_sv);

#endif
