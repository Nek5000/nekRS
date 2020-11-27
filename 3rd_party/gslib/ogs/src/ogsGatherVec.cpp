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

#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "ogsInterface.h"

#include "gatherVec.tpp"

void ogsGatherVec_add(void *gv, void *v, const int k, const size_t Nbytes, const char *type, ogs_t *ogs);
void ogsGatherVec_mul(void *gv, void *v, const int k, const size_t Nbytes, const char *type, ogs_t *ogs);
void ogsGatherVec_min(void *gv, void *v, const int k, const size_t Nbytes, const char *type, ogs_t *ogs);
void ogsGatherVec_max(void *gv, void *v, const int k, const size_t Nbytes, const char *type, ogs_t *ogs);

void ogsGatherVec(occa::memory o_gv, 
               occa::memory o_v,
               const int k, 
               const char *type, 
               const char *op, 
               ogs_t *ogs){
  ogsGatherVecStart (o_gv, o_v, k, type, op, ogs);
  ogsGatherVecFinish(o_gv, o_v, k, type, op, ogs);
}

void ogsGatherVecStart(occa::memory o_gv, 
                    occa::memory o_v, 
                    const int k,
                    const char *type, 
                    const char *op, 
                    ogs_t *ogs){
  size_t Nbytes;
  if (!strcmp(type, "float")) 
    Nbytes = sizeof(float);
  else if (!strcmp(type, "double")) 
    Nbytes = sizeof(double);
  else if (!strcmp(type, "int")) 
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int")) 
    Nbytes = sizeof(long long int);

  if (ogs->NhaloGather) {
    if (ogs::o_haloBuf.size() < ogs->NhaloGather*Nbytes*k) {
      if (ogs::o_haloBuf.size()) ogs::o_haloBuf.free();
      //      ogs::o_haloBuf = ogs->device.mappedAlloc(ogs->NhaloGather*Nbytes*k);
      //      ogs::haloBuf = ogs::o_haloBuf.getMappedPointer();
      ogs::haloBuf = ogsHostMallocPinned(ogs->device, ogs->NhaloGather*Nbytes*k, NULL, ogs::o_haloBuf, ogs::h_haloBuf);
    }
  }

  // gather halo nodes on device
  if (ogs->NhaloGather) {
    occaGatherVec(ogs->NhaloGather, k, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, o_v, ogs::o_haloBuf);
    
    ogs->device.finish();
    ogs->device.setStream(ogs::dataStream);
    ogs::o_haloBuf.copyTo(ogs::haloBuf, ogs->NhaloGather*Nbytes*k, 0, "async: true");
    ogs->device.setStream(ogs::defaultStream);
  }
}


void ogsGatherVecFinish(occa::memory o_gv, 
                     occa::memory o_v, 
                     const int k,
                     const char *type, 
                     const char *op, 
                     ogs_t *ogs){
  size_t Nbytes;
  if (!strcmp(type, "float")) 
    Nbytes = sizeof(float);
  else if (!strcmp(type, "double")) 
    Nbytes = sizeof(double);
  else if (!strcmp(type, "int")) 
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int")) 
    Nbytes = sizeof(long long int);

  if(ogs->NlocalGather) {
    occaGatherVec(ogs->NlocalGather, k, ogs->o_localGatherOffsets, ogs->o_localGatherIds, type, op, o_v, o_gv);
  }

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);
    ogs->device.finish();

    // MPI based gather using libgs
    ogsHostGatherVec(ogs::haloBuf, k, type, op, ogs->haloGshNonSym);

    // copy totally gather halo data back from HOST to DEVICE
    if (ogs->NownedHalo)
      o_gv.copyFrom(ogs::haloBuf, ogs->NownedHalo*Nbytes*k, 
                              ogs->NlocalGather*Nbytes*k, "async: true");

    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);
  }
}

void ogsGatherVec(void *gv, 
               void *v,
               const int k, 
               const char *type, 
               const char *op, 
               ogs_t *ogs){
  size_t Nbytes;
  if (!strcmp(type, "float")) 
    Nbytes = sizeof(float);
  else if (!strcmp(type, "double")) 
    Nbytes = sizeof(double);
  else if (!strcmp(type, "int")) 
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int")) 
    Nbytes = sizeof(long long int);

  if (ogs->NhaloGather) {
    if (ogs::hostBufSize < ogs->NhaloGather*Nbytes*k) {
      if (ogs::hostBufSize) free(ogs::hostBuf);
      ogs::hostBuf = (void *) malloc(ogs->NhaloGather*Nbytes*k);
    }
  }

  if (!strcmp(op, "add")) 
    ogsGatherVec_add(gv, v, k, Nbytes, type, ogs);
  else if (!strcmp(op, "mul")) 
    ogsGatherVec_mul(gv, v, k, Nbytes, type, ogs);
  else if (!strcmp(op, "min")) 
    ogsGatherVec_min(gv, v, k, Nbytes, type, ogs);
  else if (!strcmp(op, "max")) 
    ogsGatherVec_max(gv, v, k, Nbytes, type, ogs);
}

void ogsGatherVec_add(void *gv, void *v, const int k, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gatherVec_add<float>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gatherVec_add<double>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gatherVec_add<int>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gatherVec_add<long long int>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    ogsHostGatherVec(ogs::hostBuf, k, type, ogsAdd, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      memcpy((char*)gv+ogs->NlocalGather*Nbytes*k, ogs::hostBuf, ogs->NownedHalo*Nbytes*k);
  }

  if (!strcmp(type, "float")) 
    gatherVec_add<float>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gatherVec_add<double>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gatherVec_add<int>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gatherVec_add<long long int>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void ogsGatherVec_mul(void *gv, void *v, const int k, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gatherVec_mul<float>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gatherVec_mul<double>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gatherVec_mul<int>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gatherVec_mul<long long int>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    ogsHostGatherVec(ogs::hostBuf, k, type, ogsMul, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      memcpy((char*)gv+ogs->NlocalGather*Nbytes*k, ogs::hostBuf, ogs->NownedHalo*Nbytes*k);
  }

  if (!strcmp(type, "float")) 
    gatherVec_mul<float>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gatherVec_mul<double>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gatherVec_mul<int>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gatherVec_mul<long long int>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void ogsGatherVec_min(void *gv, void *v, const int k, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gatherVec_min<float>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gatherVec_min<double>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gatherVec_min<int>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gatherVec_min<long long int>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    ogsHostGatherVec(ogs::hostBuf, k, type, ogsMin, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      memcpy((char*)gv+ogs->NlocalGather*Nbytes*k, ogs::hostBuf, ogs->NownedHalo*Nbytes*k);
  }

  if (!strcmp(type, "float")) 
    gatherVec_min<float>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gatherVec_min<double>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gatherVec_min<int>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gatherVec_min<long long int>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void ogsGatherVec_max(void *gv, void *v, const int k, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gatherVec_max<float>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gatherVec_max<double>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gatherVec_max<int>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gatherVec_max<long long int>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    ogsHostGatherVec(ogs::hostBuf, k, type, ogsMax, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      memcpy((char*)gv+ogs->NlocalGather*Nbytes*k, ogs::hostBuf, ogs->NownedHalo*Nbytes*k);
  }

  if (!strcmp(type, "float")) 
    gatherVec_max<float>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gatherVec_max<double>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gatherVec_max<int>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gatherVec_max<long long int>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void occaGatherVec(const  dlong Ngather,
                const int Nentries,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_gv) {
  
  if      ((!strcmp(type, "float"))&&(!strcmp(op, "add"))) 
    ogs::gatherVecKernel_floatAdd(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "mul"))) 
    ogs::gatherVecKernel_floatMul(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "min"))) 
    ogs::gatherVecKernel_floatMin(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "max"))) 
    ogs::gatherVecKernel_floatMax(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "add"))) 
    ogs::gatherVecKernel_doubleAdd(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "mul"))) 
    ogs::gatherVecKernel_doubleMul(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "min"))) 
    ogs::gatherVecKernel_doubleMin(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "max"))) 
    ogs::gatherVecKernel_doubleMax(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "add"))) 
    ogs::gatherVecKernel_intAdd(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "mul"))) 
    ogs::gatherVecKernel_intMul(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "min"))) 
    ogs::gatherVecKernel_intMin(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "max"))) 
    ogs::gatherVecKernel_intMax(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "add"))) 
    ogs::gatherVecKernel_longAdd(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "mul"))) 
    ogs::gatherVecKernel_longMul(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "min"))) 
    ogs::gatherVecKernel_longMin(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "max"))) 
    ogs::gatherVecKernel_longMax(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v, o_gv);
}
