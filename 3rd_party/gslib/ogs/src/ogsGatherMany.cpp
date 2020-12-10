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

#include "gatherMany.tpp"

void ogsGatherMany_add(void *gv, void *v, const int k, const dlong gstride, const dlong stride, const size_t Nbytes, const char *type, ogs_t *ogs);
void ogsGatherMany_mul(void *gv, void *v, const int k, const dlong gstride, const dlong stride, const size_t Nbytes, const char *type, ogs_t *ogs);
void ogsGatherMany_min(void *gv, void *v, const int k, const dlong gstride, const dlong stride, const size_t Nbytes, const char *type, ogs_t *ogs);
void ogsGatherMany_max(void *gv, void *v, const int k, const dlong gstride, const dlong stride, const size_t Nbytes, const char *type, ogs_t *ogs);

void ogsGatherMany(occa::memory o_gv, 
               occa::memory o_v,
               const int k,
               const dlong gstride,
               const dlong stride, 
               const char *type, 
               const char *op, 
               ogs_t *ogs){
  ogsGatherManyStart (o_gv, o_v, k, gstride, stride, type, op, ogs);
  ogsGatherManyFinish(o_gv, o_v, k, gstride, stride, type, op, ogs);
}

void ogsGatherManyStart(occa::memory o_gv, 
                    occa::memory o_v, 
                    const int k,
                    const dlong gstride,
                    const dlong stride, 
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
    occaGatherMany(ogs->NhaloGather, k, stride, ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, o_v, ogs::o_haloBuf);
    
    ogs->device.finish();
    ogs->device.setStream(ogs::dataStream);
    ogs::o_haloBuf.copyTo(ogs::haloBuf, ogs->NhaloGather*Nbytes*k, 0, "async: true");
    ogs->device.setStream(ogs::defaultStream);
  }
}


void ogsGatherManyFinish(occa::memory o_gv, 
                     occa::memory o_v, 
                     const int k,
                     const dlong gstride,
                     const dlong stride, 
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
    occaGatherMany(ogs->NlocalGather, k, stride, gstride, ogs->o_localGatherOffsets, ogs->o_localGatherIds, type, op, o_v, o_gv);
  }

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);
    ogs->device.finish();

    void* H[k];
    for (int i=0;i<k;i++) H[i] = (char*)ogs::haloBuf + i*ogs->NhaloGather*Nbytes;

    // MPI based gather using libgs
    ogsHostGatherMany(H, k, type, op, ogs->haloGshNonSym);

    // copy totally gather halo data back from HOST to DEVICE
    if (ogs->NownedHalo)
      for (int i=0;i<k;i++)
        o_gv.copyFrom((char*)ogs::haloBuf+ogs->NhaloGather*Nbytes*i, 
                      ogs->NownedHalo*Nbytes, 
                      ogs->NlocalGather*Nbytes*i, "async: true");

    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);
  }
}

void ogsGatherMany(void *gv, 
               void *v,
               const int k,
               const dlong gstride,
               const dlong stride,  
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
    ogsGatherMany_add(gv, v, k, gstride, stride, Nbytes, type, ogs);
  else if (!strcmp(op, "mul")) 
    ogsGatherMany_mul(gv, v, k, gstride, stride, Nbytes, type, ogs);
  else if (!strcmp(op, "min")) 
    ogsGatherMany_min(gv, v, k, gstride, stride, Nbytes, type, ogs);
  else if (!strcmp(op, "max")) 
    ogsGatherMany_max(gv, v, k, gstride, stride, Nbytes, type, ogs);
}

void ogsGatherMany_add(void *gv, void *v, const int k, const dlong gstride, const dlong stride, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gatherMany_add<float>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gatherMany_add<double>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gatherMany_add<int>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gatherMany_add<long long int>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    void* H[k];
    for (int i=0;i<k;i++) H[i] = (char*)ogs::hostBuf + i*ogs->NhaloGather*Nbytes;

    ogsHostGatherMany(H, k, type, ogsAdd, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      for (int i=0;i<k;i++)
        memcpy((char*)gv+ogs->NlocalGather*Nbytes*i, 
               (char*)ogs::hostBuf+ogs->NhaloGather*Nbytes*i, 
               ogs->NownedHalo*Nbytes);
  }

  if (!strcmp(type, "float")) 
    gatherMany_add<float>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gatherMany_add<double>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gatherMany_add<int>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gatherMany_add<long long int>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void ogsGatherMany_mul(void *gv, void *v, const int k, const dlong gstride, const dlong stride, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gatherMany_mul<float>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gatherMany_mul<double>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gatherMany_mul<int>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gatherMany_mul<long long int>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    void* H[k];
    for (int i=0;i<k;i++) H[i] = (char*)ogs::hostBuf + i*ogs->NhaloGather*Nbytes;

    ogsHostGatherMany(H, k, type, ogsAdd, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      for (int i=0;i<k;i++)
        memcpy((char*)gv+ogs->NlocalGather*Nbytes*i, 
               (char*)ogs::hostBuf+ogs->NhaloGather*Nbytes*i, 
               ogs->NownedHalo*Nbytes);
  }

  if (!strcmp(type, "float")) 
    gatherMany_mul<float>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gatherMany_mul<double>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gatherMany_mul<int>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gatherMany_mul<long long int>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void ogsGatherMany_min(void *gv, void *v, const int k, const dlong gstride, const dlong stride, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gatherMany_min<float>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gatherMany_min<double>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gatherMany_min<int>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gatherMany_min<long long int>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    void* H[k];
    for (int i=0;i<k;i++) H[i] = (char*)ogs::hostBuf + i*ogs->NhaloGather*Nbytes;

    ogsHostGatherMany(H, k, type, ogsAdd, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      for (int i=0;i<k;i++)
        memcpy((char*)gv+ogs->NlocalGather*Nbytes*i, 
               (char*)ogs::hostBuf+ogs->NhaloGather*Nbytes*i, 
               ogs->NownedHalo*Nbytes);
  }

  if (!strcmp(type, "float")) 
    gatherMany_min<float>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gatherMany_min<double>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gatherMany_min<int>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gatherMany_min<long long int>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void ogsGatherMany_max(void *gv, void *v, const int k, const dlong gstride, const dlong stride, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gatherMany_max<float>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gatherMany_max<double>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gatherMany_max<int>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gatherMany_max<long long int>(ogs->NhaloGather, k, stride, ogs->NhaloGather, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    void* H[k];
    for (int i=0;i<k;i++) H[i] = (char*)ogs::hostBuf + i*ogs->NhaloGather*Nbytes;

    ogsHostGatherMany(H, k, type, ogsAdd, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      for (int i=0;i<k;i++)
        memcpy((char*)gv+ogs->NlocalGather*Nbytes*i, 
               (char*)ogs::hostBuf+ogs->NhaloGather*Nbytes*i, 
               ogs->NownedHalo*Nbytes);
  }

  if (!strcmp(type, "float")) 
    gatherMany_max<float>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gatherMany_max<double>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gatherMany_max<int>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gatherMany_max<long long int>(ogs->NlocalGather, k, stride, gstride, 
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void occaGatherMany(const  dlong Ngather,
                const int Nentries,
                const dlong stride,
                const dlong gstride,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_gv) {
  
  if      ((!strcmp(type, "float"))&&(!strcmp(op, "add"))) 
    ogs::gatherManyKernel_floatAdd(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "mul"))) 
    ogs::gatherManyKernel_floatMul(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "min"))) 
    ogs::gatherManyKernel_floatMin(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "max"))) 
    ogs::gatherManyKernel_floatMax(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "add"))) 
    ogs::gatherManyKernel_doubleAdd(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "mul"))) 
    ogs::gatherManyKernel_doubleMul(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "min"))) 
    ogs::gatherManyKernel_doubleMin(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "max"))) 
    ogs::gatherManyKernel_doubleMax(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "add"))) 
    ogs::gatherManyKernel_intAdd(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "mul"))) 
    ogs::gatherManyKernel_intMul(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "min"))) 
    ogs::gatherManyKernel_intMin(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "max"))) 
    ogs::gatherManyKernel_intMax(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "add"))) 
    ogs::gatherManyKernel_longAdd(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "mul"))) 
    ogs::gatherManyKernel_longMul(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "min"))) 
    ogs::gatherManyKernel_longMin(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "max"))) 
    ogs::gatherManyKernel_longMax(Ngather, Nentries, stride, gstride, o_gatherStarts, o_gatherIds, o_v, o_gv);
}
