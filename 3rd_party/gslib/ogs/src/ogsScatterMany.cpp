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

#include "scatterMany.tpp"

void ogsScatterMany_op(void *gv, void *v, 
                        const int k, 
                        const dlong sstride,  
                        const dlong stride,  
                        const size_t Nbytes, 
                        const char *type, 
                        ogs_t *ogs);

void ogsScatterMany(occa::memory o_sv, 
               occa::memory o_v, 
               const int k,
               const dlong sstride,
               const dlong stride,
               const char *type, 
               const char *op, 
               ogs_t *ogs){
  ogsScatterManyStart (o_sv, o_v, k, sstride, stride, type, op, ogs);
  ogsScatterManyFinish(o_sv, o_v, k, sstride, stride, type, op, ogs);
}

void ogsScatterManyStart(occa::memory o_sv, 
                    occa::memory o_v, 
                    const int k,
                    const dlong sstride,
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

  if(ogs->NlocalGather) {
    occaScatterMany(ogs->NlocalGather, k, stride, sstride, ogs->o_localGatherOffsets, ogs->o_localGatherIds, type, op, o_v, o_sv);
  }

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);

    if (ogs->NownedHalo) {
      for (int i=0;i<k;i++) 
        o_v.copyTo((char*)ogs::haloBuf+ogs->NhaloGather*Nbytes*i, 
                    ogs->NownedHalo*Nbytes, ogs->NlocalGather*Nbytes*i, "async: true");
    }

    ogs->device.setStream(ogs::defaultStream);
  }
}


void ogsScatterManyFinish(occa::memory o_sv, 
                     occa::memory o_v, 
                     const int k,
                     const dlong sstride,
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
    ogs->device.setStream(ogs::dataStream);
    ogs->device.finish();

    void* H[k];
    for (int i=0;i<k;i++) H[i] = (char*)ogs::haloBuf + i*ogs->NhaloGather*Nbytes;

    // MPI based scatter using gslib
    ogsHostScatterMany(H, k, type, op, ogs->haloGshNonSym);

    // copy totally scattered halo data back from HOST to DEVICE
    ogs::o_haloBuf.copyFrom(ogs::haloBuf, ogs->NhaloGather*Nbytes*k, 0, "async: true");

    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);

    occaScatterMany(ogs->NhaloGather, k, ogs->NhaloGather, sstride, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, ogs::o_haloBuf, o_sv);
  }
}

void ogsScatterMany(void *sv, 
               void *v, 
               const int k,
               const dlong sstride,
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

  ogsScatterMany_op(sv, v, k, sstride, stride, Nbytes, type, ogs);
}

void ogsScatterMany_op(void *sv, void *v, 
                       const int k,
                       const dlong sstride,
                       const dlong stride,
                       const size_t Nbytes, 
                       const char *type, 
                       ogs_t *ogs){

  if (!strcmp(type, "float")) 
    scatterMany<float>(ogs->NlocalGather,  k, stride, sstride,
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)sv);
  else if (!strcmp(type, "double")) 
    scatterMany<double>(ogs->NlocalGather, k, stride, sstride,
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)sv);
  else if (!strcmp(type, "int")) 
    scatterMany<int>(ogs->NlocalGather, k, stride, sstride, 
                      ogs->localGatherOffsets, 
                      ogs->localGatherIds, (int*)v, (int*)sv);
  else if (!strcmp(type, "long long int")) 
    scatterMany<long long int>(ogs->NlocalGather, k, stride, sstride,
                      ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)sv);

  if (ogs->NhaloGather) {
    if (ogs->NownedHalo)
      for (int i=0;i<k;i++)
        memcpy((char*)ogs::hostBuf + ogs->NhaloGather*Nbytes*i, 
               (char*) v+ogs->NlocalGather*Nbytes*i, 
               ogs->NownedHalo*Nbytes);

    
    void* H[k];
    for (int i=0;i<k;i++) H[i] = (char*)ogs::hostBuf + i*ogs->NhaloGather*Nbytes;

    // MPI based scatter using gslib
    ogsHostScatterMany(H, k, type, ogsAdd, ogs->haloGshNonSym);
  }

  if (!strcmp(type, "float")) 
    scatterMany<float>(ogs->NhaloGather, k, ogs->NhaloGather, sstride,
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)ogs::hostBuf, (float*)sv);
  else if (!strcmp(type, "double")) 
    scatterMany<double>(ogs->NhaloGather, k, ogs->NhaloGather, sstride, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)ogs::hostBuf, (double*)sv);
  else if (!strcmp(type, "int")) 
    scatterMany<int>(ogs->NhaloGather, k, ogs->NhaloGather, sstride, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)ogs::hostBuf, (int*)sv);
  else if (!strcmp(type, "long long int")) 
    scatterMany<long long int>(ogs->NhaloGather, k, ogs->NhaloGather, sstride, 
                      ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)ogs::hostBuf, (long long int*)sv);
}

void occaScatterMany(const  dlong Nscatter,
                const int Nentries,
                const dlong stride,
                const dlong sstride,
                occa::memory o_scatterStarts,
                occa::memory o_scatterIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_sv) {
  
  if      (!strcmp(type, "float")) 
    ogs::scatterManyKernel_float(Nscatter, Nentries, stride, sstride, o_scatterStarts, o_scatterIds, o_v, o_sv);
  else if (!strcmp(type, "double")) 
    ogs::scatterManyKernel_double(Nscatter, Nentries, stride, sstride, o_scatterStarts, o_scatterIds, o_v, o_sv);
  else if (!strcmp(type, "int")) 
    ogs::scatterManyKernel_int(Nscatter, Nentries, stride, sstride, o_scatterStarts, o_scatterIds, o_v, o_sv);
  else if (!strcmp(type, "long long int")) 
    ogs::scatterManyKernel_long(Nscatter, Nentries, stride, sstride, o_scatterStarts, o_scatterIds, o_v, o_sv);
}
