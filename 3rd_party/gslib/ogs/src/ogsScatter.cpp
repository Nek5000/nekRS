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

#include "scatter.tpp"

void ogsScatter_op(void *gv, void *v, const size_t Nbytes, const char *type, ogs_t *ogs);

void ogsScatter(occa::memory o_sv, 
               occa::memory o_v, 
               const char *type, 
               const char *op, 
               ogs_t *ogs){
  ogsScatterStart (o_sv, o_v, type, op, ogs);
  ogsScatterFinish(o_sv, o_v, type, op, ogs);
}

void ogsScatterStart(occa::memory o_sv, 
                    occa::memory o_v, 
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
    if (ogs::o_haloBuf.size() < ogs->NhaloGather*Nbytes) {
      if (ogs::o_haloBuf.size()) ogs::o_haloBuf.free();
      //      ogs::o_haloBuf = ogs->device.mappedAlloc(ogs->NhaloGather*Nbytes);
      //      ogs::haloBuf = ogs::o_haloBuf.getMappedPointer();

      ogs::haloBuf = ogsHostMallocPinned(ogs->device, ogs->NhaloGather*Nbytes, NULL, ogs::o_haloBuf, ogs::h_haloBuf);
    }
  }

  if(ogs->NlocalGather) {
    occaScatter(ogs->NlocalGather, ogs->o_localGatherOffsets, ogs->o_localGatherIds, type, op, o_v, o_sv);
  }

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);

    if (ogs->NownedHalo)
      o_v.copyTo(ogs::haloBuf, ogs->NownedHalo*Nbytes, 
                              ogs->NlocalGather*Nbytes, "async: true");

    ogs->device.setStream(ogs::defaultStream);
  }
}


void ogsScatterFinish(occa::memory o_sv, 
                     occa::memory o_v, 
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

    ogsHostTic(ogs->comm, 1);
    // MPI based scatter using gslib
    ogsHostScatter(ogs::haloBuf, type, op, ogs->haloGshNonSym);
    ogsHostToc();

    // copy totally scattered halo data back from HOST to DEVICE
    ogs::o_haloBuf.copyFrom(ogs::haloBuf, ogs->NhaloGather*Nbytes, 0, "async: true");

    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);

    occaScatter(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, ogs::o_haloBuf, o_sv);
  }
}

void ogsScatter(void *sv, 
               void *v, 
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
    if (ogs::hostBufSize < ogs->NhaloGather*Nbytes) {
      if (ogs::hostBufSize) free(ogs::hostBuf);
      ogs::hostBuf = (void *) malloc(ogs->NhaloGather*Nbytes);
    }
  }

  ogsScatter_op(sv, v, Nbytes, type, ogs);
}

void ogsScatter_op(void *sv, void *v, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    scatter<float>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)sv);
  else if (!strcmp(type, "double")) 
    scatter<double>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)sv);
  else if (!strcmp(type, "int")) 
    scatter<int>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)sv);
  else if (!strcmp(type, "long long int")) 
    scatter<long long int>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)sv);

  if (ogs->NhaloGather) {
    if (ogs->NownedHalo)
      memcpy(ogs::hostBuf, (char*) v+ogs->NlocalGather*Nbytes, ogs->NownedHalo*Nbytes);

    // MPI based scatter using gslib
    ogsHostScatter(ogs::hostBuf, type, ogsAdd, ogs->haloGshNonSym);
  }

  if (!strcmp(type, "float")) 
    scatter<float>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)ogs::hostBuf, (float*)sv);
  else if (!strcmp(type, "double")) 
    scatter<double>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)ogs::hostBuf, (double*)sv);
  else if (!strcmp(type, "int")) 
    scatter<int>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)ogs::hostBuf, (int*)sv);
  else if (!strcmp(type, "long long int")) 
    scatter<long long int>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)ogs::hostBuf, (long long int*)sv);
}


void occaScatter(const  dlong Nscatter,
                occa::memory o_scatterStarts,
                occa::memory o_scatterIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_sv) {
  
  if      (!strcmp(type, "float")) 
    ogs::scatterKernel_float(Nscatter, o_scatterStarts, o_scatterIds, o_v, o_sv);
  else if (!strcmp(type, "double")) 
    ogs::scatterKernel_double(Nscatter, o_scatterStarts, o_scatterIds, o_v, o_sv);
  else if (!strcmp(type, "int")) 
    ogs::scatterKernel_int(Nscatter, o_scatterStarts, o_scatterIds, o_v, o_sv);
  else if (!strcmp(type, "long long int")) 
    ogs::scatterKernel_long(Nscatter, o_scatterStarts, o_scatterIds, o_v, o_sv);
}
