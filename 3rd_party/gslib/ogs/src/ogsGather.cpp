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

#include "gather.tpp"

static double startTime = 0.0;
static double stopTime = 0.0;
static double hostTime = 0.0;
static double deviceTime = 0.0;
occa::streamTag startTag;

void ogsHostTic(MPI_Comm comm, int ifSync){
#ifdef OGS_ENABLE_TIMER
  if(ifSync) MPI_Barrier(comm);
  startTime = MPI_Wtime();
#endif
}

void ogsHostToc(){
#ifdef OGS_ENABLE_TIMER
  stopTime = MPI_Wtime();
  hostTime += (stopTime - startTime);
#endif
}

double ogsTime(bool reportHostTime){
  if(reportHostTime)
    return hostTime;
  else
    return deviceTime;
}

void ogsResetTime(){
  deviceTime = 0.0;
  hostTime = 0.0;
}

void ogsGather_add(void *gv, void *v, const size_t Nbytes, const char *type, ogs_t *ogs);
void ogsGather_mul(void *gv, void *v, const size_t Nbytes, const char *type, ogs_t *ogs);
void ogsGather_min(void *gv, void *v, const size_t Nbytes, const char *type, ogs_t *ogs);
void ogsGather_max(void *gv, void *v, const size_t Nbytes, const char *type, ogs_t *ogs);

void ogsGather(occa::memory o_gv, 
               occa::memory o_v, 
               const char *type, 
               const char *op, 
               ogs_t *ogs){
  ogsGatherStart (o_gv, o_v, type, op, ogs);
  ogsGatherFinish(o_gv, o_v, type, op, ogs);
}

void ogsGatherStart(occa::memory o_gv, 
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

  // gather halo nodes on device
  if (ogs->NhaloGather) {
    occaGather(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, o_v, ogs::o_haloBuf);
    
    ogs->device.finish();
    ogs->device.setStream(ogs::dataStream);
    ogs::o_haloBuf.copyTo(ogs::haloBuf, ogs->NhaloGather*Nbytes, 0, "async: true");
    ogs->device.setStream(ogs::defaultStream);
  }
}


void ogsGatherFinish(occa::memory o_gv, 
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

  if(ogs->NlocalGather) {
    occaGather(ogs->NlocalGather, ogs->o_localGatherOffsets, ogs->o_localGatherIds, type, op, o_v, o_gv);
  }

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);
    ogs->device.finish();

#ifdef OGS_ENABLE_TIMER
  ogsHostTic(ogs->comm, 1);
#endif
    // MPI based gather using libgs
    ogsHostGather(ogs::haloBuf, type, op, ogs->haloGshNonSym);
#ifdef OGS_ENABLE_TIMER
  ogsHostToc();
#endif

    // copy totally gather halo data back from HOST to DEVICE
    if (ogs->NownedHalo)
      o_gv.copyFrom(ogs::haloBuf, ogs->NownedHalo*Nbytes, 
                              ogs->NlocalGather*Nbytes, "async: true");

    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);
  }
}

void ogsGather(void *gv, 
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

  if (!strcmp(op, "add")) 
    ogsGather_add(gv, v, Nbytes, type, ogs);
  else if (!strcmp(op, "mul")) 
    ogsGather_mul(gv, v, Nbytes, type, ogs);
  else if (!strcmp(op, "min")) 
    ogsGather_min(gv, v, Nbytes, type, ogs);
  else if (!strcmp(op, "max")) 
    ogsGather_max(gv, v, Nbytes, type, ogs);
}

void ogsGather_add(void *gv, void *v, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gather_add<float>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gather_add<double>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gather_add<int>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gather_add<long long int>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    // MPI based scatter using gslib
    ogsHostGather(ogs::hostBuf, type, ogsAdd, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      memcpy((char*)gv+ogs->NlocalGather*Nbytes, ogs::hostBuf, ogs->NownedHalo*Nbytes);
  }

  if (!strcmp(type, "float")) 
    gather_add<float>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gather_add<double>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gather_add<int>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gather_add<long long int>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void ogsGather_mul(void *gv, void *v, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gather_mul<float>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gather_mul<double>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gather_mul<int>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gather_mul<long long int>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    // MPI based scatter using gslib
    ogsHostGather(ogs::hostBuf, type, ogsMul, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      memcpy((char*)gv+ogs->NlocalGather*Nbytes, ogs::hostBuf, ogs->NownedHalo*Nbytes);
  }

  if (!strcmp(type, "float")) 
    gather_mul<float>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gather_mul<double>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gather_mul<int>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gather_mul<long long int>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void ogsGather_min(void *gv, void *v, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gather_min<float>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gather_min<double>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gather_min<int>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gather_min<long long int>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    // MPI based scatter using gslib
    ogsHostGather(ogs::hostBuf, type, ogsMin, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      memcpy((char*)gv+ogs->NlocalGather*Nbytes, ogs::hostBuf, ogs->NownedHalo*Nbytes);
  }

  if (!strcmp(type, "float")) 
    gather_min<float>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gather_min<double>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gather_min<int>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gather_min<long long int>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void ogsGather_max(void *gv, void *v, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    gather_max<float>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)v, (float*)ogs::hostBuf);
  else if (!strcmp(type, "double")) 
    gather_max<double>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)v, (double*)ogs::hostBuf);
  else if (!strcmp(type, "int")) 
    gather_max<int>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)v, (int*)ogs::hostBuf);
  else if (!strcmp(type, "long long int")) 
    gather_max<long long int>(ogs->NhaloGather, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)v, (long long int*)ogs::hostBuf);

  if (ogs->NhaloGather) {
    // MPI based scatter using gslib
    ogsHostGather(ogs::hostBuf, type, ogsMax, ogs->haloGshNonSym);
    
    if (ogs->NownedHalo)
      memcpy((char*)gv+ogs->NlocalGather*Nbytes, ogs::hostBuf, ogs->NownedHalo*Nbytes);
  }

  if (!strcmp(type, "float")) 
    gather_max<float>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)gv);
  else if (!strcmp(type, "double")) 
    gather_max<double>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)gv);
  else if (!strcmp(type, "int")) 
    gather_max<int>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)gv);
  else if (!strcmp(type, "long long int")) 
    gather_max<long long int>(ogs->NlocalGather, ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)gv);
}

void occaGather(const  dlong Ngather,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_gv) {
  
  if      ((!strcmp(type, "float"))&&(!strcmp(op, "add"))) 
    ogs::gatherKernel_floatAdd(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "mul"))) 
    ogs::gatherKernel_floatMul(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "min"))) 
    ogs::gatherKernel_floatMin(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "max"))) 
    ogs::gatherKernel_floatMax(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "add"))) 
    ogs::gatherKernel_doubleAdd(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "mul"))) 
    ogs::gatherKernel_doubleMul(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "min"))) 
    ogs::gatherKernel_doubleMin(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "max"))) 
    ogs::gatherKernel_doubleMax(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "add"))) 
    ogs::gatherKernel_intAdd(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "mul"))) 
    ogs::gatherKernel_intMul(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "min"))) 
    ogs::gatherKernel_intMin(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "max"))) 
    ogs::gatherKernel_intMax(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "add"))) 
    ogs::gatherKernel_longAdd(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "mul"))) 
    ogs::gatherKernel_longMul(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "min"))) 
    ogs::gatherKernel_longMin(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "max"))) 
    ogs::gatherKernel_longMax(Ngather, o_gatherStarts, o_gatherIds, o_v, o_gv);
}
