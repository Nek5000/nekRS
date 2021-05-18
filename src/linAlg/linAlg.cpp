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

#include "linAlg.hpp"
#include "platform.hpp"

linAlg_t* linAlg_t::singleton = nullptr;

linAlg_t*
linAlg_t::getInstance()
{
  if(!singleton)
    singleton = new linAlg_t();
  return singleton;
}
linAlg_t::linAlg_t() {
  blocksize = BLOCKSIZE;
  serial = platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP";
  comm = platform->comm.mpiComm;
  setup();
}
void linAlg_t::reallocScratch(const dlong Nbytes)
{
  device_t& device = platform->device;
  if(h_scratch.size()) h_scratch.free();
  if(o_scratch.size()) o_scratch.free();
  //pinned scratch buffer
  {
    occa::properties props = kernelInfo;
    props["host"] = true;
    h_scratch = device.malloc(Nbytes, props);
    scratch = (dfloat*) h_scratch.ptr();
  }
  o_scratch = device.malloc(Nbytes);
}
void linAlg_t::setup() {

  device_t& device = platform->device;
  int rank;
  MPI_Comm_rank(comm, &rank);

  occa::properties kernelInfo = platform->kernelInfo;

  reallocScratch(blocksize * sizeof(dfloat));

  string oklDir;
  oklDir.assign(getenv("NEKRS_INSTALL_DIR"));
  oklDir += "/okl/linAlg/";

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("loading linAlg kernels ... "); fflush(stdout);

  occa::properties kernelInfoNoOkl = platform->kernelInfo;
  kernelInfoNoOkl["okl/enabled"] = false;

  {
      if (fillKernel.isInitialized()==false)
        fillKernel = device.buildKernel(oklDir + 
                                        "linAlgFill.okl",
                                        "fill",
                                        kernelInfo);
      if (addKernel.isInitialized()==false)
        addKernel = device.buildKernel(oklDir + 
                                        "linAlgAdd.okl",
                                        "add",
                                        kernelInfo);
      if (scaleKernel.isInitialized()==false)
        scaleKernel = device.buildKernel(oklDir + 
                                        "linAlgScale.okl",
                                        "scale",
                                        kernelInfo);
      if (scaleManyKernel.isInitialized()==false)
        scaleManyKernel = device.buildKernel(oklDir + 
                                        "linAlgScale.okl",
                                        "scaleMany",
                                        kernelInfo);
      if (axpbyKernel.isInitialized()==false){
        if(serial){
          axpbyKernel = device.buildKernel(oklDir + 
                                           "linAlgAXPBY.c",
                                           "axpby",
                                           kernelInfoNoOkl);
        } else {
          axpbyKernel = device.buildKernel(oklDir + 
                                           "linAlgAXPBY.okl",
                                           "axpby",
                                           kernelInfo);
        }
      }
      if (axpbyManyKernel.isInitialized()==false){
        if(serial){
          axpbyManyKernel = device.buildKernel(oklDir + 
                                           "linAlgAXPBY.c",
                                           "axpbyMany",
                                           kernelInfoNoOkl);
        } else {
          axpbyManyKernel = device.buildKernel(oklDir + 
                                           "linAlgAXPBY.okl",
                                           "axpbyMany",
                                           kernelInfo);
        }
      }
      if (axpbyzKernel.isInitialized()==false)
        axpbyzKernel = device.buildKernel(oklDir + 
                                          "linAlgAXPBY.okl",
                                          "axpbyz",
                                          kernelInfo);
      if (axpbyzManyKernel.isInitialized()==false)
        axpbyzManyKernel = device.buildKernel(oklDir + 
                                          "linAlgAXPBY.okl",
                                          "axpbyzMany",
                                          kernelInfo);
      if (axmyKernel.isInitialized()==false){
        if(serial){
          axmyKernel = device.buildKernel(oklDir + 
                                          "linAlgAXMY.c",
                                          "axmy",
                                          kernelInfoNoOkl);
        } else {
          axmyKernel = device.buildKernel(oklDir + 
                                          "linAlgAXMY.okl",
                                          "axmy",
                                          kernelInfo);
        }
      }
      if (axmyManyKernel.isInitialized()==false){
        if(serial){
          axmyManyKernel = device.buildKernel(oklDir + 
                                          "linAlgAXMY.c",
                                          "axmyMany",
                                          kernelInfoNoOkl);
        } else {
          axmyManyKernel = device.buildKernel(oklDir + 
                                          "linAlgAXMY.okl",
                                          "axmyMany",
                                          kernelInfo);
        }
      }
      if (axmyVectorKernel.isInitialized()==false){
        if(serial){
          axmyVectorKernel = device.buildKernel(oklDir + 
                                          "linAlgAXMY.c",
                                          "axmyVector",
                                          kernelInfoNoOkl);
        } else {
          axmyVectorKernel = device.buildKernel(oklDir + 
                                          "linAlgAXMY.okl",
                                          "axmyVector",
                                          kernelInfo);
        }
      }
      if (axmyzKernel.isInitialized()==false)
        axmyzKernel = device.buildKernel(oklDir + 
                                         "linAlgAXMY.okl",
                                         "axmyz",
                                         kernelInfo);
      if (axmyzManyKernel.isInitialized()==false)
        axmyzManyKernel = device.buildKernel(oklDir + 
                                         "linAlgAXMY.okl",
                                         "axmyzMany",
                                         kernelInfo);
      if (adyKernel.isInitialized()==false)
        adyKernel = device.buildKernel(oklDir + 
                                        "linAlgAXDY.okl",
                                        "ady",
                                        kernelInfo);
      if (adyManyKernel.isInitialized()==false)
        adyManyKernel = device.buildKernel(oklDir + 
                                        "linAlgAXDY.okl",
                                        "adyMany",
                                        kernelInfo);
      if (axdyKernel.isInitialized()==false)
        axdyKernel = device.buildKernel(oklDir + 
                                        "linAlgAXDY.okl",
                                        "axdy",
                                        kernelInfo);
      if (aydxKernel.isInitialized()==false)
        aydxKernel = device.buildKernel(oklDir + 
                                        "linAlgAXDY.okl",
                                        "aydx",
                                        kernelInfo);
      if (aydxManyKernel.isInitialized()==false)
        aydxManyKernel = device.buildKernel(oklDir + 
                                        "linAlgAXDY.okl",
                                        "aydxMany",
                                        kernelInfo);
      if (axmyzKernel.isInitialized()==false)
        axmyzKernel = device.buildKernel(oklDir + 
                                         "linAlgAXDY.okl",
                                         "axdyz",
                                         kernelInfo);
      if (sumKernel.isInitialized()==false)
        sumKernel = device.buildKernel(oklDir + 
                                        "linAlgSum.okl",
                                        "sum",
                                        kernelInfo);
      if (minKernel.isInitialized()==false)
        minKernel = device.buildKernel(oklDir + 
                                        "linAlgMin.okl",
                                        "min",
                                        kernelInfo);
       if (maxKernel.isInitialized()==false)
        maxKernel = device.buildKernel(oklDir + 
                                        "linAlgMax.okl",
                                        "max",
                                        kernelInfo);
      if (norm2Kernel.isInitialized()==false)
        norm2Kernel = device.buildKernel(oklDir + 
                                        "linAlgNorm2.okl",
                                        "norm2",
                                        kernelInfo);
      if (weightedNorm2Kernel.isInitialized()==false){
        if(serial){
          weightedNorm2Kernel = device.buildKernel(oklDir + 
                                          "linAlgWeightedNorm2.c",
                                          "weightedNorm2",
                                          kernelInfoNoOkl);
        } else {
          weightedNorm2Kernel = device.buildKernel(oklDir + 
                                          "linAlgWeightedNorm2.okl",
                                          "weightedNorm2",
                                          kernelInfo);
        }
      }
      if (weightedNorm2ManyKernel.isInitialized()==false){
        if(serial){
          weightedNorm2ManyKernel = device.buildKernel(oklDir + 
                                          "linAlgWeightedNorm2.c",
                                          "weightedNorm2Many",
                                          kernelInfoNoOkl);
        } else {
          weightedNorm2ManyKernel = device.buildKernel(oklDir + 
                                          "linAlgWeightedNorm2.okl",
                                          "weightedNorm2Many",
                                          kernelInfo);
        }
      }
      if (innerProdKernel.isInitialized()==false)
        innerProdKernel = device.buildKernel(oklDir + 
                                        "linAlgInnerProd.okl",
                                        "innerProd",
                                        kernelInfo);
      if (weightedInnerProdKernel.isInitialized()==false){
        if(serial){
          weightedInnerProdKernel = device.buildKernel(oklDir + 
                                          "linAlgWeightedInnerProd.c",
                                          "weightedInnerProd",
                                          kernelInfoNoOkl);
        } else {
          weightedInnerProdKernel = device.buildKernel(oklDir + 
                                          "linAlgWeightedInnerProd.okl",
                                          "weightedInnerProd",
                                          kernelInfo);
        }
      }
      if (weightedInnerProdManyKernel.isInitialized()==false){
        if(serial){
          weightedInnerProdManyKernel = device.buildKernel(oklDir + 
                                          "linAlgWeightedInnerProd.c",
                                          "weightedInnerProdMany",
                                          kernelInfoNoOkl);
        } else {
          weightedInnerProdManyKernel = device.buildKernel(oklDir + 
                                          "linAlgWeightedInnerProd.okl",
                                          "weightedInnerProdMany",
                                          kernelInfo);
        }
      }
      if (weightedInnerProdMultiKernel.isInitialized()==false)
        weightedInnerProdMultiKernel = device.buildKernel(oklDir + 
                                        "linAlgWeightedInnerProd.okl",
                                        "weightedInnerProdMulti",
                                        kernelInfo);
  }

  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel); fflush(stdout);
}

linAlg_t::~linAlg_t() {
  fillKernel.free();
  addKernel.free();
  scaleKernel.free();
  scaleManyKernel.free();
  axpbyKernel.free();
  axpbyManyKernel.free();
  axpbyzKernel.free();
  axpbyzManyKernel.free();
  axmyKernel.free();
  axmyManyKernel.free();
  axmyVectorKernel.free();
  axmyzKernel.free();
  axmyzManyKernel.free();
  axdyKernel.free();
  aydxKernel.free();
  aydxManyKernel.free();
  adyKernel.free();
  adyManyKernel.free();
  axdyzKernel.free();
  sumKernel.free();
  minKernel.free();
  maxKernel.free();
  norm2Kernel.free();
  weightedNorm2Kernel.free();
  weightedNorm2ManyKernel.free();
  innerProdKernel.free();
  weightedInnerProdKernel.free();
  weightedInnerProdManyKernel.free();
  weightedInnerProdMultiKernel.free();
}

/*********************/
/* vector operations */
/*********************/

// o_a[n] = alpha
void linAlg_t::fill(const dlong N, const dfloat alpha, occa::memory& o_a) {
  fillKernel(N, alpha, o_a);
}

// o_a[n] += alpha
void linAlg_t::add(const dlong N, const dfloat alpha, occa::memory& o_a, const dlong offset) {
  addKernel(N, offset, alpha, o_a);
}

// o_a[n] *= alpha
void linAlg_t::scale(const dlong N, const dfloat alpha, occa::memory& o_a)  {
  scaleKernel(N, alpha, o_a);
}
void linAlg_t::scaleMany(const dlong N, const dlong Nfields, const dlong fieldOffset, const dfloat alpha, occa::memory& o_a, const dlong offset)  {
  scaleManyKernel(N, Nfields, fieldOffset, offset, alpha, o_a);
}

// o_y[n] = beta*o_y[n] + alpha*o_x[n]
void linAlg_t::axpby(const dlong N, const dfloat alpha, occa::memory& o_x,
                    const dfloat beta,  occa::memory& o_y, const dlong xOffset, const dlong yOffset) {
  axpbyKernel(N, xOffset, yOffset, alpha, o_x, beta, o_y);
}
void linAlg_t::axpbyMany(const dlong N, const dlong Nfields, const dlong offset, const dfloat alpha, occa::memory& o_x,
                    const dfloat beta,  occa::memory& o_y) {
  axpbyManyKernel(N, Nfields, offset, alpha, o_x, beta, o_y);
}

// o_z[n] = beta*o_y[n] + alpha*o_x[n]
void linAlg_t::axpbyz(const dlong N, const dfloat alpha, occa::memory& o_x,
                      const dfloat beta, occa::memory& o_y, occa::memory& o_z) {
  axpbyzKernel(N, alpha, o_x, beta, o_y, o_z);
}
void linAlg_t::axpbyzMany(const dlong N, const dlong Nfields, const dlong fieldOffset, const dfloat alpha, occa::memory& o_x,
                      const dfloat beta, occa::memory& o_y, occa::memory& o_z) {
  axpbyzManyKernel(N, Nfields, fieldOffset, alpha, o_x, beta, o_y, o_z);
}

// o_y[n] = alpha*o_x[n]*o_y[n]
void linAlg_t::axmy(const dlong N, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y) {
  axmyKernel(N, alpha, o_x, o_y);
}
void linAlg_t::axmyMany(const dlong N, const dlong Nfields, const dlong offset,
                   const dlong mode, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y) {
  axmyManyKernel(N, Nfields, offset, mode, alpha, o_x, o_y);
}
void linAlg_t::axmyVector(const dlong N, const dlong offset,
                   const dlong mode, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y) {
  axmyVectorKernel(N, offset, mode, alpha, o_x, o_y);
}

// o_z[n] = alpha*o_x[n]*o_y[n]
void linAlg_t::axmyz(const dlong N, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y, occa::memory& o_z) {
  axmyzKernel(N, alpha, o_x, o_y, o_z);
}
void linAlg_t::axmyzMany(const dlong N, const dlong Nfields, const dlong offset, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y, occa::memory& o_z) {
  axmyzManyKernel(N, Nfields, offset, alpha, o_x, o_y, o_z);
}

// o_y[n] = alpha*o_x[n]/o_y[n]
void linAlg_t::axdy(const dlong N, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y) {
  axdyKernel(N, alpha, o_x, o_y);
}
void linAlg_t::aydx(const dlong N, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y) {
  aydxKernel(N, alpha, o_x, o_y);
}
void linAlg_t::aydxMany(const dlong N, const dlong Nfields, const dlong fieldOffset,
                   const dlong mode, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y) {
  aydxManyKernel(N, Nfields, fieldOffset, mode, alpha, o_x, o_y);
}
// o_y[n] = alpha/o_y[n]
void linAlg_t::ady(const dlong N, const dfloat alpha,
                   occa::memory& o_y) {
  adyKernel(N, alpha, o_y);
}
void linAlg_t::adyMany(const dlong N, const dlong Nfields, const dlong offset, const dfloat alpha,
                   occa::memory& o_y) {
  adyManyKernel(N, Nfields, offset, alpha, o_y);
}

// o_z[n] = alpha*o_x[n]/o_y[n]
void linAlg_t::axdyz(const dlong N, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y, occa::memory& o_z) {
  axdyzKernel(N, alpha, o_x, o_y, o_z);
}

// \sum o_a
dfloat linAlg_t::sum(const dlong N, occa::memory& o_a, MPI_Comm _comm, const dlong offset) {
  int Nblock = (N+blocksize-1)/blocksize;
  const dlong Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  sumKernel(Nblock, N, offset, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nbytes);

  dfloat sum = 0;
  for(dlong n=0;n<Nblock;++n){
    sum += scratch[n];
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DFLOAT, MPI_SUM, _comm);

  return sum;
}

// \min o_a
dfloat linAlg_t::min(const dlong N, occa::memory& o_a, MPI_Comm _comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  const dlong Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  minKernel(Nblock, N, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nbytes);

  dfloat min = scratch[0];
  for(dlong n=1;n<Nblock;++n){
    min = (scratch[n] < min) ? scratch[n]:min;
  }

  MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DFLOAT, MPI_MIN, _comm);

  return min;
}


// \max o_a
dfloat linAlg_t::max(const dlong N, occa::memory& o_a, MPI_Comm _comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  const dlong Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  maxKernel(Nblock, N, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nbytes);

  dfloat max = scratch[0];
  for(dlong n=1;n<Nblock;++n){
    max = (scratch[n] > max) ? scratch[n]:max;
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DFLOAT, MPI_MAX, _comm);

  return max;
}

// ||o_a||_2
/*
dfloat linAlg_t::norm2(const dlong N, occa::memory& o_a, MPI_Comm _comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  norm2Kernel(Nblock, N, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(dfloat));

  dfloat norm = 0;
  for(dlong n=0;n<Nblock;++n){
    norm += scratch[n];
  }

  dfloat globalnorm = 0;
  MPI_Allreduce(&norm, &globalnorm, 1, MPI_DFLOAT, MPI_SUM, _comm);

  return sqrt(globalnorm);
}
*/

// o_x.o_y
dfloat linAlg_t::innerProd(const dlong N, occa::memory& o_x, occa::memory& o_y,
                           MPI_Comm _comm, const dlong offset) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const dlong Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  innerProdKernel(Nblock, N, offset, o_x, o_y, o_scratch);

  o_scratch.copyTo(scratch, Nbytes);

  dfloat dot = 0;
  for(dlong n=0;n<Nblock;++n){
    dot += scratch[n];
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_DFLOAT, MPI_SUM, _comm);

#ifdef ENABLE_TIMER
  platform->timer.toc("dotp");
#endif
  return dot;
}

// o_w.o_x.o_y
dfloat linAlg_t::weightedInnerProd(const dlong N, occa::memory& o_w,
                                   occa::memory& o_x, occa::memory& o_y,
                                   MPI_Comm _comm) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const dlong Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  weightedInnerProdKernel(Nblock, N, o_w, o_x, o_y, o_scratch);


  dfloat dot = 0;

  if(serial){
    dot = *((dfloat*) o_scratch.ptr());
  } else {
    o_scratch.copyTo(scratch, Nbytes);
    for(dlong n=0;n<Nblock;++n){
      dot += scratch[n];
    }
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_DFLOAT, MPI_SUM, _comm);

#ifdef ENABLE_TIMER
  platform->timer.toc("dotp");
#endif
  return dot;
}
void linAlg_t::weightedInnerProdMulti(const dlong N, 
                                   const dlong NVec,
                                   const dlong Nfields,
                                   const dlong fieldOffset,
                                   occa::memory& o_w,
                                   occa::memory& o_x, occa::memory& o_y,
                                   MPI_Comm _comm, dfloat* result, const dlong offset) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const dlong Nbytes = NVec * Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  weightedInnerProdMultiKernel(Nblock, N, Nfields, fieldOffset, NVec, offset, o_w, o_x, o_y, o_scratch);

  o_scratch.copyTo(scratch, Nbytes);

  for(int field = 0; field < NVec; ++field){
    dfloat dot = 0;
    for(dlong n=0;n<Nblock;++n){
      dot += scratch[n + field * Nblock];
    }
    result[field] = dot;
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, result, NVec, MPI_DFLOAT, MPI_SUM, _comm);
#ifdef ENABLE_TIMER
  platform->timer.toc("dotp");
#endif
}
dfloat linAlg_t::weightedInnerProdMany(const dlong N, 
                                   const dlong Nfields,
                                   const dlong fieldOffset,
                                   occa::memory& o_w,
                                   occa::memory& o_x, occa::memory& o_y,
                                   MPI_Comm _comm) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const dlong Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  weightedInnerProdManyKernel(Nblock, N, Nfields, fieldOffset, o_w, o_x, o_y, o_scratch);

  dfloat dot = 0;

  if(serial){
    dot = *((dfloat*) o_scratch.ptr());
  } else {
    o_scratch.copyTo(scratch, Nbytes);
    for(dlong n=0;n<Nblock;++n){
      dot += scratch[n];
    }
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_DFLOAT, MPI_SUM, _comm);
#ifdef ENABLE_TIMER
  platform->timer.toc("dotp");
#endif

  return dot;
}

// ||o_a||_w2
dfloat linAlg_t::weightedNorm2(const dlong N, occa::memory& o_w,
                               occa::memory& o_a, MPI_Comm _comm) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const dlong Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  weightedNorm2Kernel(Nblock, N, o_w, o_a, o_scratch);



  dfloat norm = 0;
  if(serial){
    norm = *((dfloat*) o_scratch.ptr());
  } else {
    o_scratch.copyTo(scratch, Nbytes);
    for(dlong n=0;n<Nblock;++n){
      norm += scratch[n];
    }
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);
#ifdef ENABLE_TIMER
  platform->timer.toc("dotp");
#endif

  return sqrt(norm);
}
dfloat linAlg_t::weightedNorm2Many(const dlong N,
                                   const dlong Nfields,
                                   const dlong fieldOffset,
                                   occa::memory& o_w,
                               occa::memory& o_a, MPI_Comm _comm) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const dlong Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  weightedNorm2ManyKernel(Nblock, N, Nfields, fieldOffset, o_w, o_a, o_scratch);


  dfloat norm = 0;
  if(serial){
    norm = *((dfloat*) o_scratch.ptr());
  } else {
    o_scratch.copyTo(scratch, Nbytes);
    for(dlong n=0;n<Nblock;++n){
      norm += scratch[n];
    }
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);

#ifdef ENABLE_TIMER
  platform->timer.toc("dotp");
#endif
  return sqrt(norm);
}
