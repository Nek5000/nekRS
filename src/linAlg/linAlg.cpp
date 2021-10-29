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
void linAlg_t::reallocScratch(const size_t Nbytes)
{
  device_t& device = platform->device;
  if(h_scratch.size()) h_scratch.free();
  if(o_scratch.size()) o_scratch.free();
  //pinned scratch buffer
  {
    h_scratch = device.mallocHost(Nbytes);
    scratch = (dfloat*) h_scratch.ptr();
  }
  o_scratch = device.malloc(Nbytes);
}
void linAlg_t::setup() {

  auto& kernels = platform->kernels;
  int rank;
  MPI_Comm_rank(comm, &rank);

  occa::properties kernelInfo = platform->kernelInfo;

  reallocScratch(blocksize * sizeof(dfloat));

  std::string oklDir;
  oklDir.assign(getenv("NEKRS_INSTALL_DIR"));
  oklDir += "/okl/linAlg/";

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  {
    fillKernel = kernels.getKernel("fill");
    absKernel = kernels.getKernel("vabs");
    addKernel = kernels.getKernel("add");
    scaleKernel = kernels.getKernel("scale");
    scaleManyKernel = kernels.getKernel("scaleMany");
    axpbyKernel = kernels.getKernel("axpby");
    axpbyManyKernel = kernels.getKernel("axpbyMany");
    axpbyzKernel = kernels.getKernel("axpbyz");
    axpbyzManyKernel = kernels.getKernel("axpbyzMany");
    axmyKernel = kernels.getKernel("axmy");
    axmyManyKernel = kernels.getKernel("axmyMany");
    axmyVectorKernel = kernels.getKernel("axmyVector");
    axmyzKernel = kernels.getKernel("axmyz");
    axmyzManyKernel = kernels.getKernel("axmyzMany");
    adyKernel = kernels.getKernel("ady");
    adyManyKernel = kernels.getKernel("adyMany");
    axdyKernel = kernels.getKernel("axdy");
    aydxKernel = kernels.getKernel("aydx");
    aydxManyKernel = kernels.getKernel("aydxMany");
    axdyzKernel = kernels.getKernel("axdyz");
    sumKernel = kernels.getKernel("sum");
    sumManyKernel = kernels.getKernel("sumMany");
    minKernel = kernels.getKernel("min");
    maxKernel = kernels.getKernel("max");
    norm2Kernel = kernels.getKernel("norm2");
    norm2ManyKernel = kernels.getKernel("norm2Many");
    norm1Kernel = kernels.getKernel("norm1");
    norm1ManyKernel = kernels.getKernel("norm1Many");
    weightedNorm1Kernel = kernels.getKernel("weightedNorm1");
    weightedNorm1ManyKernel = kernels.getKernel("weightedNorm1Many");
    weightedNorm2Kernel = kernels.getKernel("weightedNorm2");
    weightedNorm2ManyKernel = kernels.getKernel("weightedNorm2Many");
    innerProdKernel = kernels.getKernel("innerProd");
    weightedInnerProdKernel = kernels.getKernel("weightedInnerProd");
    weightedInnerProdManyKernel = kernels.getKernel("weightedInnerProdMany");
    weightedInnerProdMultiKernel = kernels.getKernel("weightedInnerProdMulti");
  }
}

linAlg_t::~linAlg_t() {
  fillKernel.free();
  absKernel.free();
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
  norm2ManyKernel.free();
  norm1Kernel.free();
  norm1ManyKernel.free();
  weightedNorm1Kernel.free();
  weightedNorm1ManyKernel.free();
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

// o_a[n] = abs(o_a[n])
void linAlg_t::abs(const dlong N,  occa::memory& o_a) {
  absKernel(N, o_a);
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
  const size_t Nbytes = Nblock * sizeof(dfloat);
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
dfloat linAlg_t::sumMany(const dlong N, const dlong Nfields, const dlong fieldOffset, occa::memory& o_a, MPI_Comm _comm){
  int Nblock = (N+blocksize-1)/blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  sumManyKernel(Nblock, N, Nfields, fieldOffset, o_a, o_scratch);

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
  const size_t Nbytes = Nblock * sizeof(dfloat);
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
  const size_t Nbytes = Nblock * sizeof(dfloat);
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
dfloat linAlg_t::norm2(const dlong N, occa::memory& o_x, MPI_Comm _comm) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  norm2Kernel(Nblock, N, o_x, o_scratch);

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
dfloat linAlg_t::norm2Many(const dlong N, const dlong Nfields, const dlong fieldOffset, occa::memory& o_x, MPI_Comm _comm) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  norm2ManyKernel(Nblock, N, Nfields, fieldOffset, o_x, o_scratch);
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
// ||o_a||_1
dfloat linAlg_t::norm1(const dlong N, occa::memory& o_x, MPI_Comm _comm) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  norm1Kernel(Nblock, N, o_x, o_scratch);
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

  return norm;
}
dfloat linAlg_t::norm1Many(const dlong N, const dlong Nfields, const dlong fieldOffset, occa::memory& o_x, MPI_Comm _comm) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  norm1ManyKernel(Nblock, N, Nfields, fieldOffset, o_x, o_scratch);

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
  return norm;
}

// o_x.o_y
dfloat linAlg_t::innerProd(const dlong N, occa::memory& o_x, occa::memory& o_y,
                           MPI_Comm _comm, const dlong offset) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
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
  const size_t Nbytes = Nblock * sizeof(dfloat);
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
  const size_t Nbytes = NVec * Nblock * sizeof(dfloat);
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
  const size_t Nbytes = Nblock * sizeof(dfloat);
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
  const size_t Nbytes = Nblock * sizeof(dfloat);
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
  const size_t Nbytes = Nblock * sizeof(dfloat);
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

// ||o_a||_w1
dfloat linAlg_t::weightedNorm1(const dlong N, occa::memory& o_w,
                               occa::memory& o_a, MPI_Comm _comm) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  weightedNorm1Kernel(Nblock, N, o_w, o_a, o_scratch);



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

  return norm;
}
dfloat linAlg_t::weightedNorm1Many(const dlong N,
                                   const dlong Nfields,
                                   const dlong fieldOffset,
                                   occa::memory& o_w,
                               occa::memory& o_a, MPI_Comm _comm) {
#ifdef ENABLE_TIMER
  platform->timer.tic("dotp",1);
#endif
  int Nblock = (N+blocksize-1)/blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if(o_scratch.size() < Nbytes) reallocScratch(Nbytes);

  weightedNorm1ManyKernel(Nblock, N, Nfields, fieldOffset, o_w, o_a, o_scratch);


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
  return norm;
}
