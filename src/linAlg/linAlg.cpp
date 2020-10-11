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

void linAlg_t::setup() {

  int rank;
  MPI_Comm_rank(comm, &rank);

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)blocksize;

  //pinned scratch buffer
  {
    occa::properties props = kernelInfo;
    props["mapped"] = true;
    h_scratch = device.malloc(blocksize*sizeof(dfloat), props);
    scratch = (dfloat*) h_scratch.ptr(props);
  }
  o_scratch = device.malloc(blocksize*sizeof(dfloat));

  string oklDir;
  oklDir.assign(getenv("NEKRS_INSTALL_DIR"));
  oklDir += "/okl/linAlg/";

  for (int r = 0; r < 2; r++) {
    if ((r == 0 && rank == 0) || (r == 1 && rank > 0)) {
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
      if (axpbyKernel.isInitialized()==false)
        axpbyKernel = device.buildKernel(oklDir + 
                                         "linAlgAXPBY.okl",
                                         "axpby",
                                         kernelInfo);
      if (axpbyzKernel.isInitialized()==false)
        axpbyzKernel = device.buildKernel(oklDir + 
                                          "linAlgAXPBY.okl",
                                          "axpbyz",
                                          kernelInfo);
      if (axmyKernel.isInitialized()==false)
        axmyKernel = device.buildKernel(oklDir + 
                                        "linAlgAXMY.okl",
                                        "axmy",
                                        kernelInfo);
      if (axmyzKernel.isInitialized()==false)
        axmyzKernel = device.buildKernel(oklDir + 
                                         "linAlgAXMY.okl",
                                         "axmyz",
                                         kernelInfo);
      if (axdyKernel.isInitialized()==false)
        axdyKernel = device.buildKernel(oklDir + 
                                        "linAlgAXDY.okl",
                                        "axdy",
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
      if (weightedNorm2Kernel.isInitialized()==false)
        weightedNorm2Kernel = device.buildKernel(oklDir + 
                                        "linAlgWeightedNorm2.okl",
                                        "weightedNorm2",
                                        kernelInfo);
      if (innerProdKernel.isInitialized()==false)
        innerProdKernel = device.buildKernel(oklDir + 
                                        "linAlgInnerProd.okl",
                                        "innerProd",
                                        kernelInfo);
      if (weightedInnerProdKernel.isInitialized()==false)
        weightedInnerProdKernel = device.buildKernel(oklDir + 
                                        "linAlgWeightedInnerProd.okl",
                                        "weightedInnerProd",
                                        kernelInfo);
    }
    MPI_Barrier(comm);
  }
}

linAlg_t::~linAlg_t() {
  fillKernel.free();
  addKernel.free();
  scaleKernel.free();
  axpbyKernel.free();
  axpbyzKernel.free();
  axmyKernel.free();
  axmyzKernel.free();
  axdyKernel.free();
  axdyzKernel.free();
  sumKernel.free();
  minKernel.free();
  maxKernel.free();
  norm2Kernel.free();
  weightedNorm2Kernel.free();
  innerProdKernel.free();
  weightedInnerProdKernel.free();
}

/*********************/
/* vector operations */
/*********************/

// o_a[n] = alpha
void linAlg_t::fill(const dlong N, const dfloat alpha, occa::memory& o_a) {
  fillKernel(N, alpha, o_a);
}

// o_a[n] += alpha
void linAlg_t::add(const dlong N, const dfloat alpha, occa::memory& o_a) {
  addKernel(N, alpha, o_a);
}

// o_a[n] *= alpha
void linAlg_t::scale(const dlong N, const dfloat alpha, occa::memory& o_a)  {
  scaleKernel(N, alpha, o_a);
}

// o_y[n] = beta*o_y[n] + alpha*o_x[n]
void linAlg_t::axpby(const dlong N, const dfloat alpha, occa::memory& o_x,
                    const dfloat beta,  occa::memory& o_y) {
  axpbyKernel(N, alpha, o_x, beta, o_y);
}

// o_z[n] = beta*o_y[n] + alpha*o_x[n]
void linAlg_t::axpbyz(const dlong N, const dfloat alpha, occa::memory& o_x,
                      const dfloat beta, occa::memory& o_y, occa::memory& o_z) {
  axpbyzKernel(N, alpha, o_x, beta, o_y, o_z);
}

// o_y[n] = alpha*o_x[n]*o_y[n]
void linAlg_t::axmy(const dlong N, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y) {
  axmyKernel(N, alpha, o_x, o_y);
}

// o_z[n] = alpha*o_x[n]*o_y[n]
void linAlg_t::axmyz(const dlong N, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y, occa::memory& o_z) {
  axmyzKernel(N, alpha, o_x, o_y, o_z);
}

// o_y[n] = alpha*o_x[n]/o_y[n]
void linAlg_t::axdy(const dlong N, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y) {
  axdyKernel(N, alpha, o_x, o_y);
}

// o_z[n] = alpha*o_x[n]/o_y[n]
void linAlg_t::axdyz(const dlong N, const dfloat alpha,
                   occa::memory& o_x, occa::memory& o_y, occa::memory& o_z) {
  axdyzKernel(N, alpha, o_x, o_y, o_z);
}

// \sum o_a
dfloat linAlg_t::sum(const dlong N, occa::memory& o_a, MPI_Comm _comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  sumKernel(Nblock, N, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(dfloat));

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
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  minKernel(Nblock, N, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(dfloat));

  dfloat min = 9e30;
  for(dlong n=0;n<Nblock;++n){
    min = (scratch[n] < min) ? scratch[n]:min;
  }

  MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DFLOAT, MPI_MIN, _comm);

  return min;
}

// \max o_a
dfloat linAlg_t::max(const dlong N, occa::memory& o_a, MPI_Comm _comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  maxKernel(Nblock, N, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(dfloat));

  dfloat max = -9e30;
  for(dlong n=0;n<Nblock;++n){
    max = (scratch[n] < max) ? scratch[n]:max;
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
                           MPI_Comm _comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  innerProdKernel(Nblock, N, o_x, o_y, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(dfloat));

  dfloat dot = 0;
  for(dlong n=0;n<Nblock;++n){
    dot += scratch[n];
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_DFLOAT, MPI_SUM, _comm);

  return dot;
}

// o_w.o_x.o_y
dfloat linAlg_t::weightedInnerProd(const dlong N, occa::memory& o_w,
                                   occa::memory& o_x, occa::memory& o_y,
                                   MPI_Comm _comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  weightedInnerProdKernel(Nblock, N, o_w, o_x, o_y, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(dfloat));

  dfloat dot = 0;
  for(dlong n=0;n<Nblock;++n){
    dot += scratch[n];
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_DFLOAT, MPI_SUM, _comm);

  return dot;
}

// ||o_a||_w2
dfloat linAlg_t::weightedNorm2(const dlong N, occa::memory& o_w,
                               occa::memory& o_a, MPI_Comm _comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  weightedNorm2Kernel(Nblock, N, o_w, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(dfloat));

  dfloat norm = 0;
  for(dlong n=0;n<Nblock;++n){
    norm += scratch[n];
  }

  if (_comm != MPI_COMM_NULL) 
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);

  return sqrt(norm);
}
