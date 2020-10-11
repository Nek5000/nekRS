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

#ifndef LINALG_HPP
#define LINALG_HPP

#include "mpi.h"
#include "occa.hpp"
#include "types.h"

using std::string;

class linAlg_t {
private:
  occa::device device;
  occa::properties kernelInfo;
  MPI_Comm comm;
  int blocksize;

  //scratch space for reductions
  dfloat* scratch;
  occa::memory h_scratch;
  occa::memory o_scratch;

  void setup();
public:
  linAlg_t(occa::device& _device, occa::properties*& _kernelInfo, MPI_Comm& _comm) {
    blocksize = 256; 
    device = _device;
    kernelInfo = *(_kernelInfo);
    comm = _comm;
    setup();
  }

  ~linAlg_t();

  /*********************/
  /* vector operations */
  /*********************/

  // o_a[n] = alpha
  void fill(const dlong N, const dfloat alpha, occa::memory& o_a);

  // o_a[n] *= alpha
  void scale(const dlong N, const dfloat alpha, occa::memory& o_a);

  // o_a[n] += alpha
  void add(const dlong N, const dfloat alpha, occa::memory& o_a);

  // o_y[n] = beta*o_y[n] + alpha*o_x[n]
  void axpby(const dlong N, const dfloat alpha, occa::memory& o_x,
                            const dfloat beta,  occa::memory& o_y);

  // o_z[n] = beta*o_y[n] + alpha*o_x[n]
  void axpbyz(const dlong N, const dfloat alpha, occa::memory& o_x,
                             const dfloat beta,  occa::memory& o_y,
                             occa::memory& o_z);

  // o_y[n] = alpha*o_x[n]*o_y[n]
  void axmy(const dlong N, const dfloat alpha,
            occa::memory& o_x, occa::memory& o_y);

  // o_z[n] = alpha*o_x[n]*o_y[n] (new)
  void axmyz(const dlong N, const dfloat alpha,
             occa::memory& o_x, occa::memory& o_y,
             occa::memory& o_z);

  // o_y[n] = alpha*o_x[n]/o_y[n]
  void axdy(const dlong N, const dfloat alpha,
            occa::memory& o_x, occa::memory& o_y);

  // o_z[n] = alpha*o_x[n]*o_y[n]
  void axdyz(const dlong N, const dfloat alpha,
             occa::memory& o_x, occa::memory& o_y,
             occa::memory& o_z);

  // \sum o_a
  dfloat sum(const dlong N, occa::memory& o_a, MPI_Comm _comm);

  // \min o_a
  dfloat min(const dlong N, occa::memory& o_a, MPI_Comm _comm);

  // \max o_a
  dfloat max(const dlong N, occa::memory& o_a, MPI_Comm _comm);

  // ||o_a||_2
  //dfloat norm2(const dlong N, occa::memory& o_a, MPI_Comm _comm);

  // o_x.o_y
  dfloat innerProd(const dlong N, occa::memory& o_x, occa::memory& o_y,
                    MPI_Comm _comm);

  // ||o_a||_w2
  dfloat weightedNorm2(const dlong N, occa::memory& o_w, occa::memory& o_a,
                       MPI_Comm _comm);

  // o_w.o_x.o_y
  dfloat weightedInnerProd(const dlong N, occa::memory& o_w, occa::memory& o_x,
                            occa::memory& o_y, MPI_Comm _comm);

  occa::kernel fillKernel;
  occa::kernel addKernel;
  occa::kernel scaleKernel;
  occa::kernel axpbyKernel;
  occa::kernel axpbyzKernel;
  occa::kernel axmyKernel;
  occa::kernel axmyzKernel;
  occa::kernel axdyKernel;
  occa::kernel axdyzKernel;
  occa::kernel sumKernel;
  occa::kernel minKernel;
  occa::kernel maxKernel;
  occa::kernel norm2Kernel;
  occa::kernel weightedNorm2Kernel;
  occa::kernel innerProdKernel;
  occa::kernel weightedInnerProdKernel;
};

#endif
