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

#include "nrssys.hpp"

#define USE_WEIGHTED_INNER_PROD_MULTI_DEVICE 0

class linAlg_t {
private:
  occa::properties kernelInfo;
  MPI_Comm comm;
  int blocksize;
  bool serial;

  int timer = 0;

  //scratch space for reductions
  dfloat* scratch;
  occa::memory h_scratch;
  occa::memory o_scratch;

  void setup();
  void reallocScratch(const size_t Nbytes);

  void runTimers();

  ~linAlg_t();
  linAlg_t();
  static linAlg_t* singleton;
public:
  static linAlg_t* getInstance();

  void enableTimer();
  void disableTimer();

  /*********************/
  /* vector operations */
  /*********************/

  // o_a[n] = alpha
  void fill(const dlong N, const dfloat alpha, occa::memory& o_a);
  void pfill(const dlong N, const pfloat alpha, occa::memory& o_a);

  // o_a[n] = abs(o_a[n])
  void abs(const dlong N, occa::memory& o_a);

  // o_a[n] *= alpha
  void scale(const dlong N, const dfloat alpha, occa::memory& o_a);
  void scaleMany(const dlong N, const dlong Nfields, const dlong fieldOffset, const dfloat alpha, occa::memory& o_a, const dlong offset = 0);

  // o_a[n] += alpha
  void add(const dlong N, const dfloat alpha, occa::memory& o_a, const dlong offset = 0);

  // o_y[n] = beta*o_y[n] + alpha*o_x[n]
  void axpby(const dlong N, const dfloat alpha, occa::memory& o_x,
                            const dfloat beta,  occa::memory& o_y,
                            const dlong xOffset = 0, const dlong yOffset = 0);
  void paxpby(const dlong N, const pfloat alpha, occa::memory& o_x,
                             const pfloat beta,  occa::memory& o_y,
                             const dlong xOffset = 0, const dlong yOffset = 0);
  void axpbyMany(const dlong N, const dlong Nfields, const dlong offset, const dfloat alpha, occa::memory& o_x,
                            const dfloat beta,  occa::memory& o_y);
  void paxpbyMany(const dlong N, const dlong Nfields, const dlong offset, const pfloat alpha, occa::memory& o_x,
                            const pfloat beta,  occa::memory& o_y);

  // o_z[n] = beta*o_y[n] + alpha*o_x[n]
  void axpbyz(const dlong N, const dfloat alpha, occa::memory& o_x,
                             const dfloat beta,  occa::memory& o_y,
                             occa::memory& o_z);
  void axpbyzMany(const dlong N, const dlong Nfields, const dlong offset, const dfloat alpha, occa::memory& o_x,
                             const dfloat beta,  occa::memory& o_y,
                             occa::memory& o_z);

  // o_y[n] = alpha*o_x[n]*o_y[n]
  void axmy(const dlong N, const dfloat alpha,
            occa::memory& o_x, occa::memory& o_y);
  void paxmy(const dlong N, const pfloat alpha,
            occa::memory& o_x, occa::memory& o_y);

  // mode 1:
  // o_y[n,fld] = alpha*o_x[n,fld]*o_y[n,fld]
  // mode 0:
  // o_y[n,fld] = alpha*o_x[n]*o_y[n,fld]
  void axmyMany(const dlong N, 
            const dlong Nfields,
            const dlong offset, const dlong mode,
            const dfloat alpha,
            occa::memory& o_x, occa::memory& o_y);
  void axmyVector(const dlong N, 
            const dlong offset, const dlong mode,
            const dfloat alpha,
            occa::memory& o_x, occa::memory& o_y);

  // o_z[n] = alpha*o_x[n]*o_y[n] (new)
  void axmyz(const dlong N, const dfloat alpha,
             occa::memory& o_x, occa::memory& o_y,
             occa::memory& o_z);
  void paxmyz(const dlong N, const pfloat alpha,
             occa::memory& o_x, occa::memory& o_y,
             occa::memory& o_z);
  void axmyzMany(const dlong N, const dlong Nfields, const dlong offset, const dfloat alpha,
             occa::memory& o_x, occa::memory& o_y,
             occa::memory& o_z);
  void paxmyzMany(const dlong N, const dlong Nfields, const dlong offset, const pfloat alpha,
             occa::memory& o_x, occa::memory& o_y,
             occa::memory& o_z);

  // o_y[n] = alpha/o_y[n]
  void ady(const dlong N, const dfloat alpha,
            occa::memory& o_y);
  void adyMany(const dlong N, const dlong Nfields, const dlong offset, const dfloat alpha,
            occa::memory& o_y);
  void padyMany(const dlong N, const dlong Nfields, const dlong offset, const pfloat alpha,
            occa::memory& o_y);
  // o_y[n] = alpha*o_x[n]/o_y[n]
  void axdy(const dlong N, const dfloat alpha,
            occa::memory& o_x, occa::memory& o_y);

  // o_z[n] = alpha*o_x[n]*o_y[n]
  void axdyz(const dlong N, const dfloat alpha,
             occa::memory& o_x, occa::memory& o_y,
             occa::memory& o_z);
  // o_y[n] = alpha*o_y[n]/o_x[n]
  void aydx(const dlong N, const dfloat alpha,
    occa::memory& o_x, occa::memory& o_y);
  void aydxMany(const dlong N, const dlong Nfields, const dlong fieldOffset,
    const dlong mode, const dfloat alpha,
    occa::memory& o_x, occa::memory& o_y);

  // \sum o_a
  dfloat sum(const dlong N, occa::memory& o_a, MPI_Comm _comm, const dlong offset = 0);
  dfloat sumMany(const dlong N, const dlong Nfields, const dlong fieldOffset, occa::memory& o_a, MPI_Comm _comm);

  // \min o_a
  dfloat min(const dlong N, occa::memory& o_a, MPI_Comm _comm);

  // \max o_a
  dfloat max(const dlong N, occa::memory& o_a, MPI_Comm _comm);

  // ||o_a||_\infty
  dfloat amax(const dlong N, occa::memory &o_a, MPI_Comm _comm);
  dfloat
  amaxMany(const dlong N, const dlong Nfields, const dlong fieldOffset, occa::memory &o_a, MPI_Comm _comm);

  // ||o_a||_2
  dfloat norm2(const dlong N, occa::memory& o_a, MPI_Comm _comm);
  dfloat norm2Many(const dlong N, const dlong Nfields, const dlong fieldOffset, occa::memory& o_a, MPI_Comm _comm);

  // ||o_a||_1
  dfloat norm1(const dlong N, occa::memory& o_a, MPI_Comm _comm);
  dfloat norm1Many(const dlong N, const dlong Nfields, const dlong fieldOffset, occa::memory& o_a, MPI_Comm _comm);

  // o_x.o_y
  dfloat innerProd(const dlong N, occa::memory& o_x, occa::memory& o_y,
                    MPI_Comm _comm, const dlong offset = 0);

  // ||o_a||_w1
  dfloat weightedNorm1(const dlong N, occa::memory& o_w, occa::memory& o_a,
                       MPI_Comm _comm);
  dfloat weightedNorm1Many(const dlong N,
                           const dlong Nfields,
                           const dlong fieldOffset,
                           occa::memory& o_w, occa::memory& o_a,
                           MPI_Comm _comm);
  // ||o_a||_w2
  dfloat weightedNorm2(const dlong N, occa::memory& o_w, occa::memory& o_a,
                       MPI_Comm _comm);
  dfloat weightedNorm2Many(const dlong N,
                           const dlong Nfields,
                           const dlong fieldOffset,
                           occa::memory& o_w, occa::memory& o_a,
                           MPI_Comm _comm);

  // o_w.o_x.o_y
  dfloat weightedInnerProd(const dlong N, occa::memory& o_w, occa::memory& o_x,
                            occa::memory& o_y, MPI_Comm _comm);
  void weightedInnerProdMulti(const dlong N, const dlong NVec, const dlong Nfields, 
                              const dlong fieldOffset, occa::memory& o_w, occa::memory& o_x,
                              occa::memory& o_y, MPI_Comm _comm,
                              dfloat* result, const dlong offset = 0);
  void weightedInnerProdMulti(const dlong N, const dlong NVec, const dlong Nfields, 
                              const dlong fieldOffset, occa::memory& o_w, occa::memory& o_x,
                              occa::memory& o_y, MPI_Comm _comm,
                              occa::memory& o_result, const dlong offset = 0);

  dfloat weightedInnerProdMany(const dlong N,
                               const dlong Nfields, const dlong fieldOffset, occa::memory& o_w, occa::memory& o_x,
                            occa::memory& o_y, MPI_Comm _comm);

  // z = x \cross y
  void crossProduct(const dlong N,
                    const dlong fieldOffset,
                    occa::memory &o_x,
                    occa::memory &o_y,
                    occa::memory &o_z);

  void unitVector(const dlong N, const dlong fieldOffset, occa::memory &o_v);

  occa::kernel fillKernel;
  occa::kernel pfillKernel;

  occa::kernel absKernel;
  occa::kernel addKernel;
  occa::kernel scaleKernel;
  occa::kernel scaleManyKernel;
  occa::kernel axpbyKernel;
  occa::kernel paxpbyKernel;
  occa::kernel axpbyManyKernel;
  occa::kernel paxpbyManyKernel;
  occa::kernel axpbyzKernel;
  occa::kernel paxpbyzKernel;
  occa::kernel axpbyzManyKernel;
  occa::kernel axmyKernel;
  occa::kernel paxmyKernel;
  occa::kernel axmyManyKernel;
  occa::kernel axmyVectorKernel;
  occa::kernel axmyzKernel;
  occa::kernel paxmyzKernel;
  occa::kernel axmyzManyKernel;
  occa::kernel paxmyzManyKernel;
  occa::kernel axdyKernel;
  occa::kernel aydxKernel;
  occa::kernel aydxManyKernel;
  occa::kernel adyKernel;
  occa::kernel adyManyKernel;
  occa::kernel padyManyKernel;
  occa::kernel axdyzKernel;
  occa::kernel sumKernel;
  occa::kernel sumManyKernel;
  occa::kernel sumFieldKernel;
  occa::kernel minKernel;
  occa::kernel maxKernel;
  occa::kernel amaxKernel;
  occa::kernel amaxManyKernel;
  occa::kernel norm2Kernel;
  occa::kernel norm2ManyKernel;
  occa::kernel norm1Kernel;
  occa::kernel norm1ManyKernel;
  occa::kernel weightedNorm1Kernel;
  occa::kernel weightedNorm1ManyKernel;
  occa::kernel weightedNorm2Kernel;
  occa::kernel weightedNorm2ManyKernel;
  occa::kernel innerProdKernel;
  occa::kernel weightedInnerProdKernel;
  occa::kernel weightedInnerProdManyKernel;
  occa::kernel weightedInnerProdMultiKernel;
  occa::kernel weightedInnerProdMultiDeviceKernel;
  occa::kernel crossProductKernel;
  occa::kernel unitVectorKernel;
};

#endif
