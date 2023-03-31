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
#include "re2Reader.hpp"

linAlg_t *linAlg_t::singleton = nullptr;

void linAlg_t::runTimers()
{
  int nelgt, nelgv;
  const std::string meshFile = platform->options.getArgs("MESH FILE");
  re2::nelg(meshFile, nelgt, nelgv, platform->comm.mpiComm);
  const int nel = nelgv / platform->comm.mpiCommSize;

  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);

  const auto fields = 1;
  const auto Nlocal = nel * (N + 1) * (N + 1) * (N + 1);
  auto o_weight = platform->device.malloc(Nlocal * sizeof(dfloat));
  auto o_r = platform->device.malloc(Nlocal * sizeof(dfloat));
  auto o_z = platform->device.malloc(Nlocal * sizeof(dfloat));

  const auto Nrep = 20;

  {
    // warm-up
    weightedInnerProdMany(Nlocal, fields, 1, o_weight, o_r, o_z, platform->comm.mpiComm);

    platform->device.finish();
    MPI_Barrier(platform->comm.mpiComm);
    const auto tStart = MPI_Wtime();
    for (int i = 0; i < Nrep; i++) {
      weightedInnerProdMany(Nlocal, fields, 1, o_weight, o_r, o_z, platform->comm.mpiComm);
    }
    platform->device.finish();
    const auto elapsed = (MPI_Wtime() - tStart) / Nrep;
    auto elapsedMax = 0.0;
    MPI_Allreduce(&elapsed, &elapsedMax, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0)
      printf("wdotp: %.3es  ", elapsedMax);
  }

  if (platform->comm.mpiCommSize > 1) {
    platform->device.finish();
    MPI_Barrier(platform->comm.mpiComm);
    const auto tStart = MPI_Wtime();
    for (int i = 0; i < Nrep; i++) {
      weightedInnerProdMany(Nlocal, fields, 1, o_weight, o_r, o_z, MPI_COMM_SELF);
    }
    platform->device.finish();
    const auto elapsed = (MPI_Wtime() - tStart) / Nrep;
    auto elapsedMax = 0.0;
    MPI_Allreduce(&elapsed, &elapsedMax, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0)
      printf("(local: %.3es)\n", elapsedMax);
  }
  else {
    if (platform->comm.mpiRank == 0)
      printf("\n");
  }

  if (platform->comm.mpiRank == 0)
    std::cout << std::endl;

  o_weight.free();
  o_r.free();
  o_z.free();
}

linAlg_t *linAlg_t::getInstance()
{
  if (!singleton)
    singleton = new linAlg_t();
  return singleton;
}
linAlg_t::linAlg_t()
{
  blocksize = BLOCKSIZE;
  serial = platform->serial;
  comm = platform->comm.mpiComm;
  timer = 0;

  if (platform->comm.mpiRank == 0)
    std::cout << "initializing linAlg ...\n";

  setup();
  runTimers();

  if (platform->options.compareArgs("ENABLE LINALG TIMER", "TRUE"))
    timer = 1;
}
void linAlg_t::enableTimer() { timer = 1; }
void linAlg_t::disableTimer() { timer = 0; }
void linAlg_t::reallocScratch(const size_t Nbytes)
{
  device_t &device = platform->device;
  if (h_scratch.size())
    h_scratch.free();
  if (o_scratch.size())
    o_scratch.free();
  // pinned scratch buffer
  {
    h_scratch = device.mallocHost(Nbytes);
    scratch = (dfloat *)h_scratch.ptr();
  }
  o_scratch = device.malloc(Nbytes);
}
void linAlg_t::setup()
{
  auto &kernels = platform->kernels;
  int rank;
  MPI_Comm_rank(comm, &rank);

  occa::properties kernelInfo = platform->kernelInfo;

  reallocScratch(blocksize * sizeof(dfloat));

  std::string oklDir(getenv("NEKRS_KERNEL_DIR"));
  oklDir += "/linAlg/";

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  {
    fillKernel = kernels.get("fill");
    pfillKernel = kernels.get("pfill");
    absKernel = kernels.get("vabs");
    addKernel = kernels.get("add");
    scaleKernel = kernels.get("scale");
    scaleManyKernel = kernels.get("scaleMany");
    axpbyKernel = kernels.get("axpby");
    paxpbyKernel = kernels.get("paxpby");
    axpbyManyKernel = kernels.get("axpbyMany");
    paxpbyManyKernel = kernels.get("paxpbyMany");
    axpbyzKernel = kernels.get("axpbyz");
    axpbyzManyKernel = kernels.get("axpbyzMany");
    axmyKernel = kernels.get("axmy");
    paxmyKernel = kernels.get("paxmy");
    axmyManyKernel = kernels.get("axmyMany");
    axmyVectorKernel = kernels.get("axmyVector");
    axmyzKernel = kernels.get("axmyz");
    paxmyzKernel = kernels.get("paxmyz");
    axmyzManyKernel = kernels.get("axmyzMany");
    paxmyzManyKernel = kernels.get("paxmyzMany");
    adyKernel = kernels.get("ady");
    adyManyKernel = kernels.get("adyMany");
    padyManyKernel = kernels.get("padyMany");
    axdyKernel = kernels.get("axdy");
    aydxKernel = kernels.get("aydx");
    aydxManyKernel = kernels.get("aydxMany");
    axdyzKernel = kernels.get("axdyz");
    sumKernel = kernels.get("sum");
    sumManyKernel = kernels.get("sumMany");
    minKernel = kernels.get("min");
    maxKernel = kernels.get("max");
    amaxKernel = kernels.get("amax");
    amaxManyKernel = kernels.get("amaxMany");
    norm2Kernel = kernels.get("norm2");
    norm2ManyKernel = kernels.get("norm2Many");
    norm1Kernel = kernels.get("norm1");
    norm1ManyKernel = kernels.get("norm1Many");
    weightedNorm1Kernel = kernels.get("weightedNorm1");
    weightedNorm1ManyKernel = kernels.get("weightedNorm1Many");
    weightedNorm2Kernel = kernels.get("weightedNorm2");
    weightedNorm2ManyKernel = kernels.get("weightedNorm2Many");
    innerProdKernel = kernels.get("innerProd");
    weightedInnerProdKernel = kernels.get("weightedInnerProd");
    weightedInnerProdManyKernel = kernels.get("weightedInnerProdMany");
    weightedInnerProdMultiKernel = kernels.get("weightedInnerProdMulti");
    weightedInnerProdMultiDeviceKernel = kernels.get("weightedInnerProdMultiDevice");
    crossProductKernel = kernels.get("crossProduct");
    unitVectorKernel = kernels.get("unitVector");
    entrywiseMagKernel = kernels.get("entrywiseMag");
  }
}

linAlg_t::~linAlg_t()
{
  fillKernel.free();
  pfillKernel.free();
  absKernel.free();
  addKernel.free();
  scaleKernel.free();
  scaleManyKernel.free();
  axpbyKernel.free();
  paxpbyKernel.free();
  axpbyManyKernel.free();
  paxpbyManyKernel.free();
  axpbyzKernel.free();
  axpbyzManyKernel.free();
  axmyKernel.free();
  paxmyKernel.free();
  axmyManyKernel.free();
  axmyVectorKernel.free();
  axmyzKernel.free();
  paxmyzKernel.free();
  axmyzManyKernel.free();
  axdyKernel.free();
  aydxKernel.free();
  aydxManyKernel.free();
  adyKernel.free();
  adyManyKernel.free();
  padyManyKernel.free();
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
void linAlg_t::fill(const dlong N, const dfloat alpha, occa::memory &o_a) { fillKernel(N, alpha, o_a); }

// o_a[n] = alpha
void linAlg_t::pfill(const dlong N, const pfloat alpha, occa::memory &o_a) { pfillKernel(N, alpha, o_a); }

// o_a[n] = abs(o_a[n])
void linAlg_t::abs(const dlong N, occa::memory &o_a) { absKernel(N, o_a); }

// o_a[n] += alpha
void linAlg_t::add(const dlong N, const dfloat alpha, occa::memory &o_a, const dlong offset)
{
  addKernel(N, offset, alpha, o_a);
}

// o_a[n] *= alpha
void linAlg_t::scale(const dlong N, const dfloat alpha, occa::memory &o_a) { scaleKernel(N, alpha, o_a); }
void linAlg_t::scaleMany(const dlong N,
                         const dlong Nfields,
                         const dlong fieldOffset,
                         const dfloat alpha,
                         occa::memory &o_a,
                         const dlong offset)
{
  scaleManyKernel(N, Nfields, fieldOffset, offset, alpha, o_a);
}

// o_y[n] = beta*o_y[n] + alpha*o_x[n]
void linAlg_t::axpby(const dlong N,
                     const dfloat alpha,
                     occa::memory &o_x,
                     const dfloat beta,
                     occa::memory &o_y,
                     const dlong xOffset,
                     const dlong yOffset)
{
  axpbyKernel(N, xOffset, yOffset, alpha, o_x, beta, o_y);
  platform->flopCounter->add("axpby", 3 * static_cast<double>(N));
}

// o_y[n] = beta*o_y[n] + alpha*o_x[n]
void linAlg_t::paxpby(const dlong N,
                      const pfloat alpha,
                      occa::memory &o_x,
                      const pfloat beta,
                      occa::memory &o_y,
                      const dlong xOffset,
                      const dlong yOffset)
{
  paxpbyKernel(N, xOffset, yOffset, alpha, o_x, beta, o_y);
  platform->flopCounter->add("axpby", 0.5 * 3 * static_cast<double>(N));
}

void linAlg_t::axpbyMany(const dlong N,
                         const dlong Nfields,
                         const dlong offset,
                         const dfloat alpha,
                         occa::memory &o_x,
                         const dfloat beta,
                         occa::memory &o_y)
{
  axpbyManyKernel(N, Nfields, offset, alpha, o_x, beta, o_y);
  platform->flopCounter->add("axpbyMany", 3 * static_cast<double>(N) * Nfields);
}

void linAlg_t::paxpbyMany(const dlong N,
                          const dlong Nfields,
                          const dlong offset,
                          const pfloat alpha,
                          occa::memory &o_x,
                          const pfloat beta,
                          occa::memory &o_y)
{
  paxpbyManyKernel(N, Nfields, offset, alpha, o_x, beta, o_y);
  platform->flopCounter->add("axpbyMany", 0.5 * 3 * static_cast<double>(N) * Nfields);
}

// o_z[n] = beta*o_y[n] + alpha*o_x[n]
void linAlg_t::axpbyz(const dlong N,
                      const dfloat alpha,
                      occa::memory &o_x,
                      const dfloat beta,
                      occa::memory &o_y,
                      occa::memory &o_z)
{
  axpbyzKernel(N, alpha, o_x, beta, o_y, o_z);
  platform->flopCounter->add("axpbyz", 3 * static_cast<double>(N));
}
void linAlg_t::axpbyzMany(const dlong N,
                          const dlong Nfields,
                          const dlong fieldOffset,
                          const dfloat alpha,
                          occa::memory &o_x,
                          const dfloat beta,
                          occa::memory &o_y,
                          occa::memory &o_z)
{
  axpbyzManyKernel(N, Nfields, fieldOffset, alpha, o_x, beta, o_y, o_z);
  platform->flopCounter->add("axpbyzMany", 3 * static_cast<double>(N) * Nfields);
}

// o_y[n] = alpha*o_x[n]*o_y[n]
void linAlg_t::axmy(const dlong N, const dfloat alpha, occa::memory &o_x, occa::memory &o_y)
{
  axmyKernel(N, alpha, o_x, o_y);
}
void linAlg_t::paxmy(const dlong N, const pfloat alpha, occa::memory &o_x, occa::memory &o_y)
{
  paxmyKernel(N, alpha, o_x, o_y);
}

void linAlg_t::axmyMany(const dlong N,
                        const dlong Nfields,
                        const dlong offset,
                        const dlong mode,
                        const dfloat alpha,
                        occa::memory &o_x,
                        occa::memory &o_y)
{
  axmyManyKernel(N, Nfields, offset, mode, alpha, o_x, o_y);
}
void linAlg_t::axmyVector(const dlong N,
                          const dlong offset,
                          const dlong mode,
                          const dfloat alpha,
                          occa::memory &o_x,
                          occa::memory &o_y)
{
  axmyVectorKernel(N, offset, mode, alpha, o_x, o_y);
}

// o_z[n] = alpha*o_x[n]*o_y[n]
void linAlg_t::axmyz(const dlong N,
                     const dfloat alpha,
                     occa::memory &o_x,
                     occa::memory &o_y,
                     occa::memory &o_z)
{
  axmyzKernel(N, alpha, o_x, o_y, o_z);
}
void linAlg_t::paxmyz(const dlong N,
                      const pfloat alpha,
                      occa::memory &o_x,
                      occa::memory &o_y,
                      occa::memory &o_z)
{
  paxmyzKernel(N, alpha, o_x, o_y, o_z);
}
void linAlg_t::axmyzMany(const dlong N,
                         const dlong Nfields,
                         const dlong offset,
                         const dfloat alpha,
                         occa::memory &o_x,
                         occa::memory &o_y,
                         occa::memory &o_z)
{
  axmyzManyKernel(N, Nfields, offset, alpha, o_x, o_y, o_z);
}

void linAlg_t::paxmyzMany(const dlong N,
                          const dlong Nfields,
                          const dlong offset,
                          const pfloat alpha,
                          occa::memory &o_x,
                          occa::memory &o_y,
                          occa::memory &o_z)
{
  paxmyzManyKernel(N, Nfields, offset, alpha, o_x, o_y, o_z);
}

// o_y[n] = alpha*o_x[n]/o_y[n]
void linAlg_t::axdy(const dlong N, const dfloat alpha, occa::memory &o_x, occa::memory &o_y)
{
  axdyKernel(N, alpha, o_x, o_y);
}
void linAlg_t::aydx(const dlong N, const dfloat alpha, occa::memory &o_x, occa::memory &o_y)
{
  aydxKernel(N, alpha, o_x, o_y);
}
void linAlg_t::aydxMany(const dlong N,
                        const dlong Nfields,
                        const dlong fieldOffset,
                        const dlong mode,
                        const dfloat alpha,
                        occa::memory &o_x,
                        occa::memory &o_y)
{
  aydxManyKernel(N, Nfields, fieldOffset, mode, alpha, o_x, o_y);
}
// o_y[n] = alpha/o_y[n]
void linAlg_t::ady(const dlong N, const dfloat alpha, occa::memory &o_y) { adyKernel(N, alpha, o_y); }
void linAlg_t::adyMany(const dlong N,
                       const dlong Nfields,
                       const dlong offset,
                       const dfloat alpha,
                       occa::memory &o_y)
{
  adyManyKernel(N, Nfields, offset, alpha, o_y);
}
void linAlg_t::padyMany(const dlong N,
                        const dlong Nfields,
                        const dlong offset,
                        const pfloat alpha,
                        occa::memory &o_y)
{
  padyManyKernel(N, Nfields, offset, alpha, o_y);
}

// o_z[n] = alpha*o_x[n]/o_y[n]
void linAlg_t::axdyz(const dlong N,
                     const dfloat alpha,
                     occa::memory &o_x,
                     occa::memory &o_y,
                     occa::memory &o_z)
{
  axdyzKernel(N, alpha, o_x, o_y, o_z);
}

// \sum o_a
dfloat linAlg_t::sum(const dlong N, occa::memory &o_a, MPI_Comm _comm, const dlong offset)
{
  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  if (N > 1) {
    sumKernel(Nblock, N, offset, o_a, o_scratch);
    o_scratch.copyTo(scratch, Nbytes);
  }
  else {
    o_a.copyTo(scratch, Nbytes);
  }

  dfloat sum = 0;
  for (dlong n = 0; n < Nblock; ++n) {
    sum += scratch[n];
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DFLOAT, MPI_SUM, _comm);

  return sum;
}
dfloat linAlg_t::sumMany(const dlong N,
                         const dlong Nfields,
                         const dlong fieldOffset,
                         occa::memory &o_a,
                         MPI_Comm _comm)
{
  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  if (N > 1 || Nfields > 1) {
    sumManyKernel(Nblock, N, Nfields, fieldOffset, o_a, o_scratch);

    o_scratch.copyTo(scratch, Nbytes);
  }
  else {
    o_a.copyTo(scratch, Nbytes);
  }

  dfloat sum = 0;
  for (dlong n = 0; n < Nblock; ++n) {
    sum += scratch[n];
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DFLOAT, MPI_SUM, _comm);

  return sum;
}

// \min o_a
dfloat linAlg_t::min(const dlong N, occa::memory &o_a, MPI_Comm _comm)
{
  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  if (N > 1) {
    minKernel(Nblock, N, o_a, o_scratch);

    o_scratch.copyTo(scratch, Nbytes);
  }
  else {
    o_a.copyTo(scratch, Nbytes);
  }

  dfloat min = scratch[0];
  for (dlong n = 1; n < Nblock; ++n) {
    min = (scratch[n] < min) ? scratch[n] : min;
  }

  MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DFLOAT, MPI_MIN, _comm);

  return min;
}

// \max o_a
dfloat linAlg_t::max(const dlong N, occa::memory &o_a, MPI_Comm _comm)
{
  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  if (N > 1) {
    maxKernel(Nblock, N, o_a, o_scratch);

    o_scratch.copyTo(scratch, Nbytes);
  }
  else {
    o_a.copyTo(scratch, Nbytes);
  }

  dfloat max = scratch[0];
  for (dlong n = 1; n < Nblock; ++n) {
    max = (scratch[n] > max) ? scratch[n] : max;
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DFLOAT, MPI_MAX, _comm);

  return max;
}

// ||o_a||_\infty
dfloat linAlg_t::amax(const dlong N, occa::memory &o_a, MPI_Comm _comm)
{
  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  if (N > 1) {
    amaxKernel(Nblock, N, o_a, o_scratch);

    o_scratch.copyTo(scratch, Nbytes);
  }
  else {
    o_a.copyTo(scratch, Nbytes);
  }

  dfloat max = scratch[0];
  for (dlong n = 1; n < Nblock; ++n) {
    max = (scratch[n] > max) ? scratch[n] : max;
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DFLOAT, MPI_MAX, _comm);

  return max;
}

dfloat linAlg_t::amaxMany(const dlong N,
                          const dlong Nfields,
                          const dlong fieldOffset,
                          occa::memory &o_x,
                          MPI_Comm _comm)
{
  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat max = 0;
  if (N > 1 || Nfields > 1) {
    amaxManyKernel(Nblock, N, Nfields, fieldOffset, o_x, o_scratch);
    if (serial) {
      max = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        max = std::max(max, scratch[n]);
      }
    }
  }
  else {
    dfloat x;
    o_x.copyTo(&x, Nbytes);
    max = std::abs(x);
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DFLOAT, MPI_MAX, _comm);

  return max;
}

// ||o_a||_2
dfloat linAlg_t::norm2(const dlong N, occa::memory &o_x, MPI_Comm _comm)
{
  if (timer)
    platform->timer.tic("dotp", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat norm = 0;
  if (N > 1) {
    norm2Kernel(Nblock, N, o_x, o_scratch);

    if (serial) {
      norm = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        norm += scratch[n];
      }
    }
  }
  else {
    dfloat x;
    o_x.copyTo(&x, Nbytes);
    norm = x * x;
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  return sqrt(norm);
}
dfloat linAlg_t::norm2Many(const dlong N,
                           const dlong Nfields,
                           const dlong fieldOffset,
                           occa::memory &o_x,
                           MPI_Comm _comm)
{
  if (timer)
    platform->timer.tic("dotp", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat norm = 0;
  if (N > 1 || Nfields > 1) {
    norm2ManyKernel(Nblock, N, Nfields, fieldOffset, o_x, o_scratch);
    if (serial) {
      norm = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        norm += scratch[n];
      }
    }
  }
  else {
    dfloat x;
    o_x.copyTo(&x, Nbytes);
    norm = x * x;
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  return sqrt(norm);
}
// ||o_a||_1
dfloat linAlg_t::norm1(const dlong N, occa::memory &o_x, MPI_Comm _comm)
{
  if (timer)
    platform->timer.tic("dotp", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat norm = 0;
  if (N > 1) {
    norm1Kernel(Nblock, N, o_x, o_scratch);
    if (serial) {
      norm = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        norm += scratch[n];
      }
    }
  }
  else {
    dfloat x;
    o_x.copyTo(&x, Nbytes);
    norm = std::abs(x);
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  return norm;
}
dfloat linAlg_t::norm1Many(const dlong N,
                           const dlong Nfields,
                           const dlong fieldOffset,
                           occa::memory &o_x,
                           MPI_Comm _comm)
{
  if (timer)
    platform->timer.tic("dotp", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat norm = 0;
  if (N > 1 || Nfields > 1) {
    norm1ManyKernel(Nblock, N, Nfields, fieldOffset, o_x, o_scratch);

    if (serial) {
      norm = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        norm += scratch[n];
      }
    }
  }
  else {
    dfloat x;
    o_x.copyTo(&x, Nbytes);
    norm = std::abs(x);
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  return norm;
}

// o_x.o_y
dfloat
linAlg_t::innerProd(const dlong N, occa::memory &o_x, occa::memory &o_y, MPI_Comm _comm, const dlong offset)
{

  if (timer)
    platform->timer.tic("dotp", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat dot = 0;
  if (N > 1) {
    innerProdKernel(Nblock, N, offset, o_x, o_y, o_scratch);

    if (serial) {
      dot = *((dfloat *)o_scratch.ptr());
    }
    else {

      o_scratch.copyTo(scratch, Nbytes);

      for (dlong n = 0; n < Nblock; ++n) {
        dot += scratch[n];
      }
    }
  }
  else {
    dfloat x, y;
    o_x.copyTo(&x, Nbytes);
    o_y.copyTo(&y, Nbytes);
    dot = x * y;
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  return dot;
}

// o_w.o_x.o_y
dfloat linAlg_t::weightedInnerProd(const dlong N,
                                   occa::memory &o_w,
                                   occa::memory &o_x,
                                   occa::memory &o_y,
                                   MPI_Comm _comm)
{

  if (timer)
    platform->timer.tic("dotp", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat dot = 0;
  if (N > 1) {
    weightedInnerProdKernel(Nblock, N, o_w, o_x, o_y, o_scratch);

    if (serial) {
      dot = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        dot += scratch[n];
      }
    }
  }
  else {
    dfloat w, x, y;
    o_w.copyTo(&w, Nbytes);
    o_x.copyTo(&x, Nbytes);
    o_y.copyTo(&y, Nbytes);
    dot = w * x * y;
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  platform->flopCounter->add("weightedInnerProd", 3 * static_cast<double>(N));
  return dot;
}
void linAlg_t::weightedInnerProdMulti(const dlong N,
                                      const dlong NVec,
                                      const dlong Nfields,
                                      const dlong fieldOffset,
                                      occa::memory &o_w,
                                      occa::memory &o_x,
                                      occa::memory &o_y,
                                      MPI_Comm _comm,
                                      dfloat *result,
                                      const dlong offset)
{
  if (timer)
    platform->timer.tic("dotpMulti", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = NVec * Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  if (N > 1 || NVec > 1 || Nfields > 1) {
    weightedInnerProdMultiKernel(Nblock, N, Nfields, fieldOffset, NVec, offset, o_w, o_x, o_y, o_scratch);

    o_scratch.copyTo(scratch, Nbytes);

    for (int field = 0; field < NVec; ++field) {
      dfloat dot = 0;
      for (dlong n = 0; n < Nblock; ++n) {
        dot += scratch[n + field * Nblock];
      }
      result[field] = dot;
    }
  }
  else {
    dfloat w, x, y;
    o_w.copyTo(&w, Nbytes);
    o_x.copyTo(&x, Nbytes);
    o_y.copyTo(&y, Nbytes);
    result[0] = w * x * y;
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, result, NVec, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotpMulti");

  platform->flopCounter->add("weightedInnerProdMulti", NVec * static_cast<double>(N) * (2 * Nfields + 1));
}

void linAlg_t::weightedInnerProdMulti(const dlong N,
                                      const dlong NVec,
                                      const dlong Nfields,
                                      const dlong fieldOffset,
                                      occa::memory &o_w,
                                      occa::memory &o_x,
                                      occa::memory &o_y,
                                      MPI_Comm _comm,
                                      occa::memory &o_result,
                                      const dlong offset)
{
  if (timer)
    platform->timer.tic("dotpMulti", 1);

  const int Nblock = (N + blocksize - 1) / blocksize;

  if (N > 1 || NVec > 1 || Nfields > 1)
    weightedInnerProdMultiDeviceKernel(Nblock,
                                       N,
                                       Nfields,
                                       fieldOffset,
                                       NVec,
                                       offset,
                                       o_w,
                                       o_x,
                                       o_y,
                                       o_result);

  if (_comm != MPI_COMM_SELF) {
    platform->device.finish();
    MPI_Allreduce(MPI_IN_PLACE, o_result.ptr(), NVec, MPI_DFLOAT, MPI_SUM, _comm);
  }

  if (timer)
    platform->timer.toc("dotpMulti");

  platform->flopCounter->add("weightedInnerProdMulti", NVec * static_cast<double>(N) * (2 * Nfields + 1));
}

dfloat linAlg_t::weightedInnerProdMany(const dlong N,
                                       const dlong Nfields,
                                       const dlong fieldOffset,
                                       occa::memory &o_w,
                                       occa::memory &o_x,
                                       occa::memory &o_y,
                                       MPI_Comm _comm)
{
  if (timer)
    platform->timer.tic("dotp", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat dot = 0;
  if (N > 1 || Nfields > 1) {
    weightedInnerProdManyKernel(Nblock, N, Nfields, fieldOffset, o_w, o_x, o_y, o_scratch);

    if (serial) {
      dot = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        dot += scratch[n];
      }
    }
  }
  else {
    dfloat w, x, y;
    o_w.copyTo(&w, Nbytes);
    o_x.copyTo(&x, Nbytes);
    o_y.copyTo(&y, Nbytes);
    dot = w * x * y;
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  platform->flopCounter->add("weightedInnerProdMany", 3 * static_cast<double>(N) * Nfields);

  return dot;
}

// ||o_a||_w2
dfloat linAlg_t::weightedNorm2(const dlong N, occa::memory &o_w, occa::memory &o_a, MPI_Comm _comm)
{
  if (timer)
    platform->timer.tic("dotp", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat norm = 0;
  if (N > 1) {
    weightedNorm2Kernel(Nblock, N, o_w, o_a, o_scratch);

    if (serial) {
      norm = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        norm += scratch[n];
      }
    }
  }
  else {
    dfloat w, a;
    o_w.copyTo(&w, Nbytes);
    o_a.copyTo(&a, Nbytes);
    norm = w * a * a;
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  platform->flopCounter->add("weightedNorm2", 3 * static_cast<double>(N));

  return sqrt(norm);
}
dfloat linAlg_t::weightedNorm2Many(const dlong N,
                                   const dlong Nfields,
                                   const dlong fieldOffset,
                                   occa::memory &o_w,
                                   occa::memory &o_a,
                                   MPI_Comm _comm)
{
  if (timer)
    platform->timer.tic("dotp", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat norm = 0;
  if (N > 1 || Nfields > 1) {
    weightedNorm2ManyKernel(Nblock, N, Nfields, fieldOffset, o_w, o_a, o_scratch);

    if (serial) {
      norm = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        norm += scratch[n];
      }
    }
  }
  else {
    dfloat w, a;
    o_w.copyTo(&w, Nbytes);
    o_a.copyTo(&a, Nbytes);
    norm = w * a * a;
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  platform->flopCounter->add("weightedNorm2Many", 3 * static_cast<double>(N) * Nfields);
  return sqrt(norm);
}

// ||o_a||_w1
dfloat linAlg_t::weightedNorm1(const dlong N, occa::memory &o_w, occa::memory &o_a, MPI_Comm _comm)
{
  if (timer)
    platform->timer.tic("dotp", 1);

  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat norm = 0;
  if (N > 1) {
    weightedNorm1Kernel(Nblock, N, o_w, o_a, o_scratch);

    if (serial) {
      norm = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        norm += scratch[n];
      }
    }
  }
  else {
    dfloat w, a;
    o_w.copyTo(&w, Nbytes);
    o_a.copyTo(&a, Nbytes);
    norm = std::abs(w * a);
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  return norm;
}
dfloat linAlg_t::weightedNorm1Many(const dlong N,
                                   const dlong Nfields,
                                   const dlong fieldOffset,
                                   occa::memory &o_w,
                                   occa::memory &o_a,
                                   MPI_Comm _comm)
{
  if (timer)
    platform->timer.tic("dotp", 1);
  int Nblock = (N + blocksize - 1) / blocksize;
  const size_t Nbytes = Nblock * sizeof(dfloat);
  if (o_scratch.size() < Nbytes)
    reallocScratch(Nbytes);

  dfloat norm = 0;
  if (N > 1 || Nfields > 1) {
    weightedNorm1ManyKernel(Nblock, N, Nfields, fieldOffset, o_w, o_a, o_scratch);

    if (serial) {
      norm = *((dfloat *)o_scratch.ptr());
    }
    else {
      o_scratch.copyTo(scratch, Nbytes);
      for (dlong n = 0; n < Nblock; ++n) {
        norm += scratch[n];
      }
    }
  }
  else {
    dfloat w, a;
    o_w.copyTo(&w, Nbytes);
    o_a.copyTo(&a, Nbytes);
    norm = std::abs(w * a);
  }

  if (_comm != MPI_COMM_SELF)
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, _comm);

  if (timer)
    platform->timer.toc("dotp");

  return norm;
}

void linAlg_t::crossProduct(const dlong N,
                            const dlong fieldOffset,
                            occa::memory &o_x,
                            occa::memory &o_y,
                            occa::memory &o_z)
{
  crossProductKernel(N, fieldOffset, o_x, o_y, o_z);
}

void linAlg_t::unitVector(const dlong N, const dlong fieldOffset, occa::memory &o_v)
{
  unitVectorKernel(N, fieldOffset, o_v);
}

void linAlg_t::entrywiseMag(const dlong N,
                  const dlong Nfields,
                  const dlong fieldOffset,
                  occa::memory &o_a,
                  occa::memory &o_b)
{
  entrywiseMagKernel(N, Nfields, fieldOffset, o_a, o_b);
}