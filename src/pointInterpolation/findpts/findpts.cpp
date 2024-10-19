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

// for platform
#include "platform.hpp"
#include "tuple_for_each.hpp"

#include <cstdlib>
#include <vector>
#include "ogstypes.h"
#include "findpts.hpp"
#include "gslib.h"

// local data structures to switch between run-time/compile-time sizes
struct evalSrcPt_t {
  double r[findpts::dim];
  int index, proc, el;
};

template <int N> struct evalOutPt_t {
  double out[N];
  int index, proc;
};

extern "C" {
struct hash_data_3 {
  ulong hash_n;
  struct dbl_range bnd[findpts::dim];
  double fac[findpts::dim];
  uint *offset;
};

struct findpts_dummy_ms_data {
  unsigned int *nsid;
  double *distfint;
};

struct findpts_data_3 {
  struct crystal cr;
  struct findpts_local_data_3 local;
  struct hash_data_3 hash;
  struct array savpt;
  struct findpts_dummy_ms_data fdms;
  uint fevsetup;
};

auto *gslibFindptsSetup(MPI_Comm mpi_comm,
                        const dfloat *const _elx[findpts::dim],
                        const dlong n[findpts::dim],
                        const dlong nel,
                        const dlong m[findpts::dim],
                        const dfloat bbox_tol,
                        const dlong local_hash_size,
                        const dlong global_hash_size,
                        const dlong npt_max,
                        const dfloat newt_tol,
                        const dlong nsid,
                        const dfloat *const _distfint)
{

  bool useMultiSessionSupport = _distfint != nullptr;

  struct comm gs_comm;
  comm_init(&gs_comm, mpi_comm);

  const unsigned int n_[findpts::dim] = {(unsigned int)n[0], (unsigned int)n[1], (unsigned int)n[2]};
  const unsigned int m_[findpts::dim] = {(unsigned int)m[0], (unsigned int)m[1], (unsigned int)m[2]};

  const auto Nlocal = nel * n[0] * n[1] * n[2];
  std::vector<double> x(Nlocal);
  std::vector<double> y(Nlocal);
  std::vector<double> z(Nlocal);
  const double *elx[3] = {x.data(), y.data(), z.data()};
  {
    auto _x = _elx[0];
    auto _y = _elx[1];
    auto _z = _elx[2];
    for (int i = 0; i < Nlocal; i++) {
      x[i] = _x[i];
      y[i] = _y[i];
      z[i] = _z[i];
    }
  }

  void *h;
  if (useMultiSessionSupport) {
    std::vector<double> distfint(Nlocal);
    for (int i = 0; i < Nlocal; i++) {
      distfint[i] = _distfint[i];
    }

    uint unsid = nsid;
    h = findptsms_setup_3(&gs_comm,
                          elx,
                          n_,
                          (uint)nel,
                          m_,
                          (double)bbox_tol,
                          (uint)local_hash_size,
                          (uint)global_hash_size,
                          (unsigned)npt_max,
                          (double)newt_tol,
                          &unsid,
                          distfint.data());
  } else {
    h = findpts_setup_3(&gs_comm,
                        elx,
                        n_,
                        (uint)nel,
                        m_,
                        (double)bbox_tol,
                        (uint)local_hash_size,
                        (uint)global_hash_size,
                        (unsigned)npt_max,
                        (double)newt_tol);
  }

  comm_free(&gs_comm);
  return h;
}
} // extern "C"

namespace findpts
{

void findpts_t::findptsLocal(int *const code,
                             int *const el,
                             int *const elsid,
                             dfloat *const r,
                             dfloat *const dist2,
                             dfloat *const disti,
                             const dfloat *const x,
                             const dfloat *const y,
                             const dfloat *const z,
                             const dlong *const sess,
                             const int sessionIdMatch,
                             const int pn)
{
  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.tic(timerName + "findptsLocal");
  }

  if (pn == 0) {
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "findptsLocal");
    }

    // provide 0-time tic/toc to allow global reduction later in timer reporting
    if (timerLevel != TimerLevel::None) {
      platform->timer.tic(timerName + "findptsLocal::localKernel");
      platform->timer.toc(timerName + "findptsLocal::localKernel");
    }
    return;
  }

  auto o_code = platform->deviceMemoryPool.reserve<dlong>(pn);
  auto o_el = platform->deviceMemoryPool.reserve<dlong>(pn);
  auto o_elsid = platform->deviceMemoryPool.reserve<dlong>(pn);
  auto o_r = platform->deviceMemoryPool.reserve<dfloat>(dim * pn);
  auto o_dist2 = platform->deviceMemoryPool.reserve<dfloat>(pn);
  auto o_disti = platform->deviceMemoryPool.reserve<dfloat>(pn);

  auto o_sess = platform->deviceMemoryPool.reserve<dlong>(pn);
  auto o_xint = platform->deviceMemoryPool.reserve<dfloat>(pn);
  auto o_yint = platform->deviceMemoryPool.reserve<dfloat>(pn);
  auto o_zint = platform->deviceMemoryPool.reserve<dfloat>(pn);

  o_xint.copyFrom(x);
  o_yint.copyFrom(y);
  o_zint.copyFrom(z);

  if (useMultiSessionSupport) {
    o_sess.copyFrom(sess);
  }

  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "findptsLocal::localKernel");
  }
  this->localKernel(static_cast<dlong>(useMultiSessionSupport),
                    sessionId,
                    sessionIdMatch,
                    pn,
                    this->tol,
                    o_xint,
                    o_yint,
                    o_zint,
                    o_sess,
                    this->o_x,
                    this->o_y,
                    this->o_z,
                    this->o_distfint,
                    this->o_wtend_x,
                    this->o_wtend_y,
                    this->o_wtend_z,
                    this->o_c,
                    this->o_A,
                    this->o_min,
                    this->o_max,
                    this->hash_n,
                    this->o_hashMin,
                    this->o_hashFac,
                    this->o_offset,
                    o_code,
                    o_el,
                    o_elsid,
                    o_r,
                    o_dist2,
                    o_disti);
  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "findptsLocal::localKernel");
  }

  if (pn > 0) {
    o_code.copyTo(code);
    o_el.copyTo(el);
    o_elsid.copyTo(elsid);
    o_r.copyTo(r);
    o_dist2.copyTo(dist2);
    o_disti.copyTo(disti);
  }

  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.toc(timerName + "findptsLocal");
  }
}

void findpts_t::findptsLocal(int *const code,
                             int *const el,
                             int *const elsid,
                             dfloat *const r,
                             dfloat *const dist2,
                             dfloat *const disti,
                             occa::memory o_xint,
                             occa::memory o_yint,
                             occa::memory o_zint,
                             occa::memory o_sess,
                             const int sessionIdMatch,
                             const int pn)
{
  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.tic(timerName + "findptsLocal");
  }

  if (pn == 0) {
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "findptsLocal");
    }
    // provide 0-time tic/toc to allow global reduction later in timer reporting
    if (timerLevel != TimerLevel::None) {
      platform->timer.tic(timerName + "findptsLocal::localKernel");
      platform->timer.toc(timerName + "findptsLocal::localKernel");
    }
    return;
  }

  auto o_code = platform->deviceMemoryPool.reserve<dlong>(pn);
  auto o_el = platform->deviceMemoryPool.reserve<dlong>(pn);
  auto o_elsid = platform->deviceMemoryPool.reserve<dlong>(pn);
  auto o_r = platform->deviceMemoryPool.reserve<dfloat>(dim * pn);
  auto o_dist2 = platform->deviceMemoryPool.reserve<dfloat>(pn);
  auto o_disti = platform->deviceMemoryPool.reserve<dfloat>(pn);

  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "findptsLocal::localKernel");
  }

  this->localKernel(static_cast<dlong>(useMultiSessionSupport),
                    sessionId,
                    sessionIdMatch,
                    pn,
                    this->tol,
                    o_xint,
                    o_yint,
                    o_zint,
                    o_sess,
                    this->o_x,
                    this->o_y,
                    this->o_z,
                    this->o_distfint,
                    this->o_wtend_x,
                    this->o_wtend_y,
                    this->o_wtend_z,
                    this->o_c,
                    this->o_A,
                    this->o_min,
                    this->o_max,
                    this->hash_n,
                    this->o_hashMin,
                    this->o_hashFac,
                    this->o_offset,
                    o_code,
                    o_el,
                    o_elsid,
                    o_r,
                    o_dist2,
                    o_disti);

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "findptsLocal::localKernel");
  }

  if (pn > 0) {
    o_code.copyTo(code);
    o_el.copyTo(el);
    o_elsid.copyTo(elsid);
    o_r.copyTo(r);
    o_dist2.copyTo(dist2);
    o_disti.copyTo(disti);
  }

  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.toc(timerName + "findptsLocal");
  }
}

template <typename OutputType>
void findpts_t::findptsEvalImpl(occa::memory &o_out,
                                dlong findPtsDataOffset,
                                data_t *findPtsData,
                                const dlong npt,
                                const int nFields,
                                const dlong inputOffset,
                                const dlong outputOffset,
                                const occa::memory &o_in,
                                hashData_t &hash,
                                crystal &cr)
{
  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.tic(timerName + "findptsEvalImpl");
  }

  static std::vector<dfloat> out;
  if (out.size() < nFields * outputOffset) {
    constexpr int growthFactor = 2;
    out.resize(growthFactor * nFields * outputOffset);
  }

  // evaluate local points
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "findptsEvalImpl::eval local points");
  }
  if (npt > 0) {
  platform->device.occaDevice().setStream(localEvalStream);

    this->localEvalMaskKernel(npt,
                              nFields,
                              inputOffset,
                              outputOffset,
                              this->rank,
                              this->o_proc + findPtsDataOffset,
                              this->o_code + findPtsDataOffset,
                              this->o_el + findPtsDataOffset,
                              this->o_r + dim * findPtsDataOffset,
                              o_in,
                              o_out);

    o_out.copyTo(out.data(), nFields * outputOffset, 0, "async: true");
  platform->device.occaDevice().setStream(defaultStream);

  }
  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "findptsEvalImpl::eval local points");
  }

  // transfer non-local (found on a remote rank) points
  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.tic(timerName + "findptsEvalImpl::copy data to target");
  }

  if (findPtsData->updateCache) {
    const int *code = findPtsData->code_base + findPtsDataOffset;
    const int *proc = findPtsData->proc_base + findPtsDataOffset;
    const int *el = findPtsData->el_base + findPtsDataOffset;
    const dfloat *r = findPtsData->r_base + dim * findPtsDataOffset;

    int numSend = 0;
    for (int index = 0; index < npt; ++index) {
      numSend += (code[index] != CODE_NOT_FOUND && proc[index] != this->rank);
    }

    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.tic(timerName + "findptsEvalImpl::copy data to target::sarray_transfer");
    }

    static struct array src;
    array_reserve(evalSrcPt_t, &src, numSend);
    auto spt = (evalSrcPt_t *)src.ptr;

    int cnt = 0;
    for (int index = 0; index < npt; ++index) {
      if (code[index] != CODE_NOT_FOUND && proc[index] != this->rank) {
        for (int d = 0; d < dim; ++d) {
          spt[cnt].r[d] = r[index * dim + d];
        }
        spt[cnt].index = index;
        spt[cnt].proc = proc[index];
        spt[cnt].el = el[index];

        cnt++;
      }
    }
    src.n = cnt;

    sarray_transfer(evalSrcPt_t, &src, proc, 1, &cr);

    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "findptsEvalImpl::copy data to target::sarray_transfer");
    }

    // update cache
    {
      const auto n = src.n;
      findPtsData->cache.index.resize(n);
      findPtsData->cache.proc.resize(n);

      findPtsData->cache.o_el.free();
      findPtsData->cache.o_r.free();
      findPtsData->cache.o_el = platform->deviceMemoryPool.reserve<dlong>(n);
      findPtsData->cache.o_r = platform->deviceMemoryPool.reserve<dfloat>(dim * n);

      auto r = platform->memoryPool.reserve<dfloat>(findPtsData->cache.o_r.size());
      auto rPtr = r.ptr<dfloat>();
      auto el = platform->memoryPool.reserve<dlong>(findPtsData->cache.o_el.size());
      auto elPtr = el.ptr<dlong>();

      auto spt = (evalSrcPt_t *)src.ptr;
      for (int i = 0; i < n; i++) {
        for (int d = 0; d < dim; ++d) {
          rPtr[dim * i + d] = spt[i].r[d];
        }
        elPtr[i] = spt[i].el;
        findPtsData->cache.index[i] = spt[i].index;
        findPtsData->cache.proc[i] = spt[i].proc;
      }
      findPtsData->cache.o_r.copyFrom(r);
      findPtsData->cache.o_el.copyFrom(el);
      findPtsData->updateCache = false;

    }
  }

  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.toc(timerName + "findptsEvalImpl::copy data to target");
  }

  // evaluate non-local points
  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.tic(timerName + "findptsEvalImpl::eval non-local points");
  }

  static struct array outpt;
  {
    const dlong n = findPtsData->cache.index.size();
    array_reserve(OutputType, &outpt, n);
    outpt.n = n; 

    auto o_tmp = platform->deviceMemoryPool.reserve<dfloat>(nFields * n);
    const dlong offset = o_tmp.size() / nFields;

    if (timerLevel != TimerLevel::None) {
      platform->timer.tic(timerName + "findptsEvalImpl::eval non-local points::kernel");
    }

    if (n) { 
      this->localEvalKernel(n, 
                            nFields, 
                            inputOffset, 
                            offset, 
                            findPtsData->cache.o_el, 
                            findPtsData->cache.o_r, 
                            o_in, 
                            o_tmp);
    }
 
    if (timerLevel != TimerLevel::None) {
      platform->timer.toc(timerName + "findptsEvalImpl::eval non-local points::kernel");
    }
 
    auto tmp = platform->memoryPool.reserve<dfloat>(o_tmp.size());
    o_tmp.copyTo(tmp);

    auto opt = (OutputType *)outpt.ptr;
    auto tmpPtr = tmp.ptr<dfloat>();
    for (dlong i = 0; i < n; i++) { 
      for (int field = 0; field < nFields; ++field) {
        opt[i].out[field] = tmpPtr[i + field * offset];
      }
      opt[i].index = findPtsData->cache.index[i];
      opt[i].proc = findPtsData->cache.proc[i];
    }

    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.tic(timerName + "findptsEvalImpl::eval non-local points::sarray_transfer");
    }

    sarray_transfer(OutputType, &outpt, proc, 1, &cr);

    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "findptsEvalImpl::eval non-local points::sarray_transfer");
    }
  }

  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.toc(timerName + "findptsEvalImpl::eval non-local points");
  }

  // copy to user buffer
  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.tic(timerName + "findptsEvalImpl::copy results");
  }

  {
    platform->device.occaDevice().setStream(localEvalStream);
    platform->device.finish();
    platform->device.occaDevice().setStream(defaultStream);

    auto opt = (OutputType *)outpt.ptr;
    for (int i = 0; i < outpt.n; i++) {
      for (int field = 0; field < nFields; ++field) {
        out[opt[i].index + outputOffset * field] = opt[i].out[field];
      }
    }

    if (outputOffset) {
      o_out.copyFrom(out.data(), nFields * outputOffset);
    }
  }

  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.toc(timerName + "findptsEvalImpl::copy results");
    platform->timer.toc(timerName + "findptsEvalImpl");
  }
}

extern "C" {
uint hash_opt_size_3(struct findpts_local_hash_data_3 *p,
                     const struct obbox_3 *const obb,
                     const uint nel,
                     const uint max_size);
}

dlong getHashSize(const struct findpts_data_3 *fd, dlong nel, dlong max_hash_size)
{
  const findpts_local_data_3 *fd_local = &fd->local;
  auto hash_data_copy = fd_local->hd;
  return hash_opt_size_3(&hash_data_copy, fd_local->obb, nel, max_hash_size);
}

findpts_t::findpts_t(MPI_Comm comm,
                     const dfloat *const x,
                     const dfloat *const y,
                     const dfloat *const z,
                     const dlong Nq,
                     const dlong Nelements,
                     const dlong m,
                     const dfloat bbox_tol,
                     const dlong local_hash_size,
                     const dlong global_hash_size,
                     const dlong npt_max,
                     const dfloat newt_tol)
{
  findpts_t(comm,
            x,
            y,
            z,
            Nq,
            Nelements,
            m,
            bbox_tol,
            local_hash_size,
            global_hash_size,
            npt_max,
            newt_tol,
            0,
            nullptr);
}

findpts_t::findpts_t(MPI_Comm comm,
                     const dfloat *const x,
                     const dfloat *const y,
                     const dfloat *const z,
                     const dlong Nq,
                     const dlong Nelements,
                     const dlong m,
                     const dfloat bbox_tol,
                     const dlong local_hash_size,
                     const dlong global_hash_size,
                     const dlong npt_max,
                     const dfloat newt_tol,
                     const dlong sessionId_,
                     const dfloat *const distfint)
{
  sessionId = sessionId_;
  useMultiSessionSupport = distfint != nullptr;

  const dlong Nlocal = Nq * Nq * Nq * Nelements;

  const dfloat *elx[dim] = {x, y, z};
  const int n[dim] = {Nq, Nq, Nq};
  const int ms[dim] = {m, m, m};

  defaultStream = platform->device.occaDevice().getStream();
  localEvalStream = platform->device.occaDevice().createStream();

  if (platform->options.compareArgs("ENABLE FINDPTS DETAILED TIMER", "TRUE")) {
    this->timerLevel = TimerLevel::Detailed;
  }

  this->_findptsData = gslibFindptsSetup(comm,
                                         elx,
                                         n,
                                         Nelements,
                                         ms,
                                         bbox_tol,
                                         local_hash_size,
                                         global_hash_size,
                                         npt_max,
                                         newt_tol,
                                         sessionId_,
                                         distfint);

  auto findptsData = (findpts_data_3 *)this->_findptsData;

  this->comm = comm;
  MPI_Comm_rank(comm, &this->rank);

  this->tol = findptsData->local.tol;
  this->hash = &findptsData->hash;
  this->cr = &findptsData->cr;

  if (x != nullptr) {
    this->o_x = platform->device.malloc<dfloat>(Nlocal);
    this->o_y = platform->device.malloc<dfloat>(Nlocal);
    this->o_z = platform->device.malloc<dfloat>(Nlocal);
    if (useMultiSessionSupport) {
      this->o_distfint = platform->device.malloc<dfloat>(Nlocal);
    }

    this->o_x.copyFrom(x, Nlocal);
    this->o_y.copyFrom(y, Nlocal);
    this->o_z.copyFrom(z, Nlocal);
    if (useMultiSessionSupport) {
      this->o_distfint.copyFrom(distfint, Nlocal);
    }
    std::vector<dfloat> c(dim * Nelements, 0.0);
    std::vector<dfloat> A(dim * dim * Nelements, 0.0);
    std::vector<dfloat> minBound(dim * Nelements, 0.0);
    std::vector<dfloat> maxBound(dim * Nelements, 0.0);

    for (int e = 0; e < Nelements; ++e) {
      auto box = findptsData->local.obb[e];

      c[dim * e + 0] = box.c0[0];
      c[dim * e + 1] = box.c0[1];
      c[dim * e + 2] = box.c0[2];

      minBound[dim * e + 0] = box.x[0].min;
      minBound[dim * e + 1] = box.x[1].min;
      minBound[dim * e + 2] = box.x[2].min;

      maxBound[dim * e + 0] = box.x[0].max;
      maxBound[dim * e + 1] = box.x[1].max;
      maxBound[dim * e + 2] = box.x[2].max;

      for (int i = 0; i < 9; ++i) {
        A[9 * e + i] = box.A[i];
      }
    }

    this->o_c = platform->device.malloc<dfloat>(c.size());
    this->o_A = platform->device.malloc<dfloat>(A.size());
    this->o_min = platform->device.malloc<dfloat>(minBound.size());
    this->o_max = platform->device.malloc<dfloat>(maxBound.size());

    this->o_c.copyFrom(c.data(), c.size());
    this->o_A.copyFrom(A.data(), A.size());
    this->o_min.copyFrom(minBound.data(), minBound.size());
    this->o_max.copyFrom(maxBound.data(), maxBound.size());
  }

  auto hash = findptsData->local.hd;
  dfloat hashMin[dim];
  dfloat hashFac[dim];
  for (int d = 0; d < dim; ++d) {
    hashMin[d] = hash.bnd[d].min;
    hashFac[d] = hash.fac[d];
  }
  this->hash_n = hash.hash_n;
  this->o_hashMin = platform->device.malloc<dfloat>(dim);
  this->o_hashFac = platform->device.malloc<dfloat>(dim);
  this->o_hashMin.copyFrom(hashMin, dim);
  this->o_hashFac.copyFrom(hashFac, dim);

  std::string orderSuffix = "_" + std::to_string(Nq - 1);

  this->localEvalKernel = platform->kernelRequests.load("findptsLocalEval" + orderSuffix);
  this->localEvalMaskKernel = platform->kernelRequests.load("findptsLocalEvalMask" + orderSuffix);
  this->localKernel = platform->kernelRequests.load("findptsLocal" + orderSuffix);

  this->o_wtend_x = platform->device.malloc<dfloat>(6 * Nq);
  this->o_wtend_y = platform->device.malloc<dfloat>(6 * Nq);
  this->o_wtend_z = platform->device.malloc<dfloat>(6 * Nq);
  this->o_wtend_x.copyFrom(findptsData->local.fed.wtend[0], 6 * Nq);
  this->o_wtend_y.copyFrom(findptsData->local.fed.wtend[1], 6 * Nq);
  this->o_wtend_z.copyFrom(findptsData->local.fed.wtend[2], 6 * Nq);

  const auto hd_d_size = getHashSize(findptsData, Nelements, local_hash_size);

  std::vector<dlong> offsets(hd_d_size, 0);
  for (dlong i = 0; i < hd_d_size; ++i) {
    offsets[i] = findptsData->local.hd.offset[i];
  }
  this->o_offset = platform->device.malloc<dlong>(offsets.size());
  this->o_offset.copyFrom(offsets.data(), offsets.size());
}

findpts_t::~findpts_t()
{
  auto *findptsData = (findpts_data_3 *)this->_findptsData;
  if (useMultiSessionSupport) {
    findptsms_free_3(findptsData);
  } else {
    findpts_free_3(findptsData);
  }
}

static slong lfloor(dfloat x)
{
  return floor(x);
}

static ulong hash_index_aux(dfloat low, dfloat fac, ulong n, dfloat x)
{
  const slong i = lfloor((x - low) * fac);
  return i < 0 ? 0 : (n - 1 < (ulong)i ? n - 1 : (ulong)i);
}

static ulong hash_index_3(const hashData_t *p, const dfloat x[dim])
{
  const ulong n = p->hash_n;
  return (hash_index_aux(p->bnd[2].min, p->fac[2], n, x[2]) * n +
          hash_index_aux(p->bnd[1].min, p->fac[1], n, x[1])) *
             n +
         hash_index_aux(p->bnd[0].min, p->fac[0], n, x[0]);
}

struct srcPt_t {
  dfloat x[dim];
  int index, proc, sessionId;
};

struct outPt_t {
  dfloat r[dim], dist2, disti;
  int index, code, el, proc, elsid;
};

void findpts_t::find(data_t *const findPtsData,
                     const occa::memory &o_xintIn,
                     const occa::memory &o_yintIn,
                     const occa::memory &o_zintIn,
                     const dlong npt)
{
  occa::memory o_NULL;
  this->find(findPtsData, o_xintIn, o_yintIn, o_zintIn, o_NULL, 0, npt);
}

void findpts_t::find(data_t *const findPtsData,
                     const occa::memory &o_xint,
                     const occa::memory &o_yint,
                     const occa::memory &o_zint,
                     const occa::memory &o_session,
                     const dlong sessionIdMatch,
                     const dlong npt)
{
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "find");
  }

  static std::vector<dfloat> x_base;
  static std::vector<dfloat> y_base;
  static std::vector<dfloat> z_base;
  static std::vector<int> session;
  static std::vector<dfloat> disti;
  static std::vector<int> elsid;

  static std::vector<int> codeArr;
  static std::vector<int> elArr;
  static std::vector<dfloat> rArr;
  static std::vector<dfloat> dist2Arr;
  static std::vector<dfloat> x0;
  static std::vector<dfloat> x1;
  static std::vector<dfloat> x2;
  static std::vector<int> sessArr;
  static std::vector<dfloat> distiArr;
  static std::vector<int> elsidArr;

  if (x_base.size() < npt) {
    constexpr int growthFactor = 2;
    x_base.resize(growthFactor * npt);
    y_base.resize(growthFactor * npt);
    z_base.resize(growthFactor * npt);
    session.resize(growthFactor * npt);
    disti.resize(growthFactor * npt);
    elsid.resize(growthFactor * npt);
  }

  platform->timer.tic(timerName + "find::initial copy op");
  if (npt) {
    o_xint.copyTo(x_base.data(), npt);
    o_yint.copyTo(y_base.data(), npt);
    o_zint.copyTo(z_base.data(), npt);
    if (useMultiSessionSupport) {
      o_session.copyTo(session.data(), npt);
    } else {
      std::fill(session.begin(), session.end(), 0);
    }
  }
  platform->timer.toc(timerName + "find::initial copy op");

  int *const code_base = findPtsData->code_base;
  int *const proc_base = findPtsData->proc_base;
  int *const el_base = findPtsData->el_base;
  dfloat *const r_base = findPtsData->r_base;
  dfloat *const dist2_base = findPtsData->dist2_base;
  hashData_t &hash = *this->hash;
  crystal &cr = *this->cr;
  const int np = cr.comm.np, id = cr.comm.id;
  struct array hash_pt, srcPt_t, outPt_t;

  /* look locally first */
  const auto timerNameSave = timerName;
  timerName = timerName + "find::";

  findptsLocal(code_base,
               el_base,
               elsid.data(),
               r_base,
               dist2_base,
               disti.data(),
               o_xint,
               o_yint,
               o_zint,
               o_session,
               sessionIdMatch,
               npt);
  timerName = timerNameSave;

  /* send unfound and border points to global hash cells */
  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.tic(timerName + "find::unfound");
  }
  {
    int *code = code_base, *proc = proc_base;
    const dfloat *xp[dim];
    struct srcPt_t *pt;

    xp[0] = x_base.data();
    xp[1] = y_base.data();
    xp[2] = z_base.data();

    array_init(struct srcPt_t, &hash_pt, npt);
    pt = (struct srcPt_t *)hash_pt.ptr;

    dfloat x[dim];

    for (int index = 0; index < npt; ++index) {
      for (int d = 0; d < dim; ++d) {
        x[d] = *xp[d];
      }
      *proc = id;
      if ((*code != CODE_INTERNAL) || useMultiSessionSupport) {
        const auto hi = hash_index_3(&hash, x);
        for (int d = 0; d < dim; ++d) {
          pt->x[d] = x[d];
        }
        pt->index = index;
        pt->proc = hi % np;
        pt->sessionId = session[index];
        ++pt;
      }
      for (int d = 0; d < dim; ++d) {
        xp[d]++;
      }
      code++;
      proc++;
    }
    hash_pt.n = pt - (struct srcPt_t *)hash_pt.ptr;
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.tic(timerName + "find::unfound::sarray_transfer");
    }
    sarray_transfer(struct srcPt_t, &hash_pt, proc, 1, &cr);
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "find::unfound::sarray_transfer");
    }
  }

  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.toc(timerName + "find::unfound");
    platform->timer.tic(timerName + "find::send unfound");
  }

  /* look up points in hash cells, route to possible procs */
  {
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.tic(timerName + "find::send unfound::compute hash");
    }
    const unsigned int *const hash_offset = hash.offset;
    int count = 0, *proc, *proc_p;
    const struct srcPt_t *p = (struct srcPt_t *)hash_pt.ptr, *const pe = p + hash_pt.n;
    struct srcPt_t *q;
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.tic(timerName + "find::send unfound::compute hash::get count");
    }
    for (; p != pe; ++p) {
      const int hi = hash_index_3(&hash, p->x) / np;
      const int i = hash_offset[hi], ie = hash_offset[hi + 1];
      count += ie - i;
    }
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "find::send unfound::compute hash::get count");
      platform->timer.tic(timerName + "find::send unfound::compute hash::do allocations");
    }
    proc = tmalloc(int, count);
    proc_p = proc;
    array_init(struct srcPt_t, &srcPt_t, count), q = (struct srcPt_t *)srcPt_t.ptr;
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "find::send unfound::compute hash::do allocations");
      platform->timer.tic(timerName + "find::send unfound::compute hash::doubly nested loop");
    }
    p = (struct srcPt_t *)hash_pt.ptr;
    for (; p != pe; ++p) {
      const int hi = hash_index_3(&hash, p->x) / np;
      int i = hash_offset[hi];
      const int ie = hash_offset[hi + 1];
      for (; i != ie; ++i) {
        const int pp = hash_offset[i];
        if (pp == p->proc) {
          continue; /* don't send back to source proc */
        }
        *proc_p++ = pp;
        *q++ = *p;
      }
    }
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "find::send unfound::compute hash::doubly nested loop");
    }
    array_free(&hash_pt);
    srcPt_t.n = proc_p - proc;
#ifdef DIAGNOSTICS
    printf("(proc %u) hashed; routing %u/%u\n", id, (int)srcPt_t.n, count);
#endif
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "find::send unfound::compute hash");
      platform->timer.tic(timerName + "find::send unfound::sarray_transfer_ext");
    }

    sarray_transfer_ext(struct srcPt_t, &srcPt_t, reinterpret_cast<unsigned int *>(proc), sizeof(int), &cr);

    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "find::send unfound::sarray_transfer_ext");
    }
    free(proc);
  }

  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.toc(timerName + "find::send unfound");
    platform->timer.tic(timerName + "find::send back");
  }

  /* look for other procs' points, send back */
  {
    int n = srcPt_t.n;
    const struct srcPt_t *spt;
    struct outPt_t *opt;
    array_init(struct outPt_t, &outPt_t, n);
    outPt_t.n = n;
    spt = (struct srcPt_t *)srcPt_t.ptr;
    opt = (struct outPt_t *)outPt_t.ptr;
    for (; n; --n, ++spt, ++opt) {
      opt->index = spt->index;
      opt->proc = spt->proc;
    }
    spt = (struct srcPt_t *)srcPt_t.ptr;
    opt = (struct outPt_t *)outPt_t.ptr;

    // resize result buffers
    if (codeArr.size() < srcPt_t.n) {
      constexpr int growthFactor = 2;
      codeArr.resize(growthFactor * srcPt_t.n);
      elArr.resize(growthFactor * srcPt_t.n);
      rArr.resize(growthFactor * dim * srcPt_t.n);
      dist2Arr.resize(growthFactor * srcPt_t.n);

      x0.resize(growthFactor * srcPt_t.n);
      x1.resize(growthFactor * srcPt_t.n);
      x2.resize(growthFactor * srcPt_t.n);
      sessArr.resize(growthFactor * srcPt_t.n);
      distiArr.resize(growthFactor * srcPt_t.n);
      elsidArr.resize(growthFactor * srcPt_t.n);
    }

    for (int point = 0; point < srcPt_t.n; ++point) {
      x0[point] = spt[point].x[0];
      x1[point] = spt[point].x[1];
      x2[point] = spt[point].x[2];
      sessArr[point] = spt[point].sessionId;
    }

    const auto timerNameSave = timerName;
    timerName = timerName + "find::send back::";
    findptsLocal(codeArr.data(),
                 elArr.data(),
                 elsidArr.data(),
                 rArr.data(),
                 dist2Arr.data(),
                 distiArr.data(),
                 x0.data(),
                 x1.data(),
                 x2.data(),
                 sessArr.data(),
                 sessionIdMatch,
                 srcPt_t.n);
    timerName = timerNameSave;

    // unpack arrays into opt
    for (int point = 0; point < srcPt_t.n; point++) {
      opt[point].code = codeArr[point];
      opt[point].el = elArr[point];
      opt[point].elsid = elsidArr[point];
      opt[point].dist2 = dist2Arr[point];
      opt[point].disti = distiArr[point];
      for (int d = 0; d < dim; ++d) {
        opt[point].r[d] = rArr[dim * point + d];
      }
    }

    array_free(&srcPt_t);

    /* group by code to eliminate unfound points */
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.tic(timerName + "find::send back::sarray_sort");
    }

    sarray_sort(struct outPt_t, opt, outPt_t.n, code, 0, &cr.data);

    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "find::send back::sarray_sort");
    }
    n = outPt_t.n;
    while (n && opt[n - 1].code == CODE_NOT_FOUND) {
      --n;
    }
    outPt_t.n = n;
#ifdef DIAGNOSTICS
    printf("(proc %u) sending back %u found points\n", id, (int)outPt_t.n);
#endif
    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.tic(timerName + "find::send back::sarray_transfer");
    }

    sarray_transfer(struct outPt_t, &outPt_t, proc, 1, &cr);

    if (timerLevel == TimerLevel::Detailed) {
      platform->timer.toc(timerName + "find::send back::sarray_transfer");
    }
  }

  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.toc(timerName + "find::send back");
    platform->timer.tic(timerName + "find::merge");
  }

  /* merge remote results with user data */
  {
    int n = outPt_t.n;
    struct outPt_t *opt;
    if (useMultiSessionSupport) {
      sarray_sort_2(struct outPt_t,
                    outPt_t.ptr,
                    outPt_t.n,
                    index,
                    0,
                    elsid,
                    0,
                    &cr.data); /* sort by pt and session of donor element */
      uint oldindex, nextindex;
      uint oldelsid, nextelsid;
      uint istart, ioriginator;
      oldindex = 0;
      istart = 0;
      ioriginator = 0;
      double asdisti, asdist2, asr[findpts::dim];
      uint ascode, aselsid, asproc, asel; /*winner of all session */

      double csdisti, csdist2, csr[findpts::dim];
      uint cscode, cselsid, csproc, csel; /*winner of current session */

      for (opt = (struct outPt_t *)outPt_t.ptr; n; --n, ++opt) {
        const uint index = opt->index;
        nextindex = n > 1 ? (opt + 1)->index : index;
        nextelsid = n > 1 ? (opt + 1)->elsid : opt->elsid;

        if (index != oldindex || n == outPt_t.n) { /* initialize overall winner for each pt */
          oldindex = index;
          asdisti = -std::numeric_limits<double>::max();
          oldelsid = 0;
          istart = 0;
          ioriginator = 0;
          if (code_base[index] != CODE_NOT_FOUND) {
            ioriginator = 1;
          }
        }
        if (opt->elsid != oldelsid || istart == 0) { /* initialize winner for current session */
          istart = 1;
          oldelsid = opt->elsid;
          if (ioriginator == 1 && elsid[index] == opt->elsid) { /* if the originating session found the pt */
            csdisti = disti[index];
            csdist2 = dist2_base[index];
            cscode = code_base[index];
            cselsid = elsid[index];
            csproc = proc_base[index];
            csel = el_base[index];
            for (int d = 0; d < findpts::dim; ++d) {
              csr[d] = r_base[d];
            }
            ioriginator = 0;
          } else {
            cscode = CODE_NOT_FOUND;
          }
        }

        if (cscode != CODE_INTERNAL) {
          if (cscode == CODE_NOT_FOUND || opt->code == CODE_INTERNAL || opt->dist2 < csdist2) {
            csdisti = opt->disti;
            csdist2 = opt->dist2;
            cscode = opt->code;
            cselsid = opt->elsid;
            csproc = opt->proc;
            csel = opt->el;
            for (int d = 0; d < findpts::dim; ++d) {
              csr[d] = opt->r[d];
            }
          }
        }

        if (n == 1 || opt->elsid != nextelsid || index != nextindex) { /* update overall winner */
          if (csdisti >= asdisti) {
            asdisti = csdisti;
            asdist2 = csdist2;
            ascode = cscode;
            aselsid = cselsid;
            asproc = csproc;
            asel = csel;
            for (int d = 0; d < findpts::dim; ++d) {
              asr[d] = csr[d];
            }
          }
          if (index != nextindex || n == 1) { /* copy overall winner to output array */
            if (!(ioriginator == 1 && asdisti < disti[index])) {
              disti[index] = asdisti;
              dist2_base[index] = asdist2;
              code_base[index] = ascode;
              elsid[index] = aselsid;
              proc_base[index] = asproc;
              el_base[index] = asel;
              for (int d = 0; d < findpts::dim; ++d) {
                r_base[dim * index + d] = asr[d];
              }
            }
          }
        }
      }
    } else {
      for (opt = (struct outPt_t *)outPt_t.ptr; n; --n, ++opt) {
        const int index = opt->index;
        if (code_base[index] == CODE_INTERNAL) {
          continue;
        }
        if (code_base[index] == CODE_NOT_FOUND || opt->code == CODE_INTERNAL ||
            opt->dist2 < dist2_base[index]) {
          for (int d = 0; d < dim; ++d) {
            r_base[dim * index + d] = opt->r[d];
          }
          dist2_base[index] = opt->dist2;
          proc_base[index] = opt->proc;
          el_base[index] = opt->el;
          code_base[index] = opt->code;
        }
      }
    }
    array_free(&outPt_t);
  }
  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.toc(timerName + "find::merge");
    platform->timer.tic(timerName + "find::copy to device");
  }

  if (npt) {
    if (o_code.size() < npt) {
      if (o_code.size()) {
        o_code.free();
        o_el.free();
        o_r.free();
        o_proc.free();
      }
      o_code = platform->device.malloc<dlong>(npt);
      o_el = platform->device.malloc<dlong>(npt);
      o_r = platform->device.malloc<dfloat>(npt * dim);
      o_proc = platform->device.malloc<dlong>(npt);
    }

    o_code.copyFrom(code_base, npt);
    o_el.copyFrom(el_base, npt);
    o_r.copyFrom(r_base, npt * dim);
    o_proc.copyFrom(proc_base, npt);
  }
  if (timerLevel == TimerLevel::Detailed) {
    platform->timer.toc(timerName + "find::copy to device");
  }
  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "find");
  }
}

void findpts_t::find(data_t *const findPtsData,
                     const dfloat *const x_base,
                     const dfloat *const y_base,
                     const dfloat *const z_base,
                     const dlong npt)
{
  this->find(findPtsData, x_base, y_base, z_base, nullptr, 0, npt);
}

void findpts_t::find(data_t *const findPtsData,
                     const dfloat *const x_base,
                     const dfloat *const y_base,
                     const dfloat *const z_base,
                     const dlong *const session,
                     const dlong sessionIdMatch,
                     const dlong npt)
{
  occa::memory o_xint, o_yint, o_zint, o_session;
  if (npt > 0) {
    o_xint = platform->deviceMemoryPool.reserve<dfloat>(npt);
    o_yint = platform->deviceMemoryPool.reserve<dfloat>(npt);
    o_zint = platform->deviceMemoryPool.reserve<dfloat>(npt);
    o_session = platform->deviceMemoryPool.reserve<dlong>(npt);
  }

  o_xint.copyFrom(x_base, npt);
  o_yint.copyFrom(y_base, npt);
  o_zint.copyFrom(z_base, npt);
  if (session) {
    o_session.copyFrom(session, npt);
  }

  this->find(findPtsData, o_xint, o_yint, o_zint, o_session, sessionIdMatch, npt);

  if (npt > 0) {
    o_xint.free();
    o_yint.free();
    o_zint.free();
    o_session.free();
  }
}

void findpts_t::eval(const dlong npt, const occa::memory &o_in, data_t *findPtsData, occa::memory &o_out)
{
  this->eval(npt, 0, 1, 0, npt, o_in, findPtsData, o_out);
}

void findpts_t::eval(const dlong npt,
                     const dlong findPtsOffset,
                     const dlong nFields,
                     const dlong inputOffset,
                     const dlong outputOffset,
                     const occa::memory &o_in,
                     data_t *findPtsData,
                     occa::memory &o_out)
{
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic(timerName + "eval");
  }

  const auto timerNameSave = timerName;
  timerName = timerName + "eval::";

  auto fieldSizesTuple = n_tuple<int, findpts_t::maxFields>{};
  tuple_for_each(fieldSizesTuple, [&](auto T) {
    if (nFields != decltype(T)::value) {
      return;
    }
    findptsEvalImpl<evalOutPt_t<decltype(T)::value>>(o_out,
                                                     findPtsOffset,
                                                     findPtsData,
                                                     npt,
                                                     nFields,
                                                     inputOffset,
                                                     outputOffset,
                                                     o_in,
                                                     *this->hash,
                                                     *this->cr);
  });

  nekrsCheck(nFields < 1 || nFields > findpts_t::maxFields,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "Error: nFields = %d is not supported. nFields must be between 1 and %d.",
             nFields,
             findpts_t::maxFields);

  timerName = timerNameSave;

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc(timerName + "eval");
  }
}

crystal *findpts_t::crystalRouter()
{
  return this->cr;
}

} // namespace findpts
