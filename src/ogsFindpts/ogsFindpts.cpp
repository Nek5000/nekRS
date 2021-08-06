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

#include <cassert>
#include <cstdlib>
#include "ogstypes.h"
#include "ogs_FINDPTS.hpp"
#include "ogsInterface_FINDPTS.h"
#include "ogsKernels_FINDPTS.hpp"
#include "gslib.h"

// need to access internals of findpts_data structs
struct hash_data_2 {
  ulong hash_n;
  struct dbl_range bnd[2];
  double fac[2];
  uint *offset;
};
struct findpts_data_2 {
  struct crystal cr;
  struct findpts_local_data_2 local;
  struct hash_data_2 hash;
};
struct hash_data_3 {
  ulong hash_n;
  struct dbl_range bnd[3];
  double fac[3];
  uint *offset;
};
struct findpts_data_3 {
  struct crystal cr;
  struct findpts_local_data_3 local;
  struct hash_data_3 hash;
};
extern "C" {
uint findpts_local_hash_opt_size_2(struct findpts_local_hash_data_2 *p,
                               const struct obbox_2 *const obb, const uint nel,
                               const uint max_size);
uint findpts_local_hash_opt_size_3(struct findpts_local_hash_data_3 *p,
                               const struct obbox_3 *const obb, const uint nel,
                               const uint max_size);
}

// Internal routine to copy findpts_local_data to the GPU
static occa::memory ogsFindptsCopyData_2(const struct findpts_data_2 *fd,
                                         dlong nel, dlong max_hash_size,
                                         occa::device device)
{
  const findpts_local_data_2 *fd_local = &fd->local;

  //create device allocation
  dlong lag_size[2];
  for(dlong d=0;d<2;++d) lag_size[d] = gll_lag_size(fd_local->fed.n[d]);
  struct findpts_local_hash_data_2 hash_data_copy = fd_local->hd;
  dlong hd_d_size  = findpts_local_hash_opt_size_2(&hash_data_copy, fd_local->obb, nel, max_hash_size);
  dlong elx_d_size = nel*fd_local->ntot;
  dlong alloc_size =  sizeof(findpts_local_data_2)
                    + sizeof(double)*elx_d_size*2 // elx
                    + sizeof(obbox_2)*nel // obb
                    + sizeof(uint)*hd_d_size // hd.offset
                    + sizeof(double)*(fd_local->fed.n[0]+fd_local->fed.n[1]) // fed.z
                    + sizeof(double)*(lag_size[0]+lag_size[1]) // fed.lag_data
                    + sizeof(double)*(6*fd_local->fed.n[0]+6*fd_local->fed.n[1]) // fed.wtend
                    + 8*6; // room for alignment padding
                    // other structures are unused
  occa::memory dev_copy = device.malloc(alloc_size, occa::dtype::byte);
  char *dev_copy_ptr = dev_copy.ptr<char>();

  // create and fill host buffer
  char *host_copy = new char[alloc_size];
  struct findpts_local_data_2 *fd_local_copy = (struct findpts_local_data_2*)host_copy;
  *fd_local_copy = *fd_local;
  dlong working_offset = sizeof(findpts_local_data_2);

  #define SET_FIELD(field, type, size) do {\
        dlong align = 8; \
        working_offset = ((working_offset + align - 1)/align)*align; \
        fd_local_copy->field = (type*)(dev_copy_ptr+working_offset); \
        memcpy(host_copy+working_offset, fd_local->field, size*sizeof(type)); \
        working_offset += size*sizeof(type); \
      } while(0)

  SET_FIELD(elx[0], double, elx_d_size);
  SET_FIELD(elx[1], double, elx_d_size);
  SET_FIELD(obb, struct obbox_2, nel);
  SET_FIELD(hd.offset, uint, hd_d_size);
  SET_FIELD(fed.z[0], double, fd_local->fed.n[0]);
  SET_FIELD(fed.z[1], double, fd_local->fed.n[1]);
  SET_FIELD(fed.lag_data[0], double, lag_size[0]);
  SET_FIELD(fed.lag_data[1], double, lag_size[1]);
  SET_FIELD(fed.wtend[0], double, 6*fd_local->fed.n[0]);
  SET_FIELD(fed.wtend[1], double, 6*fd_local->fed.n[1]);

  #undef SET_FIELD

  // copy buffer to device
  dev_copy.copyFrom(host_copy);
  delete[] host_copy;
  return dev_copy;
}
static occa::memory ogsFindptsCopyData_3(const struct findpts_data_3 *fd,
                                         dlong nel, dlong max_hash_size,
                                         occa::device device)
{
  const findpts_local_data_3 *fd_local = &fd->local;

  //create device allocation
  dlong lag_size[3];
  for(dlong d=0;d<3;++d) lag_size[d] = gll_lag_size(fd_local->fed.n[d]);
  struct findpts_local_hash_data_3 hash_data_copy = fd_local->hd;
  dlong hd_d_size  = findpts_local_hash_opt_size_3(&hash_data_copy, fd_local->obb, nel, max_hash_size);
  dlong elx_d_size = nel*fd_local->ntot;
  dlong alloc_size =  sizeof(findpts_local_data_3)
                    + sizeof(double)*elx_d_size*3 // elx
                    + sizeof(obbox_3)*nel
                    + sizeof(uint)*hd_d_size // hd.offset
                    + sizeof(double)*(fd_local->fed.n[0]+fd_local->fed.n[1]+fd_local->fed.n[2]) // fed.z
                    + sizeof(double)*(lag_size[0]+lag_size[1]+lag_size[2]) // fed.lag_data
                    + sizeof(double)*(6*fd_local->fed.n[0]+6*fd_local->fed.n[1]+6*fd_local->fed.n[2]) // fed.wtend
                    + 8*6; // room for alignment padding
                    // other structures are unused
  occa::memory dev_copy = device.malloc(alloc_size, occa::dtype::byte);
  char *dev_copy_ptr = dev_copy.ptr<char>();

  // create and fill host buffer
  char *host_copy = new char[alloc_size];
  struct findpts_local_data_3 *fd_local_copy = (struct findpts_local_data_3*)host_copy;
  *fd_local_copy = *fd_local;
  dlong working_offset = sizeof(findpts_local_data_3);

  // ensures all fields are aligned to 8 bytes
  #define SET_FIELD(field, type, size) do {\
        dlong align = 8; \
        working_offset = ((working_offset + align - 1)/align)*align; \
        fd_local_copy->field = (type*)(dev_copy_ptr+working_offset); \
        memcpy(host_copy+working_offset, fd_local->field, size*sizeof(type)); \
        working_offset += size*sizeof(type); \
      } while(0)

  SET_FIELD(elx[0], double, elx_d_size);
  SET_FIELD(elx[1], double, elx_d_size);
  SET_FIELD(elx[2], double, elx_d_size);
  SET_FIELD(obb, struct obbox_3, nel);
  SET_FIELD(hd.offset, uint, hd_d_size);
  SET_FIELD(fed.z[0], double, fd_local->fed.n[0]);
  SET_FIELD(fed.z[1], double, fd_local->fed.n[1]);
  SET_FIELD(fed.z[2], double, fd_local->fed.n[2]);
  SET_FIELD(fed.lag_data[0], double, lag_size[0]);
  SET_FIELD(fed.lag_data[1], double, lag_size[1]);
  SET_FIELD(fed.lag_data[2], double, lag_size[2]);
  SET_FIELD(fed.wtend[0], double, 6*fd_local->fed.n[0]);
  SET_FIELD(fed.wtend[1], double, 6*fd_local->fed.n[1]);
  SET_FIELD(fed.wtend[2], double, 6*fd_local->fed.n[2]);

  #undef SET_FIELD

  // copy buffer to device
  dev_copy.copyFrom(host_copy);
  delete[] host_copy;
  return dev_copy;
}

ogs_findpts_t* ogsFindptsSetup(
  const dlong D, MPI_Comm comm,
  const dfloat* const elx[],
  const dlong n[], const dlong nel,
  const dlong m[], const dfloat bbox_tol,
  const hlong local_hash_size, const hlong global_hash_size,
  const dlong npt_max, const dfloat newt_tol,
  occa::device* device) {
  // elx, n, m have length D

  void* findpts_data;
  if (D == 2) {
    findpts_data = ogsHostFindptsSetup_2(comm, elx, n, nel, m, bbox_tol,
                                         local_hash_size, global_hash_size,
                                         npt_max, newt_tol);
  } else if (D == 3) {
    findpts_data = ogsHostFindptsSetup_3(comm, elx, n, nel, m, bbox_tol,
                                         local_hash_size, global_hash_size,
                                         npt_max, newt_tol);
  } else {
    assert(false);
  }

  ogs_findpts_t* ogs_handle = new ogs_findpts_t();
  ogs_handle->D = D;
  ogs_handle->findpts_data = findpts_data;

  if (device != nullptr) {
    ogs_handle->device = device;
    std::pair<occa::kernel, occa::kernel> kernels = ogs::initFindptsKernel(comm, *device, D, n);
    ogs_handle->local_eval_kernel = std::get<0>(kernels);
    ogs_handle->local_kernel = std::get<1>(kernels);

    // Need to copy findpts data to the
    if (D == 2) {
      ogs_handle->o_fd_local
          = ogsFindptsCopyData_2((struct findpts_data_2*)ogs_handle->findpts_data, nel, local_hash_size, *device);
    } else {
      ogs_handle->o_fd_local
          = ogsFindptsCopyData_3((struct findpts_data_3*)ogs_handle->findpts_data, nel, local_hash_size, *device);
    }
  } else {
    ogs_handle->device = nullptr;
  }

  return ogs_handle;
}


void ogsFindptsFree(ogs_findpts_t* fd) {
  if (fd->D == 2) {
    ogsHostFindptsFree_2((findpts_data_2*)fd->findpts_data);
  } else {
    ogsHostFindptsFree_3((findpts_data_3*)fd->findpts_data);
  }
  if (fd->device != nullptr) {
    // Use OCCA's reference counting to free memory and kernel objects
    fd->local_eval_kernel = occa::kernel();
    fd->local_kernel = occa::kernel();
    fd->o_fd_local = occa::memory();
  }
  delete fd;
}

void ogsFindpts(    dlong*  const  code_base  , const dlong  code_stride,
                    dlong*  const  proc_base  , const dlong  proc_stride,
                    dlong*  const    el_base  , const dlong    el_stride,
                    dfloat* const     r_base  , const dlong     r_stride,
                    dfloat* const dist2_base  , const dlong dist2_stride,
              const dfloat* const     x_base[], const dlong     x_stride[],
              const dlong npt, ogs_findpts_t* const fd,
              const bool use_device) {
  // x_base, x_stride have length D

  if (use_device) {
    if (fd->D == 2) {
      ogsDevFindpts_2 ( code_base,  code_stride,
                        proc_base,  proc_stride,
                          el_base,    el_stride,
                           r_base,     r_stride,
                       dist2_base, dist2_stride,
                           x_base,     x_stride,
                       npt, (findpts_data_2*)fd->findpts_data, fd);
    } else {
      ogsDevFindpts_3 ( code_base,  code_stride,
                        proc_base,  proc_stride,
                          el_base,    el_stride,
                           r_base,     r_stride,
                       dist2_base, dist2_stride,
                           x_base,     x_stride,
                       npt, (findpts_data_3*)fd->findpts_data, fd);
    }
  } else {
    if (fd->D == 2) {
      ogsHostFindpts_2( code_base,  code_stride,
                        proc_base,  proc_stride,
                          el_base,    el_stride,
                           r_base,     r_stride,
                       dist2_base, dist2_stride,
                           x_base,     x_stride,
                       npt, (findpts_data_2*)fd->findpts_data);
    } else {
      ogsHostFindpts_3( code_base,  code_stride,
                        proc_base,  proc_stride,
                          el_base,    el_stride,
                           r_base,     r_stride,
                       dist2_base, dist2_stride,
                           x_base,     x_stride,
                       npt, (findpts_data_3*)fd->findpts_data);
    }
  }
}

void ogsFindptsEval(
        dfloat* const  out_base, const dlong  out_stride,
  const dlong*  const code_base, const dlong code_stride,
  const dlong*  const proc_base, const dlong proc_stride,
  const dlong*  const   el_base, const dlong   el_stride,
  const dfloat* const    r_base, const dlong    r_stride,
  const dlong npt, const dfloat* const in, ogs_findpts_t *const fd) {

  if (fd->D == 2) {
    ogsHostFindptsEval_2( out_base,  out_stride,
                         code_base, code_stride,
                         proc_base, proc_stride,
                           el_base,   el_stride,
                            r_base,    r_stride,
                         npt, in, (findpts_data_2*)fd->findpts_data);
  } else {
    ogsHostFindptsEval_3( out_base,  out_stride,
                         code_base, code_stride,
                         proc_base, proc_stride,
                           el_base,   el_stride,
                            r_base,    r_stride,
                         npt, in, (findpts_data_3*)fd->findpts_data);
  }
}


void ogsFindptsEval(
        dfloat* const  out_base, const dlong  out_stride,
  const dlong*  const code_base, const dlong code_stride,
  const dlong*  const proc_base, const dlong proc_stride,
  const dlong*  const   el_base, const dlong   el_stride,
  const dfloat* const    r_base, const dlong    r_stride,
  const dlong npt, occa::memory d_in, ogs_findpts_t* const fd) {

  if (fd->D == 2) {
    ogsDevFindptsEval_2( out_base,  out_stride,
                        code_base, code_stride,
                        proc_base, proc_stride,
                          el_base,   el_stride,
                           r_base,    r_stride,
                        npt, &d_in, (findpts_data_2*)fd->findpts_data, fd);
  } else {
    ogsDevFindptsEval_3( out_base,  out_stride,
                        code_base, code_stride,
                        proc_base, proc_stride,
                          el_base,   el_stride,
                           r_base,    r_stride,
                        npt, &d_in, (findpts_data_3*)fd->findpts_data, fd);
  }
}

void ogsFindptsLocalEval(
        dfloat* const  out_base, const dlong  out_stride,
  const dlong*  const   el_base, const dlong   el_stride,
  const dfloat* const    r_base, const dlong    r_stride,
  const dlong npt, const dfloat* const in, ogs_findpts_t* const fd) {

  if (fd->D == 2) {
    ogsHostFindptsLocalEval_2(out_base,  out_stride,
                               el_base,   el_stride,
                                r_base,    r_stride,
                              npt, in, (findpts_data_2*)fd->findpts_data);
  } else {
    ogsHostFindptsLocalEval_3(out_base,  out_stride,
                               el_base,   el_stride,
                                r_base,    r_stride,
                              npt, in, (findpts_data_3*)fd->findpts_data);
  }
}

void ogsFindptsLocalEval(
  occa::memory  out_base, const dlong  out_stride,
  occa::memory   el_base, const dlong   el_stride,
  occa::memory    r_base, const dlong    r_stride,
  const dlong npt, occa::memory d_in, ogs_findpts_t* const fd) {

  if (fd->D == 2) {
    ogsDevFindptsLocalEval_2(&out_base,  out_stride,
                             & el_base,   el_stride,
                             &  r_base,    r_stride,
                             npt, &d_in, (findpts_data_2*)fd->findpts_data, fd);
  } else {
    ogsDevFindptsLocalEval_3(&out_base,  out_stride,
                             & el_base,   el_stride,
                             &  r_base,    r_stride,
                             npt, &d_in, (findpts_data_3*)fd->findpts_data, fd);
  }
}

crystal* ogsCrystalRouter(ogs_findpts_t *const fd){
  if (fd->D == 2) {
    return &((findpts_data_2*))(fd->findpts_data))->cr;
  } else {
    return &((findpts_data_3*))(fd->findpts_data))->cr;
  }
}
