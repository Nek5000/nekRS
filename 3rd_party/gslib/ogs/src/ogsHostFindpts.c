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

/* compile with C compiler (not C++) */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "gslib.h"
#include "ogs_findpts.h"

#include "ogstypes.h"

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

struct findpts_data_2 *ogsHostFindptsSetup_2(
  MPI_Comm mpi_comm,
  const dfloat *const elx[2],
  const dlong n[2], const dlong nel,
  const dlong m[2], const dfloat bbox_tol,
  const hlong local_hash_size, const hlong global_hash_size,
  const dlong npt_max, const dfloat newt_tol) {

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  struct comm gs_comm;
  comm_init(&gs_comm, mpi_comm);

  return findpts_setup_2(&gs_comm, elx, n, nel, m, bbox_tol,
                         local_hash_size, global_hash_size, npt_max, newt_tol);
}

struct findpts_data_3 *ogsHostFindptsSetup_3(
  MPI_Comm mpi_comm,
  const dfloat *const elx[3],
  const dlong n[3], const dlong nel,
  const dlong m[3], const dfloat bbox_tol,
  const hlong local_hash_size, const hlong global_hash_size,
  const dlong npt_max, const dfloat newt_tol) {

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  struct comm gs_comm;
  comm_init(&gs_comm, mpi_comm);

  return findpts_setup_3(&gs_comm, elx, n, nel, m, bbox_tol,
                         local_hash_size, global_hash_size, npt_max, newt_tol);
}


void ogsHostFindptsFree_2(struct findpts_data_2 *fd) {
  findpts_free_2(fd);
}
void ogsHostFindptsFree_3(struct findpts_data_3 *fd) {
  findpts_free_3(fd);
}

void ogsHostFindptsLagData_2(struct findpts_data_2 *const fd,
                             dfloat **lag_data, dlong *lag_data_size) {
  for (int i = 0; i < 2; ++i) {
    lag_data[i] = fd->local.fed.lag_data[i];
    lag_data_size[i] = gll_lag_size(fd->local.fed.n[i]);
  }
}
void ogsHostFindptsLagData_3(struct findpts_data_3 *const fd,
                             dfloat **lag_data, dlong *lag_data_size) {
  for (int i = 0; i < 3; ++i) {
    lag_data[i] = fd->local.fed.lag_data[i];
    lag_data_size[i] = gll_lag_size(fd->local.fed.n[i]);
  }
}

void ogsHostFindpts_2(    dlong  *const  code_base   , const dlong  code_stride   ,
                          dlong  *const  proc_base   , const dlong  proc_stride   ,
                          dlong  *const    el_base   , const dlong    el_stride   ,
                          dfloat *const     r_base   , const dlong     r_stride   ,
                          dfloat *const dist2_base   , const dlong dist2_stride   ,
                    const dfloat *const     x_base[2], const dlong     x_stride[2],
                    const dfloat npt, struct findpts_data_2 *const fd) {

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  findpts_2( code_base,  code_stride,
             proc_base,  proc_stride,
               el_base,    el_stride,
                r_base,     r_stride,
            dist2_base, dist2_stride,
                x_base,     x_stride,
            npt, fd);
}

void ogsHostFindpts_3(    dlong  *const  code_base   , const dlong  code_stride   ,
                          dlong  *const  proc_base   , const dlong  proc_stride   ,
                          dlong  *const    el_base   , const dlong    el_stride   ,
                          dfloat *const     r_base   , const dlong     r_stride   ,
                          dfloat *const dist2_base   , const dlong dist2_stride   ,
                    const dfloat *const     x_base[3], const dlong     x_stride[3],
                    const dfloat npt, struct findpts_data_3 *const fd) {

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  findpts_3( code_base,  code_stride,
             proc_base,  proc_stride,
               el_base,    el_stride,
                r_base,     r_stride,
            dist2_base, dist2_stride,
                x_base,     x_stride,
            npt, fd);
}

void ogsHostFindptsEval_2(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const code_base, const dlong code_stride,
  const dlong  *const proc_base, const dlong proc_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, const dfloat *const in, struct findpts_data_2 *const fd) {

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  findpts_eval_2( out_base,  out_stride,
                 code_base, code_stride,
                 proc_base, proc_stride,
                   el_base,   el_stride,
                    r_base,    r_stride,
                 npt, in, fd);
}

void ogsHostFindptsEval_3(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const code_base, const dlong code_stride,
  const dlong  *const proc_base, const dlong proc_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, const dfloat *const in, struct findpts_data_3 *const fd) {

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  findpts_eval_3( out_base,  out_stride,
                 code_base, code_stride,
                 proc_base, proc_stride,
                   el_base,   el_stride,
                    r_base,    r_stride,
                 npt, in, fd);
}
