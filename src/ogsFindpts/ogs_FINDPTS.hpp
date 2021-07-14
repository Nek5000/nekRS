#ifndef OGS_FINDPTS_HPP
#define OGS_FINDPTS_HPP

#include "ogs.hpp"

typedef struct {
  int D;
  void *findpts_data;
  occa::device *device;
  occa::kernel local_eval_kernel;
  occa::kernel local_kernel;
  occa::memory d_fd_local;
} ogs_findpts_t;

ogs_findpts_t *ogsFindptsSetup(
  const dlong D, MPI_Comm comm,
  const dfloat *const elx[],
  const dlong n[], const dlong nel,
  const dlong m[], const dfloat bbox_tol,
  const hlong local_hash_size, const hlong global_hash_size,
  const dlong npt_max, const dfloat newt_tol,
  occa::device *device = nullptr);
void ogsFindptsFree(ogs_findpts_t *fd);
void ogsFindpts(    dlong  *const  code_base  , const dlong  code_stride,
                    dlong  *const  proc_base  , const dlong  proc_stride,
                    dlong  *const    el_base  , const dlong    el_stride,
                    dfloat *const     r_base  , const dlong     r_stride,
                    dfloat *const dist2_base  , const dlong dist2_stride,
              const dfloat *const     x_base[], const dlong     x_stride[],
              const dlong npt, ogs_findpts_t *const fd,
              const bool use_device=true);
void ogsFindptsEval(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const code_base, const dlong code_stride,
  const dlong  *const proc_base, const dlong proc_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, const dfloat *const in, ogs_findpts_t *const fd);

void ogsFindptsEval(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const code_base, const dlong code_stride,
  const dlong  *const proc_base, const dlong proc_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, occa::memory d_in, ogs_findpts_t *const fd);

void ogsFindptsLocalEval(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, const dfloat *const in, ogs_findpts_t *const fd);

void ogsFindptsLocalEval(
  occa::memory out_base, const dlong  out_stride,
  occa::memory  el_base, const dlong   el_stride,
  occa::memory   r_base, const dlong    r_stride,
  const dlong npt, occa::memory d_in, ogs_findpts_t *const fd);



#endif
