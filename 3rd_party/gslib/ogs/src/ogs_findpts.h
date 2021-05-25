#ifndef OGS_FINDPTS_H
#define OGS_FINDPTS_H

#if !defined(MEM_H) || !defined(FINDPTS_H) || !defined(FINDPTS_LOCAL_H) || !defined(FINDPTS_EL_H) || !defined(OBBOX_H)
#warning "ogs_findpts.h" requires "mem.h", "findpts.h", "findpts_local.h", "findpts_el.h", "obbox.h"
#endif

void ogs_findpts_local_eval_internal_2(
       double    *const out_base, const unsigned out_stride,
  const unsigned *const  el_base, const unsigned  el_stride,
  const double   *const   r_base, const unsigned   r_stride,
  const unsigned pn, const void *const in, const unsigned in_stride,
  unsigned *const n, double *const lag_data[2], unsigned lag_data_size[2]);

void ogs_findpts_local_eval_internal_3(
        double    *const out_base, const unsigned out_stride,
  const unsigned *const  el_base, const unsigned  el_stride,
  const double   *const   r_base, const unsigned   r_stride,
  const unsigned pn, const void *const in, const unsigned in_stride,
  unsigned *const n, double *const lag_data[3], unsigned lag_data_size[3]);

void ogs_findpts_local_eval_2(
        double *const out_base, const unsigned out_stride,
  const uint   *const  el_base, const unsigned  el_stride,
  const double *const   r_base, const unsigned   r_stride,
  const uint npt,
  const void *const in, struct findpts_local_data_2 *const fd);

void ogs_findpts_local_eval_3(
        double *const out_base, const unsigned out_stride,
  const uint   *const  el_base, const unsigned  el_stride,
  const double *const   r_base, const unsigned   r_stride,
  const uint npt,
  const void *const in, struct findpts_local_data_3 *const fd);

void ogs_findpts_eval_2(
        double *const  out_base, const unsigned  out_stride,
  const uint   *const code_base, const unsigned code_stride,
  const uint   *const proc_base, const unsigned proc_stride,
  const uint   *const   el_base, const unsigned   el_stride,
  const double *const    r_base, const unsigned    r_stride,
  const uint npt,
  const void *const in, struct findpts_data_2 *const fd);

void ogs_findpts_eval_3(
        double *const  out_base, const unsigned  out_stride,
  const uint   *const code_base, const unsigned code_stride,
  const uint   *const proc_base, const unsigned proc_stride,
  const uint   *const   el_base, const unsigned   el_stride,
  const double *const    r_base, const unsigned    r_stride,
  const uint npt,
  const void *const in, struct findpts_data_3 *const fd);

#endif
