#ifndef OGS_FINDPTS_H
#define OGS_FINDPTS_H

#if !defined(MEM_H) || !defined(FINDPTS_H) || !defined(FINDPTS_LOCAL_H) || !defined(FINDPTS_EL_H) || !defined(OBBOX_H)
#warning "ogs_findpts.h" requires "mem.h", "findpts.h", "findpts_local.h", "findpts_el.h", "obbox.h"
#endif

struct eval_src_pt_3 { double r[3]; uint index, proc, el; };
struct eval_src_pt_2 { double r[2]; uint index, proc, el; };
struct eval_out_pt_3 { double out; uint index, proc; };
struct eval_out_pt_2 { double out; uint index, proc; };

void ogs_findpts_local_eval_internal_2(
  struct eval_out_pt_2 *opt, const struct eval_src_pt_2 *spt,
  const uint pn, const void *const in,
  struct findpts_local_data_2 *const gs_fd, const void *const ogs_fd_void);

void ogs_findpts_local_eval_internal_3(
  struct eval_out_pt_3 *opt, const struct eval_src_pt_3 *spt,
  const uint pn, const void *const in,
  struct findpts_local_data_3 *const gs_fd, const void *const ogs_fd_void);

void ogs_findpts_local_eval_2(
        double *const  out_base, const unsigned  out_stride,
  const uint   *const   el_base, const unsigned   el_stride,
  const double *const    r_base, const unsigned    r_stride,
  const uint pn, const void *const in,
  struct findpts_local_data_2 *const gs_fd, const void *const ogs_fd_void);

void ogs_findpts_local_eval_3(
        double *const  out_base, const unsigned  out_stride,
  const uint   *const   el_base, const unsigned   el_stride,
  const double *const    r_base, const unsigned    r_stride,
  const uint pn, const void *const in,
  struct findpts_local_data_3 *const gs_fd, const void *const ogs_fd_void);

void ogs_findpts_eval_2(
        double *const  out_base, const unsigned  out_stride,
  const uint   *const code_base, const unsigned code_stride,
  const uint   *const proc_base, const unsigned proc_stride,
  const uint   *const   el_base, const unsigned   el_stride,
  const double *const    r_base, const unsigned    r_stride,
  const uint npt,
  const void *const in, struct findpts_data_2 *const fd,
  const void *const ogs_fd);

void ogs_findpts_eval_3(
        double *const  out_base, const unsigned  out_stride,
  const uint   *const code_base, const unsigned code_stride,
  const uint   *const proc_base, const unsigned proc_stride,
  const uint   *const   el_base, const unsigned   el_stride,
  const double *const    r_base, const unsigned    r_stride,
  const uint npt,
  const void *const in, struct findpts_data_3 *const fd,
  const void *const ogs_fd);

#endif
