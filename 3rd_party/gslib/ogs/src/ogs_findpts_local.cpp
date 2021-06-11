
#include "ogstypes.h"
#include "ogs.hpp"

#define   AT(T,var,i)   \
        (T*)(      (char*)var##_base   +(i)*var##_stride   )
#define  CAT(T,var,i) \
  (const T*)((const char*)var##_base   +(i)*var##_stride   )
#define CATD(T,var,i,d) \
  (const T*)((const char*)var##_base[d]+(i)*var##_stride[d])

extern "C" {

  // GSLIB structs
  // depends on the assertion in ogsHostFindptsSetup that dlong == uint
  struct eval_src_pt_3 { double r[3]; dlong index, proc, el; };
  struct eval_src_pt_2 { double r[2]; dlong index, proc, el; };
  struct eval_out_pt_3 { double out; dlong index, proc; };
  struct eval_out_pt_2 { double out; dlong index, proc; };

void ogs_findpts_local_eval_internal_2(
  eval_out_pt_2 *opt, const eval_src_pt_2 *spt,
  const dlong pn, const void *const in, const dlong in_stride,
  dlong *const n, dfloat *const lag_data[2], dlong lag_data_size[2],
  const void *const fd_void)
{
  if (pn == 0) return;

  struct ogs_findpts_t *fd = (struct ogs_findpts_t*)fd_void;

  const dlong nr=n[0],ns=n[1];
  occa::device device = *fd->device;

  occa::memory workspace = device.malloc((sizeof(struct eval_out_pt_2)+sizeof(struct eval_src_pt_2))*pn,
                                         occa::dtype::byte);
  occa::memory d_out_pt = workspace; workspace += sizeof(struct eval_out_pt_2)*pn;
  occa::memory d_src_pt = workspace;
  d_src_pt.copyFrom(spt, sizeof(struct eval_src_pt_2)*pn);

  occa::memory d_out_out_base   = d_out_pt + offsetof(struct eval_out_pt_2, out);
  occa::memory d_out_proc_base  = d_out_pt + offsetof(struct eval_out_pt_2, proc);
  occa::memory d_out_index_base = d_out_pt + offsetof(struct eval_out_pt_2, index);
  occa::memory d_src_el_base    = d_src_pt + offsetof(struct eval_src_pt_2, el);
  occa::memory d_src_r_base     = d_src_pt + offsetof(struct eval_src_pt_2, r);
  occa::memory d_src_proc_base  = d_src_pt + offsetof(struct eval_src_pt_2, proc);
  occa::memory d_src_index_base = d_src_pt + offsetof(struct eval_src_pt_2, index);

  dlong out_stride = sizeof(struct eval_out_pt_2);
  dlong src_stride = sizeof(struct eval_src_pt_2);

  occa::memory d_in = *(occa::memory*)in;

  fd->local_eval_kernel(d_out_out_base,   out_stride,
                        d_out_proc_base,  out_stride,
                        d_out_index_base, out_stride,
                        d_src_el_base,    src_stride,
                        d_src_r_base,     src_stride,
                        d_src_proc_base,  src_stride,
                        d_src_index_base, src_stride,
                        pn, d_in, in_stride,
                        fd->lag_data[0], fd->lag_data[1]);
  d_out_pt.copyTo(opt, sizeof(struct eval_out_pt_2)*pn);
}



void ogs_findpts_local_eval_internal_3(
  struct eval_out_pt_3 *opt, const struct eval_src_pt_3 *spt,
  const dlong pn, const void *const in, const dlong in_stride,
  dlong *const n, dfloat *const lag_data[3], dlong lag_data_size[3],
  const void *const fd_void)
{
  if (pn == 0) return;

  ogs_findpts_t *fd = (ogs_findpts_t*)fd_void;

  const dlong nr=n[0],ns=n[1],nt=n[2];
  occa::device device = *fd->device;

  occa::memory workspace = device.malloc((sizeof(struct eval_out_pt_3)+sizeof(struct eval_src_pt_3))*pn,
                                         occa::dtype::byte);
  occa::memory d_out_pt = workspace; workspace += sizeof(struct eval_out_pt_3)*pn;
  occa::memory d_src_pt = workspace;
  d_src_pt.copyFrom(spt, sizeof(struct eval_src_pt_3)*pn);

  occa::memory d_out_out_base   = d_out_pt + offsetof(struct eval_out_pt_3, out);
  occa::memory d_out_proc_base  = d_out_pt + offsetof(struct eval_out_pt_3, proc);
  occa::memory d_out_index_base = d_out_pt + offsetof(struct eval_out_pt_3, index);
  occa::memory d_src_el_base    = d_src_pt + offsetof(struct eval_src_pt_3, el);
  occa::memory d_src_r_base     = d_src_pt + offsetof(struct eval_src_pt_3, r);
  occa::memory d_src_proc_base  = d_src_pt + offsetof(struct eval_src_pt_3, proc);
  occa::memory d_src_index_base = d_src_pt + offsetof(struct eval_src_pt_3, index);

  dlong out_stride = sizeof(struct eval_out_pt_3);
  dlong src_stride = sizeof(struct eval_src_pt_3);

  occa::memory d_in = *(occa::memory*)in;

  fd->local_eval_kernel(d_out_out_base,   out_stride,
                        d_out_proc_base,  out_stride,
                        d_out_index_base, out_stride,
                        d_src_el_base,    src_stride,
                        d_src_r_base,     src_stride,
                        d_src_proc_base,  src_stride,
                        d_src_index_base, src_stride,
                        pn, d_in, in_stride,
                        fd->lag_data[0], fd->lag_data[1], fd->lag_data[2]);

  d_out_pt.copyTo(opt, sizeof(struct eval_out_pt_3)*pn);
}

}
