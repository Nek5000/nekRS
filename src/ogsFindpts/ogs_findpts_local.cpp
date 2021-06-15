
#include <type_traits>
#include "ogstypes.h"
#include "ogs_FINDPTS.hpp"

extern "C" {

#include "gslib.h"
#include "ogs_findpts.h"

#define   AT(T,var,i)   \
        (T*)(      (char*)var##_base   +(i)*var##_stride   )
#define  CAT(T,var,i) \
  (const T*)((const char*)var##_base   +(i)*var##_stride   )
#define CATD(T,var,i,d) \
  (const T*)((const char*)var##_base[d]+(i)*var##_stride[d])

static_assert(std::is_same<dfloat, double>::value, "OGS dfloat is not compatible with GSLIB double");
static_assert(sizeof(dlong) == sizeof(uint), "OGS dlong is not compatible with GSLIB uint");

void ogs_findpts_local_eval_internal_2(
    struct eval_out_pt_2 *opt, const struct eval_src_pt_2 *spt,
    const uint pn, const void *const in,
    struct findpts_local_data_2 *const gs_fd, const void *const ogs_fd_void)
{
  if (pn == 0) return;

  ogs_findpts_t *ogs_fd = (ogs_findpts_t*)ogs_fd_void;
  occa::device device = *ogs_fd->device;
  occa::memory d_in = *(occa::memory*)in;

  occa::memory workspace = device.malloc((sizeof(struct eval_out_pt_2)+sizeof(struct eval_src_pt_2))*pn,
                                         occa::dtype::byte);
  occa::memory d_out_pt = workspace; workspace += sizeof(struct eval_out_pt_2)*pn;
  occa::memory d_src_pt = workspace;
  d_src_pt.copyFrom(spt, sizeof(struct eval_src_pt_2)*pn);

  occa::memory d_out_base   = d_out_pt + offsetof(struct eval_out_pt_2, out);
  occa::memory d_el_base    = d_src_pt + offsetof(struct eval_src_pt_2, el);
  occa::memory d_r_base     = d_src_pt + offsetof(struct eval_src_pt_2, r);

  dlong out_stride = sizeof(struct eval_out_pt_2);
  dlong src_stride = sizeof(struct eval_src_pt_2);

  ogs_fd->local_eval_kernel(d_out_base,   out_stride,
                            d_el_base,    src_stride,
                            d_r_base,     src_stride,
                            pn, d_in, gs_fd->ntot,
                            ogs_fd->lag_data[0], ogs_fd->lag_data[1]);

  d_out_pt.copyTo(opt, sizeof(struct eval_out_pt_2)*pn);
}



void ogs_findpts_local_eval_internal_3(
    struct eval_out_pt_3 *opt, const struct eval_src_pt_3 *spt,
    const uint pn, const void *const in,
    struct findpts_local_data_3 *const gs_fd, const void *const ogs_fd_void)
{
  if (pn == 0) return;

  ogs_findpts_t *ogs_fd = (ogs_findpts_t*)ogs_fd_void;
  occa::device device = *ogs_fd->device;
  occa::memory d_in = *(occa::memory*)in;

  occa::memory workspace = device.malloc((sizeof(struct eval_out_pt_3)+sizeof(struct eval_src_pt_3))*pn,
                                         occa::dtype::byte);
  occa::memory d_out_pt = workspace; workspace += sizeof(struct eval_out_pt_3)*pn;
  occa::memory d_src_pt = workspace;
  d_src_pt.copyFrom(spt, sizeof(struct eval_src_pt_3)*pn);

  occa::memory d_out_base   = d_out_pt + offsetof(struct eval_out_pt_3, out);
  occa::memory d_el_base    = d_src_pt + offsetof(struct eval_src_pt_3, el);
  occa::memory d_r_base     = d_src_pt + offsetof(struct eval_src_pt_3, r);

  dlong out_stride = sizeof(struct eval_out_pt_3);
  dlong src_stride = sizeof(struct eval_src_pt_3);

  ogs_fd->local_eval_kernel(d_out_base,   out_stride,
                            d_el_base,    src_stride,
                            d_r_base,     src_stride,
                            pn, d_in, gs_fd->ntot,
                            ogs_fd->lag_data[0], ogs_fd->lag_data[1], ogs_fd->lag_data[2]);

  d_out_pt.copyTo(opt, sizeof(struct eval_out_pt_3)*pn);
}

void ogs_findpts_local_eval_2(
          dfloat * const out_base, const uint out_stride,
    const uint   * const el_base,  const uint el_stride,
    const dfloat * const r_base,   const uint r_stride,
    const uint pn, const void *const in,
    struct findpts_local_data_2 *const gs_fd, const void *const ogs_fd_void)
{
  if (pn == 0) return;

  ogs_findpts_t *ogs_fd = (ogs_findpts_t*)ogs_fd_void;
  occa::device device = *ogs_fd->device;
  occa::memory d_in = *(occa::memory*)in;

  occa::memory workspace = device.malloc((out_stride+el_stride+r_stride)*pn,
                                         occa::dtype::byte);
  occa::memory d_out_base = workspace; workspace += out_stride*pn;
  occa::memory d_el_base  = workspace; workspace += el_stride*pn;
  occa::memory d_r_base   = workspace; workspace += r_stride*pn;
  d_el_base.copyFrom(el_base, el_stride*pn);
  d_r_base .copyFrom(r_base,  r_stride*pn);
  // don't need to copy out to device if it will be entirely overwritten
  if (out_stride != sizeof(dfloat)) {
    d_out_base.copyFrom(out_base, out_stride*pn);
  }

  ogs_fd->local_eval_kernel(d_out_base,   out_stride,
                            d_el_base,    el_stride,
                            d_r_base,     r_stride,
                            pn, d_in, gs_fd->ntot,
                            ogs_fd->lag_data[0], ogs_fd->lag_data[1]);

  d_out_base.copyTo(out_base, out_stride*pn);
}

void ogs_findpts_local_eval_3(
          dfloat * const out_base, const uint out_stride,
    const uint   * const el_base,  const uint el_stride,
    const dfloat * const r_base,   const uint r_stride,
    const uint pn, const void *const in,
    struct findpts_local_data_3 *const gs_fd, const void *const ogs_fd_void)
{
  if (pn == 0) return;

  ogs_findpts_t *ogs_fd = (ogs_findpts_t*)ogs_fd_void;
  occa::device device = *ogs_fd->device;
  occa::memory d_in = *(occa::memory*)in;

  occa::memory workspace = device.malloc((out_stride+el_stride+r_stride)*pn,
                                         occa::dtype::byte);
  occa::memory d_out_base = workspace; workspace += out_stride*pn;
  occa::memory d_el_base  = workspace; workspace += el_stride*pn;
  occa::memory d_r_base   = workspace; workspace += r_stride*pn;
  d_el_base.copyFrom(el_base, el_stride*pn);
  d_r_base .copyFrom(r_base,  r_stride*pn);
  // don't need to copy out to device if it will be entirely overwritten
  if (out_stride != sizeof(dfloat)) {
    d_out_base.copyFrom(out_base, out_stride*pn);
  }

  ogs_fd->local_eval_kernel(d_out_base,   out_stride,
                            d_el_base,    el_stride,
                            d_r_base,     r_stride,
                            pn, d_in, gs_fd->ntot,
                            ogs_fd->lag_data[0], ogs_fd->lag_data[1], ogs_fd->lag_data[2]);

  d_out_base.copyTo(out_base, out_stride*pn);
}

}
