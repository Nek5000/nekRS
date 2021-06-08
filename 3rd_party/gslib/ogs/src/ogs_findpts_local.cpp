
// work around the lack of restrict in C++
#define restrict

#include <cassert>
#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "platform.hpp"

#define   AT(T,var,i)   \
        (T*)(      (char*)var##_base   +(i)*var##_stride   )
#define  CAT(T,var,i) \
  (const T*)((const char*)var##_base   +(i)*var##_stride   )
#define CATD(T,var,i,d) \
  (const T*)((const char*)var##_base[d]+(i)*var##_stride[d])

extern "C" {

  struct eval_src_pt_3 { double r[3]; uint index, proc, el; };
  struct eval_src_pt_2 { double r[2]; uint index, proc, el; };
  struct eval_out_pt_3 { double out; uint index, proc; };
  struct eval_out_pt_2 { double out; uint index, proc; };

void ogs_findpts_local_eval_internal_2(
  struct eval_out_pt_2 *opt, const struct eval_src_pt_2 *spt,
  const unsigned pn, const void *const in, const unsigned in_stride,
  unsigned *const n, double *const lag_data[2], unsigned lag_data_size[2])
{
  if (pn == 0) return;

  const unsigned nr=n[0],ns=n[1];
  occa::device device = platform->device;

  assert(nr <= MAX_GLL_N);
  assert(ns <= MAX_GLL_N);

  occa::json malloc_props;

  occa::memory workspace = device.malloc((sizeof(struct eval_out_pt_2)+sizeof(struct eval_src_pt_2))*pn,
                                         occa::dtype::byte,
                                         malloc_props);
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

  unsigned int out_stride = sizeof(struct eval_out_pt_2);
  unsigned int src_stride = sizeof(struct eval_src_pt_2);


  occa::memory d_lag_data_0 = device.malloc(lag_data_size[0],
                                            occa::dtype::double_,
                                            malloc_props);
  d_lag_data_0.copyFrom(lag_data[0]);
  // reuse lag_data_0 if all directions have the same degree
  occa::memory d_lag_data_1;
  if (nr == ns) {
    d_lag_data_1 = d_lag_data_0;
  } else {
    d_lag_data_1 = device.malloc(lag_data_size[1],
                                            occa::dtype::double_,
                                            malloc_props);
    d_lag_data_1.copyFrom(lag_data[1]);
  }


  occa::memory d_in = *(occa::memory*)in;

  ogs::findpts_local_eval_2(d_out_out_base,   out_stride,
                            d_out_proc_base,  out_stride,
                            d_out_index_base, out_stride,
                            d_src_el_base,    src_stride,
                            d_src_r_base,     src_stride,
                            d_src_proc_base,  src_stride,
                            d_src_index_base, src_stride,
                            pn, d_in, in_stride,
                            nr, ns, d_lag_data_0, d_lag_data_1);
  d_out_pt.copyTo(opt, sizeof(struct eval_out_pt_2)*pn);
}



void ogs_findpts_local_eval_internal_3(
  struct eval_out_pt_3 *opt, const struct eval_src_pt_3 *spt,
  const unsigned pn, const void *const in, const unsigned in_stride,
  unsigned *const n, double *const lag_data[3], unsigned lag_data_size[3])
{
  if (pn == 0) return;

  const unsigned nr=n[0],ns=n[1],nt=n[2];
  occa::device device = platform->device;

  assert(nr <= MAX_GLL_N);
  assert(ns <= MAX_GLL_N);
  assert(nt <= MAX_GLL_N);

  occa::json malloc_props;

  occa::memory workspace = device.malloc((sizeof(struct eval_out_pt_3)+sizeof(struct eval_src_pt_3))*pn,
                                         occa::dtype::byte,
                                         malloc_props);
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

  unsigned int out_stride = sizeof(struct eval_out_pt_3);
  unsigned int src_stride = sizeof(struct eval_src_pt_3);

  occa::memory d_lag_data_0 = device.malloc(lag_data_size[0],
                                            occa::dtype::double_,
                                            malloc_props);
  d_lag_data_0.copyFrom(lag_data[0]);
  // reuse lag_data_0 if all directions have the same degree
  occa::memory d_lag_data_1, d_lag_data_2;
  if (nr == ns) {
    d_lag_data_1 = d_lag_data_0;
  } else {
    d_lag_data_1 = device.malloc(lag_data_size[1], occa::dtype::double_, malloc_props);
    d_lag_data_1.copyFrom(lag_data[1]);
  }
  if (nr == nt) {
    d_lag_data_2 = d_lag_data_0;
  } else {
    d_lag_data_2 = device.malloc(lag_data_size[2], occa::dtype::double_, malloc_props);
    d_lag_data_2.copyFrom(lag_data[2]);
  }

  occa::memory d_in = *(occa::memory*)in;

  ogs::findpts_local_eval_3(d_out_out_base,   out_stride,
                            d_out_proc_base,  out_stride,
                            d_out_index_base, out_stride,
                            d_src_el_base,    src_stride,
                            d_src_r_base,     src_stride,
                            d_src_proc_base,  src_stride,
                            d_src_index_base, src_stride,
                            pn, d_in, in_stride,
                            nr, ns, nt, d_lag_data_0, d_lag_data_1, d_lag_data_2);

  d_out_pt.copyTo(opt, sizeof(struct eval_out_pt_3)*pn);
}

}
