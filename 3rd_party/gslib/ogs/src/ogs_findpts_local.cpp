
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

void ogs_findpts_local_eval_internal_2(
        double   *const out_base, const unsigned out_stride,
  const unsigned *const  el_base, const unsigned  el_stride,
  const double   *const   r_base, const unsigned   r_stride,
  const unsigned pn, const void *const in, const unsigned in_stride,
  unsigned *const n, double *const lag_data[2], unsigned lag_data_size[2])
{
  if (pn == 0) return;

  const unsigned nr=n[0],ns=n[1];
  occa::device device = platform->device;

  assert(nr <= MAX_GLL_N);
  assert(ns <= MAX_GLL_N);

  occa::json malloc_props;

  occa::memory workspace = device.malloc((out_stride+el_stride+r_stride)*pn,
                                         occa::dtype::byte,
                                         malloc_props);
  occa::memory d_out_base = workspace; workspace += out_stride*pn;
  occa::memory  d_el_base = workspace; workspace += el_stride*pn;
  occa::memory   d_r_base = workspace;
  d_el_base.copyFrom(el_base, el_stride*pn);
   d_r_base.copyFrom( r_base,  r_stride*pn);

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

  ogs::findpts_local_eval_2(d_out_base, out_stride,
                             d_el_base,  el_stride,
                              d_r_base,   r_stride,
                            pn, d_in, in_stride,
                            nr, ns, d_lag_data_0, d_lag_data_1);
  d_out_base.copyTo(out_base, out_stride*pn);
}



void ogs_findpts_local_eval_internal_3(
        double   *const out_base, const unsigned out_stride,
  const unsigned *const  el_base, const unsigned  el_stride,
  const double   *const   r_base, const unsigned   r_stride,
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

  occa::memory workspace = device.malloc((out_stride+el_stride+r_stride)*pn,
                                         occa::dtype::byte,
                                         malloc_props);
  occa::memory d_out_base = workspace; workspace += out_stride*pn;
  occa::memory  d_el_base = workspace; workspace += el_stride*pn;
  occa::memory   d_r_base = workspace;
  d_el_base.copyFrom(el_base, el_stride*pn);
   d_r_base.copyFrom( r_base,  r_stride*pn);

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

  ogs::findpts_local_eval_3(d_out_base, out_stride,
                             d_el_base,  el_stride,
                              d_r_base,   r_stride,
                            pn, d_in, in_stride,
                            nr, ns, nt, d_lag_data_0, d_lag_data_1, d_lag_data_2);

  d_out_base.copyTo(out_base, out_stride*pn);
}

}
