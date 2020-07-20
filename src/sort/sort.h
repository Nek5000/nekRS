#ifndef _SORT_H_
#define _SORT_H_

#include <genmap-impl.h>

typedef enum{
  bin_sort      =0,
  hypercube_sort=1
} sort_algo;

typedef struct{
  int nfields;
  gs_dom t[3];
  uint offset[3];

  struct array *a;
  size_t unit_size,align;

  int balance;
  sort_algo algo;

  buffer buf;
} sort_data_private;
typedef sort_data_private* sort_data;
//
// parallel_bin_sort
//
int parallel_sort_private(sort_data data,struct comm *c);
int parallel_bin_sort(sort_data data,struct comm *c);
//
// parallel_hypercube_sort
//
typedef struct{
  sort_data data;
  int nprobes;
  double *probes;
  ulong *probe_cnt;
} hypercube_sort_data_private;
typedef hypercube_sort_data_private* hypercube_sort_data;

int parallel_hypercube_sort(hypercube_sort_data data,struct comm *c);
//
// Uniform parallel sort
//
#define parallel_sort(T,A,off,type,c) do {\
  sort_data_private sd;\
  sd.unit_size=sizeof(T);\
  sd.align=ALIGNOF(T);\
  sd.nfields=1;\
  sd.t[0]=type;\
  sd.offset[0]=off;\
  sd.a=A;\
  sd.algo=bin_sort;\
  sd.balance=1;\
  buffer_init(&sd.buf,1024);\
  parallel_sort_private(&sd,c);\
  buffer_free(&sd.buf);\
} while (0)

#define parallel_sort_2(T,A,off1,t1,off2,t2,c) do {\
  sort_data_private sd;\
  sd.unit_size=sizeof(T);\
  sd.align=ALIGNOF(T);\
  sd.nfields=2;\
  sd.t[0]=t1;\
  sd.offset[0]=off1;\
  sd.t[1]=t2;\
  sd.offset[1]=off2;\
  sd.a=A;\
  sd.algo=bin_sort;\
  sd.balance=1;\
  buffer_init(&sd.buf,1024);\
  parallel_sort_private(&sd,c);\
  buffer_free(&sd.buf);\
} while (0)

#endif
