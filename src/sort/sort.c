#include <sort-impl.h>
#include <float.h>
#include <genmap-impl.h>  //FIXME - include genmap-statistics

int set_dest(uint *proc,uint np,ulong start,uint size,ulong nelem)
{
  uint psize=nelem/np,i;
  if(psize==0){
    for(i=0;i<size;i++) proc[i]=start+i;
    return 0;
  }

  uint nrem =nelem-np*psize;
  uint lower=nrem*(psize+1);

  uint id1,id2;
  if(start<lower) id1=start/(psize+1);
  else id1=nrem+(start-lower)/psize;

  if((start+size)<lower) id2=(start+size)/(psize+1);
  else id2=nrem+(start+size-lower)/psize;

  i=0;
  while(id1<=id2 && i<size){
    ulong s=id1*psize+min(id1,nrem);
    ulong e=(id1+1)*psize+min(id1+1,nrem);
    e=min(e,nelem);
    while(s<=start+i && start+i<e && i<size) proc[i++]=id1;
    id1++;
  }

  return 0;
}

int load_balance(struct array *a,size_t size,struct comm *c,
    struct crystal *cr)
{
  slong in=a->n,out[2][1],buf[2][1];
  comm_scan(out,c,gs_long,gs_add,&in,1,buf);
  ulong start=out[0][0];
  ulong nelem=out[1][0];

  uint *proc; GenmapCalloc(a->n,&proc);

  set_dest(proc,c->np,start,a->n,nelem);
  sarray_transfer_ext_(a,size,proc,sizeof(uint),cr);

  GenmapFree(proc);

  return 0;
}

double get_scalar(struct array *a,uint i,uint offset,uint usize,
    gs_dom type)
{
  char* v=(char*)a->ptr+i*usize+offset;
  double data;

  switch(type){
    case gs_int:
      data=*((uint*)v);
      break;
    case gs_long:
      data=*((ulong*)v);
      break;
    case gs_double:
      data=*((double*)v);
      break;
    default:
      break;
  }

  return data;
}

void get_extrema(void *extrema_,struct sort *data,
  uint field,struct comm* c)
{
  struct array *a=data->a;
  uint usize     =data->unit_size;
  uint offset    =data->offset[field];
  gs_dom t       =data->t[field];

  double *extrema=(double *)extrema_;

  sint size=a->n;
  if(size==0){
    extrema[0]=-DBL_MAX;
    extrema[1]=-DBL_MAX;
  } else {
    extrema[0]=get_scalar(a,0     ,offset,usize,t)*-1;
    extrema[1]=get_scalar(a,size-1,offset,usize,t);
  }

  double buf[2];
  comm_allreduce(c,gs_double,gs_max,extrema,2,buf);
  extrema[0]*=-1;
}

int parallel_sort_private(struct sort *data,struct comm *c){
  struct comm dup; comm_dup(&dup,c);

  int balance =data->balance;
  sort_algo algo=data->algo;

  struct array *a=data->a;
  size_t usize   =data->unit_size;

  struct hypercube hdata;

  switch(algo){
    case bin_sort:
      parallel_bin_sort(data,c);
      break;
    case hypercube_sort:
      hdata.data=data; hdata.probes=NULL; hdata.probe_cnt=NULL;
      parallel_hypercube_sort(&hdata,&dup);
      free(hdata.probes); free(hdata.probe_cnt);
      break;
    default:
      break;
  }

  metric_tic(c,LOADBALANCE0);
  if(balance){
    metric_tic(c,LOADBALANCE1);
    struct crystal cr; crystal_init(&cr,c);
    load_balance(a,usize,c,&cr);
    metric_toc(c,LOADBALANCE1);
    metric_tic(c,LOCALSORT);
    sort_local(data);
    metric_toc(c,LOCALSORT);
    crystal_free(&cr);
  }
  metric_toc(c,LOADBALANCE0);

  comm_free(&dup);

  return 0;
}
