#include <sort-impl.h>
#include <math.h>

int init_probes(struct hypercube *data,struct comm *c)
{
  /* get input data */
  struct sort *input=data->data;

  /* Allocate space for probes and counts */
  int nprobes=data->nprobes=3;
  if(!data->probes   ) GenmapMalloc(nprobes,&data->probes   );
  if(!data->probe_cnt) GenmapMalloc(nprobes,&data->probe_cnt);

  double extrema[2];
  get_extrema((void*)extrema,data->data,0,c);
  double range=extrema[1]-extrema[0];
  double delta=range/(nprobes-1);

  data->probes[0]=extrema[0];
  data->probes[1]=extrema[0]+delta;
  data->probes[2]=extrema[1];

  return 0;
}

int update_probe_counts(struct hypercube *data,struct comm *c)
{
  struct sort *input=data->data;
  uint offset  =input->offset[0];
  gs_dom t=input->t[0];

  uint nprobes =data->nprobes;

  uint i;
  for(i=0;i<nprobes;i++) data->probe_cnt[i]=0;

  struct array *a=input->a;
  uint e;
  for(e=0;e<a->n;e++){
    double val_e=get_scalar(a,e,offset,input->unit_size,t);
    for(i=0;i<nprobes;i++)
      if(val_e<data->probes[i]) data->probe_cnt[i]++;
  }

  ulong buf[3];
  comm_allreduce(c,gs_long,gs_add,data->probe_cnt,nprobes,buf);

  return 0;
}

int update_probes(slong nelem,double *probes,ulong *probe_cnt,
    uint threshold)
{
  slong expected=nelem/2;
  if(llabs(expected-probe_cnt[1])<threshold) return 0;

  if(probe_cnt[1]<expected) probes[0]=probes[1];
  else probes[2]=probes[1];

  probes[1]=probes[0]+(probes[2]-probes[0])/2;

  return 0;
}

int transfer_elem(struct hypercube *data,struct comm *c)
{
  struct sort *input=data->data;
  struct array *a=input->a;
  uint usize     =input->unit_size;
  uint offset    =input->offset[0];
  gs_dom t  =input->t[0];

  uint size      =a->n;

  uint e,lown=0,uppern=0;
  for(e=0;e<size;e++){
    double val=get_scalar(a,e,offset,usize,t);
    if(val<data->probes[1]) lown++;
    else uppern++;
  }

  slong out[2][2],in[2],buf[2][2];
  in[0]=lown,in[1]=uppern;
  comm_scan(out,c,gs_long,gs_add,in,2,buf);

  ulong lstart=out[0][0],ustart=out[0][1];
  ulong lelem =out[1][0],uelem =out[1][1];

  uint np=c->np,lnp=np/2,unp=np-lnp;

  uint *proc; GenmapCalloc(size,&proc);
  set_dest(proc     ,lnp,lstart,lown  ,lelem);
  set_dest(proc+lown,unp,ustart,uppern,uelem);

  for(e=lown;e<size;e++) proc[e]+=lnp;

  struct crystal cr; crystal_init(&cr,c);
  sarray_transfer_ext_(a,usize,proc,sizeof(uint),&cr);
  crystal_free(&cr);

  GenmapFree(proc);

  return 0;
}

int parallel_hypercube_sort(struct hypercube *data,struct comm *c)
{
  struct sort *input=data->data;
  struct array *a   =input->a;
  gs_dom t          =input->t[0];
  uint offset       =input->offset[0];

  sint size=c->np,rank=c->id;

  slong out[2][1],buf[2][1],in=a->n;
  comm_scan(out,c,gs_long,gs_add,&in,1,buf);
  slong start=out[0][0];
  slong nelem=out[1][0];

  uint threshold=(nelem/(10*size));
  if(threshold<2) threshold=2;

  metric_tic(c,LOCALSORT);
  sort_local(data->data);
  metric_toc(c,LOCALSORT);

  if(size==1) return 0;

  metric_tic(c,UPDATEPROBE);
  init_probes        (data,c);
  update_probe_counts(data,c);
  metric_toc(c,UPDATEPROBE);

  int max_iter=log2((data->probes[2]-data->probes[0])/GENMAP_TOL),iter=0;
  while(llabs(nelem/2-data->probe_cnt[1])>threshold && iter++<max_iter){
    metric_tic(c,UPDATEPROBE);
    update_probes(nelem,data->probes,data->probe_cnt,threshold);
    update_probe_counts(data,c);
    metric_toc(c,UPDATEPROBE);
  }
  metric_tic(c,RCBTRANSFER);
  transfer_elem(data,c);
  metric_toc(c,RCBTRANSFER);

  // split the communicator
  struct comm nc;
  sint lower=(rank<size/2)?1:0;
#if defined(MPI)
  MPI_Comm nc_; MPI_Comm_split(c->c,lower,rank,&nc_);
  comm_init(&nc,nc_);
  MPI_Comm_free(&nc_);
#else
  comm_init(&nc,1);
#endif

  // TODO: Keep load balancing after each split
  parallel_hypercube_sort(data,&nc);
  comm_free(&nc);

  return 0;
}
