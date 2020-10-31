#include <sort.h>
#include <float.h>
#include <parRSB.h>
#include <genmap-impl.h>

void get_axis_len(double *length,struct array *a,struct comm *c,int ndim)
{
  double min[MAXDIM],max[MAXDIM];
  sint i;
  for(i=0;i<ndim;i++) min[i]=DBL_MAX,max[i]=-DBL_MAX;

  sint nel=a->n;
  elm_rcb* elems=a->ptr;
  for(i=0;i<nel;i++){
    if(elems[i].coord[0]<min[0]) min[0]=elems[i].coord[0];
    if(elems[i].coord[1]<min[1]) min[1]=elems[i].coord[1];
    if(ndim==3)
      if(elems[i].coord[2]<min[2]) min[2]=elems[i].coord[2];

    if(elems[i].coord[0]>max[0]) max[0]=elems[i].coord[0];
    if(elems[i].coord[1]>max[1]) max[1]=elems[i].coord[1];
    if(ndim==3)
      if(elems[i].coord[2]>max[2]) max[2]=elems[i].coord[2];
  }

  double buf[MAXDIM];
  comm_allreduce(c,gs_double,gs_min,min,MAXDIM,buf);
  comm_allreduce(c,gs_double,gs_max,max,MAXDIM,buf);

  for(i=0;i<ndim;i++)
    length[i]=max[i]-min[i];
}

int parRCB(struct comm *ci,struct array *a,int ndim){
  struct comm c; comm_dup(&c,ci);
  sint rank=c.id;
  sint size=c.np;

  double length[MAXDIM];

  while(size>1){
    metric_tic(&c,AXISLEN);
    get_axis_len(length,a,&c,ndim);
    metric_toc(&c,AXISLEN);

    int axis1=0,d;
    for(d=1;d<ndim;d++)
      if(length[d]>length[axis1]) axis1=d;

    metric_tic(&c,PARSORT);
    switch(axis1){
      case 0:
        parallel_sort(elm_rcb,a,coord[0],gs_double,0,1,&c);
        break;
      case 1:
        parallel_sort(elm_rcb,a,coord[1],gs_double,0,1,&c);
        break;
      case 2:
        parallel_sort(elm_rcb,a,coord[2],gs_double,0,1,&c);
        break;
      default:
        break;
    }
    metric_toc(&c,PARSORT);

    int p=(size+1)/2;
    int bin=(rank>=p);

    metric_tic(&c,COMMSPLIT);
    comm_ext old=c.c;
#ifdef MPI
    MPI_Comm new; MPI_Comm_split(old,bin,rank,&new);
    comm_free(&c); comm_init(&c,new);
    MPI_Comm_free(&new);
#endif
    metric_toc(&c,COMMSPLIT);

    metric_push_level();

    rank=c.id;
    size=c.np;
  }

  comm_free(&c);
  return 0;
}

void get_axis_len_local(double *length,elm_rcb *elems,uint nel,int ndim)
{
  double min[MAXDIM],max[MAXDIM];
  sint i;
  for(i=0;i<ndim;i++) min[i]=DBL_MAX,max[i]=-DBL_MAX;

  for(i=0;i<nel;i++){
    if(elems[i].coord[0]<min[0]) min[0]=elems[i].coord[0];
    if(elems[i].coord[1]<min[1]) min[1]=elems[i].coord[1];
    if(ndim==3)
      if(elems[i].coord[2]<min[2]) min[2]=elems[i].coord[2];

    if(elems[i].coord[0]>max[0]) max[0]=elems[i].coord[0];
    if(elems[i].coord[1]>max[1]) max[1]=elems[i].coord[1];
    if(ndim==3)
      if(elems[i].coord[2]>max[2]) max[2]=elems[i].coord[2];
  }

  for(i=0;i<ndim;i++)
    length[i]=max[i]-min[i];
}

void rcb_local(struct array *a,uint start,uint end,int ndim,buffer *buf){
  assert(end>=start);

  uint size=end-start;
  if(size<=2) return;

  double length[3]; elm_rcb *st=a->ptr; st+=start;
  get_axis_len_local(length,st,size,ndim);

  int axis=0;
  if(fabs(length[axis])<fabs(length[1])) axis=1;
  if(ndim==3)
    if(fabs(length[axis])<fabs(length[2])) axis=2;

  switch(axis){
    case 0:
      sarray_sort(elm_rcb,st,size,coord[0],3,buf);
      break;
    case 1:
      sarray_sort(elm_rcb,st,size,coord[1],3,buf);
      break;
    case 2:
      sarray_sort(elm_rcb,st,size,coord[2],3,buf);
      break;
    default:
      break;
  }

  uint mid=(start+end)/2;
  rcb_local(a,start,mid,ndim,buf);
  rcb_local(a,mid  ,end,ndim,buf);
}
