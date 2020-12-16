#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include <genmap-impl.h>
#include <parRSB.h>

void fparRCB_partMesh(int *part,int *seq,double *coord,int *nel,int *nv,
  int *options,int *comm,int *err)
{
  *err = 1;
  comm_ext c = MPI_Comm_f2c(*comm);
  *err=parRCB_partMesh(part,seq,coord,*nel,*nv,options,c);
}

// coord = [nel,nv,ndim]
int parRCB_partMesh(int *part,int *seq,double *coord,int nel,int nv,
  int *options,MPI_Comm comm)
{
  struct comm c; comm_init(&c,comm);
  int rank=c.id,size=c.np;

  if(rank==0)
    printf("running RCB ... ");
  fflush(stdout);

  comm_barrier(&c);
  double time0=comm_time();

  /* Load balance input data */
  slong out[2][1],buf[2][1],in=nel;
  comm_scan(out,&c,gs_long,gs_add,&in,1,buf);
  slong nelg_start=out[0][0];
  slong nelg=out[1][0];

  GenmapLong nstar=nelg/size;
  if(nstar==0) nstar=1;

  struct array eList;
  array_init(struct rcb_element,&eList,nel);

  struct rcb_element data;
  data.type=GENMAP_RCB_ELEMENT;
  data.origin=rank;

  int ndim=(nv==8)?3:2;
  int e,n,v;
  for(e=0;e<nel;++e){
    data.globalId=nelg_start+(e+1);

    GenmapLong eg=data.globalId;
    data.proc=(int) ((eg-1)/nstar);
    if(eg>size*nstar) data.proc= (eg%size)-1;

    data.coord[0]=data.coord[1]=data.coord[2]=0.0;
    for(v=0;v<nv;v++){

      for(n=0;n<ndim;n++)
        data.coord[n]+=coord[e*ndim*nv+v*ndim+n];
    }
    for(n=0;n<ndim;n++)
      data.coord[n]/=nv;

    array_cat(struct rcb_element,&eList,&data,1);
  }
  assert(eList.n==nel);

  struct crystal cr; crystal_init(&cr,&c);
  sarray_transfer(struct rcb_element,&eList,proc,1,&cr);
  nel=eList.n;
  struct rcb_element *e_ptr=eList.ptr;

  buffer bfr; buffer_init(&bfr,1024);
  sarray_sort(struct rcb_element,eList.ptr,(unsigned)nel,globalId,1,&bfr);

  double time1=comm_time();
  comm_barrier(&c);
  double time2=comm_time();

  /* Run RSB now */
  struct comm comm_rcb;
  comm_ext old=c.c;
#ifdef MPI
  MPI_Comm new; MPI_Comm_split(old,nel>0,rank,&new);
  comm_init(&comm_rcb,new); MPI_Comm_free(&new);
#else
  comm_init(&comm_rcb,1);
#endif

  if(nel>0){
    metric_init();

    rcb(&comm_rcb,&eList,ndim);

    // Do a local RCB if seq!=NULL
    if(seq!=NULL){
      rcb_local(&eList,0,eList.n,ndim,&bfr);

      e_ptr=eList.ptr;
      for(e=0; e<eList.n; e++)
        e_ptr[e].seq=e;
    }

    metric_finalize();
  }

  /* Restore original input */
  sarray_transfer(struct rcb_element,&eList,origin,1,&cr);
  nel=eList.n;
  sarray_sort(struct rcb_element,eList.ptr,eList.n,globalId,1,&bfr);

  e_ptr=eList.ptr;
  for(e=0;e<nel;e++)
    part[e]=e_ptr[e].origin;

  if(seq!=NULL)
    for(e=0;e<nel;e++)
      seq[e]=e_ptr[e].seq;

  comm_barrier(&c);
  double time=comm_time()-time0;

  /* Report time and finish */
  if(c.id==0)
    printf(" finished in %g s\n",time);
  fflush(stdout);

  array_free(&eList);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&comm_rcb);
  comm_free(&c);

  return 0;
}
