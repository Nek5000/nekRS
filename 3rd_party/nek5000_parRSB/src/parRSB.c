#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include <genmap-impl.h>
#include <parRSB.h>

void fparRSB_partMesh(int *part,long long *vtx,double *coord,int *nel,
  int *nve,int *options,int *comm,int *err)
{
  *err = 1;
  comm_ext c = MPI_Comm_f2c(*comm);
  *err = parRSB_partMesh(part, vtx, coord, *nel, *nve, options, c);
}

/*
 * part = [nel], out,
 * vtx = [nel x nve], in,
 * coord = [nel x nve x ndim], in,
 * nel = in,
 * nve = in,
 * options = null or [options[0]], in */
int parRSB_partMesh(int *part,long long *vtx,double *coord,int nel,
  int nve,int *options,MPI_Comm comm)
{
  struct comm c; comm_init(&c,comm);
  int rank=c.id,size=c.np;

  if(rank==0)
    printf("running RSB ...");
  fflush(stdout);

  comm_barrier(&c);
  double time0=comm_time();

  /* Load balance input data */
  slong out[2][1],buf[3][1],in=nel;
  comm_scan(out,&c,gs_long,gs_add,&in,1,buf);
  slong nelg_start=out[0][0];
  slong nelg=out[1][0];

  GenmapLong nstar=nelg/size;
  if(nstar==0) nstar=1;

  struct array eList;
  array_init(struct rsb_element,&eList,nel);

  struct rsb_element data;
  data.type=GENMAP_RSB_ELEMENT;
  data.origin=rank;

  int ndim=(nve==8)?3:2;
  int e,n,v;
  for(e=0; e<nel; ++e){
    data.globalId=nelg_start+(e+1);

    GenmapLong eg=data.globalId;
    data.proc=(int) ((eg-1)/nstar);
    if(eg>size*nstar) data.proc= (eg%size)-1;

    data.coord[0]=data.coord[1]=data.coord[2]=0.0;
    for(n=0; n<nve; ++n) {
      data.vertices[n]=vtx[e*nve+n];
      for(v=0;v<ndim;v++)
        data.coord[v]+=coord[e*ndim*nve+n*ndim+v];
    }
    for(v=0;v<ndim;v++)
      data.coord[v]/=nve;

    array_cat(struct rsb_element,&eList,&data,1);
  }
  assert(eList.n==nel);

  struct crystal cr; crystal_init(&cr,&c);
  sarray_transfer(struct rsb_element,&eList,proc,1,&cr);
  nel=eList.n;
  struct rsb_element *e_ptr=eList.ptr;

  buffer bfr; buffer_init(&bfr,1024);
  sarray_sort(struct rsb_element,eList.ptr,(unsigned)nel,globalId,1,&bfr);

  double time1=comm_time();
  comm_barrier(&c);
  double time2=comm_time();

  /* Run RSB now */
  struct comm comm_rsb;
  comm_ext old=c.c;
#ifdef MPI
  MPI_Comm new; MPI_Comm_split(old,nel>0,rank,&new);
  comm_init(&comm_rsb,new); MPI_Comm_free(&new);
#else
  comm_init(&comm_rsb,1);
#endif

  int verbose_level = 0;
  int print_stat = 0;
  if(options!=NULL){
    if(options[0]>0)
      verbose_level=options[1];
    if(options[0]>1)
      print_stat=options[2];
  }

  if(nel>0) {
    metric_init();

    genmap_handle h;
    genmap_init(&h, comm_rsb.c, options);

    GenmapSetArrayElements(h,&eList);
    GenmapScan(h, GenmapGetGlobalComm(h));
    GenmapSetNVertices(h, nve);

    GenmapLong nelg=GenmapGetNGlobalElements(h);
    GenmapInt id=GenmapCommRank(GenmapGetGlobalComm(h));
    GenmapInt size_=GenmapCommSize(GenmapGetGlobalComm(h));
    if((GenmapLong)size_>nelg){
      if(id==0)
        printf("Total number of elements is smaller than the "
          "number of processors.\n"
          "Run with smaller number of processors.\n");
      return 1;
    }

    genmap_rsb(h,h->verbose_level>1);

    genmap_finalize(h);

    if(print_stat>0)
      metric_print(&c);
    metric_finalize();
  }

  double time3=comm_time();
  comm_barrier(&c);
  double time4=comm_time();

  /* Restore original input */
  sarray_transfer(struct rsb_element,&eList,origin,1,&cr);
  nel=eList.n;
  sarray_sort(struct rsb_element,eList.ptr,nel,globalId,1,&bfr);

  e_ptr=eList.ptr;
  for(e=0;e<nel;e++)
    part[e]=e_ptr[e].origin;

  double time5=comm_time();
  comm_barrier(&c);
  double time=comm_time()-time0;

  if(rank==0)
    printf(" finished in %g s\n",time);

  /* Report time and finish */
  if(print_stat>0){
    double min[3],max[3],sum[3];
    min[0]=max[0]=sum[0]=time1-time0;
    min[1]=max[1]=sum[1]=time3-time2;
    min[2]=max[2]=sum[2]=time5-time4;
    comm_allreduce(&c,gs_double,gs_min,min,3,buf); // min
    comm_allreduce(&c,gs_double,gs_max,max,3,buf); // max
    comm_allreduce(&c,gs_double,gs_add,sum,3,buf); // sum
    if(rank==0){
      printf("LOADBALANCE : %g/%g/%g\n",min[0],max[0],sum[0]/c.np);
      printf("RSB         : %g/%g/%g\n",min[1],max[1],sum[1]/c.np);
      printf("RESTORE     : %g/%g/%g\n",min[2],max[2],sum[2]/c.np);
    }
    fflush(stdout);
  }

  array_free(&eList);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&comm_rsb);
  comm_free(&c);

  return 0;
}
