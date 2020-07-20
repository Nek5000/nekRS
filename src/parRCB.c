#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include <sort.h>
#include <parRSB.h>

void fparRCB_partMesh(int *part,double *vtx,int *nel,int *nv,
  int *options,int *comm,int *err)
{
  *err = 1;

  comm_ext c; c = MPI_Comm_f2c(*comm);
  *err=parRCB_partMesh(part,vtx,*nel,*nv,options,c);
}

int parRCB_partMesh(int *part,double *vtx,int nel,int nv,
  int *options,MPI_Comm comm)
{
  struct comm c; comm_init(&c,comm);
  int rank=c.id,size=c.np;

  /* load balance input data */
  slong out[2][1],buf[2][1];
  slong nell=nel;
  comm_scan(out,&c,gs_long,gs_add,&nell,1,&buf);
  slong nelg_start=out[0][0];
  slong nelg      =out[1][0];

  struct array a; array_init(elm_rcb,&a,nel);
  elm_rcb *data=a.ptr;

  int ndim=(nv==8)?3:2;

  int e, n;
  for(e=0;e<nel;++e){
    data[e].id=nelg_start+(e+1);
    data[e].orig=rank;
    for(int n=0;n<ndim;n++)
      data[e].coord[n]=vtx[e*ndim+n];
  }
  a.n=nel;

  //TODO: load balance

  struct comm rcb;
  comm_ext old=c.c;
#ifdef MPI
  MPI_Comm new; MPI_Comm_split(old,nel>0,rank,&new);
  comm_init(&rcb,new); MPI_Comm_free(&new);
#else
  comm_init(&rcb,1);
#endif

  if(nel>0){
    comm_barrier(&rcb);
    double time=comm_time();

    parRCB(&rcb,&a,ndim);

    comm_barrier(&rcb);
    time=comm_time()-time;

    if(c.id==0)
      printf(" finished in %lfs\n",time);
    fflush(stdout);
  }
  comm_free(&rcb);

  /* restore original input */
  struct crystal cr; crystal_init(&cr,&c);
  sarray_transfer(elm_rcb,&a,orig,1,&cr);
  crystal_free(&cr);

  comm_free(&c);

  assert(a.n==nel);

  buffer b; buffer_init(&b,1024);
  sarray_sort(elm_rcb,a.ptr,a.n,id,1,&b);
  buffer_free(&b);

  data=a.ptr;
  for(int e=0;e<nel;e++) part[e]=data[e].orig;

  return 0;
}
