#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include <genmap-impl.h>
#include <parRSB.h>

#define MAXNV 8 /* maximum number of vertices per element */
typedef struct {
  int proc;
  GenmapLong id;
  int part;
  long long vtx[MAXNV];
} elm_data;

void fparRSB_partMesh(int *part,long long *vtx,int *nel,int *nve,
  int *options,int *comm,int *err)
{
  *err = 1;
  GenmapCommExternal c;
  c = MPI_Comm_f2c(*comm);
  *err = parRSB_partMesh(part, vtx, *nel, *nve, options, c);
}

int parRSB_partMesh(int *part,long long *vtx,int nel,int nve,
  int *options,MPI_Comm comm)
{
  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  /* load balance input data */
  GenmapLong nelg;
  GenmapLong nell = nel;
  MPI_Allreduce(&nell,&nelg,1,MPI_LONG_LONG_INT,MPI_SUM,comm);
  GenmapLong nstar = nelg/size;
  if(nstar == 0) nstar = 1;

  GenmapLong nelg_start;
  MPI_Scan(&nell,&nelg_start,1,MPI_LONG_LONG_INT,MPI_SUM,comm);
  nelg_start-=nel;

  struct array eList;
  elm_data *data;

  array_init(elm_data,&eList,nel), eList.n=nel;
  int e, n;
  for(data=eList.ptr,e=0; e<nel; ++e) {
    data[e].id=nelg_start+(e+1);
    const GenmapLong eg=data[e].id;

    data[e].proc=(int) ((eg-1)/nstar);
    if(eg>size*nstar) data[e].proc= (eg%size)-1; 

    for(n=0; n<nve; ++n) {
      data[e].vtx[n]=vtx[e*nve+n];
    }
  }

  struct comm c; comm_init(&c,comm);
  struct crystal cr; crystal_init(&cr,&c);
  buffer buf; buffer_init(&buf, 1024);

  sarray_transfer(elm_data,&eList,proc,1,&cr);
  data=eList.ptr;
  nel =eList.n;

  sarray_sort(elm_data,data,(unsigned int)nel,id,TYPE_LONG,&buf);

  MPI_Comm commRSB;
  MPI_Comm_split(comm, nel>0, rank, &commRSB);

  if(nel>0) {
    double time0 = comm_time();
    GenmapHandle h;
    GenmapInit(&h, commRSB);

    if(options[0] != 0) {
      h->dbgLevel = options[1];
      h->printStat = options[2];
    }

    GenmapSetNLocalElements(h, (GenmapInt)nel);
    GenmapScan(h, GenmapGetGlobalComm(h));
    GenmapSetNVertices(h, nve);

    GenmapLong nelg = GenmapGetNGlobalElements(h);
    GenmapInt id = GenmapCommRank(GenmapGetGlobalComm(h));
    GenmapInt size_ = GenmapCommSize(GenmapGetGlobalComm(h));
    if((GenmapLong)size_ > nelg) {
      if(id == 0)
        printf("Total number of elements is smaller than the "
          "number of processors.\n"
          "Run with smaller number of processors.\n");
      return 1;
    }
    GenmapElements e = GenmapGetElements(h);
    GenmapLong start = GenmapGetLocalStartIndex(h);

    GenmapInt i, j;
    for(i = 0; i < nel; i++) {
      e[i].origin = id;
      for(j = 0; j < nve; j++) {
        e[i].vertices[j] = data[i].vtx[j];
      }
    }

    GenmapRSB(h,h->dbgLevel>1);

    e = GenmapGetElements(h);
    for(j = 0; j < GenmapGetNLocalElements(h); j++) {
      e[j].proc = GenmapCommRank(GenmapGetGlobalComm(h));
    }

    GenmapCrystalInit(h, GenmapGetGlobalComm(h));
    GenmapCrystalTransfer(h, GENMAP_ORIGIN);
    GenmapCrystalFinalize(h);

    assert(GenmapGetNLocalElements(h) == nel);

    e = GenmapGetElements(h);
    sarray_sort(struct GenmapElement_private,e,(unsigned int)nel,globalId,
                TYPE_LONG, &buf);

    for(i = 0; i < nel; i++) {
      data[i].part = e[i].proc;
    }

    if(id == 0 && h->dbgLevel > 0)
      printf(" finished in %lfs\n", comm_time() - time0);

    GenmapFinalize(h);
    fflush(stdout);
  }

  MPI_Comm_free(&commRSB);

  /* restore original input */
  sarray_transfer(elm_data,&eList,proc,0,&cr);
  data=eList.ptr;
  nel =eList.n;
  sarray_sort(elm_data,data,(unsigned int)nel,id,TYPE_LONG,&buf);
  MPI_Barrier(comm);

  for(e = 0; e < nel; e++) {
    part[e]=data[e].part;
  }

  array_free(&eList);
  buffer_free(&buf);
  crystal_free(&cr);
  comm_free(&c);

  return 0;
}
