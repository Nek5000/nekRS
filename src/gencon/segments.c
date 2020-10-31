#include <stdlib.h>
#include <math.h>

#include <gencon-impl.h>
#include <sort.h>

int findLocalSegments(Mesh mesh,int i,GenmapScalar tolSquared){
  Point pts=mesh->elements.ptr;
  sint npts=mesh->elements.n;

  sint j;
  for(j=0;j<npts-1;j++){
    GenmapScalar d =sqrDiff(pts[j].x[i],pts[j+1].x[i]);
    GenmapScalar dx=min(pts[j].dx,pts[j+1].dx)*tolSquared;
    if(d>dx) pts[j+1].ifSegment=1;
  }
}

int mergeSegments(Mesh mesh,struct comm *c,int i,GenmapScalar tolSquared)
{
  uint nPoints=mesh->elements.n;
  Point points=mesh->elements.ptr;

  sint rank=c->id,size=c->np;

  struct Point_private lastp=points[nPoints-1];
  lastp.proc=(rank+1)%size;

  struct array arr;
  array_init(struct Point_private,&arr,1);
  array_cat(struct Point_private,&arr,&lastp,1);

  struct crystal cr; crystal_init(&cr,c);
  sarray_transfer(struct Point_private,&arr,proc,1,&cr);
  crystal_free(&cr);

  uint n=arr.n; assert(n==1);
  lastp=((struct Point_private*)arr.ptr)[0];

  if(rank>0){
    GenmapScalar d=sqrDiff(lastp.x[i],points->x[i]);
    GenmapScalar dx=min(lastp.dx,points->dx)*tolSquared;
    if(d>dx) points->ifSegment=1;
  }

  array_free(&arr);

  return 0;
}

int findSegments(Mesh mesh,struct comm *c,GenmapScalar tol){
  int nDim=mesh->nDim,nVertex=mesh->nVertex;
  GenmapScalar tolSquared=tol*tol;

  parallel_sort(struct Point_private,&mesh->elements,x[0],
    genmap_gs_scalar,bin_sort,0,c);

  uint nPoints=mesh->elements.n;
  Point points=mesh->elements.ptr;

  buffer buf; buffer_init(&buf,1024);
  if(nDim==3)
    sarray_sort_2(struct Point_private,points,nPoints,x[1],3,x[2],3,&buf);
  else
    sarray_sort  (struct Point_private,points,nPoints,x[1],3,&buf);
  buffer_free(&buf);

  //TODO: load balance

  sint bin=1;
  if(nPoints==0)
    bin=0;

  sint rank=c->id,size=c->np; comm_ext old=c->c;
  struct comm nonZeroRanks;
#ifdef MPI
  MPI_Comm new; MPI_Comm_split(old,bin,rank,&new);
  comm_init(&nonZeroRanks,new);
  MPI_Comm_free(&new);
#else
  comm_init(&nonZeroRanks,1);
#endif

  rank=nonZeroRanks.id; size=nonZeroRanks.np;

  if(bin>0){
    slong out[2][1],buff[2][1],in[1];
    in[0]=nPoints;
    comm_scan(out,&nonZeroRanks,gs_long,gs_add,in,1,buff);
    slong start=out[0][0];

#if defined(GENMAP_DEBUG)
    printf("segments: rank=%d npts=%u start=%lld\n",rank,nPoints,start);
#endif

    sint i; for(i=0;i<nPoints;i++)
      points[i].ifSegment=0,points[i].proc=rank;

    int dim;
    for(dim=0;dim<nDim;dim++){
      findLocalSegments(mesh,dim,tolSquared);
      mergeSegments(mesh,&nonZeroRanks,dim,tolSquared);

#if defined(GENMAP_DEBUG)
    sint count=0;
    for(i=0;i<nPoints;i++)
      if(points[i].ifSegment>0)
        count++;

    in[0]=count;
    comm_allreduce(&nonZeroRanks,gs_long,gs_add,in,1,buff);
    if(rank==0)
      printf("locglob: %d %lld\n",dim+1,in[0]+1);
#endif
    }
  }

  comm_free(&nonZeroRanks);

  return 0;
}
