#include "gencon-impl.h"

int setGlobalID(Mesh mesh,struct comm *c){
  uint nPoints=mesh->elements.n;
  Point points=mesh->elements.ptr;

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

  if(bin==1){
    slong count=0;
    sint i;
    for(i=0;i<nPoints;i++)
      if(points[i].ifSegment) count++;

    slong out[2][1],buff[2][1],in[1];
    in[0]=count+!rank;
    comm_scan(out,&nonZeroRanks,gs_long,gs_add,in,1,buff);
    slong start=out[0][0];
#if defined(GENMAP_DEBUG)
    printf("rank=%d start=%lld size=%lld\n",rank,start,in[0]);
#endif

    start-=(rank>0?1:0);
    count=0;
    for(i=0;i<nPoints;i++){
      if(points[i].ifSegment) count++;
      points[i].globalId=start+count;
    }
  }

  comm_free(&nonZeroRanks);

  return 0;
}

int sendBack(Mesh mesh,struct comm *c){
  struct crystal cr; crystal_init(&cr,c);
  sarray_transfer(struct Point_private,&mesh->elements,origin,0,&cr);
  crystal_free(&cr);

  buffer buf; buffer_init(&buf,1024);
  sarray_sort(struct Point_private,mesh->elements.ptr,mesh->elements.n,
    sequenceId,1,&buf);
  buffer_free(&buf);

  return 0;
}
