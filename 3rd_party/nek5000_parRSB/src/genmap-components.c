#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>
#include <parRSB.h>

/* Find the number of disconnected components */
sint is_disconnected(struct comm *c,struct gs_data *gsh,buffer *buf,
  uint nelt_,uint nv)
{
  slong nelt=nelt_;

  GenmapLong *p;
  GenmapMalloc(nelt*nv,&p);

  slong out[2][1],buff[2][1];
  comm_scan(out,c,gs_long,gs_add,&nelt,1,buff);
  slong nelg=out[1][0];

  uint e;
  for(e=0; e<nelt*nv; e++)
    p[e]=0;

  if(c->id==0){
    for(e=0; e<nv; e++)
      p[e]=1;
  }

  slong nnz1=1,nnz0,nnzb;
  uint d;
  do{
    nnz0=nnz1,nnz1=0;

    gs(p,gs_long,gs_add,0,gsh,buf);

    for(e=0; e<nelt; e++){
      for(d=0; d<nv; d++)
        if(p[e*nv+d]>0){ nnz1++; break; }

      if(d<nv)
        for(d=0; d<nv; d++) p[e*nv+d]=1;
    }

    comm_allreduce(c,gs_long,gs_add,&nnz1,1,&nnzb);
  } while(nnz1>nnz0);

  GenmapFree(p);

  return (nnz1<nelg);
}
