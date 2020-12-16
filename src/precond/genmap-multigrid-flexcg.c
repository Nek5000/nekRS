#include <math.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

int flex_cg(genmap_handle h,GenmapComm c,mgData d,GenmapVector ri,
  int maxIter,int verbose,GenmapVector x)
{
  assert(x->size==ri->size);
  assert(x->size==GenmapGetNLocalElements(h));

  uint lelt=x->size;
  GenmapLong nelg=GenmapGetNGlobalElements(h);

  GenmapVector z0,z,dz,w,p,r;
  GenmapCreateVector(&z ,lelt);
  GenmapCreateVector(&w ,lelt);
  GenmapCreateVector(&r,lelt);
  GenmapCreateVector(&p ,lelt);
  GenmapCreateVector(&z0,lelt);
  GenmapCreateVector(&dz,lelt);

#define PREC 1
#define ORTH 1

  int rank=GenmapCommRank(c);

#if PREC
  if(rank==0 && verbose) printf("Using MG Prec.\n");
#endif

  uint i;
  for(i=0; i<lelt; i++)
    x->data[i]=0.0,r->data[i]=ri->data[i];

  GenmapCopyVector(z,r);
#if ORTH
  GenmapOrthogonalizebyOneVector(c,z,nelg);
#endif

  GenmapScalar den,alpha,beta,rz0,rz1=0,rz2,rr;

  rz1=GenmapDotVector(r,z);
  GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);
  if(GenmapCommRank(c)==0 && verbose)
    printf("rz1=%lf\n",rz1);

  GenmapCopyVector(p,z);

  i=0;
  while(i<maxIter && sqrt(rz1)>GENMAP_TOL){
    GenmapLaplacian(h,c,p->data,w->data);

    den=GenmapDotVector(p,w);
    GenmapGop(c,&den,1,GENMAP_SCALAR,GENMAP_SUM);

    alpha=rz1/den;

    GenmapAxpbyVector(x,x,1.0,p, alpha);
    GenmapAxpbyVector(r,r,1.0,w,-alpha);

    GenmapCopyVector(z0,z);
#if PREC
    mg_vcycle(z->data,r->data,d);
#else
    GenmapCopyVector(z,r);
#endif
#if ORTH
    GenmapOrthogonalizebyOneVector(c,z,nelg);
#endif

    rz0=rz1;

    rz1=GenmapDotVector(r,z);
    GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

    GenmapAxpbyVector(dz,z,1.0,z0,-1.0);
    rz2=GenmapDotVector(r,dz);
    GenmapGop(c,&rz2,1,GENMAP_SCALAR,GENMAP_SUM);
    beta=rz2/rz0;

    GenmapAxpbyVector(p,z,1.0,p,beta);
    i++;

    rr=GenmapDotVector(r,r);
    GenmapGop(c,&rr,1,GENMAP_SCALAR,GENMAP_SUM);
    if(rank==0 && verbose)
      printf("i=%d rr=%1.10e\n",i,sqrt(rr));
  }

  GenmapDestroyVector(z),GenmapDestroyVector(w);
  GenmapDestroyVector(p); GenmapDestroyVector(r);
  GenmapDestroyVector(z0),GenmapDestroyVector(dz);

  return i;
}
