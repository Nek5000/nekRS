#include <math.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

// Input z should be orthogonal to 1-vector, have unit norm.
// RQI should not change z.
int rqi(GenmapHandle h,GenmapComm c,mgData d,GenmapVector z,
  int iter,int verbose,GenmapVector y)
{
  assert(z->size==y->size);

  int ppfi;
  metric_tic(&c->gsc,PROJECTPF);
  ppfi=project_pf(h,c,d,z,20,verbose,y);
  metric_toc(&c->gsc,PROJECTPF);
  metric_acc(NPROJECTPF,ppfi);

  GenmapVector err; GenmapCreateVector(&err,z->size);
  GenmapLong nelg=GenmapGetNGlobalElements(h);
  int rank=GenmapCommRank(GenmapGetGlobalComm(h));

  uint i;
  for(i=0; i<iter; i++){
    GenmapScalar norm0=GenmapDotVector(y,y);
    GenmapGop(c,&norm0,1,GENMAP_SCALAR,GENMAP_SUM);
    GenmapScalar normi0=1.0/sqrt(norm0);

    GenmapAxpbyVector(z,z,0.0,y,normi0);
    GenmapOrthogonalizebyOneVector(h,c,z,nelg);

    metric_tic(&c->gsc,PROJECTPF);
    ppfi=project_pf(h,c,d,z,100,verbose,y);
    metric_toc(&c->gsc,PROJECTPF);
    metric_acc(NPROJECTPF,ppfi);
    GenmapOrthogonalizebyOneVector(h,c,y,nelg);

    //char fname[BUFSIZ];
    //sprintf(fname,"rqi_fiedler.%03d",i);
    //GenmapFiedlerDump(fname,y,h,c);

    GenmapScalar lambda=GenmapDotVector(y,z);
    GenmapGop(c,&lambda,1,GENMAP_SCALAR,GENMAP_SUM);

    GenmapAxpbyVector(err,y,1.0,z,-lambda);
    GenmapScalar norme=GenmapDotVector(err,err);
    GenmapGop(c,&norme,1,GENMAP_SCALAR,GENMAP_SUM);
    norme=sqrt(norme);

    GenmapScalar norm1=GenmapDotVector(y,y);
    GenmapGop(c,&norm1,1,GENMAP_SCALAR,GENMAP_SUM);
    GenmapScalar normi1=1.0/sqrt(norm1);
    metric_acc(END+i,norme*normi1);
  }

  GenmapDestroyVector(err);

  return i;
}
