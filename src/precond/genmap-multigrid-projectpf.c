#include <math.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#define MM 505 

int project_pf(GenmapHandle h,GenmapComm c,mgData d,GenmapVector ri,
  int maxIter,int verbose,GenmapVector x)
{
  assert(x->size==ri->size);
  assert(x->size==GenmapGetNLocalElements(h));

  uint lelt=x->size;
  GenmapLong nelg=GenmapGetNGlobalElements(h);

  GenmapVector z0,z,dz,w,p,r;
  GenmapCreateVector(&z ,lelt);
  GenmapCreateVector(&w ,lelt);
  GenmapCreateVector(&r ,lelt);
  GenmapCreateVector(&p ,lelt);
  GenmapCreateVector(&z0,lelt);
  GenmapCreateVector(&dz,lelt);

  assert(maxIter<MM);
  double *P; GenmapCalloc(lelt*MM,&P);
  double *W; GenmapCalloc(lelt*MM,&W);

  uint i;
  for(i=0; i<lelt; i++)
    x->data[i]=0.0,r->data[i]=ri->data[i];

  metric_tic(&d->c,PRECONVCYCLE);
  mg_vcycle(z->data,r->data,d);
  metric_toc(&d->c,PRECONVCYCLE);
  GenmapOrthogonalizebyOneVector(h,c,z,nelg);
  GenmapCopyVector(p,z);

  GenmapScalar rz1=GenmapDotVector(r,z);
  GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

  GenmapScalar rr=GenmapDotVector(r,r);
  GenmapGop(c,&rr,1,GENMAP_SCALAR,GENMAP_SUM);

  GenmapScalar alpha,beta,rz0,rz2,scale;

  int rank=GenmapCommRank(c);
  if(rank==0 && verbose)
    printf("projectpf initial rr=%g rz1=%g\n",sqrt(rr),sqrt(rz1));

  double tol=sqrt(rr)*1e-7;

  i=0; uint j,k;
  while(i<maxIter){
    metric_tic(&c->gsc,LAPLACIAN);
    GenmapLaplacian(h,c,p,w);
    metric_toc(&c->gsc,LAPLACIAN);

    GenmapScalar den=GenmapDotVector(p,w);
    GenmapGop(c,&den,1,GENMAP_SCALAR,GENMAP_SUM);
    alpha=rz1/den;

    scale=1.0/sqrt(den);
    for(j=0; j<lelt; j++){
      W[i*lelt+j]=scale*w->data[j];
      P[i*lelt+j]=scale*p->data[j];
    }

    GenmapAxpbyVector(x,x,1.0,p, alpha);
    GenmapAxpbyVector(r,r,1.0,w,-alpha);

    rr=GenmapDotVector(r,r);
    GenmapGop(c,&rr,1,GENMAP_SCALAR,GENMAP_SUM);
    if(sqrt(rr)<tol)
      break;

    GenmapCopyVector(z0,z);

    GenmapScalar norm0=GenmapDotVector(z,z);
    GenmapGop(c,&norm0,1,GENMAP_SCALAR,GENMAP_SUM);

    metric_tic(&c->gsc,PRECONVCYCLE);
    mg_vcycle(z->data,r->data,d);
    metric_toc(&c->gsc,PRECONVCYCLE);

    GenmapScalar norm1=GenmapDotVector(z,z);
    GenmapGop(c,&norm1,1,GENMAP_SCALAR,GENMAP_SUM);

    if(rank==0 && verbose)
      printf("vcycle i=%02d norm0=%g norm1=%g\n",i,sqrt(norm0),
        sqrt(norm1));

    rz0=rz1;
    GenmapOrthogonalizebyOneVector(h,c,z,nelg);
    rz1=GenmapDotVector(r,z);
    GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

    GenmapAxpbyVector(dz,z,1.0,z0,-1.0);
    rz2=GenmapDotVector(r,dz);
    GenmapGop(c,&rz2,1,GENMAP_SCALAR,GENMAP_SUM);

    if(rank==0 && verbose)
      printf("projectpf i=%d rr=%g rz1=%g rz2=%g alpha=%g\n",
        i,sqrt(rr),sqrt(rz1),sqrt(rz2),alpha);

    beta=rz2/rz0;
    GenmapAxpbyVector(p,z,1.0,p,beta);

    i++;

    for(k=0; k<lelt; k++)
      P[(MM-1)*lelt+k]=0.0;

    for(j=0; j<i; j++){
      double a=0.0;
      for(k=0; k<lelt; k++)
        a+=W[j*lelt+k]*p->data[k];
      GenmapGop(c,&a,1,GENMAP_SCALAR,GENMAP_SUM);
      for(k=0; k<lelt; k++)
        P[(MM-1)*lelt+k]+=a*P[j*lelt+k];
    }

    for(k=0; k<lelt; k++)
      p->data[k]-=P[(MM-1)*lelt+k];
  }

  GenmapDestroyVector(z),GenmapDestroyVector(w);
  GenmapDestroyVector(p); GenmapDestroyVector(r);
  GenmapDestroyVector(z0),GenmapDestroyVector(dz);

  GenmapFree(P); GenmapFree(W);

  return i;
}

int project_pf_lvl(GenmapHandle h,GenmapComm c,mgData d,GenmapScalar *ri,
  int maxIter,int verbose,int lvl_start,GenmapScalar *xo)
{
  assert(lvl_start<d->nlevels-1);

  uint  lelt=d->level_off[lvl_start+1]-d->level_off[lvl_start];
  slong nelg=lelt;
  GenmapGop(c,&nelg,1,GENMAP_LONG,GENMAP_SUM);

  GenmapVector z0,z,dz,w,p,r,x;
  GenmapCreateVector(&z ,lelt);
  GenmapCreateVector(&w ,lelt);
  GenmapCreateVector(&r ,lelt);
  GenmapCreateVector(&x ,lelt);
  GenmapCreateVector(&p ,lelt);
  GenmapCreateVector(&z0,lelt);
  GenmapCreateVector(&dz,lelt);

  assert(maxIter<MM);
  double *P; GenmapCalloc(lelt*MM,&P);
  double *W; GenmapCalloc(lelt*MM,&W);

  uint i;
  for(i=0; i<lelt; i++)
    x->data[i]=0.0,r->data[i]=ri[i];

  metric_tic(&d->c,PRECONVCYCLE);
  mg_vcycle_lvl(z->data,r->data,d,lvl_start);
  metric_toc(&d->c,PRECONVCYCLE);
  GenmapOrthogonalizebyOneVector(h,c,z,nelg);
  GenmapCopyVector(p,z);

  GenmapScalar rz1=GenmapDotVector(r,z);
  GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

  GenmapScalar rr=GenmapDotVector(r,r);
  GenmapGop(c,&rr,1,GENMAP_SCALAR,GENMAP_SUM);
  double tol=sqrt(rr)*1e-7;

  GenmapScalar alpha,beta,rz0,rz2,scale;
  uint j,k;

  csr_mat M=d->levels[lvl_start]->M;
  buffer buf; buffer_init(&buf,1024);

  i=0;
  while(i<maxIter){
    //metric_tic(&c->gsc,LAPLACIAN);
    //GenmapLaplacian(h,c,p,w);
    //metric_toc(&c->gsc,LAPLACIAN);
    metric_tic(&c->gsc,PRECONAX);
    csr_mat_gather(M,M->gsh,p->data,d->buf,&buf);
    csr_mat_apply(w->data,M,d->buf);
    metric_toc(&c->gsc,PRECONAX);

    GenmapScalar den=GenmapDotVector(p,w);
    GenmapGop(c,&den,1,GENMAP_SCALAR,GENMAP_SUM);
    alpha=rz1/den;

    scale=1.0/sqrt(den);
    for(j=0; j<lelt; j++){
      W[i*lelt+j]=scale*w->data[j];
      P[i*lelt+j]=scale*p->data[j];
    }

    GenmapAxpbyVector(x,x,1.0,p, alpha);
    GenmapAxpbyVector(r,r,1.0,w,-alpha);

    rr=GenmapDotVector(r,r);
    GenmapGop(c,&rr,1,GENMAP_SCALAR,GENMAP_SUM);
    if(sqrt(rr)<tol)
      break;

    GenmapCopyVector(z0,z);

    GenmapScalar norm0=GenmapDotVector(z,z);
    GenmapGop(c,&norm0,1,GENMAP_SCALAR,GENMAP_SUM);

    metric_tic(&c->gsc,PRECONVCYCLE);
    mg_vcycle_lvl(z->data,r->data,d,lvl_start);
    metric_toc(&c->gsc,PRECONVCYCLE);

    GenmapScalar norm1=GenmapDotVector(z,z);
    GenmapGop(c,&norm1,1,GENMAP_SCALAR,GENMAP_SUM);

    rz0=rz1;
    GenmapOrthogonalizebyOneVector(h,c,z,nelg);
    rz1=GenmapDotVector(r,z);
    GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

    GenmapAxpbyVector(dz,z,1.0,z0,-1.0);
    rz2=GenmapDotVector(r,dz);
    GenmapGop(c,&rz2,1,GENMAP_SCALAR,GENMAP_SUM);

    beta=rz2/rz0;
    GenmapAxpbyVector(p,z,1.0,p,beta);

    i++;

    for(k=0; k<lelt; k++)
      P[(MM-1)*lelt+k]=0.0;

    for(j=0; j<i; j++){
      double a=0.0;
      for(k=0; k<lelt; k++)
        a+=W[j*lelt+k]*p->data[k];
      GenmapGop(c,&a,1,GENMAP_SCALAR,GENMAP_SUM);
      for(k=0; k<lelt; k++)
        P[(MM-1)*lelt+k]+=a*P[j*lelt+k];
    }

    for(k=0; k<lelt; k++)
      p->data[k]-=P[(MM-1)*lelt+k];
  }

  for(i=0; i<lelt; i++)
    xo[i]=x->data[i];

  GenmapDestroyVector(z ); GenmapDestroyVector(w);
  GenmapDestroyVector(p ); GenmapDestroyVector(r);
  GenmapDestroyVector(z0); GenmapDestroyVector(dz);
  GenmapDestroyVector(x);

  GenmapFree(P); GenmapFree(W);

  buffer_free(&buf);

  return i;
}
