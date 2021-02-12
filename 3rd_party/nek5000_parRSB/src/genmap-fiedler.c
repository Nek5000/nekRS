#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>
//
//TODO: use a separate function to generate init vector
//
int GenmapFiedlerRQI(genmap_handle h,GenmapComm c,int max_iter,int global)
{
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapVector initVec;
  GenmapCreateVector(&initVec,lelt);

  GenmapElements elements = GenmapGetElements(h);
  GenmapInt i;
  if(global>0){
#if defined(GENMAP_PAUL)
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = GenmapGetLocalStartIndex(h) + i + 1;
    }
#else
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = (GenmapScalar) elements[i].globalId;
    }
#endif
  }else{
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }

  int verbose=h->verbose_level;
  struct comm *gsc=&c->gsc;

  GenmapOrthogonalizebyOneVector(c,initVec,GenmapGetNGlobalElements(h));

  GenmapScalar norm=GenmapDotVector(initVec,initVec);
  GenmapGop(c,&norm,1,GENMAP_SCALAR,GENMAP_SUM);
  if(verbose>0 && gsc->id==0)
    printf("RQI, |init| = %g\n",sqrt(norm));

  norm=1.0/sqrt(norm);
  GenmapScaleVector(initVec,initVec,norm);

  metric_tic(gsc,LAPLACIANSETUP1);
  GenmapInitLaplacian(h,c);
  metric_toc(gsc,LAPLACIANSETUP1);

  metric_tic(gsc,PRECONSETUP);
  mgData d; mgSetup(c,c->M,&d); d->h=h;
  metric_toc(gsc,PRECONSETUP);

  GenmapVector y; GenmapCreateZerosVector(&y,lelt);
  metric_tic(gsc,RQI);
  int iter=rqi(h,c,d,initVec,max_iter,verbose,y);
  metric_acc(NRQI,iter);
  metric_toc(gsc,RQI);

  mgFree(d);

  GenmapScalar lNorm = 0;
  for(i = 0; i < lelt; i++)
    lNorm+=y->data[i]*y->data[i];
  GenmapGop(c,&lNorm,1,GENMAP_SCALAR,GENMAP_SUM);

  GenmapScaleVector(y,y,1./sqrt(lNorm));
  for(i = 0; i < lelt; i++)
    elements[i].fiedler=y->data[i];

  GenmapDestroyVector(y);

  return iter;
}

int GenmapFiedlerLanczos(genmap_handle h,GenmapComm c,int max_iter,
  int global)
{
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapVector initVec, alphaVec, betaVec;

  GenmapCreateVector(&initVec, GenmapGetNLocalElements(h));
  GenmapElements elements = GenmapGetElements(h);

  GenmapInt i;
#if defined(GENMAP_PAUL)
  if(global>0){
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = GenmapGetLocalStartIndex(h) + i + 1;
    }
  }else{
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }
#else
  if(global>0){
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = (GenmapScalar) elements[i].globalId;
    }
  }else{
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }
#endif

  GenmapCreateVector(&alphaVec,max_iter);
  GenmapCreateVector(&betaVec,max_iter-1);
  GenmapVector *q = NULL;

  GenmapOrthogonalizebyOneVector(c,initVec,GenmapGetNGlobalElements(h));
  GenmapScalar rtr = GenmapDotVector(initVec, initVec);
  GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScalar rni = 1.0 / sqrt(rtr);
  GenmapScaleVector(initVec, initVec, rni);

#if defined(GENMAP_PAUL)
  int iter=GenmapLanczosLegendary(h,c,initVec,max_iter,&q,alphaVec,betaVec);
#else
  int iter=GenmapLanczos(h,c,initVec,max_iter,&q,alphaVec,betaVec);
#endif

  GenmapVector evLanczos, evTriDiag;
  GenmapCreateVector(&evTriDiag, iter);

#if defined(GENMAP_PAUL)
  /* Use TQLI and find the minimum eigenvalue and associated vector */
  GenmapVector *eVectors, eValues;
  GenmapTQLI(h, alphaVec, betaVec, &eVectors, &eValues);

  GenmapScalar eValMin = fabs(eValues->data[0]);
  GenmapInt eValMinI = 0;
  for(i = 1; i < iter; i++) {
    if(fabs(eValues->data[i]) < eValMin) {
      eValMin = fabs(eValues->data[i]);
      eValMinI = i;
    }
  }
  GenmapCopyVector(evTriDiag, eVectors[eValMinI]);
#else
  GenmapVector init;
  GenmapCreateVector(&init, iter);
  for(int i = 0; i < iter; i++) {
    init->data[i] = i + 1.0;
  }
  GenmapScalar avg = 0.5 * iter * (1.0 + iter) / iter;
  for(int i = 0; i < iter; i++) {
    init->data[i] -= avg;
  }
  GenmapInvPowerIter(evTriDiag, alphaVec, betaVec, init, 100);
  GenmapDestroyVector(init);
#endif

  GenmapInt j;
  GenmapCreateZerosVector(&evLanczos, lelt);
  for(i = 0; i < lelt; i++) {
    for(j = 0; j < iter; j++) {
      evLanczos->data[i] += q[j]->data[i] * evTriDiag->data[j];
    }
  }

  GenmapScalar lNorm = 0;
  for(i = 0; i < lelt; i++) {
    lNorm += evLanczos->data[i] * evLanczos->data[i];
  }

  GenmapGop(c, &lNorm, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScaleVector(evLanczos, evLanczos, 1. / sqrt(lNorm));
  for(i = 0; i < lelt; i++) {
    elements[i].fiedler = evLanczos->data[i];
  }

  GenmapDestroyVector(initVec);
  GenmapDestroyVector(alphaVec);
  GenmapDestroyVector(betaVec);
  GenmapDestroyVector(evLanczos);
  GenmapDestroyVector(evTriDiag);
#if defined(GENMAP_PAUL)
  GenmapDestroyVector(eValues);
  for(i = 0; i < iter; i++) {
    GenmapDestroyVector(eVectors[i]);
  }
  GenmapFree(eVectors);
#endif

#if defined(GENMAP_PAUL)
  for(i = 0; i < iter + 1; i++) {
    GenmapDestroyVector(q[i]);
  }
#else
  for(i = 0; i < iter; i++) {
    GenmapDestroyVector(q[i]);
  }
#endif
  GenmapFree(q);

  return iter;
}

#define write_T(dest,val,T,nunits) do{\
  memcpy(dest,&(val),sizeof(T)*nunits);\
  dest+=sizeof(T)*nunits;\
} while(0)

int GenmapFiedlerDump(const char *fname,genmap_handle h,GenmapComm comm)
{
  struct comm *c=&comm->gsc;

  MPI_File file;
  int err=MPI_File_open(c->c,fname,MPI_MODE_CREATE|MPI_MODE_WRONLY,
                        MPI_INFO_NULL,&file);
  uint rank=c->id;
  if(err!=0 && rank==0){
    fprintf(stderr,"%s:%d Error opening file %s for writing.\n",__FILE__,__LINE__,fname);
    return err;
  }

  slong nelt=GenmapGetNLocalElements(h);
  slong out[2][1],buf[2][1];
  comm_scan(out,c,gs_long,gs_add,&nelt,1,buf);
  slong start=out[0][0];
  slong nelgt=out[1][0];

  int ndim=(h->nv==8)?3:2;
  uint write_size=(ndim+1)*nelt*sizeof(double);
  if(rank==0)
    write_size+=sizeof(long)+sizeof(int); // for nelgt and ndim

  char *pbuf,*pbuf0;
  pbuf=pbuf0=(char*)calloc(write_size,sizeof(char));
  if(rank==0){
    write_T(pbuf0,nelgt,long,1);
    write_T(pbuf0,ndim ,int ,1);
  }

  GenmapElements elm=GenmapGetElements(h);
  uint i;
  for(i=0; i<nelt; i++){
    write_T(pbuf0,elm[i].coord[0],double,ndim);
    write_T(pbuf0,elm[i].fiedler ,double,1   );
  }

  MPI_Status st;
  err=MPI_File_write_ordered(file,pbuf,write_size,MPI_BYTE,&st);
  if(err!=0 && rank==0){
    fprintf(stderr,"%s:%d Error opening file %s for writing.\n",__FILE__,__LINE__,fname);
    return err;
  }

  err+=MPI_File_close(&file);
  MPI_Barrier(c->c);

  free(pbuf);

  return err;
}

int GenmapVectorDump(const char *fname,GenmapScalar *y,uint size,
  struct comm *c)
{
  MPI_File file;
  int err=MPI_File_open(c->c,fname,MPI_MODE_CREATE|MPI_MODE_WRONLY,
                        MPI_INFO_NULL,&file);
  uint rank=c->id;
  if(err!=0 && rank==0){
    fprintf(stderr,"%s:%d Error opening file %s for writing.\n",__FILE__,__LINE__,fname);
    return err;
  }

  slong nelt=size;
  slong out[2][1],buf[2][1];
  comm_scan(out,c,gs_long,gs_add,&nelt,1,buf);
  slong start=out[0][0];
  slong nelgt=out[1][0];

  uint write_size=nelt*sizeof(double);
  if(rank==0)
    write_size+=sizeof(long); // nelgt

  char *pbuf,*pbuf0;
  pbuf=pbuf0=(char*)calloc(write_size,sizeof(char));

  if(rank==0){
    write_T(pbuf0,nelgt,long,1);
  }

  uint i;
  for(i=0; i<nelt; i++){
    write_T(pbuf0,y[i],double,1);
  }

  MPI_Status st;
  err=MPI_File_write_ordered(file,pbuf,write_size,MPI_BYTE,&st);
  if(err!=0 && rank==0){
    fprintf(stderr,"%s:%d Error opening file %s for writing.\n",__FILE__,__LINE__,fname);
    return err;
  }

  err+=MPI_File_close(&file);
  MPI_Barrier(c->c);

  free(pbuf);

  return err;
}

#undef write_T
