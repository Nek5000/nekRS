#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>
#include <sort.h>
#include <parRSB.h>

void genmap_rsb(genmap_handle h,int verbose){
  int max_iter=50;
  int max_pass=50;

  GenmapComm local_c=GenmapGetLocalComm(h);
  struct comm *lc=&local_c->gsc;

  GenmapComm global_c=GenmapGetGlobalComm(h);
  struct comm *gc=&global_c->gsc;

  GenmapScan(h,local_c);
  crystal_init(&h->cr,lc);

  GenmapScan(h, GenmapGetLocalComm(h));
  uint nelt=GenmapGetNLocalElements(h);
  GenmapElements e=GenmapGetElements(h);
  GenmapInt i;
  for(i=0; i<nelt; i++)
    e[i].globalId0=GenmapGetLocalStartIndex(h)+i+1;

  buffer buf; buffer_init(&buf,1024);

  int nve=h->nv;
  int ndim=(nve==8)?3:2;
  int level=0;
  int np;

  while(GenmapCommSize(GenmapGetLocalComm(h))>1){
    local_c=GenmapGetLocalComm(h);
    lc=&local_c->gsc;
    np=lc->np;

#if defined(GENMAP_PAUL)
    int global=1;
#else
    int global=(np==gc->np);
#endif

#if defined(GENMAP_RCB_PRE_STEP)
    /* Run RCB pre-step */
    metric_tic(lc,RCB);
    rcb(lc,h->elements,ndim);
    rcb_local(h->elements,0,h->elements->n,ndim,&buf);
    metric_toc(lc,RCB);
#else
    /* Sort by global id otherwise */
    parallel_sort(struct rsb_element,h->elements,globalId0,gs_long,0,1,lc);
#endif

    /* Initialize the laplacian */
    metric_tic(lc,LAPLACIANSETUP0);
    GenmapInitLaplacianWeighted(h,local_c);
    metric_toc(lc,LAPLACIANSETUP0);

    /* Run fiedler */
    metric_tic(lc,FIEDLER);
    int ipass=0,iter;
    do {
#if defined(GENMAP_LANCZOS)
      iter=GenmapFiedlerLanczos(h,local_c,max_iter,global);
#elif defined(GENMAP_RQI)
      iter=GenmapFiedlerRQI(h,local_c,max_iter,global);
#endif
      metric_acc(NFIEDLER,iter);
      global=0;
    } while(++ipass<max_pass && iter==max_iter);
    metric_toc(lc,FIEDLER);

    /* Bisect */
    metric_tic(lc,BISECT);
    parallel_sort(struct rsb_element,h->elements,fiedler,gs_double,0,1,lc);
    int bin=1;
    if(lc->id<(np+1)/2)
      bin=0;
    // FIXME: Ugly
    GenmapSplitComm(h,&local_c,bin);
    GenmapSetLocalComm(h,local_c);
    lc=&local_c->gsc;
    GenmapScan(h,local_c);
    metric_toc(lc,BISECT);

    metric_push_level();
    level++;
  }

  /* Check if Fidler converged */
  sint converged=1;
  for(i=0; i<metric_get_levels(); i++){
    int val=(int)metric_get_value(i,NFIEDLER);
    if(val>=max_pass*max_iter){
      converged=0; break;
    }
  }
  sint bfr;
  comm_allreduce(gc,gs_int,gs_min,&converged,1,&bfr);// min
  if(converged==0 && gc->id==0)
    printf("\tWARNING: Lanczos failed to converge while partitioning!\n");

  /* Check for disconnected components */
  GenmapInitLaplacianWeighted(h,global_c);

  sint discon=is_disconnected(gc,global_c->gsw,&global_c->buf,nelt,nve);
  if(discon>0 && gc->id==0)
    printf("\tWarning: There are disconnected components!\n");

  crystal_free(&h->cr);
  buffer_free(&buf);
}
