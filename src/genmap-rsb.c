#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>

void GenmapRSB(GenmapHandle h,int verbose){
  int maxIter=50;
  int npass  =50;

  GenmapInt i;
  GenmapElements e = GenmapGetElements(h);
  GenmapScan(h, GenmapGetLocalComm(h));
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    e[i].globalId =GenmapGetLocalStartIndex(h)+i+1;
    e[i].globalId0=GenmapGetLocalStartIndex(h)+i+1;
  }

  GenmapComm global_c=GenmapGetGlobalComm(h);
  GenmapComm local_c =GenmapGetLocalComm(h);
  int rank=GenmapCommRank(global_c);

  if(rank == 0 && h->dbgLevel > 0)
    printf("running RSB ");
  fflush(stdout);

  crystal_init(&(h->cr), &(h->local->gsc));
  buffer buf0 = null_buffer;

  int nve =h->nv;
  int ndim=(nve==8)?3:2;
  int level=0;

  metric_init(); // init metrics

  while(GenmapCommSize(GenmapGetLocalComm(h)) > 1){
    GenmapComm local_c=GenmapGetLocalComm(h);
    struct comm *gsc=&local_c->gsc;
    GenmapInt np=gsc->np;

    metric_tic(gsc,RSB);

#if defined(GENMAP_PAUL)
    int global=1;
#else
    int global=(np==GenmapCommSize(GenmapGetGlobalComm(h)));
#endif

    int ipass=0,iter;
    do{
      metric_tic(gsc,FIEDLER);
#if defined(GENMAP_LANCZOS)
      iter=GenmapFiedlerLanczos(h,local_c,maxIter,global);
#elif defined(GENMAP_RQI)
      iter=GenmapFiedlerRQI(h,local_c,maxIter,global);
#endif
      metric_toc(gsc,FIEDLER);
      metric_acc(NFIEDLER,iter);

      global=0;
    }while(++ipass<npass && iter==maxIter);

    GenmapBinSort(h, GENMAP_FIEDLER, &buf0);
    metric_toc(gsc,RSB);

    GenmapInt id=GenmapCommRank(local_c);
    int bin=1; if(id<(np+1)/2) bin=0;
    GenmapSplitComm(h,&local_c,bin);
    GenmapSetLocalComm(h,local_c);

#if defined(GENMAP_PAUL)
    GenmapBinSort(h,GENMAP_GLOBALID,&buf0);
#endif

    level++;
    metric_push_level();
  }

  //metric_print(&global_c->gsc);
  metric_finalize();

  crystal_free(&(h->cr));
  buffer_free(&buf0);
}
