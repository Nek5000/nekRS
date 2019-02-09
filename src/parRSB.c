#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "genmap.h"
#include "genmap-impl.h"
#include "genmap-io.h"
#include "parRSB.h"

void fparRSB_partMesh(long long *egl, long long *vl, int *negl,
                      long long *eglcon, long long *vlcon, int *neglcon,
                      int *nve, int *options, int *comm, int *err) {
  *err = 1;

  GenmapCommExternal c;
#if defined(GENMAP_MPI)
  c = MPI_Comm_f2c(*comm);
#else
  c = 0;
#endif

  *err = parRSB_partMesh(egl, vl, negl,
                         eglcon, vlcon, *neglcon,
                         *nve, options, c);
}

int parRSB_partMesh(long long *egl, long long *vl, int *negl,
                    long long *eglcon, long long *vlcon, int neglcon,
                    int nve, int *options, MPI_Comm comm) {
  GenmapHandle h;
  GenmapInit(&h, comm, "interface");

  if(options[0] != 0) {
    h->dbgLevel = options[1];
    h->printStat = options[2];
  }

  // Check if negl is large enough
  GenmapLong neglcon_ = (GenmapLong) neglcon;
  GenmapGop(GenmapGetGlobalComm(h), &neglcon_, 1, GENMAP_LONG, GENMAP_SUM);
  GenmapInt negl_max = (GenmapInt)(neglcon_ / GenmapCommSize(
                                     GenmapGetGlobalComm(h))) + 1;
  if(negl_max > *negl) {
    printf("ERROR: negl to small to hold resulting partition!\n");
    return 1;
  }

  GenmapSetNLocalElements(h, neglcon);
  array_init(struct GenmapElement_private, &h->elementArray, neglcon);
  h->elementArray.n = neglcon;

  h->npts = neglcon * nve;
  h->nv = nve;

  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt_ = GenmapGetNLocalElements(h);
  comm_scan(out, &(GenmapGetGlobalComm(h)->gsComm), genmap_gs_long, gs_add,
            &lelt_, 1,
            buf);
  GenmapSetLocalStartIndex(h, out[0][0]);
  GenmapSetNGlobalElements(h, out[1][0]);

  GenmapElements e = GenmapGetElements(h);
  GenmapInt i, j;

  for(i = 0; i < neglcon; i++) {
    e[i].globalId = eglcon[i];
    for(j = 0; j < nve; j++) {
      e[i].vertices[j] = vlcon[i * nve + j];
    }
  }

  GenmapRSB(h);

  GenmapElements elements = GenmapGetElements(h);
  GenmapInt nv = h->nv;

  *negl = GenmapGetNLocalElements(h);

  for(i = 0; i < *negl; i++) {
    egl[i] = elements[i].globalId;
    for(j = 0; j < nv; j++) {
      vl[nv * i + j] = elements[i].vertices[j];
    }
  }

  if(h->printStat > 0) GenmapPartitionQuality(h);

  GenmapFinalize(h);

  return 0;
}
