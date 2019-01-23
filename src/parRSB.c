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
                      int *nve, int *comm, int *err) {
  *err = 1;
  setbuf(stdout, NULL);

  GenmapCommExternal c;
#if defined(GENMAP_MPI)
  c = MPI_Comm_f2c(*comm);
#else
  c = 0;
#endif

  *err = parRSB_partMesh(egl, vl, negl,
                         eglcon, vlcon, *neglcon,
                         *nve, c);
}

int parRSB_partMesh(long long *egl, long long *vl, int *negl,
                    long long *eglcon, long long *vlcon, int neglcon,
                    int nve, MPI_Comm comm) {
  GenmapHandle h;
  GenmapInit(&h, comm, "interface");

  // Check if negl is large enough
  GenmapLong neglcon_ = (GenmapLong) neglcon;
  GenmapGop(h->global, &neglcon_, 1, GENMAP_LONG, GENMAP_SUM);
  GenmapInt negl_max = (GenmapInt)(neglcon_ / GenmapNp(h->global)) + 1;
  if(negl_max > *negl) {
    printf("ERROR: negl to small to hold resulting partition!\n");
    return 1;
  }

  h->header->lelt = neglcon;
  h->header->npts = neglcon * nve;
  h->header->nv = nve;
  h->header->ndim = (nve == 8) ? 3 : 2;

  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt_ = h->header->lelt;
  comm_scan(out, &(h->global->gsComm), genmap_gs_long, gs_add, &lelt_, 1,
            buf);
  h->header->start = out[0][0];
  h->header->nel = out[1][0];

  array_init(struct GenmapElement_private, &h->elementArray, neglcon);
  h->elementArray.n = neglcon;

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
  GenmapInt nv = h->header->nv;

  *negl = h->header->lelt;

  for(i = 0; i < *negl; i++) {
    egl[i] = elements[i].globalId;
    for(j = 0; j < nv; j++) {
      vl[nv * i + j] = elements[i].vertices[j];
    }
  }

  GenmapPartitionQuality(h);

  GenmapFinalize(h);

  return 0;
}
