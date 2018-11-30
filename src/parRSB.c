#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "genmap.h"
#include "genmap-impl.h"
#include "genmap-fortran-name.h"
#include "genmap-io.h"
#include "parRSB.h"

void fparRSB_partMesh(long long *egl   , long long *vl   , int *negl, 
                      long long *eglcon, long long *vlcon, int *neglcon,
                      int *nve, int *comm, int *err) 
{
  *err = 1;
  setbuf(stdout, NULL);

  GenmapCommExternal c;
#if defined(GENMAP_MPI)
  c = MPI_Comm_f2c(*comm);
#else
  c = 0;
#endif

  *err = parRSB_partMesh(egl   , vl   , negl,
                         eglcon, vlcon, neglcon,
                         nve, c);
}

int parRSB_partMesh(GenmapLong *egl, GenmapLong *vl, GenmapInt *negl,
                    GenmapLong *eglcon, GenmapLong *vlcon, GenmapInt *neglcon,
                    GenmapInt  *nve, GenmapCommExternal comm) 
{

/*
 egl     (out)    ... Array of local global element IDs
 vl      (out)    ... Array of vertex IDs for all local elements
 negl    (in/out) ... Number of elements of egl array / Number of elements

 eglcon  (in)     ... Array of global element IDs
 vlcon   (in)     ... Array vertex IDs for all elements in eglcon
 nve     (in)     ... Array of vertices for all elements in eglcon
 neglcon (in)     ... Number of elements of eglcon array
 comm    (in)     ... Communicator
 ierr    (out)    ... Error flag
*/

  GenmapHandle h;
  h->Read = NULL;

  // Check if negl is large enough
  GenmapLong neglcon_ = (GenmapLong) neglcon;
  GenmapGop(h->global, &neglcon_, 1, GENMAP_LONG, GENMAP_SUM);
  GenmapInt rank = GenmapId(h->global);
  GenmapInt size = GenmapNp(h->global);
  GenmapInt localElements = (GenmapInt)(neglcon_ / size);
  GenmapInt nrem = (GenmapInt)(neglcon_ - localElements * size);
  localElements = localElements + rank < nrem;
  if(localElements > *negl) {
    printf("parameter negl is not large enough to hold the resulting partition.\n");
    return 1;
  }

  GenmapCreateHeader(&h->header);
  h->header->lelt = *neglcon;
  h->header->npts = (*neglcon) * (*nve);
  h->header->nv = *nve;
  h->header->ndim = (*nve == 8) ? 3 : 2;

  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt_ = h->header->lelt;
  comm_scan(out, &(h->global->gsComm), genmap_gs_long, gs_add, &lelt_, 1,
            buf);
  h->header->start = out[0][0];
  h->header->nel = out[1][0];

  array_init(struct GenmapElement_private, &h->elementArray, *neglcon);
  h->elementArray.n = *neglcon;

  GenmapElements e = GenmapGetElements(h);

  for(int i = 0; i < *neglcon; i++) {
    e[i].globalId = eglcon[i];
    for(int j = 0; j<(*nve); j++) {
      e[i].vertices[j] = vlcon[i*(*nve)+j];
    }
  }

  GenmapRSB(h);

  GenmapElements elements = GenmapGetElements(h);
  GenmapInt nv = h->header->nv;

  for(int i = 0; i < h->header->lelt; i++) {
    egl[i] = elements[i].globalId;
    for(int j = 0; j < nv; j++) {
      vl[nv * i + j] = elements[i].vertices[j];
    }
  }

  *negl = h->header->lelt;

  GenmapPartitionQuality(h);

  GenmapFinalize(h);

  return 0;
}
