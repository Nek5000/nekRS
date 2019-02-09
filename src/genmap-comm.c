#include "genmap-impl.h"
//
// GenmapComm
//
int GenmapCreateComm(GenmapComm *c, GenmapCommExternal ce) {
  GenmapMalloc(1, c);
  comm_init(&(*c)->gsComm, ce);
  (*c)->verticesHandle = NULL;
  (*c)->laplacianWeights = NULL;
  buffer_init(&(*c)->buf, 1024);
  return 0;
}

int GenmapDestroyComm(GenmapComm c) {
  if(&c->buf)
    buffer_free(&c->buf);
  if(c->verticesHandle)
    gs_free(c->verticesHandle);
  if(c->laplacianWeights)
    GenmapFree(c->laplacianWeights);
  if(&c->gsComm)
    comm_free(&c->gsComm);
  GenmapFree(c);

  return 0;
}

int GenmapCommSize(GenmapComm c) {
  return (int) c->gsComm.np;
}

int GenmapCommRank(GenmapComm c) {
  return (int) c->gsComm.id;
}

GenmapComm GenmapGetLocalComm(GenmapHandle h) {
  return h->local;
}

void GenmapSetLocalComm(GenmapHandle h, GenmapComm c) {
  h->local = c;
}

GenmapComm GenmapGetGlobalComm(GenmapHandle h) {
  return h->global;
}

void GenmapSetGlobalComm(GenmapHandle h, GenmapComm c) {
  h->global = c;
}

int GenmapGop(GenmapComm c, void *v, GenmapInt size,
              GenmapDataType type, GenmapInt op) {
#ifdef GENMAP_MPI
  if(op == GENMAP_SUM) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_SUM, c->gsComm.c);
  } else if(op == GENMAP_MAX) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_MAX, c->gsComm.c);
  } else if(op == GENMAP_MIN) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_MIN, c->gsComm.c);
  }
#endif
  return 0;
}
