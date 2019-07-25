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
  buffer_free(&c->buf);
  if(c->verticesHandle)
    gs_free(c->verticesHandle);
  if(c->laplacianWeights)
    GenmapFree(c->laplacianWeights);
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
  if(op == GENMAP_SUM) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_SUM, c->gsComm.c);
  } else if(op == GENMAP_MAX) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_MAX, c->gsComm.c);
  } else if(op == GENMAP_MIN) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_MIN, c->gsComm.c);
  }
  return 0;
}

int GenmapReduce(GenmapComm c, void *out, void *in, GenmapInt size,
                 GenmapDataType type, GenmapInt op) {
  if(op == GENMAP_SUM) {
    MPI_Reduce(in, out, size, type, MPI_SUM, 0, c->gsComm.c);
  } else if(op == GENMAP_MAX) {
    MPI_Reduce(in, out, size, type, MPI_MAX, 0, c->gsComm.c);
  } else if(op == GENMAP_MIN) {
    MPI_Reduce(in, out, size, type, MPI_MIN, 0, c->gsComm.c);
  }
  return 0;
}

int GenmapBcast(GenmapComm c, void *in, GenmapInt count, GenmapDataType type) {
  return MPI_Bcast(in, count, type, 0, c->gsComm.c);
}

void GenmapSplitComm(GenmapHandle h, GenmapComm *c, int bin) {
  GenmapCommExternal local;
  int id = GenmapCommRank(*c);
  MPI_Comm_split((*c)->gsComm.c, bin, id, &local);
  GenmapCrystalFinalize(h);
  GenmapDestroyComm(*c);
  GenmapCreateComm(c, local);
  MPI_Comm_free(&local);
  GenmapCrystalInit(h, *c);
}

int GenmapCrystalInit(GenmapHandle h, GenmapComm c) {
  crystal_init(&(h->cr), &(c->gsComm));
  return 0;
}

int GenmapCrystalTransfer(GenmapHandle h, int field) {
  if(field == GENMAP_ORIGIN)
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), origin, 0,
                    &(h->cr));
  else if(field == GENMAP_PROC)
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0,
                    &(h->cr));
  return 0;
}

int GenmapCrystalFinalize(GenmapHandle h) {
  crystal_free(&(h->cr));
  return 0;
}
