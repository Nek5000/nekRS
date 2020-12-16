#include "genmap-impl.h"

int GenmapCreateComm(GenmapComm *c_,comm_ext ce){
  GenmapMalloc(1,c_);
  GenmapComm c=*c_;
  comm_init(&c->gsc, ce);

  c->gsh=NULL;
  c->M=NULL;

  c->gsw=NULL;

  buffer_init(&c->buf,1024);

  c->b=NULL;

  return 0;
}

int GenmapDestroyComm(GenmapComm c) {
  buffer_free(&c->buf);

  if(c->gsh)
    gs_free(c->gsh);
  if(c->M)
    csr_mat_free(c->M);

  if(c->gsw)
    gs_free(c->gsw);

  if(c->b)
    GenmapFree(c->b);

  comm_free(&c->gsc);
  GenmapFree(c);

  return 0;
}

int GenmapCommSize(GenmapComm c) {
  return (int) c->gsc.np;
}

int GenmapCommRank(GenmapComm c) {
  return (int) c->gsc.id;
}

GenmapComm GenmapGetLocalComm(genmap_handle h) {
  return h->local;
}

void GenmapSetLocalComm(genmap_handle h, GenmapComm c) {
  h->local = c;
}

GenmapComm GenmapGetGlobalComm(genmap_handle h) {
  return h->global;
}

void GenmapSetGlobalComm(genmap_handle h, GenmapComm c) {
  h->global = c;
}

int GenmapGop(GenmapComm c, void *v, GenmapInt size,
              GenmapDataType type, GenmapInt op) {
  if(op == GENMAP_SUM) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_SUM, c->gsc.c);
  } else if(op == GENMAP_MAX) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_MAX, c->gsc.c);
  } else if(op == GENMAP_MIN) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_MIN, c->gsc.c);
  }
  return 0;
}

int GenmapReduce(GenmapComm c, void *out, void *in, GenmapInt size,
                 GenmapDataType type, GenmapInt op) {
  if(op == GENMAP_SUM) {
    MPI_Reduce(in, out, size, type, MPI_SUM, 0, c->gsc.c);
  } else if(op == GENMAP_MAX) {
    MPI_Reduce(in, out, size, type, MPI_MAX, 0, c->gsc.c);
  } else if(op == GENMAP_MIN) {
    MPI_Reduce(in, out, size, type, MPI_MIN, 0, c->gsc.c);
  }
  return 0;
}

int GenmapBcast(GenmapComm c, void *in, GenmapInt count, GenmapDataType type) {
  return MPI_Bcast(in, count, type, 0, c->gsc.c);
}

void GenmapSplitComm(genmap_handle h, GenmapComm *c, int bin) {
  comm_ext local;
  int id = GenmapCommRank(*c);
  MPI_Comm_split((*c)->gsc.c, bin, id, &local);
  GenmapCrystalFinalize(h);
  GenmapDestroyComm(*c);
  GenmapCreateComm(c, local);
  MPI_Comm_free(&local);
  GenmapCrystalInit(h, *c);
}

int GenmapCrystalInit(genmap_handle h, GenmapComm c) {
  crystal_init(&(h->cr), &(c->gsc));
  return 0;
}

int GenmapCrystalTransfer(genmap_handle h, int field) {
  if(field == GENMAP_ORIGIN)
    sarray_transfer(struct rsb_element,h->elements,origin,0,&h->cr);
  else if(field == GENMAP_PROC)
    sarray_transfer(struct rsb_element,h->elements,proc  ,0,&h->cr);
  return 0;
}

int GenmapCrystalFinalize(genmap_handle h) {
  crystal_free(&(h->cr));
  return 0;
}
