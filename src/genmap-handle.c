#include "genmap-impl.h"

GenmapElements GenmapGetElements(genmap_handle h) {
  return (GenmapElements) h->elements->ptr;
}

GenmapInt GenmapGetNLocalElements(genmap_handle h) {
  return h->elements->n;
}

void GenmapSetArrayElements(genmap_handle h, struct array *localElements) {
  h->elements = localElements;
}

GenmapLong GenmapGetNGlobalElements(genmap_handle h) {
  return h->nel;
}

void GenmapSetNGlobalElements(genmap_handle h, GenmapLong globalElements) {
  h->nel = globalElements;
}

GenmapLong GenmapGetLocalStartIndex(genmap_handle h) {
  return h->start;
}

void GenmapSetLocalStartIndex(genmap_handle h, GenmapLong localStart) {
  h->start = localStart;
}

int GenmapGetNVertices(genmap_handle h) {
  return h->nv;
}

void GenmapSetNVertices(genmap_handle h, int nVertices) {
  h->nv = nVertices;
}

void GenmapScan(genmap_handle h, GenmapComm c) {
  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt = GenmapGetNLocalElements(h);
  comm_scan(out,&(c->gsc),gs_long_long,gs_add,&lelt,1,buf);
  GenmapSetLocalStartIndex(h, out[0][0]);
  GenmapSetNGlobalElements(h, out[1][0]);
}
