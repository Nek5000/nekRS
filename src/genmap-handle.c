#include "genmap-impl.h"
//
// GenmapHandle
//
int GenmapCreateHandle(GenmapHandle h) {
  // Datastructures
  h->global = NULL;
  h->local  = NULL;

  h->elementArray.ptr = NULL;
  h->elementArray.n = h->elementArray.max = 0;

  h->Create(h);

  return 0;
}

int GenmapDestroyHandle(GenmapHandle h) {
  GenmapFree(h);
  return 0;
}

GenmapElements GenmapGetElements(GenmapHandle h) {
  return (GenmapElements) h->elementArray.ptr;
}

void GenmapSetElements(GenmapHandle h, GenmapElements elements) {
  h->elementArray.ptr = elements;
}

GenmapInt GenmapGetNLocalElements(GenmapHandle h) {
  return h->elementArray.n;
}

void GenmapSetNLocalElements(GenmapHandle h, GenmapInt localElements) {
  array_init(struct GenmapElement_private, &h->elementArray, localElements);
  h->elementArray.n = localElements;
}

GenmapLong GenmapGetNGlobalElements(GenmapHandle h) {
  return h->nel;
}

void GenmapSetNGlobalElements(GenmapHandle h, GenmapLong globalElements) {
  h->nel = globalElements;
}

GenmapLong GenmapGetLocalStartIndex(GenmapHandle h) {
  return h->start;
}

void GenmapSetLocalStartIndex(GenmapHandle h, GenmapLong localStart) {
  h->start = localStart;
}

int GenmapGetNVertices(GenmapHandle h) {
  return h->nv;
}

void GenmapSetNVertices(GenmapHandle h, int nVertices) {
  h->nv = nVertices;
}
