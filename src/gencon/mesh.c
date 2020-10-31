#include <genmap-impl.h>
#include <gencon-impl.h>

int MeshInit(Mesh *m_,int nel,int nDim){
  GenmapMalloc(1,m_);
  Mesh m=*m_;

  m->nelt=nel;
  m->nDim=nDim;
  m->nNeighbors=nDim;
  m->nVertex=(nDim==2)?4:8;

  array_init(struct Point_private   ,&m->elements,10); m->elements.n=0;
  array_init(struct Boundary_private,&m->boundary,10); m->boundary.n=0;

  return 0;
}

Element MeshGetElements(Mesh m){
  return (Element) m->elements.ptr;
}

int MeshFree(Mesh m){
  array_free(&m->elements);
  array_free(&m->boundary);

  free(m);

  return 0;
}
