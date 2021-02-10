#include <genmap-impl.h>
#include <gencon-impl.h>

int mesh_init(Mesh *m_,int nel,int nDim){
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

int mesh_free(Mesh m){
  array_free(&m->elements);
  array_free(&m->boundary);

  free(m);

  return 0;
}

void get_vertex_ids(long long **vertex_ids_,Mesh mesh){
  int nelt=mesh->nelt;
  int nv=(mesh->nDim==3)?8:4;

  GenmapMalloc(nelt*nv,vertex_ids_);
  long long *vertex_ids=*vertex_ids_;

  Point ptr=mesh->elements.ptr;
  int e,v,count=0;
  for(e=0; e<nelt; e++){
    for(v=0; v<nv; v++){
      vertex_ids[count]=ptr[count].globalId;
      count++;
    }
  }
}

void get_vertex_coordinates(double **coords_,Mesh mesh){
  int nelt=mesh->nelt;
  int ndim=mesh->nDim;
  int nv=(ndim==3)?8:4;

  GenmapMalloc(nelt*nv*ndim,coords_);
  double *coords=*coords_;

  Point ptr=mesh->elements.ptr;
  int e,v,d,count=0;
  for(e=0; e<nelt; e++){
    for(v=0; v<nv; v++){
      for(d=0; d<ndim; d++){
        coords[count++]=ptr[e*nv+v].x[d];
      }
    }
  }
}

int get_mesh_dim(Mesh mesh){
  return mesh->nDim;
}

int get_mesh_nel(Mesh mesh){
  return mesh->nelt;
}
