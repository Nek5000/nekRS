#include <genmap-impl.h>

#define min(a,b) ((b)<(a)?(b):(a))

typedef struct{
  GenmapULong sequenceId;
  GenmapLong neighbors[8];
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

typedef struct{
  GenmapULong elementId;
} element;

struct array *GenmapFindNeighbors(GenmapHandle h,GenmapComm c)
{
  struct comm cc=c->gsc;

  sint lelt=GenmapGetNLocalElements(h);
  sint nv  =GenmapGetNVertices(h);

  GenmapScan(h,c);
  ulong elem_id   =GenmapGetLocalStartIndex(h)+1;
  ulong sequenceId=elem_id*nv;

  size_t size=lelt*nv;
  struct array vertices; array_init(vertex,&vertices,size);

  GenmapElements elems=GenmapGetElements(h);
  sint i,j;
  for(i=0;i<lelt;i++){
    for(j=0;j<nv;j++){
      vertex t={
        .elementId =elem_id,
        .sequenceId=sequenceId,
        .vertexId  =elems[i].vertices[j],
        .workProc  =elems[i].vertices[j]%cc.np
      };
      array_cat(vertex,&vertices,&t,1);
      sequenceId++;
    }
    elem_id++;
  }

  struct crystal cr; crystal_init(&cr,&cc);
  sarray_transfer(vertex,&vertices,workProc,1,&cr);
  size=vertices.n; vertex* vPtr=vertices.ptr;

  buffer buf; buffer_init(&buf,1024);
  sarray_sort(vertex,vPtr,size,vertexId,1,&buf);

  struct array a; array_init(csr_entry,&a,10);

  //FIXME: Assumes quads or hexes
  sint s=0,e; csr_entry t;
  while(s<size){
    for(e=s+1; e<size && vPtr[s].vertexId==vPtr[e].vertexId; e++);
    int nNeighbors=min(e,size)-s;
    for(i=s;i<min(e,size);i++){
      t.r=vPtr[i].elementId; t.proc=vPtr[i].workProc;
      for(j=0;j<nNeighbors;j++){
        t.c=vPtr[s+j].elementId;
        array_cat(csr_entry,&a,&t,1);
      }
    }
    s=e;
  }

  sarray_transfer(csr_entry,&a,proc,1,&cr);
  sarray_sort_2(csr_entry,a.ptr,a.n,r,1,c,1,&buf);

  struct array *nbrs=tmalloc(struct array,1);
  array_init(entry,nbrs,lelt);

  if(lelt==0){
    crystal_free(&cr);
    buffer_free(&buf);
    array_free(&vertices);
    array_free(&a);
    return nbrs;
  }

  csr_entry *aptr=a.ptr; entry *nptr=nbrs->ptr;
  entry ee,ep; ep.r=aptr->r; ep.c=aptr->c; array_cat(entry,nbrs,&ep,1);
  for(i=1; i<a.n; i++){
    ee.r=aptr[i].r,ee.c=aptr[i].c; ulong n=nbrs->n-1;
    if(ee.r!=ep.r || ee.c!=ep.c){
      array_cat(entry,nbrs,&ee,1);
      ep=ee;
    }
  }

  sarray_sort_2(entry,nbrs->ptr,nbrs->n,r,1,c,1,&buf);

  crystal_free(&cr);
  buffer_free(&buf);
  array_free(&vertices);
  array_free(&a);

  return nbrs;
}

int GenmapInitLaplacian(GenmapHandle h,GenmapComm c)
{
  struct array *entries=GenmapFindNeighbors(h,c);
  csr_mat_setup(entries,&c->gsc,&c->M);
  array_free(entries); free(entries);

  c->gsh=get_csr_top(c->M,&c->gsc);
  GenmapMalloc(c->M->row_off[c->M->rn],&c->b);

  return 0;
}

int GenmapLaplacian(GenmapHandle h,GenmapComm c,GenmapVector u,
    GenmapVector v)
{
  assert(u->size==v->size);
  csr_mat_gather(c->M,c->gsh,u->data,c->b,&c->buf);
  csr_mat_apply(v->data,c->M,c->b);

  return 0;
}

#undef min
