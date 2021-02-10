#include <genmap-impl.h>

/* Load balance input data */
void genmap_load_balance(struct array *eList,double *coord,long long *vtx,
    uint nel,int nv,struct crystal *cr,buffer *bfr) {
  slong in = nel;
  slong out[2][1], buf[2][1];
  comm_scan(out, &cr->comm, gs_long, gs_add, &in, 1, buf);
  slong nelg_start = out[0][0];
  slong nelg = out[1][0];

  int size = cr->comm.np;
  uint nstar = nelg/size;
  if(nstar == 0)
    nstar = 1;

  size_t unit_size;
  struct rcb_element *element = NULL;

  if(vtx == NULL) { // RCB
    unit_size = sizeof(struct rcb_element);
    element = calloc(1, sizeof(struct rcb_element));
    element->type=GENMAP_RCB_ELEMENT;
  } else {
    unit_size = sizeof(struct rsb_element);
    element = calloc(1, sizeof(struct rsb_element));
    element->type=GENMAP_RSB_ELEMENT;
  }

  element->origin = cr->comm.id;

  array_init_(eList, nel, unit_size, __FILE__, __LINE__);

  int ndim=(nv==8)?3:2;
  int e,n,v;
  for(e=0; e<nel; ++e){
    slong eg = element->globalId = nelg_start + e + 1;
    element->proc = (int)((eg-1)/nstar);
    if(eg > size*nstar)
      element->proc = (eg%size)-1;

    element->coord[0] = element->coord[1] = element->coord[2] = 0.0;
    for(v=0;v<nv;v++){
      for(n=0;n<ndim;n++)
        element->coord[n] += coord[e*ndim*nv + v*ndim + n];
    }
    for(n=0;n<ndim;n++)
      element->coord[n] /= nv;

    array_cat_(unit_size, eList, element, 1, __FILE__, __LINE__);
  }
  assert(eList->n==nel);

  if(vtx != NULL) {
    struct rsb_element *elements = eList->ptr;
    for(e=0; e<nel; ++e) {
      for(v=0;v<nv;v++)
        elements[e].vertices[v]=vtx[e*nv+v];
    }
  }

  sarray_transfer_(eList, unit_size, offsetof(struct rcb_element, proc), 1, cr);
  nel=eList->n;
  if(vtx != NULL)
    sarray_sort(struct rsb_element, eList->ptr, nel, globalId, 1, bfr);
  else
    sarray_sort(struct rcb_element, eList->ptr, nel, globalId, 1, bfr);

  free(element);
}
