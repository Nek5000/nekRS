#include <stdio.h>
#include <stdlib.h>

#include <genmap-impl.h>

int genmap_init(genmap_handle *h_, comm_ext ce, parRSB_options *options) {
  GenmapMalloc(1, h_);
  genmap_handle h = *h_;

  GenmapMalloc(1, &h->global);
  comm_init(h->global, ce);

  GenmapMalloc(1, &h->local);
  comm_init(h->local, ce);

  /* Weighted Laplacian */
  h->weights = NULL;
  h->gsw = NULL;
  buffer_init(&h->buf, 1024);

  /* Un-weighted Laplacian */
  h->gsh = NULL;
  h->M = NULL;
  h->b = NULL;

  h->options = options;

  return 0;
}

int genmap_finalize(genmap_handle h) {
  /* Weighted Laplacian */
  if (h->weights != NULL)
    GenmapFree(h->weights);
  if (h->gsw != NULL)
    gs_free(h->gsw);
  buffer_free(&h->buf);

  /* Un-weighted Laplacian */
  if (h->gsh)
    gs_free(h->gsh);
  if (h->M)
    csr_mat_free(h->M);
  if (h->b != NULL)
    GenmapFree(h->b);

  if (h->global != NULL) {
    comm_free(h->global);
    GenmapFree(h->global);
  }

  if (h->local != NULL) {
    comm_free(h->local);
    GenmapFree(h->local);
  }

  GenmapFree(h);

  return 0;
}

void *genmap_get_elements(genmap_handle h) {
  return (struct rsb_element *)h->elements->ptr;
}
void genmap_set_elements(genmap_handle h, struct array *elements) {
  h->elements = elements;
}

int genmap_get_nvertices(genmap_handle h) { return h->nv; }
void genmap_set_nvertices(genmap_handle h, int nv) { h->nv = nv; }

GenmapInt genmap_get_nel(genmap_handle h) { return h->elements->n; }

GenmapULong genmap_get_partition_nel(genmap_handle h) { return h->nel; }
void genmap_set_partition_nel(genmap_handle h, GenmapULong globalElements) {
  h->nel = globalElements;
}

GenmapLong genmap_get_local_start_index(genmap_handle h) { return h->start; }
void genmap_set_local_start_index(genmap_handle h, GenmapLong localStart) {
  h->start = localStart;
}

int GenmapMallocArray(size_t n, size_t unit, void *p) {
  int ierr = posix_memalign((void **)p, GENMAP_ALIGN, n * unit);
  if (ierr)
    printf("GenmapMallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  return ierr;
}

int GenmapCallocArray(size_t n, size_t unit, void *p) {
  int ierr = 0;
  *(void **)p = calloc(n, unit);
  if (n && unit && !*(void **)p) {
    ierr = 1;
    printf("GenmapCallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  }
  return ierr;
}

int GenmapReallocArray(size_t n, size_t unit, void *p) {
  int ierr = 0;
  *(void **)p = realloc(*(void **)p, n * unit);
  if (n && unit && !*(void **)p) {
    ierr = 1;
    printf("GenmapReallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  }
  return ierr;
}

int GenmapFree(void *p) {
  free(p);
  p = NULL;
  return 0;
}
