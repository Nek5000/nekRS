#include <genmap-impl.h>

#define min(a, b) ((b) < (a) ? (b) : (a))

static void genmap_find_neighbors(struct array *nbrs, genmap_handle h,
                                  struct comm *cc) {
  sint lelt = genmap_get_nel(h);
  sint nv = genmap_get_nvertices(h);

  genmap_comm_scan(h, cc);
  ulong elem_id = genmap_get_local_start_index(h) + 1;
  ulong sequenceId = elem_id * nv;

  size_t size = lelt * nv;
  struct array vertices;
  array_init(vertex, &vertices, size);

  struct rsb_element *elems = genmap_get_elements(h);
  sint i, j;
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < nv; j++) {
      vertex vrt = {.sequenceId = sequenceId,
                    .nNeighbors = 0,
                    .elementId = elem_id,
                    .vertexId = elems[i].vertices[j],
                    .workProc = elems[i].vertices[j] % cc->np};
      array_cat(vertex, &vertices, &vrt, 1);
      sequenceId++;
    }
    elem_id++;
  }
  assert(vertices.n == lelt * nv);

  struct crystal cr;
  crystal_init(&cr, cc);

  sarray_transfer(vertex, &vertices, workProc, 1, &cr);
  size = vertices.n;
  vertex *vPtr = vertices.ptr;

  sarray_sort(vertex, vPtr, size, vertexId, 1, &h->buf);

  struct array a;
  array_init(csr_entry, &a, 10);

  // FIXME: Assumes quads or hexes
  sint s = 0, e;
  csr_entry t;
  while (s < size) {
    e = s + 1;
    while (e < size && vPtr[s].vertexId == vPtr[e].vertexId)
      e++;
    int n_neighbors = min(e, size) - s;

    for (i = s; i < min(e, size); i++) {
      t.r = vPtr[i].elementId;
      t.proc = vPtr[i].workProc;
      for (j = 0; j < n_neighbors; j++) {
        t.c = vPtr[s + j].elementId;
        array_cat(csr_entry, &a, &t, 1);
      }
    }
    s = e;
  }

  sarray_transfer(csr_entry, &a, proc, 1, &cr);
  // TODO: Check if the last line is redundant
  sarray_sort_2(csr_entry, a.ptr, a.n, r, 1, c, 1, &h->buf);
  // sarray_sort(csr_entry, a.ptr, a.n, r, 1, &h->buf);

  array_init(entry, nbrs, lelt);

  if (a.n == 0) {
    crystal_free(&cr);
    array_free(&vertices);
    array_free(&a);
  }

  csr_entry *aptr = a.ptr;
  entry *nptr = nbrs->ptr;

  entry ee = {0, 0, 0, 0, 0, 0.0}, ep = {0, 0, 0, 0, 0.0};
  ep.r = aptr[0].r;
  ep.c = aptr[0].c;
  array_cat(entry, nbrs, &ep, 1);

  for (i = 1; i < a.n; i++) {
    ee.r = aptr[i].r;
    ee.c = aptr[i].c;
    if (ee.r != ep.r || ee.c != ep.c) {
      array_cat(entry, nbrs, &ee, 1);
      ep = ee;
    }
  }

  sarray_sort_2(entry, nbrs->ptr, nbrs->n, r, 1, c, 1, &h->buf);

  crystal_free(&cr);
  array_free(&vertices);
  array_free(&a);
}

int GenmapInitLaplacian(genmap_handle h, struct comm *c) {
  struct array entries;

  metric_tic(c, FINDNBRS);
  genmap_find_neighbors(&entries, h, c);
  metric_toc(c, FINDNBRS);

  metric_tic(c, CSRMATSETUP);
  csr_mat_setup(&entries, c, &h->M);
  metric_toc(c, CSRMATSETUP);

  array_free(&entries);

  metric_toc(c, CSRTOPSETUP);
  h->gsh = get_csr_top(h->M, c);
  metric_toc(c, CSRTOPSETUP);

  GenmapRealloc(h->M->row_off[h->M->rn], &h->b);

#if defined(GENMAP_DEBUG)
  int nnz = h->M->row_off[h->M->rn];
  double fro[2] = {0.0, 0.0}, buf[2];
  for (int i = 0; i < nnz; i++) {
    fro[0] += h->M->v[i];
    fro[1] += h->M->v[i] * h->M->v[i];
  }
  comm_allreduce(c, gs_double, gs_add, &fro, 2, &buf);
  if (c->gsc.id == 0)
    printf("nrom(G,'1')=%g\nnorm(G,'fro')=%g\n", fro[0], fro[1]);
#endif

  return 0;
}

int GenmapLaplacian(genmap_handle h, GenmapScalar *u, GenmapScalar *v) {
  csr_mat_gather(h->M, h->gsh, u, h->b, &h->buf);
  csr_mat_apply(v, h->M, h->b);

  return 0;
}

#undef min
