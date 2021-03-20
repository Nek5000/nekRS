#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <genmap-impl.h>
#include <parRSB.h>

struct unmarked {
  uint index;
};

/* Find the number of disconnected components */
sint get_components(sint *component, struct rsb_element *elements,
                    struct comm *c, buffer *buf, uint nelt, uint nv) {
  slong nelt_ = nelt;
  slong out[2][1], buff[2][1];
  comm_scan(out, c, gs_long, gs_add, &nelt_, 1, buff);
  slong nelg = out[1][0];
  ulong start = out[0][0];

  if (nelg == 0)
    return 0;

  GenmapLong *p, *ids;
  GenmapMalloc(nelt * nv, &p);
  GenmapMalloc(nelt * nv, &ids);

  int null_input = 0;
  if (component == NULL) {
    GenmapMalloc(nelt, &component);
    null_input = 1;
  }

  uint e;
  for (e = 0; e < nelt; e++)
    component[e] = -1;

  struct array arr;
  array_init(struct unmarked, &arr, nelt);

  struct comm cc;

  struct unmarked u;
  uint d;
  slong nnz1, nnzg, nnzg0, nnzb, nmarked = 0;
  sint count = 0;

  do {
    // Count unmarked elements
    arr.n = 0;
    for (e = 0; e < nelt; e++) {
      if (component[e] == -1) {
        u.index = e;
        array_cat(struct unmarked, &arr, &u, 1);
      }
    }

    int bin = 0;
    if (arr.n > 0)
      bin = 1;

    genmap_comm_split(c, bin, c->id, &cc);

    nnz1 = nnzg = nnzg0 = 0;

    if (bin == 1) {
      // Initialize p
      for (e = 0; e < arr.n; e++)
        for (d = 0; d < nv; d++)
          p[e * nv + d] = 0;

      // Mark the first non-marked element as seed
      struct unmarked *ptr = (struct unmarked *)arr.ptr;
      slong first = start + ptr[0].index;
      slong first_ = first;
      comm_allreduce(&cc, gs_long, gs_min, &first_, 1, buff);

      if (first_ == first) {
        for (d = 0; d < nv; d++)
          p[0 * nv + d] = 1;
      }

      // Setup gs
      for (e = 0; e < arr.n; e++)
        for (d = 0; d < nv; d++)
          ids[e * nv + d] = elements[ptr[e].index].vertices[d];

      struct gs_data *gsh = gs_setup(ids, arr.n * nv, &cc, 0, gs_pairwise, 0);

      do {
        gs(p, gs_long, gs_add, 0, gsh, buf);

        nnz1 = 0;
        for (e = 0; e < arr.n; e++) {
          for (d = 0; d < nv; d++) {
            if (p[e * nv + d] > 0) {
              nnz1++;
              component[ptr[e].index] = count;
              break;
            }
          }

          if (d < nv) { // There was one non-zero vertex in the element
            for (d = 0; d < nv; d++)
              p[e * nv + d] = 1;
          }
        }

        nnzg0 = nnzg, nnzg = nnz1;
        comm_allreduce(&cc, gs_long, gs_add, &nnzg, 1, &nnzb);
      } while (nnzg > nnzg0);

      gs_free(gsh);
    }

    comm_free(&cc);

    comm_allreduce(c, gs_long, gs_add, &nnz1, 1, &nnzb);
    nmarked += nnz1;

    count++;
  } while (nmarked < nelg);

  if (null_input == 1)
    GenmapFree(component);

  GenmapFree(p);
  GenmapFree(ids);

  return count;
}
