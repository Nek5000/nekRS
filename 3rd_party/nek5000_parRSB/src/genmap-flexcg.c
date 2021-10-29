#include <genmap-impl.h>

int flex_cg(genmap_handle h, struct comm *gsc, mgData d, genmap_vector ri,
            int maxIter, genmap_vector x) {
  assert(x->size == ri->size);
  assert(x->size == genmap_get_nel(h));

  uint lelt = x->size;
  GenmapLong nelg = genmap_get_partition_nel(h);

  genmap_vector z0, z, dz, w, p, r;
  genmap_vector_create(&z, lelt);
  genmap_vector_create(&w, lelt);
  genmap_vector_create(&r, lelt);
  genmap_vector_create(&p, lelt);
  genmap_vector_create(&z0, lelt);
  genmap_vector_create(&dz, lelt);

#define PREC 1
#define ORTH 1

  int rank = gsc->id;

  uint i;
  for (i = 0; i < lelt; i++)
    x->data[i] = 0.0, r->data[i] = ri->data[i];

  genmap_vector_copy(z, r);
#if ORTH
  genmap_vector_ortho_one(gsc, z, nelg);
#endif

  GenmapScalar den, alpha, beta, rz0, rz1, rz2, rr, buf;

  rz1 = genmap_vector_dot(r, z);
  comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, &buf);

  genmap_vector_copy(p, z);

  i = 0;
  while (i < maxIter && sqrt(rz1) > GENMAP_TOL) {
    GenmapLaplacian(h, p->data, w->data);

    den = genmap_vector_dot(p, w);
    comm_allreduce(gsc, gs_double, gs_add, &den, 1, &buf);

    alpha = rz1 / den;

    genmap_vector_axpby(x, x, 1.0, p, alpha);
    genmap_vector_axpby(r, r, 1.0, w, -alpha);

    genmap_vector_copy(z0, z);
#if PREC
    mg_vcycle(z->data, r->data, d);
#else
    genmap_vector_copy(z, r);
#endif
#if ORTH
    genmap_vector_ortho_one(gsc, z, nelg);
#endif

    rz0 = rz1;

    rz1 = genmap_vector_dot(r, z);
    comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, &buf);

    genmap_vector_axpby(dz, z, 1.0, z0, -1.0);
    rz2 = genmap_vector_dot(r, dz);
    comm_allreduce(gsc, gs_double, gs_add, &rz2, 1, &buf);
    beta = rz2 / rz0;

    genmap_vector_axpby(p, z, 1.0, p, beta);
    i++;

    rr = genmap_vector_dot(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rr, 1, &buf);
  }

  GenmapDestroyVector(z);
  GenmapDestroyVector(w);
  GenmapDestroyVector(p);
  GenmapDestroyVector(r);
  GenmapDestroyVector(z0);
  GenmapDestroyVector(dz);

  return i;
}
