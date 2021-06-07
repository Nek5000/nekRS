#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>

int GenmapLanczosLegendary(genmap_handle h, struct comm *gsc, genmap_vector f,
                           GenmapInt niter, genmap_vector **rr,
                           genmap_vector diag, genmap_vector upper) {
  assert(diag->size == niter);
  assert(diag->size == upper->size + 1);
  assert(f->size == genmap_get_nel(h));

  if (genmap_get_partition_nel(h) < niter) {
    niter = genmap_get_partition_nel(h);
    diag->size = niter;
    upper->size = niter - 1;
  }

  GenmapScalar eps = 1.e-5;
  GenmapScalar alpha, beta;
  GenmapScalar rnorm, rtol, rni, rtr, rtz1, rtz2, pap, pap_old;
  genmap_vector r, p, w;

  rtz1 = 1.0;
  pap = 0.0;
  GenmapInt lelt = genmap_get_nel(h);

  genmap_vector_create_zeros(&p, lelt);
  genmap_vector_create(&w, lelt);
  genmap_vector_create(&r, lelt);

  genmap_vector_copy(r, f);

  genmap_vector_ortho_one(gsc, r, genmap_get_partition_nel(h));
  rtr = genmap_vector_dot(r, r);
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rni);
  rnorm = sqrt(rtr);
  rtol = rnorm * eps;
  rni = 1.0 / rnorm;

  if (*rr == NULL) {
    GenmapMalloc((size_t)(niter + 1), rr);
    GenmapInt i;
    for (i = 0; i < niter + 1; ++i)
      (*rr)[i] = NULL;
  }
  genmap_vector_create(&(*rr)[0], lelt);

  genmap_vector_scale((*rr)[0], r, rni);

  int iter;
  for (iter = 0; iter < niter; iter++) {
    rtz2 = rtz1;
    rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0)
      beta = 0.0;

    genmap_vector_axpby(p, p, beta, r, 1.0);
    genmap_vector_ortho_one(gsc, p, genmap_get_partition_nel(h));

    metric_tic(gsc, WEIGHTEDLAPLACIAN);
    GenmapLaplacianWeighted(h, p->data, w->data);
    metric_tic(gsc, WEIGHTEDLAPLACIAN);

    genmap_vector_scale(w, w, -1.0);

    pap_old = pap;
    pap = genmap_vector_dot(w, p);
    comm_allreduce(gsc, gs_double, gs_add, &pap, 1, &rni);

    alpha = rtz1 / pap;
    genmap_vector_axpby(r, r, 1.0, w, -1.0 * alpha);

    rtr = genmap_vector_dot(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rni);
    rnorm = sqrt(rtr);
    rni = 1.0 / rnorm;

    genmap_vector_create(&(*rr)[iter + 1], lelt);
    genmap_vector_scale((*rr)[iter + 1], r, rni);

    if (iter == 0) {
      diag->data[iter] = pap / rtz1;
    } else {
      diag->data[iter] = (beta * beta * pap_old + pap) / rtz1;
      upper->data[iter - 1] = -beta * pap_old / sqrt(rtz2 * rtz1);
    }

    if (rnorm < rtol) {
      diag->size = iter + 1;
      upper->size = iter;
      iter = iter + 1;
      break;
    }
  }

  metric_acc(LANCZOSTOLFINAL, rnorm);
  metric_acc(LANCZOSTOLTARGET, rtol);

  GenmapDestroyVector(p);
  GenmapDestroyVector(w);
  GenmapDestroyVector(r);

  return iter;
}

int GenmapLanczos(genmap_handle h, struct comm *gsc, genmap_vector init,
                  GenmapInt iter, genmap_vector **q, genmap_vector alpha,
                  genmap_vector beta) {
  assert(alpha->size == iter);
  assert(alpha->size == beta->size + 1);
  assert(init->size == genmap_get_nel(h));

  if (genmap_get_partition_nel(h) < iter) {
    iter = genmap_get_partition_nel(h);
    alpha->size = iter;
    beta->size = iter - 1;
  }

  GenmapInt lelt = genmap_get_nel(h);

  genmap_vector q1;
  genmap_vector_create(&q1, lelt);
  genmap_vector_copy(q1, init);

  genmap_vector_ortho_one(gsc, q1, genmap_get_partition_nel(h));
  GenmapScalar normq1 = genmap_vector_dot(q1, q1), buf;
  comm_allreduce(gsc, gs_double, gs_add, &normq1, 1, &buf);
  normq1 = sqrt(normq1);
  genmap_vector_scale(q1, q1, 1. / normq1);

  /* Set q_0 and beta_0 to zero (both uses 0-indexing) */
  genmap_vector q0;
  genmap_vector_create_zeros(&q0, lelt);
  beta->data[0] = 0.;

  genmap_vector u;
  genmap_vector_create(&u, lelt);

  if (*q == NULL) {
    GenmapMalloc((size_t)iter, q);
    GenmapInt i;
    for (i = 0; i < iter; ++i)
      (*q)[i] = NULL;
  }

  GenmapScalar b = 0.0;

  int k;
  for (k = 0; k < iter; k++) {
    genmap_vector_create(&(*q)[k], lelt);
    genmap_vector_copy((*q)[k], q1);

    metric_tic(gsc, WEIGHTEDLAPLACIAN);
    GenmapLaplacianWeighted(h, q1->data, u->data);
    metric_toc(gsc, WEIGHTEDLAPLACIAN);

    alpha->data[k] = genmap_vector_dot(q1, u);
    comm_allreduce(gsc, gs_double, gs_add, &alpha->data[k], 1, &buf);

    genmap_vector_axpby(u, u, 1., q0, -b);
    genmap_vector_axpby(u, u, 1., q1, -alpha->data[k]);

    b = genmap_vector_dot(u, u);
    comm_allreduce(gsc, gs_double, gs_add, &b, 1, &buf);
    b = sqrt(b);

    if (k < iter - 1) {
      beta->data[k] = b;
      if (fabs(b) < normq1 * GENMAP_TOL) {
        beta->size = k;
        alpha->size = k + 1;
        k = k + 1;
        break;
      }

      genmap_vector_copy(q0, q1);
      genmap_vector_scale(q1, u, 1. / beta->data[k]);
    }
  }

  GenmapDestroyVector(q0);
  GenmapDestroyVector(q1);
  GenmapDestroyVector(u);

  return k;
}
