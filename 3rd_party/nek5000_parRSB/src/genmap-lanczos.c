#include "genmap-impl.h"

int GenmapLanczos(genmap_handle h, struct comm *gsc, genmap_vector f, int niter,
                  genmap_vector **rr_, genmap_vector diag,
                  genmap_vector upper) {
  assert(diag->size == niter);
  assert(diag->size == upper->size + 1);
  assert(f->size == genmap_get_nel(h));

  GenmapInt lelt = genmap_get_nel(h);
  slong in = lelt;
  slong out[2][1], buf[2][1];
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, buf);
  slong start = out[0][0];
  slong nelg = out[1][0];

  if (nelg < niter) {
    niter = nelg;
    diag->size = niter;
    upper->size = niter - 1;
  }

  GenmapMalloc(niter + 1, rr_);
  genmap_vector *rr = *rr_;
  GenmapInt i;
  for (i = 0; i < niter + 1; ++i)
    rr[i] = NULL;

  genmap_vector r, p, w;
  genmap_vector_create_zeros(&p, lelt);
  genmap_vector_create(&w, lelt);
  genmap_vector_create(&r, lelt);

  genmap_vector_copy(r, f);

  genmap_vector_ortho_one(gsc, r, nelg);
  GenmapScalar rtr = genmap_vector_dot(r, r);
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
  GenmapScalar rnorm = sqrt(rtr);
  GenmapScalar rni = 1.0 / rnorm;

  genmap_vector_create(&rr[0], lelt);
  genmap_vector_scale(rr[0], r, rni);

  GenmapScalar eps = 1.e-5;
  GenmapScalar rtol = rnorm * eps;

  GenmapScalar rtz1 = 1.0;
  GenmapScalar pap = 0.0;
  GenmapScalar alpha, beta;
  GenmapScalar rtz2, pap_old;

  int iter;
  for (iter = 0; iter < niter; iter++) {
    rtz2 = rtz1;
    rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0)
      beta = 0.0;

    /* add2s1(p,r,beta,n) */
    for (i = 0; i < lelt; i++)
      p->data[i] = beta * p->data[i] + r->data[i];
    GenmapScalar pp = genmap_vector_dot(p, p);
    comm_allreduce(gsc, gs_double, gs_add, &pp, 1, &rni);

    genmap_vector_ortho_one(gsc, p, nelg);

    metric_tic(gsc, LAPLACIAN);
    GenmapLaplacianWeighted(h, p->data, w->data);
    metric_tic(gsc, LAPLACIAN);

    GenmapScalar ww = genmap_vector_dot(w, w);
    comm_allreduce(gsc, gs_double, gs_add, &ww, 1, &rni);

    pap_old = pap;
    pap = genmap_vector_dot(w, p);
    comm_allreduce(gsc, gs_double, gs_add, &pap, 1, &rni);
//    if (gsc->id == 0)
//      printf("iter = %d beta = %lf pp = %lf pap = %lf ww = %lf\n", iter, beta,
//             pp, pap, ww);

    alpha = rtz1 / pap;
    genmap_vector_axpby(r, r, 1.0, w, -1.0 * alpha);

    rtr = genmap_vector_dot(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rni);
    rnorm = sqrt(rtr);
    rni = 1.0 / rnorm;

    genmap_vector_create(&rr[iter + 1], lelt);
    genmap_vector_scale(rr[iter + 1], r, rni);

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

  metric_acc(LANCZOS_TOL_FINAL, rnorm);
  metric_acc(LANCZOS_TOL_TARGET, rtol);

  GenmapDestroyVector(p);
  GenmapDestroyVector(w);
  GenmapDestroyVector(r);

  return iter;
}
