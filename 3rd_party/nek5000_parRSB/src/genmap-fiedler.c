#include <limits.h>
#include <math.h>

#include <genmap-impl.h>
//
// TODO: use a separate function to generate init vector
//
int GenmapFiedlerRQI(genmap_handle h, struct comm *gsc, int max_iter,
                     int global) {
  GenmapInt lelt = genmap_get_nel(h);
  genmap_vector initVec;
  genmap_vector_create(&initVec, lelt);

  struct rsb_element *elements = genmap_get_elements(h);
  GenmapInt i;
  if (global > 0) {
    for (i = 0; i < lelt; i++)
      initVec->data[i] = genmap_get_local_start_index(h) + i + 1;
  } else {
    for (i = 0; i < lelt; i++)
      initVec->data[i] = elements[i].fiedler;
  }

  genmap_vector_ortho_one(gsc, initVec, genmap_get_partition_nel(h));

  GenmapScalar norm = genmap_vector_dot(initVec, initVec), normi;
  comm_allreduce(gsc, gs_double, gs_add, &norm, 1, &normi);
  normi = 1.0 / sqrt(norm);
  genmap_vector_scale(initVec, initVec, normi);

  metric_tic(gsc, LAPLACIAN_INIT);
  GenmapInitLaplacian(h, gsc);
  metric_toc(gsc, LAPLACIAN_INIT);

  mgData d;
  mgSetup(h, gsc, h->M, &d);

  genmap_vector y;
  genmap_vector_create_zeros(&y, lelt);

  metric_tic(gsc, RQI);
  int iter = rqi(h, gsc, d, initVec, max_iter, y);
  metric_toc(gsc, RQI);
  metric_acc(RQI_NITER, iter);

  mgFree(d);

  norm = 0;
  for (i = 0; i < lelt; i++)
    norm += y->data[i] * y->data[i];
  comm_allreduce(gsc, gs_double, gs_add, &norm, 1, &normi);

  genmap_vector_scale(y, y, 1. / sqrt(norm));
  for (i = 0; i < lelt; i++)
    elements[i].fiedler = y->data[i];

  GenmapDestroyVector(y);
  GenmapDestroyVector(initVec);

  return iter;
}

int GenmapFiedlerLanczos(genmap_handle h, struct comm *gsc, int max_iter,
                         int global) {
  GenmapInt lelt = genmap_get_nel(h);

  slong in = lelt;
  slong out[2][1], buf[2][1];
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, buf);
  slong start = out[0][0];
  slong nelg = out[1][0];

  genmap_vector initVec;
  genmap_vector_create(&initVec, lelt);
  struct rsb_element *elements = genmap_get_elements(h);
  GenmapInt i;
  if (global > 0) {
    for (i = 0; i < lelt; i++) {
      initVec->data[i] = start + i + 1.0;
      if (start + i < nelg / 2)
        initVec->data[i] += 1000 * nelg;
    }
  } else {
    for (i = 0; i < lelt; i++)
      initVec->data[i] = elements[i].fiedler;
  }

  genmap_vector alphaVec, betaVec;
  genmap_vector_create(&alphaVec, max_iter);
  genmap_vector_create(&betaVec, max_iter - 1);
  genmap_vector *q = NULL;

  GenmapScalar rtr = genmap_vector_dot(initVec, initVec);
  GenmapScalar rnorm;
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rnorm);
  rnorm = sqrt(rtr);

  genmap_vector_ortho_one(gsc, initVec, nelg);
  rtr = genmap_vector_dot(initVec, initVec);
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rnorm);
  GenmapScalar rni = 1.0 / sqrt(rtr);
  genmap_vector_scale(initVec, initVec, rni);

  int iter;
  metric_tic(gsc, LANCZOS);
  iter = GenmapLanczos(h, gsc, initVec, max_iter, &q, alphaVec, betaVec);
  metric_toc(gsc, LANCZOS);
  metric_acc(LANCZOS_NITER, iter);

  genmap_vector evLanczos, evTriDiag;
  genmap_vector_create(&evTriDiag, iter);

  /* Use TQLI and find the minimum eigenvalue and associated vector */
  genmap_vector *eVectors, eValues;
  GenmapTQLI(h, alphaVec, betaVec, &eVectors, &eValues);

  GenmapScalar eValMin = fabs(eValues->data[0]);
  GenmapInt eValMinI = 0;
  for (i = 1; i < iter; i++) {
    if (fabs(eValues->data[i]) < eValMin) {
      eValMin = fabs(eValues->data[i]);
      eValMinI = i;
    }
  }
  genmap_vector_copy(evTriDiag, eVectors[eValMinI]);

  GenmapInt j;
  genmap_vector_create_zeros(&evLanczos, lelt);
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < iter; j++)
      evLanczos->data[i] += q[j]->data[i] * evTriDiag->data[j];
  }

  GenmapScalar norm = 0, normi;
  for (i = 0; i < lelt; i++)
    norm += evLanczos->data[i] * evLanczos->data[i];

  comm_allreduce(gsc, gs_double, gs_add, &norm, 1, &normi);
  norm = sqrt(norm);

  genmap_vector_scale(evLanczos, evLanczos, 1.0 / norm);
  for (i = 0; i < lelt; i++)
    elements[i].fiedler = evLanczos->data[i];

  GenmapDestroyVector(initVec);
  GenmapDestroyVector(alphaVec);
  GenmapDestroyVector(betaVec);
  GenmapDestroyVector(evLanczos);
  GenmapDestroyVector(evTriDiag);

  GenmapDestroyVector(eValues);
  for (i = 0; i < iter; i++)
    GenmapDestroyVector(eVectors[i]);
  GenmapFree(eVectors);

  for (i = 0; i < iter + 1; i++)
    GenmapDestroyVector(q[i]);
  GenmapFree(q);

  return iter;
}
