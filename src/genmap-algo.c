#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
//
// Orthogonalize by 1-vector (vector of all 1's)
//
int GenmapOrthogonalizebyOneVector(GenmapHandle h, GenmapComm c,
                                   GenmapVector q1, GenmapLong n) {
  GenmapInt i;
  GenmapScalar sum = 0.0;
  for(i = 0;  i < q1->size; i++) {
    sum += q1->data[i];
  }

  GenmapGop(c, &sum, 1, GENMAP_SCALAR, GENMAP_SUM);
  sum /= n;
  for(i = 0;  i < q1->size; i++) {
    q1->data[i] -= sum;
  }

  return 0;
}
//
// Lanczos version used in Paul's original genmap code.
//
int GenmapLanczosLegendary(GenmapHandle h, GenmapComm c, GenmapVector f,
                           GenmapInt niter, GenmapVector **rr, GenmapVector diag,
                           GenmapVector upper) {
  assert(diag->size == niter);
  assert(diag->size == upper->size + 1);
  assert(f->size == GenmapGetNLocalElements(h));

  if(GenmapGetNGlobalElements(h) < niter) {
    niter = GenmapGetNGlobalElements(h);
    diag->size = niter;
    upper->size = niter - 1;
  }

#if defined(GENMAP_PAUL)
  GenmapScalar eps = 1.e-5;
#else
  GenmapScalar eps = 1.e-12;
#endif
  GenmapScalar alpha, beta;
  GenmapScalar rnorm, rtol, rni, rtr, rtz1, rtz2, pap = 0.0, pap_old;
  GenmapVector r, p, w, weights;

  rtz1 = 1.0;

  GenmapScalar tmp;
  GenmapInt lelt = GenmapGetNLocalElements(h);

  // Store Local Laplacian weights
  GenmapCreateVector(&weights, lelt);
  GenmapCreateZerosVector(&p, lelt);
  GenmapCreateVector(&w, lelt);
  GenmapInitLaplacian(h, c, weights);

  // Create vector r orthogonalizing init in 1-norm to (1,1,1...)
  GenmapCreateVector(&r, lelt);
  GenmapCopyVector(r, f);
  GenmapOrthogonalizebyOneVector(h, c, r, GenmapGetNGlobalElements(h));
  rtr = GenmapDotVector(r, r);
  GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
  rnorm = sqrt(rtr);
  rtol = rnorm * eps;
  rni = 1.0 / rnorm;

  // Allocate memory for q-vectors
  if(*rr == NULL) {
    GenmapMalloc((size_t)(niter + 1), rr);
    GenmapInt i;
    for(i = 0; i < niter + 1; ++i)(*rr)[i] = NULL;
  }
  GenmapCreateVector(&(*rr)[0], lelt);

  GenmapScaleVector((*rr)[0], r, rni);

  int iter;
  for(iter = 0; iter < niter; iter++) {
    rtz2 = rtz1;
    rtz1 = rtr;
    beta = rtz1 / rtz2;
    if(iter == 0) beta = 0.0;

    GenmapAxpbyVector(p, p, beta, r, 1.0);
    GenmapOrthogonalizebyOneVector(h, c, p, GenmapGetNGlobalElements(h));

    // Multiplication by the laplacian
    GenmapLaplacian(h, c, p, weights, w);
    GenmapScaleVector(w, w, -1.0);

    pap_old = pap;
    pap = GenmapDotVector(w, p);
    GenmapGop(c, &pap, 1, GENMAP_SCALAR, GENMAP_SUM);
    alpha = rtz1 / pap;
    GenmapAxpbyVector(r, r, 1.0, w, -1.0 * alpha);

    rtr = GenmapDotVector(r, r);
    GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
    if(rtr < 0) {
      diag->size = iter + 1;
      upper->size = iter;
      iter = iter + 1;
      break;
    }
    rnorm = sqrt(rtr);
    rni = 1.0 / rnorm;
    GenmapCreateVector(&(*rr)[iter + 1], lelt);
    GenmapScaleVector((*rr)[iter + 1], r, rni);

    if(iter == 0) {
      diag->data[iter] = pap / rtz1;
    } else {
      diag->data[iter] = (beta * beta * pap_old + pap) / rtz1;
      upper->data[iter - 1] = -beta * pap_old / sqrt(rtz2 * rtz1);
    }

    if(rnorm < rtol)  {
      diag->size = iter + 1;
      upper->size = iter;
      iter = iter + 1;
      break;
    }
  }

  GenmapDestroyVector(weights);
  GenmapDestroyVector(p);
  GenmapDestroyVector(w);
  GenmapDestroyVector(r);

  return iter;
}
//
// Lanczos version in Introduction to Sci. Comp by Prof. Heath.
//
int GenmapLanczos(GenmapHandle h, GenmapComm c, GenmapVector init,
                  GenmapInt iter, GenmapVector **q, GenmapVector alpha,
                  GenmapVector beta) {
  assert(alpha->size == iter);
  assert(alpha->size == beta->size + 1);
  assert(init->size == GenmapGetNLocalElements(h));

  if(GenmapGetNGlobalElements(h) < iter) {
    iter = GenmapGetNGlobalElements(h);
    alpha->size = iter;
    beta->size = iter - 1;
  }

  GenmapVector q0, q1, u;
  GenmapScalar normq1 = 0., b = 0.;

  GenmapInt lelt = GenmapGetNLocalElements(h);

  // Create vector q1 orthogonalizing init in 1-norm to (1,1,1...)
  GenmapCreateVector(&q1, lelt);
  GenmapCopyVector(q1, init);
  GenmapOrthogonalizebyOneVector(h, c, q1, GenmapGetNGlobalElements(h));
  normq1 = GenmapDotVector(q1, q1);
  GenmapGop(c, &normq1, 1, GENMAP_SCALAR, GENMAP_SUM);
  normq1 = sqrt(normq1);
  GenmapScaleVector(q1, q1, 1. / normq1);

  // Create vector u
  GenmapCreateVector(&u, lelt);

  // Set q_0 and beta_0 to zero (both uses 0-indexing)
  GenmapCreateZerosVector(&q0, lelt);
  beta->data[0] = 0.;

  // Allocate memory for q-vectors
  if(*q == NULL) {
    GenmapMalloc((size_t)iter, q);
    GenmapInt i;
    for(i = 0; i < iter; ++i)
      (*q)[i] = NULL;
  }

  // Store Local Laplacian weights
  GenmapVector weights;
  GenmapCreateVector(&weights, lelt);
  GenmapInitLaplacian(h, c, weights);

  int k;
  for(k = 0; k < iter; k++) {
    // Store q1
    GenmapCreateVector(&(*q)[k], lelt);
    GenmapCopyVector((*q)[k], q1);

    // Multiplication by the laplacian
    GenmapLaplacian(h, c, q1, weights, u);

    alpha->data[k] = GenmapDotVector(q1, u);
    GenmapGop(c, &alpha->data[k], 1, GENMAP_SCALAR, GENMAP_SUM);

    GenmapAxpbyVector(u, u, 1., q0, -b);
    GenmapAxpbyVector(u, u, 1., q1, -alpha->data[k]);

    b = GenmapDotVector(u, u);
    GenmapGop(c, &b, 1, GENMAP_SCALAR, GENMAP_SUM);
    b = sqrt(b);

    if(k < iter - 1) {
      beta->data[k] = b;
      if(fabs(b) < normq1 * GENMAP_TOL) {
        beta->size = k;
        alpha->size = k + 1;
        k = k + 1;
        break;
      }

      GenmapCopyVector(q0, q1);
      GenmapScaleVector(q1, u, 1. / beta->data[k]);
    }
  }

  GenmapDestroyVector(q0);
  GenmapDestroyVector(q1);
  GenmapDestroyVector(u);
  GenmapDestroyVector(weights);

  return k;
}

void GenmapFiedlerMinMax(GenmapHandle h, GenmapScalar *min, GenmapScalar *max) {
  *min = 1; *max = -1;

  GenmapElements e = GenmapGetElements(h);
  GenmapInt i;
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    if(e[i].fiedler < *min) {
      *min = e[i].fiedler;
    }
    if(e[i].fiedler > *max) {
      *max = e[i].fiedler;
    }
  }

  GenmapGop(GenmapGetLocalComm(h), min, 1, GENMAP_SCALAR, GENMAP_MIN);
  GenmapGop(GenmapGetLocalComm(h), max, 1, GENMAP_SCALAR, GENMAP_MAX);
}

void GenmapGlobalIdMinMax(GenmapHandle h, GenmapLong *min, GenmapLong *max) {
  *min = LONG_MAX; *max = LONG_MIN;

  GenmapElements e = GenmapGetElements(h);
  GenmapInt i;
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    if(e[i].globalId < *min) {
      *min = e[i].globalId;
    }
    if(e[i].globalId > *max) {
      *max = e[i].globalId;
    }
  }

  GenmapGop(GenmapGetLocalComm(h), min, 1, GENMAP_SCALAR, GENMAP_MIN);
  GenmapGop(GenmapGetLocalComm(h), max, 1, GENMAP_SCALAR, GENMAP_MAX);
}

GenmapInt GenmapSetFiedlerBin(GenmapHandle h) {
  GenmapScalar min, max;
  GenmapFiedlerMinMax(h, &min, &max);
  GenmapScalar range = max - min;

  GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));
  GenmapInt nbins = np;
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapElements elements = GenmapGetElements(h);

  GenmapElements p, e;
  for(p = elements, e = p + lelt; p != e; p++) {
    GenmapInt id;
    for(id = 0; id < np; id++) {
      GenmapScalar start = min + (range * id) / nbins;
      GenmapScalar end = min + (range * (id + 1)) / nbins;
      if(start <= p->fiedler && p->fiedler < end) {
        p->proc = id;
        break;
      }
    }
    if(id == np) p->proc = np - 1;
  }

  return 0;
}

GenmapInt GenmapSetGlobalIdBin(GenmapHandle h) {
  GenmapLong min, max;
  GenmapGlobalIdMinMax(h, &min, &max);
  GenmapLong range = max - min;

  GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));
  GenmapInt nbins = np;
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapElements elements = GenmapGetElements(h);

  GenmapElements p, e;
  for(p = elements, e = p + lelt; p != e; p++) {
    GenmapInt id;
    for(id = 0; id < np; id++) {
      GenmapScalar start = min + (range * id) / nbins;
      GenmapScalar end = min + (range * (id + 1)) / nbins;
      if(start <= p->globalId && p->globalId < end) {
        p->proc = id;
        break;
      }
    }
    if(id == np) p->proc = np - 1;
  }

  return 0;
}

int GenmapFiedler(GenmapHandle h, GenmapComm c, int maxIter,
                  int global) {
  // 1. Do lanczos in local communicator.
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapVector initVec, alphaVec, betaVec;

  GenmapCreateVector(&initVec, GenmapGetNLocalElements(h));
  GenmapElements elements = GenmapGetElements(h);

  GenmapInt i;
#if defined(GENMAP_PAUL)
  if(global > 0) {
    for(i = 0;  i < lelt; i++) {
      if(GenmapGetLocalStartIndex(h) + i + 1  < GenmapGetNGlobalElements(h) / 2)
        initVec->data[i] = GenmapGetLocalStartIndex(h) + i + 1 + 1000. *
                           GenmapGetNGlobalElements(h);
      else
        initVec->data[i] = GenmapGetLocalStartIndex(h) + i + 1;
    }
  } else {
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }
#else
  if(global > 0) {
    initVec->data[i] = (GenmapScalar) elements[i].globalId;
  } else {
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }
#endif

  GenmapCreateVector(&alphaVec, maxIter);
  GenmapCreateVector(&betaVec, maxIter - 1);
  GenmapVector *q = NULL;

#if defined(GENMAP_PAUL)
  GenmapOrthogonalizebyOneVector(h, c, initVec, GenmapGetNGlobalElements(h));
  GenmapScalar rtr = GenmapDotVector(initVec, initVec);
  GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScalar rni = 1.0 / sqrt(rtr);
  GenmapScaleVector(initVec, initVec, rni);
  int iter = GenmapLanczosLegendary(h, c, initVec, maxIter, &q, alphaVec,
                                    betaVec);
#else
  int iter = GenmapLanczos(h, c, initVec, maxIter, &q, alphaVec, betaVec);
#endif

  GenmapVector evLanczos, evTriDiag;
  GenmapCreateVector(&evTriDiag, iter);
#if defined(GENMAP_PAUL)
  // 2. Use TQLI and find the minimum eigenvalue and associated vector
  GenmapVector *eVectors, eValues;
  GenmapTQLI(h, alphaVec, betaVec, &eVectors, &eValues);

  GenmapScalar eValMin = fabs(eValues->data[0]);
  GenmapInt eValMinI = 0;
  for(GenmapInt i = 1; i < iter; i++) {
    if(fabs(eValues->data[i]) < eValMin) {
      eValMin = fabs(eValues->data[i]);
      eValMinI = i;
    }
  }

  GenmapCopyVector(evTriDiag, eVectors[eValMinI]);
#else
  // 2. Do inverse power iteration on local communicator and find
  // local Fiedler vector.
  GenmapVector evInit;
  GenmapCreateVector(&evInit, iter);

  // Setup initial vector and orthogonalize in 1-norm to (1,1,1...)
  for(i = 0; i < iter; i++) {
    evInit->data[i] = i + 1;
  }
  GenmapOrthogonalizebyOneVector(h, c, evInit, (GenmapLong)iter);

  GenmapInvPowerIter(evTriDiag, alphaVec, betaVec, evInit, 100);
#endif

  // Multiply tri-diagonal matrix by [q1, q2, ...q_{iter}]
  GenmapInt j;
  GenmapCreateZerosVector(&evLanczos, lelt);
  for(i = 0; i < lelt; i++) {
    for(j = 0; j < iter; j++) {
      evLanczos->data[i] += q[j]->data[i] * evTriDiag->data[j];
    }
  }

  GenmapScalar lNorm = 0;
  for(i = 0; i < lelt; i++) {
    lNorm += evLanczos->data[i] * evLanczos->data[i];
  }

  GenmapGop(c, &lNorm, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScaleVector(evLanczos, evLanczos, 1. / sqrt(lNorm));
  for(i = 0; i < lelt; i++) {
    elements[i].fiedler = evLanczos->data[i];
  }

  // n. Destory the data structures
  GenmapDestroyVector(initVec);
  GenmapDestroyVector(alphaVec);
  GenmapDestroyVector(betaVec);
  GenmapDestroyVector(evLanczos);
  GenmapDestroyVector(evTriDiag);
#if defined(GENMAP_PAUL)
  GenmapDestroyVector(eValues);
  for(i = 0; i < iter; i++) {
    GenmapDestroyVector(eVectors[i]);
  }
  GenmapFree(eVectors);
#else
  GenmapDestroyVector(evInit);

  for(i = 0; i < iter; i++) {
    GenmapDestroyVector(q[i]);
  }
  GenmapFree(q);
#endif

  return iter;
}

void GenmapBinSort(GenmapHandle h, int field, buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  if(field == 0) { // Fiedler
    // sort locally according to Fiedler vector
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
                  TYPE_DOUBLE, globalId, TYPE_LONG, buf0);
    // Sort the Fiedler vector globally
    GenmapSetFiedlerBin(h);
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0, &(h->cr));
    elements = GenmapGetElements(h);
    lelt =  (GenmapInt)h->elementArray.n;
    GenmapSetNLocalElements(h, lelt);
    // sort locally again -- now we have everything sorted
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
                  TYPE_DOUBLE, globalId, TYPE_LONG, buf0);
  } else {
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                  globalId, TYPE_LONG, globalId, TYPE_LONG, buf0);
    // Sort the Fiedler vector globally
    GenmapSetGlobalIdBin(h);
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc,
                    0, &(h->cr));
    elements = GenmapGetElements(h);
    lelt = (GenmapInt)(h->elementArray.n);
    GenmapSetNLocalElements(h, lelt);
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                  globalId, TYPE_LONG, globalId, TYPE_LONG, buf0);
  }
}

void GenmapRSB(GenmapHandle h) {
  GenmapInt id = GenmapCommRank(GenmapGetLocalComm(h));
  GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapLong nel = GenmapGetNGlobalElements(h);
  GenmapLong start = GenmapGetLocalStartIndex(h);
  GenmapElements elements = GenmapGetElements(h);

  int maxIter = 50;
  int iter = maxIter;
  int npass = 50, ipass = 0;

  if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0 && h->dbgLevel > 0)
    printf("running RSB "), fflush(stdout);
#if defined(GENMAP_MPI)
  MPI_Barrier(GenmapGetGlobalComm(h)->gsComm.c);
  double t0 = MPI_Wtime();
#else
  clock_t t0 = clock();
#endif

  crystal_init(&(h->cr), &(h->local->gsComm));
  GenmapLong out[2][1], buf[2][1];

  buffer buf0 = null_buffer;
  // Calculate the global Fiedler vector, local communicator
  // must be initialized using the global communicator, we never
  // touch global communicator

  while(GenmapCommSize(GenmapGetLocalComm(h)) > 1) {
    if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0
        && h->dbgLevel > 1) printf("."), fflush(stdout);

#if defined(GENMAP_PAUL)
    int global = 1;
#else
    int global = (GenmapCommSize(GenmapGetLocalComm(h)) == GenmapCommSize(
                    GenmapGetGlobalComm(h)));
#endif

    ipass = 0;
    do {
      iter = GenmapFiedler(h, GenmapGetLocalComm(h), maxIter, global);
      ipass++;
      global = 0;
    } while(ipass < npass && iter == maxIter);

    GenmapBinSort(h, 0, &buf0);

    GenmapLong lelt_ = (GenmapLong)GenmapGetNLocalElements(h);
    comm_scan(out, &(h->local->gsComm), genmap_gs_long, gs_add, &lelt_, 1, buf);
    start = out[0][0]; GenmapSetLocalStartIndex(h, start);
    nel = out[1][0]; GenmapSetNGlobalElements(h, nel);
    id = GenmapCommRank(GenmapGetLocalComm(h));
    np = GenmapCommSize(GenmapGetLocalComm(h));
    elements = GenmapGetElements(h);

    GenmapInt bin;
    if(id < (np + 1) / 2)
      bin = 0;
    else
      bin = 1;

    GenmapInt pNel = (GenmapInt)(nel / np);
    GenmapInt nrem = (GenmapInt)(nel - pNel * np);
    GenmapInt idCount = 0;
    while(idCount * pNel + ((idCount < nrem) ? idCount : nrem) < start)
      idCount++;

    GenmapLong upLimit = idCount * pNel + ((idCount < nrem) ? idCount : nrem);
    GenmapLong downLimit = start;
    do {
      GenmapInt end = upLimit - start < lelt ? (GenmapInt)(upLimit - start) : lelt;
      GenmapInt i;
      for(i = (GenmapInt)(downLimit - start); i < end; i++)
        elements[i].proc = idCount - 1;
      downLimit = upLimit;
      idCount++;
      upLimit = idCount * pNel + ((idCount < nrem) ? idCount : nrem);
    } while(downLimit - start < lelt);

    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0, &(h->cr));
    elements = GenmapGetElements(h);
    lelt = (GenmapInt)(h->elementArray.n); GenmapSetNLocalElements(h, lelt);

    // sort locally again -- now we have everything sorted
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
                  TYPE_DOUBLE, globalId, TYPE_LONG, &buf0);

    // Now it is time to split the communicator
    GenmapCommExternal local;
#if defined(GENMAP_MPI)
    MPI_Comm_split(h->local->gsComm.c, bin, id, &local);
#else
    local = 0;
#endif
    // finalize the crystal router
    crystal_free(&(h->cr));
    GenmapDestroyComm(h->local);

    // Create new communicator
    GenmapCreateComm(&h->local, local);
    MPI_Comm_free(&local);
    crystal_init(&(h->cr), &(h->local->gsComm));

#if defined(GENMAP_PAUL)
    GenmapBinSort(h, 1, &buf0);
    lelt = GenmapGetNLocalElements(h);
#endif

    lelt_ = (GenmapLong)lelt;
    comm_scan(out, &(GenmapGetLocalComm(h)->gsComm), genmap_gs_long, gs_add, &lelt_,
              1, buf);
    start = out[0][0]; GenmapSetLocalStartIndex(h, start);
    nel = h->nel = out[1][0];
    id = GenmapCommRank(GenmapGetLocalComm(h));
    np = GenmapCommSize(GenmapGetLocalComm(h));
    elements = GenmapGetElements(h);

  }

  crystal_free(&(h->cr));
  buffer_free(&buf0);

  double time;
#if defined(GENMAP_MPI)
  MPI_Barrier(GenmapGetGlobalComm(h)->gsComm.c);
  time = MPI_Wtime() - t0;
#else
  time = ((double)clock() - t0) / CLOCKS_PER_SEC;
#endif
  if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0 && h->dbgLevel > 0)
    printf("\nfinished in %lfs\n", time), fflush(stdout);
}
