#include <genmap-impl.h>

#include <math.h>
#include <stdio.h>
#include <time.h>
//
// Algorithms
//
// Power and inverse power iterations
int GenmapPowerIter(GenmapVector eVector, GenmapVector alpha,
                    GenmapVector beta, GenmapVector init, int iter) {
  assert(alpha->size == beta->size + 1);
  assert(alpha->size == eVector->size);

  GenmapInt  i, j;
  GenmapInt n = alpha->size;

  GenmapVector x, y, t;
  GenmapCreateVector(&x, n);
  GenmapCreateVector(&t, n);
  GenmapCreateVector(&y, n);
  GenmapCopyVector(x, init);

  GenmapScalar norm = GenmapNormVector(x, -1);

  if(n == 1) {
    eVector->data[0] = alpha->data[0];
    return 0;
  } else {
    for(j = 0; j < iter; j++) {
      // y = Ax
      y->data[0] = alpha->data[0] * x->data[0] + beta->data[0] * x->data[1];
      for(i = 1; i < n - 1; i++) {
        y->data[i] = beta->data[i - 1] * x->data[i - 1] + alpha->data[i] *
                     x->data[i] +
                     beta->data[i] * x->data[i + 1];
      }
      y->data[n - 1] = beta->data[n - 2] * x->data[n - 2] + alpha->data[n - 1]
                       *
                       x->data[n - 1];

      GenmapAxpbyVector(t, x, norm, y, -1.0);
      if(GenmapNormVector(t, 2) < GENMAP_TOL && j == iter - 1) {
        GenmapCopyVector(x, y);
        break;
      } else {
        // Normalize by inf-norm(y)
        norm = GenmapNormVector(y, -1);
        GenmapScaleVector(y, y, 1.0 / norm);
        GenmapCopyVector(x, y);
      }
    }
  }

  GenmapCopyVector(eVector, y);

  GenmapDestroyVector(x);
  GenmapDestroyVector(t);
  GenmapDestroyVector(y);

  return j;
}

int GenmapInvPowerIter(GenmapVector eVector, GenmapVector alpha,
                       GenmapVector beta, GenmapVector init, GenmapInt iter) {
  assert(alpha->size == beta->size + 1);
  assert(alpha->size == eVector->size);

  GenmapInt n = alpha->size;

  GenmapVector x, y;

  if(n == 1) {
    eVector->data[0] = alpha->data[0];
    return 0;
  } else {
    GenmapCreateVector(&x, n);
    GenmapCreateVector(&y, n);

    GenmapCopyVector(x, init);
    GenmapInt j;
    for(j = 0; j < iter; j++) {
      // Ay = x
      GenmapSymTriDiagSolve(y, x, alpha, beta);

      // Normalize by inf-norm(y)
      if(j != iter - 1) {
        GenmapScalar lambda = 1.0 / GenmapNormVector(y, -1);
        GenmapScaleVector(y, y, lambda);
      }

      GenmapCopyVector(x, y);
    }
  }

  GenmapCopyVector(eVector, y);

  GenmapDestroyVector(x);
  GenmapDestroyVector(y);

  return 0;
}
//
// Linear solve for Symmetric Tridiagonal Matrix
//
int GenmapSymTriDiagSolve(GenmapVector x, GenmapVector b,
                          GenmapVector alpha,
                          GenmapVector beta) {
  assert((x->size == b->size) && (x->size == alpha->size));
  assert(alpha->size == beta->size + 1);
  assert(b->size > 0);

  GenmapInt n = b->size;

  GenmapVector diag;
  GenmapCreateVector(&diag, n);
  GenmapCopyVector(diag, alpha);

  GenmapCopyVector(x, b);

  GenmapInt i;
  for(i = 0; i < n - 1; i++) {
    GenmapScalar m = (beta->data[i] / diag->data[i]);
    x->data[i + 1] = x->data[i + 1] - m * x->data[i];
    diag->data[i + 1] = diag->data[i + 1] - m * beta->data[i];
  }

  x->data[n - 1] = x->data[n - 1] / diag->data[n - 1];

  for(i = n - 2; i >= 0; i--) {
    x->data[i] = (x->data[i] - beta->data[i] * x->data[i + 1]) /
                 diag->data[i];
  }

  GenmapDestroyVector(diag);
  return 0;
}
//
//
int GenmapOrthogonalizebyOneVector(GenmapHandle h, GenmapComm c, GenmapVector q1, GenmapLong n) {
  GenmapInt i;
  GenmapScalar sum = 0.0;
  for(i = 0;  i < q1->size; i++) {
    sum += q1->data[i];
  }

  GenmapGop(c, &sum, 1, GENMAP_SCALAR, GENMAP_SUM);
  for(i = 0;  i < q1->size; i++) {
    q1->data[i] -= sum / n;
  }

  return 0;
}

int GenmapLanczosLegendary(GenmapHandle h, GenmapComm c, GenmapVector f,
                  GenmapInt niter, GenmapVector **rr, GenmapVector diag,
                  GenmapVector upper) {
  assert(diag->size == niter);
  assert(diag->size == upper->size + 1);
  assert(f->size == h->header->lelt);

  if(h->header->nel < niter) {
    niter = h->header->nel;
    diag->size = niter;
    upper->size = niter - 1;
  }

  GenmapScalar eps = 1.e-5;
  GenmapScalar alpha, beta;
  GenmapScalar rnorm, rtol, rni, rtr, rtz1, rtz2, pap = 0.0, pap_old;
  GenmapVector r, p, w, weights;

  rtz1 = 1.0;

  GenmapInt lelt = h->header->lelt;

  // Store Local Laplacian weights
  GenmapCreateVector(&weights, lelt);
  GenmapCreateOnesVector(&p, lelt);
  GenmapCreateVector(&w, lelt);
  h->AxInit(h, c, weights);

  // Create vector r orthogonalizing init in 1-norm to (1,1,1...)
  GenmapCreateVector(&r, lelt);
  GenmapCopyVector(r, f);
  GenmapOrthogonalizebyOneVector(h, c, r, h->header->nel);
  rtr = GenmapDotVector(r, r);
  GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
  rnorm = sqrt(rtr);
  rtol = rnorm*eps;
  rni = 1.0 / rnorm;

  // Allocate memory for q-vectors
  if(*rr == NULL){
    GenmapMalloc((size_t)(niter+1), rr);
    GenmapInt i;
    for(i = 0; i < niter+1; ++i) (*rr)[i] = NULL; 
  }
  GenmapCreateVector(&(*rr)[0], lelt);

  GenmapScaleVector((*rr)[0], r, rni);

  int iter;
  for(iter = 0; iter < niter; iter++) {
    rtz2 = rtz1;
    rtz1 = rtr;
    beta = rtz1/rtz2;
    if(iter == 0) beta = 0.0;

    GenmapAxpbyVector(p, p, beta, r, 1.0);
    GenmapOrthogonalizebyOneVector(h, c, p, h->header->nel);
    // Multiplication by the laplacian
    h->Ax(h, c, p, weights, w);
    
    pap_old = pap;
    pap = GenmapDotVector(w, p);
    GenmapGop(c, &pap, 1, GENMAP_SCALAR, GENMAP_SUM);
    alpha = rtz1 / pap;
    GenmapAxpbyVector(r, r, 1.0, w, -1.0*alpha);

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
      diag->data[iter] = (beta*beta*pap_old + pap)/rtz1;
      upper->data[iter - 1] = -beta*pap_old/sqrt(rtz2*rtz1);
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

int GenmapLanczos(GenmapHandle h, GenmapComm c, GenmapVector init,
                  GenmapInt iter, GenmapVector **q, GenmapVector alpha,
                  GenmapVector beta) {
  assert(alpha->size == iter);
  assert(alpha->size == beta->size + 1);
  assert(init->size == h->header->lelt);

  if(h->header->nel < iter) {
    iter = h->header->nel;
    alpha->size = iter;
    beta->size = iter - 1;
  }

  GenmapVector q0, q1, u;
  GenmapScalar normq1 = 0., b = 0.;

  GenmapInt lelt = h->header->lelt;

  // Create vector q1 orthogonalizing init in 1-norm to (1,1,1...)
  GenmapCreateVector(&q1, lelt);
  GenmapCopyVector(q1, init);
  GenmapOrthogonalizebyOneVector(h, c, q1, h->header->nel);
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
  if(*q == NULL){
    GenmapMalloc((size_t)(iter+1), q);
    GenmapInt i;
    for(i = 0; i < iter+1; ++i) (*q)[i] = NULL; 
  }

  // Store Local Laplacian weights
  GenmapVector weights;
  GenmapCreateVector(&weights, lelt);
  h->AxInit(h, c, weights);

  int k;
  for(k = 0; k < iter; k++) {
    // Store q1
    GenmapCreateVector(&(*q)[k], lelt);
    GenmapCopyVector((*q)[k], q1);

    // Multiplication by the laplacian
    h->Ax(h, c, q1, weights, u);

    alpha->data[k] = GenmapDotVector(q1, u);
    GenmapGop(c, &alpha->data[k], 1, GENMAP_SCALAR, GENMAP_SUM);

    GenmapAxpbyVector(u, u, 1., q0, -b);
    GenmapAxpbyVector(u, u, 1., q1, -alpha->data[k]);

    b = GenmapDotVector(u, u);
    GenmapGop(c, &b, 1, GENMAP_SCALAR, GENMAP_SUM);
    b = sqrt(b);

    if(k < iter - 1) {
      beta->data[k] = b;
#if defined(GENMAP_DEBUG)
      printf("beta, k: %lf %d\n", b, k);
#endif
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

#if defined(GENMAP_DEBUG)
  printf("k: %d, iter = %d\n", k, iter);
#endif
  return k;
}

void GenmapFiedlerMinMax(GenmapHandle h, GenmapScalar *min,
                         GenmapScalar *max) {
  *min = 1; *max = -1;

  GenmapElements e = GenmapGetElements(h);
  GenmapInt i;
  for(i = 0; i < h->header->lelt; i++) {
    if(e[i].fiedler < *min) {
      *min = e[i].fiedler;
    }
    if(e[i].fiedler > *max) {
      *max = e[i].fiedler;
    }
  }

  GenmapGop(h->local, min, 1, GENMAP_SCALAR, GENMAP_MIN);
  GenmapGop(h->local, max, 1, GENMAP_SCALAR, GENMAP_MAX);
}

GenmapInt GenmapSetProcessorId(GenmapHandle h) {
  GenmapScalar min, max;
  GenmapFiedlerMinMax(h, &min, &max);
  GenmapScalar range = max - min;

  GenmapInt np = h->Np(h->local);
  GenmapInt nbins = np;
  GenmapInt lelt = h->header->lelt;
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

int GenmapFiedler(GenmapHandle h, GenmapComm c, int maxIter,
                  int global) {
  // 1. Do lanczos in local communicator.
  GenmapInt lelt = h->header->lelt;
  GenmapVector initVec, alphaVec, betaVec;

  GenmapCreateVector(&initVec, h->header->lelt);
  GenmapElements elements = GenmapGetElements(h);

  GenmapInt i;
  if(global > 0) {
    for(i = 0;  i < lelt; i++) {
      //initVec->data[i] = (GenmapScalar) elements[i].globalId;
      if(h->header->start + lelt < h->header->nel/2)
        initVec->data[i] = h->header->start + i + 1000. * (h->header->nel/2.0);
      else
        initVec->data[i] = h->header->start + i;
    }
  } else {
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }

  GenmapCreateVector(&alphaVec, maxIter);
  GenmapCreateVector(&betaVec, maxIter - 1);
  GenmapVector *q = NULL;
  // TODO: Lanczos doesn't work well for smaller matrices
  // We need to fix this
  int iter = GenmapLanczosLegendary(h, c, initVec, maxIter, &q, alphaVec, betaVec);
  //int iter = GenmapLanczos(h, c, initVec, maxIter, &q, alphaVec, betaVec);

  // 2. Do inverse power iteration on local communicator and find
  // local Fiedler vector.
  GenmapVector evLanczos, evTriDiag, evInit;
  GenmapCreateVector(&evTriDiag, iter);
  GenmapCreateVector(&evInit, iter);
  // Setup initial vector and orthogonalize in 1-norm to (1,1,1...)
  for(i = 0; i < iter; i++) {
    evInit->data[i] = i + 1;
  }
  GenmapOrthogonalizebyOneVector(h, c, evInit, (GenmapLong)iter);

  GenmapInvPowerIter(evTriDiag, alphaVec, betaVec, evInit, 500);

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
#if defined(GENMAP_DEBUG) && defined(GENMAP_MPI)
    MPI_Barrier(h->local->gsComm.c);
    for(i = 0; i < h->Np(h->local); i++) {
      if(i == h->Id(h->local)) {
        for(j = 0; j < lelt; j++)
          printf("id = "GenmapIntFormat" globalId = "GenmapLongFormat" fiedler = "GenmapScalarFormat"\n",
                 h->Id(h->global),
                 elements[j].globalId, elements[j].fiedler);
      }
      MPI_Barrier(h->local->gsComm.c);
    }
    MPI_Barrier(h->local->gsComm.c);
#endif


  // n. Destory the data structures
  GenmapDestroyVector(initVec);
  GenmapDestroyVector(alphaVec);
  GenmapDestroyVector(betaVec);
  GenmapDestroyVector(evLanczos);
  GenmapDestroyVector(evTriDiag);
  GenmapDestroyVector(evInit);
  for(i = 0; i < iter+1; i++) {
    GenmapDestroyVector(q[i]);
  }
  GenmapFree(q);

  return iter;
}

void GenmapPrimeFactors(GenmapInt n, GenmapInt *pCount,
                        GenmapInt **primes) {
  GenmapInt nLocal = n;
  GenmapInt count = 0;
  GenmapInt countMax = 10;

  GenmapMalloc((size_t)countMax, primes);

  while(nLocal > 1) {
    GenmapInt p;
    for(p = 2; p * p <= nLocal; p++) {
      if(nLocal % p == 0) {
        while(nLocal % p == 0) nLocal /= p;

        (*primes)[count] = p;
        count++;
        if(count == countMax) {
          countMax += countMax / 2 + 1;
          GenmapRealloc((size_t)countMax, primes);
        }
      }
    }
    if(nLocal > 1) {
      (*primes)[count] = nLocal;
      count++;
      nLocal = 1;
    }
  }

  *pCount = count;
}

void GenmapRSB(GenmapHandle h) {
  GenmapInt id = h->Id(h->local);
  GenmapInt np = h->Np(h->local);
  GenmapInt lelt = h->header->lelt;
  GenmapLong nel = h->header->nel;
  GenmapLong start = h->header->start;
  GenmapElements elements = GenmapGetElements(h);

  int maxIter = 5;
  maxIter = (nel / 5 > maxIter) ? nel / 5 : maxIter;
  maxIter = (maxIter > 50) ? 50 : maxIter;

  int iter = maxIter;
  int npass = 50, ipass = 0;

  if(h->Id(h->global) == 0 && h->dbgLevel > 0) 
    printf("running RSB "), fflush(stdout);
#if defined(GENMAP_MPI)
  MPI_Barrier(h->global->gsComm.c);
  double t0 = MPI_Wtime();
#else
  clock_t t0 = clock();
#endif

  // Data needed to use gslib
  struct crystal cr;
  crystal_init(&cr, &(h->local->gsComm));
  GenmapLong out[2][1], buf[2][1];

  buffer buf0 = null_buffer;
  // Calculate the global Fiedler vector, local communicator
  // must be initialized using the global communicator, we never
  // touch global communicator

  while(h->Np(h->local) > 1) {

    if(h->Id(h->global) == 0 && h->dbgLevel > 1) printf("."), fflush(stdout);

    int global = (h->Np(h->local) == h->Np(h->global));
    ipass = 0;
    do {
      iter = GenmapFiedler(h, h->local, maxIter, global);
      ipass++;
      global = 0;
    } while(ipass < npass && iter == maxIter);

    // sort locally according to Fiedler vector
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                  fiedler,
                  TYPE_DOUBLE, globalId, TYPE_INT, &buf0);

    // Sort the Fiedler vector globally
    GenmapSetProcessorId(h);
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc,
                    0, &cr);
    elements = GenmapGetElements(h);
    lelt = h->header->lelt = (GenmapInt)h->elementArray.n;
    // sort locally again -- now we have everything sorted
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                  fiedler,
                  TYPE_DOUBLE, globalId, TYPE_INT, &buf0);

    GenmapLong lelt_ = (GenmapLong)lelt;
    comm_scan(out, &(h->local->gsComm), genmap_gs_long, gs_add, &lelt_, 1,
              buf);
    start = h->header->start = out[0][0];
    nel = h->header->nel = out[1][0];
    id = h->Id(h->local);
    np = h->Np(h->local);
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

    GenmapLong upLimit = idCount * pNel + ((idCount < nrem) ? idCount :
                                           nrem);
    GenmapLong downLimit = start;
    do {
      GenmapInt end = upLimit - start < lelt ? (GenmapInt)(
                        upLimit - start) : lelt;
      GenmapInt i;
      for(i = (GenmapInt)(downLimit - start); i < end;
          i++) elements[i].proc = idCount - 1;
      downLimit = upLimit;
      idCount++;
      upLimit = idCount * pNel + ((idCount < nrem) ? idCount : nrem);
    } while(downLimit - start < lelt);

    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc,
                    0, &cr);
    elements = GenmapGetElements(h);
    lelt = h->header->lelt = (GenmapInt)(h->elementArray.n);

    // sort locally again -- now we have everything sorted
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                  fiedler,
                  TYPE_DOUBLE, globalId, TYPE_INT, &buf0);


    // Now it is time to split the communicator
    GenmapCommExternal local;
#if defined(GENMAP_MPI)
    MPI_Comm_split(h->local->gsComm.c, bin, id, &local);
#else
    local = 0;
#endif
    // finalize the crystal router
    crystal_free(&cr);
    GenmapDestroyComm(h->local);

    // Create new communicator
    GenmapCreateComm(&(h->local), local);
    MPI_Comm_free(&local);
    crystal_init(&cr, &(h->local->gsComm));

    lelt_ = (GenmapLong)lelt;
    comm_scan(out, &(h->local->gsComm), genmap_gs_long, gs_add, &lelt_, 1,
              buf);
    start = h->header->start = out[0][0];
    nel = h->header->nel = out[1][0];
    id = h->Id(h->local);
    np = h->Np(h->local);
    elements = GenmapGetElements(h);
  }

  crystal_free(&cr);
  buffer_free(&buf0);

  double time;
#if defined(GENMAP_MPI)
  MPI_Barrier(h->global->gsComm.c);
  time = MPI_Wtime() - t0;
#else
  time = ((double)clock() - t0) / CLOCKS_PER_SEC;
#endif
  if(h->Id(h->global) == 0 && h->dbgLevel > 0) 
    printf("\nfinished in %lfs\n", time), fflush(stdout);
}
