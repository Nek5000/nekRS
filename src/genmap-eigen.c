#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>

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
    x->data[i] = (x->data[i] - beta->data[i] * x->data[i + 1]) / diag->data[i];
  }

  GenmapDestroyVector(diag);
  return 0;
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
      /* Ay = x */
      GenmapSymTriDiagSolve(y, x, alpha, beta);

      /* Normalize by inf-norm(y) */
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

GenmapScalar GenmapSign(GenmapScalar a, GenmapScalar b) {
  GenmapScalar m = 1.0 ? b >= 0. : -1.0;
  return fabs(a) * m;
}

int GenmapTQLI(GenmapHandle h, GenmapVector diagonal, GenmapVector upper,
               GenmapVector **eVectors, GenmapVector *eValues) {
  assert(diagonal->size == upper->size + 1);

  GenmapInt n = diagonal->size;

  GenmapVector d, e;
  GenmapCreateVector(&d, n);
  GenmapCopyVector(d, diagonal);
  GenmapCreateVector(&e, n);
  GenmapCopyVector(e, upper);
  e->data[n - 1] = 0.0;

  /* Create the vector to store eigenvalues */
  GenmapCreateVector(eValues, n);
  /* Init to identity */
  GenmapMalloc(n, eVectors);
  GenmapInt i;
  for(i = 0; i < n; i++) {
    GenmapCreateZerosVector(&(*eVectors)[i], n);
    (*eVectors)[i]->data[i] = 1.0;
  }

  GenmapInt j, k, l, iter, m;

  for(l = 0; l < n; l++) {
    iter = 0;
    do {
      for(m = l; m < n - 1; m++) {
        GenmapScalar dd = fabs(d->data[m]) + fabs(d->data[m + 1]);
        /* Should use a tolerance for this check */
        if(fabs(e->data[m]) / dd < GENMAP_DP_TOL) break;
      }

      if(m != l) {
        if(iter++ == 30) {
          if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0)
            printf("Too many iterations.\n");
          GenmapCopyVector(*eValues, d);
          return 1;
        }

        GenmapScalar g = (d->data[l + 1] - d->data[l]) / (2.0 * e->data[l]);
        GenmapScalar r = sqrt(g * g + 1.0);

        g = d->data[m] - d->data[l] + e->data[l] / (g + GenmapSign(r, g));
        GenmapScalar s = 1.0, c = 1.0, p = 0.0;

        for(i = m - 1; i >= l; i--) {
          GenmapScalar f = s * e->data[i];
          GenmapScalar b = c * e->data[i];
#if defined(GENMAP_PAUL)
          if(fabs(f) >= fabs(g)) {
            c = g / f;
            r = sqrt(c * c + 1.0);
            e->data[i + 1] = f * r;
            s = 1.0 / r;
            c = c * s;
          } else {
            s = f / g;
            r = sqrt(s * s + 1.0);
            e->data[i + 1] = g * r;
            c = 1.0 / r;
            s = s * c;
          }
#else
          e->data[i + 1] = r = sqrt(f * f + g * g);

          if(r < GENMAP_DP_TOL) {
            d->data[i + 1] -= p;
            e->data[m] = 0.0;
            break;
          }
          s = f / r;
          c = g / r;
#endif
          g = d->data[i + 1] - p;
          r = (d->data[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d->data[i + 1] = g + p;
          g = c * r - b;
          /* Find eigenvectors */
          for(k = 0; k < n; k++) {
            f = (*eVectors)[k]->data[i + 1];
            (*eVectors)[k]->data[i + 1] = s * (*eVectors)[k]->data[i] + c * f;
            (*eVectors)[k]->data[i] = c * (*eVectors)[k]->data[i] - s * f;
          }
          /* Done with eigenvectors */
        }

        if(r < GENMAP_DP_TOL && i >= l) continue;

        d->data[l] -= p;
        e->data[l] = g;
        e->data[m] = 0.0;
      }
    } while(m != l);
  }

  /* Orthnormalize eigenvectors -- Just normalize? */
  for(i = 0; i < n; i++) {
    for(j = 0; j < i; j++) {
      GenmapScalar tmp = (*eVectors)[i]->data[j];
      (*eVectors)[i]->data[j] = (*eVectors)[j]->data[i];
      (*eVectors)[j]->data[i] = tmp;
    }
  }

  for(k = 0; k < n; k++) {
    e->data[k] = GenmapDotVector((*eVectors)[k], (*eVectors)[k]);
    if(e->data[k] > 0.0) e->data[k] = sqrt(fabs(e->data[k]));
    GenmapScalar scale = 1.0 / e->data[k];
    GenmapScaleVector((*eVectors)[k], (*eVectors)[k], scale);
  }

  GenmapCopyVector(*eValues, d);

  GenmapDestroyVector(d);
  GenmapDestroyVector(e);

  return 0;
}
