#include "metrics.h"
#include "multigrid.h"
#include "parrsb-impl.h"
#include "sort.h"

#define MM 500

extern void matrix_inverse(int N, double *A);

inline static scalar dot(scalar *y, scalar *x, uint n) {
  scalar result = 0.0;
  for (uint i = 0; i < n; i++)
    result += x[i] * y[i];

  return result;
}

inline static void ortho(scalar *q, uint lelt, ulong n, struct comm *c) {
  uint i;
  scalar sum = 0.0;
  for (i = 0; i < lelt; i++)
    sum += q[i];

  scalar buf;
  comm_allreduce(c, gs_double, gs_add, &sum, 1, &buf);
  sum /= n;

  for (i = 0; i < lelt; i++)
    q[i] -= sum;
}

struct fiedler {
  scalar fiedler;
  uint proc, seq;
  int part0, part1;
};

int power_serial(double *y, uint N, double *A, int verbose) {
  time_t t;
  srand((unsigned)time(&t));

  int i;
  scalar norm = 0.0;
  for (i = 0; i < N; i++) {
    y[i] = (rand() % 50) / 50.0;
    norm += y[i] * y[i];
  }

  scalar normi = 1.0 / sqrt(norm);
  for (i = 0; i < N; i++)
    y[i] *= normi;

  double *Ay = tcalloc(double, N);
  int j, k, l;
  scalar err = 1.0, lambda;
  for (i = 0; i < 100; i++) {
    norm = 0.0;
    for (j = 0; j < N; j++) {
      Ay[j] = 0.0;
      for (k = 0; k < N; k++) {
        Ay[j] += A[j * N + k] * y[k];
      }
      norm += Ay[j] * Ay[j];
    }

    if (i > 0)
      err = (sqrt(norm) - lambda) / lambda;
    lambda = sqrt(norm);

    normi = 1.0 / sqrt(norm);
    for (j = 0; j < N; j++)
      y[j] = Ay[j] * normi;

    if (fabs(err) < 1.e-12)
      break;
  }
  free(Ay);

  return i;
}

int inv_power_serial(double *y, uint N, double *A, int verbose) {
  double *Ainv = tcalloc(double, N *N);
  int j, k;
  for (j = 0; j < N; j++) {
    for (k = 0; k < N; k++)
      Ainv[j * N + k] = A[k * N + j];
  }

  matrix_inverse(N, Ainv);

  for (j = 0; j < N; j++) {
    for (k = 0; k < N; k++)
      A[j * N + k] = Ainv[k * N + j];
  }
  j = power_serial(y, N, Ainv, verbose);

  free(Ainv);

  return j;
}

static int project(scalar *x, uint n, scalar *b, struct laplacian *L,
                   struct mg *d, struct comm *c, int miter, double tol,
                   int null_space, int verbose, buffer *bfr) {
  slong out[2][1], buf[2][1], in = n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  ulong ng = out[1][0];

  if (ng == 0)
    return 0;

  scalar *z = (scalar *)tcalloc(scalar, 6 * n);
  scalar *w = z + n, *r = w + n, *p = r + n, *z0 = p + n, *dz = z0 + n;
  scalar *P = tcalloc(scalar, 2 * (miter + 1) * n);
  scalar *W = P + n * (miter + 1);

  uint i;
  for (i = 0; i < n; i++)
    x[i] = 0, r[i] = b[i];

  scalar rr = dot(r, r, n);
  comm_allreduce(c, gs_double, gs_add, &rr, 1, buf);
  scalar rtol = rr * tol * tol;

  for (i = 0; i < n; i++)
    z[i] = r[i];
  if (null_space)
    ortho(z, n, ng, c);
  scalar rz1 = dot(z, z, n);
  comm_allreduce(c, gs_double, gs_add, &rz1, 1, buf);

  for (i = 0; i < n; i++)
    p[i] = z[i];

  scalar alpha, beta, rzt, rz2;

  uint j, k;
  for (i = 0; i < miter; i++) {
    // mat_vec_csr(w, p, S, gsh, wrk, bfr);
    laplacian(w, L, p, bfr);

    scalar pw = dot(p, w, n);
    comm_allreduce(c, gs_double, gs_add, &pw, 1, buf);
    alpha = rz1 / pw;

    pw = 1 / sqrt(pw);
    for (j = 0; j < n; j++)
      W[i * n + j] = pw * w[j], P[i * n + j] = pw * p[j];

    for (j = 0; j < n; j++)
      x[j] += alpha * p[j], r[j] -= alpha * w[j];

    rr = dot(r, r, n);
    comm_allreduce(c, gs_double, gs_add, &rr, 1, buf);

    if (rr < rtol || sqrt(rr) < tol)
      break;

    for (j = 0; j < n; j++)
      z0[j] = z[j];

    mg_vcycle(z, r, d, c, bfr);

    rzt = rz1;
    if (null_space)
      ortho(z, n, ng, c);
    rz1 = dot(r, z, n);
    comm_allreduce(c, gs_double, gs_add, &rz1, 1, buf);

    for (j = 0; j < n; j++)
      dz[j] = z[j] - z0[j];
    rz2 = dot(r, dz, n);
    comm_allreduce(c, gs_double, gs_add, &rz2, 1, buf);

    if (c->id == 0 && verbose > 0) {
      printf("rr = %lf rtol = %lf rz0 = %lf rz1 = %lf rz2 = %lf\n", rr, rtol,
             rzt, rz1, rz2);
      fflush(stdout);
    }

    beta = rz2 / rzt;
    for (j = 0; j < n; j++)
      p[j] = z[j] + beta * p[j];

    for (k = 0; k < n; k++)
      P[miter * n + k] = 0;

    for (j = 0; j <= i; j++) {
      pw = 0;
      for (k = 0; k < n; k++)
        pw += W[j * n + k] * p[k];
      comm_allreduce(c, gs_double, gs_add, &pw, 1, buf);
      for (k = 0; k < n; k++)
        P[miter * n + k] += pw * P[j * n + k];
    }

    for (k = 0; k < n; k++)
      p[k] -= P[miter * n + k];
  }

  free(z);
  free(P);

  return i == miter ? i : i + 1;
}

// Input z should be orthogonal to 1-vector, have unit norm.
// inverse iteration should not change z.
static int inverse(scalar *y, struct array *elements, int nv, scalar *z,
                   struct comm *gsc, int miter, int mpass, double tol,
                   int factor, int sagg, int grammian, slong nelg,
                   buffer *buf) {
  metric_tic(gsc, RSB_INVERSE_SETUP);
  uint lelt = elements->n;
  struct rsb_element *elems = (struct rsb_element *)elements->ptr;
  struct laplacian *wl = laplacian_init(elems, lelt, nv, GS, gsc, buf);

  // Reserve enough memory in buffer
  size_t wrk = sizeof(ulong) * lelt + sizeof(slong) * nv * lelt;
  buffer_reserve(buf, wrk);

  slong out[2][1], bfr[2][1], in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  ulong *eid = (ulong *)buf->ptr;
  slong *vtx = (slong *)(eid + lelt);
  uint i, j, k, l;
  for (i = k = 0; i < lelt; i++) {
    eid[i] = start + i + 1;
    for (j = 0; j < nv; j++)
      vtx[k++] = elems[i].vertices[j];
  }

  // Setup LAMG preconditioner
  struct crystal cr;
  crystal_init(&cr, gsc);
  struct par_mat *L = par_csr_setup_con(lelt, eid, vtx, nv, 1, gsc, &cr, buf);
  for (uint i = 0; i < L->rn; i++) {
    for (uint j = L->adj_off[i], je = L->adj_off[i + 1]; j < je; j++)
      L->adj_val[j] /= L->diag_val[i];
    L->diag_val[i] = 1.0;
  }
  struct mg *d = mg_setup(L, factor, sagg, &cr, buf);
  crystal_free(&cr);
  metric_toc(gsc, RSB_INVERSE_SETUP);

  scalar *err = tcalloc(scalar, 2 * lelt + miter * (lelt + miter + 1));
  scalar *Z = err + lelt, *M = Z + lelt, *GZ = M + miter * lelt;
  scalar *rhs = GZ + miter * miter, *v = rhs + miter;

  metric_tic(gsc, RSB_INVERSE);
  for (i = 0; i < mpass; i++) {
    int ppfi = project(y, lelt, z, wl, d, gsc, miter, tol, 1, 0, buf);

    ortho(y, lelt, nelg, gsc);

    scalar lambda = dot(y, z, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &lambda, 1, bfr);

    for (uint j = 0; j < lelt; j++)
      err[j] = y[j] - lambda * z[j];
    scalar norme = dot(err, err, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &norme, 1, bfr);
    norme = sqrt(norme);

    scalar norm = dot(y, y, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &norm, 1, bfr);
    scalar normi = 1.0 / sqrt(norm);

    for (j = 0; j < lelt; j++)
      z[j] = y[j] * normi;

    ortho(z, lelt, nelg, gsc);

    int N = i + 1;
    if (grammian == 1) {
      // if k>1;
      //  Z(:,k)=z-Z(:,1:k-1)*(Z(:,1:k-1)'*z);
      //  Z(:,k)=Z(:,k)/norm(Z(:,k));
      // end;
      if (i > 0) {
        // rhs = Z[1:k-1,:]*z
        for (j = 0; j < i; j++) {
          rhs[j] = 0.0;
          for (l = 0; l < lelt; l++)
            rhs[j] += Z[j * lelt + l] * z[l];
        }
        // Global reduction rhs[j]
        comm_allreduce(gsc, gs_double, gs_add, rhs, i, bfr);

        // Z[k,:] = z[:] - Z[:,1:lelt]*rhs[:]
        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] = z[l];

        for (j = 0; j < i; j++) {
          for (l = 0; l < lelt; l++)
            Z[i * lelt + l] = Z[i * lelt + l] - rhs[j] * Z[j * lelt + l];
        }

        // Z[k,:]= Z[k,:]/||Z[k,:]||
        norm = 0.0;
        for (l = 0; l < lelt; l++)
          norm += Z[i * lelt + l] * Z[i * lelt + l];

        comm_allreduce(gsc, gs_double, gs_add, &norm, 1, bfr);
        norm = 1.0 / sqrt(norm);

        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] *= norm;

        // M=Z(1:k,:)*G*Z(1:k,:);
        for (j = 0; j < N; j++) {
          laplacian(GZ, wl, &Z[j * lelt], buf);
          for (k = 0; k < N; k++) {
            M[k * N + j] = 0.0;
            for (l = 0; l < lelt; l++)
              M[k * N + j] += Z[k * lelt + l] * GZ[l];
          }
        }

        // Global reduction of M
        comm_allreduce(gsc, gs_double, gs_add, M, N * N, buf->ptr);

        // Inverse power iterarion on M
        inv_power_serial(v, N, M, 0);

        for (j = 0; j < lelt; j++)
          z[j] = 0.0;

        for (j = 0; j < N; j++) {
          for (k = 0; k < lelt; k++)
            z[k] += Z[j * lelt + k] * v[j];
        }
        ortho(z, lelt, nelg, gsc);
      } else {
        // Z(k,:) = z;
        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] = z[l];
      }
    }

    if (ppfi == 1)
      break;
  }
  metric_toc(gsc, RSB_INVERSE);

  laplacian_free(wl);
  if (L) {
    par_mat_free(L);
    free(L);
  }
  mg_free(d);
  if (err)
    free(err);

  return i;
}

static double sign(scalar a, scalar b) {
  scalar m = b >= 0.0 ? 1.0 : -1.0;
  return fabs(a) * m;
}

static int tqli(scalar *eVectors, scalar *eValues, sint n, scalar *diagonal,
                scalar *upper, int id) {
  if (n == 0)
    return 0;

  scalar *d = tcalloc(scalar, 2 * n), *e = d + n;
  sint i;
  for (i = 0; i < n; i++)
    d[i] = diagonal[i];
  for (i = 0; i < n - 1; i++)
    e[i] = upper[i];
  e[n - 1] = 0.0;

  for (i = 0; i < n; i++) {
    for (uint j = 0; j < n; j++)
      eVectors[i * n + j] = 0;
    eVectors[i * n + i] = 1;
  }

  int j, k, l, iter, m;
  for (l = 0; l < n; l++) {
    iter = 0;
    do {
      for (m = l; m < n - 1; m++) {
        scalar dd = fabs(d[m]) + fabs(d[m + 1]);
        /* Should use a tolerance for this check */
        if (fabs(e[m]) / dd < SCALAR_TOL)
          break;
      }

      if (m != l) {
        if (iter++ == 30) {
          if (id == 0)
            printf("Too many iterations.\n");
          // vec_copy(*eValues, d);
          for (i = 0; i < n; i++)
            eValues[i] = d[i];
          return 1;
        }

        scalar g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        scalar r = sqrt(g * g + 1.0);

        g = d[m] - d[l] + e[l] / (g + sign(r, g));
        scalar s = 1.0, c = 1.0, p = 0.0;

        for (i = m - 1; i >= l; i--) {
          scalar f = s * e[i];
          scalar b = c * e[i];

          if (fabs(f) >= fabs(g)) {
            c = g / f;
            r = sqrt(c * c + 1.0);
            e[i + 1] = f * r;
            s = 1.0 / r;
            c = c * s;
          } else {
            s = f / g;
            r = sqrt(s * s + 1.0);
            e[i + 1] = g * r;
            c = 1.0 / r;
            s = s * c;
          }

          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          /* Find eigenvectors */
          for (k = 0; k < n; k++) {
            f = eVectors[k * n + i + 1];
            eVectors[k * n + i + 1] = s * eVectors[k * n + i] + c * f;
            eVectors[k * n + i] = c * eVectors[k * n + i] - s * f;
          }
          /* Done with eigenvectors */
        }

        if (r < SCALAR_TOL && i >= l)
          continue;

        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }

  /* Orthnormalize eigenvectors -- Just normalize? */
  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      scalar tmp = eVectors[i * n + j];
      eVectors[i * n + j] = eVectors[j * n + i];
      eVectors[j * n + i] = tmp;
    }
  }

  for (k = 0; k < n; k++) {
    e[k] = 0;
    for (uint i = 0; i < n; i++)
      e[k] += eVectors[k * n + i] * eVectors[k * n + i];
    if (e[k] > 0.0)
      e[k] = sqrt(fabs(e[k]));
    scalar scale = 1.0 / e[k];
    for (uint i = 0; i < n; i++)
      eVectors[k * n + i] *= scale;
  }

  // vec_copy(*eValues, d);
  for (i = 0; i < n; i++)
    eValues[i] = d[i];

  free(d);

  return 0;
}

static int lanczos_aux(scalar *diag, scalar *upper, scalar *rr, uint lelt,
                       ulong nelg, int niter, double tol, scalar *f,
                       struct laplacian *gl, struct comm *gsc, buffer *bfr) {
  scalar *r = tcalloc(scalar, 3 * lelt), *p = r + lelt, *w = p + lelt;
  // vec_copy(r, f);
  uint i;
  for (i = 0; i < lelt; i++)
    r[i] = f[i];

  // vec_ortho(gsc, r, nelg);
  ortho(r, lelt, nelg, gsc);

  scalar rtz1 = 1, pap = 0, alpha, beta, rtz2, pap_old;

  scalar rtr = dot(r, r, lelt), buf[2];
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
  scalar rnorm = sqrt(rtr), rtol = rnorm * tol;

  // vec_scale(rr[0], r, rni);
  scalar rni = 1.0 / rnorm;
  for (i = 0; i < lelt; i++)
    rr[0 * lelt + i] = r[i] * rni;

  int iter;
  for (iter = 0; iter < niter; iter++) {
    rtz2 = rtz1, rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0)
      beta = 0.0;

    // add2s1(p,r,beta,n)
    for (i = 0; i < lelt; i++)
      p[i] = beta * p[i] + r[i];

    scalar pp = dot(p, p, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &pp, 1, buf);

    // vec_ortho(gsc, p, nelg);
    ortho(p, lelt, nelg, gsc);

    laplacian(w, gl, p, bfr);

    scalar ww = dot(w, w, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &ww, 1, buf);

    pap_old = pap, pap = dot(w, p, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &pap, 1, buf);
    /*
    if (gsc->id == 0)
      printf("host iter = %d beta = %lf pp = %lf pap = %lf\n", iter, beta, pp,
             pap);
             */

    alpha = rtz1 / pap;
    // vec_axpby(r, r, 1.0, w, -1.0 * alpha);
    for (i = 0; i < lelt; i++)
      r[i] = r[i] - alpha * w[i];

    rtr = dot(r, r, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
    rnorm = sqrt(rtr), rni = 1.0 / rnorm;

    // vec_scale(rr[iter + 1], r, rni);
    for (i = 0; i < lelt; i++)
      rr[(iter + 1) * lelt + i] = r[i] * rni;

    if (iter == 0) {
      diag[iter] = pap / rtz1;
    } else {
      diag[iter] = (beta * beta * pap_old + pap) / rtz1;
      upper[iter - 1] = -beta * pap_old / sqrt(rtz2 * rtz1);
    }

    if (rnorm < rtol) {
      iter++;
      break;
    }
  }

  metric_acc(TOL_FNL, rnorm);
  metric_acc(TOL_TGT, rtol);

  free(r);

  return iter;
}

static int lanczos(scalar *fiedler, struct array *elements, int nv,
                   scalar *initv, struct comm *gsc, int miter, int mpass,
                   double tol, slong nelg, buffer *bfr) {
  metric_tic(gsc, RSB_LANCZOS_SETUP);
  uint lelt = elements->n;
  struct rsb_element *elems = (struct rsb_element *)elements->ptr;
  struct laplacian *wl = laplacian_init(elems, lelt, nv, GS, gsc, bfr);
  metric_toc(gsc, RSB_LANCZOS_SETUP);

  if (nelg < miter)
    miter = nelg;

  scalar *alpha = tcalloc(scalar, 2 * miter - 1), *beta = alpha + miter;
  scalar *rr = tcalloc(scalar, (miter + 1) * lelt);
  scalar *eVectors = tcalloc(scalar, miter * miter);
  scalar *eValues = tcalloc(scalar, miter);
  int iter = miter, ipass;
  for (ipass = 0; iter == miter && ipass < mpass; ipass++) {
    double t = comm_time();
    iter = lanczos_aux(alpha, beta, rr, lelt, nelg, miter, tol, initv, wl, gsc,
                       bfr);
    metric_acc(RSB_LANCZOS, comm_time() - t);

    t = comm_time();
    // Use TQLI and find the minimum eigenvalue and associated vector
    tqli(eVectors, eValues, iter, alpha, beta, gsc->id);
    scalar eValMin = fabs(eValues[0]);
    uint eValMinI = 0;
    for (uint i = 1; i < iter; i++) {
      if (fabs(eValues[i]) < eValMin) {
        eValMin = fabs(eValues[i]);
        eValMinI = i;
      }
    }

    for (uint i = 0; i < lelt; i++) {
      fiedler[i] = 0.0;
      for (uint j = 0; j < iter; j++)
        fiedler[i] += rr[j * lelt + i] * eVectors[eValMinI * iter + j];
    }
    ortho(fiedler, lelt, nelg, gsc);
    for (uint i = 0; i < lelt; i++)
      initv[i] = fiedler[i];
    metric_acc(RSB_LANCZOS_TQLI, comm_time() - t);
  }

  free(alpha), free(rr), free(eVectors), free(eValues);
  laplacian_free(wl);

  return (ipass - 1) * miter + iter;
}

int fiedler(struct array *elements, int nv, parrsb_options *opts,
            struct comm *gsc, buffer *buf, int verbose) {
  metric_tic(gsc, RSB_FIEDLER_SETUP);
  uint lelt = elements->n;
  slong out[2][1], wrk[2][1], in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], nelg = out[1][0];

  scalar *initv = tcalloc(scalar, lelt);
  for (uint i = 0; i < lelt; i++) {
    initv[i] = start + i + 1.0;
    if (start + i < nelg / 2)
      initv[i] += 1000 * nelg;
  }

  ortho(initv, lelt, nelg, gsc);
  scalar rtr = dot(initv, initv, lelt), rni;
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rni);

  rni = 1.0 / sqrt(rtr);
  for (uint i = 0; i < lelt; i++)
    initv[i] *= rni;
  metric_toc(gsc, RSB_FIEDLER_SETUP);

  metric_tic(gsc, RSB_FIEDLER_CALC);
  int iter = 0;
  scalar *f = tcalloc(scalar, lelt);
  switch (opts->rsb_algo) {
  case 0:
    iter = lanczos(f, elements, nv, initv, gsc, opts->rsb_max_iter,
                   opts->rsb_max_passes, opts->rsb_tol, nelg, buf);
    break;
  case 1:
    iter = inverse(f, elements, nv, initv, gsc, opts->rsb_max_iter,
                   opts->rsb_max_passes, opts->rsb_tol, opts->rsb_mg_factor,
                   opts->rsb_mg_sagg, opts->rsb_mg_grammian, nelg, buf);
    break;
  default:
    break;
  }
  metric_toc(gsc, RSB_FIEDLER_CALC);
  metric_acc(RSB_FIEDLER_CALC_NITER, iter);

  scalar norm = 0;
  for (uint i = 0; i < lelt; i++)
    norm += f[i] * f[i];

  scalar normi;
  comm_allreduce(gsc, gs_double, gs_add, &norm, 1, &normi);
  normi = 1.0 / sqrt(norm);

  for (uint i = 0; i < lelt; i++)
    f[i] *= normi;

  struct rsb_element *elems = (struct rsb_element *)elements->ptr;
  for (uint i = 0; i < lelt; i++)
    elems[i].fiedler = f[i];

  if (initv)
    free(initv);
  if (f)
    free(f);
  return 0;
}
