#include <genmap-impl.h>

// Input z should be orthogonal to 1-vector, have unit norm.
// RQI should not change z.
int rqi(genmap_handle h, struct comm *gsc, mgData d, genmap_vector z,
        int max_iter, genmap_vector y) {
  assert(z->size == y->size);

  uint lelt = z->size;
  genmap_vector err;
  genmap_vector_create(&err, lelt);

  int rank = genmap_comm_rank(genmap_global_comm(h));
  GenmapLong nelg = genmap_get_partition_nel(h);

  // Grammian
  GenmapScalar *Z, *GZ, *M, *rhs, *v, *buf;
  GenmapMalloc(max_iter * lelt, &Z);
  GenmapMalloc(lelt, &GZ);
  GenmapMalloc(max_iter * max_iter, &M);
  GenmapMalloc(max_iter, &rhs);
  GenmapMalloc(max_iter, &v);
  GenmapMalloc(max_iter * max_iter, &buf);

  metric_tic(gsc, PROJECT);
  int ppfi = project(h, gsc, d, z, 100, y);
  metric_toc(gsc, PROJECT);
  metric_acc(PROJECT_NITER, ppfi);

  uint i, j, k, l;
  for (i = 0; i < max_iter; i++) {
    GenmapScalar norm = genmap_vector_dot(y, y);
    comm_allreduce(gsc, gs_double, gs_add, &norm, 1, buf);
    GenmapScalar normi = 1.0 / sqrt(norm);

    genmap_vector_axpby(z, z, 0.0, y, normi);
    genmap_vector_ortho_one(gsc, z, nelg);

    int N = i + 1;

    if (h->options->rsb_grammian == 1) {
      // if k>1;
      //  Z(:,k)=z-Z(:,1:k-1)*(Z(:,1:k-1)'*z);
      //  Z(:,k)=Z(:,k)/norm(Z(:,k));
      // end;
      if (i > 0) {
        // rhs = Z[1:k-1,:]*z
        for (j = 0; j < i; j++) {
          rhs[j] = 0.0;
          for (l = 0; l < lelt; l++)
            rhs[j] += Z[j * lelt + l] * z->data[l];
        }
        // Global reduction rhs[j]
        comm_allreduce(gsc, gs_double, gs_add, rhs, i, buf);

        // Z[k,:] = z[:] - Z[:,1:lelt]*rhs[:]
        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] = z->data[l];

        for (j = 0; j < i; j++) {
          for (l = 0; l < lelt; l++)
            Z[i * lelt + l] = Z[i * lelt + l] - rhs[j] * Z[j * lelt + l];
        }

        // Z[k,:]= Z[k,:]/||Z[k,:]||
        norm = 0.0;
        for (l = 0; l < lelt; l++)
          norm += Z[i * lelt + l] * Z[i * lelt + l];

        comm_allreduce(gsc, gs_double, gs_add, &norm, 1, buf);
        norm = 1.0 / sqrt(norm);

        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] *= norm;

        // M=Z(1:k,:)*G*Z(1:k,:);
        for (j = 0; j < N; j++) {
          GenmapLaplacian(h, &Z[j * lelt], GZ);
          for (k = 0; k < N; k++) {
            M[k * N + j] = 0.0;
            for (l = 0; l < lelt; l++)
              M[k * N + j] += Z[k * lelt + l] * GZ[l];
          }
        }

        // Global reduction of M
        comm_allreduce(gsc, gs_double, gs_add, M, N * N, buf);

        // Inverse power iterarion on M
        genmap_inverse_power(v, N, M, 0);

        for (j = 0; j < lelt; j++)
          z->data[j] = 0.0;

        for (j = 0; j < N; j++) {
          for (k = 0; k < lelt; k++)
            z->data[k] += Z[j * lelt + k] * v[j];
        }
        genmap_vector_ortho_one(gsc, z, nelg);
      } else {
        // Z(k,:) = z;
        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] = z->data[l];
      }
    }

    metric_tic(gsc, PROJECT);
    ppfi = project(h, gsc, d, z, 100, y);
    metric_toc(gsc, PROJECT);
    metric_acc(PROJECT_NITER, ppfi);

    genmap_vector_ortho_one(gsc, y, nelg);

    GenmapScalar lambda = genmap_vector_dot(y, z);
    comm_allreduce(gsc, gs_double, gs_add, &lambda, 1, buf);

    genmap_vector_axpby(err, y, 1.0, z, -lambda);
    GenmapScalar norme = genmap_vector_dot(err, err);
    comm_allreduce(gsc, gs_double, gs_add, &norme, 1, buf);
    norme = sqrt(norme);

    if (ppfi == 1)
      break;
  }

  GenmapFree(Z);
  GenmapFree(GZ);
  GenmapFree(M);
  GenmapFree(rhs);
  GenmapFree(v);
  GenmapFree(buf);

  GenmapDestroyVector(err);

  return i + 1;
}
