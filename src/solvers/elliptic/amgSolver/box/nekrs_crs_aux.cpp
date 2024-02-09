#include <assert.h>

#include "nekrs_crs.hpp"

static void check_alloc_(void *ptr, const char *file, unsigned line) {
  if (ptr)
    return;
  fprintf(stderr, "check_alloc failure: %s:%d\n", file, line);
  exit(EXIT_FAILURE);
}

#define check_alloc(ptr) check_alloc_(ptr, __FILE__, __LINE__)

static void gen_crs_basis(double *b, int j_, dfloat *z, int Nq, int Np) {
  double *zr = (double *)calloc(Nq, sizeof(double));
  double *zs = (double *)calloc(Nq, sizeof(double));
  double *zt = (double *)calloc(Nq, sizeof(double));
  double *z0 = (double *)calloc(Nq, sizeof(double));
  double *z1 = (double *)calloc(Nq, sizeof(double));
  check_alloc(zr), check_alloc(zs), check_alloc(zt);
  check_alloc(z0), check_alloc(z1);

  for (int i = 0; i < Nq; i++) {
    z0[i] = 0.5 * (1 - z[i]);
    z1[i] = 0.5 * (1 + z[i]);
  }

  memcpy(zr, z0, Nq * sizeof(double));
  memcpy(zs, z0, Nq * sizeof(double));
  memcpy(zt, z0, Nq * sizeof(double));

  int jj = j_ + 1;
  if (jj % 2 == 0)
    memcpy(zr, z1, Nq * sizeof(double));
  if (jj == 3 || jj == 4 || jj == 7 || jj == 8)
    memcpy(zs, z1, Nq * sizeof(double));
  if (jj > 4)
    memcpy(zt, z1, Nq * sizeof(double));

  for (int k = 0; k < Nq; k++) {
    for (int j = 0; j < Nq; j++) {
      for (int i = 0; i < Nq; i++) {
        int n = i + Nq * j + Nq * Nq * k + j_ * Np;
        b[n] = zr[i] * zs[j] * zt[k];
      }
    }
  }

  free(zr), free(zs), free(zt), free(z0), free(z1);
}

static void get_local_crs_galerkin(double *a, int nc, mesh_t *mf,
                                   elliptic_t *ef) {
  int nelt = mf->Nelements, Np = mf->Np;
  size_t size = nelt * Np;
  size_t ncrs = nc * Np;

  double *b_ = tcalloc(double, ncrs);
  dfloat *b = tcalloc(dfloat, ncrs);
  check_alloc(b_), check_alloc(b);
  for (int j = 0; j < nc; j++)
    gen_crs_basis(b_, j, mf->gllz, mf->Nq, mf->Np);

  for (size_t i = 0; i < ncrs; i++)
    b[i] = b_[i];

  dfloat *u = tcalloc(dfloat, size);
  dfloat *w = tcalloc(dfloat, size);
  check_alloc(u), check_alloc(w);

  occa::memory o_u = platform->device.malloc(size * sizeof(dfloat), u);
  occa::memory o_w = platform->device.malloc(size * sizeof(dfloat), w);
  occa::memory o_upf = platform->device.malloc(size * sizeof(pfloat));
  occa::memory o_wpf = platform->device.malloc(size * sizeof(pfloat));

  int i, j, k, e;
  for (j = 0; j < nc; j++) {
    for (e = 0; e < nelt; e++)
      memcpy(&u[e * Np], &b[j * Np], Np * sizeof(dfloat));

    o_u.copyFrom(u);
    platform->copyDfloatToPfloatKernel(mf->Nlocal, o_u, o_upf);
    ellipticAx(ef, mf->Nelements, mf->o_elementList, o_upf, o_wpf,
               pfloatString);
    platform->copyPfloatToDfloatKernel(mf->Nlocal, o_wpf, o_w);
    o_w.copyTo(w);

    for (e = 0; e < nelt; e++) {
      for (i = 0; i < nc; i++) {
        a[i + j * nc + e * nc * nc] = 0.0;
        for (k = 0; k < Np; k++)
          a[i + j * nc + e * nc * nc] += b[k + i * Np] * w[k + e * Np];
      }
    }
  }

  free(b), free(b_), free(w), free(u);
  o_u.free(), o_w.free(), o_upf.free(), o_wpf.free();
}

static void set_mat_ij(uint *ia, uint *ja, int nc, int nelt) {
  uint i, j, e;
  for (e = 0; e < nelt; e++) {
    for (j = 0; j < nc; j++) {
      for (i = 0; i < nc; i++) {
        ia[i + j * nc + nc * nc * e] = e * nc + i;
        ja[i + j * nc + nc * nc * e] = e * nc + j;
      }
    }
  }
}

void jl_setup_aux(uint *ntot_, ulong **gids_, uint *nnz_, uint **ia_,
                  uint **ja_, double **a_, elliptic_t *elliptic,
                  elliptic_t *ellipticf) {
  mesh_t *mesh = elliptic->mesh, *meshf = ellipticf->mesh;
  assert(mesh->Nelements == meshf->Nelements);
  uint nelt = meshf->Nelements, nc = mesh->Np;

  // Set global ids: copy and apply the mask
  uint ntot = *ntot_ = nelt * nc;
  ulong *gids = *gids_ = tcalloc(ulong, ntot);
  check_alloc(gids);

  for (int j = 0; j < nelt * nc; j++)
    gids[j] = mesh->globalIds[j];

  if (elliptic->Nmasked) {
    dlong *mask_ids = (dlong *)calloc(elliptic->Nmasked, sizeof(dlong));
    elliptic->o_maskIds.copyTo(mask_ids, elliptic->Nmasked * sizeof(dlong));
    for (int n = 0; n < elliptic->Nmasked; n++)
      gids[mask_ids[n]] = 0;
    free(mask_ids);
  }

  // Set coarse matrix
  uint nnz = *nnz_ = nc * nc * nelt;
  double *a = *a_ = tcalloc(double, nnz);
  check_alloc(a);

  get_local_crs_galerkin(a, nc, meshf, ellipticf);

  uint *ia = *ia_ = tcalloc(uint, nnz), *ja = *ja_ = tcalloc(uint, nnz);
  check_alloc(ia), check_alloc(ja);
  set_mat_ij(ia, ja, nc, nelt);
}

#undef check_alloc
