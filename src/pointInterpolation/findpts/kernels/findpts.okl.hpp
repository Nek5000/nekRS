
struct findptsElementPoint_t {
  dfloat x[p_D], r[p_D], oldr[p_D], dist2, dist2p, tr;
  dlong flags;
};

struct findptsElementGFace_t {
  dfloat *x[p_D], *dxdn[p_D];
};
struct findptsElementGEdge_t {
  dfloat *x[p_D], *dxdn1[p_D], *dxdn2[p_D], *d2xdn1[p_D], *d2xdn2[p_D];
};
struct findptsElementGPT_t {
  dfloat x[p_D], jac[p_D * p_D], hes[18];
};

struct findptsElementData_t {
  dlong npt_max;
  findptsElementPoint_t *p;
  dlong n[p_D];
  dfloat *z[p_D];
  void *lag[p_D];
  dfloat *lag_data[p_D];
  dfloat *wtend[p_D];

  const dfloat *x[p_D];
  dlong side_init;
  dlong *sides;
  findptsElementGFace_t face[2 * p_D];
  findptsElementGEdge_t edge[2 * 2 * p_D];
  findptsElementGPT_t pt[1 << p_D];

  dfloat *work;
};

struct dbl_range_t {
  dfloat min, max;
};
struct obbox_t {
  dfloat c0[p_D], A[p_D * p_D];
  dbl_range_t x[p_D];
};

struct findptsLocalHashData_t {
  dlong hash_n;
  dbl_range_t bnd[p_D];
  dfloat fac[p_D];
  dlong *offset;
  dlong max;
};