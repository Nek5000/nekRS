
struct findpts_el_pt {
  dfloat x[p_D],r[p_D],oldr[p_D],dist2,dist2p,tr;
  dlong  flags;
};

#if p_D==3
struct findpts_el_gface { dfloat *x[p_D], *dxdn[p_D]; };
struct findpts_el_gedge { dfloat *x[p_D], *dxdn1[p_D], *dxdn2[p_D],
                                       *d2xdn1[p_D], *d2xdn2[p_D]; };
struct findpts_el_gpt   { dfloat x[p_D], jac[p_D*p_D], hes[18]; };
#else
struct findpts_el_gedge { dfloat *x[p_D], *dxdn1[p_D]; };
struct findpts_el_gpt   { dfloat x[p_D], jac[p_D*p_D], hes[4]; };
#endif

struct findpts_el_data {
  dlong npt_max;
  struct findpts_el_pt *p; // unused: storage for an single element's data

  dlong n[p_D];
  dfloat *z[p_D];
  void *lag[p_D]; // unused: we don't use the function pointer approach on device
  dfloat *lag_data[p_D];
  dfloat *wtend[p_D];

  const dfloat *x[p_D]; // unused: storage for an single element's data

  dlong side_init; // unused: storage for an single element's data
  dfloat *sides; // unused: storage for an single element's data
#if p_D==3
  struct findpts_el_gface face[2*p_D];   // unused: storage for an single element's data
  struct findpts_el_gedge edge[2*2*p_D]; // unused: storage for an single element's data
#else
  struct findpts_el_gedge edge[2*p_D]; // unused: storage for an single element's data
#endif
  struct findpts_el_gpt pt[1<<p_D]; // unused: storage for an single element's data

  dfloat *work; // unused: workspace for an single element's data
};

struct dbl_range { dfloat min, max; };
struct obbox { dfloat c0[p_D], A[p_D*p_D];
               struct dbl_range x[p_D]; };

struct findpts_local_hash_data {
  dlong hash_n;
  struct dbl_range bnd[p_D];
  dfloat fac[p_D];
  dlong *offset;
  dlong max;
};

struct findpts_local_data {
  dlong ntot;
  const dfloat *elx[p_D];
  struct obbox *obb;
  struct findpts_local_hash_data hd;
  struct findpts_el_data fed;
  dfloat tol;
};


#define   AT(T,var,i)   \
        (T*)(      (char*)var##_base   +(i)*var##_stride   )
#define  CAT(T,var,i) \
  (const T*)((const char*)var##_base   +(i)*var##_stride   )
#define CATD(T,var,i,d) \
  (const T*)((const char*)var##_base[d]+(i)*var##_stride[d])

