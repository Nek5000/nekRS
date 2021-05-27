#if !defined(nekrs_interp_hpp_)
#define nekrs_interp_hpp_

#include "nrs.hpp"

// Contains data for doing interpolations on a particular mesh
struct interp_data {
  nrs_t *nrs;
  double newton_tol;
  unsigned D;
  void *findpts;
};

// Does the setup to interpolate fields on the given mesh
//   nrs            ... nekRS configuration data
//   newton_tol     ... tolerance for newton solve (use 0 for default)
// returns pointer to the interpolation handle
struct interp_data* interp_setup(nrs_t *nrs, double newton_tol);

// Frees a previously setup interpolation handle
void interp_free(struct interp_data *handle);

// Interpolates the field at the given points
//   fld            ... source field(s) (dfloat[nrs->fieldOffset*nfld])
//   nfld           ... number of fields
//   x              ... array of pointers to the interpolation points (dfloat[D][n])
//   x_stride       ... array of the strides to the interpolation points (dlong[D])
//   n              ... number of points
//   iwk            ... integer working array to hold point location information (dlong[3][nmax])
//   rwk            ... real working array to hold the point local information (dlong[D+1][nmax])
//   nmax           ... leading dimension of iwk and rwk, at least as large as n
//   if_need_pts    ... wheather the interpolation points need to be located (proc,el,r,s,t)
//   handle         ... handle
//   out            ... array of pointers to the output arrays (dfloat[D][n])
//   out_stride     ... array of the strides of the output arrays (dlong[D])
void interp_nfld(dfloat *fld, dlong nfld,
                 dfloat *x[], dlong x_stride[], dlong n,
                 dlong *iwk, dfloat *rwk, dlong nmax,
                 bool if_need_pts, struct interp_data *handle,
                 dfloat *out[], dlong out_stride[]);

// Interpolates the velocity fields at the give points
//   uvw_base       ... array of pointers to the velocity output arrays (dfloat[D][n])
//   uvw_stride     ... array of the strides of the velocity output arrays (dlong[D])
//   xyz_base       ... array of pointers to the coordinate arrays (dfloat[D][n])
//   xyz_stride     ... array of the strides of the coordinate arrays (dlong[D])
//   n              ... number of points to interpolate
//   nrs            ... the NekRS data
void interp_velocity(dfloat *uvw_base[], dlong uvw_stride[],
                     dfloat *xyz_base[], dlong xyz_stride[],
                     int n, nrs_t *nrs);

#endif
