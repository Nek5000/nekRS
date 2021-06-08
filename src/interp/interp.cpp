
#include <mpi.h>
#include "nrs.hpp"
#include "platform.hpp"
#include <vector>

#include "interp.hpp"


// uncomment the following line to check the accuracy of the OCCA implementation
// of interp_velocity versus the CPU implementation
//#define check_dev_accuracy

#ifdef check_dev_accuracy
#include "nekInterfaceAdapter.hpp"
#endif



struct interp_data* interp_setup(nrs_t *nrs, double newton_tol)
{

  if (newton_tol < 5e-13) {
    newton_tol = 5e-13;
  }
  int npt_max = 128;
  int bb_tol = 0.01;

  mesh_t *mesh = nrs->meshV;

  dlong nmsh = mesh->N;
  dlong nelm = mesh->Nelements;
  dlong D = mesh->dim;

  // element geometry
  dfloat *elx[3] = {mesh->x, mesh->y, mesh->z};

  // element dimensions
  dlong n1[3] = {mesh->N+1, mesh->N+1, mesh->N+1};

  dlong m1[3] = {2*n1[0], 2*n1[1], 2*n1[2]};

  // used for # of cells in hash tables
  dlong hash_size = nelm*n1[0]*n1[1];
  if (D == 3) hash_size *= n1[2];

  MPI_Comm comm = platform_t::getInstance()->comm.mpiComm;

  ogs_findpts_t *findpts_handle = ogsFindptsSetup(D, comm, elx, n1, nelm, m1, bb_tol,
                                                  hash_size, hash_size, npt_max, newton_tol);

  struct interp_data *handle = new interp_data();
  handle->nrs = nrs;
  handle->newton_tol = newton_tol;
  handle->D = D;
  handle->findpts = findpts_handle;

  return handle;
}


void interp_free(struct interp_data *handle)
{
  ogsFindptsFree(handle->findpts);
  delete handle;
}

// Accepts any pointer-like type supported by ogsFindptsEval
// I.E. dfloat*, occa::memory (w/ a dtype of dfloat)
template<typename fld_ptr>
void interp_nfld(fld_ptr fld, dlong nfld,
                 dfloat *x[], dlong x_stride[], dlong n,
                 dlong *iwk, dfloat *rwk, dlong nmax,
                 bool if_need_pts, struct interp_data *handle,
                 dfloat *out[], dlong out_stride[])
{
  assert(n <= nmax);

  dlong  *code  = iwk;
  dlong  *proc  = iwk+2*nmax;
  dlong  *el    = iwk+nmax;
  dfloat *r     = rwk+nmax;
  dfloat *dist2 = rwk;

  dlong D = handle->D;

  // if the location of the points need to be computed
  if (if_need_pts) {
    unsigned nfail = 0;
    // findpts takes strides in terms of bytes, but interp_nfld takes strides in terms of elements
    dlong *x_stride_bytes = (dlong*)malloc(D*sizeof(dlong));
    for (int i = 0; i < D; ++i) x_stride_bytes[i] = x_stride[i]*sizeof(dfloat);
    ogsFindpts(code,  1*sizeof(dlong),
               proc,  1*sizeof(dlong),
               el,    1*sizeof(dlong),
               r,     D*sizeof(dfloat),
               dist2, 1*sizeof(dfloat),
               x,     x_stride_bytes,
               n, handle->findpts);
    free(x_stride_bytes);

    for (int in = 0; in < n; ++in) {
      if (code[in] == 1) {
        if (dist2[in] > 10*handle->newton_tol) {
          nfail += 1;
          //if (nfail < 5) write(6,'(a,1p4e15.7)')     ' WARNING: point on boundary or outside the mesh xy[z]d^2: ',     xp(in),yp(in),zp(in),rwk(in,1)
          if (nfail < 5){
            std::cerr << " WARNING: point on boundary or outside the mesh xy[z]d^2: "
                      << x[0][in*x_stride[0]] << "," << x[1][in*x_stride[1]] << ", " << x[2][in*x_stride[2]] << ", " << dist2[in] << std::endl;
          }
        }
      } else if (code[in] == 2) {
        nfail += 1;
        //if (nfail < 5) write(6,'(a,1p3e15.7)')        ' WARNING: point not within mesh xy[z]: !',        xp(in),yp(in),zp(in)
        if (nfail < 5){
          std::cerr << " WARNING: point not within mesh xy[z]d^2: "
                    << x[0][in*x_stride[0]] << "," << x[1][in*x_stride[1]] << ", " << x[2][in*x_stride[2]] << std::endl;
        }
      }
    }

    //  nn(1) = iglsum(n,1)
    //  nn(2) = iglsum(nfail,1)
    //  if(nio.eq.0) then
    //    if(nn(2).gt.0 .or. loglevel.gt.2) write(6,1) nn(1),nn(2)
    //1     format('   total number of points = ',i12,/,'   failed = '
    // &         ,i12,/,' done :: intp_nfld')
    //  endif
  }

  for (int ifld = 0; ifld < nfld; ++ifld) {
    dlong in_offset  = ifld*handle->nrs->fieldOffset;

    ogsFindptsEval(out     [ifld], out_stride[ifld]*sizeof(dfloat),
                   code,      1               *sizeof(dlong),
                   proc,      1               *sizeof(dlong),
                   el,        1               *sizeof(dlong),
                   r,         D               *sizeof(dfloat),
                   n, fld+in_offset, handle->findpts);
  }
}

// instantiations for host and occa memory
template
void interp_nfld(const dfloat *fld, dlong nfld,
                 dfloat *x[], dlong x_stride[], dlong n,
                 dlong *iwk, dfloat *rwk, dlong nmax,
                 bool if_need_pts, struct interp_data *handle,
                 dfloat *out[], dlong out_stride[]);
template
void interp_nfld(occa::memory fld, dlong nfld,
                 dfloat *x[], dlong x_stride[], dlong n,
                 dlong *iwk, dfloat *rwk, dlong nmax,
                 bool if_need_pts, struct interp_data *handle,
                 dfloat *out[], dlong out_stride[]);


void interp_velocity(dfloat *uvw_base[], dlong uvw_stride[],
                     dfloat *xyz_base[], dlong xyz_stride[],
                     int n, nrs_t *nrs)
{
  // the interp handle is cached to avoid repeated setups and frees
  static interp_data *interp_handle;
  static bool called = false;
  if (!called || interp_handle->nrs != nrs) {
    if (called) {
      interp_free(interp_handle);
    }
    called = true;
    interp_handle = interp_setup(nrs, 0);
  }

  dlong D = interp_handle->D;
  char *workspace = (char*)malloc(sizeof(dfloat)*n*(D+1) + sizeof(int)*n*3);
  dfloat *rwork = (dfloat*)workspace;
  int    *iwork = (int*)(workspace + sizeof(dfloat)*n*(D+1));

#ifdef check_dev_accuracy
  dfloat **uvw_copy = (dfloat**)malloc(D*sizeof(dfloat*));
  for (int i = 0; i < D; ++i) {
    uvw_copy[i] = (dfloat*)malloc(uvw_stride[i]*n*sizeof(dfloat));
    memcpy(uvw_copy[i], uvw_base[i], uvw_stride[i]*n*sizeof(dfloat));
  }

  nek::ocopyToNek();
#endif //check_dev_accuracy

  occa::memory o_U_dfloat = nrs->o_U.cast(occa::dtype::get<dfloat>());
  interp_nfld(o_U_dfloat, nrs->dim,
              xyz_base, xyz_stride, n,
              iwork, rwork, n, true, interp_handle,
              uvw_base, uvw_stride);

#ifdef check_dev_accuracy
  interp_nfld(nrs->U, nrs->dim,
              xyz_base, xyz_stride, n,
              iwork, rwork, n, true, interp_handle,
              uvw_copy, uvw_stride);

  // controls the tolerences for printing warnings
  const dfloat abs_tol = 2e-15;
  const dfloat rel_tol = 1;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < D; ++j) {
      dfloat out_dev = uvw_base[j][i*uvw_stride[j]];
      dfloat out_ref = uvw_copy[j][i*uvw_stride[j]];
      if (std::abs(out_ref - out_dev) > abs_tol
          || std::abs(out_ref - out_dev) > rel_tol*std::abs(out_ref)) {
        printf("WARNING: ogs_findpts_eval varied at point %d: %e != %e (diff %e)\n", i, out_ref, out_dev, out_ref-out_dev);
      }
    }
  }
  for (int i = 0; i < D; ++i) {
    free(uvw_copy[i]);
  }
  free(uvw_copy);
#endif // check_dev_accuracy

  free(workspace);
}
