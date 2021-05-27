
#include <mpi.h>
#include "nrs.hpp"
#include "platform.hpp"
#include <vector>

#include "interp.hpp"

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

  void *findpts_handle = ogsFindptsSetup(D, comm, elx, n1, nelm, m1, bb_tol,
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
  ogsFindptsFree((ogs_findpts_t*)handle->findpts);
  delete handle;
}

void interp_nfld(dfloat *fld, dlong nfld,
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

  dlong D = handle->nrs->dim;

  unsigned nfail = 0;
  if (if_need_pts) {
    // findpts takes strides in terms of bytes, but interp_nfld takes strides in terms of elements
    dlong *x_stride_bytes = (dlong*)malloc(D*sizeof(dlong));
    for (int i = 0; i < D; ++i) x_stride_bytes[i] = x_stride[i]*sizeof(dfloat);
    ogsFindpts(code,  1*sizeof(dlong),
               proc,  1*sizeof(dlong),
               el,    1*sizeof(dlong),
               r,     D*sizeof(dfloat),
               dist2, 1*sizeof(dfloat),
               x,     x_stride_bytes,
               n, (ogs_findpts_t*)handle->findpts);
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
  }

  for (int ifld = 0; ifld < nfld; ++ifld) {
     dlong in_offset  = ifld*handle->nrs->fieldOffset;

     ogsFindptsEval(out[ifld], out_stride[ifld]*sizeof(dfloat),
                    code,      1               *sizeof(dlong),
                    proc,      1               *sizeof(dlong),
                    el,        1               *sizeof(dlong),
                    r,         D               *sizeof(dfloat),
                    n, fld+in_offset, (ogs_findpts_t*)handle->findpts);
  }

//  nn(1) = iglsum(n,1)
//  nn(2) = iglsum(nfail,1)
//  if(nio.eq.0) then
//    if(nn(2).gt.0 .or. loglevel.gt.2) write(6,1) nn(1),nn(2)
//1     format('   total number of points = ',i12,/,'   failed = '
// &         ,i12,/,' done :: intp_nfld')
//  endif
}


void interp_velocity(dfloat *uvw_base[], dlong uvw_stride[],
                     dfloat *xyz_base[], dlong xyz_stride[],
                     int n, nrs_t *nrs)
{
  dlong D = nrs->dim;
  char *workspace = (char*)malloc(sizeof(dfloat)*n*(D+1) + sizeof(int)*n*3);
  dfloat *rwork = (dfloat*)workspace;
  int    *iwork = (int*)(workspace + sizeof(dfloat)*n*(D+1));

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

  interp_nfld(nrs->U, nrs->dim,
              xyz_base, xyz_stride, n,
              iwork, rwork, n, true, interp_handle,
              uvw_base, uvw_stride);

  free(workspace);
}
