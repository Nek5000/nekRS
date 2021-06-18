
#include <mpi.h>
#include "nrs.hpp"
#include "platform.hpp"
#include <vector>

#include "interp.hpp"

interp_data* interp_setup(nrs_t* nrs, double newton_tol)
{

  if (newton_tol < 5e-13) {
    newton_tol = 5e-13;
  }
  int npt_max = 128;
  int bb_tol = 0.01;

  mesh_t* mesh = nrs->meshV;

  dlong nmsh = mesh->N;
  dlong nelm = mesh->Nelements;
  dlong D = mesh->dim;

  // element geometry
  dfloat* elx[3] = {mesh->x, mesh->y, mesh->z};

  // element dimensions
  dlong n1[3] = {mesh->N+1, mesh->N+1, mesh->N+1};

  dlong m1[3] = {2*n1[0], 2*n1[1], 2*n1[2]};

  // used for # of cells in hash tables
  dlong hash_size = nelm*n1[0]*n1[1];
  if (D == 3) hash_size *= n1[2];

  MPI_Comm comm = platform_t::getInstance()->comm.mpiComm;

  ogs_findpts_t* findpts_handle = ogsFindptsSetup(D, comm, elx, n1, nelm, m1, bb_tol,
                                                  hash_size, hash_size, npt_max, newton_tol,
                                                  (occa::device*)&platform_t::getInstance()->device);

  interp_data* handle = new interp_data();
  handle->nrs = nrs;
  handle->newton_tol = newton_tol;
  handle->D = D;
  handle->findpts = findpts_handle;

  return handle;
}


void interp_free(interp_data* handle)
{
  ogsFindptsFree(handle->findpts);
  delete handle;
}

// Accepts any pointer-like type supported by ogsFindptsEval
// I.E. dfloat*, occa::memory (w/ a dtype of dfloat)
template<typename fld_ptr>
void interp_nfld(fld_ptr fld, dlong nfld,
                 dfloat* x[], dlong x_stride[], dlong n,
                 dlong* iwork, dfloat* rwork, dlong nmax,
                 bool if_need_pts, interp_data* handle,
                 dfloat* out[], dlong out_stride[],
                 bool dev_findpts)
{
  assert(n <= nmax);

  dlong*  code  = iwork;
  dlong*  proc  = iwork+2*nmax;
  dlong*  el    = iwork+nmax;
  dfloat* r     = rwork+nmax;
  dfloat* dist2 = rwork;

  dlong D = handle->D;

  // if the location of the points need to be computed
  if (if_need_pts) {
    unsigned nfail = 0;
    // findpts takes strides in terms of bytes, but interp_nfld takes strides in terms of elements
    dlong* x_stride_bytes = (dlong*)malloc(D*sizeof(dlong));
    for (int i = 0; i < D; ++i) x_stride_bytes[i] = x_stride[i]*sizeof(dfloat);
    ogsFindpts(code,  1*sizeof(dlong),
               proc,  1*sizeof(dlong),
               el,    1*sizeof(dlong),
               r,     D*sizeof(dfloat),
               dist2, 1*sizeof(dfloat),
               x,     x_stride_bytes,
               n, handle->findpts,
               dev_findpts);
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
          std::cerr << " WARNING: point not within mesh xy[z]: "
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

    ogsFindptsEval(out[ifld], out_stride[ifld]*sizeof(dfloat),
                   code,      1               *sizeof(dlong),
                   proc,      1               *sizeof(dlong),
                   el,        1               *sizeof(dlong),
                   r,         D               *sizeof(dfloat),
                   n, fld+in_offset, handle->findpts);
  }
}

// instantiations for host and occa memory
template
void interp_nfld(const dfloat* fld, dlong nfld,
                 dfloat* x[], dlong x_stride[], dlong n,
                 dlong* iwk, dfloat* rwk, dlong nmax,
                 bool if_need_pts, interp_data* handle,
                 dfloat* out[], dlong out_stride[],
                 bool dev_findpts);
template
void interp_nfld(dfloat* fld, dlong nfld,
                 dfloat* x[], dlong x_stride[], dlong n,
                 dlong* iwk, dfloat* rwk, dlong nmax,
                 bool if_need_pts, interp_data* handle,
                 dfloat* out[], dlong out_stride[],
                 bool dev_findpts);
template
void interp_nfld(occa::memory fld, dlong nfld,
                 dfloat* x[], dlong x_stride[], dlong n,
                 dlong* iwk, dfloat* rwk, dlong nmax,
                 bool if_need_pts, interp_data* handle,
                 dfloat* out[], dlong out_stride[],
                 bool dev_findpts);


void interp_velocity(dfloat *uvw_base[], dlong uvw_stride[],
                     dfloat *xyz_base[], dlong xyz_stride[],
                     int n, nrs_t *nrs, bool check_occa)
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
  char* workspace = (char*)malloc(sizeof(dfloat)*n*(D+1) + sizeof(int)*n*3);
  dfloat* rwork = (dfloat*)workspace;
  int*    iwork = (int*)(workspace + sizeof(dfloat)*n*(D+1));


  char* workspace_copy;
  dfloat* rwork_copy;
  int* iwork_copy;
  dfloat** uvw_copy;
  if (check_occa) {
    workspace_copy = (char*)malloc(sizeof(dfloat)*n*(D+1) + sizeof(int)*n*3);
    memcpy(workspace_copy, workspace, sizeof(dfloat)*n*(D+1) + sizeof(int)*n*3);
    rwork_copy = (dfloat*)workspace_copy;
    iwork_copy = (int*)(workspace_copy + sizeof(dfloat)*n*(D+1));

    uvw_copy = (dfloat**)malloc(D*sizeof(dfloat*));
    for (int i = 0; i < D; ++i) {
      uvw_copy[i] = (dfloat*)malloc(uvw_stride[i]*n*sizeof(dfloat));
      memcpy(uvw_copy[i], uvw_base[i], uvw_stride[i]*n*sizeof(dfloat));
    }

    nrs->o_U.copyTo(nrs->U);
  }

  MPI_Comm comm = platform_t::getInstance()->comm.mpiComm;
  int mpi_rank = platform_t::getInstance()->comm.mpiRank;
  int mpi_size = platform_t::getInstance()->comm.mpiCommSize;

  occa::memory o_U_dfloat = nrs->o_U.cast(occa::dtype::get<dfloat>());
  interp_nfld(o_U_dfloat, nrs->dim,
              xyz_base, xyz_stride, n,
              iwork, rwork, n, true, interp_handle,
              uvw_base, uvw_stride, true);

  if (check_occa) {
    interp_nfld(nrs->U, nrs->dim,
                xyz_base, xyz_stride, n,
                iwork_copy, rwork_copy, n, true, interp_handle,
                uvw_copy, uvw_stride, false);

    // controls the tolerences for printing warnings
    const dfloat abs_tol = 2e-15;
    const dfloat rel_tol = 1;

    dlong*   code_dev = iwork;
    dlong*   proc_dev = iwork+2*n;
    dlong*     el_dev = iwork+n;
    dfloat*     r_dev = rwork+n;
    dfloat* dist2_dev = rwork;

    dlong*   code_ref = iwork_copy;
    dlong*   proc_ref = iwork_copy+2*n;
    dlong*     el_ref = iwork_copy+n;
    dfloat*     r_ref = rwork_copy+n;
    dfloat* dist2_ref = rwork_copy;

    for(int r = 0; r<mpi_size; ++r) {
      if (r == mpi_rank) {
        for (int i = 0; i < n; ++i) {
          if (code_ref[i] != code_dev[i]
              || (code_ref[i] != 2
                  && (   proc_ref[i] != proc_dev[i]
                      ||   el_ref[i] !=   el_dev[i]
                      || (dist2_dev[i]>1e-16 && dist2_ref[i]*2 < dist2_dev[i])))) {
            printf("%d: WARNING: ogs_findpts varied at %4d: %e, %e, %e\n",
                   mpi_rank, i, xyz_base[0][i*xyz_stride[0]], xyz_base[1][i*xyz_stride[1]], xyz_base[2][i*xyz_stride[2]]);
            printf("%d:                          %1d, %3d, %4d, %e, %e, %e, %e\n",
                   mpi_rank, code_ref[i], proc_ref[i], el_ref[i], r_ref[3*i], r_ref[3*i+1], r_ref[3*i+2], dist2_ref[i]);
            printf("%d:                          %1d, %3d, %4d, %e, %e, %e, %e\n",
                   mpi_rank, code_dev[i], proc_dev[i], el_dev[i], r_dev[3*i], r_dev[3*i+1], r_dev[3*i+2], dist2_dev[i]);
          }
          for (int j = 0; j < D; ++j) {
            dfloat out_dev = uvw_base[j][i*uvw_stride[j]];
            dfloat out_ref = uvw_copy[j][i*uvw_stride[j]];
            if (std::abs(out_ref - out_dev) > abs_tol
                || std::abs(out_ref - out_dev) > rel_tol*std::abs(out_ref)) {
              printf("WARNING: ogs_findpts_eval varied at point %d: %e != %e (diff %e)\n", i, out_ref, out_dev, out_ref-out_dev);
            }
          }
        }
      }
      MPI_Barrier(comm);
    }
    for (int i = 0; i < D; ++i) {
      free(uvw_copy[i]);
    }
    free(uvw_copy);
    free(workspace_copy);
  }

  free(workspace);
}
