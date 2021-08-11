
#include <cfloat>
#include "neknek.hpp"
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "ogsKernels_FINDPTS.hpp"

static void reserveAllocation(nrs_t *nrs, dlong npt) {
  neknek_t *neknek = nrs->neknek;
  const dlong D = nrs->dim;
  occa::device &device = platform_t::getInstance()->device;

  if(neknek->val_interp == nullptr || neknek->npt != npt) {
    if(neknek->val_interp != nullptr) {
      delete [] neknek->val_interp;
      delete [] neknek->proc;
    }

    dfloat *ralloc = new dfloat[npt*(1+D+D+neknek->Nscalar)];
    neknek->val_interp  = ralloc; ralloc += (D+neknek->Nscalar)*npt;
    neknek->eval_buffer = ralloc; ralloc += npt;
    neknek->r           = ralloc; ralloc += npt*D;

    dlong *ialloc = new dlong[npt*3 + nrs->fieldOffset+1];
    neknek->proc        = ialloc; ialloc += npt;
    neknek->code        = ialloc; ialloc += npt;
    neknek->el          = ialloc; ialloc += npt;
    neknek->point_map   = ialloc; ialloc += nrs->fieldOffset+1;

    neknek->o_point_map  = device.malloc(nrs->fieldOffset+1,      occa::dtype::get<dlong>());
    if(npt > 0) {
      neknek->o_val_interp = device.malloc((D+neknek->Nscalar)*npt, occa::dtype::get<dfloat>());
    } else {
      // OCCA doesn't handle length 0 arrays well
      neknek->o_val_interp = device.malloc(1, occa::dtype::get<dfloat>());
    }

    // DEBUG:
    for(int i = 0; i < (D+neknek->Nscalar)*npt; ++i) {
      neknek->val_interp[i] = std::nan("");
    }

    neknek->npt = npt;
  }

}


static void find_interp_points(nrs_t* nrs){
  neknek_t *neknek = nrs->neknek;

  const dlong nsessions = neknek->nsessions;
  const dlong sessionID = neknek->sessionID;
  MPI_Comm global_comm = neknek->global_comm;

  const dlong D = nrs->dim;
  const mesh_t *mesh = nrs->meshV;
  const dlong nmsh = mesh->N;
  const dlong nelm = mesh->Nelements;
  const dlong nfac = mesh->Nfaces;
  const dlong nfpt = mesh->Nfp;

  occa::device &device = platform_t::getInstance()->device;

  // Setup findpts
  dfloat tol = 5e-13;
  dlong npt_max = 128;
  dfloat bb_tol = 0.01;
  dlong n1[3] = {nmsh+1, nmsh+1, nmsh+1};
  dlong nf1[3] = {2*n1[0],2*n1[1],2*n1[2]};
  dlong npt_per_elm = n1[0]*n1[1]*(D==3?n1[2]:1);
  dlong ntot = npt_per_elm*nelm;

  dfloat *elx_null[3] = {nullptr, nullptr, nullptr};
  dfloat *elx     [3] = {mesh->x, mesh->y, mesh->z};

  if (neknek->ogs_handle != nullptr) {
    ogsFindptsFree(neknek->ogs_handle);
  }

  ogs_findpts_t **ogs_handles = new ogs_findpts_t*[nsessions];
  for(dlong i = 0; i < nsessions; ++i) {
    ogs_handles[i] = ogsFindptsSetup(D, global_comm,
                                     (i == sessionID ? elx : elx_null),
                                     n1, (i == sessionID ? nelm : 0),
                                     nf1, bb_tol, ntot, ntot, npt_max, tol,
                                     &device);
  }
  neknek->ogs_handle = ogsFindptsSetup(D, global_comm, elx, n1, nelm,
                                       nf1, bb_tol, ntot, ntot, npt_max, tol,
                                       &device);

  MPI_Barrier(MPI_COMM_WORLD);

  constexpr dlong faceMap[6] = {5, 0, 1, 2, 3, 4};

  dlong num_interp_faces = 0;
  dlong *intflag = (dlong*)nek::ptr("intflag");
  for (dlong e = 0; e < nelm; ++e) {
    for (dlong f = 0; f < nfac; ++f) {
      num_interp_faces += intflag[faceMap[f]+nfac*e]!=0;
    }
  }
  dlong num_interp_bounds = num_interp_faces*nfpt;
  reserveAllocation(nrs, num_interp_bounds);

  dfloat *interp_x = new dfloat[num_interp_bounds*D];
  dlong ip = 0;
  std::fill(neknek->point_map, neknek->point_map+nrs->fieldOffset, -1);
  for (dlong e = 0; e < nelm; ++e) {
    for (dlong f = 0; f < nfac; ++f) {
      
      for(dlong m = 0; m < nfpt; ++m) {
        dlong id = nfac*nfpt*e + nfpt*f + m;
        dlong idM = mesh->vmapM[id];

        if (intflag[faceMap[f]+nfac*e]) {
          for(dlong d = 0; d<D; ++d) {
            interp_x[ip + d*num_interp_bounds] = elx[d][idM];
          }
          neknek->point_map[idM] = ip;
          ++ip;
        }
      }
    }
  }
  neknek->point_map[nrs->fieldOffset] = neknek->npt;
  neknek->o_point_map.copyFrom(neknek->point_map);

  dfloat *ralloc_temp = new dfloat[num_interp_bounds*(2+D)];
  dfloat *dist2      = ralloc_temp; ralloc_temp +=   num_interp_bounds;
  dfloat *dist2_temp = ralloc_temp; ralloc_temp +=   num_interp_bounds;
  dfloat *    r_temp = ralloc_temp; ralloc_temp += D*num_interp_bounds;

  dlong *ialloc_temp = new dlong[num_interp_bounds*3];
  dlong  * proc_temp = ialloc_temp; ialloc_temp +=   num_interp_bounds;
  dlong  * code_temp = ialloc_temp; ialloc_temp +=   num_interp_bounds;
  dlong  *   el_temp = ialloc_temp; ialloc_temp +=   num_interp_bounds;

  for (dlong i = 0; i < num_interp_bounds; ++i) {
    dist2[i] = DBL_MAX;
    neknek->code[i] = 2;
  }

  MPI_Barrier(global_comm);

  dfloat *interp_x_arr[3] = {interp_x, interp_x+num_interp_bounds, interp_x+2*num_interp_bounds};
  dfloat *interp_x_nul[3] = {nullptr, nullptr, nullptr};
  dlong   interp_x_str[3] = {1*sizeof(dfloat), 1*sizeof(dfloat), 1*sizeof(dfloat)};
  for(dlong sess = 0; sess < nsessions; ++sess) {
    ogsFindpts( code_temp,   1*sizeof(dlong),
                proc_temp,   1*sizeof(dlong),
                  el_temp,   1*sizeof(dlong),
                   r_temp,   D*sizeof(dfloat),
               dist2_temp,   1*sizeof(dfloat),
               (sess==sessionID)?interp_x_nul:interp_x_arr, interp_x_str,
               (sess==sessionID)?0:num_interp_bounds,
               ogs_handles[sess]);

    if(sess!=sessionID) {
      for(dlong i = 0; i < num_interp_bounds; ++i) {
        // TODO review this condition
        // ideally should find point farthest from a boundary
        if((code_temp[i] == 0 || code_temp[i] == 1) && dist2[i] > dist2_temp[i]){
          neknek->code[i] =  code_temp[i];
                 dist2[i] = dist2_temp[i];
          neknek->proc[i] =  proc_temp[i];
          neknek->el  [i] =    el_temp[i];
          for(dlong d = 0; d<D; ++d) {
            neknek->r[D*i+d] = r_temp[D*i+d];
          }
        }
      }
    }
  }
  // TODO add warning prints


  for(dlong i = 0; i < nsessions; ++i) {
    ogsFindptsFree(ogs_handles[i]);
  }

  delete [] interp_x;
  delete [] dist2;     //ralloc_temp;
  delete [] proc_temp; //ialloc_temp;

}

void neknek_setup(nrs_t *nrs)
{
  if(platform->options.compareArgs("BUILD ONLY", "TRUE")) {
    int maxSessions;
    platform->options.getArgs("NEKNEK MAX NUM SESSIONS", maxSessions);
    if (maxSessions >= 2) {
      occa::device device = platform_t::getInstance()->device;
      dlong n1[3] = {nrs->meshV->N+1, nrs->meshV->N+1, nrs->meshV->N+1};
      MPI_Comm comm = platform->comm.mpiComm;

      // precompile kernels
      // OCCA automatically garbage collections
      std::pair<occa::kernel, occa::kernel> kernels = ogs::initFindptsKernel(comm, device, 3, n1);
    }
    return;
  }

  neknek_t *neknek = nrs->neknek;
  if(!neknek->connected) {
    reserveAllocation(nrs, 0);
    neknek->point_map[nrs->fieldOffset] = 0;
    neknek->o_point_map.copyFrom(neknek->point_map);
    return;
  }

  const dlong nsessions = neknek->nsessions;
  MPI_Comm global_comm = neknek->global_comm;
  dlong global_rank;
  MPI_Comm_rank(global_comm, &global_rank);

  if(global_rank == 0) printf("configuring neknek with %d sessions\n", nsessions);

  dlong nfld[2];
  nfld[0] = nrs->Nscalar;
  nfld[1] = -nfld[0];
  MPI_Allreduce(MPI_IN_PLACE, nfld, 2, MPI_DLONG, MPI_MIN, global_comm);
  nfld[1] = -nfld[1];
  if (nfld[0] != nfld[1]) {
    if(global_rank == 0) {
      std::cout << "WARNING: varying numbers of scalars; only updating " << nfld[0] << std::endl;
    }
  }
  neknek->Nscalar = nfld[0];

  platform->options.getArgs("NEKNEK CORRECTOR STEPS", neknek->NcorrectorSteps);

  const dlong movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
  dlong globalMovingMesh;
  MPI_Allreduce(&movingMesh, &globalMovingMesh, 1, MPI_DLONG, MPI_MAX, global_comm);
  neknek->globalMovingMesh = globalMovingMesh;

  find_interp_points(nrs);
}

static void field_eval(nrs_t *nrs, dlong field, occa::memory in) {

  neknek_t *neknek = nrs->neknek;
  const dlong D = nrs->dim;
  dfloat *out = neknek->val_interp+field*neknek->npt;

  ogsFindptsEval(out,          1*sizeof(dfloat),
                 neknek->code, 1*sizeof(dlong),
                 neknek->proc, 1*sizeof(dlong),
                 neknek->el,   1*sizeof(dlong),
                 neknek->r,    D*sizeof(dfloat),
                 neknek->npt, in,
                 neknek->ogs_handle);

}

void neknek_update_boundary(nrs_t *nrs)
{
  neknek_t *neknek = nrs->neknek;
  if(!neknek->connected) return;

  if (neknek->globalMovingMesh) {
    find_interp_points(nrs);
  }

  const dlong D = nrs->dim;

  occa::memory o_U = nrs->o_U.cast(occa::dtype::get<dfloat>());
  for(dlong d = 0; d < D; ++d) {
    field_eval(nrs, d, o_U+d*nrs->fieldOffset);
  }
  if (neknek->Nscalar > 0) {
    occa::memory o_S = nrs->cds->o_S.cast(occa::dtype::get<dfloat>());
    for(dlong i = 0; i < neknek->Nscalar; ++i) {
      field_eval(nrs, D+i, nrs->cds->o_S+nrs->cds->fieldOffset[i]);
    }
  }
  // TODO Allow for higher order extrapolation
  // TODO figure out chk_outflow

  neknek->o_val_interp.copyFrom(neknek->val_interp);
}
                                              

