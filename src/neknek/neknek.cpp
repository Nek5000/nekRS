
#include <cfloat>
#include "neknek.hpp"
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "ogsKernels_FINDPTS.hpp"

static void reserveAllocation(nrs_t *nrs, dlong npt) {
  neknek_t *neknek = nrs->neknek;
  const dlong D = nrs->dim;
  occa::device &device = platform_t::getInstance()->device;

  if(neknek->valInterp == nullptr || neknek->npt != npt) {
    if(neknek->valInterp != nullptr) {
      delete [] neknek->valInterp;
      delete [] neknek->proc;
    }

    dfloat *ralloc = new dfloat[npt*(D+D+neknek->Nscalar)];
    neknek->valInterp  = ralloc; ralloc += (D+neknek->Nscalar)*npt;
    neknek->r           = ralloc; ralloc += npt*D;

    dlong *ialloc = new dlong[npt*3 + nrs->fieldOffset+1];
    neknek->proc        = ialloc; ialloc += npt;
    neknek->code        = ialloc; ialloc += npt;
    neknek->el          = ialloc; ialloc += npt;
    neknek->pointMap    = ialloc; ialloc += nrs->fieldOffset+1;

    neknek->o_pointMap  = device.malloc(nrs->fieldOffset+1, occa::dtype::get<dlong>());
    if(npt > 0) {
      neknek->o_valInterp = device.malloc((D+neknek->Nscalar)*npt, occa::dtype::get<dfloat>());
    } else {
      // OCCA can't give length 0 arrays to kernels
      neknek->o_valInterp = device.malloc(1, occa::dtype::get<dfloat>());
    }

    neknek->npt = npt;
  }

}


static void findInterpPoints(nrs_t* nrs){
  neknek_t *neknek = nrs->neknek;

  const dlong nsessions = neknek->nsessions;
  const dlong sessionID = neknek->sessionID;
  MPI_Comm globalComm = neknek->globalComm;

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

  if (neknek->ogsHandle != nullptr) {
    ogsFindptsFree(neknek->ogsHandle);
  }

  ogs_findpts_t **ogsHandles = new ogs_findpts_t*[nsessions];
  for(dlong i = 0; i < nsessions; ++i) {
    ogsHandles[i] = ogsFindptsSetup(D, globalComm,
                                    (i == sessionID ? elx : elx_null),
                                    n1, (i == sessionID ? nelm : 0),
                                    nf1, bb_tol, ntot, ntot, npt_max, tol,
                                    &device);
  }
  neknek->ogsHandle = ogsFindptsSetup(D, globalComm, elx, n1, nelm,
                                      nf1, bb_tol, ntot, ntot, npt_max, tol,
                                      &device);

  constexpr dlong faceMap[6] = {5, 0, 1, 2, 3, 4};

  dlong num_interp_faces = 0;
  dlong *intflag = (dlong*)nek::ptr("intflag");
  for (dlong e = 0; e < nelm; ++e) {
    for (dlong f = 0; f < nfac; ++f) {
      num_interp_faces += intflag[faceMap[f]+nfac*e]!=0;
    }
  }
  dlong npt = num_interp_faces*nfpt;
  reserveAllocation(nrs, npt);

  dfloat *interpX = new dfloat[npt*D];
  dlong ip = 0;
  std::fill(neknek->pointMap, neknek->pointMap+nrs->fieldOffset, -1);
  for (dlong e = 0; e < nelm; ++e) {
    for (dlong f = 0; f < nfac; ++f) {

      for(dlong m = 0; m < nfpt; ++m) {
        dlong id = nfac*nfpt*e + nfpt*f + m;
        dlong idM = mesh->vmapM[id];

        if (intflag[faceMap[f]+nfac*e]) {
          for(dlong d = 0; d<D; ++d) {
            interpX[ip + d*npt] = elx[d][idM];
          }
          neknek->pointMap[idM] = ip;
          ++ip;
        }
      }
    }
  }
  neknek->pointMap[nrs->fieldOffset] = neknek->npt;
  neknek->o_pointMap.copyFrom(neknek->pointMap);

  dfloat *rallocTemp = new dfloat[npt*(2+D)];
  dfloat *dist2      = rallocTemp; rallocTemp +=   npt;
  dfloat *dist2Temp  = rallocTemp; rallocTemp +=   npt;
  dfloat *    rTemp  = rallocTemp; rallocTemp += D*npt;

  dlong *iallocTemp = new dlong[npt*3];
  dlong  * procTemp = iallocTemp; iallocTemp +=   npt;
  dlong  * codeTemp = iallocTemp; iallocTemp +=   npt;
  dlong  *   elTemp = iallocTemp; iallocTemp +=   npt;

  for (dlong i = 0; i < npt; ++i) {
    dist2[i] = DBL_MAX;
    neknek->code[i] = 2;
  }

  dfloat *interpX_3[3] = {interpX, interpX+npt, interpX+2*npt};
  dfloat *nullptr_3[3] = {nullptr, nullptr, nullptr};
  dlong interpXStride[3] = {1*sizeof(dfloat), 1*sizeof(dfloat), 1*sizeof(dfloat)};
  for(dlong sess = 0; sess < nsessions; ++sess) {
    ogsFindpts( codeTemp,   1*sizeof(dlong),
                procTemp,   1*sizeof(dlong),
                  elTemp,   1*sizeof(dlong),
                   rTemp,   D*sizeof(dfloat),
               dist2Temp,   1*sizeof(dfloat),
               (sess==sessionID)?nullptr_3:interpX_3, interpXStride,
               (sess==sessionID)?0:npt,
               ogsHandles[sess]);

    if(sess!=sessionID) {
      for(dlong i = 0; i < npt; ++i) {
        // TODO review this condition
        // ideally should find point farthest from a boundary
        if((codeTemp[i] == 0 || codeTemp[i] == 1) && dist2[i] > dist2Temp[i]){
          neknek->code[i] =  codeTemp[i];
                 dist2[i] = dist2Temp[i];
          neknek->proc[i] =  procTemp[i];
          neknek->el  [i] =    elTemp[i];
          for(dlong d = 0; d<D; ++d) {
            neknek->r[D*i+d] = rTemp[D*i+d];
          }
        }
      }
    }
  }
  // TODO add warning prints


  for(dlong i = 0; i < nsessions; ++i) {
    ogsFindptsFree(ogsHandles[i]);
  }

  delete [] interpX;
  delete [] dist2;    //rallocTemp;
  delete [] procTemp; //iallocTemp;

}

void neknekSetup(nrs_t *nrs)
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
    neknek->pointMap[nrs->fieldOffset] = 0;
    neknek->o_pointMap.copyFrom(neknek->pointMap);
    return;
  }

  const dlong nsessions = neknek->nsessions;
  MPI_Comm globalComm = neknek->globalComm;
  dlong globalRank;
  MPI_Comm_rank(globalComm, &globalRank);

  if(globalRank == 0) printf("configuring neknek with %d sessions\n", nsessions);

  dlong nFields[2];
  nFields[0] = nrs->Nscalar;
  nFields[1] = -nFields[0];
  MPI_Allreduce(MPI_IN_PLACE, nFields, 2, MPI_DLONG, MPI_MIN, globalComm);
  nFields[1] = -nFields[1];
  if (nFields[0] != nFields[1]) {
    if(globalRank == 0) {
      std::cout << "WARNING: varying numbers of scalars; only updating " << nFields[0] << std::endl;
    }
  }
  neknek->Nscalar = nFields[0];

  platform->options.getArgs("NEKNEK CORRECTOR STEPS", neknek->NcorrectorSteps);

  const dlong movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
  dlong globalMovingMesh;
  MPI_Allreduce(&movingMesh, &globalMovingMesh, 1, MPI_DLONG, MPI_MAX, globalComm);
  neknek->globalMovingMesh = globalMovingMesh;

  findInterpPoints(nrs);
}

static void fieldEval(nrs_t *nrs, dlong field, occa::memory in) {

  neknek_t *neknek = nrs->neknek;
  const dlong D = nrs->dim;
  dfloat *out = neknek->valInterp+field*neknek->npt;

  ogsFindptsEval(out,          1*sizeof(dfloat),
                 neknek->code, 1*sizeof(dlong),
                 neknek->proc, 1*sizeof(dlong),
                 neknek->el,   1*sizeof(dlong),
                 neknek->r,    D*sizeof(dfloat),
                 neknek->npt, in,
                 neknek->ogsHandle);

}

void neknekUpdateBoundary(nrs_t *nrs)
{
  neknek_t *neknek = nrs->neknek;
  if(!neknek->connected) return;

  if (neknek->globalMovingMesh) {
    findInterpPoints(nrs);
  }

  const dlong D = nrs->dim;

  occa::memory o_U = nrs->o_U.cast(occa::dtype::get<dfloat>());
  for(dlong d = 0; d < D; ++d) {
    fieldEval(nrs, d, o_U+d*nrs->fieldOffset);
  }
  if (neknek->Nscalar > 0) {
    occa::memory o_S = nrs->cds->o_S.cast(occa::dtype::get<dfloat>());
    for(dlong i = 0; i < neknek->Nscalar; ++i) {
      fieldEval(nrs, D+i, nrs->cds->o_S+nrs->cds->fieldOffset[i]);
    }
  }
  // TODO Allow for higher order extrapolation
  // TODO figure out chk_outflow

  neknek->o_valInterp.copyFrom(neknek->valInterp);
}
