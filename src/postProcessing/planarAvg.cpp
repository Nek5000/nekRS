#include "nrs.hpp"
#include "platform.hpp"
#include "linAlg.hpp"
#include "nekInterfaceAdapter.hpp"
#include "postProcessing.hpp"

namespace {
inline int mod1(int i, int n) {
  if(!i) return 0;
  return (i+n-1)%n + 1;
}

void get_exyz(int& ex, int& ey, int& ez,int eg, int nelx, int nely)
{
  ex = mod1(eg, nelx);
  ey = 1 + (mod1(eg, nelx*nely) - 1)/nelx;
  ez = 1 + (eg-1)/(nelx*nely);
}

oogs_t *gtpp_gs_setup(nrs_t *nrs, int nelgx, int nelgy, int nelgz, std::string dir)
{
  mesh_t* mesh = nrs->meshV;
  const auto nelgxy = nelgx*nelgy;
  const auto nelgyz = nelgy*nelgz;
  const auto nelgzx = nelgz*nelgx;

  auto *ids = (hlong *) calloc(mesh->Nlocal, sizeof(hlong)); 

  for(int iel = 0; iel < mesh->Nelements; iel++) {
    const auto eg = nek::lglel(iel) + 1;
    int ex, ey, ez;
    const auto nx1 = mesh->Nq;
    const auto ny1 = mesh->Nq;
    const auto nz1 = mesh->Nq;

    // Enumerate points in the y-z plane
    if(dir == "x") {
       get_exyz(ex,ey,ez,eg,nelgx,nelgyz);
       const auto ex_g = ey;
       for(int k = 0; k < mesh->Nq; k++) {
         for(int j = 0; j < mesh->Nq; j++) {
           for(int i = 0; i < mesh->Nq; i++) {
             const auto id = iel*mesh->Np + k*mesh->Nq*mesh->Nq + j*mesh->Nq + i; 
             ids[id] = (j+1) + ny1*k + ny1*nz1*(ex_g-1);
           }
         }
       } 
    }

    // Enumerate points in the x-z plane
    if(dir == "y") {
       get_exyz(ex,ey,ez,eg,nelgx,nelgy);
       const auto ex_g = (ez-1)*nelgx+ex;
       for(int k = 0; k < mesh->Nq; k++) {
         for(int j = 0; j < mesh->Nq; j++) {
           for(int i = 0; i < mesh->Nq; i++) {
             const auto id = iel*mesh->Np + k*mesh->Nq*mesh->Nq + j*mesh->Nq + i; 
             ids[id] = (k+1) + nz1*i + nx1*nz1*(ex_g-1);
           }
         }
       } 
    }

    // Enumerate points in the x-y plane
    if(dir == "z") {
       get_exyz(ex,ey,ez,eg,nelgxy,1);
       const auto ex_g = ex;
       for(int k = 0; k < mesh->Nq; k++) {
         for(int j = 0; j < mesh->Nq; j++) {
           for(int i = 0; i < mesh->Nq; i++) {
             const auto id = iel*mesh->Np + k*mesh->Nq*mesh->Nq + j*mesh->Nq + i; 
             ids[id] = (i+1) + nx1*j + nx1*ny1*(ex_g-1) + 1;
           }
         }
       } 
    }
  }

#if 0
  dfloat *idsDfloat = (dfloat *) calloc(mesh->Nlocal, sizeof(dfloat));
  for(int i = 0; i < mesh->Nlocal; i++) idsDfloat[i] = ids[i];
  occa::memory o_idsDfloat = platform->device.malloc(1*mesh->Nlocal * sizeof(dfloat), idsDfloat);
  writeFld("ids", 0.0, 1, 1, &o_NULL, &o_NULL, &o_idsDfloat, 1);
#endif

  auto ogsh = ogsSetup(mesh->Nlocal, ids, platform->comm.mpiComm, 1, platform->device.occaDevice());
  free(ids);
  auto oogsh = oogs::setup(ogsh, 6, nrs->fieldOffset, ogsDfloat, NULL, OOGS_AUTO);
  return oogsh;
}

} // namespace

void postProcessing::planarAvg(nrs_t *nrs, const std::string& dir, int NELGX, int NELGY, int NELGZ, int nflds, occa::memory o_avg)
{
  mesh_t* mesh = nrs->meshV;
  static auto firstTime = 1;
  const auto fieldOffsetByte = nrs->fieldOffset * sizeof(dfloat);
  static occa::memory o_avgWeights;
  static oogs_t *oogs_x;
  static oogs_t *oogs_y;
  static oogs_t *oogs_z;

  if(firstTime) {
    o_avgWeights = platform->device.malloc(3*fieldOffsetByte);

    oogs_x = gtpp_gs_setup(nrs, NELGX, NELGY, NELGZ, "x"); 
    auto o_avgWeights_x = o_avgWeights.slice(0*fieldOffsetByte, fieldOffsetByte);
    o_avgWeights_x.copyFrom(mesh->o_LMM, mesh->Nlocal*sizeof(dfloat));
    oogs::startFinish(o_avgWeights_x, 1, mesh->Nlocal, ogsDfloat, ogsAdd, oogs_x);
    platform->linAlg->ady(mesh->Nlocal, 1, o_avgWeights_x);
    platform->linAlg->axmy(mesh->Nlocal, 1, mesh->o_LMM, o_avgWeights_x);

    oogs_y = gtpp_gs_setup(nrs, NELGX*NELGY, 1, NELGZ, "y");
    auto o_avgWeights_y = o_avgWeights.slice(1*fieldOffsetByte, fieldOffsetByte);
    o_avgWeights_y.copyFrom(mesh->o_LMM, mesh->Nlocal*sizeof(dfloat));
    oogs::startFinish(o_avgWeights_y, 1, mesh->Nlocal, ogsDfloat, ogsAdd, oogs_y);
    platform->linAlg->ady(mesh->Nlocal, 1, o_avgWeights_y);
    platform->linAlg->axmy(mesh->Nlocal, 1, mesh->o_LMM, o_avgWeights_y);

    oogs_z = gtpp_gs_setup(nrs, NELGX*NELGY, 1, NELGZ, "z");
    auto o_avgWeights_z = o_avgWeights.slice(2*fieldOffsetByte, fieldOffsetByte);
    o_avgWeights_z.copyFrom(mesh->o_LMM, mesh->Nlocal*sizeof(dfloat));
    oogs::startFinish(o_avgWeights_z, 1, mesh->Nlocal, ogsDfloat, ogsAdd, oogs_z);
    platform->linAlg->ady(mesh->Nlocal, 1, o_avgWeights_z);
    platform->linAlg->axmy(mesh->Nlocal, 1, mesh->o_LMM, o_avgWeights_z);

    firstTime = 0;
  }

  occa::memory o_wghts;
  oogs_t *gsh;
  if(dir == "x") {
    o_wghts = o_avgWeights.slice(0*fieldOffsetByte, fieldOffsetByte);
    gsh = oogs_x;
  } else if(dir == "y") {
    o_wghts = o_avgWeights.slice(1*fieldOffsetByte, fieldOffsetByte);
    gsh = oogs_y;
  } else if(dir == "z") {
    o_wghts = o_avgWeights.slice(2*fieldOffsetByte, fieldOffsetByte);
    gsh = oogs_z;
  } else {
    if (platform->comm.mpiRank == 0) printf("ERROR in planarAvg: Unknown direction!");
    ABORT(EXIT_FAILURE);
  }
  for(int ifld = 0; ifld < nflds; ifld++) {
    auto o_wrk = o_avg.slice(ifld*fieldOffsetByte, fieldOffsetByte);
    platform->linAlg->axmy(mesh->Nlocal, 1, o_wghts, o_wrk);
  } 
  oogs::startFinish(o_avg, nflds, nrs->fieldOffset, ogsDfloat, ogsAdd, gsh);
} 
