#include "nrs.hpp"
#include "platform.hpp"
#include "linAlg.hpp"
#include "nekInterfaceAdapter.hpp"
#include "postProcessing.hpp"
#include <numeric>
#include <algorithm>

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
  writeFld("id" + dir, 0.0, 1, 1, &o_NULL, &o_NULL, &o_idsDfloat, 1);
#endif

  auto ogsh = ogsSetup(mesh->Nlocal, ids, platform->comm.mpiComm, 1, platform->device.occaDevice());
  free(ids);
  auto oogsh = oogs::setup(ogsh, 6, nrs->fieldOffset, ogsDfloat, NULL, OOGS_AUTO);
  return oogsh;
}

void fusedPlanarAvg(nrs_t *nrs, const std::string & direction, int NELGX, int NELGY, int NELGZ, int nflds, occa::memory o_avg)
{
  static_assert(std::is_same<dlong,int>::value, "dlong != int");

  static bool issueWarning = true;

  if (!platform->device.deviceAtomic) {
    if (platform->comm.mpiRank == 0 && issueWarning) {
      std::cout << "Device atomics are not supported!\n";
      std::cout << "Relying on slower, separate planar averaging operations.\n";
    }
    issueWarning = false;

    const auto firstDir = direction.substr(0, 1);
    const auto secondDir = direction.substr(1);
    postProcessing::planarAvg(nrs, firstDir, NELGX, NELGY, NELGZ, nflds, o_avg);
    postProcessing::planarAvg(nrs, secondDir, NELGX, NELGY, NELGZ, nflds, o_avg);

    return;
  }

  const auto * mesh = nrs->meshV;
  static occa::memory o_locToGlobE;

  static occa::memory o_scratch;
  static occa::memory h_scratch;
  static dfloat * scratch;

  int elemDir = -1;
  if(direction == "xy" || direction == "yx"){
    elemDir = NELGZ;
  }
  if(direction == "xz" || direction == "zx"){
    elemDir = NELGY;
  }
  if(direction == "yz" || direction == "zy"){
    elemDir = NELGX;
  }

  const auto Nwords = nflds * mesh->Nq * elemDir;
  const auto Nbytes = Nwords * sizeof(dfloat);

  if(o_scratch.size() < Nbytes){
    if(o_scratch.size()) o_scratch.free();
    if(h_scratch.size()) h_scratch.free();
    {
      h_scratch = platform->device.mallocHost(Nbytes);
      scratch = (dfloat *)h_scratch.ptr();
    }
    o_scratch = platform->device.malloc(Nbytes);
  }

  if(o_locToGlobE.size() == 0){
    std::vector<dlong> globalElement(mesh->Nelements, 0);
    for(int element = 0; element < mesh->Nelements; ++element){
      const auto ge = nek::lglel(element);
      globalElement[element] = ge;
    }

    o_locToGlobE = platform->device.malloc(mesh->Nelements * sizeof(dlong), globalElement.data());
  }

  auto upperCaseDir = direction;
  std::transform(upperCaseDir.begin(), upperCaseDir.end(), upperCaseDir.begin(), ::toupper);
  std::sort(upperCaseDir.begin(), upperCaseDir.end());

  auto gatherPlanarValuesKernel = platform->kernels.get("gatherPlanarValues" + upperCaseDir);
  auto scatterPlanarValuesKernel = platform->kernels.get("scatterPlanarValues" + upperCaseDir);

  platform->linAlg->fill(Nwords, 0.0, o_scratch);

  gatherPlanarValuesKernel(mesh->Nelements,
    nflds,
    nrs->fieldOffset,
    NELGX,
    NELGY,
    NELGZ,
    o_locToGlobE,
    o_avg,
    o_scratch);

  platform->comm.allreduce(o_scratch, Nwords, comm_t::type::dfloat, comm_t::op::sum, platform->comm.mpiComm);

  scatterPlanarValuesKernel(mesh->Nelements,
    nflds,
    nrs->fieldOffset,
    NELGX,
    NELGY,
    NELGZ,
    o_locToGlobE,
    o_scratch,
    o_avg);

}

} // namespace


void postProcessing::planarAvg(nrs_t *nrs, const std::string& dir, int NELGX, int NELGY, int NELGZ, int nflds, occa::memory o_avg)
{
  mesh_t* mesh = nrs->meshV;
  const auto fieldOffsetByte = nrs->fieldOffset * sizeof(dfloat);

  static occa::memory o_avgWeight_x;
  static occa::memory o_avgWeight_y;
  static occa::memory o_avgWeight_z;

  static occa::memory o_avgWeight_xy;
  static occa::memory o_avgWeight_xz;
  static occa::memory o_avgWeight_yz;

  static oogs_t *oogs_x = nullptr;
  static oogs_t *oogs_y = nullptr;
  static oogs_t *oogs_z = nullptr;

  occa::memory o_wghts;
  oogs_t *gsh = nullptr;

  if(dir == "x") {
    o_wghts = o_avgWeight_x;
    gsh = oogs_x;
  } else if(dir == "y") {
    o_wghts = o_avgWeight_y;
    gsh = oogs_y;
  } else if(dir == "z") {
    o_wghts = o_avgWeight_z;
    gsh = oogs_z;
  } else if(dir == "xy" || dir == "yx") {
    o_wghts = o_avgWeight_xy;
  } else if(dir == "xz" || dir == "zx") {
    o_wghts = o_avgWeight_xz;
  } else if(dir == "yz" || dir == "zy") {
    o_wghts = o_avgWeight_yz;
  } else {
    if (platform->comm.mpiRank == 0) printf("ERROR in planarAvg: Unknown direction!");
    ABORT(EXIT_FAILURE);
  }

  if(!gsh && o_wghts.size() == 0) {

    if(dir == "x"){
      oogs_x = gtpp_gs_setup(nrs, NELGX, NELGY, NELGZ, "x"); 
      gsh = oogs_x;
      o_avgWeight_x = platform->device.malloc(fieldOffsetByte);
      o_wghts = o_avgWeight_x;
    }
    else if(dir == "y"){
      oogs_y = gtpp_gs_setup(nrs, NELGX, NELGY, NELGZ, "y");
      gsh = oogs_y;
      o_avgWeight_y = platform->device.malloc(fieldOffsetByte);
      o_wghts = o_avgWeight_y;
    }
    else if(dir == "z"){
      oogs_z = gtpp_gs_setup(nrs, NELGX*NELGY, 1, NELGZ, "z");
      gsh = oogs_z;
      o_avgWeight_z = platform->device.malloc(fieldOffsetByte);
      o_wghts = o_avgWeight_z;
    }
    else if(dir == "xy" || dir == "yx"){
      o_avgWeight_xy = platform->device.malloc(fieldOffsetByte);
      o_wghts = o_avgWeight_xy;
    }
    else if(dir == "xz" || dir == "zx"){
      o_avgWeight_xz = platform->device.malloc(fieldOffsetByte);
      o_wghts = o_avgWeight_xz;
    }
    else if(dir == "yz" || dir == "zy"){
      o_avgWeight_yz = platform->device.malloc(fieldOffsetByte);
      o_wghts = o_avgWeight_yz;
    }
    else{
      if (platform->comm.mpiRank == 0) printf("ERROR in planarAvg: Unknown direction!");
      ABORT(EXIT_FAILURE);
    }

    o_wghts.copyFrom(mesh->o_LMM, mesh->Nlocal*sizeof(dfloat));
    if(dir.length() > 1){
      fusedPlanarAvg(nrs, dir, NELGX, NELGY, NELGZ, 1, o_wghts);
    } else {
      oogs::startFinish(o_wghts, 1, mesh->Nlocal, ogsDfloat, ogsAdd, gsh);
    }
    platform->linAlg->ady(mesh->Nlocal, 1, o_wghts);
    platform->linAlg->axmy(mesh->Nlocal, 1, mesh->o_LMM, o_wghts);

  }

  for(int ifld = 0; ifld < nflds; ifld++) {
    auto o_wrk = o_avg.slice(ifld*fieldOffsetByte, fieldOffsetByte);
    platform->linAlg->axmy(mesh->Nlocal, 1, o_wghts, o_wrk);
  } 

  if(dir.length() > 1){
    fusedPlanarAvg(nrs, dir, NELGX, NELGY, NELGZ, nflds, o_avg);
  } else {
    oogs::startFinish(o_avg, nflds, nrs->fieldOffset, ogsDfloat, ogsAdd, gsh);
  }
} 
