#include "platform.hpp"
#include "linAlg.hpp"
#include "mesh.h"
#include "nekInterfaceAdapter.hpp"
#include <numeric>
#include <algorithm>

namespace
{
inline int mod1(int i, int n)
{
  if (!i) {
    return 0;
  }
  return (i + n - 1) % n + 1;
}

void get_exyz(int &ex, int &ey, int &ez, int eg, int nelx, int nely)
{
  ex = mod1(eg, nelx);
  ey = 1 + (mod1(eg, nelx * nely) - 1) / nelx;
  ez = 1 + (eg - 1) / (nelx * nely);
}

oogs_t *gtpp_gs_setup(mesh_t *mesh, int nelgx, int nelgy, int nelgz, std::string dir)
{
  const auto nelgxy = nelgx * nelgy;
  const auto nelgyz = nelgy * nelgz;
  const auto nelgzx = nelgz * nelgx;

  auto *ids = (hlong *)calloc(mesh->Nlocal, sizeof(hlong));

  for (int iel = 0; iel < mesh->Nelements; iel++) {
    const auto eg = nek::localElementIdToGlobal(iel) + 1;
    int ex, ey, ez;
    const auto nx1 = mesh->Nq;
    const auto ny1 = mesh->Nq;
    const auto nz1 = mesh->Nq;

    // Enumerate points in the y-z plane
    if (dir == "x") {
      get_exyz(ex, ey, ez, eg, nelgx, nelgyz);
      const auto ex_g = ey;
      for (int k = 0; k < mesh->Nq; k++) {
        for (int j = 0; j < mesh->Nq; j++) {
          for (int i = 0; i < mesh->Nq; i++) {
            const auto id = iel * mesh->Np + k * mesh->Nq * mesh->Nq + j * mesh->Nq + i;
            ids[id] = (j + 1) + ny1 * k + ny1 * nz1 * (ex_g - 1);
          }
        }
      }
    }

    // Enumerate points in the x-z plane
    if (dir == "y") {
      get_exyz(ex, ey, ez, eg, nelgx, nelgy);
      const auto ex_g = (ez - 1) * nelgx + ex;
      for (int k = 0; k < mesh->Nq; k++) {
        for (int j = 0; j < mesh->Nq; j++) {
          for (int i = 0; i < mesh->Nq; i++) {
            const auto id = iel * mesh->Np + k * mesh->Nq * mesh->Nq + j * mesh->Nq + i;
            ids[id] = (k + 1) + nz1 * i + nx1 * nz1 * (ex_g - 1);
          }
        }
      }
    }

    // Enumerate points in the x-y plane
    if (dir == "z") {
      get_exyz(ex, ey, ez, eg, nelgxy, 1);
      const auto ex_g = ex;
      for (int k = 0; k < mesh->Nq; k++) {
        for (int j = 0; j < mesh->Nq; j++) {
          for (int i = 0; i < mesh->Nq; i++) {
            const auto id = iel * mesh->Np + k * mesh->Nq * mesh->Nq + j * mesh->Nq + i;
            ids[id] = (i + 1) + nx1 * j + nx1 * ny1 * (ex_g - 1) + 1;
          }
        }
      }
    }
  }

  auto ogsh = ogsSetup(mesh->Nlocal, ids, platform->comm.mpiComm, 1, platform->device.occaDevice());
  free(ids);
  auto oogsh = oogs::setup(ogsh, 6, /* dummy */ mesh->Nlocal, ogsDfloat, NULL, OOGS_AUTO);
  return oogsh;
}

void fusedPlanarAvg(mesh_t *mesh,
                    const std::string &direction,
                    int NELGX,
                    int NELGY,
                    int NELGZ,
                    int nflds,
                    dlong fldOffset,
                    occa::memory o_avg)
{
  static bool issueWarning = true;

  if (!platform->device.deviceAtomic) {
    if (platform->comm.mpiRank == 0 && issueWarning) {
      std::cout << "Device atomics are not supported!\n";
      std::cout << "Relying on slower, separate planar averaging operations.\n";
    }
    issueWarning = false;

    const auto firstDir = direction.substr(0, 1);
    const auto secondDir = direction.substr(1);
    planarAvg(mesh, firstDir, NELGX, NELGY, NELGZ, nflds, fldOffset, o_avg);
    planarAvg(mesh, secondDir, NELGX, NELGY, NELGZ, nflds, fldOffset, o_avg);

    return;
  }

  static occa::memory o_locToGlobE;

  int elemDir = -1;
  if (direction == "xy" || direction == "yx") {
    elemDir = NELGZ;
  }
  if (direction == "xz" || direction == "zx") {
    elemDir = NELGY;
  }
  if (direction == "yz" || direction == "zy") {
    elemDir = NELGX;
  }

  const auto Nlocal = nflds * mesh->Nq * elemDir;
  auto o_scratch = platform->deviceMemoryPool.reserve<dfloat>(Nlocal);

  if (o_locToGlobE.byte_size() == 0) {
    std::vector<dlong> globalElement(mesh->Nelements, 0);
    for (int element = 0; element < mesh->Nelements; ++element) {
      const auto ge = nek::localElementIdToGlobal(element);
      globalElement[element] = ge;
    }

    o_locToGlobE = platform->device.malloc<dlong>(mesh->Nelements);
    o_locToGlobE.copyFrom(globalElement.data());
  }

  auto upperCaseDir = direction;
  std::transform(upperCaseDir.begin(), upperCaseDir.end(), upperCaseDir.begin(), ::toupper);
  std::sort(upperCaseDir.begin(), upperCaseDir.end());

  auto gatherPlanarValuesKernel = platform->kernelRequests.load("gatherPlanarValues" + upperCaseDir);
  auto scatterPlanarValuesKernel = platform->kernelRequests.load("scatterPlanarValues" + upperCaseDir);

  platform->linAlg->fill(Nlocal, 0.0, o_scratch);

  gatherPlanarValuesKernel(mesh->Nelements,
                           nflds,
                           fldOffset,
                           NELGX,
                           NELGY,
                           NELGZ,
                           o_locToGlobE,
                           o_avg,
                           o_scratch);

  platform->comm.allreduce(o_scratch, Nlocal, comm_t::op::sum, platform->comm.mpiComm);

  scatterPlanarValuesKernel(mesh->Nelements,
                            nflds,
                            fldOffset,
                            NELGX,
                            NELGY,
                            NELGZ,
                            o_locToGlobE,
                            o_scratch,
                            o_avg);
}

} // namespace

void planarAvg(mesh_t *mesh,
               const std::string &dir,
               int NELGX,
               int NELGY,
               int NELGZ,
               int nflds,
               dlong fldOffset,
               occa::memory &o_avg)
{
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

  if (dir == "x") {
    o_wghts = o_avgWeight_x;
    gsh = oogs_x;
  } else if (dir == "y") {
    o_wghts = o_avgWeight_y;
    gsh = oogs_y;
  } else if (dir == "z") {
    o_wghts = o_avgWeight_z;
    gsh = oogs_z;
  } else if (dir == "xy" || dir == "yx") {
    o_wghts = o_avgWeight_xy;
  } else if (dir == "xz" || dir == "zx") {
    o_wghts = o_avgWeight_xz;
  } else if (dir == "yz" || dir == "zy") {
    o_wghts = o_avgWeight_yz;
  } else {
    nekrsAbort(platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "Unknown direction!");
  }

  if (!gsh && o_wghts.byte_size() == 0) {

    if (dir == "x") {
      oogs_x = gtpp_gs_setup(mesh, NELGX, NELGY, NELGZ, "x");
      gsh = oogs_x;
      o_avgWeight_x = platform->device.malloc<dfloat>(mesh->Nlocal);
      o_wghts = o_avgWeight_x;
    } else if (dir == "y") {
      oogs_y = gtpp_gs_setup(mesh, NELGX, NELGY, NELGZ, "y");
      gsh = oogs_y;
      o_avgWeight_y = platform->device.malloc<dfloat>(mesh->Nlocal);
      o_wghts = o_avgWeight_y;
    } else if (dir == "z") {
      oogs_z = gtpp_gs_setup(mesh, NELGX * NELGY, 1, NELGZ, "z");
      gsh = oogs_z;
      o_avgWeight_z = platform->device.malloc<dfloat>(mesh->Nlocal);
      o_wghts = o_avgWeight_z;
    } else if (dir == "xy" || dir == "yx") {
      o_avgWeight_xy = platform->device.malloc<dfloat>(mesh->Nlocal);
      o_wghts = o_avgWeight_xy;
    } else if (dir == "xz" || dir == "zx") {
      o_avgWeight_xz = platform->device.malloc<dfloat>(mesh->Nlocal);
      o_wghts = o_avgWeight_xz;
    } else if (dir == "yz" || dir == "zy") {
      o_avgWeight_yz = platform->device.malloc<dfloat>(mesh->Nlocal);
      o_wghts = o_avgWeight_yz;
    } else {
      nekrsAbort(platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "Unknown direction!");
    }

    o_wghts.copyFrom(mesh->o_LMM);
    if (dir.length() > 1) {
      fusedPlanarAvg(mesh, dir, NELGX, NELGY, NELGZ, 1, fldOffset, o_wghts);
    } else {
      oogs::startFinish(o_wghts, 1, mesh->Nlocal, ogsDfloat, ogsAdd, gsh);
    }
    platform->linAlg->ady(mesh->Nlocal, 1, o_wghts);
    platform->linAlg->axmy(mesh->Nlocal, 1, mesh->o_LMM, o_wghts);
  }

  for (int ifld = 0; ifld < nflds; ifld++) {
    auto o_wrk = o_avg.slice(ifld * fldOffset, fldOffset);
    platform->linAlg->axmy(mesh->Nlocal, 1, o_wghts, o_wrk);
  }

  if (dir.length() > 1) {
    fusedPlanarAvg(mesh, dir, NELGX, NELGY, NELGZ, nflds, fldOffset, o_avg);
  } else {
    oogs::startFinish(o_avg, nflds, fldOffset, ogsDfloat, ogsAdd, gsh);
  }
}
