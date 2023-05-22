#include <limits>
#include "nrs.hpp"
#include "linAlg.hpp"

occa::memory cdsSolve(const int is, cds_t* cds, dfloat time, int stage)
{
  std::string sid = scalarDigitStr(is);

  platform->timer.tic("scalar rhs", 1);  
  mesh_t* mesh = cds->mesh[0];
  if(is) {
    mesh = cds->meshV;
  }

  platform->o_mempool.slice1.copyFrom(cds->o_BF, cds->fieldOffset[is] * sizeof(dfloat), 0,  
                                      cds->fieldOffsetScan[is] * sizeof(dfloat));

  cds->neumannBCKernel(mesh->Nelements,
                       1,
                       mesh->o_sgeo,
                       mesh->o_vmapM,
                       mesh->o_EToB,
                       is,
                       time,
                       cds->fieldOffset[is],
                       cds->EToBOffset,
                       mesh->o_x,
                       mesh->o_y,
                       mesh->o_z,
                       cds->o_Ue,
                       cds->o_S,
                       cds->o_EToB,
                       cds->o_diff,
                       cds->o_rho,
                       *(cds->o_usrwrk),
                       cds->o_BF);

  platform->timer.toc("scalar rhs");

  const occa::memory &o_S0 =
      (platform->options.compareArgs("SCALAR" + sid + " INITIAL GUESS", "EXTRAPOLATION") && stage == 1)
          ? cds->o_Se.slice(cds->fieldOffsetScan[is] * sizeof(dfloat), cds->fieldOffset[is] * sizeof(dfloat))
          : cds->o_S.slice(cds->fieldOffsetScan[is] * sizeof(dfloat), cds->fieldOffset[is] * sizeof(dfloat));
  platform->o_mempool.slice0.copyFrom(o_S0, mesh->Nlocal * sizeof(dfloat));
  auto o_BF_i = cds->o_BF.slice(cds->fieldOffsetScan[is] * sizeof(dfloat), cds->fieldOffset[is] * sizeof(dfloat));
  ellipticSolve(cds->solver[is], o_BF_i, platform->o_mempool.slice0);

  return platform->o_mempool.slice0;
}


