#include <limits>
#include "nrs.hpp"
#include "linAlg.hpp"

occa::memory cdsSolve(const int is, cds_t* cds, double time, int stage)
{
  std::string sid = scalarDigitStr(is);

  mesh_t* mesh = cds->mesh[0];
  if(is) {
    mesh = cds->meshV;
  }

  platform->timer.tic("scalar rhs", 1);

  auto o_rhs = platform->o_memPool.reserve<dfloat>(cds->fieldOffset[is]);
  o_rhs.copyFrom(cds->o_BF, cds->fieldOffset[is], 0, cds->fieldOffsetScan[is]);

  cds->neumannBCKernel(mesh->Nelements,
                       1,
                       mesh->o_sgeo,
                       mesh->o_vmapM,
                       mesh->o_EToB,
                       is,
                       time,
                       cds->fieldOffset[is],
                       0,
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
                       o_rhs);

  platform->timer.toc("scalar rhs");

  auto o_S = [&]()
  {
     auto o_S0 = platform->o_memPool.reserve<dfloat>(cds->fieldOffset[is]);
     if (platform->options.compareArgs("SCALAR" + sid + " INITIAL GUESS", "EXTRAPOLATION") && stage == 1)
       o_S0.copyFrom(cds->o_Se, cds->fieldOffset[is], 0, cds->fieldOffsetScan[is]);
     else
       o_S0.copyFrom(cds->o_S, cds->fieldOffset[is], 0, cds->fieldOffsetScan[is]);
     
     return o_S0;
  }();

  ellipticSolve(cds->solver[is], o_rhs, o_S);

  return o_S;
}


