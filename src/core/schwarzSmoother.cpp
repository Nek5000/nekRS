#include "schwarzSmoother.hpp"
#include <util.hpp>
#include "nekInterfaceAdapter.hpp"
#include <vector>
void reconfigurePressureSolver(elliptic_t* pSolver){
  auto Nelements = mesh->Nelements;
  auto precon = pSolver->precon;
  auto levels = precon->parAlmond->levels;
  for(int level = 0; level < precon->parAlmond->numLevels; ++level){
    MGLevel* curr_level = dynamic_cast<MGLevel*>(levels[level]);
    if(curr_level){
      auto elliptic = curr_level->elliptic;
      auto N = elliptic->mesh->Nq;
      auto N2 = (N+2)*(N+2);
      auto N3 = N2*(N+2);
      std::vector<dfloat> Sx(N2*Nelements);
      std::vector<dfloat> Sy(N2*Nelements);
      std::vector<dfloat> Sz(N2*Nelements);
      std::vector<dfloat> D(N3*Nelements);
      std::vector<dfloat> wts(N*N*4*3*Nelements);
      get_nek_operators(Sx.data(), Sy.data(), Sz.data(), D.data(), wts.data(), level);
      curr_level->fdm_op->build(Sx.data(), Sy.data(), Sz.data(), D.data(), wts.data());
    }
  }
}
}
