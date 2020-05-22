#include "schwarzSmoother.hpp"
#include "nekInterfaceAdapter.hpp"
#include <vector>
void reconfigurePressureSolver(elliptic_t* pSolver){
  const dlong Nelements = pSolver->mesh->Nelements;
  parAlmond::multigridLevel** mglevels = pSolver->precon->parAlmond->levels;
  for(int level = 0; level < pSolver->nLevels; ++level){
    MGLevel* currLevel = dynamic_cast<MGLevel*>(mglevels[level]);
    if(currLevel){
      const int N = pSolver->levels[level]+1;
      const int N2 = (N+2)*(N+2);
      const int N3 = N2*(N+2);
      std::vector<dfloat> Sx(N2*Nelements);
      std::vector<dfloat> Sy(N2*Nelements);
      std::vector<dfloat> Sz(N2*Nelements);
      std::vector<dfloat> D(N3*Nelements);
      std::vector<dfloat> wts(N*N*4*3*Nelements);
      get_nek_operators(Sx.data(), Sy.data(), Sz.data(), D.data(), wts.data(), 3-level);
      currLevel->build(Sx.data(), Sy.data(), Sz.data(), D.data(), wts.data());
    }
  }
}