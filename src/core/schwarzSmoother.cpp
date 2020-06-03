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

      std::vector<pfloat> casted_Sx(N2*Nelements);
      std::vector<pfloat> casted_Sy(N2*Nelements);
      std::vector<pfloat> casted_Sz(N2*Nelements);
      std::vector<pfloat> casted_D(N3*Nelements);
      std::vector<pfloat> casted_wts(N*N*4*3*Nelements);
      nek_schwarzOperators(Sx.data(), Sy.data(), Sz.data(), D.data(), wts.data(), pSolver->nLevels-level);
      for(dlong i = 0; i < N2*Nelements; ++i){
        casted_Sx[i] = static_cast<pfloat>(Sx[i]);
        casted_Sy[i] = static_cast<pfloat>(Sy[i]);
        casted_Sz[i] = static_cast<pfloat>(Sz[i]);
      }
      for(dlong i = 0; i < N3*Nelements; ++i){
        casted_D[i] = static_cast<pfloat>(D[i]);
      }
      for(dlong i = 0; i < N*N*4*3*Nelements; ++i){
        casted_wts[i] = static_cast<pfloat>(wts[i]);
      }

      currLevel->build(casted_Sx.data(), casted_Sy.data(), casted_Sz.data(), casted_D.data(), casted_wts.data());
    }
  }
}
