#include <elliptic.h>
#include "ellipticPrecon.h"
#include "ellipticMultiGrid.h"

void
ellipticMultiGridUpdateLambda(elliptic_t* elliptic)
{
  mesh_t* mesh = elliptic->mesh;
  precon_t* precon = (precon_t*) elliptic->precon;
  MGSolver_t::multigridLevel** levels = precon->MGSolver->levels;
  const int numMGLevels = elliptic->nLevels;
  for(int levelIndex = 0; levelIndex < numMGLevels; levelIndex++){
    auto mgLevel = dynamic_cast<pMGLevel*>(levels[levelIndex]);

    if(levelIndex > 0){
      auto prevLevel = dynamic_cast<pMGLevel*>(levels[levelIndex-1]);
      elliptic_t* ellipticFine = prevLevel->elliptic;
      elliptic_t* ellipticCoarse = mgLevel->elliptic;
      const int Nfq = ellipticFine->mesh->Nq;
      const int Ncq = ellipticCoarse->mesh->Nq;

      precon_t* precon = (precon_t*) ellipticCoarse->precon;
      precon->coarsenKernel(ellipticCoarse->mesh->Nelements, ellipticCoarse->o_interp, 
                            ellipticFine->o_lambda, ellipticCoarse->o_lambda);
    }

    if(elliptic->options.compareArgs("MULTIGRID DOWNWARD SMOOTHER","JACOBI") ||
       elliptic->options.compareArgs("MULTIGRID UPWARD SMOOTHER","JACOBI") ||
       elliptic->options.compareArgs("MULTIGRID SMOOTHER", "DAMPEDJACOBI"))
    {
      const bool coarsestLevel = levelIndex == numMGLevels-1;
      if(!coarsestLevel || elliptic->options.compareArgs("MULTIGRID COARSE SOLVE", "FALSE"))
        ellipticUpdateJacobi(mgLevel->elliptic,mgLevel->o_invDiagA);
    }
    
  }
}
