#include <elliptic.h>
#include "ellipticPrecon.h"
#include "ellipticMultiGrid.h"

void
ellipticMultiGridUpdateLambda(elliptic_t* elliptic)
{
  precon_t *precon = elliptic->precon;
  MGSolver_t::multigridLevel** levels = precon->MGSolver->levels;
  {
    auto level = dynamic_cast<pMGLevel*>(levels[0]);
    elliptic_t* ellipticFine = level->elliptic;
    platform->copyDfloatToPfloatKernel(elliptic->mesh->Nlocal, elliptic->o_lambda0, ellipticFine->o_lambda0);
    if(!ellipticFine->poisson)
      platform->copyDfloatToPfloatKernel(elliptic->mesh->Nlocal, elliptic->o_lambda1, ellipticFine->o_lambda1);
  }

  for(int levelIndex = 1; levelIndex < elliptic->nLevels; levelIndex++) {
    auto mgLevel = dynamic_cast<pMGLevel*>(levels[levelIndex]);
    elliptic_t* ellipticCoarse = mgLevel->elliptic;

    auto prevLevel = dynamic_cast<pMGLevel*>(levels[levelIndex-1]);
    elliptic_t* ellipticFine = prevLevel->elliptic;

    precon_t *precon = ellipticCoarse->precon;
    precon->coarsenKernel(ellipticCoarse->mesh->Nelements, 
                          ellipticCoarse->o_interp, 
                          ellipticFine->o_lambda0, ellipticCoarse->o_lambda0);
    if(!elliptic->poisson)
      precon->coarsenKernel(ellipticCoarse->mesh->Nelements, 
                            ellipticCoarse->o_interp, 
                            ellipticFine->o_lambda1, ellipticCoarse->o_lambda1);
 
  }
}
