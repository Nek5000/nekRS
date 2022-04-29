#include <elliptic.h>

void
ellipticMultiGridUpdateLambda(elliptic_t* elliptic)
{
  mesh_t* mesh = elliptic->mesh;
  precon_t* precon = elliptic->precon;
  parAlmond::multigridLevel** levels = precon->parAlmond->levels;
  const int numMGLevels = elliptic->nLevels;
  for(int levelIndex = 0; levelIndex < numMGLevels; levelIndex++){
    auto mgLevel = dynamic_cast<MGLevel*>(levels[levelIndex]);

    if(levelIndex == 0){
      elliptic_t* ellipticFine = mgLevel->elliptic;
      ellipticFine->copyDfloatToPfloatKernel(2 * mesh->Nelements * mesh->Np,
        elliptic->o_lambda,
        ellipticFine->o_lambdaPfloat);
    }
    else {
      auto prevLevel = dynamic_cast<MGLevel*>(levels[levelIndex-1]);
      elliptic_t* ellipticFine = prevLevel->elliptic;
      elliptic_t* ellipticCoarse = mgLevel->elliptic;
      const int Nfq = ellipticFine->mesh->Nq;
      const int Ncq = ellipticCoarse->mesh->Nq;
      ellipticCoarse->copyPfloatToDPfloatKernel(2 * ellipticFine->mesh->Nelements * ellipticFine->mesh->Np,
        ellipticFine->o_lambdaPfloat,
        ellipticFine->o_lambda);

      ellipticCoarse->precon->coarsenKernel(2 * ellipticCoarse->mesh->Nelements, ellipticCoarse->o_interp, ellipticFine->o_lambda, ellipticCoarse->o_lambda);

      ellipticCoarse->copyDfloatToPfloatKernel(2 * ellipticCoarse->mesh->Nelements * ellipticCoarse->mesh->Np,
        ellipticCoarse->o_lambda,
        ellipticCoarse->o_lambdaPfloat);
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
