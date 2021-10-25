#include <elliptic.h>

bool updateDiagonal(setupAide& options)
{
  return options.compareArgs("MULTIGRID DOWNWARD SMOOTHER","JACOBI") ||
         options.compareArgs("MULTIGRID UPWARD SMOOTHER","JACOBI") ||
         options.compareArgs("MULTIGRID SMOOTHER", "DAMPEDJACOBI");
}

void
ellipticUpdateLambda(elliptic_t* elliptic)
{
  const bool updateDiag = updateDiagonal(elliptic->options);
  mesh_t* mesh = elliptic->mesh;
  precon_t* precon = elliptic->precon;
  parAlmond::multigridLevel** levels = precon->parAlmond->levels;
  const int numMGLevels = elliptic->nLevels;
  for(int levelIndex = 0; levelIndex < numMGLevels; levelIndex++){
    auto mgLevel = dynamic_cast<MGLevel*>(levels[levelIndex]);

    if(updateDiag){
      const bool coarsestLevel = levelIndex == numMGLevels-1;
      if(!coarsestLevel || elliptic->options.compareArgs("MULTIGRID COARSE SOLVE", "FALSE"))
        ellipticUpdateJacobi(mgLevel->elliptic,mgLevel->o_invDiagA);
    }
    
    if(levelIndex == 0){
      elliptic_t* ellipticFine = mgLevel->elliptic;
      ellipticFine->copyDfloatToPfloatKernel(mesh->Nelements * mesh->Np,
        elliptic->o_lambda,
        ellipticFine->o_lambdaPfloat);
    }
    else {
      auto prevLevel = dynamic_cast<MGLevel*>(levels[levelIndex-1]);
      elliptic_t* ellipticFine = prevLevel->elliptic;
      elliptic_t* ellipticCoarse = mgLevel->elliptic;
      const int Nfq = ellipticFine->mesh->Nq;
      const int Ncq = ellipticCoarse->mesh->Nq;
      occa::memory o_lambdaFine = elliptic->o_p;
      ellipticCoarse->copyPfloatToDPfloatKernel(ellipticFine->mesh->Nelements * ellipticFine->mesh->Np,
        ellipticFine->o_lambdaPfloat,
        o_lambdaFine);

      ellipticCoarse->precon->coarsenKernel(ellipticCoarse->mesh->Nelements, mgLevel->o_interp, o_lambdaFine, ellipticCoarse->o_lambda);

      ellipticCoarse->copyDfloatToPfloatKernel(ellipticCoarse->mesh->Nelements * ellipticCoarse->mesh->Np,
        ellipticCoarse->o_lambda,
        ellipticCoarse->o_lambdaPfloat);
    }
  }
}