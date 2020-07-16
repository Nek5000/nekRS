/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "elliptic.h"
#include <timer.hpp>

int ellipticSolve(elliptic_t *elliptic, dfloat tol,
                  occa::memory &o_r, occa::memory &o_x){
  
  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int Niter = 0;
  int maxIter = 1000; 
  
  if(elliptic->var_coeff && options.compareArgs("PRECONDITIONER", "JACOBI"))
    ellipticUpdateJacobi(elliptic);

  if(elliptic->allNeumann) {// zero mean of RHS
    ellipticZeroMean(elliptic, o_r);
  }
  
  options.getArgs("MAXIMUM ITERATIONS", maxIter);

  options.getArgs("SOLVER TOLERANCE", tol);

  if(options.compareArgs("RESIDUAL PROJECTION","ENABLED")){
    timer::tic("preSolveProjection");
    elliptic->residualProjection->preSolveProjection(o_r);
    timer::toc("preSolveProjection");
  }
  
  if(!options.compareArgs("KRYLOV SOLVER", "NONBLOCKING"))
    Niter = pcg (elliptic, o_r, o_x, tol, maxIter);
  else{
    if(!options.compareArgs("KRYLOV SOLVER", "FLEXIBLE")){
      Niter = nbpcg (elliptic, o_r, o_x, tol, maxIter);
    }
    else{      
      Niter = nbfpcg (elliptic, o_r, o_x, tol, maxIter);
    }
  }

  if(elliptic->allNeumann){ // zero mean of RHS
    ellipticZeroMean(elliptic, o_x);
  }

  if(options.compareArgs("RESIDUAL PROJECTION","ENABLED")){
    timer::tic("postSolveProjection");
    elliptic->residualProjection->postSolveProjection(o_x);
    timer::toc("postSolveProjection");
  }
  
  return Niter;

}
