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

int pcg(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_x,
        const dfloat tol, const int MAXIT)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int fixedIterationCountFlag = 0;
  int enableGatherScatters = 1;
  int flexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");
  int verbose = options.compareArgs("VERBOSE", "TRUE");

  if(options.compareArgs("FIXED ITERATION COUNT", "TRUE"))
    fixedIterationCountFlag = 1;

  // register scalars
  dfloat rdotz1 = 0;
  dfloat rdotz2 = 0;

  // now initialized
  dfloat alpha = 0, beta = 0, pAp = 0;

  /*aux variables */
  occa::memory &o_p  = elliptic->o_p;
  occa::memory &o_z  = elliptic->o_z;
  occa::memory &o_Ap = elliptic->o_Ap;
  occa::memory &o_weight = elliptic->o_invDegree;

  elliptic->fillKernel(elliptic->Nfields*elliptic->Ntotal, 0.0, o_p);

  pAp = 0;
  rdotz1 = 1;

  const dfloat rdotr0 = ellipticWeightedNorm2(elliptic, o_weight, o_r) * elliptic->resNormFactor;
  const dfloat TOL =  tol * tol;
  if (verbose && (mesh->rank == 0))
    printf("CG: initial res norm %12.12f WE NEED TO GET TO %12.12f \n",
           sqrt(rdotr0), sqrt(TOL));

  //if (rdotr0 <= TOL && !fixedIterationCountFlag) return 0;
  if (rdotr0 == 0) return 0;

  int iter;
  for(iter = 1; iter <= MAXIT; ++iter) {
    // z = Precon^{-1} r
    ellipticPreconditioner(elliptic, o_r, o_z);

    rdotz2 = rdotz1;

    // r.z
    rdotz1 = ellipticWeightedInnerProduct(elliptic, o_weight, o_r, o_z);
    if(flexible) {
      dfloat zdotAp = ellipticWeightedInnerProduct(elliptic, o_weight, o_z, o_Ap);
      beta = -alpha * zdotAp / rdotz2;
    } else {
      beta = (iter == 1) ? 0:rdotz1 / rdotz2;
    }

    // p = z + beta*p
    ellipticScaledAdd(elliptic, 1.f, o_z, beta, o_p);

    // A*p
    ellipticOperator(elliptic, o_p, o_Ap, dfloatString);

    // dot(p,A*p)
    pAp =  ellipticWeightedInnerProduct(elliptic, o_weight, o_p, o_Ap);
    alpha = rdotz1 / pAp;

    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)
    dfloat rdotr =
      ellipticUpdatePCG(elliptic, o_p, o_Ap, alpha, o_x, o_r) * elliptic->resNormFactor;

    if (verbose && (mesh->rank == 0)) {
      if(rdotr < 0)
        printf("WARNING CG: rdotr = %17.15lf\n", rdotr);

      printf("CG: it %d r norm %12.12le alpha = %le \n", iter, sqrt(rdotr), alpha);
    }

    if(rdotr <= TOL && !fixedIterationCountFlag) break;
  }
  return iter;
}
