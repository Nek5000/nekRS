/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"

namespace parAlmond {

void solver_t::pcg(const int maxIt, const dfloat tol){

  const dlong m = levels[0]->Nrows;
  const dlong n = levels[0]->Ncols;

  ktype = PCG;

  // use parAlmond's buffers
  dfloat *r = levels[0]->rhs;
  dfloat *z = levels[0]->x;

  // initial residual
  dfloat rdotr0 = vectorInnerProd(m, r, r, levels[0]->comm);

  dfloat *x  = (dfloat *) calloc(n,sizeof(dfloat));
  dfloat *Ap = (dfloat *) calloc(n,sizeof(dfloat));
  dfloat *p  = (dfloat *) calloc(n,sizeof(dfloat));

  //sanity check
  if (rdotr0<=(tol*tol)) {
    memcpy(levels[0]->x, x, m*sizeof(dfloat));
    free(x); free(p); free(Ap);
    return;
  }

  // Precondition, z = M^{-1}*r
  if(ctype==KCYCLE) {
    this->kcycle(0);
  } else if(ctype==VCYCLE) {
    this->vcycle(0);
  }
  memcpy(p, z, m*sizeof(dfloat));

  dfloat rdotz0 = vectorInnerProd(m, r, z, levels[0]->comm);

  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  dfloat alpha, beta, pAp;

  int Niter = 0;
  while(rdotr0>(tol*tol)){
    //   Ap = A*p;
    levels[0]->Ax(p, Ap);

    dfloat pAp = vectorInnerProd(m, p, Ap, levels[0]->comm);

    alpha = rdotz0/pAp;

    // update solution
    //    x = x + alpha * p;
    vectorAdd(m, alpha, p, 1.0, x);

    // update residual
    // r = r - alpha * Ap;
    vectorAdd(m, -alpha, Ap, 1.0, r);

    dfloat rdotr1 = vectorInnerProd(m, r, r, levels[0]->comm);

    if(rdotr1 < tol*tol) {
      rdotr0 = rdotr1;
      break;
    }

    // Precondition, z = M^{-1}*r
    if(ctype==KCYCLE) {
      this->kcycle(0);
    } else if(ctype==VCYCLE) {
      this->vcycle(0);
    }

    dfloat rdotz1 = vectorInnerProd(m, r, z, levels[0]->comm);

    if(ctype==KCYCLE) {
      // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
      dfloat zdotAp = vectorInnerProd(m, z, Ap, levels[0]->comm);
      beta = -alpha*zdotAp/rdotz0;
    } else {
      beta = rdotz1/rdotz0;
    }

    // p = z + beta*p
    vectorAdd(m, 1.0, z, beta, p);

    // switch rdotz0 <= rdotz1
    rdotz0 = rdotz1;

    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;

    Niter++;

    printf("Almond PCG iter %d, res = %g\n", Niter, sqrt(rdotr0));

    if(Niter==maxIt) break;
  }

  //copy result back to parAlmond's x storage
  memcpy(levels[0]->x, x, m*sizeof(dfloat));
  free(x); free(p); free(Ap);
}

void solver_t::device_pcg(const int maxIt, const dfloat tol){

  const dlong m = levels[0]->Nrows;
  const dlong n = levels[0]->Ncols;

  ktype = PCG;

  // use parAlmond's buffers
  occa::memory &o_r = levels[0]->o_rhs;
  occa::memory &o_z = levels[0]->o_x;

  // initial residual
  dfloat rdotr0 = vectorInnerProd(m, o_r, o_r, levels[0]->comm);

  //  dfloat TOL =  mymax(tol*tol*normB,tol*tol);
  
  occa::memory o_x  = device.malloc(n*sizeof(dfloat),levels[0]->x);
  occa::memory o_Ap = device.malloc(n*sizeof(dfloat),levels[0]->x);
  occa::memory o_p  = device.malloc(n*sizeof(dfloat),levels[0]->x);

  //    x = 0;
  vectorSet(m, 0.0, o_x);

  //sanity check
  if (rdotr0<=(tol*tol)) {
    levels[0]->o_x.copyFrom(o_x);
    printf("Almond PCG iter %d, res = %g\n", 0, sqrt(rdotr0));
    o_x.free(); o_p.free(); o_Ap.free();
    return;
  }

  // Precondition, z = M^{-1}*r
  if(ctype==KCYCLE) {
    this->device_kcycle(0);
  } else if(ctype==VCYCLE) {
    this->device_vcycle(0);
  }
  o_p.copyFrom(o_z);

  dfloat rdotz0 = vectorInnerProd(m, o_r, o_z, levels[0]->comm);

  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  dfloat alpha, beta, pAp;

  int Niter = 0;
  while(rdotr0>(tol*tol)){
    //   Ap = A*p;
    levels[0]->Ax(o_p, o_Ap);

    dfloat pAp = vectorInnerProd(m, o_p, o_Ap, levels[0]->comm);

    alpha = rdotz0/pAp;

    // update solution
    //    x = x + alpha * p;
    vectorAdd(m, alpha, o_p, 1.0, o_x);

    // update residual
    // r = r - alpha * Ap;
    vectorAdd(m, -alpha, o_Ap, 1.0, o_r);

    dfloat rdotr1 = vectorInnerProd(m, o_r, o_r, levels[0]->comm);

    if(rdotr1 < tol*tol) {
      rdotr0 = rdotr1;
      break;
    }

    // Precondition, z = M^{-1}*r
    if(ctype==KCYCLE) {
      this->device_kcycle(0);
    } else if(ctype==VCYCLE) {
      this->device_vcycle(0);
    }

    dfloat rdotz1 = vectorInnerProd(m, o_r, o_z, levels[0]->comm);

    if(ctype==KCYCLE) {
      // flexible pcg beta = (z.(-alpha*Ap))/zdotz0
      dfloat zdotAp = vectorInnerProd(m, o_z, o_Ap, levels[0]->comm);
      beta = -alpha*zdotAp/rdotz0;
    } else if(ctype==VCYCLE) {
      beta = rdotz1/rdotz0;
    }

    // p = z + beta*p
    vectorAdd(m, 1.0, o_z, beta, o_p);

    // switch rdotz0 <= rdotz1
    rdotz0 = rdotz1;

    // switch rdotz0,rdotr0 <= rdotz1,rdotr1
    rdotr0 = rdotr1;

    Niter++;

    //printf("Almond PCG iter %d, res = %g\n", Niter, sqrt(rdotr0));

    if(Niter==maxIt) break;
  }

  //copy result back to parAlmond's x storage
  levels[0]->o_x.copyFrom(o_x);

  printf("Almond PCG iter %d, res = %g\n", Niter, sqrt(rdotr0));

  o_x.free(); o_p.free(); o_Ap.free();
}

} //namespace parAlmond
