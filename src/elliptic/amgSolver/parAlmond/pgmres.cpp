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

void gmresUpdate(dlong Nrows,
                 dfloat *x,
                 dfloat **V,
                 dfloat *H,
                 dfloat *s,
                 int Niter,
                 int maxIt){

  dfloat *y = (dfloat *) calloc(Niter, sizeof(dfloat));

  for(int k=Niter-1; k>=0; --k){
    y[k] = s[k];

    for(int m=k+1; m<maxIt; ++m)
      y[k] -= H[k + m*(maxIt+1)]*y[m];

    y[k] /= H[k + k*(maxIt+1)];
  }

  for(int j=0; j<Niter; ++j){
    vectorAdd(Nrows, y[j], V[j], 1.0, x);
  }

  free(y);
}

void gmresUpdate(dlong Nrows,
                 occa::memory o_x,
                 occa::memory *o_V,
                 dfloat *H,
                 dfloat *s,
                 int Niter,
                 int maxIt){

  dfloat *y = (dfloat *) calloc(Niter, sizeof(dfloat));

  for(int k=Niter-1; k>=0; --k){
    y[k] = s[k];

    for(int m=k+1; m<Niter; ++m)
      y[k] -= H[k + m*(maxIt+1)]*y[m];

    y[k] /= H[k + k*(maxIt+1)];
  }

  for(int j=0; j<Niter; ++j){
    vectorAdd(Nrows, y[j], o_V[j], 1.0, o_x);
  }

  free(y);
}

void solver_t::pgmres(const int maxIt,
                      const dfloat tol){

  const dlong m = levels[0]->Nrows;
  const dlong n = levels[0]->Ncols;

  ktype = GMRES;

  // use parAlmond's buffers
  dfloat *r = levels[0]->rhs;
  dfloat *z = levels[0]->x;

  // initial residual
  dfloat nb = sqrt(vectorInnerProd(m, r, r, levels[0]->comm));

  //    x = 0;
  dfloat *x  = (dfloat *) calloc(n,sizeof(dfloat));
  vectorSet(m, 0.0, x);

  //sanity check
  if (nb<=tol) {
    memcpy(levels[0]->x, x, m*sizeof(dfloat));
    free(x);
    return;
  }

  // M r = b - A*x0
  if(ctype==KCYCLE) {
    this->kcycle(0);
  } else if(ctype==VCYCLE) {
    this->vcycle(0);
  }
  memcpy(r, z, m*sizeof(dfloat));

  dfloat nr = sqrt(vectorInnerProd(m, r, r, levels[0]->comm));

  dfloat *s = (dfloat *) calloc(maxIt+1, sizeof(dfloat));
  s[0] = nr;

  dfloat **V = (dfloat **) calloc(maxIt,sizeof(dfloat *));

  for(int i=0; i<maxIt; ++i){
    V[i] = (dfloat *) calloc(m, sizeof(dfloat)); //TODO this is way too much memory if maxit is large
  }

  // V(:,0) = r/nr
  vectorAdd(m, (1./nr), r,  0., V[0]);

  dfloat *H = (dfloat*) calloc((maxIt+1)*(maxIt+1), sizeof(dfloat));
  dfloat *J = (dfloat*) calloc(4*maxIt, sizeof(dfloat));

  dfloat *Av = (dfloat *) calloc(m, sizeof(dfloat));
  dfloat *w  = (dfloat *) calloc(m, sizeof(dfloat));

  int Niter=0;

  for(int i=0; i<maxIt; i++){

    Niter = i+1;
    // Av = A*V(:.i)
    levels[0]->Ax(V[i], Av);

    // M w = A vi
    memcpy(r, Av, m*sizeof(dfloat));

    if(ctype==KCYCLE) {
      this->kcycle(0);
    } else if(ctype==VCYCLE) {
      this->vcycle(0);
    }
    memcpy(w, z, m*sizeof(dfloat));

    for(int k=0; k<=i; ++k){
      dfloat hki = vectorInnerProd(m, w, V[k], levels[0]->comm);

      // w = w - hki*V[k]
      vectorAdd(m, -hki, V[k], 1.0, w);

      // H(k,i) = hki
      H[k + i*(maxIt+1)] = hki;
    }

    dfloat wdotw = vectorInnerProd(m, w, w, levels[0]->comm);

    H[i+1 + i*(maxIt+1)] = sqrt(wdotw);

    for(int k=0; k<i; ++k){
      dfloat h1 = H[k +     i*(maxIt+1)];
      dfloat h2 = H[k + 1 + i*(maxIt+1)];

      H[k +     i*(maxIt+1)] = J[4*k    ]*h1 + J[4*k + 2]*h2;
      H[k + 1 + i*(maxIt+1)] = J[4*k + 1]*h1 + J[4*k + 3]*h2;
    }

    dfloat h1 = H[i + i*(maxIt+1)];
    dfloat h2 = H[i + 1 + i*(maxIt+1)];
    dfloat hr = sqrt(h1*h1 + h2*h2);

    H[i   +  i*(maxIt+1)] = hr;
    H[i+1 +  i*(maxIt+1)] = 0.;

    dfloat ct = h1/hr;
    dfloat st = h2/hr;
    J[4*i    ] =  ct;     J[4*i + 2] = st;
    J[4*i + 1] = -st;     J[4*i + 3] = ct;

    dfloat s1 = s[i];
    dfloat s2 = s[i+1];

    s[i  ] =  ct*s1 + st*s2;
    s[i+1] = -st*s1 + ct*s2;

    if(fabs(s[i+1]) < tol) break;

    if(i < maxIt-1){
      dfloat wdotw = vectorInnerProd(m, w, w, levels[0]->comm);
      dfloat nw = sqrt(wdotw);

      // V(:,i+1) = w/nw
      vectorAdd(m,1./nw, w, 0.0, V[i+1]);
    }
  }

  gmresUpdate(m, x, V, H, s, Niter, maxIt);

  //copy result back to parAlmond's x storage
  memcpy(levels[0]->x, x, m*sizeof(dfloat));

  free(x);
  free(s); free(V);
  free(H); free(J);
  free(Av); free(w);

  if(Niter == maxIt)
    printf("gmres did not converge in given number of iterations\n");
}

void solver_t::device_pgmres(const int maxIt,
                             const dfloat tol){

  const dlong m = levels[0]->Nrows;
  const dlong n = levels[0]->Ncols;

  // use parAlmond's buffers
  occa::memory &o_r = levels[0]->o_rhs;
  occa::memory &o_z = levels[0]->o_x;

  // initial residual
  dfloat nb = sqrt(vectorInnerProd(m, o_r, o_r, levels[0]->comm));

  occa::memory  o_x = device.malloc(n*sizeof(dfloat), levels[0]->x);
  occa::memory  o_Av= device.malloc(n*sizeof(dfloat), levels[0]->x);
  occa::memory  o_w = device.malloc(n*sizeof(dfloat), levels[0]->x);

  //sanity check
  if (nb<=tol) {
    levels[0]->o_x.copyFrom(o_x);
    printf("Almond PGMRES iter %d, res = %g\n", 0, nb);
    o_x.free(); o_Av.free(); o_w.free();
    return;
  }

  // M r = b - A*x0
  if(ctype==KCYCLE) {
    this->device_kcycle(0);
  } else if(ctype==VCYCLE) {
    this->device_vcycle(0);
  }
  o_r.copyFrom(o_z);


  dfloat nr = sqrt(vectorInnerProd(m, o_r, o_r, levels[0]->comm));

  dfloat *s = (dfloat *) calloc(maxIt+1, sizeof(dfloat));
  s[0] = nr;

  occa::memory *o_V = (occa::memory *) calloc(maxIt, sizeof(occa::memory));
  for(int i=0; i<maxIt; ++i){
    o_V[i] = device.malloc(n*sizeof(dfloat), levels[0]->x);
  }

  // V(:,0) = r/nr
  vectorAdd(m, (1./nr), o_r, 0., o_V[0]);

  dfloat *H = (dfloat *) calloc((maxIt+1)*(maxIt+1), sizeof(dfloat));
  dfloat *J = (dfloat *) calloc(4*maxIt, sizeof(dfloat));

  int Niter = 0;

  int i;
  for(i=0; i<maxIt; i++){

    Niter = i+1;

    // r = A*V(:.i)
    levels[0]->Ax(o_V[i], o_r);

    // M w = A vi
    if(ctype==KCYCLE) {
      this->device_kcycle(0);
    } else if(ctype==VCYCLE) {
      this->device_vcycle(0);
    }

    for(int k=0; k<=i; ++k){
      dfloat hki = vectorInnerProd(m, o_z, o_V[k], levels[0]->comm);

      // w = w - hki*V[k]
      vectorAdd(m, -hki, o_V[k], 1.0, o_z);

      // H(k,i) = hki
      H[k + i*(maxIt+1)] = hki;
    }

    dfloat nw = sqrt(vectorInnerProd(m, o_z, o_z, levels[0]->comm));
    H[i+1 + i*(maxIt+1)] = nw;

    for(int k=0; k<i; ++k){
      dfloat h1 = H[k +     i*(maxIt+1)];
      dfloat h2 = H[k + 1 + i*(maxIt+1)];

      H[k +     i*(maxIt+1)] = J[4*k    ]*h1 + J[4*k + 2]*h2;
      H[k + 1 + i*(maxIt+1)] = J[4*k + 1]*h1 + J[4*k + 3]*h2;
    }

    dfloat h1 = H[i + i*(maxIt+1)];
    dfloat h2 = H[i + 1 + i*(maxIt+1)];
    dfloat hr = sqrt(h1*h1 + h2*h2);

    H[i   +  i*(maxIt+1)] = hr;
    H[i+1 +  i*(maxIt+1)] = 0;

    dfloat ct = h1/hr;
    dfloat st = h2/hr;
    J[4*i    ] =  ct;     J[4*i + 2] = st;
    J[4*i + 1] = -st;     J[4*i + 3] = ct;

    dfloat s1 = s[i];
    dfloat s2 = s[i+1];

    s[i  ] =  ct*s1 + st*s2;
    s[i+1] = -st*s1 + ct*s2;

    if(fabs(s[i+1]) < tol) break;

    if(i < maxIt-1){
      // V(:,i+1) = w/nw
      vectorAdd(m, 1./nw, o_z, 0.0, o_V[i+1]);
    }
  }

  gmresUpdate(m, o_x, o_V, H, s, Niter, maxIt);

  //copy result back to parAlmond's x storage
  levels[0]->o_x.copyFrom(o_x);

  printf("Almond PGMRES iter %d, res = %g\n", Niter, fabs(s[i+1]));

  if(Niter == maxIt)
    printf("gmres did not converge in given number of iterations \n");

  for(int i=0; i<maxIt; ++i)
    o_V[i].free();
  free((void*)o_V);

  free(s);
  free(H); free(J);

  o_Av.free();
  o_w.free();
  o_x.free();
}

} //namepsace parAlmond