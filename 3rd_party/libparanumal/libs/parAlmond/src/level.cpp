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

multigridLevel::multigridLevel(dlong N, dlong M, KrylovType ktype_, MPI_Comm comm_):
  Nrows(N), Ncols(M), ktype(ktype_) {
  comm = comm_;
}

multigridLevel::~multigridLevel() {

  if (x  ) free(x  );
  if (rhs) free(rhs);
  if (res) free(res);

  if (ck) free(ck);
  if (vk) free(vk);
  if (wk) free(wk);

  if (weight) free(weight);

  if (o_x.size()) o_x.free();
  if (o_rhs.size()) o_rhs.free();
  if (o_res.size()) o_res.free();

  if (o_ck.size()) o_ck.free();
  if (o_vk.size()) o_vk.free();
  if (o_wk.size()) o_wk.free();

  if (o_weight.size()) o_weight.free();

}

void multigridLevel::kcycleOp1(dfloat *alpha1, dfloat *rho1,
                               dfloat *norm_rhs, dfloat *norm_rhstilde) {

  //ck = x
  memcpy(ck, x, Nrows*sizeof(dfloat));

  // vk = A*ck
  this->Ax(ck,vk);

  dfloat rho[3];

  if(ktype == PCG)
    kcycleCombinedOp1(Nrows, rho, ck, rhs, vk, weight, weighted, comm);

  if(ktype == GMRES)
    kcycleCombinedOp1(Nrows, rho, vk, rhs, vk, weight, weighted, comm);

  *alpha1 = rho[0];
  *rho1   = rho[1];
  *norm_rhs = sqrt(rho[2]);

  const dfloat a = -(*alpha1)/(*rho1);

  // rhs = rhs - (alpha1/rho1)*vk
  *norm_rhstilde = sqrt(vectorAddInnerProd(Nrows, a, vk, 1.0, rhs, o_weight, weighted,comm));
}

void multigridLevel::kcycleOp2(const dfloat alpha1, const dfloat rho1) {

  // w = A*x
  this->Ax(x,wk);

  dfloat rho[3];

  if(ktype == PCG)
    kcycleCombinedOp2(Nrows,rho, x, vk, wk, rhs, weight, weighted, comm);

  if(ktype == GMRES)
    kcycleCombinedOp2(Nrows,rho, wk, vk, wk, rhs, weight, weighted, comm);

  const dfloat gamma  = rho[0];
  const dfloat beta   = rho[1];
  const dfloat alpha2 = rho[2];

  if(fabs(rho1) > (dfloat) 1e-20){

    const dfloat rho2 = beta - gamma*gamma/rho1;

    if(fabs(rho2) > (dfloat) 1e-20){
      // x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ck + (alpha2/rho2)*dk
      const dfloat a = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
      const dfloat b = alpha2/rho2;

      vectorAdd(Nrows, a, ck, b, x);
    }
  }
}

void multigridLevel::device_kcycleOp1(dfloat *alpha1, dfloat *rho1,
                               dfloat *norm_rhs, dfloat *norm_rhstilde) {

  //ck = x
  o_ck.copyFrom(o_x, Nrows*sizeof(dfloat));

  // vk = A*ck
  this->Ax(o_ck,o_vk);

  dfloat rho[3];

  if(ktype == PCG)
    kcycleCombinedOp1(Nrows, rho, o_ck, o_rhs, o_vk, o_weight, weighted, comm);

  if(ktype == GMRES)
    kcycleCombinedOp1(Nrows, rho, o_vk, o_rhs, o_vk, o_weight, weighted, comm);

  *alpha1 = rho[0];
  *rho1   = rho[1];
  *norm_rhs = sqrt(rho[2]);

  const dfloat a = -(*alpha1)/(*rho1);

  // rhs = rhs - (alpha1/rho1)*vk
  *norm_rhstilde = sqrt(vectorAddInnerProd(Nrows, a, o_vk, 1.0, o_rhs, o_weight, weighted,comm));
}

void multigridLevel::device_kcycleOp2(const dfloat alpha1, const dfloat rho1) {

  // w = A*x
  this->Ax(o_x,o_wk);

  dfloat rho[3];

  if(ktype == PCG)
    kcycleCombinedOp2(Nrows,rho, o_x, o_vk, o_wk, o_rhs, o_weight, weighted, comm);

  if(ktype == GMRES)
    kcycleCombinedOp2(Nrows,rho, o_wk, o_vk, o_wk, o_rhs, o_weight, weighted, comm);

  const dfloat gamma  = rho[0];
  const dfloat beta   = rho[1];
  const dfloat alpha2 = rho[2];

  if(fabs(rho1) > (dfloat) 1e-20){

    const dfloat rho2 = beta - gamma*gamma/rho1;

    if(fabs(rho2) > (dfloat) 1e-20){
      // x = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ck + (alpha2/rho2)*dk
      const dfloat a = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
      const dfloat b = alpha2/rho2;

      vectorAdd(Nrows, a, o_ck, b, o_x);
    }
  }
}

}