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




//
extern "C" void ellipticAxVarHex3D(const dlong & Nelements,
                        const dlong & offset,
                        const dfloat* __restrict__ ggeo,
                        const dfloat* __restrict__ D,
                        const dfloat* __restrict__ S,
                        const dfloat* __restrict__ lambda,
                        const dfloat* __restrict__ q,
                        dfloat* __restrict__ Aq )
{
  dfloat s_q[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqr[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqs[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqt[p_Nq][p_Nq][p_Nq];

  #pragma omp parallel for private(s_q, s_Gqr, s_Gqs, s_Gqt)
  for(dlong e = 0; e < Nelements; ++e) {
    const dlong element = e;

#pragma unroll
    for(int k = 0; k < p_Nq; k++)
#pragma unroll
      for(int j = 0; j < p_Nq; ++j)
#pragma unroll
        for(int i = 0; i < p_Nq; ++i) {
          const dlong base = i + j * p_Nq + k * p_Nq * p_Nq + element * p_Np;
          const dfloat qbase = q[base];
          s_q[k][j][i] = qbase;
        }

#pragma unroll
    for(int k = 0; k < p_Nq; ++k)
#pragma unroll
      for(int j = 0; j < p_Nq; ++j)
#pragma unroll
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_G00 = ggeo[gbase + p_G00ID * p_Np];
          const dfloat r_G01 = ggeo[gbase + p_G01ID * p_Np];
          const dfloat r_G11 = ggeo[gbase + p_G11ID * p_Np];
          const dfloat r_G12 = ggeo[gbase + p_G12ID * p_Np];
          const dfloat r_G02 = ggeo[gbase + p_G02ID * p_Np];
          const dfloat r_G22 = ggeo[gbase + p_G22ID * p_Np];

          const dlong id = element * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_lam0 = lambda[id + 0 * offset];

          dfloat qr = 0.f;
          dfloat qs = 0.f;
          dfloat qt = 0.f;

#pragma unroll
          for(int m = 0; m < p_Nq; m++){
            qr += S[m*p_Nq + i] * s_q[k][j][m];
            qs += S[m*p_Nq + j] * s_q[k][m][i];
            qt += S[m*p_Nq + k] * s_q[m][j][i];
          }

          dfloat Gqr = r_G00 * qr;
          Gqr += r_G01 * qs;
          Gqr += r_G02 * qt;

          dfloat Gqs = r_G01 * qr;
          Gqs += r_G11 * qs;
          Gqs += r_G12 * qt;

          dfloat Gqt = r_G02 * qr;
          Gqt += r_G12 * qs;
          Gqt += r_G22 * qt;

          s_Gqr[k][j][i] = r_lam0 * Gqr;
          s_Gqs[k][j][i] = r_lam0 * Gqs;
          s_Gqt[k][j][i] = r_lam0 * Gqt;
        }

#pragma unroll
    for(int k = 0; k < p_Nq; k++)
#pragma unroll
      for(int j = 0; j < p_Nq; ++j)
#pragma unroll
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_GwJ = ggeo[gbase + p_GWJID * p_Np];

          const dlong id = element * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_lam1 = lambda[id + 1 * offset];

          dfloat r_Aq = r_GwJ * r_lam1 * s_q[k][j][i];
          dfloat r_Aqr = 0, r_Aqs = 0, r_Aqt = 0;

#pragma unroll
          for(int m = 0; m < p_Nq; m++){
            r_Aqr += D[m*p_Nq+i] * s_Gqr[k][j][m];
            r_Aqs += D[m*p_Nq+j] * s_Gqs[k][m][i];
            r_Aqt += D[m*p_Nq+k] * s_Gqt[m][j][i];
          }

          Aq[id] = r_Aqr + r_Aqs + r_Aqt + r_Aq;
        }
  }
}
