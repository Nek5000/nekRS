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

extern "C" void FUNC(ellipticPartialAxHex3D)(const dlong & Nelements,
                     const dlong & offset,
                     const dlong & loffset,
                     const dlong* __restrict__ elementList,
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

  dfloat s_D[p_Nq][p_Nq];
  dfloat s_S[p_Nq][p_Nq];

  for(int j = 0; j < p_Nq; ++j)
    for(int i = 0; i < p_Nq; ++i) {
      s_D[j][i] = D[j * p_Nq + i];
      s_S[j][i] = S[j * p_Nq + i];
    }

#ifdef __NEKRS__OMP__
  #pragma omp parallel for private(s_q, s_Gqr, s_Gqs, s_Gqt)
#endif
  for(dlong e = 0; e < Nelements; ++e) {
    const dlong element = elementList[e];

    for(int k = 0; k < p_Nq; k++)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong base = i + j * p_Nq + k * p_Nq * p_Nq + element * p_Np;
          s_q[k][j][i] = q[base];
        }

    for(int k = 0; k < p_Nq; ++k)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_G00 = ggeo[gbase + p_G00ID * p_Np];
          const dfloat r_G01 = ggeo[gbase + p_G01ID * p_Np];
          const dfloat r_G11 = ggeo[gbase + p_G11ID * p_Np];
          const dfloat r_G12 = ggeo[gbase + p_G12ID * p_Np];
          const dfloat r_G02 = ggeo[gbase + p_G02ID * p_Np];
          const dfloat r_G22 = ggeo[gbase + p_G22ID * p_Np];

          dfloat qr = 0.f;
          dfloat qs = 0.f;
          dfloat qt = 0.f;

          for(int m = 0; m < p_Nq; m++) {
            qr += s_D[i][m] * s_q[k][j][m];
            qs += s_D[j][m] * s_q[k][m][i];
            qt += s_D[k][m] * s_q[m][j][i];
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

          s_Gqr[k][j][i] = Gqr;
          s_Gqs[k][j][i] = Gqs;
          s_Gqt[k][j][i] = Gqt;
        }

    for(int k = 0; k < p_Nq; k++)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;

          dfloat r_Aq = 0;
#ifndef p_poisson
          r_Aq = ggeo[gbase + p_GWJID * p_Np] * lambda[1*loffset] * s_q[k][j][i];
#endif
          dfloat r_Aqr = 0, r_Aqs = 0, r_Aqt = 0;

          for(int m = 0; m < p_Nq; m++)
            r_Aqr += s_S[i][m] * s_Gqr[k][j][m];
          for(int m = 0; m < p_Nq; m++)
            r_Aqs += s_S[j][m] * s_Gqs[k][m][i];
          for(int m = 0; m < p_Nq; m++)
            r_Aqt += s_S[k][m] * s_Gqt[m][j][i];

          const dlong id = element * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          Aq[id] = lambda[0*loffset]*(r_Aqr + r_Aqs + r_Aqt) + r_Aq;
        }
  }
}
