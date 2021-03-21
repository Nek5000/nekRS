#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <occa.hpp>

using namespace std;
using namespace occa;
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

extern "C" void subCycleStrongCubatureVolumeHex3D(const int & Nelements,
                                                  const int * __restrict__ elementList,
                                                  const dfloat * __restrict__ cubD,
                                                  const dfloat * __restrict__ cubInterpT,
                                                  const dfloat * __restrict__ cubProjectT,
                                                  const int & offset,
                                                  const int & cubatureOffset,
                                                  const int & NUoffset,
                                                  const dfloat * __restrict__ invLumpedMassMatrix,
                                                  const dfloat * __restrict__ BdivW,
                                                  const dfloat & c0,
                                                  const dfloat & c1,
                                                  const dfloat & c2,
                                                  const dfloat * __restrict__ conv,
                                                  const dfloat * __restrict__ Ud,
                                                  dfloat * __restrict__ NU) {
  // (phi, U.grad Ud)
  dfloat r_c[3] = {c0, c1,c2};
  dfloat s_cubD[p_cubNq][p_cubNq];
  dfloat s_cubInterpT[p_Nq][p_cubNq];
  dfloat s_cubProjectT[p_cubNq][p_Nq];
  dfloat s_U[p_cubNq][p_cubNq];
  dfloat s_V[p_cubNq][p_cubNq];
  dfloat s_W[p_cubNq][p_cubNq];
  dfloat s_Ud[p_cubNq][p_cubNq];
  dfloat s_Vd[p_cubNq][p_cubNq];
  dfloat s_Wd[p_cubNq][p_cubNq];
  dfloat s_Ud1[p_Nq][p_cubNq];
  dfloat s_Vd1[p_Nq][p_cubNq];
  dfloat s_Wd1[p_Nq][p_cubNq];
  dfloat r_U2[256][p_cubNq], r_V2[256][p_cubNq], r_W2[256][p_cubNq];
  dfloat r_Ud[256][p_cubNq], r_Vd[256][p_cubNq], r_Wd[256][p_cubNq];
  for (int e = 0; e < Nelements; ++e) {
    int _occa_exclusive_index = 0;
    const int element = elementList[e];
    #pragma unroll
    for (int j = 0; j < p_cubNq; ++j) {
      #pragma unroll
      for (int i = 0; i < p_cubNq; ++i) {
        const int id = i + j * p_cubNq;
        if (id < p_Nq * p_cubNq) {
          s_cubInterpT[j][i] = cubInterpT[id];
          s_cubProjectT[j][i] = cubProjectT[id];
        }
        s_cubD[j][i] = cubD[id];
        #pragma unroll
        for (int k = 0; k < p_cubNq; ++k) {
          r_Ud[_occa_exclusive_index][k] = 0;
          r_Vd[_occa_exclusive_index][k] = 0;
          r_Wd[_occa_exclusive_index][k] = 0;
        }
        ++_occa_exclusive_index;
      }
    }
    #pragma unroll
    for (int c = 0; c < p_Nq; ++c) {
      #pragma unroll
      for (int b = 0; b < p_Nq; ++b) {
        #pragma unroll
        for (int a = 0; a < p_Nq; ++a) {
            // this can be improved
            const int id = element * p_Np + c * p_Nq * p_Nq + b * p_Nq + a;
            s_Ud[b][a] = Ud[id + 0 * offset];
            s_Vd[b][a] = Ud[id + 1 * offset];
            s_Wd[b][a] = Ud[id + 2 * offset];
        }
      }

      // interpolate in 'r'
      #pragma unroll
      for (int b = 0; b < p_Nq; ++b) {
        #pragma unroll
        for (int i = 0; i < p_cubNq; ++i) {
            dfloat Ud1 = 0, Vd1 = 0, Wd1 = 0;
            #pragma unroll
            for (int a = 0; a < p_Nq; ++a) {
              dfloat Iia = s_cubInterpT[a][i];
              Ud1 += Iia * s_Ud[b][a];
              Vd1 += Iia * s_Vd[b][a];
              Wd1 += Iia * s_Wd[b][a];
            }
            s_Ud1[b][i] = Ud1;
            s_Vd1[b][i] = Vd1;
            s_Wd1[b][i] = Wd1;
        }
      }

      // interpolate in 's'
      _occa_exclusive_index = 0;
      #pragma unroll
      for (int j = 0; j < p_cubNq; ++j) {
        #pragma unroll
        for (int i = 0; i < p_cubNq; ++i) {
          dfloat Ud2 = 0, Vd2 = 0, Wd2 = 0;

          // interpolate in b
          #pragma unroll
          for (int b = 0; b < p_Nq; ++b) {
            dfloat Ijb = s_cubInterpT[b][j];
            Ud2 += Ijb * s_Ud1[b][i];
            Vd2 += Ijb * s_Vd1[b][i];
            Wd2 += Ijb * s_Wd1[b][i];
          }

          // interpolate in c progressively
          #pragma unroll
          for (int k = 0; k < p_cubNq; ++k) {
            dfloat Ikc = s_cubInterpT[c][k];
            r_Ud[_occa_exclusive_index][k] += Ikc * Ud2;
            r_Vd[_occa_exclusive_index][k] += Ikc * Vd2;
            r_Wd[_occa_exclusive_index][k] += Ikc * Wd2;
          }
          ++_occa_exclusive_index;
        }
      }
    }
    #pragma unroll p_cubNq
    for (int k = 0; k < p_cubNq; ++k) {
      ;
      _occa_exclusive_index = 0;
      #pragma unroll
      for (int j = 0; j < p_cubNq; ++j) {
        #pragma unroll
        for (int i = 0; i < p_cubNq; ++i) {
          s_Ud[j][i] = r_Ud[_occa_exclusive_index][k];
          s_Vd[j][i] = r_Vd[_occa_exclusive_index][k];
          s_Wd[j][i] = r_Wd[_occa_exclusive_index][k];
          ++_occa_exclusive_index;
        }
      }
      ;
      _occa_exclusive_index = 0;
      #pragma unroll
      for (int j = 0; j < p_cubNq; ++j) {
        #pragma unroll
        for (int i = 0; i < p_cubNq; ++i) {
          dfloat Udr = 0, Uds = 0, Udt = 0;
          dfloat Vdr = 0, Vds = 0, Vdt = 0;
          dfloat Wdr = 0, Wds = 0, Wdt = 0;
          #pragma unroll
          for (int n = 0; n < p_cubNq; ++n) {
            dfloat Din = s_cubD[i][n];
            Udr += Din * s_Ud[j][n];
            Vdr += Din * s_Vd[j][n];
            Wdr += Din * s_Wd[j][n];
          }
          #pragma unroll
          for (int n = 0; n < p_cubNq; ++n) {
            dfloat Djn = s_cubD[j][n];
            Uds += Djn * s_Ud[n][i];
            Vds += Djn * s_Vd[n][i];
            Wds += Djn * s_Wd[n][i];
          }
          #pragma unroll
          for (int n = 0; n < p_cubNq; ++n) {
            dfloat Dkn = s_cubD[k][n];
            Udt += Dkn * r_Ud[_occa_exclusive_index][n];
            Vdt += Dkn * r_Vd[_occa_exclusive_index][n];
            Wdt += Dkn * r_Wd[_occa_exclusive_index][n];
          }
          dfloat Uhat = 0.0;
          dfloat Vhat = 0.0;
          dfloat What = 0.0;
          const int id = element * p_cubNp + k * p_cubNq * p_cubNq + j * p_cubNq + i;
          #pragma unroll
          for (int s = 0; s < p_nEXT; ++s) {
            const int s_offset = s * p_NVfields * cubatureOffset;
            const dfloat coeff = r_c[s];
            Uhat += coeff * conv[id + 0 * cubatureOffset + s_offset];
            Vhat += coeff * conv[id + 1 * cubatureOffset + s_offset];
            What += coeff * conv[id + 2 * cubatureOffset + s_offset];
          }

          // U*dUdx + V*dUdy + W*dUdz = (U*(drdx*dUdr+dsdx*dUds+dtdx*dUdt) + V*(drdy*dUdr ..))

          // I_f^t*(J_f*C_f^t)*G_f*\hat{D}_f*I_f*u
          r_U2[_occa_exclusive_index][k] = Uhat * Udr + Vhat * Uds + What * Udt;
          r_V2[_occa_exclusive_index][k] = Uhat * Vdr + Vhat * Vds + What * Vdt;
          r_W2[_occa_exclusive_index][k] = Uhat * Wdr + Vhat * Wds + What * Wdt;
          ++_occa_exclusive_index;
        }
      }
    }

    // now project back in t
    #pragma unroll
    for (int c = 0; c < p_Nq; ++c) {
      ;
      _occa_exclusive_index = 0;
      #pragma unroll
      for (int j = 0; j < p_cubNq; ++j) {
        #pragma unroll
        for (int i = 0; i < p_cubNq; ++i) {
          dfloat rhsU = 0, rhsV = 0, rhsW = 0;
          #pragma unroll
          for (int k = 0; k < p_cubNq; ++k) {
            dfloat Ikc = s_cubInterpT[c][k];
            rhsU += Ikc * r_U2[_occa_exclusive_index][k];
            rhsV += Ikc * r_V2[_occa_exclusive_index][k];
            rhsW += Ikc * r_W2[_occa_exclusive_index][k];
          }
          s_U[j][i] = rhsU;
          s_V[j][i] = rhsV;
          s_W[j][i] = rhsW;
          ++_occa_exclusive_index;
        }
      }
      #pragma unroll
      for (int b = 0; b < p_Nq; ++b) {
        #pragma unroll
        for (int i = 0; i < p_cubNq; ++i) {
            dfloat rhsU = 0, rhsV = 0, rhsW = 0;
            #pragma unroll
            for (int j = 0; j < p_cubNq; ++j) {
              dfloat Ijb = s_cubInterpT[b][j];
              rhsU += Ijb * s_U[j][i];
              rhsV += Ijb * s_V[j][i];
              rhsW += Ijb * s_W[j][i];
            }
            s_Ud[b][i] = rhsU;
            s_Vd[b][i] = rhsV;
            s_Wd[b][i] = rhsW;
        }
      }
      #pragma unroll
      for (int b = 0; b < p_Nq; ++b) {
        #pragma unroll
        for (int a = 0; a < p_Nq; ++a) {
            dfloat rhsU = 0, rhsV = 0, rhsW = 0;
            #pragma unroll
            for (int i = 0; i < p_cubNq; ++i) {
              dfloat Iia = s_cubInterpT[a][i];
              rhsU += Iia * s_Ud[b][i];
              rhsV += Iia * s_Vd[b][i];
              rhsW += Iia * s_Wd[b][i];
            }
            const int id = element * p_Np + c * p_Nq * p_Nq + b * p_Nq + a;
            dfloat invLMM = p_MovingMesh ? 0.0 : invLumpedMassMatrix[id];
            dfloat bdivw = 0.0;
            if (p_MovingMesh) {
              #pragma unroll
              for (int s = 0; s < p_nEXT; s++) {
                const dfloat coeff = r_c[s];
                invLMM += coeff * invLumpedMassMatrix[id + s * offset];
                bdivw += coeff * BdivW[id + s * offset];
              }
            }
            NU[id + 0 * offset + NUoffset] = (rhsU - bdivw * Ud[id + 0 * offset]) * invLMM;
            NU[id + 1 * offset + NUoffset] = (rhsV - bdivw * Ud[id + 1 * offset]) * invLMM;
            NU[id + 2 * offset + NUoffset] = (rhsW - bdivw * Ud[id + 2 * offset]) * invLMM;
        }
      }
    }
  }
}