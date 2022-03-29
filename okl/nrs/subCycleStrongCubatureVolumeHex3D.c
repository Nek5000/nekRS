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

extern "C" void FUNC(subCycleStrongCubatureVolumeHex3D)(const int & Nelements,
                                                  const int * __restrict__ elementList,
                                                  const dfloat * __restrict__ cubD,
                                                  const dfloat * __restrict__ cubInterpT,
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
  dfloat s_U[p_cubNq][p_cubNq];
  dfloat s_V[p_cubNq][p_cubNq];
  dfloat s_W[p_cubNq][p_cubNq];
  dfloat s_Ud[p_cubNq][p_cubNq];
  dfloat s_Vd[p_cubNq][p_cubNq];
  dfloat s_Wd[p_cubNq][p_cubNq];
  dfloat s_Ud1[p_Nq][p_cubNq];
  dfloat s_Vd1[p_Nq][p_cubNq];
  dfloat s_Wd1[p_Nq][p_cubNq];
  dfloat r_U2[p_cubNq][p_cubNq][p_cubNq], r_V2[p_cubNq][p_cubNq][p_cubNq], r_W2[p_cubNq][p_cubNq][p_cubNq];
  dfloat r_Ud[p_cubNq][p_cubNq][p_cubNq], r_Vd[p_cubNq][p_cubNq][p_cubNq], r_Wd[p_cubNq][p_cubNq][p_cubNq];
  for (int j = 0; j < p_cubNq; ++j) {
    for (int i = 0; i < p_cubNq; ++i) {
      const int id = i + j * p_cubNq;
      if (id < p_Nq * p_cubNq) {
        s_cubInterpT[j][i] = cubInterpT[id];
      }
      s_cubD[j][i] = cubD[id];
    }
  }

#ifdef __NEKRS__OMP__
  #pragma omp parallel for private(s_U, s_V, s_W, s_Ud, s_Vd, s_Wd, s_Ud1, s_Vd1, s_Wd1, r_U2, r_V2, r_W2, r_Ud, r_Vd, r_Wd)
#endif
  for (int e = 0; e < Nelements; ++e) {
    const int element = elementList[e];
    for (int j = 0; j < p_cubNq; ++j) {
      for (int i = 0; i < p_cubNq; ++i) {
        for (int k = 0; k < p_cubNq; ++k) {
          r_Ud[j][i][k] = 0;
          r_Vd[j][i][k] = 0;
          r_Wd[j][i][k] = 0;
        }
      }
    }
    for (int c = 0; c < p_Nq; ++c) {
      for (int b = 0; b < p_Nq; ++b) {
        for (int a = 0; a < p_Nq; ++a) {
            const int id = element * p_Np + c * p_Nq * p_Nq + b * p_Nq + a;
            s_Ud[b][a] = Ud[id + 0 * offset];
            s_Vd[b][a] = Ud[id + 1 * offset];
            s_Wd[b][a] = Ud[id + 2 * offset];
        }
      }

      // interpolate in 'r'
      for (int b = 0; b < p_Nq; ++b) {
        for (int i = 0; i < p_cubNq; ++i) {
            dfloat Ud1 = 0, Vd1 = 0, Wd1 = 0;
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
      for (int j = 0; j < p_cubNq; ++j) {
        for (int i = 0; i < p_cubNq; ++i) {
          dfloat Ud2 = 0, Vd2 = 0, Wd2 = 0;

          // interpolate in b
          for (int b = 0; b < p_Nq; ++b) {
            dfloat Ijb = s_cubInterpT[b][j];
            Ud2 += Ijb * s_Ud1[b][i];
            Vd2 += Ijb * s_Vd1[b][i];
            Wd2 += Ijb * s_Wd1[b][i];
          }

          // interpolate in c progressively
          for (int k = 0; k < p_cubNq; ++k) {
            dfloat Ikc = s_cubInterpT[c][k];
            r_Ud[j][i][k] += Ikc * Ud2;
            r_Vd[j][i][k] += Ikc * Vd2;
            r_Wd[j][i][k] += Ikc * Wd2;
          }
        }
      }
    }

    // Uhat * dr
    for (int j = 0; j < p_cubNq; ++j) {
      for (int k = 0; k < p_cubNq; ++k) {
        for (int i = 0; i < p_cubNq; ++i) {
          dfloat Udr = 0;
          dfloat Vdr = 0;
          dfloat Wdr = 0;
          for (int n = 0; n < p_cubNq; ++n) {
            dfloat Din = s_cubD[i][n];
            Udr += Din * r_Ud[j][n][k];
            Vdr += Din * r_Vd[j][n][k];
            Wdr += Din * r_Wd[j][n][k];
          }
          dfloat Uhat = 0.0;
          const int id = element * p_cubNp + k * p_cubNq * p_cubNq + j * p_cubNq + i;
          for (int s = 0; s < p_nEXT; ++s) {
            const int s_offset = s * p_NVfields * cubatureOffset;
            const dfloat coeff = r_c[s];
            Uhat += coeff * conv[id + 0 * cubatureOffset + s_offset];
          }

          // U*dUdx + V*dUdy + W*dUdz = (U*(drdx*dUdr+dsdx*dUds+dtdx*dUdt) + V*(drdy*dUdr ..))

          // I_f^t*(J_f*C_f^t)*G_f*\hat{D}_f*I_f*u
          r_U2[j][i][k] = Uhat * Udr;
          r_V2[j][i][k] = Uhat * Vdr;
          r_W2[j][i][k] = Uhat * Wdr;
        }
      }
    }
    // Vhat * ds
    for (int j = 0; j < p_cubNq; ++j) {
      for (int i = 0; i < p_cubNq; ++i) {
        for (int k = 0; k < p_cubNq; ++k) {
          dfloat Uds = 0;
          dfloat Vds = 0;
          dfloat Wds = 0;
          for (int n = 0; n < p_cubNq; ++n) {
            dfloat Djn = s_cubD[j][n];
            Uds += Djn * r_Ud[n][i][k];
            Vds += Djn * r_Vd[n][i][k];
            Wds += Djn * r_Wd[n][i][k];
          }
          dfloat Vhat = 0.0;
          const int id = element * p_cubNp + k * p_cubNq * p_cubNq + j * p_cubNq + i;
          for (int s = 0; s < p_nEXT; ++s) {
            const int s_offset = s * p_NVfields * cubatureOffset;
            const dfloat coeff = r_c[s];
            Vhat += coeff * conv[id + 1 * cubatureOffset + s_offset];
          }

          // U*dUdx + V*dUdy + W*dUdz = (U*(drdx*dUdr+dsdx*dUds+dtdx*dUdt) + V*(drdy*dUdr ..))

          // I_f^t*(J_f*C_f^t)*G_f*\hat{D}_f*I_f*u
          r_U2[j][i][k] += Vhat * Uds;
          r_V2[j][i][k] += Vhat * Vds;
          r_W2[j][i][k] += Vhat * Wds;
        }
      }
    }
    // What * dt
    for (int j = 0; j < p_cubNq; ++j) {
      for (int k = 0; k < p_cubNq; ++k) {
        for (int i = 0; i < p_cubNq; ++i) {
          dfloat Udt = 0;
          dfloat Vdt = 0;
          dfloat Wdt = 0;
          for (int n = 0; n < p_cubNq; ++n) {
            dfloat Dkn = s_cubD[k][n];
            Udt += Dkn * r_Ud[j][i][n];
            Vdt += Dkn * r_Vd[j][i][n];
            Wdt += Dkn * r_Wd[j][i][n];
          }
          dfloat What = 0.0;
          const int id = element * p_cubNp + k * p_cubNq * p_cubNq + j * p_cubNq + i;
          for (int s = 0; s < p_nEXT; ++s) {
            const int s_offset = s * p_NVfields * cubatureOffset;
            const dfloat coeff = r_c[s];
            What += coeff * conv[id + 2 * cubatureOffset + s_offset];
          }

          // U*dUdx + V*dUdy + W*dUdz = (U*(drdx*dUdr+dsdx*dUds+dtdx*dUdt) + V*(drdy*dUdr ..))

          // I_f^t*(J_f*C_f^t)*G_f*\hat{D}_f*I_f*u
          r_U2[j][i][k] += What * Udt;
          r_V2[j][i][k] += What * Vdt;
          r_W2[j][i][k] += What * Wdt;
        }
      }
    }

    // now project back in t
    for (int c = 0; c < p_Nq; ++c) {
      for (int j = 0; j < p_cubNq; ++j) {
        for (int i = 0; i < p_cubNq; ++i) {
          dfloat rhsU = 0, rhsV = 0, rhsW = 0;
          for (int k = 0; k < p_cubNq; ++k) {
            dfloat Ikc = s_cubInterpT[c][k];
            rhsU += Ikc * r_U2[j][i][k];
            rhsV += Ikc * r_V2[j][i][k];
            rhsW += Ikc * r_W2[j][i][k];
          }
          s_U[j][i] = rhsU;
          s_V[j][i] = rhsV;
          s_W[j][i] = rhsW;
        }
      }
      for (int b = 0; b < p_Nq; ++b) {
        for (int i = 0; i < p_cubNq; ++i) {
            dfloat rhsU = 0, rhsV = 0, rhsW = 0;
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
      for (int b = 0; b < p_Nq; ++b) {
        for (int a = 0; a < p_Nq; ++a) {
            dfloat rhsU = 0, rhsV = 0, rhsW = 0;
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
