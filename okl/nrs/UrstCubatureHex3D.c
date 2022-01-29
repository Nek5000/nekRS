extern "C" void FUNC(UrstCubatureHex3D)(const dlong &Nelements,
                                        const dfloat *__restrict__ cubvgeo,
                                        const dfloat *__restrict__ cubInterpT,
                                        const dlong &offset,
                                        const dlong &cubatureOffset,
                                        const dfloat *__restrict__ U,
                                        const dfloat *__restrict__ W,
                                        dfloat *__restrict__ result)
{
  dfloat s_cubInterpT[p_Nq][p_cubNq];
  dfloat s_U[p_Nq][p_Nq];
  dfloat s_V[p_Nq][p_Nq];
  dfloat s_W[p_Nq][p_Nq];
  dfloat s_U1[p_Nq][p_cubNq];
  dfloat s_V1[p_Nq][p_cubNq];
  dfloat s_W1[p_Nq][p_cubNq];
  dfloat r_U[p_cubNq][p_cubNq][p_cubNq], r_V[p_cubNq][p_cubNq][p_cubNq], r_W[p_cubNq][p_cubNq][p_cubNq];

  for (dlong j = 0; j < p_cubNq; ++j) {
    for (dlong i = 0; i < p_cubNq; ++i) {
      const dlong id = i + j * p_cubNq;
      if (id < p_Nq * p_cubNq) {
        s_cubInterpT[j][i] = cubInterpT[id];
      }
    }
  }

  for (dlong element = 0; element < Nelements; ++element) {
    for (dlong c = 0; c < p_Nq; ++c) {
      for (dlong b = 0; b < p_cubNq; ++b) {
        for (dlong a = 0; a < p_cubNq; ++a) {
          if (c == 0) {
            for (dlong k = 0; k < p_cubNq; ++k) {
              r_U[b][a][k] = 0.0;
              r_V[b][a][k] = 0.0;
              r_W[b][a][k] = 0.0;
            }
          }
          if (a < p_Nq && b < p_Nq) {
            const dlong id = element * p_Np + c * p_Nq * p_Nq + b * p_Nq + a;
            dfloat Ue = U[id + 0 * offset];
            dfloat Ve = U[id + 1 * offset];
            dfloat We = U[id + 2 * offset];
            if (p_relative) {
              Ue -= W[id + 0 * offset];
              Ve -= W[id + 1 * offset];
              We -= W[id + 2 * offset];
            }
            s_U[b][a] = Ue;
            s_V[b][a] = Ve;
            s_W[b][a] = We;
          }
        }
      }

      // dlongerpolate in 'r'
      for (dlong b = 0; b < p_cubNq; ++b) {
        for (dlong i = 0; i < p_cubNq; ++i) {
          if (b < p_Nq) {
            dfloat U1 = 0, V1 = 0, W1 = 0;
            for (dlong a = 0; a < p_Nq; ++a) {
              dfloat Iia = s_cubInterpT[a][i];
              U1 += Iia * s_U[b][a];
              V1 += Iia * s_V[b][a];
              W1 += Iia * s_W[b][a];
            }
            s_U1[b][i] = U1;
            s_V1[b][i] = V1;
            s_W1[b][i] = W1;
          }
        }
      }

      // dlongerpolate in 's'
      for (dlong j = 0; j < p_cubNq; ++j) {
        for (dlong i = 0; i < p_cubNq; ++i) {
          dfloat U2 = 0, V2 = 0, W2 = 0;

          // dlongerpolate in b
          for (dlong b = 0; b < p_Nq; ++b) {
            dfloat Ijb = s_cubInterpT[b][j];
            U2 += Ijb * s_U1[b][i];
            V2 += Ijb * s_V1[b][i];
            W2 += Ijb * s_W1[b][i];
          }

          // dlongerpolate in c progressively
          for (dlong k = 0; k < p_cubNq; ++k) {
            dfloat Ikc = s_cubInterpT[c][k];
            r_U[j][i][k] += Ikc * U2;
            r_V[j][i][k] += Ikc * V2;
            r_W[j][i][k] += Ikc * W2;
          }
        }
      }
    }
    for (dlong k = 0; k < p_cubNq; ++k) {
      for (dlong j = 0; j < p_cubNq; ++j) {
        for (dlong i = 0; i < p_cubNq; ++i) {
          const dlong gid = element * p_cubNp * p_Nvgeo + k * p_cubNq * p_cubNq + j * p_cubNq + i;

          const dfloat drdx = cubvgeo[gid + p_RXID * p_cubNp];
          const dfloat drdy = cubvgeo[gid + p_RYID * p_cubNp];
          const dfloat drdz = cubvgeo[gid + p_RZID * p_cubNp];
          const dfloat dsdx = cubvgeo[gid + p_SXID * p_cubNp];
          const dfloat dsdy = cubvgeo[gid + p_SYID * p_cubNp];
          const dfloat dsdz = cubvgeo[gid + p_SZID * p_cubNp];
          const dfloat dtdx = cubvgeo[gid + p_TXID * p_cubNp];
          const dfloat dtdy = cubvgeo[gid + p_TYID * p_cubNp];
          const dfloat dtdz = cubvgeo[gid + p_TZID * p_cubNp];
          const dfloat JW = cubvgeo[gid + p_JWID * p_cubNp];

          const dfloat Un = r_U[j][i][k];
          const dfloat Vn = r_V[j][i][k];
          const dfloat Wn = r_W[j][i][k];
          const dlong id = element * p_cubNp + k * p_cubNq * p_cubNq + j * p_cubNq + i;
          result[id + 0 * cubatureOffset] = JW * (Un * drdx + Vn * drdy + Wn * drdz);
          result[id + 1 * cubatureOffset] = JW * (Un * dsdx + Vn * dsdy + Wn * dsdz);
          result[id + 2 * cubatureOffset] = JW * (Un * dtdx + Vn * dtdy + Wn * dtdz);
        }
      }
    }
  }
}