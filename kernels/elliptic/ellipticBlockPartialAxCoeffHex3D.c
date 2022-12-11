extern "C" void FUNC(ellipticBlockPartialAxCoeffHex3D)(const dlong & Nelements,
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
  dfloat s_q[3][p_Nq][p_Nq][p_Nq];
  dfloat s_Gqr[3][p_Nq][p_Nq][p_Nq];
  dfloat s_Gqs[3][p_Nq][p_Nq][p_Nq];
  dfloat s_Gqt[3][p_Nq][p_Nq][p_Nq];

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
          s_q[0][k][j][i] = q[base + 0 * offset];
          s_q[1][k][j][i] = q[base + 1 * offset];
          s_q[2][k][j][i] = q[base + 2 * offset];
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

          const dlong id      = element * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;

          const dfloat r_lam00 = lambda[id + 0 * offset + 0 * loffset];
          const dfloat r_lam10 = lambda[id + 0 * offset + 1 * loffset];
          const dfloat r_lam20 = lambda[id + 0 * offset + 2 * loffset];

          dfloat qr0 = 0.f, qr1 = 0.f, qr2 = 0.f;
          dfloat qs0 = 0.f, qs1 = 0.f, qs2 = 0.f;
          dfloat qt0 = 0.f, qt1 = 0.f, qt2 = 0.f;

          for(int m = 0; m < p_Nq; m++) {
            qr0 += s_S[m][i] * s_q[0][k][j][m];
            qs0 += s_S[m][j] * s_q[0][k][m][i];
            qt0 += s_S[m][k] * s_q[0][m][j][i];
            //
            qr1 += s_S[m][i] * s_q[1][k][j][m];
            qs1 += s_S[m][j] * s_q[1][k][m][i];
            qt1 += s_S[m][k] * s_q[1][m][j][i];

            qr2 += s_S[m][i] * s_q[2][k][j][m];
            qs2 += s_S[m][j] * s_q[2][k][m][i];
            qt2 += s_S[m][k] * s_q[2][m][j][i];
          }

          dfloat Gqr0 = r_G00 * qr0 + r_G01 * qs0 + r_G02 * qt0;
          dfloat Gqs0 = r_G01 * qr0 + r_G11 * qs0 + r_G12 * qt0;
          dfloat Gqt0 = r_G02 * qr0 + r_G12 * qs0 + r_G22 * qt0;

          dfloat Gqr1 = r_G00 * qr1 + r_G01 * qs1 + r_G02 * qt1;
          dfloat Gqs1 = r_G01 * qr1 + r_G11 * qs1 + r_G12 * qt1;
          dfloat Gqt1 = r_G02 * qr1 + r_G12 * qs1 + r_G22 * qt1;

          dfloat Gqr2 = r_G00 * qr2 + r_G01 * qs2 + r_G02 * qt2;
          dfloat Gqs2 = r_G01 * qr2 + r_G11 * qs2 + r_G12 * qt2;
          dfloat Gqt2 = r_G02 * qr2 + r_G12 * qs2 + r_G22 * qt2;

          s_Gqr[0][k][j][i] = r_lam00 * Gqr0;
          s_Gqs[0][k][j][i] = r_lam00 * Gqs0;
          s_Gqt[0][k][j][i] = r_lam00 * Gqt0;

          s_Gqr[1][k][j][i] = r_lam10 * Gqr1;
          s_Gqs[1][k][j][i] = r_lam10 * Gqs1;
          s_Gqt[1][k][j][i] = r_lam10 * Gqt1;

          s_Gqr[2][k][j][i] = r_lam20 * Gqr2;
          s_Gqs[2][k][j][i] = r_lam20 * Gqs2;
          s_Gqt[2][k][j][i] = r_lam20 * Gqt2;
        }

    for(int k = 0; k < p_Nq; k++)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase  = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_GwJ = ggeo[gbase + p_GWJID * p_Np];

          const dlong id = element * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;

          const dfloat r_lam01 = lambda[id + 1 * offset + 0 * loffset];
          const dfloat r_lam11 = lambda[id + 1 * offset + 1 * loffset];
          const dfloat r_lam21 = lambda[id + 1 * offset + 2 * loffset];

          dfloat r_Aq0 = r_GwJ * r_lam01 * s_q[0][k][j][i];
          dfloat r_Aq1 = r_GwJ * r_lam11 * s_q[1][k][j][i];
          dfloat r_Aq2 = r_GwJ * r_lam21 * s_q[2][k][j][i];

          dfloat r_Aqr0 = 0.f, r_Aqs0 = 0.f, r_Aqt0 = 0.f;
          dfloat r_Aqr1 = 0.f, r_Aqs1 = 0.f, r_Aqt1 = 0.f;
          dfloat r_Aqr2 = 0.f, r_Aqs2 = 0.f, r_Aqt2 = 0.f;

          for(int m = 0; m < p_Nq; m++) {
            r_Aqr0 += s_D[m][i] * s_Gqr[0][k][j][m];
            r_Aqr1 += s_D[m][i] * s_Gqr[1][k][j][m];
            r_Aqr2 += s_D[m][i] * s_Gqr[2][k][j][m];
          }

          for(int m = 0; m < p_Nq; m++) {
            r_Aqs0 += s_D[m][j] * s_Gqs[0][k][m][i];
            r_Aqs1 += s_D[m][j] * s_Gqs[1][k][m][i];
            r_Aqs2 += s_D[m][j] * s_Gqs[2][k][m][i];
          }

          for(int m = 0; m < p_Nq; m++) {
            r_Aqt0 += s_D[m][k] * s_Gqt[0][m][j][i];
            r_Aqt1 += s_D[m][k] * s_Gqt[1][m][j][i];
            r_Aqt2 += s_D[m][k] * s_Gqt[2][m][j][i];
          }

          Aq[id + 0 * offset] = r_Aqr0 + r_Aqs0 + r_Aqt0 + r_Aq0;
          Aq[id + 1 * offset] = r_Aqr1 + r_Aqs1 + r_Aqt1 + r_Aq1;
          Aq[id + 2 * offset] = r_Aqr2 + r_Aqs2 + r_Aqt2 + r_Aq2;
        }
  }
}
