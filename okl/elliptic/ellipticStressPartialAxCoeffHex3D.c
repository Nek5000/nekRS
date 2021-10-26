extern "C" void FUNC(ellipticStressPartialAxCoeffHex3D)(const dlong &Nelements,
                              const dlong &offset,
                              const dlong &loffset,
                              const dlong* __restrict__ elementList,
                              const dfloat* __restrict__ vgeo,
                              const dfloat* __restrict__ D,
                              const dfloat* __restrict__ S,
                              const dfloat* __restrict__ lambda,
                              const dfloat* __restrict__ q,
                              dfloat* __restrict__ Aq)
{
  dfloat s_D[p_Nq][p_Nq];

  dfloat s_U[p_Nq][p_Nq][p_Nq];
  dfloat s_V[p_Nq][p_Nq][p_Nq];
  dfloat s_W[p_Nq][p_Nq][p_Nq];

  dfloat s_SUr[p_Nq][p_Nq][p_Nq];
  dfloat s_SUs[p_Nq][p_Nq][p_Nq];
  dfloat s_SUt[p_Nq][p_Nq][p_Nq];

  dfloat s_SVr[p_Nq][p_Nq][p_Nq];
  dfloat s_SVs[p_Nq][p_Nq][p_Nq];
  dfloat s_SVt[p_Nq][p_Nq][p_Nq];

  dfloat s_SWr[p_Nq][p_Nq][p_Nq];
  dfloat s_SWs[p_Nq][p_Nq][p_Nq];
  dfloat s_SWt[p_Nq][p_Nq][p_Nq];

  for(int j = 0; j < p_Nq; ++j)
    for(int i = 0; i < p_Nq; ++i)
      s_D[j][i] = D[j * p_Nq + i];

#ifdef __NEKRS__OMP__
  #pragma omp parallel for private(s_U, s_V, s_W, s_SUr, s_SUs, s_SUt, s_SVr, s_SVs, s_SVt, s_SWr, s_SWs, s_SWt)
#endif
  for(dlong elem = 0; elem < Nelements; ++elem) {
    dlong e = elementList[elem];
    for(int k = 0; k < p_Nq; ++k)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong id = e * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          s_U[k][j][i] = q[id + 0 * offset];
          s_V[k][j][i] = q[id + 1 * offset];
          s_W[k][j][i] = q[id + 2 * offset];
        }

    // loop over slabs
    for(int k = 0; k < p_Nq; ++k)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gid = i + j * p_Nq + k * p_Nq * p_Nq + e * p_Np * p_Nvgeo;
          const dfloat rx = vgeo[gid + p_RXID * p_Np];
          const dfloat ry = vgeo[gid + p_RYID * p_Np];
          const dfloat rz = vgeo[gid + p_RZID * p_Np];

          const dfloat sx = vgeo[gid + p_SXID * p_Np];
          const dfloat sy = vgeo[gid + p_SYID * p_Np];
          const dfloat sz = vgeo[gid + p_SZID * p_Np];

          const dfloat tx = vgeo[gid + p_TXID * p_Np];
          const dfloat ty = vgeo[gid + p_TYID * p_Np];
          const dfloat tz = vgeo[gid + p_TZID * p_Np];

          const dfloat JW = vgeo[gid + p_JWID * p_Np];

          // compute 1D derivatives
          dfloat ur = 0.f, us = 0.f, ut = 0.f;
          dfloat vr = 0.f, vs = 0.f, vt = 0.f;
          dfloat wr = 0.f, ws = 0.f, wt = 0.f;
          for(int m = 0; m < p_Nq; ++m) {
            const dfloat Dim = s_D[i][m]; // Dr
            const dfloat Djm = s_D[j][m]; // Ds
            const dfloat Dkm = s_D[k][m]; // Dt

            ur += Dim * s_U[k][j][m];
            us += Djm * s_U[k][m][i];
            ut += Dkm * s_U[m][j][i];
            //
            vr += Dim * s_V[k][j][m];
            vs += Djm * s_V[k][m][i];
            vt += Dkm * s_V[m][j][i];
            //
            wr += Dim * s_W[k][j][m];
            ws += Djm * s_W[k][m][i];
            wt += Dkm * s_W[m][j][i];
          }

          const dlong id = e * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat u_lam0 = lambda[id + 0 * offset + 0 * loffset];
          const dfloat v_lam0 = lambda[id + 0 * offset + 1 * loffset];
          const dfloat w_lam0 = lambda[id + 0 * offset + 2 * loffset];

          const dfloat dudx = rx * ur + sx * us + tx * ut;
          const dfloat dudy = ry * ur + sy * us + ty * ut;
          const dfloat dudz = rz * ur + sz * us + tz * ut;

          const dfloat dvdx = rx * vr + sx * vs + tx * vt;
          const dfloat dvdy = ry * vr + sy * vs + ty * vt;
          const dfloat dvdz = rz * vr + sz * vs + tz * vt;

          const dfloat dwdx = rx * wr + sx * ws + tx * wt;
          const dfloat dwdy = ry * wr + sy * ws + ty * wt;
          const dfloat dwdz = rz * wr + sz * ws + tz * wt;

          const dfloat s11 = u_lam0 * JW * (dudx + dudx);
          const dfloat s12 = u_lam0 * JW * (dudy + dvdx);
          const dfloat s13 = u_lam0 * JW * (dudz + dwdx);

          const dfloat s21 = v_lam0 * JW * (dvdx + dudy);
          const dfloat s22 = v_lam0 * JW * (dvdy + dvdy);
          const dfloat s23 = v_lam0 * JW * (dvdz + dwdy);

          const dfloat s31 = w_lam0 * JW * (dwdx + dudz);
          const dfloat s32 = w_lam0 * JW * (dwdy + dvdz);
          const dfloat s33 = w_lam0 * JW * (dwdz + dwdz);

          s_SUr[k][j][i] =  rx * s11 + ry * s12 + rz * s13;
          s_SUs[k][j][i] =  sx * s11 + sy * s12 + sz * s13;
          s_SUt[k][j][i] =  tx * s11 + ty * s12 + tz * s13;
          //
          s_SVr[k][j][i] =  rx * s21 + ry * s22 + rz * s23;
          s_SVs[k][j][i] =  sx * s21 + sy * s22 + sz * s23;
          s_SVt[k][j][i] =  tx * s21 + ty * s22 + tz * s23;
          //
          s_SWr[k][j][i] =  rx * s31 + ry * s32 + rz * s33;
          s_SWs[k][j][i] =  sx * s31 + sy * s32 + sz * s33;
          s_SWt[k][j][i] =  tx * s31 + ty * s32 + tz * s33;
        }

// loop over slabs
    for(int k = 0; k < p_Nq; ++k)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          dfloat r_Au = 0.f, r_Av = 0.f, r_Aw = 0.f;
          for(int m = 0; m < p_Nq; m++) {
            const dfloat Dim = s_D[m][i]; // Dr'
            const dfloat Djm = s_D[m][j]; // Ds'
            const dfloat Dkm = s_D[m][k]; // Dt'

            r_Au += Dim * s_SUr[k][j][m];
            r_Au += Djm * s_SUs[k][m][i];
            r_Au += Dkm * s_SUt[m][j][i];

            r_Av += Dim * s_SVr[k][j][m];
            r_Av += Djm * s_SVs[k][m][i];
            r_Av += Dkm * s_SVt[m][j][i];

            r_Aw += Dim * s_SWr[k][j][m];
            r_Aw += Djm * s_SWs[k][m][i];
            r_Aw += Dkm * s_SWt[m][j][i];
          }
          const dlong id      = e * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat u_lam1 = lambda[id + 1 * offset + 0 * loffset];
          const dfloat v_lam1 = lambda[id + 1 * offset + 1 * loffset];
          const dfloat w_lam1 = lambda[id + 1 * offset + 2 * loffset];

          const dlong gid = i + j * p_Nq + k * p_Nq * p_Nq + e * p_Np * p_Nvgeo;
          const dfloat JW = vgeo[gid + p_JWID * p_Np];
          // store in register
          Aq[id + 0 * offset] =  r_Au + u_lam1 * JW * s_U[k][j][i];
          Aq[id + 1 * offset] =  r_Av + v_lam1 * JW * s_V[k][j][i];
          Aq[id + 2 * offset] =  r_Aw + w_lam1 * JW * s_W[k][j][i];
        }
  }
}
