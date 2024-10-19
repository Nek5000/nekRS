// "pre" region from Algorithm 5, https://arxiv.org/pdf/2205.08909.pdf
extern "C" void FUNC(combinedPCGPreMatVec)(const dlong &N,
                                           const dlong &updateX,
                                           const dlong &preco,
                                           const dlong &fieldOffset,
                                           const dfloat &alphakm1,
                                           const dfloat &alphakm2,
                                           const dfloat &betakm1,
                                           const dfloat &betakm2,
                                           const dfloat &alphaDivBetakm2, // avoid division in x-update
                                           const dfloat *__restrict__ Minv,
                                           const dfloat *__restrict__ v,
                                           dfloat *__restrict__ p,
                                           dfloat *__restrict__ x,
                                           dfloat *__restrict__ r)
{
  if (updateX) {
    for (int fld = 0; fld < p_Nfields; fld++) {
      for (dlong n = 0; n < N; ++n) {
        const dlong id = n + fld * fieldOffset;
        const dfloat M = preco ? Minv[id] : 1.0;
        const dfloat pkm1 = p[id];
        const dfloat rkm1 = r[id];
        const dfloat vkm1 = v[id];
        x[id] += alphakm1 * pkm1 + alphaDivBetakm2 * (pkm1 - M * rkm1);
        const dfloat rk = rkm1 - alphakm1 * vkm1;
        r[id] = rk;
        p[id] = M * rk + betakm1 * pkm1;
      }
    }
  } else {
    for (int fld = 0; fld < p_Nfields; fld++) {
      for (dlong n = 0; n < N; ++n) {
        const dlong id = n + fld * fieldOffset;
        const dfloat M = preco ? Minv[id] : 1.0;
        const dfloat pkm1 = p[id];
        const dfloat rkm1 = r[id];
        const dfloat vkm1 = v[id];
        const dfloat rk = rkm1 - alphakm1 * vkm1;
        r[id] = rk;
        p[id] = M * rk + betakm1 * pkm1;
      }
    }
  }
}
