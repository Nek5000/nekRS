// "post" region from Algorithm 5, https://arxiv.org/pdf/2205.08909.pdf
extern "C" void FUNC(combinedPCGPostMatVec)(const dlong &N,
                                            const dlong &fieldOffset,
                                            const dlong &preco,
                                            const dfloat *__restrict__ resWeight,
                                            const dfloat *__restrict__ weights,
                                            const dfloat *__restrict__ Minv,
                                            const dfloat *__restrict__ v,
                                            const dfloat *__restrict__ p,
                                            const dfloat *__restrict__ r,
                                            dfloat *__restrict__ reduction)
{
  dfloat sums[p_nReduction];
  for (int i = 0; i < p_nReduction; ++i) {
    sums[i] = 0.0;
  }

  for (int id = 0; id < N; ++id) {
    const dfloat wt = weights[id];
    const dfloat resWt = resWeight[id];

    for (int fld = 0; fld < p_Nfields; ++fld) {
      const dlong n = id + fld * fieldOffset;
      const dfloat M = preco ? Minv[n] : 1.0;
      const dfloat pk = p[n];
      const dfloat rk = r[n];
      const dfloat vk = v[n];
      const dfloat Mvk = M * vk;
      const dfloat Mrk = M * rk;

      sums[p_gamma] += rk * rk * resWt;
      sums[p_a] += pk * vk * wt;
      sums[p_b] += rk * vk * wt;
      sums[p_c] += vk * vk * wt;
      sums[p_d] += rk * Mrk * wt;
      sums[p_e] += rk * Mvk * wt;
      sums[p_f] += vk * Mvk * wt;
    }
  }
  for (int i = 0; i < p_nReduction; ++i) {
    reduction[i] = sums[i];
  }
}
