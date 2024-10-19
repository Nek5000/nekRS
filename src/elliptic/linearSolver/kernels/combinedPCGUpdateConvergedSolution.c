// converged solution update from Algorithm 5, https://arxiv.org/pdf/2205.08909.pdf
extern "C" void
FUNC(combinedPCGUpdateConvergedSolution)(const dlong &N,
                                         const dlong &singleVectorUpdate,
                                         const dlong &preco,
                                         const dlong &fieldOffset,
                                         const dfloat &alphak,
                                         const dfloat &alphakm1,
                                         const dfloat &betakm1,
                                         const dfloat &alphaDivBetakm1, // avoid division in x-update
                                         const dfloat *__restrict__ Minv,
                                         const dfloat *__restrict__ p,
                                         const dfloat *__restrict__ r,
                                         dfloat *__restrict__ x)
{
  if (singleVectorUpdate) {
    for (int fld = 0; fld < p_Nfields; fld++) {
      for (dlong n = 0; n < N; ++n) {
        const dlong id = n + fld * fieldOffset;
        const dfloat pk = p[id];
        x[id] += alphak * pk;
      }
    }
  } else {
    for (int fld = 0; fld < p_Nfields; fld++) {
      for (dlong n = 0; n < N; ++n) {
        const dlong id = n + fld * fieldOffset;
        const dfloat pk = p[id];
        const dfloat rk = r[id];
        const dfloat M = preco ? Minv[id] : 1.0;
        x[id] += alphak * pk + alphaDivBetakm1 * (pk - M * rk);
      }
    }
  }
}
