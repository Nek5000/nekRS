// This kernel is needed as it used for mixed-precision Jacobi preconditioning 
extern "C" void FUNC(axmyzManyPfloat)(const dlong & N,
                   const dlong & Nfields,
                   const dlong & offset,
                   const pfloat & alpha,
                   const dfloat * __restrict__ x,
                   const pfloat * __restrict__ y,
                   dfloat * __restrict__ z){

#ifdef __NEKRS__OMP__ 
  #pragma omp parallel for
#endif
  for(dlong n=0;n<N;++n){
    for(int fld = 0; fld < Nfields; ++fld){
      const int id = n + fld *offset;
      z[id] = dfloat(alpha*x[id]*y[id]);
    }
  }
}
