extern "C" void FUNC(scaledAdd) (const dlong & N,
  const pfloat & alpha,
  const  pfloat *  __restrict__ x,
  const pfloat & beta,
  pfloat *  __restrict__ y){
  
#ifdef __NEKRS__OMP__
  #pragma omp parallel for
#endif
  for(dlong n = 0; n < N; ++n){
      y[n] = alpha*x[n] + beta*y[n];
  }
}