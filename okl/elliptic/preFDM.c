

extern "C" void FUNC(preFDM)(const dlong& Nelements,
                    const pfloat* __restrict__ u,
                    pfloat* __restrict__ work1)
{
  #define getIdx(k,j,i,e) ((k)*p_Nq_e*p_Nq_e+(j)*p_Nq_e+(i)+(e)*p_Nq_e*p_Nq_e*p_Nq_e)
  #define getIdx2(k,j,i,e) ((k-1)*p_Nq*p_Nq+(j-1)*p_Nq+(i-1)+(e)*p_Nq*p_Nq*p_Nq)
  #define sWork1(k,j,i,e) (work1[(getIdx(k,j,i,e))])
  #define uArr(k,j,i,e) (u[(getIdx2(k,j,i,e))])

#ifdef __NEKRS__OMP__
  #pragma omp parallel for
#endif
  for (dlong elem = 0; elem < Nelements; elem++) {
    #pragma unroll 
    for(int k = 0; k < p_Nq_e; ++k){
      #pragma unroll 
      for(int j = 0; j < p_Nq_e; ++j){
        #pragma unroll 
        for(int i = 0; i < p_Nq_e; ++i){
          const bool iBound = i>=1 && i <(p_Nq_e-1);
          const bool jBound = j>=1 && j <(p_Nq_e-1);
          const bool kBound = k>=1 && k <(p_Nq_e-1);
          if(iBound && jBound && kBound){
            const dlong elem_offset = elem * p_Nq * p_Nq * p_Nq;
            const dlong idx = i + j * p_Nq + k * p_Nq * p_Nq + elem_offset;
            sWork1(k,j,i,elem) = uArr(k,j,i,elem);
          } else {
            sWork1(k,j,i,elem) = 0.0;
          }
        }
      }
    }


    #pragma unroll 
    for(int j = 1; j < p_Nq_e-1; ++j){
      #pragma unroll 
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 2;
          sWork1(l1,j,k,elem) = uArr(l2,j,k,elem);
      }
    }
    #pragma unroll 
    for(int j = 1; j < p_Nq_e-1; ++j){
      #pragma unroll 
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 2;
          sWork1(p_Nq_e - l1 - 1,j,k,elem) = uArr(p_Nq_e - l2 - 1,j,k,elem);
      }
    }


    #pragma unroll 
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll 
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 2;
          sWork1(i,l1,k,elem) = uArr(i,l2,k,elem);
      }
    }
    #pragma unroll 
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll 
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 2;
          sWork1(i,p_Nq_e - l1 - 1,k,elem) = uArr(i,p_Nq_e - l2 - 1,k,elem);
      }
    }
    #pragma unroll 
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll 
      for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 0;
          const int l2 = 2;
          sWork1(i,j,l1,elem) = uArr(i,j,l2,elem);
      }
    }
    #pragma unroll 
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll 
      for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 0;
          const int l2 = 2;
          sWork1(i,j,p_Nq_e - l1 - 1,elem) = uArr(i,j,p_Nq_e - l2 - 1,elem);
      }
    }
  }
  #undef getIdx
  #undef getIdx2
  #undef sWork1
  #undef uArr
}
