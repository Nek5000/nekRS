

extern "C" void FUNC(postFDM)(const dlong& Nelements,
                     pfloat* __restrict__ my_work1,
                     pfloat* __restrict__ my_work2,
                     pfloat* __restrict__ Su,
                     const pfloat* __restrict__ wts)
{
  pfloat work1[p_Nq_e][p_Nq_e][p_Nq_e];
  pfloat work2[p_Nq_e][p_Nq_e][p_Nq_e];
#ifdef __NEKRS__OMP__
  #pragma omp parallel for private(work1, work2)
#endif
  for (dlong elem = 0; elem < Nelements; ++elem) {
    #pragma unroll
    for(int k = 0; k < p_Nq_e; ++k){
      #pragma unroll
      for(int j = 0; j < p_Nq_e; ++j){
        #pragma unroll
        for(int i = 0; i < p_Nq_e; ++i) {
          const dlong elem_offset = elem * p_Nq_e * p_Nq_e * p_Nq_e;
          const dlong idx = i + j * p_Nq_e + k * p_Nq_e * p_Nq_e + elem_offset;
          work1[k][j][i] = my_work2[idx];
          work2[k][j][i] = my_work1[idx];
        }
      }
    }
    #pragma unroll
    for(int j = 1; j < p_Nq_e-1; ++j){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 0;
          work1[l1][j][k] = work1[l1][j][k] - work2[l2][j][k];
      }
    }
    #pragma unroll
    for(int j = 1; j < p_Nq_e-1; ++j){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 0;
          work1[p_Nq_e - l1 - 1][j][k] = work1[p_Nq_e - l1 - 1][j][k] -
                                         work2[p_Nq_e - l2 - 1][j][k];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 0;
          work1[i][l1][k] = work1[i][l1][k] - work2[i][l2][k];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 0;
          work1[i][p_Nq_e - l1 - 1][k] = work1[i][p_Nq_e - l1 - 1][k] -
                                         work2[i][p_Nq_e - l2 - 1][k];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 0;
          const int l2 = 0;
          work1[i][j][l1] = work1[i][j][l1] - work2[i][j][l2];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 0;
          const int l2 = 0;
          work1[i][j][p_Nq_e - l1 - 1] = work1[i][j][p_Nq_e - l1 - 1] -
                                         work2[i][j][p_Nq_e - l2 - 1];
      }
    }
    #pragma unroll
    for(int j = 1; j < p_Nq_e-1; ++j){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 2;
          const int l2 = 0;
          work1[l1][j][k] = work1[l1][j][k] + work1[l2][j][k];
      }
    }
    #pragma unroll
    for(int j = 1; j < p_Nq_e-1; ++j){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 2;
          const int l2 = 0;
          work1[p_Nq_e - l1 - 1][j][k] = work1[p_Nq_e - l1 - 1][j][k] +
                                         work1[p_Nq_e - l2 - 1][j][k];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 2;
          const int l2 = 0;
          work1[i][l1][k] = work1[i][l1][k] + work1[i][l2][k];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 2;
          const int l2 = 0;
          work1[i][p_Nq_e - l1 - 1][k] = work1[i][p_Nq_e - l1 - 1][k] +
                                         work1[i][p_Nq_e - l2 - 1][k];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 2;
          const int l2 = 0;
          work1[i][j][l1] = work1[i][j][l1] + work1[i][j][l2];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 2;
          const int l2 = 0;
          work1[i][j][p_Nq_e - l1 - 1] = work1[i][j][p_Nq_e - l1 - 1] +
                                         work1[i][j][p_Nq_e - l2 - 1];
      }
    }
    #pragma unroll
    for(int k = 0; k < p_Nq; ++k){
      #pragma unroll
      for(int j = 0; j < p_Nq; ++j){
        #pragma unroll
        for(int i = 0; i < p_Nq; ++i){
            const dlong elem_offset = elem * p_Nq * p_Nq * p_Nq;
            const dlong idx = i + j * p_Nq + k * p_Nq * p_Nq + elem_offset;
            Su[idx] = work1[k + 1][j + 1][i + 1] * wts[idx];
        }
      }
    }
  }
}
