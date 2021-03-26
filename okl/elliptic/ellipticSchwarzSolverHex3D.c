extern "C" void preFDM(const dlong& Nelements,
                    const pfloat* __restrict__ u,
                    pfloat* __restrict__ work1)
{
  #define getIdx(k,j,i,e) ((k)*p_Nq_e*p_Nq_e+(j)*p_Nq_e+(i)+(e)*p_Nq_e*p_Nq_e*p_Nq_e)
  #define getIdx2(k,j,i,e) ((k-1)*p_Nq*p_Nq+(j-1)*p_Nq+(i-1)+(e)*p_Nq*p_Nq*p_Nq)
  #define sWork1(k,j,i,e) (work1[(getIdx(k,j,i,e))])
  #define uArr(k,j,i,e) (u[(getIdx2(k,j,i,e))])
  #pragma omp parallel for
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

extern "C" void postFDM(const dlong& Nelements,
                     pfloat* __restrict__ my_work1,
                     pfloat* __restrict__ my_work2,
                     pfloat* __restrict__ Su,
                     const pfloat* __restrict__ wts)
{
  pfloat work1[p_Nq_e][p_Nq_e][p_Nq_e];
  pfloat work2[p_Nq_e][p_Nq_e][p_Nq_e];
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

extern "C" void fusedFDM(
  const dlong& Nelements,
  const dlong& localNelements,
  const dlong* __restrict__  elementList,
  pfloat* __restrict__ Su,
  const pfloat* __restrict__ S_x,
  const pfloat* __restrict__ S_y,
  const pfloat* __restrict__ S_z,
  const pfloat* __restrict__ inv_L,
#if p_restrict
  const dfloat* __restrict__ wts,
#endif
  pfloat* __restrict__ u
  )
{
#define getIdx(k,j,i,e) ((k)*p_Nq_e*p_Nq_e+(j)*p_Nq_e+(i)+(e)*p_Nq_e*p_Nq_e*p_Nq_e)
#define work1(k,j,i,e) (u[(getIdx(k,j,i,e))])
  pfloat S_x_e[p_Nq_e][p_Nq_e];
  pfloat S_y_e[p_Nq_e][p_Nq_e];
  pfloat S_z_e[p_Nq_e][p_Nq_e];
  pfloat S_x_eT[p_Nq_e][p_Nq_e];
  pfloat S_y_eT[p_Nq_e][p_Nq_e];
  pfloat S_z_eT[p_Nq_e][p_Nq_e];
  pfloat tmp[p_Nq_e][p_Nq_e][p_Nq_e];
  pfloat work2[p_Nq_e][p_Nq_e][p_Nq_e];

  for (dlong my_elem = 0; my_elem < Nelements; ++my_elem) {
    const dlong element = my_elem;
    const dlong elem = element;
    #pragma unroll
    for(int j = 1; j < p_Nq_e-1; ++j){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 2;
          work1(l1,j,k,elem) = work1(l1,j,k,elem) - work1(l2,j,k,elem);
      }
    }
    #pragma unroll
    for(int j = 1; j < p_Nq_e-1; ++j){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 2;
          work1(p_Nq_e - l1 - 1,j,k,elem) = work1(p_Nq_e - l1 - 1,j,k,elem) -
                                         work1(p_Nq_e - l2 - 1,j,k,elem);
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 2;
          work1(i,l1,k,elem) = work1(i,l1,k,elem) - work1(i,l2,k,elem);
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 2;
          work1(i,p_Nq_e - l1 - 1,k,elem) = work1(i,p_Nq_e - l1 - 1,k,elem) -
                                         work1(i,p_Nq_e - l2 - 1,k,elem);
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 0;
          const int l2 = 2;
          work1(i,j,l1,elem) = work1(i,j,l1,elem) - work1(i,j,l2,elem);
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 0;
          const int l2 = 2;
          work1(i,j,p_Nq_e - l1 - 1,elem) = work1(i,j,p_Nq_e - l1 - 1,elem) -
                                         work1(i,j,p_Nq_e - l2 - 1,elem);
      }
    }
    #pragma unroll
    for (int i = 0; i < p_Nq_e; i++){
      #pragma unroll
      for (int j = 0; j < p_Nq_e; j++) {
        const int ij = j + i * p_Nq_e;
        S_x_e[i][j] = S_x[ij + element * p_Nq_e * p_Nq_e];
        S_y_e[i][j] = S_y[ij + element * p_Nq_e * p_Nq_e];
        S_z_e[i][j] = S_z[ij + element * p_Nq_e * p_Nq_e];
        S_x_eT[j][i] = S_x_e[i][j];
        S_y_eT[j][i] = S_y_e[i][j];
        S_z_eT[j][i] = S_z_e[i][j];
      }
    }
    #pragma unroll
    for (int k = 0; k < p_Nq_e; k++) {
      #pragma unroll
      for (int j = 0; j < p_Nq_e; j++) {
        #pragma unroll
        for (int i = 0; i < p_Nq_e; i++) {
          pfloat value = 0.0;
          #pragma unroll
          for (int l = 0; l < p_Nq_e; l++)
            value += S_x_e[l][j] * work1(k,i,l,elem);
          work2[k][i][j] = value;
        }
      }
    }
    #pragma unroll
    for (int k = 0; k < p_Nq_e; k++) {
      #pragma unroll
      for (int j = 0; j < p_Nq_e; j++) {
        #pragma unroll
        for (int i = 0; i < p_Nq_e; i++) {
          pfloat value = 0.0;
          #pragma unroll
          for (int l = 0; l < p_Nq_e; l++)
            value += S_y_e[l][j] * work2[k][l][i];
          //work1(j,i,k,elem) = value;
          tmp[j][k][i] = value;
        }
      }
    }

    #pragma unroll
    for (int k = 0; k < p_Nq_e; k++) {
      #pragma unroll
      for (int j = 0; j < p_Nq_e; j++) {
        #pragma unroll
        for (int i = 0; i < p_Nq_e; i++) {
          const int v = i + j * p_Nq_e + k * p_Nq_e * p_Nq_e;
          pfloat value = 0.0;
          #pragma unroll
          for (int l = 0; l < p_Nq_e; l++)
            value += S_z_e[l][k] * tmp[j][l][i];
          work2[k][i][j] = value * inv_L[v + element * p_Nq_e * p_Nq_e * p_Nq_e];
        }
      }
    }

    #pragma unroll
    for (int k = 0; k < p_Nq_e; k++) {
      #pragma unroll
      for (int j = 0; j < p_Nq_e; j++) {
        #pragma unroll
        for (int i = 0; i < p_Nq_e; i++) {
          pfloat value = 0.0;
          #pragma unroll
          for (int l = 0; l < p_Nq_e; l++)
            value += S_x_eT[l][i] * work2[k][l][j];
          tmp[k][j][i] = value;
        }
      }
    }

    #pragma unroll
    for (int k = 0; k < p_Nq_e; k++) {
      #pragma unroll
      for (int j = 0; j < p_Nq_e; j++) {
        #pragma unroll
        for (int i = 0; i < p_Nq_e; i++) {
          pfloat value = 0.0;
          #pragma unroll
          for (int l = 0; l < p_Nq_e; l++)
            value += S_y_eT[l][j] * tmp[k][l][i];
          work2[j][k][i] = value;
        }
      }
    }
    #pragma unroll
    for (int k = 0; k < p_Nq_e; k++) {
      #pragma unroll
      for (int j = 0; j < p_Nq_e; j++) {
        #pragma unroll
        for (int i = 0; i < p_Nq_e; i++) {
          pfloat value = 0.0;
          #pragma unroll
          for (int l = 0; l < p_Nq_e; l++)
            value += S_z_eT[l][k] * work2[j][l][i];

#if (!p_restrict)
          const dlong elem_offset = element * p_Nq_e * p_Nq_e * p_Nq_e;
          const int v = i + j * p_Nq_e + k * p_Nq_e * p_Nq_e + elem_offset;
          Su[v] = value;
#endif
          tmp[k][j][i] = value;
        }
      }
    }
#if (!p_restrict)
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
    #pragma unroll
    for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 0;
          const int l2 = 0;
          work2[l1][j][k] = tmp[l2][j][k];
          work2[p_Nq_e - l1 - 1][j][k] = tmp[p_Nq_e - l2 - 1][j][k];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 0;
          work2[i][l1][k] = tmp[i][l2][k];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int k = 1; k < p_Nq_e-1; ++k){
          const int l1 = 0;
          const int l2 = 0;
          work2[i][p_Nq_e - l1 - 1][k] = tmp[i][p_Nq_e - l2 - 1][k];
      }
    }

    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 0;
          const int l2 = 0;
          work2[i][j][l1] = tmp[i][j][l2];
      }
    }
    #pragma unroll
    for(int i = 1; i < p_Nq_e-1; ++i){
      #pragma unroll
      for(int j = 1; j < p_Nq_e-1; ++j){
          const int l1 = 0;
          const int l2 = 0;
          work2[i][j][p_Nq_e - l1 - 1] = tmp[i][j][p_Nq_e - l2 - 1];
      }
    }

    #pragma unroll
    for(int k = 0; k < p_Nq_e; ++k){
      #pragma unroll
      for(int j = 0; j < p_Nq_e; ++j){
        #pragma unroll
        for(int i = 0; i < p_Nq_e; ++i) {
          const dlong elem_offset = element * p_Nq_e * p_Nq_e * p_Nq_e;
          const dlong idx = i + j * p_Nq_e + k * p_Nq_e * p_Nq_e + elem_offset;
          u[idx] = work2[k][j][i];
        }
      }
    }
#else  /* if (!p_restrict) */
    #pragma unroll
    for(int k = 0; k < p_Nq; ++k){
      #pragma unroll
      for(int j = 0; j < p_Nq; ++j){
        #pragma unroll
        for(int i = 0; i < p_Nq; ++i){
            const dlong elem_offset = element * p_Nq * p_Nq * p_Nq;
            const dlong idx = i + j * p_Nq + k * p_Nq * p_Nq + elem_offset;
            Su[idx] = tmp[k + 1][j + 1][i + 1] * wts[idx];
        }
      }
    }

#endif
  }
#undef getIdx
#undef work1
}
