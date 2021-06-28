

extern "C" void fusedFDM(
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

  for (dlong my_elem = 0; my_elem < localNelements; ++my_elem) {
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
