/*

   The MIT License (MIT)

   Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

 */

extern "C" void FUNC(ellipticPreconProlongateHex3D)(const dlong& Nelements,
                                               const pfloat* __restrict__  R,
                                               const pfloat* __restrict__  qc,
                                               pfloat* __restrict__  qN)
{
  dfloat r_q[p_NqCoarse][p_NqFine][p_NqFine];

  dfloat s_q[p_NqCoarse][p_NqCoarse];
  dfloat s_Pq[p_NqFine][p_NqCoarse];

  dfloat s_R[p_NqCoarse][p_NqFine];
  for(int j = 0; j < p_NqCoarse; ++j){
    for(int i = 0; i < p_NqFine; ++i) {
      int t = i + j * p_NqFine;
      const dfloat r = R[t];
      s_R[j][i] = r;
    }
  }
#ifdef __NEKRS__OMP__
  #pragma omp parallel for private(s_Pq, r_q, s_q)
#endif
  for(dlong e = 0; e < Nelements; ++e) {


    for(int j = 0; j < p_NqCoarse; ++j)
      for(int i = 0; i < p_NqFine; ++i) {
        const int t = i + j * p_NqFine;
        if(t < p_NqCoarse * p_NqCoarse) {
          #pragma unroll
          for(int k = 0; k < p_NqFine; ++k)
            r_q[j][i][k] = 0;

          for(int k = 0; k < p_NqCoarse; ++k) {
            const int id = t + k * p_NqCoarse * p_NqCoarse + e * p_NpCoarse;
            const dfloat tmp = qc[id];

            #pragma unroll
            for(int m = 0; m < p_NqFine; ++m)
              r_q[j][i][m] += s_R[k][m] * tmp;
          }
        }
      }

    for(int k = 0; k < p_NqFine; ++k) {

      for(int j = 0; j < p_NqCoarse; ++j)
        for(int i = 0; i < p_NqFine; ++i) {
          const int t = i + j * p_NqFine;
          if(t < p_NqCoarse * p_NqCoarse) {
            const int ti = t % p_NqCoarse;
            const int tj = t / p_NqCoarse;

            s_q[tj][ti] = r_q[j][i][k];
          }
        }

      for(int j = 0; j < p_NqCoarse; ++j)
        for(int i = 0; i < p_NqFine; ++i) {
          const int t = i + j * p_NqFine;

          if(t < p_NqCoarse * p_NqFine) {
            const int ti = t % p_NqCoarse;
            const int tj = t / p_NqCoarse;

            dfloat res = 0;

            #pragma unroll
            for(int m = 0; m < p_NqCoarse; ++m)
              res += s_R[m][tj] * s_q[m][ti];

            s_Pq[tj][ti] = res;
          }
        }

      for(int j = 0; j < p_NqFine; ++j)
        for(int i = 0; i < p_NqFine; ++i) {
          dfloat res = 0;

          #pragma unroll
          for(int m = 0; m < p_NqCoarse; ++m)
            res += s_R[m][i] * s_Pq[j][m];

          const int id = i + j * p_NqFine + k * p_NqFine * p_NqFine + e * p_NpFine;
          qN[id] += res;
        }
    }
  }
}
