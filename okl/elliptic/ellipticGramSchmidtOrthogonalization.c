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

// w = w - \sum_i^k (w^T v_i) v_i
// w = w - \sum_i^k y_i v_i
extern "C" void FUNC(gramSchmidtOrthogonalization)(const dlong & Nblock, const dlong & N,
                                    const dlong & offset,
                                    const dlong & gmresSize,
                                    const dfloat * __restrict__ weights,
                                    const dfloat*  __restrict__ y,
                                    const dfloat*  __restrict__ V,
                                    dfloat*  __restrict__ w,
                                    dfloat*  __restrict__ reduction)
{
  for(int j = 0; j < gmresSize; ++j){
    const dfloat yj = y[j];
    #pragma unroll
    for(int fld = 0; fld < p_Nfields; fld++){
      for(dlong n = 0; n < N; ++n) {
        const dfloat Vnj = V[n + fld * offset + j * offset * p_Nfields];
        w[n + fld * offset] -= yj * Vnj;
      }
    }
  }
  dfloat sum = 0.0;
  #pragma unroll
  for(int fld = 0; fld < p_Nfields; fld++){
    for(dlong n = 0 ; n < N; ++n){
      const dfloat weight = weights[n];
      const dfloat w_curr = w[n + fld * offset];
      sum += w_curr * w_curr * weight;
    }
  }
  reduction[0] = sum;
}
