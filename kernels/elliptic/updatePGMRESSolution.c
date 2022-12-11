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

// x = x + Zy
extern "C" void FUNC(updatePGMRESSolution)(const dlong & N,
                                    const dlong & offset,
                                    const dlong & gmresSize,
                                    const dfloat* __restrict__ y,
                                    const dfloat*  __restrict__ Z,
                                    dfloat*  __restrict__ x)
{
#ifdef __NEKRS__OMP__
  #pragma omp parallel for collapse(3)
#endif
  for(int j = 0; j < gmresSize; ++j){
    #pragma unroll
    for(int fld = 0; fld < p_Nfields; ++fld){
      for(int n = 0 ; n < N; ++n){
        const dfloat yj = y[j];
        const dfloat Znj = Z[n + fld * offset + j * offset * p_Nfields];
        x[n + fld * offset] += Znj * yj;
      }
    }
  }
}
