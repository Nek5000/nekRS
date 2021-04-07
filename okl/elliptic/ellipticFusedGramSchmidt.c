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

extern "C" void fusedGramSchmidt(const dlong & Nblocks, const dlong & N,
                                    const dlong & offset,
                                    const dfloat & h,
                                    const dfloat* __restrict__ weights,
                                    const dfloat* __restrict__ vk,
                                    const dfloat* __restrict__ vkp1,
                                    dfloat* __restrict__ w,
                                    dfloat* __restrict__ reduction)
{
  dfloat rdotr = 0.0;
  #pragma omp parallel for collapse(2)
  for(int fld = 0 ; fld < p_eNfields; fld++){
    for(int id = 0 ; id < N; ++id){
      const dfloat w_curr = w[id + fld * offset];
      const dfloat v_curr = vk[id + fld * offset];
      const dfloat w_new = w_curr - h * v_curr;
      if(!p_lastIter){
        rdotr += w_new * vkp1[id + fld * offset] * weights[id];
      } else {
        rdotr += w_new * w_new * weights[id];
      }
      w[id + fld * offset] = w_new;
    }
  }

  reduction[0] = rdotr;
}
