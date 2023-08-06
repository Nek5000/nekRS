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

// r = b - Ax
// (r,r)
extern "C" void FUNC(fusedResidualAndNorm)(const dlong & Nblocks, const dlong & N,
                                    const dlong & offset,
                                    const dfloat* __restrict__ weights,
                                    const dfloat* __restrict__ b_vec,
                                    const dfloat* __restrict__ Ax,
                                    dfloat* __restrict__ r,
                                    dfloat* __restrict__ reduction)
{
  dfloat rdotr = 0.0;
  for(int fld = 0 ; fld < p_Nfields; ++fld){
    for(int id = 0 ; id < N; ++id){
      const dfloat rnew = b_vec[id + fld * offset] - Ax[id + fld * offset];
      r[id + fld * offset] = rnew;
      rdotr += rnew * rnew * weights[id];
    }
  }
  reduction[0] = rdotr;
}
