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

extern "C"
void ellipticAxHex3D(const dlong & Nelements,
                     const dfloat* __restrict__ ggeo,
                     const dfloat* __restrict__ D,
                     const dfloat* __restrict__ S,
                     const dfloat & lambda,
                     const dfloat* __restrict__ q,
                     dfloat* __restrict__ Aq )
{
  dfloat s_q[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqr[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqs[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqt[p_Nq][p_Nq][p_Nq];

  dfloat s_D[p_Nq][p_Nq];
  dfloat s_S[p_Nq][p_Nq];

  for(int j = 0; j < p_Nq; ++j)
    for(int i = 0; i < p_Nq; ++i) {
      s_D[j][i] = D[j * p_Nq + i];
      s_S[j][i] = S[j * p_Nq + i];
    }

  for(dlong e = 0; e < Nelements; ++e) {
    const dlong element = e;

    for(int k = 0; k < p_Nq; k++)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong base = i + j * p_Nq + k * p_Nq * p_Nq + element * p_Np;
          const dfloat qbase = q[base];
          s_q[k][j][i] = qbase;
        }

    for(int k = 0; k < p_Nq; ++k)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_G00 = ggeo[gbase + p_G00ID * p_Np];
          const dfloat r_G01 = ggeo[gbase + p_G01ID * p_Np];
          const dfloat r_G11 = ggeo[gbase + p_G11ID * p_Np];
          const dfloat r_G12 = ggeo[gbase + p_G12ID * p_Np];
          const dfloat r_G02 = ggeo[gbase + p_G02ID * p_Np];
          const dfloat r_G22 = ggeo[gbase + p_G22ID * p_Np];

          dfloat qr = 0.f;
          dfloat qs = 0.f;
          dfloat qt = 0.f;

          for(int m = 0; m < p_Nq; m++) {
            qr += s_S[m][i] * s_q[k][j][m];
            qs += s_S[m][j] * s_q[k][m][i];
            qt += s_S[m][k] * s_q[m][j][i];
          }

          dfloat Gqr = r_G00 * qr;
          Gqr += r_G01 * qs;
          Gqr += r_G02 * qt;

          dfloat Gqs = r_G01 * qr;
          Gqs += r_G11 * qs;
          Gqs += r_G12 * qt;

          dfloat Gqt = r_G02 * qr;
          Gqt += r_G12 * qs;
          Gqt += r_G22 * qt;

          s_Gqr[k][j][i] = Gqr;
          s_Gqs[k][j][i] = Gqs;
          s_Gqt[k][j][i] = Gqt;
        }

    for(int k = 0; k < p_Nq; k++)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_GwJ = ggeo[gbase + p_GWJID * p_Np];

          dfloat r_Aq = r_GwJ * lambda * s_q[k][j][i];
          dfloat r_Aqr = 0, r_Aqs = 0, r_Aqt = 0;

          for(int m = 0; m < p_Nq; m++)
            r_Aqr += s_D[m][i] * s_Gqr[k][j][m];
          for(int m = 0; m < p_Nq; m++)
            r_Aqs += s_D[m][j] * s_Gqs[k][m][i];
          for(int m = 0; m < p_Nq; m++)
            r_Aqt += s_D[m][k] * s_Gqt[m][j][i];

          const dlong id = element * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          Aq[id] = r_Aqr + r_Aqs + r_Aqt + r_Aq;
        }
  }
}
extern "C"
void ellipticAxVarHex3D(const dlong & Nelements,
                        const dlong & offset,
                        const dfloat* __restrict__ ggeo,
                        const dfloat* __restrict__ D,
                        const dfloat* __restrict__ S,
                        const dfloat* __restrict__ lambda,
                        const dfloat* __restrict__ q,
                        dfloat* __restrict__ Aq )
{
  dfloat s_q[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqr[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqs[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqt[p_Nq][p_Nq][p_Nq];

  dfloat s_D[p_Nq][p_Nq];
  dfloat s_S[p_Nq][p_Nq];

  for(int j = 0; j < p_Nq; ++j)
    for(int i = 0; i < p_Nq; ++i) {
      s_D[j][i] = D[j * p_Nq + i];
      s_S[j][i] = S[j * p_Nq + i];
    }


  for(dlong e = 0; e < Nelements; ++e) {
    const dlong element = e;

    for(int k = 0; k < p_Nq; k++)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong base = i + j * p_Nq + k * p_Nq * p_Nq + element * p_Np;
          const dfloat qbase = q[base];
          s_q[k][j][i] = qbase;
        }

    for(int k = 0; k < p_Nq; ++k)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_G00 = ggeo[gbase + p_G00ID * p_Np];
          const dfloat r_G01 = ggeo[gbase + p_G01ID * p_Np];
          const dfloat r_G11 = ggeo[gbase + p_G11ID * p_Np];
          const dfloat r_G12 = ggeo[gbase + p_G12ID * p_Np];
          const dfloat r_G02 = ggeo[gbase + p_G02ID * p_Np];
          const dfloat r_G22 = ggeo[gbase + p_G22ID * p_Np];

          const dlong id = element * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_lam0 = lambda[id + 0 * offset];

          dfloat qr = 0.f;
          dfloat qs = 0.f;
          dfloat qt = 0.f;

          for(int m = 0; m < p_Nq; m++) {
            qr += s_S[m][i] * s_q[k][j][m];
            qs += s_S[m][j] * s_q[k][m][i];
            qt += s_S[m][k] * s_q[m][j][i];
          }

          dfloat Gqr = r_G00 * qr;
          Gqr += r_G01 * qs;
          Gqr += r_G02 * qt;

          dfloat Gqs = r_G01 * qr;
          Gqs += r_G11 * qs;
          Gqs += r_G12 * qt;

          dfloat Gqt = r_G02 * qr;
          Gqt += r_G12 * qs;
          Gqt += r_G22 * qt;

          s_Gqr[k][j][i] = r_lam0 * Gqr;
          s_Gqs[k][j][i] = r_lam0 * Gqs;
          s_Gqt[k][j][i] = r_lam0 * Gqt;
        }

    for(int k = 0; k < p_Nq; k++)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_GwJ = ggeo[gbase + p_GWJID * p_Np];

          const dlong id = element * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_lam1 = lambda[id + 1 * offset];

          dfloat r_Aq = r_GwJ * r_lam1 * s_q[k][j][i];
          dfloat r_Aqr = 0, r_Aqs = 0, r_Aqt = 0;

          for(int m = 0; m < p_Nq; m++)
            r_Aqr += s_D[m][i] * s_Gqr[k][j][m];
          for(int m = 0; m < p_Nq; m++)
            r_Aqs += s_D[m][j] * s_Gqs[k][m][i];
          for(int m = 0; m < p_Nq; m++)
            r_Aqt += s_D[m][k] * s_Gqt[m][j][i];

          Aq[id] = r_Aqr + r_Aqs + r_Aqt + r_Aq;
        }
  }
}

extern "C"
void ellipticBlockAxVarHex3D_N3(const dlong & Nelements,
                                const dlong & offset,
                                const dlong & loffset,
                                const dfloat* __restrict__ ggeo,
                                const dfloat* __restrict__ D,
                                const dfloat* __restrict__ S,
                                const dfloat* __restrict__ lambda,
                                const dfloat* __restrict__ q,
                                dfloat* __restrict__ Aq )
{
  dfloat s_q[3][p_Nq][p_Nq][p_Nq];
  dfloat s_Gqr[3][p_Nq][p_Nq][p_Nq];
  dfloat s_Gqs[3][p_Nq][p_Nq][p_Nq];
  dfloat s_Gqt[3][p_Nq][p_Nq][p_Nq];

  dfloat s_D[p_Nq][p_Nq];
  dfloat s_S[p_Nq][p_Nq];

  for(int j = 0; j < p_Nq; ++j)
    for(int i = 0; i < p_Nq; ++i) {
      s_D[j][i] = D[j * p_Nq + i];
      s_S[j][i] = S[j * p_Nq + i];
    }

  for(dlong e = 0; e < Nelements; ++e) {
    const dlong element = e;

    for(int k = 0; k < p_Nq; k++)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong base = i + j * p_Nq + k * p_Nq * p_Nq + element * p_Np;
          s_q[0][k][j][i] = q[base + 0 * offset];
          s_q[1][k][j][i] = q[base + 1 * offset];
          s_q[2][k][j][i] = q[base + 2 * offset];
        }

    for(int k = 0; k < p_Nq; ++k)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_G00 = ggeo[gbase + p_G00ID * p_Np];
          const dfloat r_G01 = ggeo[gbase + p_G01ID * p_Np];
          const dfloat r_G11 = ggeo[gbase + p_G11ID * p_Np];
          const dfloat r_G12 = ggeo[gbase + p_G12ID * p_Np];
          const dfloat r_G02 = ggeo[gbase + p_G02ID * p_Np];
          const dfloat r_G22 = ggeo[gbase + p_G22ID * p_Np];

          const dlong id      = element * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;

          const dfloat r_lam00 = lambda[id + 0 * offset + 0 * loffset];
          const dfloat r_lam10 = lambda[id + 0 * offset + 1 * loffset];
          const dfloat r_lam20 = lambda[id + 0 * offset + 2 * loffset];

          dfloat qr0 = 0.f, qr1 = 0.f, qr2 = 0.f;
          dfloat qs0 = 0.f, qs1 = 0.f, qs2 = 0.f;
          dfloat qt0 = 0.f, qt1 = 0.f, qt2 = 0.f;

          for(int m = 0; m < p_Nq; m++) {
            qr0 += s_S[m][i] * s_q[0][k][j][m];
            qs0 += s_S[m][j] * s_q[0][k][m][i];
            qt0 += s_S[m][k] * s_q[0][m][j][i];
            //
            qr1 += s_S[m][i] * s_q[1][k][j][m];
            qs1 += s_S[m][j] * s_q[1][k][m][i];
            qt1 += s_S[m][k] * s_q[1][m][j][i];

            qr2 += s_S[m][i] * s_q[2][k][j][m];
            qs2 += s_S[m][j] * s_q[2][k][m][i];
            qt2 += s_S[m][k] * s_q[2][m][j][i];
          }

          dfloat Gqr0 = r_G00 * qr0 + r_G01 * qs0 + r_G02 * qt0;
          dfloat Gqs0 = r_G01 * qr0 + r_G11 * qs0 + r_G12 * qt0;
          dfloat Gqt0 = r_G02 * qr0 + r_G12 * qs0 + r_G22 * qt0;

          dfloat Gqr1 = r_G00 * qr1 + r_G01 * qs1 + r_G02 * qt1;
          dfloat Gqs1 = r_G01 * qr1 + r_G11 * qs1 + r_G12 * qt1;
          dfloat Gqt1 = r_G02 * qr1 + r_G12 * qs1 + r_G22 * qt1;

          dfloat Gqr2 = r_G00 * qr2 + r_G01 * qs2 + r_G02 * qt2;
          dfloat Gqs2 = r_G01 * qr2 + r_G11 * qs2 + r_G12 * qt2;
          dfloat Gqt2 = r_G02 * qr2 + r_G12 * qs2 + r_G22 * qt2;

          s_Gqr[0][k][j][i] = r_lam00 * Gqr0;
          s_Gqs[0][k][j][i] = r_lam00 * Gqs0;
          s_Gqt[0][k][j][i] = r_lam00 * Gqt0;

          s_Gqr[1][k][j][i] = r_lam10 * Gqr1;
          s_Gqs[1][k][j][i] = r_lam10 * Gqs1;
          s_Gqt[1][k][j][i] = r_lam10 * Gqt1;

          s_Gqr[2][k][j][i] = r_lam20 * Gqr2;
          s_Gqs[2][k][j][i] = r_lam20 * Gqs2;
          s_Gqt[2][k][j][i] = r_lam20 * Gqt2;
        }

    for(int k = 0; k < p_Nq; k++)
      for(int j = 0; j < p_Nq; ++j)
        for(int i = 0; i < p_Nq; ++i) {
          const dlong gbase  = element * p_Nggeo * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;
          const dfloat r_GwJ = ggeo[gbase + p_GWJID * p_Np];

          const dlong id = element * p_Np + k * p_Nq * p_Nq + j * p_Nq + i;

          const dfloat r_lam01 = lambda[id + 1 * offset + 0 * loffset];
          const dfloat r_lam11 = lambda[id + 1 * offset + 1 * loffset];
          const dfloat r_lam21 = lambda[id + 1 * offset + 2 * loffset];

          dfloat r_Aq0 = r_GwJ * r_lam01 * s_q[0][k][j][i];
          dfloat r_Aq1 = r_GwJ * r_lam11 * s_q[1][k][j][i];
          dfloat r_Aq2 = r_GwJ * r_lam21 * s_q[2][k][j][i];

          dfloat r_Aqr0 = 0.f, r_Aqs0 = 0.f, r_Aqt0 = 0.f;
          dfloat r_Aqr1 = 0.f, r_Aqs1 = 0.f, r_Aqt1 = 0.f;
          dfloat r_Aqr2 = 0.f, r_Aqs2 = 0.f, r_Aqt2 = 0.f;

          for(int m = 0; m < p_Nq; m++) {
            r_Aqr0 += s_D[m][i] * s_Gqr[0][k][j][m];
            r_Aqr1 += s_D[m][i] * s_Gqr[1][k][j][m];
            r_Aqr2 += s_D[m][i] * s_Gqr[2][k][j][m];
          }

          for(int m = 0; m < p_Nq; m++) {
            r_Aqs0 += s_D[m][j] * s_Gqs[0][k][m][i];
            r_Aqs1 += s_D[m][j] * s_Gqs[1][k][m][i];
            r_Aqs2 += s_D[m][j] * s_Gqs[2][k][m][i];
          }

          for(int m = 0; m < p_Nq; m++) {
            r_Aqt0 += s_D[m][k] * s_Gqt[0][m][j][i];
            r_Aqt1 += s_D[m][k] * s_Gqt[1][m][j][i];
            r_Aqt2 += s_D[m][k] * s_Gqt[2][m][j][i];
          }

          Aq[id + 0 * offset] = r_Aqr0 + r_Aqs0 + r_Aqt0 + r_Aq0;
          Aq[id + 1 * offset] = r_Aqr1 + r_Aqs1 + r_Aqt1 + r_Aq1;
          Aq[id + 2 * offset] = r_Aqr2 + r_Aqs2 + r_Aqt2 + r_Aq2;
        }
  }
}

//
extern "C"
void ellipticStressAxVarHex3D(const dlong &Nelements,
                                 const dlong &offset,
                                 const dlong &loffset,
                                 const dfloat * __restrict__ vgeo,
                                 const dfloat * __restrict__ D,
                                 const dfloat * __restrict__ S,
                                 const dfloat * __restrict__ lambda,
                                 const dfloat * __restrict__ q,
                                 dfloat * __restrict__ Aq){

    dfloat s_D[p_Nq][p_Nq];

    dfloat s_U[p_Nq][p_Nq][p_Nq];
    dfloat s_V[p_Nq][p_Nq][p_Nq];
    dfloat s_W[p_Nq][p_Nq][p_Nq];

    dfloat s_SUr[p_Nq][p_Nq][p_Nq];
    dfloat s_SUs[p_Nq][p_Nq][p_Nq];
    dfloat s_SUt[p_Nq][p_Nq][p_Nq];

    dfloat s_SVr[p_Nq][p_Nq][p_Nq];
    dfloat s_SVs[p_Nq][p_Nq][p_Nq];
    dfloat s_SVt[p_Nq][p_Nq][p_Nq];

    dfloat s_SWr[p_Nq][p_Nq][p_Nq];
    dfloat s_SWs[p_Nq][p_Nq][p_Nq];
    dfloat s_SWt[p_Nq][p_Nq][p_Nq];

    for(int j=0;j<p_Nq;++j){
      for(int i=0;i<p_Nq;++i){
      s_D[j][i] = D[j*p_Nq+i];
    }
  }
    

for(dlong e=0; e<Nelements; ++e){

    for(int k=0;k<p_Nq;++k){ 
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){
            const dlong id = e*p_Np+k*p_Nq*p_Nq+j*p_Nq+i;
            s_U[k][j][i] = q[id + 0*offset];
            s_V[k][j][i] = q[id + 1*offset];
            s_W[k][j][i] = q[id + 2*offset];
        }
      }
    }
    


    // loop over slabs
     for(int k=0;k<p_Nq;++k){ 
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){   
          const dlong gid = i + j*p_Nq + k*p_Nq*p_Nq + e*p_Np*p_Nvgeo;
          const dfloat rx = vgeo[gid + p_RXID*p_Np];
          const dfloat ry = vgeo[gid + p_RYID*p_Np];
          const dfloat rz = vgeo[gid + p_RZID*p_Np];
          
          const dfloat sx = vgeo[gid + p_SXID*p_Np];
          const dfloat sy = vgeo[gid + p_SYID*p_Np];
          const dfloat sz = vgeo[gid + p_SZID*p_Np];
          
          const dfloat tx = vgeo[gid + p_TXID*p_Np];
          const dfloat ty = vgeo[gid + p_TYID*p_Np];
          const dfloat tz = vgeo[gid + p_TZID*p_Np];
          
          const dfloat JW = vgeo[gid + p_JWID*p_Np];

          // compute 1D derivatives
          dfloat ur = 0.f, us = 0.f, ut = 0.f;
          dfloat vr = 0.f, vs = 0.f, vt = 0.f;
          dfloat wr = 0.f, ws = 0.f, wt = 0.f;
          for(int m=0;m<p_Nq;++m){
            const dfloat Dim = s_D[i][m]; // Dr
            const dfloat Djm = s_D[j][m]; // Ds
            const dfloat Dkm = s_D[k][m]; // Dt

            ur += Dim*s_U[k][j][m];
            us += Djm*s_U[k][m][i];
            ut += Dkm*s_U[m][j][i];
            //
            vr += Dim*s_V[k][j][m];
            vs += Djm*s_V[k][m][i];
            vt += Dkm*s_V[m][j][i];
            //
            wr += Dim*s_W[k][j][m];
            ws += Djm*s_W[k][m][i];
            wt += Dkm*s_W[m][j][i];
          }

          const dlong id = e*p_Np + k*p_Nq*p_Nq + j*p_Nq + i;  
          const dfloat u_lam0 = lambda[id + 0*offset + 0*loffset]; 
          // const dfloat u_lam1 = lambda[id + 1*offset + 0*loffset];
          const dfloat v_lam0 = lambda[id + 0*offset + 1*loffset]; 
          // const dfloat v_lam1 = lambda[id + 1*offset + 1*loffset];
          const dfloat w_lam0 = lambda[id + 0*offset + 2*loffset]; 
          // const dfloat w_lam1 = lambda[id + 1*offset + 2*loffset];

         
          const dfloat dudx = rx*ur + sx*us + tx*ut; 
          const dfloat dudy = ry*ur + sy*us + ty*ut; 
          const dfloat dudz = rz*ur + sz*us + tz*ut; 

          const dfloat dvdx = rx*vr + sx*vs + tx*vt; 
          const dfloat dvdy = ry*vr + sy*vs + ty*vt; 
          const dfloat dvdz = rz*vr + sz*vs + tz*vt; 

          const dfloat dwdx = rx*wr + sx*ws + tx*wt; 
          const dfloat dwdy = ry*wr + sy*ws + ty*wt; 
          const dfloat dwdz = rz*wr + sz*ws + tz*wt; 

          const dfloat s11 = u_lam0*JW*(dudx + dudx); 
          const dfloat s12 = u_lam0*JW*(dudy + dvdx); 
          const dfloat s13 = u_lam0*JW*(dudz + dwdx); 

          const dfloat s21 = v_lam0*JW*(dvdx + dudy); 
          const dfloat s22 = v_lam0*JW*(dvdy + dvdy); 
          const dfloat s23 = v_lam0*JW*(dvdz + dwdy); 

          const dfloat s31 = w_lam0*JW*(dwdx + dudz); 
          const dfloat s32 = w_lam0*JW*(dwdy + dvdz); 
          const dfloat s33 = w_lam0*JW*(dwdz + dwdz); 

          s_SUr[k][j][i] =  rx*s11 + ry*s12 + rz*s13;
          s_SUs[k][j][i] =  sx*s11 + sy*s12 + sz*s13;
          s_SUt[k][j][i] =  tx*s11 + ty*s12 + tz*s13;
          //
          s_SVr[k][j][i] =  rx*s21 + ry*s22 + rz*s23;
          s_SVs[k][j][i] =  sx*s21 + sy*s22 + sz*s23;
          s_SVt[k][j][i] =  tx*s21 + ty*s22 + tz*s23;
          //
          s_SWr[k][j][i] =  rx*s31 + ry*s32 + rz*s33;
          s_SWs[k][j][i] =  sx*s31 + sy*s32 + sz*s33;
          s_SWt[k][j][i] =  tx*s31 + ty*s32 + tz*s33;

         
        }
      }
    }


// loop over slabs
    for(int k=0;k<p_Nq;++k){ 
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){
          dfloat r_Au = 0.f, r_Av = 0.f, r_Aw = 0.f;
          for(int m = 0; m < p_Nq; m++) {
            const dfloat Dim = s_D[m][i]; // Dr'
            const dfloat Djm = s_D[m][j]; // Ds'
            const dfloat Dkm = s_D[m][k]; // Dt'

            r_Au += Dim*s_SUr[k][j][m];
            r_Au += Djm*s_SUs[k][m][i];
            r_Au += Dkm*s_SUt[m][j][i];

            r_Av += Dim*s_SVr[k][j][m];
            r_Av += Djm*s_SVs[k][m][i];
            r_Av += Dkm*s_SVt[m][j][i];

            r_Aw += Dim*s_SWr[k][j][m];
            r_Aw += Djm*s_SWs[k][m][i];
            r_Aw += Dkm*s_SWt[m][j][i];
          }
          const dlong id      = e*p_Np +k*p_Nq*p_Nq+ j*p_Nq + i;
          const dfloat u_lam1 = lambda[id + 1*offset + 0*loffset];
          const dfloat v_lam1 = lambda[id + 1*offset + 1*loffset];
          const dfloat w_lam1 = lambda[id + 1*offset + 2*loffset];
         
          const dlong gid = i + j*p_Nq + k*p_Nq*p_Nq + e*p_Np*p_Nvgeo;
          const dfloat JW = vgeo[gid + p_JWID*p_Np];
           // store in register
          Aq[id+0*offset] =  r_Au + u_lam1*JW*s_U[k][j][i]; 
          Aq[id+1*offset] =  r_Av + v_lam1*JW*s_V[k][j][i];
          Aq[id+2*offset] =  r_Aw + w_lam1*JW*s_W[k][j][i];
        }
      }
    }
  }
}

