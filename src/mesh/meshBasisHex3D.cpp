/*

   The MIT License (MIT)

   Copyright (c) 2020 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "nrssys.hpp"
#include "mesh.h"
#include "mesh3D.h"

// ------------------------------------------------------------------------
// HEX 3D NODES
// ------------------------------------------------------------------------
void NodesHex3D(int _N, dfloat* _r, dfloat* _s, dfloat* _t)
{
  int _Nq = _N + 1;

  dfloat* r1D = (dfloat*) malloc(_Nq * sizeof(dfloat));
  JacobiGLL(_N, r1D); //Gauss-Legendre-Lobatto nodes

  //Tensor product
  for (int k = 0; k < _Nq; k++)
    for (int j = 0; j < _Nq; j++)
      for (int i = 0; i < _Nq; i++) {
        _r[i + j * _Nq + k * _Nq * _Nq] = r1D[i];
        _s[i + j * _Nq + k * _Nq * _Nq] = r1D[j];
        _t[i + j * _Nq + k * _Nq * _Nq] = r1D[k];
      }

  free(r1D);
}

void FaceNodesHex3D(int _N, dfloat* _r, dfloat* _s, dfloat* _t, int* _faceNodes)
{
  int _Nq = _N + 1;
  int _Nfp = _Nq * _Nq;
  int _Np = _Nq * _Nq * _Nq;

  int cnt[6];
  for (int i = 0; i < 6; i++) cnt[i] = 0;

  dfloat deps = 1.;
  while((1. + deps) > 1.)
    deps *= 0.5;

  const dfloat NODETOL = 1000. * deps;

  for (int n = 0; n < _Np; n++) {
    if(fabs(_t[n] + 1) < NODETOL)
      _faceNodes[0 * _Nfp + (cnt[0]++)] = n;
    if(fabs(_s[n] + 1) < NODETOL)
      _faceNodes[1 * _Nfp + (cnt[1]++)] = n;
    if(fabs(_r[n] - 1) < NODETOL)
      _faceNodes[2 * _Nfp + (cnt[2]++)] = n;
    if(fabs(_s[n] - 1) < NODETOL)
      _faceNodes[3 * _Nfp + (cnt[3]++)] = n;
    if(fabs(_r[n] + 1) < NODETOL)
      _faceNodes[4 * _Nfp + (cnt[4]++)] = n;
    if(fabs(_t[n] - 1) < NODETOL)
      _faceNodes[5 * _Nfp + (cnt[5]++)] = n;
  }
}

#if 0

void mesh_t::VertexNodesHex3D(int _N, dfloat* _r, dfloat* _s, dfloat* _t, int* _vertexNodes)
{
  int _Nq = _N + 1;
  int _Np = _Nq * _Nq * _Nq;

  dfloat deps = 1.;
  while((1. + deps) > 1.)
    deps *= 0.5;

  const dfloat NODETOL = 1000. * deps;

  for(int n = 0; n < _Np; ++n) {
    if( (_r[n] + 1) * (_r[n] + 1) + (_s[n] + 1) * (_s[n] + 1) + (_t[n] + 1) * (_t[n] + 1) < NODETOL)
      _vertexNodes[0] = n;
    if( (_r[n] - 1) * (_r[n] - 1) + (_s[n] + 1) * (_s[n] + 1) + (_t[n] + 1) * (_t[n] + 1) < NODETOL)
      _vertexNodes[1] = n;
    if( (_r[n] - 1) * (_r[n] - 1) + (_s[n] - 1) * (_s[n] - 1) + (_t[n] + 1) * (_t[n] + 1) < NODETOL)
      _vertexNodes[2] = n;
    if( (_r[n] + 1) * (_r[n] + 1) + (_s[n] - 1) * (_s[n] - 1) + (_t[n] + 1) * (_t[n] + 1) < NODETOL)
      _vertexNodes[3] = n;
    if( (_r[n] + 1) * (_r[n] + 1) + (_s[n] + 1) * (_s[n] + 1) + (_t[n] - 1) * (_t[n] - 1) < NODETOL)
      _vertexNodes[4] = n;
    if( (_r[n] - 1) * (_r[n] - 1) + (_s[n] + 1) * (_s[n] + 1) + (_t[n] - 1) * (_t[n] - 1) < NODETOL)
      _vertexNodes[5] = n;
    if( (_r[n] - 1) * (_r[n] - 1) + (_s[n] - 1) * (_s[n] - 1) + (_t[n] - 1) * (_t[n] - 1) < NODETOL)
      _vertexNodes[6] = n;
    if( (_r[n] + 1) * (_r[n] + 1) + (_s[n] - 1) * (_s[n] - 1) + (_t[n] - 1) * (_t[n] - 1) < NODETOL)
      _vertexNodes[7] = n;
  }
}

void mesh_t::EquispacedNodesHex3D(int _N, dfloat* _r, dfloat* _s, dfloat* _t)
{
  int _Nq = _N + 1;

  //Equispaced 1D nodes
  dfloat* r1D = (dfloat*) malloc(_Nq * sizeof(dfloat));
  dfloat dr = 2.0 / _N;
  for (int i = 0; i < _Nq; i++) r1D[i] = -1.0 + i * dr;

  //Tensor product
  for (int k = 0; k < _Nq; k++)
    for (int j = 0; j < _Nq; j++)
      for (int i = 0; i < _Nq; i++) {
        _r[i + j * _Nq + k * _Nq * _Nq] = r1D[i];
        _s[i + j * _Nq + k * _Nq * _Nq] = r1D[j];
        _t[i + j * _Nq + k * _Nq * _Nq] = r1D[k];
      }

  free(r1D);
}

void mesh_t::EquispacedEToVHex3D(int _N, int* _EToV)
{
  int _Nq = _N + 1;
  int _Nverts = 4;

  //Tensor product
  int cnt = 0;
  for (int k = 0; k < _N; k++)
    for (int j = 0; j < _N; j++)
      for (int i = 0; i < _N; i++) {
        //tet 1 (0,3,2,7)
        _EToV[cnt * _Nverts + 0] = i  + (j  ) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 1] = i + 1 + (j + 1) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 2] = i  + (j + 1) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 3] = i + 1 + (j + 1) * _Nq + (k + 1) * _Nq * _Nq;
        cnt++;

        //tet 2 (0,1,3,7)
        _EToV[cnt * _Nverts + 0] = i  + (j  ) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 1] = i + 1 + (j  ) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 2] = i + 1 + (j + 1) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 3] = i + 1 + (j + 1) * _Nq + (k + 1) * _Nq * _Nq;
        cnt++;

        //tet 3 (0,2,6,7)
        _EToV[cnt * _Nverts + 0] = i  + (j  ) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 1] = i  + (j + 1) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 2] = i  + (j + 1) * _Nq + (k + 1) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 3] = i + 1 + (j + 1) * _Nq + (k + 1) * _Nq * _Nq;
        cnt++;

        //tet 4 (0,6,4,7)
        _EToV[cnt * _Nverts + 0] = i  + (j  ) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 1] = i  + (j + 1) * _Nq + (k + 1) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 2] = i  + (j  ) * _Nq + (k + 1) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 3] = i + 1 + (j + 1) * _Nq + (k + 1) * _Nq * _Nq;
        cnt++;

        //tet 5 (0,5,1,7)
        _EToV[cnt * _Nverts + 0] = i  + (j  ) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 1] = i + 1 + (j  ) * _Nq + (k + 1) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 2] = i + 1 + (j  ) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 3] = i + 1 + (j + 1) * _Nq + (k + 1) * _Nq * _Nq;
        cnt++;

        //tet 6 (0,4,5,7)
        _EToV[cnt * _Nverts + 0] = i  + (j  ) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 1] = i  + (j  ) * _Nq + (k + 1) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 2] = i + 1 + (j  ) * _Nq + (k + 1) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 3] = i + 1 + (j + 1) * _Nq + (k + 1) * _Nq * _Nq;
        cnt++;
      }
}

void mesh_t::SEMFEMEToVHex3D(int _N, int* _EToV)
{
  int _Nq = _N + 1;
  int _Nverts = 8;

  //Tensor product
  int cnt = 0;
  for (int k = 0; k < _N; k++)
    for (int j = 0; j < _N; j++)
      for (int i = 0; i < _N; i++) {
        _EToV[cnt * _Nverts + 0] = i  + (j  ) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 1] = i + 1 + (j  ) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 2] = i + 1 + (j + 1) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 3] = i  + (j + 1) * _Nq + (k  ) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 4] = i  + (j  ) * _Nq + (k + 1) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 5] = i + 1 + (j  ) * _Nq + (k + 1) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 6] = i + 1 + (j + 1) * _Nq + (k + 1) * _Nq * _Nq;
        _EToV[cnt * _Nverts + 7] = i  + (j + 1) * _Nq + (k + 1) * _Nq * _Nq;
        cnt++;
      }
}

// ------------------------------------------------------------------------
// ORTHONORMAL BASIS POLYNOMIALS
// ------------------------------------------------------------------------
void mesh_t::OrthonormalBasisHex3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat* P)
{
  *P = JacobiP(a,0,0,i) * JacobiP(b,0,0,j) * JacobiP(c,0,0,k);
}

void mesh_t::GradOrthonormalBasisHex3D(dfloat a,
                                       dfloat b,
                                       dfloat c,
                                       int i,
                                       int j,
                                       int k,
                                       dfloat* Pr,
                                       dfloat* Ps,
                                       dfloat* Pt)
{
  *Pr = GradJacobiP(a,0,0,i) * JacobiP(b,0,0,j) * JacobiP(c,0,0,k);
  *Ps = JacobiP(a,0,0,i) * GradJacobiP(b,0,0,j) * JacobiP(c,0,0,k);
  *Pt = JacobiP(a,0,0,i) * JacobiP(b,0,0,j) * GradJacobiP(c,0,0,k);
}

// ------------------------------------------------------------------------
// 2D VANDERMONDE MATRICES
// ------------------------------------------------------------------------

void mesh_t::VandermondeHex3D(int _N, int Npoints, dfloat* _r, dfloat* _s, dfloat* _t, dfloat* V)
{
  int _Nq = _N + 1;
  int _Np = _Nq * _Nq * _Nq;

  for(int n = 0; n < Npoints; n++)
    for(int k = 0; k < _Nq; k++)
      for(int j = 0; j < _Nq; j++)
        for(int i = 0; i < _Nq; i++) {
          int id = n * _Np + i + j * _Nq + k * _Nq * _Nq;
          OrthonormalBasisHex3D(_r[n], _s[n], _t[n], i, j, k, V + id);
        }
}

void mesh_t::GradVandermondeHex3D(int _N,
                                  int Npoints,
                                  dfloat* _r,
                                  dfloat* _s,
                                  dfloat* _t,
                                  dfloat* Vr,
                                  dfloat* Vs,
                                  dfloat* Vt)
{
  int _Nq = _N + 1;
  int _Np = _Nq * _Nq * _Nq;

  for(int n = 0; n < Npoints; n++)
    for(int k = 0; k < _Nq; k++)
      for(int j = 0; j < _Nq; j++)
        for(int i = 0; i < _Nq; i++) {
          int id = n * _Np + i + j * _Nq + k * _Nq * _Nq;
          GradOrthonormalBasisHex3D(_r[n], _s[n], _t[n], i, j, k, Vr + id, Vs + id, Vt + id);
        }
}

// ------------------------------------------------------------------------
// 2D OPERATOR MATRICES
// ------------------------------------------------------------------------
void mesh_t::MassMatrixHex3D(int _Np, dfloat* V, dfloat* _MM)
{
  // masMatrix = inv(V')*inv(V) = inv(V*V')
  for(int n = 0; n < _Np; ++n)
    for(int m = 0; m < _Np; ++m) {
      dfloat res = 0;
      for(int i = 0; i < _Np; ++i)
        res += V[n * _Np + i] * V[m * _Np + i];
      _MM[n * _Np + m] = res;
    }
  matrixInverse(_Np, _MM);
}

void mesh_t::LumpedMassMatrixHex3D(int _N, dfloat* _gllw, dfloat* _MM)
{
  int _Nq = _N + 1;
  int _Np = _Nq * _Nq * _Nq;

  // LumpedMassMatrix = gllw \ctimes gllw \ctimes gllw
  for(int k = 0; k < _Nq; ++k)
    for(int n = 0; n < _Nq; ++n)
      for(int m = 0; m < _Nq; ++m) {
        int id = n + m * _Nq + k * _Nq * _Nq;
        _MM[id + id * _Np] = _gllw[n] * _gllw[m] * _gllw[k];
      }
}

void mesh_t::invLumpedMassMatrixHex3D(int _N, dfloat* _gllw, dfloat* _invMM)
{
  int _Nq = _N + 1;
  int _Np = _Nq * _Nq * _Nq;

  // invLumpedMassMatrix = invgllw \ctimes invgllw
  for(int k = 0; k < _Nq; ++k)
    for(int n = 0; n < _Nq; ++n)
      for(int m = 0; m < _Nq; ++m) {
        int id = n + m * _Nq + k * _Nq * _Nq;
        _invMM[id + id * _Np] = 1.0 / (_gllw[n] * _gllw[m] * _gllw[k]);
      }
}

void mesh_t::DmatrixHex3D(int _N, int Npoints, dfloat* _r, dfloat* _s, dfloat* _t,
                          dfloat* _Dr, dfloat* _Ds, dfloat* _Dt)
{
  int _Nq = _N + 1;
  int _Np = _Nq * _Nq * _Nq;

  dfloat* V  = (dfloat*) calloc(Npoints * _Np, sizeof(dfloat));
  dfloat* Vr = (dfloat*) calloc(Npoints * _Np, sizeof(dfloat));
  dfloat* Vs = (dfloat*) calloc(Npoints * _Np, sizeof(dfloat));
  dfloat* Vt = (dfloat*) calloc(Npoints * _Np, sizeof(dfloat));

  VandermondeHex3D(_N, Npoints, _r, _s, _t, V);
  GradVandermondeHex3D(_N, Npoints, _r, _s, _t, Vr, Vs, Vt);

  //Dr = Vr/V, Ds = Vs/V, Dt = Vt/V
  matrixRightSolve(_Np, _Np, Vr, _Np, _Np, V, _Dr);
  matrixRightSolve(_Np, _Np, Vs, _Np, _Np, V, _Ds);
  matrixRightSolve(_Np, _Np, Vt, _Np, _Np, V, _Dt);

  free(V);
  free(Vr);
  free(Vs);
  free(Vt);
}

void mesh_t::InterpolationMatrixHex3D(int _N,
                                      int NpointsIn, dfloat* rIn, dfloat* sIn, dfloat* tIn,
                                      int NpointsOut, dfloat* rOut, dfloat* sOut, dfloat* tOut,
                                      dfloat* I)
{
  int _Nq = _N + 1;
  int _Np = _Nq * _Nq * _Nq;

  // need NpointsIn = _Np
  nrsCheck(NpointsIn != _Np, MPI_COMM_SELF, EXIT_FAILURE, 
           "Invalid Interplation operator requested", "");

  dfloat* VIn = (dfloat*) malloc(NpointsIn * _Np * sizeof(dfloat));
  dfloat* VOut = (dfloat*) malloc(NpointsOut * _Np * sizeof(dfloat));

  VandermondeHex3D(_N, NpointsIn,   rIn, sIn, tIn, VIn);
  VandermondeHex3D(_N, NpointsOut, rOut, sOut, tOut, VOut);

  matrixRightSolve(NpointsOut, _Np, VOut, NpointsIn, _Np, VIn, I);

  free(VIn);
  free(VOut);
}

#endif
