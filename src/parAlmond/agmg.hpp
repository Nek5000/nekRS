/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#ifndef PARALMOND_AGMGLEVEL_HPP
#define PARALMOND_AGMGLEVEL_HPP

namespace parAlmond {

class agmgLevel: public multigridLevel {

public:
  parCSR   *A,   *P,   *R;
  parHYB *o_A, *o_P, *o_R;

  SmoothType stype;
  dfloat lambda, lambda1, lambda0; //smoothing params

  int ChebyshevIterations;

  bool gatherLevel;
  ogs_t *ogs;
  dfloat *Gx, *Sx;
  occa::memory o_Sx, o_Gx;

  agmgLevel(parCSR *AA, KrylovType Ktype);
  agmgLevel(parCSR *AA, parCSR *PP, parCSR *RR, KrylovType Ktype);
  ~agmgLevel();

  void Ax(dfloat        *x, dfloat        *Ax);
  void Ax(occa::memory o_x, occa::memory o_Ax);

  void smooth(dfloat        *rhs, dfloat        *x, bool x_is_zero);
  void smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero);

  void residual(dfloat        *rhs, dfloat        *x, dfloat        *res);
  void residual(occa::memory o_rhs, occa::memory o_x, occa::memory o_res);

  void coarsen(dfloat        *x, dfloat        *Cx);
  void coarsen(occa::memory o_x, occa::memory o_Cx);

  void prolongate(dfloat        *x, dfloat        *Px);
  void prolongate(occa::memory o_x, occa::memory o_Px);

  void smoothJacobi(dfloat *r, dfloat *x, const bool x_is_zero);
  void smoothDampedJacobi(dfloat *r, dfloat *x, const bool x_is_zero);
  void smoothChebyshev(dfloat *r, dfloat *x, const bool x_is_zero);

  void smoothJacobi(occa::memory o_r, occa::memory o_x, bool x_is_zero);
  void smoothDampedJacobi(occa::memory o_r, occa::memory o_x, bool x_is_zero);
  void smoothChebyshev(occa::memory o_r, occa::memory o_x, bool x_is_zero);

  void Report();
};


agmgLevel *coarsenAgmgLevel(agmgLevel *level, KrylovType ktype, setupAide options);

parCSR* strongGraph(parCSR *A);

void formAggregates(parCSR *A, parCSR *C,
                     hlong* FineToCoarse,
                     hlong* globalAggStarts,
                     setupAide options);

parCSR *constructProlongation(parCSR *A, hlong *FineToCoarse,
                            hlong *globalAggStarts, dfloat **nullCoarseA);

parCSR *transpose(parCSR *A);

parCSR *galerkinProd(parCSR *A, parCSR *P);



void setupAgmgSmoother(agmgLevel *level, SmoothType s, int ChebIterations);

void allocateAgmgVectors(agmgLevel *level, int k, int numLevels, CycleType ctype);

void syncAgmgToDevice(agmgLevel *level, int k, int numLevels, CycleType ctype);

}

#endif
