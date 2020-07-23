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

typedef struct
{
  hlong row;
  hlong col;
  int ownerRank;
  dfloat val;
}nonZero_t;

typedef struct
{
  long long int preconBytes;

  ogs_t* ogs;
  ogs_t* FEMogs;

  dfloat* zP;
  occa::memory o_zP;

  dfloat* xG, * rhsG;
  occa::memory o_xG, o_rhsG;

  occa::memory o_Gr;
  occa::memory o_Gz;
  occa::memory o_Sr;

  occa::memory o_vmapPP;
  occa::memory o_faceNodesP;

  occa::memory o_oasForward;
  occa::memory o_oasBack;
  occa::memory o_oasDiagInvOp;

  occa::memory o_oasForwardDg;
  occa::memory o_oasBackDg;
  occa::memory o_oasDiagInvOpDg;
  occa::memory o_invDegreeDGP;

  occa::memory o_oasForwardDgT;
  occa::memory o_oasBackDgT;

  occa::kernel restrictKernel;

  occa::kernel coarsenKernel;
  occa::kernel prolongateKernel;

  occa::kernel overlappingPatchKernel;
  occa::kernel exactPatchSolverKernel;
  occa::kernel approxPatchSolverKernel;
  occa::kernel exactFacePatchSolverKernel;
  occa::kernel approxFacePatchSolverKernel;
  occa::kernel exactBlockJacobiSolverKernel;
  occa::kernel approxBlockJacobiSolverKernel;
  occa::kernel patchGatherKernel;
  occa::kernel facePatchGatherKernel;
  occa::kernel CGLocalPatchKernel;

  occa::memory o_rFEM;
  occa::memory o_zFEM;
  occa::memory o_GrFEM;
  occa::memory o_GzFEM;

  occa::kernel SEMFEMInterpKernel;
  occa::kernel SEMFEMAnterpKernel;

  ogs_t* ogsP, * ogsDg;

  occa::memory o_diagA;
  occa::memory o_invDiagA;
  occa::memory o_invAP;
  occa::memory o_invDegreeAP;
  occa::memory o_patchesIndex;

  // coarse grid basis for preconditioning
  occa::memory o_V1, o_Vr1, o_Vs1, o_Vt1;
  occa::memory o_r1, o_z1;
  dfloat* r1, * z1;

  void* xxt;

  occa::memory o_coarseInvDegree;

  int coarseNp;
  hlong coarseTotal;
  hlong* coarseOffsets;
  dfloat* B, * tmp2;
  occa::memory* o_B, o_tmp2;
  void* xxt2;
  parAlmond::solver_t* parAlmond;

  // block Jacobi precon
  occa::memory o_invMM;
  occa::kernel blockJacobiKernel;
  occa::kernel partialblockJacobiKernel;

  //dummy almond level to store the OAS smoothing op
  // agmgLevel *OASLevel;
  // void **OASsmoothArgs;

  //SEMFEM variables
  mesh_t* femMesh;

  // Overlapping Additive Schwarz variables
  void* ellipticOneRing;

  dlong NoneRingSendTotal;
  hlong* oneRingSendList;
  hlong* NoneRingSend;

  void* oneRingSendBuffer;

  MPI_Request* oneRingSendRequests;

  hlong NoneRingRecvTotal;
  hlong* NoneRingRecv;
  void* oneRingRecvBuffer;
  MPI_Request* oneRingRecvRequests;

  void* ellipticOasCoarse;
  occa::memory o_oasRestrictionMatrix; // Y
  occa::memory o_oasProlongationMatrix;// Y

  occa::memory o_oasCoarseTmp; // Y
  occa::memory o_oasFineTmp;   // Y

  occa::memory o_oneRingSendList;  // Y
  occa::memory o_oneRingSendBuffer;  //Y

  occa::memory o_oneRingRecvBuffer; // Y

  occa::kernel oasRestrictionKernel;
  occa::kernel oasProlongationKernel;

  ogs_t* oasOgs;

  bool additive;
} precon_t;
