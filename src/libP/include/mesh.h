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

#ifndef MESH_H
#define MESH_H 1

#include <unistd.h>

#include "mpi.h"
#include <math.h>
#include <stdlib.h>
#include <occa.hpp>

#include "types.h"
#include "ogs.hpp"

#include "timer.h"

#include "setupAide.hpp"

#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12

extern "C" { // Start C linkage
typedef struct {

  MPI_Comm comm;
  int rank, size; // MPI rank and size (process count)
  
  int dim;
  int Nverts, Nfaces, NfaceVertices;

  hlong Nnodes;
  dfloat *EX; // coordinates of vertices for each element
  dfloat *EY;
  dfloat *EZ;

  dlong Nelements;
  hlong *EToV; // element-to-vertex connectivity
  dlong *EToE; // element-to-element connectivity
  int   *EToF; // element-to-(local)face connectivity
  int   *EToP; // element-to-partition/process connectivity
  int   *EToB; // element-to-boundary condition type

  hlong *elementInfo; //type of element

  // boundary faces
  hlong NboundaryFaces; // number of boundary faces
  hlong *boundaryInfo; // list of boundary faces (type, vertex-1, vertex-2, vertex-3)

  // MPI halo exchange info
  dlong  totalHaloPairs;  // number of elements to be sent in halo exchange
  dlong *haloElementList; // sorted list of elements to be sent in halo exchange
  int *NhaloPairs;      // number of elements worth of data to send/recv
  int  NhaloMessages;     // number of messages to send

  dlong *haloGetNodeIds; // volume node ids of outgoing halo nodes
  dlong *haloPutNodeIds; // volume node ids of incoming halo nodes

  void *haloSendRequests;
  void *haloRecvRequests;

  dlong NinternalElements; // number of elements that can update without halo exchange
  dlong NnotInternalElements; // number of elements that cannot update without halo exchange

  // CG gather-scatter info
  hlong *globalIds;
  hlong *maskedGlobalIds;
  void *gsh, *hostGsh; // gslib struct pointer
  ogs_t *ogs; //occa gs pointer

  // list of elements that are needed for global gather-scatter
  dlong NglobalGatherElements;
  dlong *globalGatherElementList;
  occa::memory o_globalGatherElementList;

  // list of elements that are not needed for global gather-scatter
  dlong NlocalGatherElements;
  dlong *localGatherElementList;
  occa::memory o_localGatherElementList;

  //list of fair pairs
  dlong NfacePairs;
  dlong *EToFPairs;
  dlong *FPairsToE;
  int *FPairsToF;

  // NBN: streams / command queues
  occa::stream stream0, stream1;

  // volumeGeometricFactors;
  dlong Nvgeo;
  dfloat *vgeo;

  // second order volume geometric factors
  dlong Nggeo;
  dfloat *ggeo;

  // volume node info
  int N, Np;
  dfloat *r, *s, *t;    // coordinates of local nodes
  dfloat *Dr, *Ds, *Dt; // collocation differentiation matrices
  dfloat *Dmatrices;
  dfloat *MM, *invMM;           // reference mass matrix
  dfloat *LMM, *invLMM;
  dfloat *Srr,*Srs, *Srt; //element stiffness matrices
  dfloat *Ssr,*Sss, *Sst;
  dfloat *Str,*Sts, *Stt;
  dfloat *Smatrices;
  int maxNnzPerRow;
  dfloat *x, *y, *z;    // coordinates of physical nodes
  
  dfloat sphereRadius;  // for Quad3D 
  
  dfloat volume;

  // indices of vertex nodes
  int *vertexNodes;

  // quad specific quantity
  int Nq, NqP, NpP;
  
  dfloat *D; // 1D differentiation matrix (for tensor-product)
  dfloat *DW; // weak 1D differentiation matrix (for tensor-product)
  dfloat *gllz; // 1D GLL quadrature nodes
  dfloat *gllw; // 1D GLL quadrature weights

  int gjNq;
  dfloat *gjr,*gjw; // 1D nodes and weights for Gauss Jacobi quadature
  dfloat *gjI,*gjD; // 1D GLL to Gauss node interpolation and differentiation matrices
  dfloat *gjD2;     // 1D GJ to GJ node differentiation

  // transform to/from eigenmodes of 1D laplacian (with built in weighting)
  dfloat *oasForward;
  dfloat *oasBack;
  dfloat *oasDiagOp;

  // transform to/from eigenmode of IPDG 1D laplacian
  dfloat *oasForwardDg;
  dfloat *oasBackDg;
  dfloat *oasDiagOpDg;

  //rotated node ids
  int *rmapP;

  //reference patch inverse (for OAS precon)
  dfloat *invAP;

  // face node info
  int Nfp;        // number of nodes per face
  int *faceNodes; // list of element reference interpolation nodes on element faces
  dlong *vmapM;     // list of volume nodes that are face nodes
  dlong *vmapP;     // list of volume nodes that are paired with face nodes
  dlong *mapP;     // list of surface nodes that are paired with -ve surface  nodes
  int *faceVertices; // list of mesh vertices on each face

  dfloat *LIFT; // lift matrix
  dfloat *FMM;  // Face Mass Matrix
  dfloat *sMT; // surface mass (MM*LIFT)^T

  dlong   Nsgeo;
  dfloat *sgeo;

  // field info for PDE solver
  int Nfields;
  dfloat *q;    // solution data array
  dfloat *fQM, *fQP; //solution trace arrays
  dfloat *rhsq, *rhsq2, *rhsq3; // right hand side data array
  dfloat *resq; // residual data array (for LSERK time-stepping)

  dfloat Lambda2; // square of penalty paramater used in constructing q^*

  // cubature
  int cubNp, cubNfp, cubNq;
  dfloat *cubr, *cubs, *cubt, *cubw; // coordinates and weights of local cubature nodes
  dfloat *cubx, *cuby, *cubz;    // coordinates of physical nodes
  dfloat *cubInterp; // interpolate from W&B to cubature nodes
  dfloat *cubProject; // projection matrix from cubature nodes to W&B nodes
  dfloat *cubD;       // 1D differentiation matrix
  dfloat *cubDiffInterp;     // 1D weak differentiation matrix
  dfloat *cubDW;     // 1D weak differentiation matrix
  dfloat *cubDrW;    // 'r' weak differentiation matrix
  dfloat *cubDsW;    // 's' weak differentiation matrix
  dfloat *cubDtW;    // 't' weak differentiation matrix
  dfloat *cubDWmatrices;

  dfloat *cubvgeo;  //volume geometric data at cubature points
  dfloat *cubsgeo;  //surface geometric data at cubature points
  dfloat *cubggeo;  //second type volume geometric data at cubature points
  
  // c2 at cubature points (for wadg)
  dfloat *c2;

  //source injection
  dfloat *sourceq;
  dfloat sourceX0, sourceY0, sourceZ0, sourceT0, sourceC2, sourceFreq;
  int sourceNelements;
  dlong *MRABsourceNelements;
  dlong *sourceElements;

  // surface integration node info
  int    intNfp;    // number of integration nodes on each face
  dfloat *intInterp; // interp from surface node to integration nodes
  dfloat *intLIFT;   // lift from surface integration nodes to W&B volume nodes
  dfloat *intx, *inty, *intz; // coordinates of suface integration nodes

  // Bernstein-Bezier info
  dfloat *VB, *invVB; // Bernstein Vandermonde matrices
  dfloat *BBMM;
  dfloat *invVB1D, *invVB2D;
  int *D0ids, *D1ids, *D2ids, *D3ids; // Bernstein deriv matrix indices
  dfloat *Dvals; // Bernstein deriv matrix values
  int *D0Tids, *D1Tids, *D2Tids, *D3Tids; // Bernstein transpose deriv matrix indices
  dfloat *DTvals; // Bernstein transpose deriv matrix values
  dfloat *VBq, *PBq; // cubature interpolation/projection matrices
  int *L0ids; // L0 matrix ids
  dfloat *L0vals; // L0 values (L0 tridiagonal in 2D)
  int *ELids; // lift reduction matrix indices
  dfloat *ELvals; // lift reduction matrix values
  int max_EL_nnz; // max number of non-zeros per row of EL
  int *BBRaiseids; //Bernstein elevate matrix indices
  dfloat *BBRaiseVals; //Bernstein elevate matrix values
  dfloat *BBLower; //Berstein projection matrix.

  //degree raising and lowering interpolation matrices
  dfloat *interpRaise;
  dfloat *interpLower;

  //sparse basis info
  dfloat *sparseV, *invSparseV;
  dfloat *sparseMM;
  int* FaceModes;
  int SparseNnzPerRow;
  int SparseNnzPerRowNonPadded;
  int *sparseStackedNZ;
  dfloat *sparseSrrT;
  dfloat *sparseSrsT;
  dfloat *sparseSssT;
  int *Ind;

  dlong *mmapM, *mmapP; 
  int   *mmapS;
  dfloat *mapSgn;

  // time stepping info
  dfloat dt; // time step
  dfloat startTime ; // Start Time
  dfloat finalTime; // final time to run acoustics to
  int   NtimeSteps;// number of time steps
  int   errorStep; // number of steps between error calculations
  int   Nrk;
  dfloat rka[5], rkb[5], rkc[6]; // AK: deprecated

  // MRAB,SAAB coefficients
  dfloat mrab[3], mrabb[3], saab[3], saabexp; // AK: deprecated 
  int MRABNlevels;
  int *MRABlevel;
  dlong *MRABNelements, *MRABNhaloElements;
  dlong **MRABelementIds, **MRABhaloIds;
  int *MRABshiftIndex;

  dlong *MRABpmlNelements, *MRABpmlNhaloElements;
  dlong **MRABpmlElementIds, **MRABpmlIds;
  dlong **MRABpmlHaloElementIds, **MRABpmlHaloIds;

  dlong pmlNelements, nonPmlNelements;
  dlong *nonPmlElementIds, *pmlElementIds, *pmlIds;  
  int shiftIndex;

  dfloat dtfactor; //Delete later for script run 
  dfloat maxErrorBoltzmann;

  dfloat *errtmp;
  dfloat rkC[7], rkA[7*7], rkE[7];

  occa::memory o_rkq, o_rkrhsq, o_rkerr; // deprecated, AK.
  occa::memory o_errtmp;
  occa::memory o_rkA, o_rkE;

  // ploting info for generating field vtu
  int    plotNverts;    // number of vertices for each plot element
  int    plotNp;        // number of plot nodes per element
  int    plotNelements; // number of "plot elements" per element
  int   *plotEToV;      // triangulation of plot nodes
  dfloat *plotR, *plotS, *plotT; // coordinates of plot nodes in reference element
  dfloat *plotInterp;    // warp & blend to plot node interpolation matrix

  int *contourEToV;
  dfloat *contourVX, *contourVY, *contourVZ;
  dfloat *contourInterp, *contourInterp1, *contourFilter; 

  //SEMFEM data
  int NpFEM, NelFEM;
  int *FEMEToV;
  dfloat *rFEM, *sFEM, *tFEM;
  dfloat *SEMFEMInterp;

  occa::memory o_SEMFEMInterp;
  occa::memory o_SEMFEMAnterp;

  // Boltzmann specific stuff
  dfloat RT, sqrtRT, tauInv, Ma, Re; // Deprecated: AK

  // pml stuff
  int    pmlNfields;
  //  dlong    pmlNelements; // deprecated
  dlong   *pmlElementList; // deprecated

  int Ntscale; // Will be removed, for time accuracy test
  
  dfloat *invTau; // deprecated in Boltzmann


  // Probe Data
  int probeN, probeNTotal; 
  dfloat *probeR, *probeS, *probeT;
  // dfloat *probeX, *probeY, *probeZ;  
  dlong *probeElementIds, *probeIds;  
  dfloat *probeI; 

  // occa stuff
  occa::device device;

  occa::stream defaultStream;
  occa::stream dataStream;
  occa::stream computeStream;

  occa::memory o_q, o_rhsq, o_resq, o_fQM, o_fQP;

  occa::memory o_Dr, o_Ds, o_Dt, o_LIFT, o_MM, o_invMM, o_MMPfloat;
  occa::memory o_DrT, o_DsT, o_DtT, o_LIFTT;
  occa::memory o_LMM, o_invLMM;
  occa::memory o_Dmatrices;
  occa::memory o_DmatricesPfloat;
  occa::memory o_FMMT;
  occa::memory o_sMT;

  occa::memory o_D; // tensor product differentiation matrix (for Hexes)
  occa::memory o_DW; // tensor product differentiation matrix (for Hexes)
  occa::memory o_SrrT, o_SrsT, o_SrtT; //element stiffness matrices
  occa::memory o_SsrT, o_SssT, o_SstT;
  occa::memory o_Srr, o_Srs, o_Srt, o_Sss, o_Sst, o_Stt; // for char4-based kernels
  occa::memory o_Smatrices;
  occa::memory o_SmatricesPfloat;
  occa::memory o_IndT, o_IndTchar;
  occa::memory o_India, o_Indja;
  occa::memory o_StrT, o_StsT, o_SttT;
  occa::memory o_Ind; // for sparse index storage

  occa::memory o_vgeo, o_sgeo;
  occa::memory o_vmapM, o_vmapP, o_mapP;

  occa::memory o_rmapP;

  occa::memory o_EToE, o_EToF, o_EToB, o_x, o_y, o_z;

  occa::memory o_EToFPairs, o_FPairsToE, o_FPairsToF;

  // cubature (for wadg)
  occa::memory o_intLIFTT, o_intInterpT, o_intx, o_inty, o_intz;
  occa::memory o_cubDWT, o_cubD;
  occa::memory o_cubDrWT, o_cubDsWT, o_cubDtWT, o_cubDiffInterpT;
  occa::memory o_cubDWmatrices;
  occa::memory o_cubInterpT, o_cubProjectT;
  occa::memory o_invMc; // for comparison: inverses of weighted mass matrices

  occa::memory o_cubvgeo, o_cubsgeo, o_cubggeo;

  occa::memory o_c2;

  //MRAB element lists
  occa::memory *o_MRABelementIds;
  occa::memory *o_MRABhaloIds;
  occa::memory *o_MRABpmlElementIds;
  occa::memory *o_MRABpmlIds;
  occa::memory *o_MRABpmlHaloElementIds;
  occa::memory *o_MRABpmlHaloIds;


  // DG halo exchange info
  occa::memory o_haloElementList;
  occa::memory o_haloBuffer;
  occa::memory o_haloGetNodeIds;
  occa::memory o_haloPutNodeIds;
  
  occa::memory o_internalElementIds;
  occa::memory o_notInternalElementIds;

  // Bernstein-Bezier occa arrays
  occa::memory o_BBMM;
  occa::memory o_D0ids, o_D1ids, o_D2ids, o_D3ids, o_Dvals; // Bernstein deriv matrix indices
  occa::memory o_packedDids; // char4 packed increments (D1ids-D0ids)

  occa::memory o_invVB1DT, o_invVB2DT;
  occa::memory o_VBq, o_PBq; // cubature interpolation/projection matrices
  occa::memory o_L0ids, o_L0vals, o_ELids, o_ELvals;

  /* sparse basis occa arrays */
  occa::memory o_sparseStackedNZ;
  occa::memory o_sparseSrrT;
  occa::memory o_sparseSrsT;
  occa::memory o_sparseSssT;
  occa::memory o_mapSgn;

  // pml vars
  occa::memory o_sigmax, o_sigmay, o_sigmaz; // AK: deprecated


  occa::memory o_pmlElementIds;
  occa::memory o_nonPmlElementIds;
  occa::memory o_pmlIds;

  occa::memory o_pmlElementList;
  
  occa::memory o_ggeo; // second order geometric factors
  occa::memory o_ggeoPfloat; // second order geometric factors
  occa::memory o_projectL2; // local weights for projection.

  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel updateKernel;
  occa::kernel traceUpdateKernel;
  occa::kernel haloExtractKernel;
  occa::kernel partialSurfaceKernel;
  occa::kernel haloGetKernel;
  occa::kernel haloPutKernel;

  // Just for test will be deleted after temporal testsAK
  occa::kernel RKupdateKernel;
  occa::kernel RKpmlUpdateKernel;


  occa::kernel gatherKernel;
  occa::kernel scatterKernel;
  occa::kernel gatherScatterKernel;

  occa::kernel getKernel;
  occa::kernel putKernel;

  occa::kernel sumKernel;
  occa::kernel addScalarKernel;

  occa::kernel AxKernel;
  occa::kernel innerProductKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotDivideKernel;

  occa::kernel gradientKernel;
  occa::kernel ipdgKernel;

  occa::kernel maskKernel;
  occa::kernel maskPfloatKernel;

  // Boltzmann Specific Kernels
  occa::kernel relaxationKernel;
  occa::kernel pmlRelaxationKernel;
  
}mesh_t;

// serial sort
void mysort(hlong *data, int N, const char *order);

// sort entries in an array in parallel
void parallelSort(int size, int rank, MPI_Comm comm,
		  int N, void *vv, size_t sz,
		  int (*compare)(const void *, const void *),
		  void (*match)(void *, void *)
		  );

#define mymax(a,b) (((a)>(b))?(a):(b))
#define mymin(a,b) (((a)<(b))?(a):(b))

/* hash function */
unsigned int hash(const unsigned int value) ;

/* dimension independent mesh operations */
void meshConnect(mesh_t *mesh);

/* build parallel face connectivity */
void meshParallelConnect(mesh_t *mesh);

/* build global connectivity in parallel */
void meshParallelConnectNodes(mesh_t *mesh, int isTmesh, int nrsBuildOnly);

void meshHaloSetup(mesh_t *mesh);

/* extract whole elements for the halo exchange */
void meshHaloExtract(mesh_t *mesh, size_t Nbytes, void *sourceBuffer, void *haloBuffer);

void meshHaloExchange(mesh_t *mesh,
    size_t Nbytes,         // message size per element
    void *sourceBuffer,
    void *sendBuffer,    // temporary buffer
    void *recvBuffer);

void meshHaloExchangeStart(mesh_t *mesh,
    size_t Nbytes,       // message size per element
    void *sendBuffer,    // temporary buffer
    void *recvBuffer);


void meshHaloExchangeFinish(mesh_t *mesh);

void meshHaloExchangeBlocking(mesh_t *mesh,
			     size_t Nbytes,       // message size per element
			     void *sendBuffer,    // temporary buffer
			      void *recvBuffer);

// print out parallel partition i
void meshPartitionStatistics(mesh_t *mesh);

// build element-boundary connectivity
void meshConnectBoundary(mesh_t *mesh);

void meshParallelGatherScatterSetup(mesh_t *mesh,
                                      dlong N,
                                      hlong *globalIds,
                                      MPI_Comm &comm,
                                      int verbose);

// generic mesh setup
mesh_t *meshSetup(char *filename, int N, setupAide &options);
void meshFree(mesh_t*);

void occaTimerTic(occa::device device,std::string name);
void occaTimerToc(occa::device device,std::string name);

extern "C"
{
  void * xxtSetup(uint num_local_rows,
      void* row_ids,
      uint nnz,
      void*   A_i,
      void*   A_j,
      void* A_vals,
      int null_space,
      const char* inttype,
      const char* floattype);

  void xxtSolve(void* x,
      void* A,
      void* rhs);

  void xxtFree(void* A) ;
}

extern "C"
{
  void dgesv_ ( int     *N, int     *NRHS, double  *A,
                int     *LDA,
                int     *IPIV, 
                double  *B,
                int     *LDB,
                int     *INFO );

  // void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, 
  //             const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B, 
  //             const int *LDB, double *BETA, double *C, const int *LDC);

   void dgemm_ (char *, char *, int *, int *, int *,
         const dfloat *, const dfloat * __restrict, int *,
         const dfloat * __restrict, int *,
         const dfloat *, dfloat * __restrict, int *);

  void sgesv_(int *N, int *NRHS,float  *A, int *LDA, int *IPIV, float  *B, int *LDB,int *INFO);

  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
              double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
  
  double dlange_(char *NORM, int *M, int *N, double *A, int *LDA, double *WORK);
  void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM,
                double *RCOND, double *WORK, int *IWORK, int *INFO );
}

void readDfloatArray(FILE *fp, const char *label, dfloat **A, int *Nrows, int* Ncols);
void readIntArray   (FILE *fp, const char *label, int **A   , int *Nrows, int* Ncols);

void meshApplyElementMatrix(mesh_t *mesh, dfloat *A, dfloat *q, dfloat *Aq);
void meshApplyVectorElementMatrix(mesh_t *mesh, int Nfield, const dlong offset, dfloat *A, dfloat *q, dfloat *Aq);

void meshRecursiveSpectralBisectionPartition(mesh_t *mesh);

void matrixInverse(int N, dfloat *A);
dfloat matrixConditionNumber(int N, dfloat *A);

#if 0
void *occaHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem);
#else
void *occaHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem, occa::memory &h_mem);
#endif

void matrixRightSolve(int NrowsA, int NcolsA, dfloat *A, int NrowsB, int NcolsB, dfloat *B, dfloat *C);
void matrixEig(int N, dfloat *A, dfloat *VR, dfloat *WR, dfloat *WI);
void matrixTranspose(const int M, const int N,
                     const dfloat  *A, const int LDA,
                           dfloat *AT, const int LDAT);

// 1D mesh basis functions
void Nodes1D(int _N, dfloat *_r);
void EquispacedNodes1D(int _N, dfloat *_r);
void OrthonormalBasis1D(dfloat a, int i, dfloat *P);
void GradOrthonormalBasis1D(dfloat a, int i, dfloat *Pr);
void Vandermonde1D(int _N, int Npoints, dfloat *_r, dfloat *V);
void GradVandermonde1D(int _N, int Npoints, dfloat *_r, dfloat *Vr);
void MassMatrix1D(int _Np, dfloat *V, dfloat *_MM);
void Dmatrix1D(int _N, int NpointsIn, dfloat *_rIn,
                               int NpointsOut, dfloat *_rOut, dfloat *_Dr);
void DWmatrix1D(int _N, dfloat *_D, dfloat *_DT);

void InterpolationMatrix1D(int _N,
                               int NpointsIn, dfloat *rIn,
                               int NpointsOut, dfloat *rOut,
                               dfloat *I);
void DegreeRaiseMatrix1D(int Nc, int Nf, dfloat *P);
void CubatureWeakDmatrix1D(int _Nq, int _cubNq,
                                     dfloat *_cubProject, dfloat *_cubD, dfloat *_cubPDT);
dfloat JacobiP(dfloat a, dfloat alpha, dfloat beta, int _N);
dfloat GradJacobiP(dfloat a, dfloat alpha, dfloat beta, int _N);
void JacobiGLL(int _N, dfloat *_x, dfloat *_w = nullptr);
void JacobiGQ(dfloat alpha, dfloat beta, int _N, dfloat *_x, dfloat *_w);
} // end C Linkage
#endif

