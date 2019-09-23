#if !defined(nekrs_ins_hpp_)
#define nekrs_ins_hpp_

#include "mesh2D.h"
#include "mesh3D.h"
#include "elliptic.h"
#include "cds.h"

extern "C" { // Begin C Linkage
typedef struct {

  int dim, elementType;

  mesh_t *mesh;
  elliptic_t *uSolver;
  elliptic_t *vSolver;
  elliptic_t *wSolver;
  elliptic_t *pSolver;
  cds_t      *cds;
  
  int Nscalar; 
  setupAide options;
  setupAide vOptions, pOptions; 	

  // INS SOLVER OCCA VARIABLES
  dfloat nu, Re;
  int NVfields, NTfields;
  dlong fieldOffset;
  dlong Ntotal;

  int Nblock;

  dfloat dt, dti;          // time step
  dfloat time;
  int tstep;
  dfloat g0, ig0, lambda;      // helmhotz solver -lap(u) + lamda u
  dfloat startTime;   
  dfloat finalTime;   

  int temporalOrder;
  int ExplicitOrder; 
  int   NtimeSteps;  // number of time steps 
  int   Nstages;     
  int   outputStep;
  int   isOutputStep;
  int   outputForceStep; 

  int ARKswitch;
  
  int NiterU, NiterV, NiterW, NiterP;

  //solver tolerances
  dfloat presTOL, velTOL;

  dfloat idt, inu; // hold some inverses
  
  dfloat *U, *P;
  dfloat *NU, *LU, *GP;
  dfloat *GU;   
  dfloat *rhsU, *rhsV, *rhsW, *rhsP;   
  dfloat *rkU, *rkP, *PI;
  dfloat *rkNU, *rkLU, *rkGP;

  dfloat *FU; // Additional source terms for explicit contribution
  
  dfloat *Vort, *Div;

  //RK Subcycle Data
  int SNrk;
  dfloat *Srka, *Srkb, *Srkc; 

  //ARK data
  int Nrk;
  dfloat *rkC;
  dfloat *erkA, *irkA, *prkA;
  dfloat *erkB, *irkB, *prkB;
  dfloat *erkE, *irkE, *prkE;
  int embeddedRKFlag;

  //EXTBDF data
  dfloat *extbdfA, *extbdfB, *extbdfC;
  dfloat *extC;

  int *VmapB, *PmapB;
  occa::memory o_VmapB, o_PmapB;

  //halo data
  dfloat *vSendBuffer;
  dfloat *vRecvBuffer;
  dfloat *pSendBuffer;
  dfloat *pRecvBuffer;
  dfloat * velocityHaloGatherTmp;

  occa::memory o_vSendBuffer,h_vSendBuffer;
  occa::memory o_vRecvBuffer,h_vRecvBuffer;
  occa::memory o_pSendBuffer,h_pSendBuffer;
  occa::memory o_pRecvBuffer,h_pRecvBuffer;
  occa::memory o_gatherTmpPinned, h_gatherTmpPinned;

  int Nsubsteps;  
  dfloat *Ud, *Ue, *resU, *rhsUd, sdt;
  occa::memory o_Ud, o_Ue, o_resU, o_rhsUd;

  dfloat *cU, *cUd;
  occa::memory o_cU, o_cUd;

  dfloat *Wrk;
  occa::memory o_Wrk; 

  // Cfl related
  int computedDh; 
  occa::memory o_idH; // i.e. inverse of 1D Gll Spacing for quad and Hex

  int readRestartFile,writeRestartFile, restartedFromFile;

  // Filter Stabilization Matrix
  int filterNc; // filter cut modes i.e. below is not touched
  dfloat *filterM, filterS; 
  occa::memory o_filterMT; // transpose of filter matrix 
  occa::kernel VFilterKernel; // Relaxation-Term based filtering
  occa::kernel SFilterKernel; // Relaxation-Term based filtering

  occa::kernel scaledAddKernel;
  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel ;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel subCycleExtKernel;

  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;

  occa::memory o_invLumpedMassMatrix;
  
  occa::kernel constrainKernel;
  
  occa::memory o_U, o_P;
  occa::memory o_rhsU, o_rhsV, o_rhsW, o_rhsP; 

  occa::memory o_NU, o_LU, o_GP, o_NC;
  occa::memory o_GU;

  occa::memory o_FU; 

  occa::memory o_UH, o_VH, o_WH;
  occa::memory o_rkU, o_rkP, o_PI;
  occa::memory o_rkNU, o_rkLU, o_rkGP;

  occa::memory o_Vort, o_Div;

  occa::memory o_vHaloBuffer, o_pHaloBuffer; 
  occa::memory o_velocityHaloGatherTmp;

  occa::kernel haloGetKernel;
  occa::kernel haloPutKernel;
  
  //ARK data
  occa::memory o_rkC;
  occa::memory o_erkA, o_irkA, o_prkA;
  occa::memory o_erkB, o_irkB, o_prkB;
  occa::memory o_erkE, o_irkE, o_prkE;

  //EXTBDF data
  occa::memory o_extbdfA, o_extbdfB, o_extbdfC;
  occa::memory o_extC;

  occa::kernel velocityHaloExtractKernel;
  occa::kernel velocityHaloScatterKernel;
  occa::kernel pressureHaloExtractKernel;
  occa::kernel pressureHaloScatterKernel;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel advectionCubatureVolumeKernel;
  occa::kernel advectionCubatureSurfaceKernel;

  occa::kernel advectionStrongVolumeKernel;
  occa::kernel advectionStrongCubatureVolumeKernel;
  
  occa::kernel diffusionKernel;
  occa::kernel velocityGradientKernel;

  occa::kernel gradientVolumeKernel;
  occa::kernel gradientSurfaceKernel;

  occa::kernel divergenceVolumeKernel;
  occa::kernel divergenceSurfaceKernel;

  occa::kernel divergenceStrongVolumeKernel;
  
  occa::kernel pressureRhsKernel;
  occa::kernel pressureRhsBCKernel;
  occa::kernel pressureAddBCKernel;
  occa::kernel pressurePenaltyKernel;
  occa::kernel pressureUpdateKernel;

  occa::kernel velocityRhsKernel;
  occa::kernel velocityRhsBCKernel;
  occa::kernel velocityAddBCKernel;
  occa::kernel velocityUpdateKernel;  
  
  occa::kernel setScalarKernel; 

  occa::kernel cflKernel; 
  occa::kernel maxKernel; 

  int TOMBO;  
  occa::kernel pressureAxKernel; 
  occa::kernel velocityAxKernel; 
  occa::kernel curlKernel; 
  occa::kernel curlBKernel; // needed for 2D
  occa::kernel invMassMatrixKernel; 
  occa::kernel massMatrixKernel; 

  int *EToB;
  occa::memory o_EToB;

  occa::memory o_InvM;
  occa::properties *kernelInfo;

}ins_t;

} // end C Linkage

#endif
