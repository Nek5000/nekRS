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
  mesh_t *meshT;

  elliptic_t *uSolver;
  elliptic_t *vSolver;
  elliptic_t *wSolver;
  elliptic_t *uvwSolver;
  elliptic_t *pSolver;
  cds_t      *cds;

  dlong ellipticWrkOffset; 
  
  int Nscalar; 
  setupAide options;
  setupAide vOptions, pOptions; 	

  // INS SOLVER OCCA VARIABLES
  int NVfields, NTfields;
  dlong fieldOffset;
  dlong Nlocal, Ntotal;

  int Nblock;

  dfloat dt, idt;
  dfloat time;
  int tstep;
  dfloat g0, ig0;
  dfloat startTime;   
  dfloat finalTime;   

  int   cht;

  int   temporalOrder;
  int   ExplicitOrder; 
  int   NtimeSteps;  // number of time steps 
  int   Nstages;     
  int   outputStep;
  int   isOutputStep;
  int   outputForceStep; 

  int NiterU, NiterV, NiterW, NiterP;
  dfloat presTOL, velTOL;

  dfloat *U, *P;
  dfloat *BF, *FU;
  dfloat *PI;

  //RK Subcycle Data
  int SNrk;
  dfloat *Srka, *Srkb, *Srkc;
  occa::memory o_Srka, o_Srkb; 

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

  int *VmapB;
  occa::memory o_VmapB;

  dlong *elementInfo;
  occa::memory o_elementInfo;

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
  
  occa::memory o_wrk0, o_wrk1, o_wrk2, o_wrk3, o_wrk4, o_wrk5, o_wrk6,
               o_wrk9, o_wrk12, o_wrk15; 

  int Nsubsteps;  
  dfloat *Ue, sdt;
  occa::memory o_Ue;

  int lowMach;
  dfloat *div;
  occa::memory o_div;

  dfloat rho, mue;
  occa::memory o_rho, o_mue;

  dfloat *usrwrk;
  occa::memory o_usrwrk; 

  occa::memory o_idH; // i.e. inverse of 1D Gll Spacing for quad and Hex

  int readRestartFile,writeRestartFile, restartedFromFile;

  int filterNc; // filter cut modes i.e. below is not touched
  dfloat *filterM, filterS; 
  occa::memory o_filterMT; // transpose of filter matrix 
  occa::kernel filterRTKernel; // Relaxation-Term based filtering

  occa::kernel qtlKernel;
  occa::kernel pressureAddQtlKernel;

  occa::kernel pqKernel;
  occa::kernel ncKernel;

  occa::kernel scalarScaledAddKernel;
  occa::kernel scaledAddKernel;
  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel ;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel velocityExtKernel;

  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;

  occa::kernel constrainKernel;
  
  occa::memory o_U, o_P;

  occa::memory o_BF;
  occa::memory o_FU; 

  int var_coeff;
  dfloat *prop, *ellipticCoeff;
  occa::memory o_prop, o_ellipticCoeff;

  occa::memory o_UH;
  occa::memory o_PI;

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

  occa::kernel divergenceVolumeKernel;
  occa::kernel divergenceSurfaceKernel;

  occa::kernel divergenceStrongVolumeKernel;
  occa::kernel sumMakefKernel; 
  occa::kernel pressureRhsKernel;
  occa::kernel pressureAddBCKernel;
  occa::kernel pressurePenaltyKernel;
  occa::kernel pressureUpdateKernel;

  occa::kernel velocityRhsKernel;
  occa::kernel velocityRhsBCKernel;
  occa::kernel velocityAddBCKernel;
  
  occa::kernel setScalarKernel; 

  occa::kernel cflKernel; 
  occa::kernel maxKernel; 

  occa::kernel setEllipticCoeffKernel;
  occa::kernel setEllipticCoeffPressureKernel;

  int TOMBO;  
  occa::kernel AxKernel; 
  occa::kernel curlKernel; 
  occa::kernel invMassMatrixKernel; 
  occa::kernel massMatrixKernel; 
  
  int *EToB;
  occa::memory o_EToB;

  occa::memory o_InvM;
  occa::properties *kernelInfo;

}ins_t;

} // end C Linkage

#endif
