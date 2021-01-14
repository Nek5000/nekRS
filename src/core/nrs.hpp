#if !defined(nekrs_nekrs_hpp_)
#define nekrs_nekrs_hpp_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <getopt.h>

#include "nrssys.hpp"
#include "mesh3D.h"
#include "elliptic.h"
#include "cds.hpp"
#include "linAlg.hpp"
#include "timer.hpp"
#include "inipp.hpp"

typedef struct
{
  int dim, elementType;

  mesh_t* mesh;
  mesh_t* meshT;

  elliptic_t* uSolver;
  elliptic_t* vSolver;
  elliptic_t* wSolver;
  elliptic_t* uvwSolver;
  elliptic_t* pSolver;

  cds_t* cds;

  oogs_t* gsh;

  linAlg_t* linAlg;

  dlong ellipticWrkOffset;

  int flow;

  int Nscalar;
  setupAide options;
  setupAide vOptions, pOptions;

  inipp::Ini<char> *par;

  int NVfields, NTfields;
  dlong fieldOffset;
  dlong Nlocal, Ntotal;

  int Nblock;

  dfloat dt[3], idt;
  int tstep;
  int lastStep;
  dfloat g0, ig0;

  int cht;

  int temporalOrder;
  int ExplicitOrder;
  int Nstages;
  int isOutputStep;
  int outputForceStep;

  int NiterU, NiterV, NiterW, NiterP;

  dfloat* U, * P;
  dfloat* BF, * FU;

  //RK Subcycle Data
  int SNrk;
  dfloat* Srka, * Srkb, * Srkc;
  occa::memory o_Srka, o_Srkb;

  //ARK data
  int Nrk;
  dfloat* rkC;

  //EXTBDF data
  dfloat* extbdfA, * extbdfB, * extbdfC;
  dfloat* extC;

  int* VmapB;
  occa::memory o_VmapB;

  occa::memory o_wrk0, o_wrk1, o_wrk2, o_wrk3, o_wrk4, o_wrk5, o_wrk6, o_wrk7,
               o_wrk9, o_wrk12, o_wrk15;

  int Nsubsteps;
  dfloat* Ue, sdt;
  occa::memory o_Ue;

  dfloat* div;
  occa::memory o_div;

  dfloat rho, mue;
  occa::memory o_rho, o_mue;

  dfloat* usrwrk;
  occa::memory o_usrwrk;

  occa::memory o_idH; // i.e. inverse of 1D Gll Spacing for quad and Hex

  int filterNc; // filter cut modes i.e. below is not touched
  dfloat* filterM, filterS;
  occa::memory o_filterMT; // transpose of filter matrix
  occa::kernel filterRTKernel; // Relaxation-Term based filtering

  occa::kernel qtlKernel;
  occa::kernel pressureAddQtlKernel;
  occa::kernel pressureStressKernel;

  occa::kernel PQKernel;
  occa::kernel mueDivKernel;
  occa::kernel dotMultiplyKernel;

  occa::kernel scalarScaledAddKernel;
  occa::kernel scaledAddKernel;
  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel extrapolateKernel;

  occa::kernel wgradientVolumeKernel;

  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;

  occa::kernel constrainKernel;

  occa::memory o_U, o_P;

  occa::memory o_BF;
  occa::memory o_FU;

  int var_coeff;
  dfloat* prop, * ellipticCoeff;
  occa::memory o_prop, o_ellipticCoeff;

  occa::memory o_UH;

  occa::memory o_vHaloBuffer, o_pHaloBuffer;
  occa::memory o_velocityHaloGatherTmp;

  occa::kernel haloGetKernel;
  occa::kernel haloPutKernel;

  //ARK data
  occa::memory o_rkC;

  //EXTBDF data
  occa::memory o_extbdfA, o_extbdfB, o_extbdfC;
  occa::memory o_extC;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionCubatureVolumeKernel;

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
  occa::kernel pressureDirichletBCKernel;
  occa::kernel pressurePenaltyKernel;
  occa::kernel pressureUpdateKernel;

  occa::kernel velocityRhsKernel;
  occa::kernel velocityNeumannBCKernel;
  occa::kernel velocityDirichletBCKernel;

  occa::kernel fillKernel;

  occa::kernel cflKernel;
  occa::kernel maxKernel;

  occa::kernel setEllipticCoeffKernel;
  occa::kernel setEllipticCoeffPressureKernel;

  occa::kernel pressureAxKernel;
  occa::kernel curlKernel;
  occa::kernel invMassMatrixKernel;
  occa::kernel massMatrixKernel;

  occa::kernel maskCopyKernel;

  int* EToB;
  occa::memory o_EToB;

  occa::properties* kernelInfo;
} nrs_t;


#include "io.hpp"

occa::device occaDeviceConfig(setupAide &options, MPI_Comm comm);

// std::to_string might be not accurate enough
static string to_string_f(double a)
{
  stringstream s;
  s << std::scientific << a;
  return s.str();
}

static std::vector<std::string> serializeString(const std::string sin)
{
  std::vector<std::string> slist;
  string s(sin);
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  std::stringstream ss;
  ss.str(s);
  while( ss.good() ) {
    std::string substr;
    std::getline(ss, substr, ',');
    slist.push_back(substr);
  }
  return slist;
}

#endif
