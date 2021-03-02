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
#include "platform.hpp"

struct nrs_t
{
  int dim, elementType;

  mesh_t* meshV;
  mesh_t* meshT;
  linAlg_t* linAlg;

  elliptic_t* uSolver;
  elliptic_t* vSolver;
  elliptic_t* wSolver;
  elliptic_t* uvwSolver;
  elliptic_t* pSolver;

  cds_t* cds;

  oogs_t* gsh;

  dlong ellipticWrkOffset;

  int flow;

  int Nscalar;
  setupAide options;
  setupAide vOptions, pOptions;

  inipp::Ini<char> *par;

  int NVfields, NTfields;
  dlong fieldOffset;
  dlong Ntotal;


  dfloat dt[3], idt;
  dfloat p0th[3] = {0.0, 0.0, 0.0};
  dfloat dp0thdt;
  int tstep;
  int lastStep;
  dfloat g0, ig0;

  int cht;

  int nEXT;
  int nBDF;
  int isOutputStep;
  int outputForceStep;

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
  dfloat* coeffEXT, * coeffBDF, * coeffSubEXT;
  dfloat* extC;

  int* VmapB;
  occa::memory o_VmapB;

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
  occa::kernel advectMeshVelocityKernel;

  occa::kernel pressureAddQtlKernel;
  occa::kernel pressureStressKernel;

  occa::kernel PQKernel;
  occa::kernel mueDivKernel;

  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel extrapolateKernel;
  occa::kernel subCycleRKKernel;
  occa::kernel subCycleExtrapolateFieldKernel;
  occa::kernel subCycleExtrapolateScalarKernel;

  occa::kernel wgradientVolumeKernel;

  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;

  occa::kernel constrainKernel;

  occa::memory o_U, o_P;

  occa::memory o_BF;
  occa::memory o_FU;

  dfloat* wrk;

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
  occa::memory o_coeffEXT, o_coeffBDF, o_coeffSubEXT;
  occa::memory o_extC;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionCubatureVolumeKernel;

  occa::kernel advectionStrongVolumeKernel;
  occa::kernel advectionStrongCubatureVolumeKernel;

  occa::kernel diffusionKernel;
  occa::kernel velocityGradientKernel;

  occa::kernel gradientVolumeKernel;

  occa::kernel wDivergenceVolumeKernel;
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

  occa::kernel cflKernel;

  occa::kernel setEllipticCoeffKernel;
  occa::kernel setEllipticCoeffPressureKernel;

  occa::kernel pressureAxKernel;
  occa::kernel curlKernel;
  occa::kernel maskCopyKernel;

  int* EToB;
  occa::memory o_EToB;

  occa::properties* kernelInfo;
};


#include "io.hpp"

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
