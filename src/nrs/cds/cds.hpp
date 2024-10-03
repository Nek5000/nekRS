#ifndef CDS_H
#define CDS_H

#include "platform.hpp"
#include "solver.hpp"
#include "elliptic.hpp"
#include "cvode.hpp"

struct cdsConfig_t {

  int Nscalar;
  mesh_t *meshT;
  mesh_t *meshV;
  dfloat *g0 = nullptr;
  dfloat *dt = nullptr;
  int nBDF;
  dfloat *coeffBDF = nullptr;
  occa::memory o_coeffBDF;
  int nEXT;
  dfloat *coeffEXT = nullptr;
  occa::memory o_coeffEXT;
  deviceMemory<dfloat> *o_usrwrk;
  dlong vFieldOffset;
  dlong vCubatureOffset;
  dlong fieldOffset;
  int Nsubsteps;
  occa::memory o_U;
  occa::memory o_Ue;
  occa::memory o_Urst;
  occa::memory o_relUrst;

  bool dpdt = false;
  dfloat *dp0thdt = nullptr;
  dfloat *alpha0Ref = nullptr;

};

class cds_t : public solver_t
{

  using userSource_t = std::function<void(double)>;
  using userProperties_t = std::function<void(double)>;
  using userImplicitLinearTerm_t = std::function<occa::memory(double, int)>;

public:
  static constexpr double targetTimeBenchmark{0.2};

  cds_t(cdsConfig_t &cfg);
  void solve(double time, int stage);
  void makeNLT(int is, double time, int tstep, occa::memory &o_Subcycling);
  occa::memory advectionSubcyling(int nEXT, double time, int scalarIdx);

  void saveSolutionState();
  void restoreSolutionState();

  void applyAVM();

  userSource_t userSource = nullptr;
  userProperties_t userProperties = nullptr;
  userImplicitLinearTerm_t userImplicitLinearTerm = nullptr;

  std::vector<mesh_t *> mesh;
  std::vector<dlong> fieldOffset;
  std::vector<dlong> fieldOffsetScan;
  occa::memory o_fieldOffsetScan;
  dlong fieldOffsetSum;
  mesh_t *meshV;
  std::vector<elliptic *> solver;
  cvode_t *cvode = nullptr;

  bool anyCvodeSolver = false;
  bool anyEllipticSolver = false;

  bool cht = false;

  int NVfields = 3;
  int NSfields = 0;

  oogs_t *gsh = nullptr;
  oogs_t *gshT = nullptr;
  QQt *qqt = nullptr;
  QQt *qqtT = nullptr;

  dlong vFieldOffset;
  dlong vCubatureOffset;
  dfloat *dt = nullptr;
  dfloat *g0 = nullptr;

  int nEXT;
  int nBDF;

  std::vector<int> compute;
  std::vector<int> cvodeSolve;
  occa::memory o_compute;
  occa::memory o_cvodeSolve;

  std::vector<dfloat> filterS;
  occa::memory o_applyFilterRT;
  occa::memory o_filterS;
  occa::memory o_filterRT;
  int applyFilter = 0;

  int Nsubsteps = 0;

  bool dpdt = false;
  dfloat *dp0thdt = nullptr;
  dfloat *alpha0Ref = nullptr;

  int *EToB = nullptr;
  occa::memory o_EToB;
  dlong EToBOffset;

  deviceMemory<dfloat> *o_usrwrk;

  occa::memory o_U;
  occa::memory o_Ue;
  occa::memory o_relUrst;
  occa::memory o_Urst;

  occa::memory o_S, o_Se;
  occa::memory o_prop;
  occa::memory o_rho, o_diff;
  occa::memory o_NLT, o_JwF;

  dfloat *coeffEXT, *coeffBDF;
  occa::memory o_coeffEXT, o_coeffBDF;

  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;
  occa::kernel filterRTKernel;
  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel advectionCubatureVolumeKernel;
  occa::kernel advectionCubatureSurfaceKernel;
  occa::kernel strongAdvectionVolumeKernel;
  occa::kernel strongAdvectionCubatureVolumeKernel;
  occa::kernel advectMeshVelocityKernel;
  occa::kernel neumannBCKernel;
  occa::kernel dirichletBCKernel;
  occa::kernel maskCopyKernel;
  occa::kernel maskCopy2Kernel;

  occa::properties *kernelInfo;

private:
  occa::memory o_Ssave;
  occa::memory o_NLTsave;
  occa::memory o_Spropsave;
};

#endif
