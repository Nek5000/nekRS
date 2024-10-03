#include "advectionSubCycling.hpp"
#include "nrs.hpp"

#include "cds.hpp"
#include "lowPassFilter.hpp"
#include "avm.hpp"
#include "bcMap.hpp"

static void advectionFlops(mesh_t *mesh, int Nfields)
{
  const auto cubNq = mesh->cubNq;
  const auto cubNp = mesh->cubNp;
  const auto Nq = mesh->Nq;
  const auto Np = mesh->Np;
  const auto Nelements = mesh->Nelements;
  double flopCount = 0.0; // per elem basis
  if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
    flopCount += 4. * Nq * (cubNp + cubNq * cubNq * Nq + cubNq * Nq * Nq); // interpolation
    flopCount += 6. * cubNp * cubNq;                                       // apply Dcub
    flopCount += 5 * cubNp; // compute advection term on cubature mesh
    flopCount += mesh->Np;  // weight by inv. mass matrix
  } else {
    flopCount += 8 * (Np * Nq + Np);
  }

  flopCount *= Nelements;
  flopCount *= Nfields;

  platform->flopCounter->add("advection", flopCount);
}

occa::memory cds_t::advectionSubcyling(int nEXT, double time, int scalarIdx)
{
  const auto movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  const auto mesh = this->mesh[scalarIdx];
  const auto gshV = this->gsh;

  const auto nFields = 1;
  const auto fieldOffset = this->fieldOffset[scalarIdx];
  const auto fieldOffsetSum = this->fieldOffsetSum;
  const auto meshOffset =
      this->fieldOffset[0]; // used for mesh->o_invLMM and mesh->o_divU in case of moving mesh

  auto o_U = this->o_S.slice(this->fieldOffsetScan[scalarIdx], fieldOffset);

  auto kernel = (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
                      ? this->subCycleStrongCubatureVolumeKernel
                      : this->subCycleStrongVolumeKernel;

  return advectionSubcyclingRK(mesh,
                               this->meshV,
                               time,
                               this->dt,
                               this->Nsubsteps,
                               this->coeffBDF,
                               nEXT,
                               nFields,
                               kernel,
                               gshV,
                               meshOffset,
                               fieldOffset,
                               this->vCubatureOffset,
                               fieldOffsetSum,
                               (movingMesh) ? this->o_relUrst : this->o_Urst,
                               o_U);
}

cds_t::cds_t(cdsConfig_t &cfg)
{
  auto &options = platform->options;
  const std::string section = "cds-";
  platform_t *platform = platform_t::getInstance();

  this->NSfields = cfg.Nscalar;

  this->mesh.resize(this->NSfields);
  this->fieldOffset.resize(this->NSfields);
  this->fieldOffsetScan.resize(this->NSfields);
  this->solver.resize(this->NSfields);
  this->compute.resize(this->NSfields);
  this->cvodeSolve.resize(this->NSfields);
  this->filterS.resize(this->NSfields);

  this->NVfields = cfg.meshV->dim;
  this->g0 = cfg.g0;
  this->dt = cfg.dt;
  this->nBDF = cfg.nBDF;
  this->coeffBDF = cfg.coeffBDF;
  this->o_coeffBDF = cfg.o_coeffBDF;
  this->nEXT = cfg.nEXT;
  this->coeffEXT = cfg.coeffEXT;
  this->o_coeffEXT = cfg.o_coeffEXT;
  this->o_usrwrk = cfg.o_usrwrk;
  this->vFieldOffset = cfg.fieldOffset;
  this->Nsubsteps = cfg.Nsubsteps;
  this->vCubatureOffset = cfg.vCubatureOffset;
  this->fieldOffset[0] = cfg.fieldOffset;
  this->o_U = cfg.o_U;
  this->o_Ue = cfg.o_Ue;
  this->o_Urst = cfg.o_Urst;
  this->o_relUrst = cfg.o_relUrst;
  this->mesh[0] = cfg.meshT;
  this->meshV = cfg.meshV;

  this->dpdt = cfg.dpdt;
  this->dp0thdt = cfg.dp0thdt;
  this->alpha0Ref = cfg.alpha0Ref;

  this->cht = (cfg.meshT != cfg.meshV) ? true : false;

  this->fieldOffsetScan[0] = 0;
  dlong sum = this->fieldOffset[0];
  for (int s = 1; s < this->NSfields; ++s) {
    this->fieldOffset[s] = this->fieldOffset[0]; // same for all scalars
    this->fieldOffsetScan[s] = sum;
    sum += this->fieldOffset[s];
    this->mesh[s] = this->meshV;
  }
  this->fieldOffsetSum = sum;

  this->o_fieldOffsetScan = platform->device.malloc<dlong>(this->NSfields, this->fieldOffsetScan.data());

  this->gsh = oogs::setup(this->meshV->ogs, 1, 0, ogsDfloat, NULL, OOGS_AUTO);
  this->qqt = new QQt(this->gsh);

  this->gshT = (this->cht)
                   ? oogs::setup(this->mesh[0]->ogs, 1, 0, ogsDfloat, NULL, OOGS_AUTO)
                   : this->gsh;
  this->qqtT = new QQt(this->gshT);

  this->o_prop = platform->device.malloc<dfloat>(2 * this->fieldOffsetSum);
  this->o_diff = this->o_prop.slice(0 * this->fieldOffsetSum);
  this->o_rho = this->o_prop.slice(1 * this->fieldOffsetSum);

  for (int is = 0; is < this->NSfields; is++) {
    const std::string sid = scalarDigitStr(is);

    if (options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
      continue;
    }

    auto mesh = (is) ? this->meshV : this->mesh[0]; // only first scalar can be a CHT mesh

    dfloat diff = 1;
    dfloat rho = 1;
    options.getArgs("SCALAR" + sid + " DIFFUSIVITY", diff);
    options.getArgs("SCALAR" + sid + " DENSITY", rho);

    auto o_diff = this->o_diff + this->fieldOffsetScan[is];
    auto o_rho = this->o_rho + this->fieldOffsetScan[is];
    platform->linAlg->fill(mesh->Nlocal, diff, o_diff);
    platform->linAlg->fill(mesh->Nlocal, rho, o_rho);
  }

  this->anyCvodeSolver = false;
  this->anyEllipticSolver = false;

  this->EToBOffset = this->mesh[0]->Nelements * this->mesh[0]->Nfaces;

  this->EToB = (int *)calloc(this->EToBOffset * this->NSfields, sizeof(int));

  for (int is = 0; is < this->NSfields; is++) {
    std::string sid = scalarDigitStr(is);

    this->compute[is] = 1;
    if (options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
      this->compute[is] = 0;
      this->cvodeSolve[is] = 0;
      continue;
    }

    this->cvodeSolve[is] = options.compareArgs("SCALAR" + sid + " SOLVER", "CVODE");
    this->anyCvodeSolver |= this->cvodeSolve[is];
    this->anyEllipticSolver |= (!this->cvodeSolve[is] && this->compute[is]);

    auto mesh = (is) ? this->meshV : this->mesh[0]; // only first scalar can be a CHT mesh

    int cnt = 0;
    for (int e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        this->EToB[cnt + this->EToBOffset * is] = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], "scalar" + sid);
        cnt++;
      }
    }
  }
  this->o_EToB = platform->device.malloc<int>(this->EToBOffset * this->NSfields, this->EToB);

  this->o_compute = platform->device.malloc<dlong>(this->NSfields, this->compute.data());
  this->o_cvodeSolve = platform->device.malloc<dlong>(this->NSfields, this->cvodeSolve.data());

  int nFieldsAlloc = this->anyEllipticSolver ? std::max(this->nBDF, this->nEXT) : 1;
  this->o_S = platform->device.malloc<dfloat>(nFieldsAlloc * this->fieldOffsetSum);

  nFieldsAlloc = this->anyEllipticSolver ? this->nEXT : 1;
  this->o_NLT = platform->device.malloc<dfloat>(nFieldsAlloc * this->fieldOffsetSum);

  if (this->anyEllipticSolver) {
    this->o_Se = platform->device.malloc<dfloat>(this->fieldOffsetSum);
    this->o_JwF = platform->device.malloc<dfloat>(this->fieldOffsetSum);
  }

  bool scalarFilteringEnabled = false;
  bool avmEnabled = false;
  {
    for (int is = 0; is < this->NSfields; is++) {
      std::string sid = scalarDigitStr(is);

      if (!options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "NONE")) {
        scalarFilteringEnabled = true;
      }

      if (options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "AVM_AVERAGED_MODAL_DECAY")) {
        avmEnabled = true;
      }
    }
  }

  this->applyFilter = 0;

  if (scalarFilteringEnabled) {

    std::vector<dlong> applyFilterRT(this->NSfields, 0);
    const dlong Nmodes = this->mesh[0]->N + 1;
    this->o_filterRT = platform->device.malloc<dfloat>(this->NSfields * Nmodes * Nmodes);
    this->o_filterS = platform->device.malloc<dfloat>(this->NSfields);
    this->o_applyFilterRT = platform->device.malloc<dlong>(this->NSfields);
    for (int is = 0; is < this->NSfields; is++) {
      std::string sid = scalarDigitStr(is);

      if (options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "NONE")) {
        continue;
      }
      if (!this->compute[is]) {
        continue;
      }

      if (options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "HPFRT")) {
        int filterNc = -1;
        options.getArgs("SCALAR" + sid + " HPFRT MODES", filterNc);
        dfloat filterS;
        options.getArgs("SCALAR" + sid + " HPFRT STRENGTH", filterS);
        filterS = -1.0 * fabs(filterS);
        this->filterS[is] = filterS;

        this->o_filterRT.copyFrom(lowPassFilterSetup(this->mesh[is], filterNc),
                                  Nmodes * Nmodes,
                                  is * Nmodes * Nmodes);

        applyFilterRT[is] = 1;
        this->applyFilter = 1;
      }
    }

    this->o_filterS.copyFrom(this->filterS.data(), this->NSfields);
    this->o_applyFilterRT.copyFrom(applyFilterRT.data(), this->NSfields);

    if (avmEnabled) {
      avm::setup(meshV, gsh);
    }
  }

  std::string kernelName;
  const std::string suffix = "Hex3D";
  {
    kernelName = "strongAdvectionVolume" + suffix;
    this->strongAdvectionVolumeKernel = platform->kernelRequests.load(section + kernelName);

    if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
      kernelName = "strongAdvectionCubatureVolume" + suffix;
      this->strongAdvectionCubatureVolumeKernel = platform->kernelRequests.load(section + kernelName);
    }

    kernelName = "advectMeshVelocity" + suffix;
    this->advectMeshVelocityKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "maskCopy";
    this->maskCopyKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "maskCopy2";
    this->maskCopy2Kernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "neumannBC" + suffix;
    this->neumannBCKernel = platform->kernelRequests.load(section + kernelName);
    kernelName = "dirichletBC";
    this->dirichletBCKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "filterRT" + suffix;
    this->filterRTKernel = platform->kernelRequests.load("core-" + kernelName);

    if (this->Nsubsteps) {
      if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
        kernelName = "subCycleStrongCubatureVolume" + suffix;
        this->subCycleStrongCubatureVolumeKernel = platform->kernelRequests.load(section + kernelName);
      }
      kernelName = "subCycleStrongVolume" + suffix;
      this->subCycleStrongVolumeKernel = platform->kernelRequests.load(section + kernelName);
    }
  }

  if (this->anyCvodeSolver) {
    this->cvode = new cvode_t(this);
  }
}

void cds_t::makeNLT(int is, double time, int tstep, occa::memory &o_Usubcycling)
{
    const std::string sid = scalarDigitStr(is);

    auto mesh = (is) ? this->meshV : this->mesh[0];
    const dlong isOffset = this->fieldOffsetScan[is];

    if (platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "HPFRT")) {
      const auto fieldOffset = this->fieldOffset[0]; // same for all scalars
      this->filterRTKernel(this->meshV->Nelements,
                           is,
                           1,
                           this->o_fieldOffsetScan,
                           this->o_applyFilterRT,
                           this->o_filterRT,
                           this->o_filterS,
                           this->o_rho,
                           this->o_S,
                           this->o_NLT);

      double flops = 6 * mesh->Np * mesh->Nq + 4 * mesh->Np;
      flops *= static_cast<double>(mesh->Nelements);
      platform->flopCounter->add("scalarFilterRT", flops);
    }

    const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
    if (movingMesh && !this->Nsubsteps) {
      this->advectMeshVelocityKernel(this->meshV->Nelements,
                                     mesh->o_vgeo,
                                     mesh->o_D,
                                     isOffset,
                                     this->vFieldOffset,
                                     this->o_rho,
                                     mesh->o_U,
                                     this->o_S,
                                     this->o_NLT);
      double flops = 18 * mesh->Np * mesh->Nq + 21 * mesh->Np;
      flops *= static_cast<double>(mesh->Nelements);
      platform->flopCounter->add("scalar advectMeshVelocity", flops);
    }

    if (platform->options.compareArgs("ADVECTION", "TRUE")) {
      if (this->Nsubsteps) {
        o_Usubcycling = this->advectionSubcyling(std::min(tstep, this->nEXT), time, is);
      } else {
        if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
          this->strongAdvectionCubatureVolumeKernel(this->meshV->Nelements,
                                                    1,
                                                    0, /* weighted */
                                                    0, /* sharedRho */
                                                    mesh->o_vgeo,
                                                    mesh->o_cubDiffInterpT,
                                                    mesh->o_cubInterpT,
                                                    mesh->o_cubProjectT,
                                                    this->o_compute + is,
                                                    this->o_fieldOffsetScan + is,
                                                    this->vFieldOffset,
                                                    this->vCubatureOffset,
                                                    this->o_S,
                                                    this->o_Urst,
                                                    this->o_rho,
                                                    this->o_NLT);
        } else {
          this->strongAdvectionVolumeKernel(this->meshV->Nelements,
                                            1,
                                            0, /* weighted */
                                            mesh->o_vgeo,
                                            mesh->o_D,
                                            this->o_compute + is,
                                            this->o_fieldOffsetScan + is,
                                            this->vFieldOffset,
                                            this->o_S,
                                            this->o_Urst,
                                            this->o_rho,
                                            this->o_NLT);
        }
        advectionFlops(this->mesh[0], 1);
      }
    }
}

void cds_t::saveSolutionState()
{
  if (!o_Ssave.isInitialized()) {
    o_Ssave = platform->device.malloc<dfloat>(o_S.length());
    o_NLTsave = platform->device.malloc<dfloat>(o_NLT.length());
    o_Spropsave = platform->device.malloc<dfloat>(o_prop.length());
  }

  o_Ssave.copyFrom(o_S, o_S.length());
  o_NLTsave.copyFrom(o_NLT, o_NLT.length());
  o_Spropsave.copyFrom(o_prop, o_prop.length());
}

void cds_t::restoreSolutionState()
{
  o_Ssave.copyTo(o_S, o_S.length());
  o_NLTsave.copyTo(o_NLT, o_NLT.length());
  o_Spropsave.copyTo(o_prop, o_prop.length());
}

void cds_t::applyAVM()
{
  auto verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  auto mesh = this->meshV; // assumes mesh is the same for all scalars
  static std::vector<occa::memory> o_diff0(NSfields);

  static std::vector<occa::memory> o_nuAVM;
  static auto initialized = false;
  if (!initialized) {
    for (int is = 0; is < NSfields; is++) {
      const auto sid = scalarDigitStr(is);

      if (platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "AVM_AVERAGED_MODAL_DECAY")) {
        nekrsCheck(mesh->N < 5,
                   platform->comm.mpiComm,
                   EXIT_FAILURE,
                   "%s\n",
                   "AVM requires polynomialOrder >= 5!");

        o_diff0[is] = platform->device.malloc<dfloat>(fieldOffset[is]);
        o_diff0[is].copyFrom(o_diff, fieldOffset[is], 0, fieldOffsetScan[is]);
      }
    }
    initialized = true;
  }

  for (int scalarIndex = 0; scalarIndex < NSfields; scalarIndex++) {
    const auto sid = scalarDigitStr(scalarIndex);
 
    if (!platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "AVM_AVERAGED_MODAL_DECAY"))
      continue; 

    // restore inital viscosity
    o_diff.copyFrom(o_diff0.at(scalarIndex),
                    fieldOffset[scalarIndex],
                    fieldOffsetScan[scalarIndex]);
 
    dfloat kappa = 1.0;
    platform->options.getArgs("SCALAR" + sid + " REGULARIZATION AVM ACTIVATION WIDTH", kappa);
 
    dfloat logS0 = 2.0; // threshold smoothness exponent (activate for logSk > logS0 - kappa) 
    platform->options.getArgs("SCALAR" + sid + " REGULARIZATION AVM DECAY THRESHOLD", logS0);
 
    dfloat scalingCoeff = 1.0;
    platform->options.getArgs("SCALAR" + sid + " REGULARIZATION AVM SCALING COEFF", scalingCoeff);

    dfloat absTol = 0;
    platform->options.getArgs("SCALAR" + sid + " REGULARIZATION AVM ABSOLUTE TOL", absTol);

    const bool makeCont = platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION AVM C0", "TRUE");
 
    auto o_Si = o_S.slice(fieldOffsetScan[scalarIndex], mesh->Nlocal);
    auto o_eps = avm::viscosity(vFieldOffset, o_U, o_Si, absTol, scalingCoeff, logS0, kappa, makeCont);
 
    if (verbose) {
      const dfloat maxEps = platform->linAlg->max(mesh->Nlocal, o_eps, platform->comm.mpiComm);
      const dfloat minEps = platform->linAlg->min(mesh->Nlocal, o_eps, platform->comm.mpiComm);
      occa::memory o_S_slice = o_diff + fieldOffsetScan[scalarIndex];
      const dfloat maxDiff = platform->linAlg->max(mesh->Nlocal, o_S_slice, platform->comm.mpiComm);
      const dfloat minDiff = platform->linAlg->min(mesh->Nlocal, o_S_slice, platform->comm.mpiComm);
 
      if (platform->comm.mpiRank == 0) {
        printf("applying a min/max artificial viscosity of (%f,%f) to scalar%s with min/max visc (%f,%f)\n",
               minEps,
               maxEps,
               sid.c_str(),
               minDiff,
               maxDiff);
      }
    }
 
    platform->linAlg->axpby(mesh->Nlocal, 1.0, o_eps, 1.0, o_diff, 0, fieldOffsetScan[scalarIndex]);
 
    if (verbose) {
      occa::memory o_S_slice = o_diff + fieldOffsetScan[scalarIndex];
      const dfloat maxDiff = platform->linAlg->max(mesh->Nlocal, o_S_slice, platform->comm.mpiComm);
      const dfloat minDiff = platform->linAlg->min(mesh->Nlocal, o_S_slice, platform->comm.mpiComm);
 
      if (platform->comm.mpiRank == 0) {
        printf("scalar%s now has a min/max visc: (%f,%f)\n", sid.c_str(), minDiff, maxDiff);
      }
    }
  }
}
