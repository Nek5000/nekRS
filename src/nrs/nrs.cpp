#include "nrs.hpp"
#include "bdry.hpp"
#include "bcMap.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "lowPassFilter.hpp"
#include "avm.hpp"
#include "re2Reader.hpp"
#include "cds.hpp"
#include "advectionSubCycling.hpp"
#include "elliptic.hpp"
#include "createEToBV.hpp"
#include "iofldFactory.hpp"

static void computeDivUErr(nrs_t *nrs, dfloat &divUErrVolAvg, dfloat &divUErrL2)
{
  auto mesh = nrs->mesh;

  auto o_divErr = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);

  nrs->divergenceVolumeKernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, nrs->fieldOffset, nrs->o_U, o_divErr);

  double flops = 18 * (mesh->Np * mesh->Nq + mesh->Np);
  flops *= static_cast<double>(mesh->Nelements);

  platform->flopCounter->add("divergenceVolumeKernel", flops);

  oogs::startFinish(o_divErr, 1, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_invLMM, o_divErr);

  platform->linAlg->axpby(mesh->Nlocal, 1.0, nrs->o_div, -1.0, o_divErr);
  divUErrL2 = platform->linAlg->weightedNorm2(mesh->Nlocal, mesh->o_LMM, o_divErr, platform->comm.mpiComm) /
              sqrt(mesh->volume);

  divUErrVolAvg =
      platform->linAlg->innerProd(mesh->Nlocal, mesh->o_LMM, o_divErr, platform->comm.mpiComm) / mesh->volume;

  divUErrVolAvg = std::abs(divUErrVolAvg);
}

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

static void assignKernels(nrs_t *nrs)
{
  std::string kernelName;
  const std::string suffix = "Hex3D";
  {
    const std::string section = "nrs-";

    kernelName = "computeFieldDotNormal";
    nrs->computeFieldDotNormalKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "computeFaceCentroid";
    nrs->computeFaceCentroidKernel = platform->kernelRequests.load(section + kernelName);

    {
      kernelName = "strongAdvectionVolume" + suffix;
      nrs->strongAdvectionVolumeKernel = platform->kernelRequests.load(section + kernelName);

      if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
        kernelName = "strongAdvectionCubatureVolume" + suffix;
        nrs->strongAdvectionCubatureVolumeKernel = platform->kernelRequests.load(section + kernelName);
      }
    }

    kernelName = "SijOij" + suffix;
    nrs->SijOijKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "curl" + suffix;
    nrs->curlKernel = platform->kernelRequests.load("core-" + kernelName);

    kernelName = "gradientVolume" + suffix;
    nrs->gradientVolumeKernel = platform->kernelRequests.load("core-" + kernelName);

    kernelName = "wGradientVolume" + suffix;
    nrs->wgradientVolumeKernel = platform->kernelRequests.load("core-" + kernelName);

    kernelName = "wDivergenceVolume" + suffix;
    nrs->wDivergenceVolumeKernel = platform->kernelRequests.load("core-" + kernelName);

    kernelName = "divergenceVolume" + suffix;
    nrs->divergenceVolumeKernel = platform->kernelRequests.load("core-" + kernelName);

    kernelName = "vectorFilterRT" + suffix;
    nrs->filterRTKernel = platform->kernelRequests.load("core-" + kernelName);

    kernelName = "sumMakef";
    nrs->sumMakefKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "divergenceSurface" + suffix;
    nrs->divergenceSurfaceKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "advectMeshVelocity" + suffix;
    nrs->advectMeshVelocityKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "pressureRhs" + suffix;
    nrs->pressureRhsKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "pressureStress" + suffix;
    nrs->pressureStressKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "pressureDirichletBC" + suffix;
    nrs->pressureDirichletBCKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "velocityRhs" + suffix;
    nrs->velocityRhsKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "averageNormalBcType";
    nrs->averageNormalBcTypeKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "fixZeroNormalMask";
    nrs->fixZeroNormalMaskKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "applyZeroNormalMask";
    nrs->applyZeroNormalMaskKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "initializeZeroNormalMask";
    nrs->initializeZeroNormalMaskKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "velocityDirichletBC" + suffix;
    nrs->velocityDirichletBCKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "velocityNeumannBC" + suffix;
    nrs->velocityNeumannBCKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "UrstCubature" + suffix;
    nrs->UrstCubatureKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "Urst" + suffix;
    nrs->UrstKernel = platform->kernelRequests.load(section + kernelName);

    if (nrs->Nsubsteps) {
      if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
        kernelName = "subCycleStrongCubatureVolume" + suffix;
        nrs->subCycleStrongCubatureVolumeKernel = platform->kernelRequests.load(section + kernelName);
      }
      kernelName = "subCycleStrongVolume" + suffix;
      nrs->subCycleStrongVolumeKernel = platform->kernelRequests.load(section + kernelName);
    }

    kernelName = "extrapolate";
    nrs->extrapolateKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "maskCopy";
    nrs->maskCopyKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "maskCopy2";
    nrs->maskCopy2Kernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "cfl" + suffix;
    nrs->cflKernel = platform->kernelRequests.load(section + kernelName);

    kernelName = "pressureAddQtl";
    nrs->pressureAddQtlKernel = platform->kernelRequests.load(section + kernelName);
  }
}

static void setupEllipticSolvers(nrs_t *nrs)
{
  auto createEToB = [](std::string field, mesh_t *mesh) {
    std::vector<int> EToB(mesh->Nelements * mesh->Nfaces);
    for (dlong e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        const int bID = mesh->EToB[f + e * mesh->Nfaces];
        EToB[f + e * mesh->Nfaces] = bcMap::ellipticType(bID, field);
      }
    }
    return EToB;
  };

  if (nrs->Nscalar) {
    cds_t *cds = nrs->cds;

    for (int is = 0; is < cds->NSfields; is++) {
      std::string sid = scalarDigitStr(is);

      if (!cds->compute[is]) {
        continue;
      }

      auto mesh = cds->mesh[is];

      const auto solverName = cds->cvodeSolve[is] ? "CVODE" : "ELLIPTIC";
      if (platform->comm.mpiRank == 0) {
        std::cout << "================= " << solverName << " SETUP SCALAR" << sid << " ===============\n";
      }

      const int nbrBIDs = bcMap::size("scalar" + sid);
      for (int bID = 1; bID <= nbrBIDs; bID++) {
        std::string bcTypeText(bcMap::text(bID, "scalar" + sid));
        if (platform->comm.mpiRank == 0 && bcTypeText.size()) {
          printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str());
        }
      }

      if (cds->cvodeSolve[is]) {
        continue;
      }

      auto EToB = createEToB("scalar" + sid, mesh);

      auto o_rho_i = cds->o_rho.slice(cds->fieldOffsetScan[is], mesh->Nlocal);
      auto o_lambda0 = cds->o_diff.slice(cds->fieldOffsetScan[is], mesh->Nlocal);
      auto o_lambda1 = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
      platform->linAlg->axpby(mesh->Nlocal, *cds->g0 / cds->dt[0], o_rho_i, 0.0, o_lambda1);

      cds->solver[is] = new elliptic("scalar" + sid, mesh, nrs->fieldOffset, EToB, o_lambda0, o_lambda1);
    }
  }

  if (nrs->flow) {
    auto mesh = nrs->mesh;

    if (platform->comm.mpiRank == 0) {
      printf("================ ELLIPTIC SETUP VELOCITY ================\n");
    }

    if (platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE")) {
      platform->options.setArgs("VELOCITY NFIELDS", std::to_string(mesh->dim));
    }

    bool unalignedBoundary = bcMap::unalignedMixedBoundary("velocity");
    nekrsCheck(unalignedBoundary && !platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE"),
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s\n",
               "SHL or unaligned SYM boundaries require solver = pcg+block");

    for (int bID = 1; bID <= bcMap::size("velocity"); bID++) {
      std::string bcTypeText(bcMap::text(bID, "velocity"));
      if (platform->comm.mpiRank == 0 && bcTypeText.size()) {
        printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str());
      }
    }

    auto o_lambda0 = nrs->o_mue;
    auto o_lambda1 = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
    platform->linAlg->axpby(mesh->Nlocal, nrs->g0 / nrs->dt[0], nrs->o_rho, 0.0, o_lambda1);

    auto EToBx = createEToB("x-velocity", mesh);
    auto EToBy = createEToB("y-velocity", mesh);
    auto EToBz = createEToB("z-velocity", mesh);

    if (platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE")) {
      std::vector<int> EToB;
      EToB.insert(std::end(EToB), std::begin(EToBx), std::end(EToBx));
      EToB.insert(std::end(EToB), std::begin(EToBy), std::end(EToBy));
      EToB.insert(std::end(EToB), std::begin(EToBz), std::end(EToBz));

      nrs->uvwSolver = new elliptic("velocity", mesh, nrs->fieldOffset, EToB, o_lambda0, o_lambda1);

      if (unalignedBoundary) {
        nrs->o_zeroNormalMaskVelocity =
            platform->device.malloc<dfloat>(nrs->uvwSolver->Nfields() * nrs->uvwSolver->fieldOffset());
        nrs->o_EToBVVelocity = platform->device.malloc<int>(nrs->mesh->Nlocal);
        createEToBV(nrs->mesh, nrs->uvwSolver->EToB(), nrs->o_EToBVVelocity);
        auto o_EToB = platform->device.malloc<int>(mesh->Nelements * mesh->Nfaces * nrs->uvwSolver->Nfields(),
                                                   nrs->uvwSolver->EToB().data());
        createZeroNormalMask(nrs, mesh, o_EToB, nrs->o_EToBVVelocity, nrs->o_zeroNormalMaskVelocity);

        auto applyZeroNormalMaskLambda =
            [nrs, mesh](dlong Nelements, const occa::memory &o_elementList, occa::memory &o_x) {
              applyZeroNormalMask(nrs,
                                  mesh,
                                  Nelements,
                                  o_elementList,
                                  nrs->uvwSolver->o_EToB(),
                                  nrs->o_zeroNormalMaskVelocity,
                                  o_x);
            };
        nrs->uvwSolver->applyZeroNormalMask(applyZeroNormalMaskLambda);
      }
    } else {
      nrs->uSolver = new elliptic("velocity", mesh, nrs->fieldOffset, EToBx, o_lambda0, o_lambda1);
      nrs->vSolver = new elliptic("velocity", mesh, nrs->fieldOffset, EToBy, o_lambda0, o_lambda1);
      nrs->wSolver = new elliptic("velocity", mesh, nrs->fieldOffset, EToBz, o_lambda0, o_lambda1);
    }
  }

  if (nrs->flow) {
    auto mesh = nrs->mesh;

    if (platform->comm.mpiRank == 0) {
      printf("================ ELLIPTIC SETUP PRESSURE ================\n");
    }

    auto o_lambda0 = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
    platform->linAlg->adyz(mesh->Nlocal, 1.0, nrs->o_rho, o_lambda0);

    auto EToB = createEToB("pressure", mesh);

    nrs->pSolver = new elliptic("pressure", mesh, nrs->fieldOffset, EToB, o_lambda0, o_NULL);
    if (nrs->cds) {
      nrs->cds->dpdt = nrs->pSolver->nullSpace();
    }
  }

  if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
    auto mesh = (nrs->cht) ? nrs->cds->mesh[0] : nrs->mesh;

    if (platform->comm.mpiRank == 0) {
      printf("================ ELLIPTIC SETUP MESH ================\n");
    }

    if (platform->options.compareArgs("MESH BLOCK SOLVER", "TRUE")) {
      platform->options.setArgs("MESH NFIELDS", std::to_string(mesh->dim));
    }

    std::vector<int> EToB;
    auto EToBx = createEToB("x-mesh", mesh);
    EToB.insert(std::end(EToB), std::begin(EToBx), std::end(EToBx));
    auto EToBy = createEToB("y-mesh", mesh);
    EToB.insert(std::end(EToB), std::begin(EToBy), std::end(EToBy));
    auto EToBz = createEToB("z-mesh", mesh);
    EToB.insert(std::end(EToB), std::begin(EToBz), std::end(EToBz));

    auto o_lambda0 = nrs->o_meshMue;

    nrs->meshSolver = new elliptic("mesh", mesh, nrs->fieldOffset, EToB, o_lambda0, o_NULL);

    bool unalignedBoundary = bcMap::unalignedMixedBoundary("mesh");
    if (unalignedBoundary) {
      nrs->o_zeroNormalMaskMeshVelocity =
          platform->device.malloc<dfloat>(nrs->meshSolver->Nfields() * nrs->meshSolver->fieldOffset());
      nrs->o_EToBVMeshVelocity = platform->device.malloc<int>(mesh->Nlocal);
      auto o_EToB = platform->device.malloc<int>(mesh->Nelements * mesh->Nfaces * nrs->meshSolver->Nfields(),
                                                 nrs->meshSolver->EToB().data());
      createEToBV(mesh, nrs->meshSolver->EToB(), nrs->o_EToBVMeshVelocity);
      createZeroNormalMask(nrs, mesh, o_EToB, nrs->o_EToBVMeshVelocity, nrs->o_zeroNormalMaskMeshVelocity);
      auto applyZeroNormalMaskLambda =
          [nrs, mesh](dlong Nelements, const occa::memory &o_elementList, occa::memory &o_x) {
            applyZeroNormalMask(nrs,
                                mesh,
                                Nelements,
                                o_elementList,
                                nrs->meshSolver->o_EToB(),
                                nrs->o_zeroNormalMaskMeshVelocity,
                                o_x);
          };
      nrs->meshSolver->applyZeroNormalMask(applyZeroNormalMaskLambda);
    }
  }
}

int nrs_t::numberActiveFields()
{
  int Nscalar = 0;
  platform->options.getArgs("NUMBER OF SCALARS", Nscalar);

  int fields = 0;
  if (!platform->options.compareArgs("VELOCITY SOLVER", "NONE")) {
    fields++;
  }
  for (int is = 0; is < Nscalar; ++is) {
    std::string sid = scalarDigitStr(is);
    if (!platform->options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
      fields++;
    }
  }
  return fields;
}

occa::memory nrs_t::advectionSubcycling(int nEXT, double time)
{
  const auto movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  const auto gsh = this->gsh;

  const auto nFields = this->NVfields;
  const auto fieldOffset = this->fieldOffset;
  const auto fieldOffsetSum = nFields * fieldOffset;
  const auto meshOffset =
      this->fieldOffset; // used for mesh->o_invLMM and mesh->o_divU in case of moving mesh

  auto o_U = this->o_U;

  auto o_Urst = (movingMesh) ? this->o_relUrst : this->o_Urst;
  auto opKernel = (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
                      ? this->subCycleStrongCubatureVolumeKernel
                      : this->subCycleStrongVolumeKernel;

  return advectionSubcyclingRK(mesh,
                               mesh,
                               time,
                               this->dt,
                               this->Nsubsteps,
                               this->coeffBDF,
                               nEXT,
                               nFields,
                               opKernel,
                               gsh,
                               meshOffset,
                               fieldOffset,
                               this->cubatureOffset,
                               fieldOffsetSum,
                               o_Urst,
                               o_U);
}

void nrs_t::printMinMax()
{
  if (platform->comm.mpiRank == 0) {
    printf("================= INITIAL CONDITION ====================\n");
  }

  {
    auto mesh = (cht) ? cds->mesh[0] : this->mesh;
    auto o_x = mesh->o_x;
    auto o_y = mesh->o_y;
    auto o_z = mesh->o_z;

    const auto xMin = platform->linAlg->min(mesh->Nlocal, o_x, platform->comm.mpiComm);
    const auto yMin = platform->linAlg->min(mesh->Nlocal, o_y, platform->comm.mpiComm);
    const auto zMin = platform->linAlg->min(mesh->Nlocal, o_z, platform->comm.mpiComm);
    const auto xMax = platform->linAlg->max(mesh->Nlocal, o_x, platform->comm.mpiComm);
    const auto yMax = platform->linAlg->max(mesh->Nlocal, o_y, platform->comm.mpiComm);
    const auto zMax = platform->linAlg->max(mesh->Nlocal, o_z, platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("XYZ   min/max: %g %g  %g %g  %g %g\n", xMin, xMax, yMin, yMax, zMin, zMax);
    }
  }

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    auto mesh = (cht) ? cds->mesh[0] : this->mesh;
    auto o_ux = mesh->o_U + 0 * this->fieldOffset;
    auto o_uy = mesh->o_U + 1 * this->fieldOffset;
    auto o_uz = mesh->o_U + 2 * this->fieldOffset;
    const auto uxMin = platform->linAlg->min(mesh->Nlocal, o_ux, platform->comm.mpiComm);
    const auto uyMin = platform->linAlg->min(mesh->Nlocal, o_uy, platform->comm.mpiComm);
    const auto uzMin = platform->linAlg->min(mesh->Nlocal, o_uz, platform->comm.mpiComm);
    const auto uxMax = platform->linAlg->max(mesh->Nlocal, o_ux, platform->comm.mpiComm);
    const auto uyMax = platform->linAlg->max(mesh->Nlocal, o_uy, platform->comm.mpiComm);
    const auto uzMax = platform->linAlg->max(mesh->Nlocal, o_uz, platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("UMSH  min/max: %g %g  %g %g  %g %g\n", uxMin, uxMax, uyMin, uyMax, uzMin, uzMax);
    }
  }

  {
    auto o_ux = this->o_U + 0 * this->fieldOffset;
    auto o_uy = this->o_U + 1 * this->fieldOffset;
    auto o_uz = this->o_U + 2 * this->fieldOffset;
    const auto uxMin = platform->linAlg->min(mesh->Nlocal, o_ux, platform->comm.mpiComm);
    const auto uyMin = platform->linAlg->min(mesh->Nlocal, o_uy, platform->comm.mpiComm);
    const auto uzMin = platform->linAlg->min(mesh->Nlocal, o_uz, platform->comm.mpiComm);
    const auto uxMax = platform->linAlg->max(mesh->Nlocal, o_ux, platform->comm.mpiComm);
    const auto uyMax = platform->linAlg->max(mesh->Nlocal, o_uy, platform->comm.mpiComm);
    const auto uzMax = platform->linAlg->max(mesh->Nlocal, o_uz, platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("U     min/max: %g %g  %g %g  %g %g\n", uxMin, uxMax, uyMin, uyMax, uzMin, uzMax);
    }
  }

  {
    const auto prMin = platform->linAlg->min(mesh->Nlocal, this->o_P, platform->comm.mpiComm);
    const auto prMax = platform->linAlg->max(mesh->Nlocal, this->o_P, platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("P     min/max: %g %g\n", prMin, prMax);
    }
  }

  if (this->Nscalar) {
    auto cds = this->cds;
    if (platform->comm.mpiRank == 0) {
      printf("S     min/max:");
    }

    int cnt = 0;
    for (int is = 0; is < cds->NSfields; is++) {
      cnt++;

      auto mesh = cds->mesh[is];

      auto o_si = this->cds->o_S + this->cds->fieldOffsetScan[is];
      const auto siMin = platform->linAlg->min(mesh->Nlocal, o_si, platform->comm.mpiComm);
      const auto siMax = platform->linAlg->max(mesh->Nlocal, o_si, platform->comm.mpiComm);
      if (platform->comm.mpiRank == 0) {
        if (cnt > 1) {
          printf("  ");
        } else {
          printf(" ");
        }
        printf("%g %g", siMin, siMax);
      }
    }
    if (platform->comm.mpiRank == 0) {
      printf("\n");
    }
  }
}

void nrsSetDefaultSettings(setupAide *options)
{
  options->setArgs("CONSTANT FLOW RATE", "FALSE");

  options->setArgs("NUMBER OF SCALARS", "0");

  options->setArgs("BDF ORDER", "2");
  options->setArgs("EXT ORDER", "3");
  options->setArgs("SUBCYCLING STEPS", "0");

  options->setArgs("ADVECTION", "TRUE");
  options->setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");

  options->setArgs("MIN ADJUST DT RATIO", "0.5");
  options->setArgs("MAX ADJUST DT RATIO", "1.5");

  options->setArgs("MESH SOLVER", "NONE");
  options->setArgs("MOVING MESH", "FALSE");

  options->setArgs("PRESSURE VISCOUS TERMS", "TRUE");

  options->setArgs("VARIABLE DT", "FALSE");
}

nrs_t::nrs_t()
{
  platform->solver = this;

  int result = 0;
  MPI_Comm_compare(platform->comm.mpiCommParent, platform->comm.mpiComm, &result);
  this->multiSession = result;

  platform->options.getArgs("MESH DIMENSION", this->NVfields);
  platform->options.getArgs("ELEMENT TYPE", this->elementType);

}

void nrs_t::init()
{
  this->flow = 1;
  if (platform->options.compareArgs("VELOCITY SOLVER", "NONE")) {
    this->flow = 0;
  }

  if (this->flow) {
    if (platform->options.compareArgs("VELOCITY STRESSFORMULATION", "TRUE")) {
      platform->options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    }
  }

  this->cht = [&]() {
    int nelgt, nelgv;
    const std::string meshFile = platform->options.getArgs("MESH FILE");
    re2::nelg(meshFile, nelgt, nelgv, platform->comm.mpiComm);

    nekrsCheck(nelgt != nelgv && platform->options.compareArgs("MOVING MESH", "TRUE"),
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s\n",
               "Conjugate heat transfer not supported in a moving mesh!");

    return (nelgt > nelgv) ? 1 : 0;
  }();

  nekrsCheck((cht || platform->options.compareArgs("LOWMACH", "TRUE")) &&
                 !platform->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE"),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "conjugate heat transfer or lowMach requires a TEMPERATURE field");

  platform->options.getArgs("SUBCYCLING STEPS", this->Nsubsteps);

  // g0 * dt required for Helmholtz coefficient (eigenvalues for Chebyshev during setup)
  this->g0 = 1;
  this->dt[0] = 1e-2;
  platform->options.getArgs("DT", this->dt[0]);

  platform->options.getArgs("BDF ORDER", this->nBDF);
  platform->options.getArgs("EXT ORDER", this->nEXT);
  if (this->Nsubsteps) {
    this->nEXT = this->nBDF;
  }

  nekrsCheck(this->nEXT < this->nBDF,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "EXT order needs to be >= BDF order!");

  platform->options.getArgs("NUMBER OF SCALARS", this->Nscalar);

  int N, cubN;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  platform->options.getArgs("CUBATURE POLYNOMIAL DEGREE", cubN);

  nek::setup(numberActiveFields());

  auto getMesh = [&]() {
    auto [meshT, meshV] = createMesh(platform->comm.mpiComm, N, cubN, this->cht, platform->kernelInfo);
    if (!cht) {
      meshV = meshT;
    }

    auto offset = meshV->Np * (meshV->Nelements);
    offset = std::max(offset, meshT->Np * (meshT->Nelements));

    auto cubOffset = offset;
    if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
      cubOffset = std::max(cubOffset, meshV->Nelements * meshV->cubNp);
    }

    offset = alignStride<dfloat>(offset);
    cubOffset = alignStride<dfloat>(cubOffset);

    this->fieldOffset = offset;
    this->cubatureOffset = cubOffset;

    meshT->fieldOffset = this->fieldOffset;
    meshV->fieldOffset = this->fieldOffset;

    return std::pair{meshT, meshV};
  }();

  this->mesh = getMesh.second;
  this->meshV = this->mesh;
  auto meshT = getMesh.first;

  auto verifyBC = [&]() {
    auto fields = nrsFieldsToSolve(platform->options);

    for (const auto &field : fields) {
      auto msh = (this->cht && (field == "scalar00" || field == "mesh")) ? meshT : mesh;
      nekrsCheck(msh->Nbid != bcMap::size(field),
                 platform->comm.mpiComm,
                 EXIT_FAILURE,
                 "Size of %s boundaryTypeMap (%d) does not match number of boundary IDs in mesh (%d)!\n",
                 field.c_str(),
                 bcMap::size(field),
                 msh->Nbid);
    }

    bcMap::checkBoundaryAlignment(mesh);
  };

  verifyBC();

  this->coeffEXT = (dfloat *)calloc(this->nEXT, sizeof(dfloat));
  this->coeffBDF = (dfloat *)calloc(this->nBDF, sizeof(dfloat));

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    auto mesh = meshT;
    const int nBDF = std::max(this->nBDF, this->nEXT);
    auto o_tmp = platform->device.malloc<dfloat>(mesh->Nlocal);
    o_tmp.copyFrom(mesh->o_LMM, mesh->Nlocal);
    mesh->o_LMM.free();
    mesh->o_LMM = platform->device.malloc<dfloat>(this->fieldOffset * nBDF);
    mesh->o_LMM.copyFrom(o_tmp, mesh->Nlocal);

    o_tmp.copyFrom(mesh->o_invLMM, mesh->Nlocal);
    mesh->o_invLMM.free();
    mesh->o_invLMM = platform->device.malloc<dfloat>(this->fieldOffset * nBDF);
    mesh->o_invLMM.copyFrom(o_tmp, mesh->Nlocal);

    const int nAB = std::max(this->nEXT, mesh->nAB);
    mesh->o_U = platform->device.malloc<dfloat>(this->NVfields * nAB * this->fieldOffset);
    mesh->o_Ue = platform->device.malloc<dfloat>(this->NVfields * nAB * this->fieldOffset);
    if (this->Nsubsteps) {
      mesh->o_divU = platform->device.malloc<dfloat>(this->fieldOffset * nAB);
    }
  }

  {
    const dlong Nstates = this->Nsubsteps ? std::max(this->nBDF, this->nEXT) : 1;
    bool useCVODE = platform->options.compareArgs("CVODE", "TRUE");
    if ((useCVODE || this->Nsubsteps) && platform->options.compareArgs("MOVING MESH", "TRUE")) {
      this->o_relUrst = platform->device.malloc<dfloat>(Nstates * this->NVfields * this->cubatureOffset);
    }
    if (!this->Nsubsteps || platform->options.compareArgs("MOVING MESH", "FALSE")) {
      this->o_Urst = platform->device.malloc<dfloat>(Nstates * this->NVfields * this->cubatureOffset);
    }
  }

  this->o_U =
      platform->device.malloc<dfloat>(this->NVfields * std::max(this->nBDF, this->nEXT) * this->fieldOffset);

  this->o_Ue = platform->device.malloc<dfloat>(this->NVfields * this->fieldOffset);

  this->o_P = platform->device.malloc<dfloat>(this->fieldOffset);

  this->o_JwF = platform->device.malloc<dfloat>(this->NVfields * this->fieldOffset);
  this->o_NLT = platform->device.malloc<dfloat>(this->NVfields * this->nEXT * this->fieldOffset);

  int nProperties = 2;
  if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
    nProperties = 4;
  }

  this->o_prop = platform->device.malloc<dfloat>(nProperties * this->fieldOffset);
  this->o_mue = this->o_prop.slice(0 * this->fieldOffset);
  this->o_rho = this->o_prop.slice(1 * this->fieldOffset);
  if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
    this->o_meshMue = this->o_prop.slice(2 * this->fieldOffset);
    this->o_meshRho = this->o_prop.slice(3 * this->fieldOffset);
  }

  dfloat mue = 1;
  dfloat rho = 1;
  platform->options.getArgs("VISCOSITY", mue);
  platform->options.getArgs("DENSITY", rho);

  platform->linAlg->fill(mesh->Nlocal, mue, this->o_mue);
  platform->linAlg->fill(mesh->Nlocal, rho, this->o_rho);
  if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
    auto o_mue = this->o_prop + 2 * this->fieldOffset;
    auto o_rho = this->o_prop + 3 * this->fieldOffset;
    platform->linAlg->fill(meshT->Nlocal, 1.0, o_mue);
    platform->linAlg->fill(meshT->Nlocal, 0.0, o_rho);
  }

  if (platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")) {
    this->o_Uc = platform->device.malloc<dfloat>(this->NVfields * this->fieldOffset);
    this->o_Pc = platform->device.malloc<dfloat>(mesh->Nlocal);
    this->o_prevProp = platform->device.malloc<dfloat>(2 * this->fieldOffset);
    this->o_prevProp.copyFrom(this->o_prop, this->o_prevProp.length());
  }

  this->o_div = platform->device.malloc<dfloat>(this->fieldOffset);

  this->o_coeffEXT = platform->device.malloc<dfloat>(this->nEXT, this->coeffEXT);
  this->o_coeffBDF = platform->device.malloc<dfloat>(this->nBDF, this->coeffBDF);

  this->gsh = oogs::setup(mesh->ogs, this->NVfields, this->fieldOffset, ogsDfloat, NULL, OOGS_AUTO);
  this->qqt = new QQt(this->gsh);
  this->qqtT = this->qqt;

  if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
    this->gshMesh = oogs::setup(meshT->ogs, this->NVfields, this->fieldOffset, ogsDfloat, NULL, OOGS_AUTO);
  }

  auto EToB = [&](const std::string &field, mesh_t *mesh) {
    std::vector<int> EToB(mesh->Nelements * mesh->Nfaces);
    int cnt = 0;
    for (int e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        EToB[cnt] = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], field);
        cnt++;
      }
    }
    auto o_EToB = platform->device.malloc<int>(EToB.size());
    o_EToB.copyFrom(EToB.data());
    return o_EToB;
  };

  if (this->flow) {
    this->o_EToB = EToB("velocity", mesh);
  }

  if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
    this->o_EToBMeshVelocity = EToB("mesh", meshT);
  }

  if (platform->options.compareArgs("VELOCITY REGULARIZATION METHOD", "HPFRT")) {
    int nModes = -1;
    dfloat strength = 1.0;
    platform->options.getArgs("VELOCITY HPFRT STRENGTH", strength);
    platform->options.getArgs("VELOCITY HPFRT MODES", nModes);
    this->filterS = -std::abs(strength);
    this->o_filterRT = lowPassFilterSetup(mesh, nModes);
  }

  assignKernels(this);

  if (this->Nscalar) {
    cdsConfig_t cfg;

    cfg.Nscalar = this->Nscalar;
    cfg.g0 = &this->g0;
    cfg.dt = this->dt;
    cfg.nBDF = this->nBDF;
    cfg.coeffBDF = this->coeffBDF;
    cfg.o_coeffBDF = this->o_coeffBDF;
    cfg.nEXT = this->nEXT;
    cfg.coeffEXT = this->coeffEXT;
    cfg.o_coeffEXT = this->o_coeffEXT;
    cfg.o_usrwrk = &this->o_usrwrk;
    cfg.vFieldOffset = this->fieldOffset;
    cfg.Nsubsteps = this->Nsubsteps;
    cfg.vCubatureOffset = this->cubatureOffset;
    cfg.fieldOffset = this->fieldOffset;
    cfg.o_U = this->o_U;
    cfg.o_Ue = this->o_Ue;
    cfg.o_Urst = this->o_Urst;
    cfg.o_relUrst = this->o_relUrst;
    cfg.meshT = meshT;
    cfg.meshV = mesh;
    cfg.dp0thdt = &this->dp0thdt;
    cfg.alpha0Ref = &this->alpha0Ref;

    this->cds = new cds_t(cfg);

    this->qqtT = this->cds->qqtT;

    if (this->cds->cvode) {
      this->cds->cvode->setEvaluateProperties(
          std::bind(&nrs_t::evaluateProperties, this, std::placeholders::_1));
      this->cds->cvode->setEvaluateDivergence(
          std::bind(&nrs_t::evaluateDivergence, this, std::placeholders::_1));
    }
  }

  this->setIC();

  if (this->cds) {
    if (this->userScalarSource) {
      cds->userSource = this->userScalarSource;
    }
    if (this->userProperties) {
      cds->userProperties = this->userProperties;
    }
    if (this->userScalarImplicitLinearTerm) {
      this->cds->userImplicitLinearTerm = this->userScalarImplicitLinearTerm;
    }
  }

  double startTime;
  platform->options.getArgs("START TIME", startTime);
  this->evaluateProperties(startTime);

  // CVODE can only be initialized once the initial condition is known
  if (cds) {
    if (cds->cvode) {
      cds->cvode->initialize();
    }
  }

  printMeshMetrics(meshT);
  if (mesh != meshT) {
    printMeshMetrics(mesh);
  }

  setupEllipticSolvers(this);
}

void nrs_t::restartFromFile(const std::string &restartStr)
{
  auto options = serializeString(restartStr, '+');
  const auto fileName = options[0];
  options.erase(options.begin());

  if (platform->comm.mpiRank == 0) {
    if (options.size()) {
      std::cout << "restart options: ";
    }
    for (const auto &element : options) {
      std::cout << element << "  ";
    }
    std::cout << std::endl;
  }

  auto requestedStep = [&]() {
    auto it = std::find_if(options.begin(), options.end(), [](const std::string &s) {
      return s.find("step") != std::string::npos;
    });

    std::string val;
    if (it != options.end()) {
      val = serializeString(*it, '=').at(1);
      options.erase(it);
    }
    return (val.empty()) ? -1 : std::stoi(val);
  }();

  auto requestedTime = [&]() {
    auto it = std::find_if(options.begin(), options.end(), [](const std::string &s) {
      return s.find("time") != std::string::npos;
    });

    std::string val;
    if (it != options.end()) {
      val = serializeString(*it, '=').at(1);
      options.erase(it);
    }
    return val;
  }();

  auto pointInterpolation = [&]() {
    auto it = std::find_if(options.begin(), options.end(), [](const std::string &s) {
      return s.find("int") != std::string::npos;
    });

    auto found = false;
    if (it != options.end()) {
      found = true;
      options.erase(it);
    }
    return found;
  }();

  const auto requestedFields = [&]() {
    std::vector<std::string> flds;
    for (const auto &entry : {"x", "u", "p", "t", "s"}) {
      auto it = std::find_if(options.begin(), options.end(), [entry](const std::string &s) {
        std::string ss = s;
        lowerCase(ss);
        return ss.find(entry) != std::string::npos;
      });
      if (it != options.end()) {
        std::string s = *it;
        lowerCase(s);
        std::cout << "requested field: " << s << std::endl;
        flds.push_back(s);
      }
    }
    return flds;
  }();

  auto fileNameEndsWithBp = [&]() {
    const std::string suffix = ".bp";
    if (fileName.size() >= suffix.size()) {
      return fileName.compare(fileName.size() - suffix.size(), suffix.size(), suffix) == 0;
    }
    return false;
  }();
  auto iofld = iofldFactory::create((fileNameEndsWithBp) ? "adios" : "");
  iofld->open((cht) ? cds->mesh[0] : mesh, iofld::mode::read, fileName, requestedStep);

  const auto avaiableFields = iofld->availableVariables();
  if (platform->comm.mpiRank == 0 && platform->verbose) {
    for (const auto &entry : avaiableFields) {
      std::cout << " found variable " << entry << std::endl;
    }
  }

  double time = -1;
  iofld->addVariable("time", time);
  if (platform->options.compareArgs("LOWMACH", "TRUE")) {
    iofld->addVariable("p0th", p0th[0]);
  }

  auto checkOption = [&](const std::string &name) {
    if (requestedFields.size() == 0) {
      return true; // nothing specfied -> assign all
    }
    if (std::find(requestedFields.begin(), requestedFields.end(), name) != requestedFields.end()) {
      return true;
    }
    return false;
  };

  auto isAvailable = [&](const std::string &name) {
    return std::find(avaiableFields.begin(), avaiableFields.end(), name) != avaiableFields.end();
  };

  if (checkOption("x") && isAvailable("mesh")) {
    std::vector<occa::memory> o_iofldX;
    auto mesh = (cht) ? cds->mesh[0] : this->mesh;
    o_iofldX.push_back(mesh->o_x);
    o_iofldX.push_back(mesh->o_y);
    o_iofldX.push_back(mesh->o_z);
    iofld->addVariable("mesh", o_iofldX);
  }

  if (checkOption("u") && isAvailable("velocity")) {
    std::vector<occa::memory> o_iofldU;
    o_iofldU.push_back(o_U.slice(0 * fieldOffset, mesh->Nlocal));
    o_iofldU.push_back(o_U.slice(1 * fieldOffset, mesh->Nlocal));
    o_iofldU.push_back(o_U.slice(2 * fieldOffset, mesh->Nlocal));
    iofld->addVariable("velocity", o_iofldU);
  }

  if (checkOption("p") && isAvailable("pressure")) {
    std::vector<occa::memory> o_iofldP = {o_P.slice(0, mesh->Nlocal)};
    iofld->addVariable("pressure", o_iofldP);
  }

  if (Nscalar) {
    std::vector<occa::memory> o_iofldT;
    if (checkOption("t") && isAvailable("temperature")) {
      auto mesh = (cht) ? cds->mesh[0] : this->mesh;
      o_iofldT.push_back(cds->o_S.slice(0, mesh->Nlocal));
      iofld->addVariable("temperature", o_iofldT);
    }

    const auto scalarStart = (o_iofldT.size()) ? 1 : 0;
    for (int i = scalarStart; i < Nscalar; i++) {
      const auto sid = scalarDigitStr(i - scalarStart);
      if (checkOption("s" + sid) && isAvailable("scalar" + sid)) {
        auto o_Si = cds->o_S.slice(cds->fieldOffsetScan[i], mesh->Nlocal);
        std::vector<occa::memory> o_iofldSi = {o_Si};
        iofld->addVariable("scalar" + sid, o_iofldSi);
      }
    }
  }

  if (pointInterpolation) {
    iofld->readAttribute("interpolate", "true");
  }

  iofld->process();
  iofld->close();

  platform->options.setArgs("START TIME", (requestedTime.size()) ? requestedTime : to_string_f(time));
}

void nrs_t::setIC()
{
  getICFromNek();

  if (!platform->options.getArgs("RESTART FILE NAME").empty()) {
    restartFromFile(platform->options.getArgs("RESTART FILE NAME"));
  }

  if (platform->comm.mpiRank == 0) {
    std::cout << "calling UDF_Setup ... \n" << std::flush;
  }
  udf.setup();
  if (platform->comm.mpiRank == 0) {
    std::cout << "done\n" << std::flush;
  }

  if (cht) {
    cds->mesh[0]->update();
  }
  mesh->update();

  auto projC0 = [&](oogs_t *gsh, mesh_t *mesh, int nFields, dlong fieldOffset, occa::memory &o_in) {
    platform->linAlg->axmyMany(mesh->Nlocal, nFields, fieldOffset, 0, 1.0, mesh->o_LMM, o_in);
    oogs::startFinish(o_in, nFields, fieldOffset, ogsDfloat, ogsAdd, gsh);
    platform->linAlg->axmyMany(mesh->Nlocal, nFields, fieldOffset, 0, 1.0, mesh->o_invLMM, o_in);
  };

  projC0(gsh, mesh, NVfields, fieldOffset, o_U);

  projC0(gsh, mesh, 1, fieldOffset, o_P);

  if (Nscalar) {
    for (int s = 0; s < Nscalar; ++s) {
      const std::string sid = scalarDigitStr(s);
      if (platform->options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
        continue;
      }
      auto gsh = (s == 0) ? cds->gshT : cds->gsh;
      auto o_Si = cds->o_S + cds->fieldOffsetScan[s];
      projC0(gsh, cds->mesh[s], 1, cds->fieldOffset[s], o_Si);
    }
  }

  double startTime;
  platform->options.getArgs("START TIME", startTime);
  copyToNek(startTime, 0, true); // ensure both codes are in sync

  nekrsCheck(platform->options.compareArgs("LOWMACH", "TRUE") && p0th[0] <= 1e-6,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Unreasonable p0th value %g!",
             p0th[0]);
}

void nrs_t::printRunStat(int step)
{
  const int rank = platform->comm.mpiRank;
  auto comm_ = platform->comm.mpiComm;

  platform->timer.set("velocity proj",
                      platform->timer.query("velocity proj pre", "DEVICE:MAX") +
                          platform->timer.query("velocity proj post", "DEVICE:MAX"),
                      platform->timer.count("velocity proj pre"));

  platform->timer.set("pressure proj",
                      platform->timer.query("pressure proj pre", "DEVICE:MAX") +
                          platform->timer.query("pressure proj post", "DEVICE:MAX"),
                      platform->timer.count("pressure proj pre"));

  platform->timer.set("scalar proj",
                      platform->timer.query("scalar proj pre", "DEVICE:MAX") +
                          platform->timer.query("scalar proj post", "DEVICE:MAX"),
                      platform->timer.count("scalar proj pre"));

  platform->timer.set("mesh proj",
                      platform->timer.query("mesh proj pre", "DEVICE:MAX") +
                          platform->timer.query("mesh proj post", "DEVICE:MAX"),
                      platform->timer.count("mesh proj pre"));

  double gsTime = ogsTime(/* reportHostTime */ true);
  MPI_Allreduce(MPI_IN_PLACE, &gsTime, 1, MPI_DOUBLE, MPI_MAX, comm_);

  const double tElapsedTime = platform->timer.query("elapsed", "DEVICE:MAX");

  if (rank == 0) {
    std::cout << "\n>>> runtime statistics (step= " << step << "  totalElapsed= " << tElapsedTime << "s"
              << "):\n";
  }

  std::cout.setf(std::ios::scientific);
  int outPrecisionSave = std::cout.precision();
  std::cout.precision(5);

  if (rank == 0) {
    std::cout << "name                    "
              << "time          "
              << "abs%  "
              << "rel%  "
              << "calls\n";
  }

  const double tElapsedTimeSolve = platform->timer.query("elapsedStepSum", "DEVICE:MAX");
  platform->timer.printStatSetElapsedTimeSolve(tElapsedTimeSolve);
  const double tSetup = platform->timer.query("setup", "DEVICE:MAX");

  const double tMinSolveStep = platform->timer.query("minSolveStep", "DEVICE:MAX");
  const double tMaxSolveStep = platform->timer.query("maxSolveStep", "DEVICE:MAX");

  const double tScalarCvode = platform->timer.query("cvode_t::solve", "DEVICE:MAX");

  bool printFlops = !platform->options.compareArgs("PRESSURE PRECONDITIONER", "SEMFEM") && tScalarCvode < 0;

  const double flops =
      platform->flopCounter->get(platform->comm.mpiComm) / (tElapsedTimeSolve * platform->comm.mpiCommSize);

  platform->timer.printStatEntry("  solve                 ", tElapsedTimeSolve, tElapsedTimeSolve);
  if (tElapsedTimeSolve > 0 && rank == 0) {
    std::cout << "    min                 " << tMinSolveStep << "s\n";
    std::cout << "    max                 " << tMaxSolveStep << "s\n";
    if (printFlops) {
      std::cout << "    flops/rank          " << flops << "\n";
    }
  }

  auto lpmLocalKernelPredicate = [](const std::string &tag) {
    return tag.find("lpm_t::") != std::string::npos && tag.find("localKernel") != std::string::npos;
  };

  auto lpmLocalEvalKernelPredicate = [](const std::string &tag) {
    return tag.find("lpm_t::") != std::string::npos && tag.find("localEvalKernel") != std::string::npos;
  };

  auto neknekLocalKernelPredicate = [](const std::string &tag) {
    return tag.find("neknek_t::") != std::string::npos && tag.find("localKernel") != std::string::npos;
  };

  auto neknekLocalEvalKernelPredicate = [](const std::string &tag) {
    return tag.find("neknek_t::") != std::string::npos && tag.find("localEvalKernel") != std::string::npos;
  };

  platform->timer.printStatEntry("    checkpointing       ",
                                 "checkpointing",
                                 "DEVICE:MAX",
                                 tElapsedTimeSolve);
  platform->timer.printStatEntry("    udfExecuteStep      ",
                                 "udfExecuteStep",
                                 "DEVICE:MAX",
                                 tElapsedTimeSolve);
  const double tudf = platform->timer.query("udfExecuteStep", "DEVICE:MAX");
  platform->timer.printStatEntry("      lpm integrate     ", "lpm_t::integrate", "DEVICE:MAX", tudf);
  const double tlpm = platform->timer.query("lpm_t::integrate", "DEVICE:MAX");
  platform->timer.printStatEntry("        userRHS         ", "lpm_t::integrate::userRHS", "DEVICE:MAX", tlpm);
  const double tParticleRHS = platform->timer.query("lpm_t::integrate::userRHS", "DEVICE:MAX");
  platform->timer.printStatEntry("          interpolate   ",
                                 "lpm_t::integrate::userRHS::interpolate",
                                 "DEVICE:MAX",
                                 tParticleRHS);
  const double tInterpPart = platform->timer.query("lpm_t::integrate::userRHS::interpolate", "DEVICE:MAX");
  auto [tLocalKernel, nLocalKernel] =
      platform->timer.sumAllMatchingTags(lpmLocalEvalKernelPredicate, "DEVICE:MAX");
  platform->timer.printStatEntry("            eval kernel ", tLocalKernel, nLocalKernel, tInterpPart);
  platform->timer.printStatEntry("        findpts         ", "lpm_t::integrate::find", "DEVICE:MAX", tlpm);
  const double tFindPart = platform->timer.query("lpm_t::integrate::find", "DEVICE:MAX");
  auto [tFindKernel, nFindKernel] = platform->timer.sumAllMatchingTags(lpmLocalKernelPredicate, "DEVICE:MAX");
  platform->timer.printStatEntry("          find kernel   ", tFindKernel, nFindKernel, tFindPart);
  platform->timer.printStatEntry("        delete          ", "lpm_t::deleteParticles", "DEVICE:MAX", tlpm);
  platform->timer.printStatEntry("      lpm add           ", "lpm_t::addParticles", "DEVICE:MAX", tudf);
  platform->timer.printStatEntry("      lpm write         ", "lpm_t::write", "DEVICE:MAX", tudf);

  const double tDiv = platform->timer.query("udfDiv", "DEVICE:MAX");
  platform->timer.printStatEntry("    udfDiv              ", "udfDiv", "DEVICE:MAX", tElapsedTimeSolve);

  const double tMakef = platform->timer.query("makef", "DEVICE:MAX");
  platform->timer.printStatEntry("    makef               ", "makef", "DEVICE:MAX", tElapsedTimeSolve);
  platform->timer.printStatEntry("      udfUEqnSource     ", "udfUEqnSource", "DEVICE:MAX", tMakef);

  const double tMakeq = platform->timer.query("makeq", "DEVICE:MAX");
  platform->timer.printStatEntry("    makeq               ", "makeq", "DEVICE:MAX", tElapsedTimeSolve);
  platform->timer.printStatEntry("      udfSEqnSource     ", "udfSEqnSource", "DEVICE:MAX", tMakeq);

  platform->timer.printStatEntry("    udfProperties       ",
                                 "udfProperties",
                                 "DEVICE:MAX",
                                 tElapsedTimeSolve);

  platform->timer.printStatEntry("    meshUpdate          ", "meshUpdate", "DEVICE:MAX", tElapsedTimeSolve);
  const double tMesh = platform->timer.query("meshSolve", "DEVICE:MAX");
  platform->timer.printStatEntry("    meshSolve           ", "meshSolve", "DEVICE:MAX", tElapsedTimeSolve);
  platform->timer.printStatEntry("      preconditioner    ", "mesh preconditioner", "DEVICE:MAX", tMesh);
  platform->timer.printStatEntry("      initial guess     ", "mesh proj", "DEVICE:MAX", tMesh);

  const double tNekNek = platform->timer.query("neknek update boundary", "DEVICE:MAX");
  platform->timer.printStatEntry("    neknek              ",
                                 "neknek update boundary",
                                 "DEVICE:MAX",
                                 tElapsedTimeSolve);
  platform->timer.printStatEntry("      sync              ", "neknek sync", "DEVICE:MAX", tNekNek);
  platform->timer.printStatEntry("      exchange          ", "neknek exchange", "DEVICE:MAX", tNekNek);
  const double tExchange = platform->timer.query("neknek exchange", "DEVICE:MAX");
  std::tie(tLocalKernel, nLocalKernel) =
      platform->timer.sumAllMatchingTags(neknekLocalEvalKernelPredicate, "DEVICE:MAX");
  platform->timer.printStatEntry("        eval kernel     ", tLocalKernel, nLocalKernel, tExchange);
  platform->timer.printStatEntry("      findpts           ",
                                 "neknek updateInterpPoints",
                                 "DEVICE:MAX",
                                 tNekNek);
  const double tFindpts = platform->timer.query("neknek updateInterpPoints", "DEVICE:MAX");

  if (tFindpts > 0.0) {
    std::tie(tFindKernel, nFindKernel) =
        platform->timer.sumAllMatchingTags(neknekLocalKernelPredicate, "DEVICE:MAX");
    platform->timer.printStatEntry("        find kernel     ", tFindKernel, nFindKernel, tFindpts);
  }

  const double tVelocity = platform->timer.query("velocitySolve", "DEVICE:MAX");
  platform->timer.printStatEntry("    velocitySolve       ",
                                 "velocitySolve",
                                 "DEVICE:MAX",
                                 tElapsedTimeSolve);
  platform->timer.printStatEntry("      rhs               ", "velocity rhs", "DEVICE:MAX", tVelocity);
  platform->timer.printStatEntry("      preconditioner    ",
                                 "velocity preconditioner",
                                 "DEVICE:MAX",
                                 tVelocity);
  platform->timer.printStatEntry("      initial guess     ", "velocity proj", "DEVICE:MAX", tVelocity);

  const double tPressure = platform->timer.query("pressureSolve", "DEVICE:MAX");
  platform->timer.printStatEntry("    pressureSolve       ",
                                 "pressureSolve",
                                 "DEVICE:MAX",
                                 tElapsedTimeSolve);
  platform->timer.printStatEntry("      rhs               ", "pressure rhs", "DEVICE:MAX", tPressure);

  const double tPressurePreco = platform->timer.query("pressure preconditioner", "DEVICE:MAX");
  platform->timer.printStatEntry("      preconditioner    ",
                                 "pressure preconditioner",
                                 "DEVICE:MAX",
                                 tPressure);

  auto tags = platform->timer.tags();
  for (int i = 15; i > 0; i--) {
    const std::string tag = "pressure preconditioner smoother N=" + std::to_string(i);
    if (std::find(tags.begin(), tags.end(), tag) == tags.end()) {
      continue;
    }
    platform->timer.printStatEntry("        pMG smoother    ", tag, "DEVICE:MAX", tPressurePreco);
  }

  platform->timer.printStatEntry("        coarse grid     ", "coarseSolve", "DEVICE:MAX", tPressurePreco);
  platform->timer.printStatEntry("      initial guess     ", "pressure proj", "DEVICE:MAX", tPressure);

  int nScalar = 0;
  platform->options.getArgs("NUMBER OF SCALARS", nScalar);

  const double tScalar = platform->timer.query("scalarSolve", "DEVICE:MAX");
  platform->timer.printStatEntry("    scalarSolve         ", "scalarSolve", "DEVICE:MAX", tElapsedTimeSolve);
  platform->timer.printStatEntry("      rhs               ", "scalar rhs", "DEVICE:MAX", tScalar);

  auto cvodeMakeQPredicate = [](const std::string &tag) {
    bool match = tag.find("cvode_t::") != std::string::npos && tag.find("makeq") != std::string::npos;
    // ensure children of the timer aren't doubly counted
    return match && tag.find("makeq::") == std::string::npos;
  };
  auto [tMakeqCvode, nMakeqCvode] = platform->timer.sumAllMatchingTags(cvodeMakeQPredicate, "DEVICE:MAX");

  auto cvodeUdfSEqnSourcePredicate = [](const std::string &tag) {
    bool match = tag.find("cvode_t::") != std::string::npos && tag.find("udfSEqnSource") != std::string::npos;
    // ensure children of the timer aren't doubly counted
    return match && tag.find("udfSEqnSource::") == std::string::npos;
  };
  auto [tSEqnSourceCvode, nSEqnSourceCvode] =
      platform->timer.sumAllMatchingTags(cvodeUdfSEqnSourcePredicate, "DEVICE:MAX");

  auto cvodeLocalPointSourcePredicate = [](const std::string &tag) {
    bool match = tag.find("cvode_t::") != std::string::npos && tag.find("pointSource") != std::string::npos;
    // ensure children of the timer aren't doubly counted
    return match;
  };
  auto [tLocalPointSource, nLocalPointSource] =
      platform->timer.sumAllMatchingTags(cvodeLocalPointSourcePredicate, "DEVICE:MAX");

  auto cvodePropertiesPredicate = [](const std::string &tag) {
    bool match =
        tag.find("cvode_t::") != std::string::npos && tag.find("evaluateProperties") != std::string::npos;
    // ensure children of the timer aren't doubly counted
    return match && tag.find("evaluateProperties::") == std::string::npos;
  };
  auto [tPropCvode, nPropCvode] = platform->timer.sumAllMatchingTags(cvodePropertiesPredicate, "DEVICE:MAX");
  platform->timer.printStatEntry("    scalarSolveCvode    ",
                                 "cvode_t::solve",
                                 "DEVICE:MAX",
                                 tElapsedTimeSolve);
  platform->timer.printStatEntry("      makeq             ", tMakeqCvode, nMakeqCvode, tScalarCvode);
  platform->timer.printStatEntry("        udfSEqnSource   ", tSEqnSourceCvode, nSEqnSourceCvode, tMakeqCvode);
  platform->timer.printStatEntry("      local pt src      ",
                                 tLocalPointSource,
                                 nLocalPointSource,
                                 tScalarCvode);
  platform->timer.printStatEntry("      udfProperties     ", tPropCvode, nPropCvode, tScalarCvode);

  auto precoTimeScalars = 0.0;
  auto precoCallsScalars = 0.0;
  for (int is = 0; is < nScalar; is++) {
    std::string sid = scalarDigitStr(is);
    precoTimeScalars += platform->timer.query("scalar" + sid + " preconditioner", "DEVICE:MAX");
    precoCallsScalars += platform->timer.count("scalar" + sid + " preconditioner");
  }
  platform->timer.set("scalar preconditioner", precoTimeScalars, precoCallsScalars);

  platform->timer.printStatEntry("      preconditioner    ", "scalar preconditioner", "DEVICE:MAX", tScalar);
  platform->timer.printStatEntry("      initial guess     ", "scalar proj", "DEVICE:MAX", tScalar);

  platform->timer.printStatEntry("    gsMPI               ", gsTime, tElapsedTimeSolve);

  platform->timer.printStatEntry("    dotp                ", "dotp", "DEVICE:MAX", tElapsedTimeSolve);

  platform->timer.printStatEntry("    dotp multi          ", "dotpMulti", "DEVICE:MAX", tElapsedTimeSolve);

  if (platform->comm.mpiRank == 0) {
    std::cout << std::endl;
  }
  platform->device.printMemoryUsage(platform->comm.mpiComm);
  if (platform->comm.mpiRank == 0) {
    std::cout << std::endl;
  }

  std::cout.unsetf(std::ios::scientific);
  std::cout.precision(outPrecisionSave);
}

void nrs_t::finalize()
{
  if (this->uSolver) {
    delete this->uSolver;
  }
  if (this->vSolver) {
    delete this->vSolver;
  }
  if (this->wSolver) {
    delete this->wSolver;
  }
  if (this->uvwSolver) {
    delete this->uvwSolver;
  }
  if (this->pSolver) {
    delete this->pSolver;
  }
  for (int is = 0; is < this->Nscalar; is++) {
    if (this->cds->solver[is]) {
      delete this->cds->solver[is];
    }
  }
  if (this->cds) {
    if (this->cds->cvode) {
      delete this->cds->cvode;
    }
  }
  if (this->meshSolver) {
    delete this->meshSolver;
  }
  checkpointWriter.reset();
}

void nrs_t::makeNLT(double time, int tstep, occa::memory &o_Usubcycling)
{
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  if (platform->options.compareArgs("VELOCITY REGULARIZATION METHOD", "HPFRT")) {
    this->filterRTKernel(mesh->Nelements,
                         this->o_filterRT,
                         this->filterS,
                         this->fieldOffset,
                         this->o_U,
                         this->o_NLT);
    double flops = 24 * mesh->Np * mesh->Nq + 3 * mesh->Np;
    flops *= static_cast<double>(mesh->Nelements);
    platform->flopCounter->add("velocityFilterRT", flops);
  }

  if (movingMesh && !this->Nsubsteps) {
    this->advectMeshVelocityKernel(mesh->Nelements,
                                   mesh->o_vgeo,
                                   mesh->o_D,
                                   this->fieldOffset,
                                   mesh->o_U,
                                   this->o_U,
                                   this->o_NLT);
    double flops = 54 * mesh->Np * mesh->Nq + 63 * mesh->Np;
    flops *= static_cast<double>(mesh->Nelements);
    platform->flopCounter->add("velocity advectMeshVelocity", flops);
  }

  if (platform->options.compareArgs("ADVECTION", "TRUE")) {
    if (this->Nsubsteps) {
      o_Usubcycling = this->advectionSubcycling(std::min(tstep, this->nEXT), time);
    } else {
      auto o_adv = platform->deviceMemoryPool.reserve<dfloat>(this->NVfields * this->fieldOffset);

      if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
        this->strongAdvectionCubatureVolumeKernel(mesh->Nelements,
                                                  mesh->o_vgeo,
                                                  mesh->o_cubDiffInterpT,
                                                  mesh->o_cubInterpT,
                                                  mesh->o_cubProjectT,
                                                  this->fieldOffset,
                                                  this->cubatureOffset,
                                                  this->o_U,
                                                  this->o_Urst,
                                                  o_adv);
      } else {
        this->strongAdvectionVolumeKernel(mesh->Nelements,
                                          mesh->o_vgeo,
                                          mesh->o_D,
                                          this->fieldOffset,
                                          this->o_U,
                                          this->o_Urst,
                                          o_adv);
      }

      platform->linAlg->axpby(this->NVfields * this->fieldOffset, -1.0, o_adv, 1.0, this->o_NLT);

      advectionFlops(this->mesh, this->NVfields);
    }
  }
}

void nrs_t::printStepInfo(double time, int tstep, bool printStepInfo, bool printVerboseInfo)
{
  cds_t *cds = this->cds;

  const double elapsedStep = platform->timer.query("elapsedStep", "DEVICE:MAX");
  const double elapsedStepSum = platform->timer.query("elapsedStepSum", "DEVICE:MAX");
  bool verboseInfo = platform->options.compareArgs("VERBOSE SOLVER INFO", "TRUE");
  const auto cfl = this->computeCFL();
  dfloat divUErrVolAvg, divUErrL2;

  if (verboseInfo) {
    computeDivUErr(this, divUErrVolAvg, divUErrL2);
  }

  auto printSolverInfo = [](elliptic *solver, const std::string &name) {
    const auto [prevProjVecs, nProjVecs] = solver->projectionCounters();
    if (nProjVecs > 0) {
      if (prevProjVecs > 0) {
        printf("%-10s: resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
               std::string("proj" + name).c_str(),
               solver->initialResidual(),
               solver->initialGuessResidual(),
               solver->initialResidual() / solver->initialGuessResidual(),
               prevProjVecs,
               nProjVecs);
      }
    }
    printf("%-10s: iter %03d  resNorm0 %.2e  resNorm %.2e\n",
           name.c_str(),
           solver->Niter(),
           solver->initialGuessResidual(),
           solver->finalResidual());
  };

  if (platform->comm.mpiRank == 0) {
    if (verboseInfo && printVerboseInfo) {
      bool cvodePrinted = false;
      for (int is = 0; is < this->Nscalar; is++) {
        if (cds->compute[is] && !cds->cvodeSolve[is]) {
          const auto sid = scalarDigitStr(is);
          printSolverInfo(cds->solver[is], "S" + sid);
        } else if (cds->cvodeSolve[is] && !cvodePrinted) {
          this->cds->cvode->printInfo(true);
          cvodePrinted = true;
        }
      }

      if (this->neknek) {
        printf("%-10s: sync %.2e  exchange %.2e\n", "neknek", this->neknek->tSync(), this->neknek->tExch());
      }

      if (this->flow) {
        printSolverInfo(this->pSolver, "P");

        if (this->uvwSolver) {
          printSolverInfo(this->uvwSolver, "UVW");

        } else {
          printSolverInfo(this->uSolver, "U");
          printSolverInfo(this->vSolver, "V");
          printSolverInfo(this->wSolver, "W");
        }
        printf("%-10s: %.2e  %.2e\n", "divUErr", divUErrVolAvg, divUErrL2);
      }

      if (this->meshSolver) {
        printSolverInfo(this->meshSolver, "MSH");
      }
    }

    if (platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")) {
      this->flowRatePrintInfo(verboseInfo && printVerboseInfo);
    }

    const auto printTimers = printStepInfo && this->timeStepConverged;

    if (printStepInfo) {
      printf("step= %d  t= %.8e  dt=%.1e  C= %.3f", tstep, time, this->dt[0], cfl);
      if (!printTimers) {
        std::cout << std::endl;
      }
    }

    if (printTimers) {
      printf("  elapsedStep= %.2es  elapsedStepSum= %.5es\n", elapsedStep, elapsedStepSum);
    }
  }

  if (this->cds) {
    if (this->cds->cvode) {
      this->cds->cvode->resetCounters();
    }
  }

  bool largeCFLCheck = (cfl > 30) && this->numberActiveFields();

  nekrsCheck(largeCFLCheck || std::isnan(cfl) || std::isinf(cfl),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "Unreasonable CFL!");
}

void nrs_t::writeCheckpoint(double t, int step, bool enforceOutXYZ, bool enforceFP64, int N_, bool uniform)
{
  if (!checkpointWriter) {
    checkpointWriter = iofldFactory::create();
  }

  const auto outXYZ =
      (enforceOutXYZ) ? true : platform->options.compareArgs("CHECKPOINT OUTPUT MESH", "TRUE");

  if (!checkpointWriter->isInitialized()) {
    auto visMesh = (cht) ? cds->mesh[0] : mesh;
    checkpointWriter->open(visMesh, iofld::mode::write, platform->options.getArgs("CASENAME"));

    if (platform->options.compareArgs("LOWMACH", "TRUE")) {
      checkpointWriter->addVariable("p0th", p0th[0]);
    }

    if (platform->options.compareArgs("VELOCITY CHECKPOINTING", "TRUE")) {
      std::vector<occa::memory> o_V;
      for (int i = 0; i < mesh->dim; i++) {
        o_V.push_back(o_U.slice(i * fieldOffset, visMesh->Nlocal));
      }
      checkpointWriter->addVariable("velocity", o_V);
    }

    if (platform->options.compareArgs("PRESSURE CHECKPOINTING", "TRUE")) {
      auto o_p = std::vector<occa::memory>{o_P.slice(0, visMesh->Nlocal)};
      checkpointWriter->addVariable("pressure", o_p);
    }

    for (int i = 0; i < Nscalar; i++) {
      if (platform->options.compareArgs("SCALAR" + scalarDigitStr(i) + " CHECKPOINTING", "TRUE")) {
        const auto temperatureExists = platform->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE");
        std::vector<occa::memory> o_Si = {cds->o_S.slice(cds->fieldOffsetScan[i], visMesh->Nlocal)};
        if (i == 0 && temperatureExists) {
          checkpointWriter->addVariable("temperature", o_Si);
        } else {
          const auto is = (temperatureExists) ? i - 1 : i;
          checkpointWriter->addVariable("scalar" + scalarDigitStr(is), o_Si);
        }
      }
    }
  }

  const auto Nfld = [&]() {
    int N;
    platform->options.getArgs("POLYNOMIAL DEGREE", N);
    return (N_) ? N_ : N;
  }();
  checkpointWriter->writeAttribute("polynomialOrder", std::to_string(Nfld));

  auto FP64 = platform->options.compareArgs("CHECKPOINT PRECISION", "FP64");
  if (enforceFP64) {
    FP64 = true;
  }
  checkpointWriter->writeAttribute("precision", (FP64) ? "64" : "32");
  checkpointWriter->writeAttribute("uniform", (uniform) ? "true" : "false");
  checkpointWriter->writeAttribute("outputMesh", (outXYZ) ? "true" : "false");

  checkpointWriter->addVariable("time", t);

  for (const auto &entry : userCheckpointFields) {
    checkpointWriter->addVariable(entry.first, entry.second);
  }

  checkpointWriter->process();
}

int nrs_t::lastStepLocalSession(double timeNew, int tstep, double elapsedTime)
{
  double endTime = -1;
  platform->options.getArgs("END TIME", endTime);

  int numSteps = -1;
  platform->options.getArgs("NUMBER TIMESTEPS", numSteps);

  int last = 0;
  if (!platform->options.getArgs("STOP AT ELAPSED TIME").empty()) {
    double maxElaspedTime;
    platform->options.getArgs("STOP AT ELAPSED TIME", maxElaspedTime);
    if (elapsedTime > 60.0 * maxElaspedTime) {
      last = 1;
    }
  } else if (endTime >= 0) {
    const double eps = 1e-10;
    last = fabs(timeNew - endTime) < eps || timeNew > endTime;
  } else {
    last = (tstep == numSteps);
  }
  return last;
}

int nrs_t::setLastStep(double timeNew, int tstep, double elapsedTime)
{
  int last = lastStepLocalSession(timeNew, tstep, elapsedTime);

  lastStep = last;
  return last;
}

void nrs_t::copyToNek(double time, int tstep, bool updateMesh)
{
  *(nekData.istep) = tstep;
  copyToNek(time, updateMesh);
}

void nrs_t::copyToNek(double time, bool updateMesh_)
{
  if (platform->comm.mpiRank == 0) {
    printf("copying solution to nek\n");
    fflush(stdout);
  }

  *(nekData.time) = time;
  *(nekData.p0th) = p0th[0];

  auto updateMesh = [&]() {
    auto mesh = (cht) ? cds->mesh[0] : this->mesh;

    auto [x, y, z] = mesh->xyzHost();
    for (int i = 0; i < mesh->Nlocal; i++) {
      nekData.xm1[i] = x[i];
      nekData.ym1[i] = y[i];
      nekData.zm1[i] = z[i];
    }
    nek::recomputeGeometry();
  };

  {
    auto U = platform->memoryPool.reserve<dfloat>(mesh->dim * fieldOffset);
    o_U.copyTo(U, U.size());
    auto vx = U.ptr<dfloat>() + 0 * fieldOffset;
    auto vy = U.ptr<dfloat>() + 1 * fieldOffset;
    auto vz = U.ptr<dfloat>() + 2 * fieldOffset;
    for (int i = 0; i < mesh->Nlocal; i++) {
      nekData.vx[i] = vx[i];
      nekData.vy[i] = vy[i];
      nekData.vz[i] = vz[i];
    }
  }

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    auto mesh = (cht) ? cds->mesh[0] : this->mesh;

    auto U = platform->memoryPool.reserve<dfloat>(mesh->dim * fieldOffset);
    mesh->o_U.copyTo(U, U.size());
    auto wx = U.ptr<dfloat>() + 0 * fieldOffset;
    auto wy = U.ptr<dfloat>() + 1 * fieldOffset;
    auto wz = U.ptr<dfloat>() + 2 * fieldOffset;
    for (int i = 0; i < mesh->Nlocal; i++) {
      nekData.wx[i] = wx[i];
      nekData.wy[i] = wy[i];
      nekData.wz[i] = wz[i];
    }
    updateMesh_ = true;
  }

  if (updateMesh_) {
    updateMesh();
  }

  {
    auto P = platform->memoryPool.reserve<dfloat>(mesh->Nlocal);
    o_P.copyTo(P, P.size());
    auto Pptr = P.ptr<dfloat>();
    for (int i = 0; i < mesh->Nlocal; i++) {
      nekData.pr[i] = Pptr[i];
    }
  }

  if (Nscalar) {
    const dlong nekFieldOffset = nekData.lelt * std::pow(nekData.nx1, nekData.ndim);
    for (int is = 0; is < Nscalar; is++) {
      auto mesh = cds->mesh[is];

      auto S = platform->memoryPool.reserve<dfloat>(mesh->Nlocal);
      cds->o_S.copyTo(S, S.size(), 0, cds->fieldOffsetScan[is]);

      auto Sptr = S.ptr<dfloat>();
      auto Ti = nekData.t + is * nekFieldOffset;
      for (int i = 0; i < mesh->Nlocal; i++) {
        Ti[i] = Sptr[i];
      }
    }
  }
}

void nrs_t::copyFromNek()
{
  double time; // dummy
  copyFromNek(time);
}

void nrs_t::copyFromNek(double &time)
{
  if (platform->comm.mpiRank == 0) {
    printf("copying solution from nek\n");
    fflush(stdout);
  }

  time = *(nekData.time);
  p0th[0] = *(nekData.p0th);

  {
    auto U = platform->memoryPool.reserve<dfloat>(mesh->dim * fieldOffset);
    auto vx = U.ptr<dfloat>() + 0 * fieldOffset;
    auto vy = U.ptr<dfloat>() + 1 * fieldOffset;
    auto vz = U.ptr<dfloat>() + 2 * fieldOffset;
    for (int i = 0; i < mesh->Nlocal; i++) {
      vx[i] = nekData.vx[i];
      vy[i] = nekData.vy[i];
      vz[i] = nekData.vz[i];
    }
    o_U.copyFrom(U, U.size());
  }

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    auto mesh = (cht) ? cds->mesh[0] : this->mesh;

    auto U = platform->memoryPool.reserve<dfloat>(mesh->dim * fieldOffset);
    auto wx = U.ptr<dfloat>() + 0 * fieldOffset;
    auto wy = U.ptr<dfloat>() + 1 * fieldOffset;
    auto wz = U.ptr<dfloat>() + 2 * fieldOffset;
    for (int i = 0; i < mesh->Nlocal; i++) {
      wx[i] = nekData.wx[i];
      wy[i] = nekData.wy[i];
      wz[i] = nekData.wz[i];
    }
    mesh->o_U.copyFrom(U, U.size());
  }

  {
    auto P = platform->memoryPool.reserve<dfloat>(o_P.size());
    auto Pptr = P.ptr<dfloat>();
    for (int i = 0; i < mesh->Nlocal; i++) {
      Pptr[i] = nekData.pr[i];
    }
    o_P.copyFrom(P, P.size());
  }

  if (Nscalar) {
    const dlong nekFieldOffset = nekData.lelt * std::pow(nekData.nx1, nekData.ndim);
    for (int is = 0; is < Nscalar; is++) {
      auto mesh = cds->mesh[is];
      auto Ti = nekData.t + is * nekFieldOffset;

      auto S = platform->memoryPool.reserve<dfloat>(mesh->Nlocal);

      auto Sptr = S.ptr<dfloat>();
      for (int i = 0; i < mesh->Nlocal; i++) {
        Sptr[i] = Ti[i];
      }
      cds->o_S.copyFrom(S, S.size(), cds->fieldOffsetScan[is], 0);
    }
  }
}

void nrs_t::getICFromNek()
{
  nek::getIC();
  copyFromNek();
}
