#include "nrs.hpp"
#include "bdry.hpp"
#include "bcMap.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "hpf.hpp"
#include "avm.hpp"
#include "re2Reader.hpp"

#include "cdsSetup.cpp"

std::vector<std::string> fieldsToSolve(setupAide &options)
{
  int Nscalar = 0;
  options.getArgs("NUMBER OF SCALARS", Nscalar);

  std::vector<std::string> fields;

  if (!options.compareArgs("MESH SOLVER", "NONE")) {
    fields.push_back("mesh");
  }

  if (!options.compareArgs("VELOCITY SOLVER", "NONE")) {
    fields.push_back("velocity");
  }

  for (int i = 0; i < Nscalar; i++) {
    const auto sid = scalarDigitStr(i);
    if (!options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
      fields.push_back("scalar" + sid);
    }
  }
  return fields;
}

void printICMinMax(nrs_t *nrs)
{
  if (platform->comm.mpiRank == 0) {
    printf("================= INITIAL CONDITION ====================\n");
  }

  {
    auto mesh = nrs->_mesh;
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
    auto mesh = nrs->_mesh;
    auto o_ux = mesh->o_U + 0 * nrs->fieldOffset * sizeof(dfloat);
    auto o_uy = mesh->o_U + 1 * nrs->fieldOffset * sizeof(dfloat);
    auto o_uz = mesh->o_U + 2 * nrs->fieldOffset * sizeof(dfloat);
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
    auto mesh = nrs->meshV;
    auto o_ux = nrs->o_U + 0 * nrs->fieldOffset * sizeof(dfloat);
    auto o_uy = nrs->o_U + 1 * nrs->fieldOffset * sizeof(dfloat);
    auto o_uz = nrs->o_U + 2 * nrs->fieldOffset * sizeof(dfloat);
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
    auto mesh = nrs->meshV;
    const auto prMin = platform->linAlg->min(mesh->Nlocal, nrs->o_P, platform->comm.mpiComm);
    const auto prMax = platform->linAlg->max(mesh->Nlocal, nrs->o_P, platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("P     min/max: %g %g\n", prMin, prMax);
    }
  }

  if (nrs->Nscalar) {
    auto cds = nrs->cds;
    if (platform->comm.mpiRank == 0) {
      printf("S     min/max:");
    }

    int cnt = 0;
    for (int is = 0; is < cds->NSfields; is++) {
      cnt++;

      mesh_t *mesh;
      (is) ? mesh = cds->meshV : mesh = cds->mesh[0]; // only first scalar can be a CHT mesh

      auto o_si = nrs->cds->o_S + nrs->cds->fieldOffsetScan[is] * sizeof(dfloat);
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

void nrsSetup(MPI_Comm comm, setupAide &options, nrs_t *nrs)
{
  platform_t *platform = platform_t::getInstance();
  device_t &device = platform->device;
  nrs->kernelInfo = new occa::properties();
  *(nrs->kernelInfo) = platform->kernelInfo;
  occa::properties &kernelInfo = *nrs->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  kernelInfo["include_paths"].asArray();

  int N, cubN;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  platform->options.getArgs("CUBATURE POLYNOMIAL DEGREE", cubN);
  platform->options.getArgs("NUMBER OF SCALARS", nrs->Nscalar);
  platform->options.getArgs("MESH DIMENSION", nrs->dim);
  platform->options.getArgs("ELEMENT TYPE", nrs->elementType);

  {
#if 1
    if (platform->device.mode() == "Serial") {
      platform->options.setArgs("ENABLE GS COMM OVERLAP", "FALSE");
    }
#endif

    if (platform->comm.mpiCommSize == 1) {
      platform->options.setArgs("ENABLE GS COMM OVERLAP", "FALSE");
    }

    if (platform->comm.mpiRank == 0 && platform->options.compareArgs("ENABLE GS COMM OVERLAP", "FALSE")) {
      std::cout << "ENABLE GS COMM OVERLAP disabled\n\n";
    }
  }

  nrs->flow = 1;
  if (platform->options.compareArgs("VELOCITY SOLVER", "NONE")) {
    nrs->flow = 0;
  }

  if (nrs->flow) {
    if (platform->options.compareArgs("VELOCITY STRESSFORMULATION", "TRUE")) {
      platform->options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    }
  }

  {
    int nelgt, nelgv;
    const std::string meshFile = options.getArgs("MESH FILE");
    re2::nelg(meshFile, nelgt, nelgv, platform->comm.mpiComm);

    nrsCheck(nelgt != nelgv && platform->options.compareArgs("MOVING MESH", "TRUE"),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "Conjugate heat transfer not supported in a moving mesh!");

    nrsCheck(nelgt != nelgv && !platform->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE"),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "Conjugate heat transfer requires a temperature field!");

    bool coupled = neknekCoupled();
    nrsCheck(nelgt != nelgv && coupled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "Conjugate heat transfer + neknek not supported!");
  }

  // init nek
  {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    std::string casename;
    platform->options.getArgs("CASENAME", casename);

    nek::setup(nrs);
    nek::setic();
    nek::userchk();
    if (platform->comm.mpiRank == 0) {
      std::cout << "\n";
    }
  }

  nrs->cht = 0;
  {
    hlong NelementsV = nekData.nelv;
    hlong NelementsT = nekData.nelt;
    MPI_Allreduce(MPI_IN_PLACE, &NelementsV, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);
    MPI_Allreduce(MPI_IN_PLACE, &NelementsT, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);
    if ((NelementsT > NelementsV) && nrs->Nscalar) {
      nrs->cht = 1;
    }

    nrsCheck(nrs->cht && (NelementsT <= NelementsV),
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "Invalid solid element partitioning");
  }

  nrs->_mesh = createMesh(comm, N, cubN, nrs->cht, kernelInfo);
  nrs->meshV = (mesh_t *)nrs->_mesh->fluid;
  mesh_t *mesh = nrs->meshV;

  // verify boundary conditions
  {
    auto fields = fieldsToSolve(options);

    for (const auto &field : fields) {
      auto msh = (nrs->cht && (field == "scalar00" || field == "mesh")) ? nrs->_mesh : mesh;
      nrsCheck(msh->Nbid != bcMap::size(field),
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "Size of %s boundaryTypeMap does not match number of boundary IDs in mesh!\n",
               field.c_str());
    }

    std::vector<mesh_t *> meshList;
    meshList.push_back(nrs->_mesh);
    if (nrs->meshV != nrs->_mesh) {
      meshList.push_back(nrs->meshV);
    }

    for (const auto &msh : meshList) {
      bcMap::checkBoundaryAlignment(msh);
      bcMap::remapUnalignedBoundaries(msh);
    }
  }

  nrs->NVfields = nrs->dim;

  platform->options.getArgs("SUBCYCLING STEPS", nrs->Nsubsteps);
  platform->options.getArgs("DT", nrs->dt[0]);

  nrs->idt = 1 / nrs->dt[0];
  nrs->g0 = 1;

  platform->options.getArgs("BDF ORDER", nrs->nBDF);
  platform->options.getArgs("EXT ORDER", nrs->nEXT);
  if (nrs->Nsubsteps) {
    nrs->nEXT = nrs->nBDF;
  }

  nrsCheck(nrs->nEXT < nrs->nBDF,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "%s\n",
           "EXT order needs to be >= BDF order!");

  nrs->coeffEXT = (dfloat *)calloc(nrs->nEXT, sizeof(dfloat));
  nrs->coeffBDF = (dfloat *)calloc(nrs->nBDF, sizeof(dfloat));

  nrs->nRK = 4;

  dfloat mue = 1;
  dfloat rho = 1;
  platform->options.getArgs("VISCOSITY", mue);
  platform->options.getArgs("DENSITY", rho);

  { // setup fieldOffset
    nrs->fieldOffset = mesh->Np * (mesh->Nelements + mesh->totalHaloPairs);
    mesh_t *meshT = nrs->_mesh;
    nrs->fieldOffset = std::max(nrs->fieldOffset, meshT->Np * (meshT->Nelements + meshT->totalHaloPairs));
    nrs->fieldOffset = alignStride<dfloat>(nrs->fieldOffset);
  }
  nrs->_mesh->fieldOffset = nrs->fieldOffset;

  { // setup cubatureOffset
    if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
      nrs->cubatureOffset = std::max(nrs->fieldOffset, mesh->Nelements * mesh->cubNp);
    } else {
      nrs->cubatureOffset = nrs->fieldOffset;
    }
    nrs->cubatureOffset = alignStride<dfloat>(nrs->cubatureOffset);
  }

  if (nrs->Nsubsteps) {
    int Sorder;
    platform->options.getArgs("SUBCYCLING TIME ORDER", Sorder);
    if (Sorder == 4 && nrs->nRK == 4) { // ERK(4,4)
      dfloat rka[4] = {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0};
      dfloat rkb[4] = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
      dfloat rkc[4] = {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0};
      nrs->coeffsfRK = (dfloat *)calloc(nrs->nRK, sizeof(dfloat));
      nrs->weightsRK = (dfloat *)calloc(nrs->nRK, sizeof(dfloat));
      nrs->nodesRK = (dfloat *)calloc(nrs->nRK, sizeof(dfloat));
      memcpy(nrs->coeffsfRK, rka, nrs->nRK * sizeof(dfloat));
      memcpy(nrs->weightsRK, rkb, nrs->nRK * sizeof(dfloat));
      memcpy(nrs->nodesRK, rkc, nrs->nRK * sizeof(dfloat));
    } else {
      nrsCheck(true, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "Unsupported subcycling scheme!");
    }
    nrs->o_coeffsfRK = device.malloc(nrs->nRK * sizeof(dfloat), nrs->coeffsfRK);
    nrs->o_weightsRK = device.malloc(nrs->nRK * sizeof(dfloat), nrs->weightsRK);
  }

  // setup mempool
  int ellipticMaxFields = 1;
  if (platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE") ||
      !platform->options.compareArgs("MESH SOLVER", "NONE")) {
    ellipticMaxFields = nrs->NVfields;
  }

  int wrkFields = 10;
  if (nrs->Nsubsteps) {
    wrkFields = 9 + 3 * nrs->NVfields;
  }
  if (options.compareArgs("MOVING MESH", "TRUE")) {
    wrkFields += nrs->NVfields;
  }

  const int mempoolNflds = std::max(wrkFields, 2 * nrs->NVfields + elliptic_t::NWorkspaceFields * ellipticMaxFields);
  platform->create_mempool(nrs->fieldOffset, mempoolNflds);

  if (options.compareArgs("MOVING MESH", "TRUE")) {
    const int nBDF = std::max(nrs->nBDF, nrs->nEXT);
    platform->o_mempool.slice0.copyFrom(mesh->o_LMM, mesh->Nlocal * sizeof(dfloat));
    mesh->o_LMM.free();
    mesh->o_LMM = platform->device.malloc(nrs->fieldOffset * nBDF, sizeof(dfloat));
    mesh->o_LMM.copyFrom(platform->o_mempool.slice0, mesh->Nlocal * sizeof(dfloat));

    platform->o_mempool.slice0.copyFrom(mesh->o_invLMM, mesh->Nlocal * sizeof(dfloat));
    mesh->o_invLMM.free();
    mesh->o_invLMM = platform->device.malloc(nrs->fieldOffset * nBDF, sizeof(dfloat));
    mesh->o_invLMM.copyFrom(platform->o_mempool.slice0, mesh->Nlocal * sizeof(dfloat));

    const int nAB = std::max(nrs->nEXT, mesh->nAB);
    mesh->U = (dfloat *)calloc(nrs->NVfields * nrs->fieldOffset * nAB, sizeof(dfloat));
    mesh->o_U = platform->device.malloc((nrs->NVfields * nAB * sizeof(dfloat)) * nrs->fieldOffset, mesh->U);
    mesh->o_Ue = platform->device.malloc((nrs->NVfields * nAB * sizeof(dfloat)) * nrs->fieldOffset);
    if (nrs->Nsubsteps) {
      mesh->o_divU = platform->device.malloc(nrs->fieldOffset * nAB, sizeof(dfloat));
    }
  }

  {
    const dlong Nstates = nrs->Nsubsteps ? std::max(nrs->nBDF, nrs->nEXT) : 1;
    bool useCVODE = platform->options.compareArgs("CVODE", "TRUE");
    if ((useCVODE || nrs->Nsubsteps) && platform->options.compareArgs("MOVING MESH", "TRUE")) {
      nrs->o_relUrst =
          platform->device.malloc((Nstates * nrs->NVfields * sizeof(dfloat)) * nrs->cubatureOffset);
    }
    if (!nrs->Nsubsteps || platform->options.compareArgs("MOVING MESH", "FALSE")) {
      nrs->o_Urst = platform->device.malloc((Nstates * nrs->NVfields * sizeof(dfloat)) * nrs->cubatureOffset);
    }
  }

  nrs->U = (dfloat *)calloc(nrs->NVfields * std::max(nrs->nBDF, nrs->nEXT) * nrs->fieldOffset, sizeof(dfloat));
  nrs->o_U = platform->device.malloc(nrs->NVfields * std::max(nrs->nBDF, nrs->nEXT) * nrs->fieldOffset * sizeof(dfloat),
                                     nrs->U);

  nrs->o_Ue = platform->device.malloc((nrs->NVfields * sizeof(dfloat)) * nrs->fieldOffset);

  nrs->P = (dfloat *)calloc(nrs->fieldOffset, sizeof(dfloat));
  nrs->o_P = platform->device.malloc(nrs->fieldOffset * sizeof(dfloat), nrs->P);

  nrs->o_BF = platform->device.malloc((nrs->NVfields * sizeof(dfloat)) * nrs->fieldOffset);
  nrs->o_FU = platform->device.malloc((nrs->NVfields * nrs->nEXT * sizeof(dfloat)) * nrs->fieldOffset);

  nrs->o_ellipticCoeff = device.malloc((2 * sizeof(dfloat)) * nrs->fieldOffset);

  int nProperties = 2;
  if (!options.compareArgs("MESH SOLVER", "NONE")) {
    nProperties = 4;
  }

  nrs->o_prop = device.malloc((nProperties * sizeof(dfloat)) * nrs->fieldOffset);
  nrs->o_mue = nrs->o_prop.slice((0 * sizeof(dfloat)) * nrs->fieldOffset);
  nrs->o_rho = nrs->o_prop.slice((1 * sizeof(dfloat)) * nrs->fieldOffset);
  if (!options.compareArgs("MESH SOLVER", "NONE")) {
    nrs->o_meshMue = nrs->o_prop.slice((2 * sizeof(dfloat)) * nrs->fieldOffset);
    nrs->o_meshRho = nrs->o_prop.slice((3 * sizeof(dfloat)) * nrs->fieldOffset);
  }

  platform->linAlg->fill(mesh->Nlocal, mue, nrs->o_mue);
  platform->linAlg->fill(mesh->Nlocal, rho, nrs->o_rho);
  if (!options.compareArgs("MESH SOLVER", "NONE")) {
    auto o_mue = nrs->o_prop + (2 * nrs->fieldOffset)*sizeof(dfloat);
    auto o_rho = nrs->o_prop + (3 * nrs->fieldOffset)*sizeof(dfloat);
    platform->linAlg->fill(mesh->Nlocal, 1.0, o_mue);
    platform->linAlg->fill(mesh->Nlocal, 0.0, o_rho);
  }

  if (platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")) {
    nrs->o_Uc = platform->device.malloc((nrs->NVfields * sizeof(dfloat)) * nrs->fieldOffset);
    nrs->o_Pc = platform->device.malloc(nrs->fieldOffset * sizeof(dfloat));
    nrs->o_prevProp = device.malloc((2 * sizeof(dfloat)) * nrs->fieldOffset);
    nrs->o_prevProp.copyFrom(nrs->o_prop, nrs->o_prevProp.size());
  }

  nrs->o_div = device.malloc(nrs->fieldOffset * sizeof(dfloat));

  nrs->o_coeffEXT = platform->device.malloc(nrs->nEXT * sizeof(dfloat), nrs->coeffEXT);
  nrs->o_coeffBDF = platform->device.malloc(nrs->nBDF * sizeof(dfloat), nrs->coeffBDF);

  nrs->gsh = oogs::setup(mesh->ogs, nrs->NVfields, nrs->fieldOffset, ogsDfloat, NULL, OOGS_AUTO);

  if (!options.compareArgs("MESH SOLVER", "NONE")) {
    mesh_t *meshT = nrs->_mesh;
    nrs->gshMesh = oogs::setup(meshT->ogs, nrs->NVfields, nrs->fieldOffset, ogsDfloat, NULL, OOGS_AUTO);
  }

  if (nrs->flow) {
    nrs->EToB = (int *)calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
    int cnt = 0;
    for (int e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        nrs->EToB[cnt] = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], "velocity");
        cnt++;
      }
    }
    nrs->o_EToB = device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int), nrs->EToB);
  }

  if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
    nrs->EToBMeshVelocity = (int *)calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
    int cnt = 0;
    for (int e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        int bc = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], "mesh");
        nrs->EToBMeshVelocity[cnt] = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], "mesh");
        cnt++;
      }
    }
    nrs->o_EToBMeshVelocity =
        device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int), nrs->EToBMeshVelocity);
  }

  if (platform->options.compareArgs("VELOCITY REGULARIZATION METHOD", "HPFRT")) {

    nrs->filterNc = -1;
    dfloat filterS;
    platform->options.getArgs("VELOCITY HPFRT STRENGTH", filterS);
    platform->options.getArgs("VELOCITY HPFRT MODES", nrs->filterNc);
    filterS = -1.0 * fabs(filterS);
    nrs->filterS = filterS;

    nrs->o_filterMT = hpfSetup(nrs->meshV, nrs->filterNc);
  }

  // build kernels
  std::string kernelName;
  const std::string suffix = "Hex3D";
  {
    const std::string section = "nrs-";
    kernelName = "nStagesSum3";
    nrs->nStagesSum3Kernel = platform->kernels.get(section + kernelName);

    kernelName = "computeFieldDotNormal";
    nrs->computeFieldDotNormalKernel = platform->kernels.get(section + kernelName);

    kernelName = "computeFaceCentroid";
    nrs->computeFaceCentroidKernel = platform->kernels.get(section + kernelName);

    {
      kernelName = "strongAdvectionVolume" + suffix;
      nrs->strongAdvectionVolumeKernel = platform->kernels.get(section + kernelName);
      kernelName = "strongAdvectionCubatureVolume" + suffix;
      nrs->strongAdvectionCubatureVolumeKernel = platform->kernels.get(section + kernelName);
    }

    kernelName = "curl" + suffix;
    nrs->curlKernel = platform->kernels.get(section + kernelName);

    kernelName = "SijOij" + suffix;
    nrs->SijOijKernel = platform->kernels.get(section + kernelName);

    kernelName = "gradientVolume" + suffix;
    nrs->gradientVolumeKernel = platform->kernels.get(section + kernelName);

    kernelName = "wGradientVolume" + suffix;
    nrs->wgradientVolumeKernel = platform->kernels.get(section + kernelName);

    {
      kernelName = "sumMakef";
      nrs->sumMakefKernel = platform->kernels.get(section + kernelName);
    }

    kernelName = "wDivergenceVolume" + suffix;
    nrs->wDivergenceVolumeKernel = platform->kernels.get(section + kernelName);
    kernelName = "divergenceVolume" + suffix;
    nrs->divergenceVolumeKernel = platform->kernels.get(section + kernelName);

    kernelName = "divergenceSurface" + suffix;
    nrs->divergenceSurfaceKernel = platform->kernels.get(section + kernelName);

    kernelName = "advectMeshVelocity" + suffix;
    nrs->advectMeshVelocityKernel = platform->kernels.get(section + kernelName);

    kernelName = "pressureRhs" + suffix;
    nrs->pressureRhsKernel = platform->kernels.get(section + kernelName);

    kernelName = "pressureStress" + suffix;
    nrs->pressureStressKernel = platform->kernels.get(section + kernelName);

    kernelName = "pressureDirichletBC" + suffix;
    nrs->pressureDirichletBCKernel = platform->kernels.get(section + kernelName);

    kernelName = "velocityRhs" + suffix;
    nrs->velocityRhsKernel = platform->kernels.get(section + kernelName);

    kernelName = "averageNormalBcType";
    nrs->averageNormalBcTypeKernel = platform->kernels.get(section + kernelName);

    kernelName = "fixZeroNormalMask";
    nrs->fixZeroNormalMaskKernel = platform->kernels.get(section + kernelName);

    kernelName = "applyZeroNormalMask";
    nrs->applyZeroNormalMaskKernel = platform->kernels.get(section + kernelName);

    kernelName = "initializeZeroNormalMask";
    nrs->initializeZeroNormalMaskKernel = platform->kernels.get(section + kernelName);

    kernelName = "velocityDirichletBC" + suffix;
    nrs->velocityDirichletBCKernel = platform->kernels.get(section + kernelName);

    kernelName = "velocityNeumannBC" + suffix;
    nrs->velocityNeumannBCKernel = platform->kernels.get(section + kernelName);

    kernelName = "UrstCubature" + suffix;
    nrs->UrstCubatureKernel = platform->kernels.get(section + kernelName);

    kernelName = "Urst" + suffix;
    nrs->UrstKernel = platform->kernels.get(section + kernelName);

    if (nrs->Nsubsteps) {
      if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
        kernelName = "subCycleStrongCubatureVolume" + suffix;
        nrs->subCycleStrongCubatureVolumeKernel = platform->kernels.get(section + kernelName);
      }
      kernelName = "subCycleStrongVolume" + suffix;
      nrs->subCycleStrongVolumeKernel = platform->kernels.get(section + kernelName);

      kernelName = "subCycleRKUpdate";
      nrs->subCycleRKUpdateKernel = platform->kernels.get(section + kernelName);
      kernelName = "subCycleRK";
      nrs->subCycleRKKernel = platform->kernels.get(section + kernelName);

      kernelName = "subCycleInitU0";
      nrs->subCycleInitU0Kernel = platform->kernels.get(section + kernelName);
    }

    kernelName = "extrapolate";
    nrs->extrapolateKernel = platform->kernels.get(section + kernelName);

    kernelName = "maskCopy";
    nrs->maskCopyKernel = platform->kernels.get(section + kernelName);

    kernelName = "maskCopy2";
    nrs->maskCopy2Kernel = platform->kernels.get(section + kernelName);

    kernelName = "mask";
    nrs->maskKernel = platform->kernels.get(section + kernelName);

    kernelName = "filterRT" + suffix;
    nrs->filterRTKernel = platform->kernels.get(section + kernelName);

    kernelName = "cfl" + suffix;
    nrs->cflKernel = platform->kernels.get(section + kernelName);

    kernelName = "pressureAddQtl";
    nrs->pressureAddQtlKernel = platform->kernels.get(section + kernelName);

    kernelName = "setEllipticCoeff";
    nrs->setEllipticCoeffKernel = platform->kernels.get(section + kernelName);
    kernelName = "setEllipticCoeffPressure";
    nrs->setEllipticCoeffPressureKernel = platform->kernels.get(section + kernelName);
  }

  if (nrs->Nscalar) {
    nrs->cds = cdsSetup(nrs, platform->options);
    if (nrs->cds->anyCvodeSolver) {
      nrs->cvode = new cvode_t(nrs);
      nrs->cds->cvode = nrs->cvode;
    }
  }

  // get IC + t0 from nek
  double startTime;
  nek::copyFromNek(startTime);
  platform->options.setArgs("START TIME", to_string_f(startTime));

  // udf setup
  if (platform->comm.mpiRank == 0) {
    printf("calling udf_setup ... ");
  }
  fflush(stdout);

  udf.setup(nrs);

  if (platform->comm.mpiRank == 0) {
    printf("done\n");
  }
  fflush(stdout);

  nrs->p0the = nrs->p0th[0];

  // in case the user modifies mesh in udf.setup
  nrs->_mesh->o_x.copyFrom(nrs->_mesh->x);
  nrs->_mesh->o_y.copyFrom(nrs->_mesh->y);
  nrs->_mesh->o_z.copyFrom(nrs->_mesh->z);
  if (nrs->meshV != nrs->_mesh) {
    nrs->meshV->update();
  }
  nrs->_mesh->update();

  // in case the user sets IC in udf.setup
  nrs->o_U.copyFrom(nrs->U);
  nrs->o_P.copyFrom(nrs->P);
  if (nrs->Nscalar) {
    nrs->cds->o_S.copyFrom(nrs->cds->S);
  }
  if (options.compareArgs("MOVING MESH", "TRUE")) {
    mesh->o_U.copyFrom(mesh->U);
  }

  // ensure both codes see the same mesh + IC
  nek::ocopyToNek(startTime, 0);

  // update props based on IC
  evaluateProperties(nrs, startTime);

  // CVODE can only be initialized once the initial condition
  // is known, however, a user may need to set function ptrs
  // on to cvode_t object.
  // Hence, the actual CVODE initialization part of cvode_t
  // is done below.
  if (nrs->cvode) {
    nrs->cvode->initialize();
  }

  if (platform->comm.mpiRank == 0) {
    std::cout << std::endl;
  }
  printMeshMetrics(nrs->_mesh);

  printICMinMax(nrs);

  // setup elliptic solver

  if (nrs->Nscalar) {
    cds_t *cds = nrs->cds;

    for (int is = 0; is < cds->NSfields; is++) {
      std::string sid = scalarDigitStr(is);

      if (!cds->compute[is]) {
        continue;
      }

      mesh_t *mesh;
      (is) ? mesh = cds->meshV : mesh = cds->mesh[0]; // only first scalar can be a CHT mesh

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

      cds->solver[is] = new elliptic_t();
      cds->solver[is]->name = "scalar" + sid;
      cds->solver[is]->Nfields = 1;
      cds->solver[is]->fieldOffset = nrs->fieldOffset;
      cds->solver[is]->mesh = mesh;

      cds->solver[is]->poisson = 0;

      cds->setEllipticCoeffKernel(mesh->Nlocal,
                                  cds->g0 * cds->idt,
                                  cds->fieldOffsetScan[is],
                                  nrs->fieldOffset,
                                  0,
                                  cds->o_diff,
                                  cds->o_rho,
                                  o_NULL,
                                  cds->o_ellipticCoeff);

      cds->solver[is]->o_lambda0 = cds->o_ellipticCoeff.slice(0 * nrs->fieldOffset * sizeof(dfloat));
      cds->solver[is]->o_lambda1 = cds->o_ellipticCoeff.slice(1 * nrs->fieldOffset * sizeof(dfloat));

      cds->solver[is]->EToB = (int *)calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
      for (dlong e = 0; e < mesh->Nelements; e++) {
        for (int f = 0; f < mesh->Nfaces; f++) {
          const int bID = mesh->EToB[f + e * mesh->Nfaces];
          cds->solver[is]->EToB[f + e * mesh->Nfaces] = bcMap::ellipticType(bID, "scalar" + sid);
        }
      }

      ellipticSolveSetup(cds->solver[is]);
    }
  }

  if (nrs->flow) {
    if (platform->comm.mpiRank == 0) {
      printf("================ ELLIPTIC SETUP VELOCITY ================\n");
    }

    nrs->uvwSolver = NULL;

    bool unalignedBoundary = bcMap::unalignedMixedBoundary("velocity");

    nrsCheck(unalignedBoundary && !options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE"),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "SHL or unaligned SYM boundaries require solver = pcg+block");

    if (platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE")) {
      nrs->uvwSolver = new elliptic_t();
    }

    for (int bID = 1; bID <= bcMap::size("velocity"); bID++) {
      std::string bcTypeText(bcMap::text(bID, "velocity"));
      if (platform->comm.mpiRank == 0 && bcTypeText.size()) {
        printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str());
      }
    }

    nrs->setEllipticCoeffKernel(mesh->Nlocal,
                                nrs->g0 * nrs->idt,
                                0 * nrs->fieldOffset,
                                nrs->fieldOffset,
                                0,
                                nrs->o_mue,
                                nrs->o_rho,
                                o_NULL,
                                nrs->o_ellipticCoeff);

    if (nrs->uvwSolver) {
      nrs->uvwSolver->name = "velocity";
      nrs->uvwSolver->stressForm = 0;
      if (options.compareArgs("VELOCITY STRESSFORMULATION", "TRUE")) {
        nrs->uvwSolver->stressForm = 1;
      }
      nrs->uvwSolver->Nfields = nrs->NVfields;
      nrs->uvwSolver->fieldOffset = nrs->fieldOffset;
      nrs->uvwSolver->mesh = mesh;
      nrs->uvwSolver->o_lambda0 = nrs->o_ellipticCoeff.slice(0 * nrs->fieldOffset * sizeof(dfloat));
      nrs->uvwSolver->o_lambda1 = nrs->o_ellipticCoeff.slice(1 * nrs->fieldOffset * sizeof(dfloat));
      nrs->uvwSolver->poisson = 0;
      nrs->uvwSolver->EToB =
          (int *)calloc(mesh->Nelements * mesh->Nfaces * nrs->uvwSolver->Nfields, sizeof(int));
      for (int fld = 0; fld < nrs->uvwSolver->Nfields; fld++) {
        std::string key;
        if (fld == 0) {
          key = "x-velocity";
        }
        if (fld == 1) {
          key = "y-velocity";
        }
        if (fld == 2) {
          key = "z-velocity";
        }
        for (dlong e = 0; e < mesh->Nelements; e++) {
          for (int f = 0; f < mesh->Nfaces; f++) {
            const int offset = fld * mesh->Nelements * mesh->Nfaces;
            const int bID = mesh->EToB[f + e * mesh->Nfaces];
            nrs->uvwSolver->EToB[f + e * mesh->Nfaces + offset] = bcMap::ellipticType(bID, key);
          }
        }
      }

      if (unalignedBoundary) {
        nrs->o_zeroNormalMaskVelocity =
            platform->device.malloc((nrs->uvwSolver->Nfields * sizeof(dfloat)) * nrs->uvwSolver->fieldOffset);
        nrs->o_EToBVVelocity = platform->device.malloc(nrs->meshV->Nlocal * sizeof(int));
        createEToBV(nrs->meshV, nrs->uvwSolver->EToB, nrs->o_EToBVVelocity);
        auto o_EToB =
            platform->device.malloc(mesh->Nelements * mesh->Nfaces * nrs->uvwSolver->Nfields * sizeof(int),
                                    nrs->uvwSolver->EToB);
        createZeroNormalMask(nrs, mesh, o_EToB, nrs->o_EToBVVelocity, nrs->o_zeroNormalMaskVelocity);

        nrs->uvwSolver->applyZeroNormalMask =
            [nrs, mesh](dlong Nelements, occa::memory &o_elementList, occa::memory &o_x) {
              applyZeroNormalMask(nrs,
                                  mesh,
                                  Nelements,
                                  o_elementList,
                                  nrs->uvwSolver->o_EToB,
                                  nrs->o_zeroNormalMaskVelocity,
                                  o_x);
            };
      }
      ellipticSolveSetup(nrs->uvwSolver);
    } else {
      nrs->uSolver = new elliptic_t();
      nrs->uSolver->name = "velocity";
      nrs->uSolver->Nfields = 1;
      nrs->uSolver->fieldOffset = nrs->fieldOffset;
      nrs->uSolver->mesh = mesh;
      nrs->uSolver->o_lambda0 = nrs->o_ellipticCoeff.slice(0 * nrs->fieldOffset * sizeof(dfloat));
      nrs->uSolver->o_lambda1 = nrs->o_ellipticCoeff.slice(1 * nrs->fieldOffset * sizeof(dfloat));
      nrs->uSolver->poisson = 0;
      nrs->uSolver->EToB = (int *)calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
      for (dlong e = 0; e < mesh->Nelements; e++) {
        for (int f = 0; f < mesh->Nfaces; f++) {
          const int bID = mesh->EToB[f + e * mesh->Nfaces];
          nrs->uSolver->EToB[f + e * mesh->Nfaces] = bcMap::ellipticType(bID, "x-velocity");
        }
      }

      ellipticSolveSetup(nrs->uSolver);

      nrs->vSolver = new elliptic_t();
      nrs->vSolver->name = "velocity";
      nrs->vSolver->Nfields = 1;
      nrs->vSolver->fieldOffset = nrs->fieldOffset;
      nrs->vSolver->mesh = mesh;
      nrs->vSolver->o_lambda0 = nrs->o_ellipticCoeff.slice(0 * nrs->fieldOffset * sizeof(dfloat));
      nrs->vSolver->o_lambda1 = nrs->o_ellipticCoeff.slice(1 * nrs->fieldOffset * sizeof(dfloat));
      nrs->vSolver->poisson = 0;
      nrs->vSolver->EToB = (int *)calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
      for (dlong e = 0; e < mesh->Nelements; e++) {
        for (int f = 0; f < mesh->Nfaces; f++) {
          const int bID = mesh->EToB[f + e * mesh->Nfaces];
          nrs->vSolver->EToB[f + e * mesh->Nfaces] = bcMap::ellipticType(bID, "y-velocity");
        }
      }

      ellipticSolveSetup(nrs->vSolver);

      nrs->wSolver = new elliptic_t();
      nrs->wSolver->name = "velocity";
      nrs->wSolver->Nfields = 1;
      nrs->wSolver->fieldOffset = nrs->fieldOffset;
      nrs->wSolver->mesh = mesh;
      nrs->wSolver->o_lambda0 = nrs->o_ellipticCoeff.slice(0 * nrs->fieldOffset * sizeof(dfloat));
      nrs->wSolver->o_lambda1 = nrs->o_ellipticCoeff.slice(1 * nrs->fieldOffset * sizeof(dfloat));
      nrs->wSolver->poisson = 0;
      nrs->wSolver->EToB = (int *)calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
      for (dlong e = 0; e < mesh->Nelements; e++) {
        for (int f = 0; f < mesh->Nfaces; f++) {
          const int bID = mesh->EToB[f + e * mesh->Nfaces];
          nrs->wSolver->EToB[f + e * mesh->Nfaces] = bcMap::ellipticType(bID, "z-velocity");
        }
      }

      ellipticSolveSetup(nrs->wSolver);
    }
  } // flow

  if (nrs->flow) {
    if (platform->comm.mpiRank == 0) {
      printf("================ ELLIPTIC SETUP PRESSURE ================\n");
    }

    nrs->pSolver = new elliptic_t();
    nrs->pSolver->name = "pressure";
    nrs->pSolver->Nfields = 1;
    nrs->pSolver->fieldOffset = nrs->fieldOffset;
    nrs->pSolver->mesh = mesh;

    nrs->pSolver->poisson = 1;

    // lambda0 = 1/rho  lambda1 = 0
    nrs->setEllipticCoeffPressureKernel(mesh->Nlocal, nrs->fieldOffset, nrs->o_rho, nrs->o_ellipticCoeff);

    nrs->pSolver->o_lambda0 = nrs->o_ellipticCoeff.slice(0 * nrs->fieldOffset * sizeof(dfloat));
    nrs->pSolver->o_lambda1 = nrs->o_ellipticCoeff.slice(1 * nrs->fieldOffset * sizeof(dfloat));

    nrs->pSolver->EToB = (int *)calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
    for (dlong e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        const int bID = mesh->EToB[f + e * mesh->Nfaces];
        nrs->pSolver->EToB[f + e * mesh->Nfaces] = bcMap::ellipticType(bID, "pressure");
      }
    }

    ellipticSolveSetup(nrs->pSolver);

  } // flow

  if (!options.compareArgs("MESH SOLVER", "NONE")) {
    mesh_t *mesh = nrs->_mesh;

    if (platform->comm.mpiRank == 0) {
      printf("================ ELLIPTIC SETUP MESH ================\n");
    }

    const int nbrBIDs = bcMap::size("mesh");
    for (int bID = 1; bID <= nbrBIDs; bID++) {
      std::string bcTypeText(bcMap::text(bID, "mesh"));
      if (platform->comm.mpiRank == 0 && bcTypeText.size()) {
        printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str());
      }
    }

    nrs->setEllipticCoeffKernel(mesh->Nlocal,
                                1.0,
                                0 * nrs->fieldOffset,
                                nrs->fieldOffset,
                                0,
                                nrs->o_meshMue,
                                nrs->o_meshRho,
                                o_NULL,
                                nrs->o_ellipticCoeff);

    nrs->meshSolver = new elliptic_t();
    nrs->meshSolver->name = "mesh";
    nrs->meshSolver->stressForm = 0;
    if (options.compareArgs("MESH STRESSFORMULATION", "TRUE")) {
      nrs->meshSolver->stressForm = 1;
    }
    nrs->meshSolver->Nfields = nrs->NVfields;
    nrs->meshSolver->fieldOffset = nrs->fieldOffset;
    nrs->meshSolver->mesh = mesh;
    nrs->meshSolver->o_lambda0 = nrs->o_ellipticCoeff.slice(0 * nrs->fieldOffset * sizeof(dfloat));
    nrs->meshSolver->o_lambda1 = nrs->o_ellipticCoeff.slice(1 * nrs->fieldOffset * sizeof(dfloat));
    nrs->meshSolver->poisson = 1;

    nrs->meshSolver->EToB =
        (int *)calloc(mesh->Nelements * mesh->Nfaces * nrs->meshSolver->Nfields, sizeof(int));
    for (int fld = 0; fld < nrs->meshSolver->Nfields; fld++) {
      std::string key;
      if (fld == 0) {
        key = "x-mesh";
      }
      if (fld == 1) {
        key = "y-mesh";
      }
      if (fld == 2) {
        key = "z-mesh";
      }
      for (dlong e = 0; e < mesh->Nelements; e++) {
        for (int f = 0; f < mesh->Nfaces; f++) {
          const int offset = fld * mesh->Nelements * mesh->Nfaces;
          const int bID = mesh->EToB[f + e * mesh->Nfaces];
          nrs->meshSolver->EToB[f + e * mesh->Nfaces + offset] = bcMap::ellipticType(bID, key);
        }
      }
    }

    bool unalignedBoundary = bcMap::unalignedMixedBoundary("mesh");
    if (unalignedBoundary) {
      nrs->o_zeroNormalMaskMeshVelocity =
          platform->device.malloc((nrs->meshSolver->Nfields * sizeof(dfloat)) * nrs->meshSolver->fieldOffset);
      nrs->o_EToBVMeshVelocity = platform->device.malloc(mesh->Nlocal * sizeof(int));
      auto o_EToB =
          platform->device.malloc(mesh->Nelements * mesh->Nfaces * nrs->meshSolver->Nfields * sizeof(int),
                                  nrs->meshSolver->EToB);
      createEToBV(mesh, nrs->meshSolver->EToB, nrs->o_EToBVMeshVelocity);
      createZeroNormalMask(nrs, mesh, o_EToB, nrs->o_EToBVMeshVelocity, nrs->o_zeroNormalMaskMeshVelocity);
      nrs->meshSolver->applyZeroNormalMask =
          [nrs, mesh](dlong Nelements, occa::memory &o_elementList, occa::memory &o_x) {
            applyZeroNormalMask(nrs,
                                mesh,
                                Nelements,
                                o_elementList,
                                nrs->meshSolver->o_EToB,
                                nrs->o_zeroNormalMaskMeshVelocity,
                                o_x);
          };
    }
    ellipticSolveSetup(nrs->meshSolver);
  }
}
