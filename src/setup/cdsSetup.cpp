cds_t *cdsSetup(nrs_t *nrs, setupAide options)
{
  const std::string section = "cds-";
  cds_t *cds = new cds_t();
  platform_t *platform = platform_t::getInstance();
  device_t &device = platform->device;

  cds->mesh[0] = nrs->_mesh;
  mesh_t *mesh = cds->mesh[0];
  cds->meshV = nrs->_mesh->fluid;
  cds->elementType = nrs->elementType;
  cds->dim = nrs->dim;
  cds->NVfields = nrs->NVfields;
  cds->NSfields = nrs->Nscalar;

  cds->coeffEXT = nrs->coeffEXT;
  cds->coeffBDF = nrs->coeffBDF;
  cds->coeffSubEXT = nrs->coeffSubEXT;
  cds->nBDF = nrs->nBDF;
  cds->nEXT = nrs->nEXT;
  cds->o_coeffEXT = nrs->o_coeffEXT;
  cds->o_coeffBDF = nrs->o_coeffBDF;
  cds->o_coeffSubEXT = nrs->o_coeffSubEXT;

  cds->o_usrwrk = &(nrs->o_usrwrk);

  cds->vFieldOffset = nrs->fieldOffset;
  cds->vCubatureOffset = nrs->cubatureOffset;
  cds->fieldOffset[0] = nrs->fieldOffset;
  cds->fieldOffsetScan[0] = 0;
  dlong sum = cds->fieldOffset[0];
  for (int s = 1; s < cds->NSfields; ++s) {
    cds->fieldOffset[s] = cds->fieldOffset[0];
    cds->fieldOffsetScan[s] = sum;
    sum += cds->fieldOffset[s];
    cds->mesh[s] = cds->mesh[0];
  }
  cds->fieldOffsetSum = sum;

  cds->gsh = nrs->gsh;
  cds->gshT = (nrs->cht) ? oogs::setup(mesh->ogs, 1, nrs->fieldOffset, ogsDfloat, NULL, OOGS_AUTO) : cds->gsh;

  cds->U = nrs->U;
  cds->S = (dfloat *)calloc(std::max(cds->nBDF, cds->nEXT) * cds->fieldOffsetSum, sizeof(dfloat));
  cds->BF = (dfloat *)calloc(cds->fieldOffsetSum, sizeof(dfloat));
  cds->FS = (dfloat *)calloc(cds->nEXT * cds->fieldOffsetSum, sizeof(dfloat));

  cds->Nsubsteps = nrs->Nsubsteps;
  if (cds->Nsubsteps) {
    cds->nRK = nrs->nRK;
    cds->coeffsfRK = nrs->coeffsfRK;
    cds->weightsRK = nrs->weightsRK;
    cds->nodesRK = nrs->nodesRK;
    cds->o_coeffsfRK = nrs->o_coeffsfRK;
    cds->o_weightsRK = nrs->o_weightsRK;
  }

  cds->dt = nrs->dt;

  cds->prop = (dfloat *)calloc(2 * cds->fieldOffsetSum, sizeof(dfloat));

  for (int is = 0; is < cds->NSfields; is++) {
    const int scalarWidth = getDigitsRepresentation(NSCALAR_MAX - 1);
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(scalarWidth) << is;
    std::string sid = ss.str();

    if (options.compareArgs("SCALAR" + sid + " SOLVER", "NONE"))
      continue;

    dfloat diff = 1;
    dfloat rho = 1;
    options.getArgs("SCALAR" + sid + " DIFFUSIVITY", diff);
    options.getArgs("SCALAR" + sid + " DENSITY", rho);

    const dlong off = cds->fieldOffsetSum;
    for (int e = 0; e < mesh->Nelements; e++)
      for (int n = 0; n < mesh->Np; n++) {
        cds->prop[0 * off + cds->fieldOffsetScan[is] + e * mesh->Np + n] = diff;
        cds->prop[1 * off + cds->fieldOffsetScan[is] + e * mesh->Np + n] = rho;
      }
  }

  cds->o_prop = device.malloc(2 * cds->fieldOffsetSum * sizeof(dfloat), cds->prop);
  cds->o_diff = cds->o_prop.slice(0 * cds->fieldOffsetSum * sizeof(dfloat));
  cds->o_rho = cds->o_prop.slice(1 * cds->fieldOffsetSum * sizeof(dfloat));

  cds->o_ellipticCoeff = nrs->o_ellipticCoeff;

  cds->o_U = nrs->o_U;
  cds->o_Ue = nrs->o_Ue;
  cds->o_S =
      platform->device.malloc(std::max(cds->nBDF, cds->nEXT) * cds->fieldOffsetSum * sizeof(dfloat), cds->S);
  cds->o_Se = platform->device.malloc(cds->fieldOffsetSum, sizeof(dfloat));
  cds->o_BF = platform->device.malloc(cds->fieldOffsetSum * sizeof(dfloat), cds->BF);
  cds->o_FS = platform->device.malloc(cds->nEXT * cds->fieldOffsetSum * sizeof(dfloat), cds->FS);

  cds->o_relUrst = nrs->o_relUrst;
  cds->o_Urst = nrs->o_Urst;

  for (int is = 0; is < cds->NSfields; is++) {
    const int scalarWidth = getDigitsRepresentation(NSCALAR_MAX - 1);
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(scalarWidth) << is;
    std::string sid = ss.str();

    cds->compute[is] = 1;
    if (options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
      cds->compute[is] = 0;
      continue;
    }

    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh[0]; // only first scalar can be a CHT mesh

    cds->options[is] = options;

    cds->options[is].setArgs("REGULARIZATION RAMP CONSTANT",
                             options.getArgs("SCALAR" + sid + " REGULARIZATION RAMP CONSTANT"));
    cds->options[is].setArgs("REGULARIZATION AVM C0",
                             options.getArgs("SCALAR" + sid + " REGULARIZATION AVM C0"));
    cds->options[is].setArgs("REGULARIZATION METHOD",
                             options.getArgs("SCALAR" + sid + " REGULARIZATION METHOD"));
    cds->options[is].setArgs("REGULARIZATION VISMAX COEFF",
                             options.getArgs("SCALAR" + sid + " REGULARIZATION VISMAX COEFF"));
    cds->options[is].setArgs("REGULARIZATION SCALING COEFF",
                             options.getArgs("SCALAR" + sid + " REGULARIZATION SCALING COEFF"));
    cds->options[is].setArgs("HPFRT STRENGTH", options.getArgs("SCALAR" + sid + " HPFRT STRENGTH"));
    cds->options[is].setArgs("HPFRT MODES", options.getArgs("SCALAR" + sid + " HPFRT MODES"));
    cds->options[is].setArgs("KRYLOV SOLVER", options.getArgs("SCALAR" + sid + " KRYLOV SOLVER"));
    cds->options[is].setArgs("PGMRES RESTART", options.getArgs("SCALAR" + sid + " PGMRES RESTART"));
    cds->options[is].setArgs("DISCRETIZATION", options.getArgs("SCALAR DISCRETIZATION"));
    cds->options[is].setArgs("BASIS", options.getArgs("SCALAR BASIS"));
    cds->options[is].setArgs("PRECONDITIONER", options.getArgs("SCALAR" + sid + " PRECONDITIONER"));
    cds->options[is].setArgs("SOLVER TOLERANCE", options.getArgs("SCALAR" + sid + " SOLVER TOLERANCE"));
    cds->options[is].setArgs("LINEAR SOLVER STOPPING CRITERION",
                             options.getArgs("SCALAR" + sid + " LINEAR SOLVER STOPPING CRITERION"));
    cds->options[is].setArgs("INITIAL GUESS", options.getArgs("SCALAR" + sid + " INITIAL GUESS"));
    cds->options[is].setArgs("RESIDUAL PROJECTION VECTORS",
                             options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION VECTORS"));
    cds->options[is].setArgs("RESIDUAL PROJECTION START",
                             options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION START"));
    cds->options[is].setArgs("MAXIMUM ITERATIONS", options.getArgs("SCALAR" + sid + " MAXIMUM ITERATIONS"));

    dfloat largeNumber = 1 << 20;
    cds->EToB[is] = (int *)calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
    int *EToB = cds->EToB[is];
    int cnt = 0;
    for (int e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        EToB[cnt] = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], "scalar" + sid);
        cnt++;
      }
    }
    cds->o_EToB[is] = device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int), EToB);
  }

  bool scalarFilteringEnabled = false;
  bool avmEnabled = false;
  {
    for (int is = 0; is < cds->NSfields; is++) {
      if (!cds->options[is].compareArgs("REGULARIZATION METHOD", "NONE"))
        scalarFilteringEnabled = true;
      if (cds->options[is].compareArgs("REGULARIZATION METHOD", "HPF_RESIDUAL"))
        avmEnabled = true;
      if (cds->options[is].compareArgs("REGULARIZATION METHOD", "HIGHEST_MODAL_DECAY"))
        avmEnabled = true;
    }
  }

  if (scalarFilteringEnabled) {
    const dlong Nmodes = cds->mesh[0]->N + 1;
    cds->o_filterMT = platform->device.malloc(cds->NSfields * Nmodes * Nmodes, sizeof(dfloat));
    for (int is = 0; is < cds->NSfields; is++) {
      if (cds->options[is].compareArgs("REGULARIZATION METHOD", "NONE"))
        continue;
      if (!cds->compute[is])
        continue;
      int filterNc = -1;
      cds->options[is].getArgs("HPFRT MODES", filterNc);
      dfloat filterS;
      cds->options[is].getArgs("HPFRT STRENGTH", filterS);
      filterS = -1.0 * fabs(filterS);
      cds->filterS[is] = filterS;

      dfloat *A = filterSetup(cds->mesh[is], filterNc);

      const dlong Nmodes = cds->mesh[is]->N + 1;
      cds->o_filterMT.copyFrom(A, Nmodes * Nmodes * sizeof(dfloat), is * Nmodes * Nmodes * sizeof(dfloat));

      free(A);
    }
  }

  if (avmEnabled)
    avm::setup(cds);

  std::string kernelName;
  const std::string suffix = "Hex3D";
  {
    kernelName = "strongAdvectionVolume" + suffix;
    cds->strongAdvectionVolumeKernel = platform->kernels.get(section + kernelName);

    kernelName = "strongAdvectionCubatureVolume" + suffix;
    cds->strongAdvectionCubatureVolumeKernel = platform->kernels.get(section + kernelName);

    kernelName = "advectMeshVelocity" + suffix;
    cds->advectMeshVelocityKernel = platform->kernels.get(section + kernelName);

    kernelName = "maskCopy";
    cds->maskCopyKernel = platform->kernels.get(section + kernelName);

    kernelName = "sumMakef";
    cds->sumMakefKernel = platform->kernels.get(section + kernelName);

    kernelName = "helmholtzBC" + suffix;
    cds->helmholtzRhsBCKernel = platform->kernels.get(section + kernelName);
    kernelName = "dirichletBC";
    cds->dirichletBCKernel = platform->kernels.get(section + kernelName);

    kernelName = "setEllipticCoeff";
    cds->setEllipticCoeffKernel = platform->kernels.get(section + kernelName);

    kernelName = "filterRT" + suffix;
    cds->filterRTKernel = platform->kernels.get(section + kernelName);

    kernelName = "nStagesSum3";
    cds->nStagesSum3Kernel = platform->kernels.get(section + kernelName);

    if (cds->Nsubsteps) {
      if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
        kernelName = "subCycleStrongCubatureVolume" + suffix;
        cds->subCycleStrongCubatureVolumeKernel = platform->kernels.get(section + kernelName);
      }
      kernelName = "subCycleStrongVolume" + suffix;
      cds->subCycleStrongVolumeKernel = platform->kernels.get(section + kernelName);

      kernelName = "subCycleRKUpdate";
      cds->subCycleRKUpdateKernel = platform->kernels.get(section + kernelName);
      kernelName = "subCycleRK";
      cds->subCycleRKKernel = platform->kernels.get(section + kernelName);

      kernelName = "subCycleInitU0";
      cds->subCycleInitU0Kernel = platform->kernels.get(section + kernelName);
    }
  }

  return cds;
}

