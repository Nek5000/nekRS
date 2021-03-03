#include "nrs.hpp"
#include "meshSetup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "filter.hpp"
#include "bcMap.hpp"
#include <vector>
#include <map>

static cds_t* cdsSetup(ins_t* ins, mesh_t* mesh, setupAide options, occa::properties &kernelInfoH);

void nrsSetup(MPI_Comm comm, setupAide &options, nrs_t *nrs)
{
  platform_t* platform = platform_t::getInstance();
  device_t& device = platform->device;
  nrs->options = options;
  nrs->kernelInfo = new occa::properties();
  *(nrs->kernelInfo) = platform->kernelInfo;
  occa::properties& kernelInfo = *nrs->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  kernelInfo["include_paths"].asArray();

  int N, cubN;
  int buildOnly = 0;
  string install_dir;
  if(nrs->options.compareArgs("BUILD ONLY", "TRUE")) buildOnly = 1;
  nrs->options.getArgs("POLYNOMIAL DEGREE", N);
  nrs->options.getArgs("CUBATURE POLYNOMIAL DEGREE", cubN);
  nrs->options.getArgs("NUMBER OF SCALARS", nrs->Nscalar);
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  nrs->options.getArgs("MESH DIMENSION", nrs->dim);
  nrs->options.getArgs("ELEMENT TYPE", nrs->elementType);

  nrs->flow = 1;
  if(nrs->options.compareArgs("VELOCITY", "FALSE")) nrs->flow = 0;
  if(nrs->options.compareArgs("VELOCITY SOLVER", "NONE")) nrs->flow = 0;

  if(nrs->flow) {
    if(nrs->options.compareArgs("STRESSFORMULATION", "TRUE"))
       nrs->options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
  }


  // jit compile + init nek
  {  
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    string casename;
    nrs->options.getArgs("CASENAME", casename);

    int err = 0;
    int npTarget = size;
    if (buildOnly) nrs->options.getArgs("NP TARGET", npTarget);
    if (rank == 0) err = buildNekInterface(casename.c_str(), mymax(5, nrs->Nscalar), N, npTarget, nrs->options);
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, comm);
    if (err) ABORT(EXIT_FAILURE);; 

    if (!buildOnly) {
      nek::setup(comm, nrs->options, nrs);
      nek::setic();
      nek::userchk();
    }
  }

  nrs->cht = 0;
  if (nekData.nelv != nekData.nelt && nrs->Nscalar) nrs->cht = 1;

  // create mesh
  if (buildOnly) {
    mesh_t* mesh = new mesh_t();
    mesh->linAlg = nrs->linAlg;
    createMeshDummy(mesh, comm, N, cubN, nrs->options,  kernelInfo);
    nrs->meshT = mesh;
    nrs->meshV= mesh;
  } else {
    mesh_t* mesh = new mesh_t();
    mesh->linAlg = nrs->linAlg;
    createMesh(mesh, comm, N, cubN, nrs->cht, nrs->options,  kernelInfo);
    nrs->meshT = mesh;
    nrs->meshV= mesh;
    if (nrs->cht) {
      mesh_t* meshV = new mesh_t();
      meshV->linAlg = nrs->linAlg;
      createMeshV(meshV, comm, N, cubN, nrs->meshT, nrs->options, kernelInfo);
      nrs->meshV= meshV;
    }
  }
  mesh_t* mesh = nrs->meshV;

  if (nrs->cht && !nrs->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE")) {
    if (platform->comm.mpiRank == 0) cout << "Conjugate heat transfer requires solving for temperature!\n"; 
    ABORT(EXIT_FAILURE);;
  } 

  { 
    dlong retVal; 
    MPI_Allreduce(&mesh->NinternalElements,&retVal,1,MPI_DLONG,MPI_MIN,platform->comm.mpiComm);
    if(platform->comm.mpiRank == 0) printf("min NinternalElements: %d (ratio: %4.2f)\n", retVal, (double)retVal/mesh->Nelements);
  }

  occa::properties kernelInfoV  = kernelInfo;
  occa::properties kernelInfoP  = kernelInfo;
  occa::properties kernelInfoS  = kernelInfo;

  nrs->NVfields = 3;
  nrs->NTfields = nrs->NVfields + 1;   // Total Velocity + Pressure

  nrs->nRK = 4;


  mesh->Nfields = 1;

  nrs->coeffEXT = (dfloat*) calloc(3, sizeof(dfloat));
  nrs->coeffBDF = (dfloat*) calloc(3, sizeof(dfloat));
  nrs->coeffSubEXT = (dfloat*) calloc(3, sizeof(dfloat));
  nrs->coeffRK = (dfloat*) calloc(3, sizeof(dfloat));

  if (nrs->options.compareArgs("TIME INTEGRATOR", "TOMBO1")) {
    nrs->nEXT = 1;
    nrs->nBDF = 1;
  } else if (nrs->options.compareArgs("TIME INTEGRATOR", "TOMBO2")) {
    nrs->nEXT = 2;
    nrs->nBDF = 2;
  } else if (nrs->options.compareArgs("TIME INTEGRATOR", "TOMBO3")) {
    nrs->nEXT = 3;
    nrs->nBDF = 3;
  }

  dfloat mue = 1;
  dfloat rho = 1;
  nrs->options.getArgs("VISCOSITY", mue);
  nrs->options.getArgs("DENSITY", rho);

  nrs->options.getArgs("SUBCYCLING STEPS",nrs->Nsubsteps);
  nrs->options.getArgs("DT", nrs->dt[0]);

  const dlong Nlocal = mesh->Np * mesh->Nelements;
  const dlong Ntotal = mesh->Np * (mesh->Nelements + mesh->totalHaloPairs);

  mesh->Nlocal = Nlocal;
  nrs->fieldOffset = Ntotal;

  // ensure that offset is large enough for v and t mesh and is properly aligned
  {
    const dlong NtotalT = nrs->meshT->Np * (nrs->meshT->Nelements + nrs->meshT->totalHaloPairs);
    nrs->fieldOffset = mymax(Ntotal, NtotalT);

    int PAGESIZE = 4096; // default is 4kB
    char* tmp;
    tmp = getenv("NEKRS_PAGE_SIZE");
    if (tmp != NULL) PAGESIZE = std::stoi(tmp);
    const int pageW = PAGESIZE / sizeof(dfloat);
    if (nrs->fieldOffset % pageW) nrs->fieldOffset = (nrs->fieldOffset / pageW + 1) * pageW;
  }
  nrs->meshT->fieldOffset = nrs->fieldOffset;
  nrs->meshV->fieldOffset = nrs->fieldOffset;

  if(nrs->Nsubsteps) {
    int Sorder;
    nrs->options.getArgs("SUBCYCLING TIME ORDER", Sorder);
    if(Sorder == 4 && nrs->nRK == 4) { // ERK(4,4)
      dfloat rka[4] = {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0};
      dfloat rkb[4] = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
      dfloat rkc[4] = {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0};
      nrs->Srka = (dfloat*) calloc(nrs->nRK, sizeof(dfloat));
      nrs->Srkb = (dfloat*) calloc(nrs->nRK, sizeof(dfloat));
      nrs->Srkc = (dfloat*) calloc(nrs->nRK, sizeof(dfloat));
      memcpy(nrs->Srka, rka, nrs->nRK * sizeof(dfloat));
      memcpy(nrs->Srkb, rkb, nrs->nRK * sizeof(dfloat));
      memcpy(nrs->Srkc, rkc, nrs->nRK * sizeof(dfloat));
    }else{
      if(platform->comm.mpiRank == 0) cout << "Unsupported subcycling scheme!\n";
      ABORT(1);
    }
    nrs->o_Srka = device.malloc(nrs->nRK * sizeof(dfloat), nrs->Srka);
    nrs->o_Srkb = device.malloc(nrs->nRK * sizeof(dfloat), nrs->Srkb);
  }

  // setup scratch space
  const int wrkNflds = 6;
  const int ellipticWrkNflds = 15;
  nrs->ellipticWrkOffset = wrkNflds * nrs->fieldOffset;

  const int scratchNflds = wrkNflds + ellipticWrkNflds;
  platform->create_mempool(nrs->fieldOffset, scratchNflds);

  platform->o_slice0  = platform->o_slice0;
  platform->o_slice1  = platform->o_slice1;
  platform->o_slice2  = platform->o_slice2;
  platform->o_slice3  = platform->o_slice3;
  platform->o_slice4  = platform->o_slice4;
  platform->o_slice5  = platform->o_slice5;
  platform->o_slice6  = platform->o_slice6;
  platform->o_slice7  = platform->o_slice7;
  platform->o_slice9  = platform->o_slice9;
  platform->o_slice12 = platform->o_slice12;
  platform->o_slice15 = platform->o_slice15;
  if(options.compareArgs("MOVING MESH", "TRUE")){
    // realloc o_LMM, o_invLMMM to be large enough
    const int nBDF = nrs->nBDF;
    {
      platform->o_mempool.copyFrom(mesh->o_LMM, mesh->Nlocal * sizeof(dfloat));
      mesh->o_LMM = platform->device.malloc(nrs->fieldOffset * nBDF * sizeof(dfloat));
      mesh->o_LMM.copyFrom(platform->o_mempool, mesh->Nlocal * sizeof(dfloat));
      platform->o_mempool.copyFrom(mesh->o_invLMM, mesh->Nlocal * sizeof(dfloat));
      mesh->o_invLMM = platform->device.malloc(nrs->fieldOffset * nBDF * sizeof(dfloat));
      mesh->o_invLMM.copyFrom(platform->o_mempool, mesh->Nlocal * sizeof(dfloat));
    }

    mesh->o_BdivW = platform->device.malloc(nrs->fieldOffset * nBDF * sizeof(dfloat), platform->mempool);

    const int nAB = mesh->nAB;
    mesh->U = (dfloat*) calloc(nrs->NVfields * nrs->fieldOffset * nAB, sizeof(dfloat));
    mesh->o_U = platform->device.malloc(nrs->NVfields * nrs->fieldOffset * nAB * sizeof(dfloat), mesh->U);
  }


  nrs->U  = (dfloat*) calloc(nrs->NVfields * nrs->nBDF * nrs->fieldOffset,sizeof(dfloat));
  nrs->Ue = (dfloat*) calloc(nrs->NVfields * nrs->fieldOffset,sizeof(dfloat));
  nrs->P  = (dfloat*) calloc(nrs->fieldOffset,sizeof(dfloat));
  nrs->BF = (dfloat*) calloc(nrs->NVfields * nrs->fieldOffset,sizeof(dfloat));
  nrs->FU = (dfloat*) calloc(nrs->NVfields * nrs->nEXT * nrs->fieldOffset,sizeof(dfloat));

  nrs->o_U  = platform->device.malloc(nrs->NVfields * nrs->nBDF * nrs->fieldOffset * sizeof(dfloat), nrs->U);
  nrs->o_Ue = platform->device.malloc(nrs->NVfields * nrs->fieldOffset * sizeof(dfloat), nrs->Ue);
  nrs->o_P  = platform->device.malloc(nrs->fieldOffset * sizeof(dfloat), nrs->P);
  nrs->o_BF = platform->device.malloc(nrs->NVfields * nrs->fieldOffset * sizeof(dfloat), nrs->BF);
  nrs->o_FU = platform->device.malloc(nrs->NVfields * nrs->nEXT * nrs->fieldOffset * sizeof(dfloat), nrs->FU);

  nrs->var_coeff = 1; // use always var coeff elliptic
  nrs->ellipticCoeff = (dfloat*) calloc(2 * nrs->fieldOffset,sizeof(dfloat));
  nrs->o_ellipticCoeff = device.malloc(2 * nrs->fieldOffset * sizeof(dfloat),
                                             nrs->ellipticCoeff);

  nrs->prop =  (dfloat*) calloc(2 * nrs->fieldOffset,sizeof(dfloat));
  for (int e = 0; e < mesh->Nelements; e++)
    for (int n = 0; n < mesh->Np; n++) {
      nrs->prop[0 * nrs->fieldOffset + e * mesh->Np + n] = mue;
      nrs->prop[1 * nrs->fieldOffset + e * mesh->Np + n] = rho;
    }
  nrs->o_prop = device.malloc(2 * nrs->fieldOffset * sizeof(dfloat), nrs->prop);
  nrs->o_mue = nrs->o_prop.slice(0 * nrs->fieldOffset * sizeof(dfloat));
  nrs->o_rho = nrs->o_prop.slice(1 * nrs->fieldOffset * sizeof(dfloat));

  nrs->div   = (dfloat*) calloc(nrs->fieldOffset,sizeof(dfloat));
  nrs->o_div = device.malloc(nrs->fieldOffset * sizeof(dfloat), nrs->div);

  dfloat rkC[4]  = {1.0, 0.0, -1.0, -2.0};
  nrs->o_coeffEXT = platform->device.malloc(3 * sizeof(dfloat), nrs->coeffEXT);
  nrs->o_coeffBDF = platform->device.malloc(3 * sizeof(dfloat), nrs->coeffBDF);
  nrs->o_coeffSubEXT = platform->device.malloc(3 * sizeof(dfloat), nrs->coeffEXT);
  nrs->o_coeffRK    = platform->device.malloc(3 * sizeof(dfloat), nrs->coeffRK);

  // define aux kernel constants
  kernelInfo["defines/" "p_eNfields"] = nrs->NVfields;
  kernelInfo["defines/" "p_NVfields"] = nrs->NVfields;
  kernelInfo["defines/" "p_nEXT"] =  nrs->nEXT;
  kernelInfo["defines/" "p_nBDF"] =  nrs->nBDF;
  if(nrs->Nsubsteps)
    kernelInfo["defines/" "p_SUBCYCLING"] =  1;
  else
    kernelInfo["defines/" "p_SUBCYCLING"] =  0;

  kernelInfo["defines/" "p_blockSize"] = BLOCKSIZE;
  //kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  const int movingMesh = nrs->options.compareArgs("MOVING MESH", "TRUE");
  kernelInfo["defines/" "p_MovingMesh"] = movingMesh;

  // jit compile udf kernels
  if (udf.loadKernels) {
    if (platform->comm.mpiRank == 0) cout << "loading udf kernels ... ";
    udf.loadKernels(nrs);
    if (platform->comm.mpiRank == 0) cout << "done" << endl;
  }

  meshParallelGatherScatterSetup(mesh, mesh->Nlocal, mesh->globalIds, platform->comm.mpiComm, 0);
  oogs_mode oogsMode = OOGS_AUTO; 
  if(nrs->options.compareArgs("THREAD MODEL", "SERIAL")) oogsMode = OOGS_DEFAULT;
  nrs->gsh = oogs::setup(mesh->ogs, nrs->NVfields, nrs->fieldOffset, ogsDfloat, NULL, oogsMode);

  linAlg_t * linAlg = nrs->linAlg;

  if(!buildOnly) {
    int err = 0;
    dlong gNelements = mesh->Nelements;
    MPI_Allreduce(MPI_IN_PLACE, &gNelements, 1, MPI_DLONG, MPI_SUM, platform->comm.mpiComm);
    const dfloat sum2 = (dfloat)gNelements * mesh->Np;
    linAlg->fillKernel(nrs->fieldOffset, 1.0, platform->o_slice0);
    ogsGatherScatter(platform->o_slice0, ogsDfloat, ogsAdd, mesh->ogs);
    linAlg->axmyKernel(Nlocal, 1.0, mesh->ogs->o_invDegree, platform->o_slice0); 
    dfloat* tmp = (dfloat*) calloc(Nlocal, sizeof(dfloat));
    platform->o_slice0.copyTo(tmp, Nlocal * sizeof(dfloat));
    dfloat sum1 = 0;
    for(int i = 0; i < Nlocal; i++) sum1 += tmp[i];
    MPI_Allreduce(MPI_IN_PLACE, &sum1, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    sum1 = abs(sum1 - sum2) / sum2;
    if(sum1 > 1e-15) {
      if(platform->comm.mpiRank == 0) printf("ogsGatherScatter test err=%g!\n", sum1);
      fflush(stdout);
      err++;
    }

    mesh->ogs->o_invDegree.copyTo(tmp, Nlocal * sizeof(dfloat));
    double* vmult = (double*) nek::ptr("vmult");
    sum1 = 0;
    for(int i = 0; i < Nlocal; i++) sum1 += abs(tmp[i] - vmult[i]);
    MPI_Allreduce(MPI_IN_PLACE, &sum1, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    if(sum1 > 1e-15) {
      if(platform->comm.mpiRank == 0) printf("multiplicity test err=%g!\n", sum1);
      fflush(stdout);
      err++;
    }

    if(err) ABORT(1);
    free(tmp);
  }

  // setup boundary mapping
  dfloat largeNumber = 1 << 20;
  nrs->VmapB = (int*) calloc(mesh->Nelements * mesh->Np,sizeof(int));
  for (int e = 0; e < mesh->Nelements; e++)
    for (int n = 0; n < mesh->Np; n++) nrs->VmapB[n + e * mesh->Np] = largeNumber;

  nrs->EToB = (int*) calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));

  int cnt = 0;
  for (int e = 0; e < mesh->Nelements; e++)
    for (int f = 0; f < mesh->Nfaces; f++) {
      int bc = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], "velocity");
      nrs->EToB[cnt] = bc;
      if (bc > 0) {
        for (int n = 0; n < mesh->Nfp; n++) {
          int fid = mesh->faceNodes[n + f * mesh->Nfp];
          nrs->VmapB[fid + e * mesh->Np] = mymin(bc,nrs->VmapB[fid + e * mesh->Np]); // Dirichlet wnrs
        }
      }
      cnt++;
    }

  ogsGatherScatter(nrs->VmapB, ogsInt, ogsMin, mesh->ogs);
  for (int n = 0; n < mesh->Nelements * mesh->Np; n++)
    if (nrs->VmapB[n] == largeNumber) nrs->VmapB[n] = 0;

  nrs->o_EToB = device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int),nrs->EToB);
  nrs->o_VmapB = device.malloc(mesh->Nelements * mesh->Np * sizeof(int), nrs->VmapB);

  if(nrs->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
    filterSetup(nrs);

  // build kernels
  string fileName, kernelName;
  const string suffix = "Hex3D";
  const string oklpath = install_dir + "/okl/core/";

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("loading ns kernels ... "); fflush(stdout);

  {

      occa::properties kernelInfoBC = kernelInfo;
      const string bcDataFile = install_dir + "/include/core/bcData.h";
      kernelInfoBC["includes"] += bcDataFile.c_str();
      string boundaryHeaderFileName;
      nrs->options.getArgs("DATA FILE", boundaryHeaderFileName);
      kernelInfoBC["includes"] += realpath(boundaryHeaderFileName.c_str(), NULL);

      fileName = oklpath + "nrsAdvection" + suffix + ".okl";
      kernelName = "nrsStrongAdvectionVolume" + suffix;
      nrs->advectionStrongVolumeKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
      kernelName = "nrsStrongAdvectionCubatureVolume" + suffix;
      nrs->advectionStrongCubatureVolumeKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsCurl" + suffix + ".okl";
      kernelName = "nrsCurl" + suffix;
      nrs->curlKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsGradient" + suffix + ".okl";
      kernelName = "nrsGradientVolume" + suffix;
      nrs->gradientVolumeKernel =  device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "nrswGradientVolume" + suffix;
      nrs->wgradientVolumeKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsSumMakef" + suffix + ".okl";
      kernelName = "nrsSumMakef" + suffix;
      nrs->sumMakefKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsDivergence" + suffix + ".okl";
      kernelName = "nrswDivergenceVolume" + suffix;
      nrs->wDivergenceVolumeKernel =
        platform->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);
      kernelName = "nrsDivergenceVolume" + suffix;
      nrs->divergenceVolumeKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      kernelName = "nrsDivergenceSurfaceTOMBO" + suffix;
      nrs->divergenceSurfaceKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      fileName = oklpath + "nrsAdvectMeshVelocity.okl";
      kernelName = "nrsAdvectMeshVelocity";
      nrs->advectMeshVelocityKernel =
        platform->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      // nrsSurfaceFlux kernel requires that p_blockSize >= p_Nq * p_Nq
      if( BLOCKSIZE < mesh->Nq * mesh->Nq ){
        if(platform->comm.mpiRank == 0)
          printf("ERROR: nrsSurfaceFlux kernel requires BLOCKSIZE >= Nq * Nq."
            "BLOCKSIZE = %d, Nq*Nq = %d\n", BLOCKSIZE, mesh->Nq * mesh->Nq);
        ABORT(EXIT_FAILURE);
      }

      fileName = oklpath + "nrsPressureRhs" + suffix + ".okl";
      kernelName = "nrsPressureRhsTOMBO" + suffix;
      nrs->pressureRhsKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsPressureStress" + suffix + ".okl";
      kernelName = "nrsPressureStress" + suffix;
      nrs->pressureStressKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsPressureBC" + suffix + ".okl";
      kernelName = "nrsPressureDirichletBC" + suffix;
      nrs->pressureDirichletBCKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      fileName = oklpath + "nrsPressureUpdate" + ".okl";
      kernelName = "nrsPressureUpdate";
      nrs->pressureUpdateKernel =  device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "nrsVelocityRhs" + suffix + ".okl";
      kernelName = "nrsVelocityRhsTOMBO" + suffix;
      nrs->velocityRhsKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsVelocityBC" + suffix + ".okl";
      kernelName = "nrsVelocityDirichletBC" + suffix;
      nrs->velocityDirichletBCKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      kernelName = "nrsVelocityNeumannBC" + suffix;
      nrs->velocityNeumannBCKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      fileName = oklpath + "nrsSubCycle" + suffix + ".okl";
      kernelName = "nrsSubCycleStrongCubatureVolume" + suffix;
      nrs->subCycleStrongCubatureVolumeKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      kernelName = "nrsSubCycleStrongVolume" + suffix;
      nrs->subCycleStrongVolumeKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsSubCycleRKUpdate" + ".okl";
      kernelName = "nrsSubCycleLSERKUpdate";
      if(nrs->nRK == 4) kernelName = "nrsSubCycleERKUpdate";
      nrs->subCycleRKUpdateKernel =
        platform->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
      kernelName = "nrsSubCycleRK";
      nrs->subCycleRKKernel =
        platform->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      kernelName = "nrsSubCycleExtrapolateField";
      nrs->subCycleExtrapolateFieldKernel =
        platform->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
      kernelName = "nrsSubCycleExtrapolateScalar";
      nrs->subCycleExtrapolateScalarKernel =
        platform->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsExtrapolate" + ".okl";
      kernelName = "nrsMultiExtrapolate";
      nrs->extrapolateKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      // ===========================================================================

      fileName = oklpath + "mask" + ".okl";

      kernelName = "maskCopy";
      nrs->maskCopyKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      // ===========================================================================

      fileName = oklpath + "nrsFilterRT" + suffix + ".okl";
      kernelName = "nrsFilterRT" + suffix;
      nrs->filterRTKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsCfl" + suffix + ".okl";
      kernelName = "nrsCfl" + suffix;
      nrs->cflKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "nrsPressureAddQtl" + ".okl";
      kernelName = "nrsPressureAddQtl";
      nrs->pressureAddQtlKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "setEllipticCoeff.okl";
      kernelName = "setEllipticCoeff";
      nrs->setEllipticCoeffKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "setEllipticCoeffPressure";
      nrs->setEllipticCoeffPressureKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "nrsPQ.okl";
      kernelName = "nrsPQ";
      nrs->PQKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "nrsMueDiv.okl";
      kernelName = "nrsMueDiv";
      nrs->mueDivKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel); fflush(stdout);

  if(nrs->Nscalar) {
    mesh_t* msh;
    (nrs->cht) ? msh = nrs->meshT : msh = nrs->meshV;
    nrs->cds = cdsSetup(nrs, msh, nrs->options, kernelInfoS);
  }

  if(!buildOnly) {
    // get IC + t0 from nek
    double startTime;
    nek::copyToNek(startTime);
    nrs->options.setArgs("START TIME", to_string_f(startTime));

    if(platform->comm.mpiRank == 0)  printf("calling udf_setup ... "); fflush(stdout);
    udf.setup(nrs);
    if(platform->comm.mpiRank == 0)  printf("done\n"); fflush(stdout);
   }

  // setup elliptic solvers

  const int nbrBIDs = bcMap::size(0);
  int NBCType = nbrBIDs + 1;

  if(nrs->Nscalar) {
    mesh_t* mesh;
    (nrs->cht) ? mesh = nrs->meshT : mesh = nrs->meshV;
    cds_t* cds = nrs->cds;

    for (int is = 0; is < cds->NSfields; is++) {
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(2) << is;
      string sid = ss.str();
 
      if(!cds->compute[is]) continue;
 
      mesh_t* mesh;
      (is) ? mesh = cds->meshV : mesh = cds->meshT[0]; // only first scalar can be a CHT mesh

      if (platform->comm.mpiRank == 0)
        cout << "================= ELLIPTIC SETUP SCALAR" << sid << " ===============\n";

      int nbrBIDs = bcMap::size(0);
      if(nrs->cht && is == 0) nbrBIDs = bcMap::size(1);
      int* sBCType = (int*) calloc(nbrBIDs + 1, sizeof(int));
 
      for (int bID = 1; bID <= nbrBIDs; bID++) {
        string bcTypeText(bcMap::text(bID, "scalar" + sid));
        if(platform->comm.mpiRank == 0) printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str());
        sBCType[bID] = bcMap::type(bID, "scalar" + sid);
      }
 
      cds->solver[is] = new elliptic_t();
      cds->solver[is]->blockSolver = 0;
      cds->solver[is]->Nfields = 1;
      cds->solver[is]->Ntotal = nrs->fieldOffset;
      cds->solver[is]->wrk = platform->mempool + nrs->ellipticWrkOffset;
      cds->solver[is]->o_wrk = platform->o_mempool.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
      cds->solver[is]->mesh = mesh;
      cds->solver[is]->dim = cds->dim;
      cds->solver[is]->elementType = cds->elementType;
      cds->solver[is]->BCType = (int*) calloc(nbrBIDs + 1,sizeof(int));
      memcpy(cds->solver[is]->BCType,sBCType,(nbrBIDs + 1) * sizeof(int));
      free(sBCType);
      cds->solver[is]->var_coeff = cds->var_coeff;
      for (int i = 0; i < 2 * nrs->fieldOffset; i++) nrs->ellipticCoeff[i] = 1;
      cds->solver[is]->lambda = cds->ellipticCoeff;
      cds->solver[is]->o_lambda = cds->o_ellipticCoeff;
      cds->solver[is]->loffset = 0;
 
      cds->solver[is]->options = cds->options[is];
      ellipticSolveSetup(cds->solver[is], kernelInfoS);
    }
  }

  if (nrs->flow) {
    if (platform->comm.mpiRank == 0) printf("================ ELLIPTIC SETUP VELOCITY ================\n");

    nrs->uvwSolver = NULL;

    if(nrs->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE"))
      nrs->uvwSolver = new elliptic_t();

    int* uvwBCType = (int*) calloc(3 * NBCType, sizeof(int));
    int* uBCType = uvwBCType + 0 * NBCType;
    int* vBCType = uvwBCType + 1 * NBCType;
    int* wBCType = uvwBCType + 2 * NBCType;
    for (int bID = 1; bID <= nbrBIDs; bID++) {
      string bcTypeText(bcMap::text(bID, "velocity"));
      if(platform->comm.mpiRank == 0) printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str());

      uBCType[bID] = bcMap::type(bID, "x-velocity");
      vBCType[bID] = bcMap::type(bID, "y-velocity");
      wBCType[bID] = bcMap::type(bID, "z-velocity");
    }

    nrs->vOptions = options;
    nrs->vOptions.setArgs("KRYLOV SOLVER",        options.getArgs("VELOCITY KRYLOV SOLVER"));
    nrs->vOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("VELOCITY SOLVER TOLERANCE"));
    nrs->vOptions.setArgs("DISCRETIZATION",       options.getArgs("VELOCITY DISCRETIZATION"));
    nrs->vOptions.setArgs("BASIS",                options.getArgs("VELOCITY BASIS"));
    nrs->vOptions.setArgs("PRECONDITIONER",       options.getArgs("VELOCITY PRECONDITIONER"));
    nrs->vOptions.setArgs("RESIDUAL PROJECTION",       options.getArgs("VELOCITY RESIDUAL PROJECTION"));
    nrs->vOptions.setArgs("RESIDUAL PROJECTION VECTORS",       options.getArgs("VELOCITY RESIDUAL PROJECTION VECTORS"));
    nrs->vOptions.setArgs("RESIDUAL PROJECTION START",       options.getArgs("VELOCITY RESIDUAL PROJECTION START"));
    nrs->vOptions.setArgs("MULTIGRID COARSENING", options.getArgs("VELOCITY MULTIGRID COARSENING"));
    nrs->vOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("VELOCITY MULTIGRID SMOOTHER"));
    nrs->vOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE",
                          options.getArgs("VELOCITY MULTIGRID CHEBYSHEV DEGREE"));
    nrs->vOptions.setArgs("PARALMOND CYCLE",      options.getArgs("VELOCITY PARALMOND CYCLE"));
    nrs->vOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("VELOCITY PARALMOND SMOOTHER"));
    nrs->vOptions.setArgs("PARALMOND PARTITION",  options.getArgs("VELOCITY PARALMOND PARTITION"));
    nrs->vOptions.setArgs("PARALMOND CHEBYSHEV DEGREE",
                          options.getArgs("VELOCITY PARALMOND CHEBYSHEV DEGREE"));
    nrs->vOptions.setArgs("PARALMOND AGGREGATION STRATEGY",
                          options.getArgs("VELOCITY PARALMOND AGGREGATION STRATEGY"));

    // coeff used by ellipticSetup to detect allNeumann
    for (int i = 0; i < 2 * nrs->fieldOffset; i++) nrs->ellipticCoeff[i] = 1;

    if(nrs->uvwSolver) {
      nrs->uvwSolver->blockSolver = 1;
      nrs->uvwSolver->stressForm = 0;
      if(options.compareArgs("STRESSFORMULATION", "TRUE"))
        nrs->uvwSolver->stressForm = 1;
      nrs->uvwSolver->Nfields = nrs->NVfields;
      nrs->uvwSolver->Ntotal = nrs->fieldOffset;
      nrs->uvwSolver->wrk = platform->mempool + nrs->ellipticWrkOffset;
      nrs->uvwSolver->o_wrk = platform->o_mempool.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
      nrs->uvwSolver->mesh = mesh;
      nrs->uvwSolver->options = nrs->vOptions;
      nrs->uvwSolver->dim = nrs->dim;
      nrs->uvwSolver->elementType = nrs->elementType;
      nrs->uvwSolver->NBCType = NBCType;
      nrs->uvwSolver->BCType = (int*) calloc(nrs->NVfields * NBCType,sizeof(int));
      memcpy(nrs->uvwSolver->BCType,uvwBCType,nrs->NVfields * NBCType * sizeof(int));
      nrs->uvwSolver->var_coeff = nrs->var_coeff;
      nrs->uvwSolver->lambda = nrs->ellipticCoeff;
      nrs->uvwSolver->o_lambda = nrs->o_ellipticCoeff;
      nrs->uvwSolver->loffset = 0; // use same ellipticCoeff for u,v and w

      ellipticSolveSetup(nrs->uvwSolver, kernelInfoV);
    } else {
      nrs->uSolver = new elliptic_t();
      nrs->uSolver->blockSolver = 0;
      nrs->uSolver->Nfields = 1;
      nrs->uSolver->Ntotal = nrs->fieldOffset;
      nrs->uSolver->wrk = platform->mempool + nrs->ellipticWrkOffset;
      nrs->uSolver->o_wrk = platform->o_mempool.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
      nrs->uSolver->mesh = mesh;
      nrs->uSolver->options = nrs->vOptions;
      nrs->uSolver->dim = nrs->dim;
      nrs->uSolver->elementType = nrs->elementType;
      nrs->uSolver->NBCType = NBCType;
      nrs->uSolver->BCType = (int*) calloc(NBCType,sizeof(int));
      memcpy(nrs->uSolver->BCType,uBCType,NBCType * sizeof(int));
      nrs->uSolver->var_coeff = nrs->var_coeff;
      nrs->uSolver->lambda = nrs->ellipticCoeff;
      nrs->uSolver->o_lambda = nrs->o_ellipticCoeff;
      nrs->uSolver->loffset = 0;

      ellipticSolveSetup(nrs->uSolver, kernelInfoV);

      nrs->vSolver = new elliptic_t();
      nrs->vSolver->blockSolver = 0;
      nrs->vSolver->Nfields = 1;
      nrs->vSolver->Ntotal = nrs->fieldOffset;
      nrs->vSolver->wrk = platform->mempool + nrs->ellipticWrkOffset;
      nrs->vSolver->o_wrk = platform->o_mempool.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
      nrs->vSolver->mesh = mesh;
      nrs->vSolver->options = nrs->vOptions;
      nrs->vSolver->dim = nrs->dim;
      nrs->vSolver->elementType = nrs->elementType;
      nrs->vSolver->NBCType = NBCType;
      nrs->vSolver->BCType = (int*) calloc(NBCType,sizeof(int));
      memcpy(nrs->vSolver->BCType,vBCType,NBCType * sizeof(int));
      nrs->vSolver->var_coeff = nrs->var_coeff;
      nrs->vSolver->lambda = nrs->ellipticCoeff;
      nrs->vSolver->o_lambda = nrs->o_ellipticCoeff;
      nrs->vSolver->loffset = 0;

      ellipticSolveSetup(nrs->vSolver, kernelInfoV);

      if (nrs->dim == 3) {
        nrs->wSolver = new elliptic_t();
        nrs->wSolver->blockSolver = 0;
        nrs->wSolver->Nfields = 1;
        nrs->wSolver->Ntotal = nrs->fieldOffset;
        nrs->wSolver->wrk = platform->mempool + nrs->ellipticWrkOffset;
        nrs->wSolver->o_wrk = platform->o_mempool.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
        nrs->wSolver->mesh = mesh;
        nrs->wSolver->options = nrs->vOptions;
        nrs->wSolver->dim = nrs->dim;
        nrs->wSolver->elementType = nrs->elementType;
        nrs->wSolver->NBCType = NBCType;
        nrs->wSolver->BCType = (int*) calloc(NBCType,sizeof(int));
        memcpy(nrs->wSolver->BCType,wBCType,NBCType * sizeof(int));
        nrs->wSolver->var_coeff = nrs->var_coeff;
        nrs->wSolver->lambda = nrs->ellipticCoeff;
        nrs->wSolver->o_lambda = nrs->o_ellipticCoeff;
        nrs->wSolver->loffset = 0;

        ellipticSolveSetup(nrs->wSolver, kernelInfoV);
      }
    }
  } // flow

  if (nrs->flow) {
    if (platform->comm.mpiRank == 0) printf("================ ELLIPTIC SETUP PRESSURE ================\n");

    int* pBCType = (int*) calloc(NBCType, sizeof(int));
    for (int bID = 1; bID <= nbrBIDs; bID++)
      pBCType[bID] = bcMap::type(bID, "pressure");

    nrs->pOptions = options;
    nrs->pOptions.setArgs("KRYLOV SOLVER",        options.getArgs("PRESSURE KRYLOV SOLVER"));
    nrs->pOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("PRESSURE SOLVER TOLERANCE"));
    nrs->pOptions.setArgs("DISCRETIZATION",       options.getArgs("PRESSURE DISCRETIZATION"));
    nrs->pOptions.setArgs("BASIS",                options.getArgs("PRESSURE BASIS"));
    nrs->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRESSURE PRECONDITIONER"));
    nrs->pOptions.setArgs("MULTIGRID COARSENING", options.getArgs("PRESSURE MULTIGRID COARSENING"));
    nrs->pOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
    nrs->pOptions.setArgs("MULTIGRID DOWNWARD SMOOTHER",
                          options.getArgs("PRESSURE MULTIGRID DOWNWARD SMOOTHER"));
    nrs->pOptions.setArgs("MULTIGRID UPWARD SMOOTHER",
                          options.getArgs("PRESSURE MULTIGRID UPWARD SMOOTHER"));
    nrs->pOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE",
                          options.getArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE"));
    nrs->pOptions.setArgs("PARALMOND CYCLE",      options.getArgs("PRESSURE PARALMOND CYCLE"));
    nrs->pOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
    nrs->pOptions.setArgs("PARALMOND PARTITION",  options.getArgs("PRESSURE PARALMOND PARTITION"));
    nrs->pOptions.setArgs("PARALMOND CHEBYSHEV DEGREE",
                          options.getArgs("PRESSURE PARALMOND CHEBYSHEV DEGREE"));
    nrs->pOptions.setArgs("PARALMOND AGGREGATION STRATEGY",
                          options.getArgs("PRESSURE PARALMOND AGGREGATION STRATEGY"));
    nrs->pOptions.setArgs("RESIDUAL PROJECTION", options.getArgs("PRESSURE RESIDUAL PROJECTION"));
    nrs->pOptions.setArgs("RESIDUAL PROJECTION VECTORS",
                          options.getArgs("PRESSURE RESIDUAL PROJECTION VECTORS"));
    nrs->pOptions.setArgs("RESIDUAL PROJECTION START",
                          options.getArgs("PRESSURE RESIDUAL PROJECTION START"));
    nrs->pOptions.setArgs("MULTIGRID VARIABLE COEFFICIENT", "FALSE");

    nrs->pSolver = new elliptic_t();
    nrs->pSolver->blockSolver = 0;
    nrs->pSolver->Nfields = 1;
    nrs->pSolver->Ntotal = nrs->fieldOffset;
    nrs->pSolver->wrk = platform->mempool + nrs->ellipticWrkOffset;
    nrs->pSolver->o_wrk = platform->o_mempool.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
    nrs->pSolver->mesh = mesh;
    nrs->pSolver->dim = nrs->dim;
    nrs->pSolver->elementType = nrs->elementType;
    nrs->pSolver->BCType = (int*) calloc(nbrBIDs + 1,sizeof(int));
    memcpy(nrs->pSolver->BCType,pBCType,(nbrBIDs + 1) * sizeof(int));
    nrs->pSolver->var_coeff = 1;
    //// coeff used by ellipticSetup to detect allNeumann
    // and coeff[0] to setup MG levels
    for (int i = 0; i < 2 * nrs->fieldOffset; i++) nrs->ellipticCoeff[i] = 0;
    nrs->pSolver->lambda = nrs->ellipticCoeff;
    nrs->pSolver->o_lambda = nrs->o_ellipticCoeff;
    nrs->pSolver->loffset = 0;

    string p_mglevels;
    if(nrs->pOptions.getArgs("MULTIGRID COARSENING", p_mglevels)) {
      std::vector<std::string> mgLevelList;
      mgLevelList = serializeString(p_mglevels);
      nrs->pSolver->nLevels = mgLevelList.size();
      nrs->pSolver->levels = (int*) calloc(nrs->pSolver->nLevels,sizeof(int));
      for(int i = 0; i < nrs->pSolver->nLevels; ++i)
        nrs->pSolver->levels[i] = std::atoi(mgLevelList.at(i).c_str());

      if(nrs->pSolver->levels[0] > mesh->N || 
         nrs->pSolver->levels[nrs->pSolver->nLevels-1] < 1) {
        if(platform->comm.mpiRank == 0) printf("ERROR: Invalid multigrid coarsening!\n");
        ABORT(EXIT_FAILURE);;
      }
      nrs->pOptions.setArgs("MULTIGRID COARSENING","CUSTOM");
    } else if(nrs->pOptions.compareArgs("MULTIGRID DOWNWARD SMOOTHER","ASM") ||
              nrs->pOptions.compareArgs("MULTIGRID DOWNWARD SMOOTHER","RAS")) {
      std::map<int,std::vector<int> > mg_level_lookup =
      {
        {1,{1}},
        {2,{2,1}},
        {3,{3,1}},
        {4,{4,2,1}},
        {5,{5,3,1}},
        {6,{6,3,1}},
        {7,{7,3,1}},
        {8,{8,5,1}},
        {9,{9,5,1}},
        {10,{10,6,1}},
        {11,{11,6,1}},
        {12,{12,7,1}},
        {13,{13,7,1}},
        {14,{14,8,1}},
        {15,{15,9,1}},
      };

      const std::vector<int>& levels = mg_level_lookup.at(mesh->Nq - 1);
      nrs->pSolver->nLevels = levels.size();
      nrs->pSolver->levels = (int*) calloc(nrs->pSolver->nLevels,sizeof(int));
      for(int i = 0; i < nrs->pSolver->nLevels; ++i)
        nrs->pSolver->levels[i] = levels.at(i);
      nrs->pOptions.setArgs("MULTIGRID COARSENING","CUSTOM");
    } else if(nrs->pOptions.compareArgs("MULTIGRID DOWNWARD SMOOTHER","JAC")) {
      std::map<int,std::vector<int> > mg_level_lookup =
      {
        {1,{1}},
        {2,{2,1}},
        {3,{3,1}},
        {4,{4,2,1}},
        {5,{5,3,1}},
        {6,{6,4,2,1}},
        {7,{7,5,3,1}},
        {8,{8,6,4,1}},
        {9,{9,7,5,1}},
        {10,{10,8,5,1}},
        {11,{11,9,5,1}},
        {12,{12,10,5,1}},
        {13,{13,11,5,1}},
        {14,{14,12,5,1}},
        {15,{15,13,5,1}},
      };

      const std::vector<int>& levels = mg_level_lookup.at(mesh->Nq - 1);
      nrs->pSolver->nLevels = levels.size();
      nrs->pSolver->levels = (int*) calloc(nrs->pSolver->nLevels,sizeof(int));
      for(int i = 0; i < nrs->pSolver->nLevels; ++i)
        nrs->pSolver->levels[i] = levels.at(i);
      nrs->pOptions.setArgs("MULTIGRID COARSENING","CUSTOM");
    }

    nrs->pSolver->options = nrs->pOptions;
    ellipticSolveSetup(nrs->pSolver, kernelInfoP);

  } // flow
}

static cds_t* cdsSetup(nrs_t* nrs, mesh_t* mesh, setupAide options, occa::properties &kernelInfoH)
{
  cds_t* cds = new cds_t();
  platform_t* platform = platform_t::getInstance();
  device_t& device = platform->device;
  cds->meshT[0] = mesh;

  string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));

  // set mesh, options
  cds->meshV       = nrs->meshV;
  cds->elementType = nrs->elementType;
  cds->dim         = nrs->dim;
  cds->NVfields    = nrs->NVfields;
  cds->NSfields    = nrs->Nscalar;

  cds->coeffEXT = nrs->coeffEXT;
  cds->coeffBDF = nrs->coeffBDF;
  cds->coeffSubEXT = nrs->coeffSubEXT;
  cds->coeffRK    = nrs->coeffRK;

  cds->nBDF       = nrs->nBDF;
  cds->nEXT       = nrs->nEXT;

  // time stepper
  dfloat rkC[4]  = {1.0, 0.0, -1.0, -2.0};
  cds->o_coeffEXT = nrs->o_coeffEXT;
  cds->o_coeffBDF = nrs->o_coeffBDF;
  cds->o_coeffSubEXT = nrs->o_coeffSubEXT;
  cds->o_coeffRK    = nrs->o_coeffRK;

  cds->o_usrwrk = &(nrs->o_usrwrk);

  dlong Nlocal = mesh->Np * mesh->Nelements;
  dlong Ntotal = mesh->Np * (mesh->Nelements + mesh->totalHaloPairs);
  mesh->Nlocal = Nlocal;
  cds->fieldOffset  = Ntotal;

  cds->vFieldOffset = nrs->fieldOffset;
  cds->fieldOffset  = nrs->fieldOffset;

  platform->o_slice0 = platform->o_slice0;
  platform->o_slice1 = platform->o_slice1;
  platform->o_slice2 = platform->o_slice2;
  platform->o_slice3 = platform->o_slice3;
  platform->o_slice4 = platform->o_slice4;
  platform->o_slice5 = platform->o_slice5;
  platform->o_slice6 = platform->o_slice6;

  cds->gsh = nrs->gsh;
  
  if(nrs->cht) {
    meshParallelGatherScatterSetup(mesh, mesh->Nlocal, mesh->globalIds, platform->comm.mpiComm, 0);
    oogs_mode oogsMode = OOGS_AUTO; 
    if(options.compareArgs("THREAD MODEL", "SERIAL")) oogsMode = OOGS_DEFAULT;
    cds->gshT = oogs::setup(mesh->ogs, 1, cds->fieldOffset, ogsDfloat, NULL, oogsMode);
  } else {
    cds->gshT = cds->gsh;
  }

  // Solution storage at interpolation nodes
  cds->U     = nrs->U; // Point to INS side Velocity
  cds->S     =
    (dfloat*) calloc(cds->NSfields * cds->nBDF * cds->fieldOffset,sizeof(dfloat));
  cds->BF    = (dfloat*) calloc(cds->NSfields * cds->fieldOffset,sizeof(dfloat));
  cds->FS    =
    (dfloat*) calloc(cds->NSfields * cds->nBDF * cds->fieldOffset,sizeof(dfloat));

  cds->Nsubsteps = nrs->Nsubsteps;
  if(cds->Nsubsteps) {
    cds->nRK   = nrs->nRK;
    cds->Srka   = nrs->Srka;
    cds->Srkb   = nrs->Srkb;
    cds->Srkc   = nrs->Srkc;
    cds->o_Srka = nrs->o_Srka;
    cds->o_Srkb = nrs->o_Srkb;
  }

  cds->dt  = nrs->dt;
  cds->sdt = nrs->sdt;

  cds->prop = (dfloat*) calloc(cds->NSfields * 2 * cds->fieldOffset,sizeof(dfloat));
  for(int is = 0; is < cds->NSfields; is++) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << is;
    string sid = ss.str();

    if(options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) continue;

    dfloat diff = 1;
    dfloat rho = 1;
    options.getArgs("SCALAR" + sid + " DIFFUSIVITY", diff);
    options.getArgs("SCALAR" + sid + " DENSITY", rho);

    const dlong off = cds->NSfields * cds->fieldOffset;
    for (int e = 0; e < mesh->Nelements; e++)
      for (int n = 0; n < mesh->Np; n++) {
        cds->prop[0 * off + is * cds->fieldOffset + e * mesh->Np + n] = diff;
        cds->prop[1 * off + is * cds->fieldOffset + e * mesh->Np + n] = rho;
      }
  }
  cds->o_prop =
    device.malloc(cds->NSfields * 2 * cds->fieldOffset * sizeof(dfloat), cds->prop);
  cds->o_diff = cds->o_prop.slice(0 * cds->NSfields * cds->fieldOffset * sizeof(dfloat));
  cds->o_rho  = cds->o_prop.slice(1 * cds->NSfields * cds->fieldOffset * sizeof(dfloat));

  cds->var_coeff = 1; // use always var coeff elliptic
  cds->ellipticCoeff   = nrs->ellipticCoeff;
  cds->o_ellipticCoeff = nrs->o_ellipticCoeff;

  cds->o_U  = nrs->o_U;
  cds->o_Ue = nrs->o_Ue;
  cds->o_S  =
    platform->device.malloc(cds->NSfields * cds->nBDF * cds->fieldOffset * sizeof(dfloat), cds->S);
  cds->o_Se =
    platform->device.malloc(cds->NSfields * cds->nBDF * cds->fieldOffset * sizeof(dfloat));
  cds->o_BF = platform->device.malloc(cds->NSfields * cds->fieldOffset * sizeof(dfloat), cds->BF);
  cds->o_FS =
    platform->device.malloc(cds->NSfields * cds->nEXT * cds->fieldOffset * sizeof(dfloat),
                        cds->FS);

  for (int is = 0; is < cds->NSfields; is++) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << is;
    string sid = ss.str();

    cds->compute[is] = 1;
    if (options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
      cds->compute[is] = 0;
      continue;
    }

    mesh_t* mesh;
    (is) ? mesh = cds->meshV : mesh = cds->meshT[0]; // only first scalar can be a CHT mesh
 
    cds->options[is] = options;

    cds->options[is].setArgs("KRYLOV SOLVER", options.getArgs("SCALAR SOLVER"));
    cds->options[is].setArgs("DISCRETIZATION", options.getArgs("SCALAR DISCRETIZATION"));
    cds->options[is].setArgs("BASIS", options.getArgs("SCALAR BASIS"));
    cds->options[is].setArgs("PRECONDITIONER", options.getArgs("SCALAR" + sid + " PRECONDITIONER"));
    cds->options[is].setArgs("SOLVER TOLERANCE",
                         options.getArgs("SCALAR" + sid +  " SOLVER TOLERANCE"));
    cds->options[is].setArgs("RESIDUAL PROJECTION",  options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION"));
    cds->options[is].setArgs("RESIDUAL PROJECTION VECTORS",  options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION VECTORS"));
    cds->options[is].setArgs("RESIDUAL PROJECTION START",  options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION START"));

    // setup boundary mapping
    dfloat largeNumber = 1 << 20;
    cds->mapB[is] = (int*) calloc(mesh->Nelements * mesh->Np,sizeof(int));
    int* mapB = cds->mapB[is];
    for (int e = 0; e < mesh->Nelements; e++)
      for (int n = 0; n < mesh->Np; n++) mapB[n + e * mesh->Np] = largeNumber;

    cds->EToB[is] = (int*) calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
    int* EToB = cds->EToB[is];

    int cnt = 0;
    for (int e = 0; e < mesh->Nelements; e++)
      for (int f = 0; f < mesh->Nfaces; f++) {
        int bc = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], "scalar" + sid);
        EToB[cnt] = bc;
        if (bc > 0) {
          for (int n = 0; n < mesh->Nfp; n++) {
            int fid = mesh->faceNodes[n + f * mesh->Nfp];
            mapB[fid + e * mesh->Np] = mymin(bc,mapB[fid + e * mesh->Np]); // Dirichlet wnrs
          }
        }
        cnt++;
      }

    ogsGatherScatter(mapB, ogsInt, ogsMin, mesh->ogs);

    for (int n = 0; n < mesh->Nelements * mesh->Np; n++)
      if (mapB[n] == largeNumber) mapB[n] = 0;

    cds->o_EToB[is] = device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int), EToB);
    cds->o_mapB[is] = device.malloc(mesh->Nelements * mesh->Np * sizeof(int), mapB);
  }

  // build kernels
  occa::properties kernelInfo = *nrs->kernelInfo;
  occa::properties kernelInfoBC = kernelInfo;
  //kernelInfo["defines/" "p_NSfields"]  = cds->NSfields;

  string fileName, kernelName;
  const string suffix = "Hex3D";
  const string oklpath = install_dir + "/okl/core/";

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("loading cds kernels ... "); fflush(stdout);

   {
      fileName = oklpath + "cdsAdvection" + suffix + ".okl";

      const string bcDataFile = install_dir + "/include/core/bcData.h";
      kernelInfoBC["includes"] += bcDataFile.c_str();
      string boundaryHeaderFileName;
      options.getArgs("DATA FILE", boundaryHeaderFileName);
      kernelInfoBC["includes"] += realpath(boundaryHeaderFileName.c_str(), NULL);

      kernelName = "cdsStrongAdvectionVolume" + suffix;
      cds->advectionStrongVolumeKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "cdsStrongAdvectionCubatureVolume" + suffix;
      cds->advectionStrongCubatureVolumeKernel =  device.buildKernel(fileName,
                                                                           kernelName,
                                                                           kernelInfo);

      // ===========================================================================
      fileName = oklpath + "cdsAdvectMeshVelocity.okl";
      kernelName = "cdsAdvectMeshVelocity";
      cds->advectMeshVelocityKernel =
        platform->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "mask.okl";
      kernelName = "maskCopy";
      cds->maskCopyKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName   = oklpath + "cdsSumMakef" + suffix + ".okl";
      kernelName = "cdsSumMakef" + suffix;
      cds->sumMakefKernel =  device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "cdsHelmholtzBC" + suffix + ".okl";
      kernelName = "cdsHelmholtzBC" + suffix;
      cds->helmholtzRhsBCKernel =  device.buildKernel(fileName, kernelName, kernelInfoBC);

      kernelName = "cdsDirichletBC";
      cds->dirichletBCKernel =  device.buildKernel(fileName, kernelName, kernelInfoBC);

      fileName = oklpath + "setEllipticCoeff.okl";
      kernelName = "setEllipticCoeff";
      cds->setEllipticCoeffKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "cdsFilterRT" + suffix + ".okl";
      kernelName = "cdsFilterRT" + suffix;
      cds->filterRTKernel =
        device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      if(cds->Nsubsteps) {
        fileName = oklpath + "cdsSubCycle" + suffix + ".okl";
        kernelName = "cdsSubCycleStrongCubatureVolume" + suffix;
        cds->subCycleStrongCubatureVolumeKernel =  device.buildKernel(fileName,
                                                                            kernelName,
                                                                            kernelInfo);

        kernelName = "cdsSubCycleStrongVolume" + suffix;
        cds->subCycleStrongVolumeKernel =
          device.buildKernel(fileName, kernelName, kernelInfo);

        fileName = oklpath + "cdsSubCycleRKUpdate.okl";
        kernelName = "cdsSubCycleLSERKUpdate";
        if(cds->nRK == 4) kernelName = "cdsSubCycleERKUpdate";
        cds->subCycleRKUpdateKernel =  platform->device.buildKernel(fileName, kernelName, kernelInfo);
        kernelName = "cdsSubCycleRK";
        cds->subCycleRKKernel =  platform->device.buildKernel(fileName, kernelName, kernelInfo);
        kernelName = "subCycleExtrapolateField";
        cds->subCycleExtrapolateFieldKernel =  platform->device.buildKernel(fileName, kernelName, kernelInfo);
        kernelName = "cdsSubCycleExtrapolateScalar";
        cds->subCycleExtrapolateScalarKernel =  platform->device.buildKernel(fileName, kernelName, kernelInfo);
      }
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel); fflush(stdout);

  return cds;
}
