#include "nrs.hpp"
#include "meshSetup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "bcMap.hpp"
#include <vector>
#include <map>
#include "filter.hpp"
#include "avm.hpp"

namespace{
cds_t* cdsSetup(nrs_t* nrs, setupAide options, occa::properties &kernelInfoBC);
}

void nrsSetup(MPI_Comm comm, setupAide &options, nrs_t *nrs)
{
  platform_t* platform = platform_t::getInstance();
  device_t& device = platform->device;
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
  if(platform->options.compareArgs("BUILD ONLY", "TRUE")) buildOnly = 1;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  platform->options.getArgs("CUBATURE POLYNOMIAL DEGREE", cubN);
  platform->options.getArgs("NUMBER OF SCALARS", nrs->Nscalar);
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  platform->options.getArgs("MESH DIMENSION", nrs->dim);
  platform->options.getArgs("ELEMENT TYPE", nrs->elementType);
  if(platform->device.mode() == "Serial")
    platform->options.setArgs("ENABLE OVERLAP", "FALSE");

  nrs->flow = 1;
  if(platform->options.compareArgs("VELOCITY", "FALSE")) nrs->flow = 0;
  if(platform->options.compareArgs("VELOCITY SOLVER", "NONE")) nrs->flow = 0;

  if(nrs->flow) {
    if(platform->options.compareArgs("STRESSFORMULATION", "TRUE"))
       platform->options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
  }


  // jit compile + init nek
  {  
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    string casename;
    platform->options.getArgs("CASENAME", casename);

    int err = 0;
    int npTarget = size;
    if (buildOnly) platform->options.getArgs("NP TARGET", npTarget);
    if (rank == 0) err = buildNekInterface(casename.c_str(), mymax(5, nrs->Nscalar), N, npTarget, platform->options);
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, comm);
    if (err) ABORT(EXIT_FAILURE);; 

    if (!buildOnly) {
      nek::setup(comm, platform->options, nrs);
      nek::setic();
      nek::userchk();
      if (platform->comm.mpiRank == 0) cout << "\n";
    }
  }

  nrs->cht = 0;
  if (nekData.nelv != nekData.nelt && nrs->Nscalar) nrs->cht = 1;
  if (nrs->cht && !platform->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE")) {
    if (platform->comm.mpiRank == 0) cout << "Conjugate heat transfer requires solving for temperature!\n"; 
    ABORT(EXIT_FAILURE);;
  } 

  nrs->_mesh = createMesh(comm, N, cubN, nrs->cht, kernelInfo);
  nrs->meshV = (mesh_t*) nrs->_mesh->fluid;
  mesh_t* mesh = nrs->meshV;

  { 
    double val = (double)mesh->NlocalGatherElements/mesh->Nelements; 
    MPI_Allreduce(MPI_IN_PLACE,&val,1,MPI_DOUBLE,MPI_MIN,platform->comm.mpiComm);
    if(platform->comm.mpiRank == 0) 
      printf("min %2.0f%% of the local elements are internal\n", 100*val);
  }

  occa::properties kernelInfoV  = kernelInfo;
  occa::properties kernelInfoP  = kernelInfo;
  occa::properties kernelInfoS  = kernelInfo;

  nrs->NVfields = 3;
  nrs->NTfields = nrs->NVfields + 1;   // Total Velocity + Pressure
  mesh->Nfields = 1;

  platform->options.getArgs("SUBCYCLING STEPS",nrs->Nsubsteps);
  platform->options.getArgs("DT", nrs->dt[0]);

  if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO1")) {
    nrs->nBDF = 1;
  } else if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO2")) {
    nrs->nBDF = 2;
  } else if (platform->options.compareArgs("TIME INTEGRATOR", "TOMBO3")) {
    nrs->nBDF = 3;
  }
  nrs->nEXT = 3;
  if(nrs->Nsubsteps) nrs->nEXT = nrs->nBDF;
  nrs->coeffEXT = (dfloat*) calloc(nrs->nEXT, sizeof(dfloat));
  nrs->coeffBDF = (dfloat*) calloc(nrs->nBDF, sizeof(dfloat));

  nrs->nRK = 4;
  nrs->coeffSubEXT = (dfloat*) calloc(3, sizeof(dfloat));

  dfloat mue = 1;
  dfloat rho = 1;
  platform->options.getArgs("VISCOSITY", mue);
  platform->options.getArgs("DENSITY", rho);

  const dlong Nlocal = mesh->Nlocal;

  { // setup fieldOffset
    nrs->fieldOffset = mesh->Np * (mesh->Nelements + mesh->totalHaloPairs);
    mesh_t* meshT = nrs->_mesh; 
    nrs->fieldOffset = mymax(nrs->fieldOffset, meshT->Np * (meshT->Nelements + meshT->totalHaloPairs));

    int PAGESIZE = 4096; // default is 4kB
    char* tmp;
    tmp = getenv("NEKRS_PAGE_SIZE");
    if (tmp != NULL) PAGESIZE = std::stoi(tmp);
    const int pageW = PAGESIZE / sizeof(dfloat);
    if (nrs->fieldOffset % pageW) nrs->fieldOffset = (nrs->fieldOffset / pageW + 1) * pageW;
  }

  nrs->_mesh->fieldOffset = nrs->fieldOffset;

  if(nrs->Nsubsteps) {
    int Sorder;
    platform->options.getArgs("SUBCYCLING TIME ORDER", Sorder);
    if(Sorder == 4 && nrs->nRK == 4) { // ERK(4,4)
      dfloat rka[4] = {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0};
      dfloat rkb[4] = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
      dfloat rkc[4] = {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0};
      nrs->coeffsfRK = (dfloat*) calloc(nrs->nRK, sizeof(dfloat));
      nrs->weightsRK = (dfloat*) calloc(nrs->nRK, sizeof(dfloat));
      nrs->nodesRK = (dfloat*) calloc(nrs->nRK, sizeof(dfloat));
      memcpy(nrs->coeffsfRK, rka, nrs->nRK * sizeof(dfloat));
      memcpy(nrs->weightsRK, rkb, nrs->nRK * sizeof(dfloat));
      memcpy(nrs->nodesRK, rkc, nrs->nRK * sizeof(dfloat));
    }else{
      if(platform->comm.mpiRank == 0) cout << "Unsupported subcycling scheme!\n";
      ABORT(1);
    }
    nrs->o_coeffsfRK = device.malloc(nrs->nRK * sizeof(dfloat), nrs->coeffsfRK);
    nrs->o_weightsRK = device.malloc(nrs->nRK * sizeof(dfloat), nrs->weightsRK);
  }

  // setup scratch space
  const int wrkNflds = 6;
  const int ellipticWrkNflds = 15;
  nrs->ellipticWrkOffset = wrkNflds * nrs->fieldOffset;

  const int scratchNflds = wrkNflds + ellipticWrkNflds;
  platform->create_mempool(nrs->fieldOffset, scratchNflds);

  if(options.compareArgs("MOVING MESH", "TRUE")){
    const int nBDF = std::max(nrs->nBDF, nrs->nEXT);
    platform->o_mempool.slice0.copyFrom(mesh->o_LMM, mesh->Nlocal * sizeof(dfloat));
    mesh->o_LMM.free();
    mesh->o_LMM = platform->device.malloc(nrs->fieldOffset * nBDF ,  sizeof(dfloat));
    mesh->o_LMM.copyFrom(platform->o_mempool.slice0, mesh->Nlocal * sizeof(dfloat));
    platform->o_mempool.slice0.copyFrom(mesh->o_invLMM, mesh->Nlocal * sizeof(dfloat));
    mesh->o_invLMM.free();
    mesh->o_invLMM = platform->device.malloc(nrs->fieldOffset * nBDF ,  sizeof(dfloat));
    mesh->o_invLMM.copyFrom(platform->o_mempool.slice0, mesh->Nlocal * sizeof(dfloat));

    const int nAB = std::max(nrs->nEXT, mesh->nAB);
    mesh->U = (dfloat*) calloc(nrs->NVfields * nrs->fieldOffset * nAB, sizeof(dfloat));
    mesh->o_U = platform->device.malloc(nrs->NVfields * nrs->fieldOffset * nAB * sizeof(dfloat), mesh->U);
    if(nrs->Nsubsteps)
      mesh->o_divU = platform->device.malloc(nrs->fieldOffset * nAB, sizeof(dfloat));
  }

  {
    dlong offset;
    if(platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
      offset = std::max(nrs->fieldOffset, mesh->Nelements * mesh->cubNp);
    else
      offset = nrs->fieldOffset;
 
    const dlong Nstates = nrs->Nsubsteps ? std::max(nrs->nBDF, nrs->nEXT) : 1;
    if(nrs->Nsubsteps && platform->options.compareArgs("MOVING MESH", "TRUE"))
      nrs->o_relUrst = platform->device.malloc(Nstates * nrs->NVfields * offset, sizeof(dfloat));
    else
      nrs->o_Urst = platform->device.malloc(Nstates * nrs->NVfields * offset, sizeof(dfloat));
  }

  nrs->U  = (dfloat*) calloc(nrs->NVfields * std::max(nrs->nBDF, nrs->nEXT) * nrs->fieldOffset,sizeof(dfloat));
  nrs->Ue = (dfloat*) calloc(nrs->NVfields * nrs->fieldOffset,sizeof(dfloat));
  nrs->P  = (dfloat*) calloc(nrs->fieldOffset,sizeof(dfloat));
  nrs->BF = (dfloat*) calloc(nrs->NVfields * nrs->fieldOffset,sizeof(dfloat));
  nrs->FU = (dfloat*) calloc(nrs->NVfields * nrs->nEXT * nrs->fieldOffset,sizeof(dfloat));

  nrs->o_U  = platform->device.malloc(nrs->NVfields * std::max(nrs->nBDF,nrs->nEXT) * nrs->fieldOffset * sizeof(dfloat), nrs->U);
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

  nrs->o_coeffEXT = platform->device.malloc(nrs->nEXT * sizeof(dfloat), nrs->coeffEXT);
  nrs->o_coeffBDF = platform->device.malloc(nrs->nBDF * sizeof(dfloat), nrs->coeffBDF);
  nrs->o_coeffSubEXT = platform->device.malloc(nrs->nEXT * sizeof(dfloat), nrs->coeffEXT);

  // define aux kernel constants
  kernelInfo["defines/" "p_eNfields"] = nrs->NVfields;
  kernelInfo["defines/" "p_NVfields"] = nrs->NVfields;

  occa::properties kernelInfoBC = *(nrs->kernelInfo);

  // jit compile udf kernels
  if (udf.loadKernels) {
    occa::properties* tmpKernelInfo = nrs->kernelInfo;
    nrs->kernelInfo = &kernelInfoBC;
    if (platform->comm.mpiRank == 0) cout << "loading udf kernels ... ";
    udf.loadKernels(nrs);
    if (platform->comm.mpiRank == 0) cout << "done" << endl;
    nrs->kernelInfo = tmpKernelInfo;
  }
  const string bcDataFile = install_dir + "/include/core/bcData.h";
  kernelInfoBC["includes"] += bcDataFile.c_str();
  string boundaryHeaderFileName;
  platform->options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfoBC["includes"] += realpath(boundaryHeaderFileName.c_str(), NULL);


  meshParallelGatherScatterSetup(mesh, mesh->Nlocal, mesh->globalIds, platform->comm.mpiComm, 0);
  oogs_mode oogsMode = OOGS_AUTO; 
  //if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP") oogsMode = OOGS_DEFAULT;
  nrs->gsh = oogs::setup(mesh->ogs, nrs->NVfields, nrs->fieldOffset, ogsDfloat, NULL, oogsMode);

  linAlg_t * linAlg = platform->linAlg;

  if(!buildOnly) {
    int err = 0;
    dlong gNelements = mesh->Nelements;
    MPI_Allreduce(MPI_IN_PLACE, &gNelements, 1, MPI_DLONG, MPI_SUM, platform->comm.mpiComm);
    const dfloat sum2 = (dfloat)gNelements * mesh->Np;
    linAlg->fillKernel(nrs->fieldOffset, 1.0, platform->o_mempool.slice0);
    ogsGatherScatter(platform->o_mempool.slice0, ogsDfloat, ogsAdd, mesh->ogs);
    linAlg->axmyKernel(Nlocal, 1.0, mesh->ogs->o_invDegree, platform->o_mempool.slice0); 
    dfloat* tmp = (dfloat*) calloc(Nlocal, sizeof(dfloat));
    platform->o_mempool.slice0.copyTo(tmp, Nlocal * sizeof(dfloat));
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

  if(platform->options.compareArgs("FILTER STABILIZATION", "RELAXATION")){

    nrs->filterNc = -1;
    dfloat filterS;
    platform->options.getArgs("HPFRT STRENGTH", filterS);
    platform->options.getArgs("HPFRT MODES", nrs->filterNc);
    filterS = -1.0 * fabs(filterS);
    nrs->filterS = filterS;

    dfloat* A = filterSetup(nrs->meshV, nrs->filterNc);

    const dlong Nmodes = nrs->meshV->N + 1;

    nrs->o_filterMT = platform->device.malloc(Nmodes * Nmodes * sizeof(dfloat), A);

    free(A);
  }

  // build kernels
  string fileName, kernelName;
  const string suffix = "Hex3D";
  const string oklpath = install_dir + "/okl/";

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("loading ns kernels ... "); fflush(stdout);

  {
      fileName = oklpath + "core/nStagesSum.okl";
      kernelName = "nStagesSum3";
      nrs->nStagesSum3Kernel =
        device.buildKernel(fileName, kernelName, platform->kernelInfo);

      {
        occa::properties prop = kernelInfo;
        prop["defines/" "p_cubNq"] = nrs->meshV->cubNq;
        prop["defines/" "p_cubNp"] = nrs->meshV->cubNp;
	fileName = oklpath + "nrs/advection" + suffix + ".okl";
        kernelName = "strongAdvectionVolume" + suffix;
        nrs->advectionStrongVolumeKernel =
          device.buildKernel(fileName, kernelName, prop);
        kernelName = "strongAdvectionCubatureVolume" + suffix;
        nrs->advectionStrongCubatureVolumeKernel =
          device.buildKernel(fileName, kernelName, prop);
      }

      fileName = oklpath + "nrs/curl" + suffix + ".okl";
      kernelName = "curl" + suffix;
      nrs->curlKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "nrs/gradient" + suffix + ".okl";
      kernelName = "gradientVolume" + suffix;
      nrs->gradientVolumeKernel =  device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "nrswGradientVolume" + suffix;
      nrs->wgradientVolumeKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      {
        occa::properties prop = kernelInfo;
        const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
        prop["defines/" "p_nEXT"] =  nrs->nEXT;
        prop["defines/" "p_nBDF"] =  nrs->nBDF;
        prop["defines/" "p_MovingMesh"] = movingMesh;
        if(nrs->Nsubsteps)
          prop["defines/" "p_SUBCYCLING"] = 1;
        else
          prop["defines/" "p_SUBCYCLING"] = 0;
          
        fileName = oklpath + "nrs/sumMakef.okl";
        kernelName = "sumMakef";
        nrs->sumMakefKernel =  device.buildKernel(fileName, kernelName, prop);
      }

      fileName = oklpath + "nrs/divergence" + suffix + ".okl";
      kernelName = "nrswDivergenceVolume" + suffix;
      nrs->wDivergenceVolumeKernel =
        platform->device.buildKernel(fileName, kernelName, kernelInfoBC);
      kernelName = "divergenceVolume" + suffix;
      nrs->divergenceVolumeKernel =
        device.buildKernel(fileName, kernelName, kernelInfoBC);

      kernelName = "divergenceSurfaceTOMBO" + suffix;
      nrs->divergenceSurfaceKernel =
        device.buildKernel(fileName, kernelName, kernelInfoBC);

      fileName = oklpath + "nrs/advectMeshVelocityHex3D.okl";
      kernelName = "advectMeshVelocityHex3D";
      nrs->advectMeshVelocityKernel =
        platform->device.buildKernel(fileName, kernelName, kernelInfo);

      // nrsSurfaceFlux kernel requires that p_blockSize >= p_Nq * p_Nq
      if( BLOCKSIZE < mesh->Nq * mesh->Nq ){
        if(platform->comm.mpiRank == 0)
          printf("ERROR: nrsSurfaceFlux kernel requires BLOCKSIZE >= Nq * Nq."
            "BLOCKSIZE = %d, Nq*Nq = %d\n", BLOCKSIZE, mesh->Nq * mesh->Nq);
        ABORT(EXIT_FAILURE);
      }

      fileName = oklpath + "nrs/pressureRhs" + suffix + ".okl";
      kernelName = "pressureRhsTOMBO" + suffix;
      nrs->pressureRhsKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "nrs/pressureStress" + suffix + ".okl";
      kernelName = "pressureStress" + suffix;
      nrs->pressureStressKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "nrs/pressureBC" + suffix + ".okl";
      kernelName = "pressureDirichletBC" + suffix;
      nrs->pressureDirichletBCKernel =
        device.buildKernel(fileName, kernelName, kernelInfoBC);

      fileName = oklpath + "nrs/pressureUpdate" + ".okl";
      kernelName = "pressureUpdate";
      nrs->pressureUpdateKernel =  device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "nrs/velocityRhs" + suffix + ".okl";
      kernelName = "velocityRhsTOMBO" + suffix;
      nrs->velocityRhsKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "nrs/velocityBC" + suffix + ".okl";
      kernelName = "velocityDirichletBC" + suffix;
      nrs->velocityDirichletBCKernel =
        device.buildKernel(fileName, kernelName, kernelInfoBC);

      kernelName = "velocityNeumannBC" + suffix;
      nrs->velocityNeumannBCKernel =
        device.buildKernel(fileName, kernelName, kernelInfoBC);

      occa::properties prop = kernelInfo;
      const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
      prop["defines/" "p_relative"] = movingMesh && nrs->Nsubsteps;
      prop["defines/" "p_cubNq"] =  nrs->meshV->cubNq;
      prop["defines/" "p_cubNp"] =  nrs->meshV->cubNp;
      fileName = oklpath + "nrs/Urst" + suffix + ".okl";
      kernelName = "UrstCubature" + suffix;
      nrs->UrstCubatureKernel =
        device.buildKernel(fileName, kernelName, prop);

      kernelName = "Urst" + suffix;
      nrs->UrstKernel =
        device.buildKernel(fileName, kernelName, prop);


      if(nrs->Nsubsteps){
        occa::properties prop = kernelInfo;
        const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
        prop["defines/" "p_MovingMesh"] = movingMesh;
        prop["defines/" "p_nEXT"] =  nrs->nEXT;
        prop["defines/" "p_nBDF"] =  nrs->nBDF;
        prop["defines/" "p_cubNq"] =  nrs->meshV->cubNq;
        prop["defines/" "p_cubNp"] =  nrs->meshV->cubNp;
	
        fileName = oklpath + "nrs/subCycle" + suffix + ".okl";
        occa::properties subCycleStrongCubatureProps = prop;
        if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP"){
          fileName = oklpath + "nrs/subCycle" + suffix + ".c";
          subCycleStrongCubatureProps["okl/enabled"] = false;
        }
        kernelName = "subCycleStrongCubatureVolume" + suffix;
        nrs->subCycleStrongCubatureVolumeKernel =
          device.buildKernel(fileName, kernelName, subCycleStrongCubatureProps);
        fileName = oklpath + "nrs/subCycle" + suffix + ".okl";
        nrs->subCycleStrongVolumeKernel =
          device.buildKernel(fileName, kernelName, prop);

        fileName = oklpath + "nrs/subCycleRKUpdate" + ".okl";
        kernelName = "subCycleLSERKUpdate";
        if(nrs->nRK == 4) kernelName = "subCycleERKUpdate";
        nrs->subCycleRKUpdateKernel =
          platform->device.buildKernel(fileName, kernelName, prop);
        kernelName = "subCycleRK";
        nrs->subCycleRKKernel =
          platform->device.buildKernel(fileName, kernelName, prop);

        kernelName = "subCycleInitU0";
        nrs->subCycleInitU0Kernel =  platform->device.buildKernel(fileName, kernelName, prop);
      }

      fileName = oklpath + "nrs/extrapolate" + ".okl";
      kernelName = "multiExtrapolate";
      nrs->extrapolateKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "core/mask" + ".okl";
      kernelName = "maskCopy";
      nrs->maskCopyKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "nrs/regularization/filterRT" + suffix + ".okl";
      kernelName = "filterRT" + suffix;
      nrs->filterRTKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      occa::properties cflProps = kernelInfo;
      cflProps["defines/ " "p_MovingMesh"] = movingMesh;
      fileName = oklpath + "nrs/cfl" + suffix + ".okl";
      kernelName = "cfl" + suffix;
      if( BLOCKSIZE < mesh->Nq * mesh->Nq ){
        if(platform->comm.mpiRank == 0)
          printf("ERROR: cfl kernel requires BLOCKSIZE >= Nq * Nq."
            "BLOCKSIZE = %d, Nq*Nq = %d\n", BLOCKSIZE, mesh->Nq * mesh->Nq);
        ABORT(EXIT_FAILURE);
      }
      nrs->cflKernel =
        device.buildKernel(fileName, kernelName, cflProps);

      fileName = oklpath + "nrs/pressureAddQtl" + ".okl";
      kernelName = "pressureAddQtl";
      nrs->pressureAddQtlKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "core/setEllipticCoeff.okl";
      kernelName = "setEllipticCoeff";
      nrs->setEllipticCoeffKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);
      kernelName = "setEllipticCoeffPressure";
      nrs->setEllipticCoeffPressureKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "nrs/mueDiv.okl";
      kernelName = "mueDiv";
      nrs->mueDivKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel); fflush(stdout);

  if(nrs->Nscalar) {
    nrs->cds = cdsSetup(nrs, platform->options, kernelInfoBC);
  }

  if(!buildOnly) {
    // get IC + t0 from nek
    double startTime;
    nek::copyFromNek(startTime);
    platform->options.setArgs("START TIME", to_string_f(startTime));

    if(platform->comm.mpiRank == 0)  printf("calling udf_setup ... "); fflush(stdout);
    udf.setup(nrs);
    if(platform->comm.mpiRank == 0)  printf("done\n"); fflush(stdout);
   }

  // setup elliptic solvers

  const int nbrBIDs = bcMap::size(0);
  int NBCType = nbrBIDs + 1;

  if(nrs->Nscalar) {
    cds_t* cds = nrs->cds;

    for (int is = 0; is < cds->NSfields; is++) {
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(2) << is;
      string sid = ss.str();
 
      if(!cds->compute[is]) continue;
 
      mesh_t* mesh;
      (is) ? mesh = cds->meshV : mesh = cds->mesh[0]; // only first scalar can be a CHT mesh

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
      cds->solver[is]->name = "scalar" + sid;
      cds->solver[is]->blockSolver = 0;
      cds->solver[is]->Nfields = 1;
      cds->solver[is]->Ntotal = nrs->fieldOffset;
      cds->solver[is]->wrk = platform->mempool.slice0 + nrs->ellipticWrkOffset;
      cds->solver[is]->o_wrk = platform->o_mempool.o_ptr.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
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

    if(platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE"))
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
    nrs->vOptions.setArgs("PGMRES RESTART",        options.getArgs("VELOCITY PGMRES RESTART"));
    nrs->vOptions.setArgs("KRYLOV SOLVER",        options.getArgs("VELOCITY KRYLOV SOLVER"));
    nrs->vOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("VELOCITY SOLVER TOLERANCE"));
    nrs->vOptions.setArgs("DISCRETIZATION",       options.getArgs("VELOCITY DISCRETIZATION"));
    nrs->vOptions.setArgs("BASIS",                options.getArgs("VELOCITY BASIS"));
    nrs->vOptions.setArgs("PRECONDITIONER",       options.getArgs("VELOCITY PRECONDITIONER"));
    nrs->vOptions.setArgs("RESIDUAL PROJECTION",       options.getArgs("VELOCITY RESIDUAL PROJECTION"));
    nrs->vOptions.setArgs("RESIDUAL PROJECTION VECTORS",       options.getArgs("VELOCITY RESIDUAL PROJECTION VECTORS"));
    nrs->vOptions.setArgs("RESIDUAL PROJECTION START",       options.getArgs("VELOCITY RESIDUAL PROJECTION START"));
    nrs->vOptions.setArgs("RESIDUAL PROJECTION METHOD",       options.getArgs("VELOCITY RESIDUAL PROJECTION METHOD"));
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
    nrs->vOptions.setArgs("MAXIMUM ITERATIONS", options.getArgs("VELOCITY MAXIMUM ITERATIONS"));

    // coeff used by ellipticSetup to detect allNeumann
    for (int i = 0; i < 2 * nrs->fieldOffset; i++) nrs->ellipticCoeff[i] = 1;

    if(nrs->uvwSolver) {
      nrs->uvwSolver->blockSolver = 1;
      nrs->uvwSolver->stressForm = 0;
      if(options.compareArgs("STRESSFORMULATION", "TRUE"))
        nrs->uvwSolver->stressForm = 1;
      nrs->uvwSolver->Nfields = nrs->NVfields;
      nrs->uvwSolver->Ntotal = nrs->fieldOffset;
      nrs->uvwSolver->wrk = platform->mempool.slice0 + nrs->ellipticWrkOffset;
      nrs->uvwSolver->o_wrk = platform->o_mempool.o_ptr.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
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
      nrs->uSolver->wrk = platform->mempool.slice0 + nrs->ellipticWrkOffset;
      nrs->uSolver->o_wrk = platform->o_mempool.o_ptr.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
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
      nrs->vSolver->wrk = platform->mempool.slice0 + nrs->ellipticWrkOffset;
      nrs->vSolver->o_wrk = platform->o_mempool.o_ptr.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
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
        nrs->wSolver->wrk = platform->mempool.slice0 + nrs->ellipticWrkOffset;
        nrs->wSolver->o_wrk = platform->o_mempool.o_ptr.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
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

    if(platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE")) {
      nrs->uvwSolver->name = "velocity";
    } else {
      nrs->uSolver->name = "x-velocity";
      nrs->vSolver->name = "y-velocity";
      nrs->wSolver->name = "v-velocity";
    }
  } // flow

  if (nrs->flow) {
    if (platform->comm.mpiRank == 0) printf("================ ELLIPTIC SETUP PRESSURE ================\n");

    int* pBCType = (int*) calloc(NBCType, sizeof(int));
    for (int bID = 1; bID <= nbrBIDs; bID++)
      pBCType[bID] = bcMap::type(bID, "pressure");

    nrs->pOptions = options;
    nrs->pOptions.setArgs("PGMRES RESTART",       options.getArgs("PRESSURE PGMRES RESTART"));
    nrs->pOptions.setArgs("KRYLOV SOLVER",        options.getArgs("PRESSURE KRYLOV SOLVER"));
    nrs->pOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("PRESSURE SOLVER TOLERANCE"));
    nrs->pOptions.setArgs("DISCRETIZATION",       options.getArgs("PRESSURE DISCRETIZATION"));
    nrs->pOptions.setArgs("BASIS",                options.getArgs("PRESSURE BASIS"));
    nrs->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRESSURE PRECONDITIONER"));
    nrs->pOptions.setArgs("SEMFEM SOLVER", options.getArgs("PRESSURE SEMFEM SOLVER"));
    nrs->pOptions.setArgs("SEMFEM SOLVER PRECISION", options.getArgs("PRESSURE SEMFEM SOLVER PRECISION"));
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
    nrs->pOptions.setArgs("RESIDUAL PROJECTION METHOD",       options.getArgs("PRESSURE RESIDUAL PROJECTION METHOD"));
    nrs->pOptions.setArgs("MULTIGRID VARIABLE COEFFICIENT", "FALSE");
    nrs->pOptions.setArgs("MAXIMUM ITERATIONS", options.getArgs("PRESSURE MAXIMUM ITERATIONS"));

    nrs->pSolver = new elliptic_t();
    nrs->pSolver->name = "pressure";
    nrs->pSolver->blockSolver = 0;
    nrs->pSolver->Nfields = 1;
    nrs->pSolver->Ntotal = nrs->fieldOffset;
    nrs->pSolver->wrk = platform->mempool.slice0 + nrs->ellipticWrkOffset;
    nrs->pSolver->o_wrk = platform->o_mempool.o_ptr.slice(nrs->ellipticWrkOffset * sizeof(dfloat));
    nrs->pSolver->mesh = mesh;
    nrs->pSolver->dim = nrs->dim;
    nrs->pSolver->elementType = nrs->elementType;
    nrs->pSolver->BCType = (int*) calloc(nbrBIDs + 1,sizeof(int));
    memcpy(nrs->pSolver->BCType,pBCType,(nbrBIDs + 1) * sizeof(int));

    nrs->pSolver->var_coeff = 1;

    // coeff used by ellipticSetup to detect allNeumann
    for (int i = 0; i < 2 * nrs->fieldOffset; i++) nrs->ellipticCoeff[i] = 0;
    nrs->pSolver->lambda = nrs->ellipticCoeff;
    nrs->pSolver->o_lambda = nrs->o_ellipticCoeff;
    nrs->pSolver->loffset = 0;

    string p_mglevels;
    if(nrs->pOptions.getArgs("MULTIGRID COARSENING", p_mglevels)) {
      std::vector<std::string> mgLevelList;
      mgLevelList = serializeString(p_mglevels,',');
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

namespace{
cds_t* cdsSetup(nrs_t* nrs, setupAide options, occa::properties& kernelInfoBC)
{
  cds_t* cds = new cds_t();
  platform_t* platform = platform_t::getInstance();
  device_t& device = platform->device;
  string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));

  cds->mesh[0]     = nrs->_mesh;
  mesh_t* mesh     = cds->mesh[0];
  cds->meshV       = nrs->_mesh->fluid;
  cds->elementType = nrs->elementType;
  cds->dim         = nrs->dim;
  cds->NVfields    = nrs->NVfields;
  cds->NSfields    = nrs->Nscalar;

  cds->coeffEXT    = nrs->coeffEXT;
  cds->coeffBDF    = nrs->coeffBDF;
  cds->coeffSubEXT = nrs->coeffSubEXT;
  cds->nBDF        = nrs->nBDF;
  cds->nEXT        = nrs->nEXT;
  cds->o_coeffEXT  = nrs->o_coeffEXT;
  cds->o_coeffBDF  = nrs->o_coeffBDF;
  cds->o_coeffSubEXT = nrs->o_coeffSubEXT;

  cds->o_usrwrk = &(nrs->o_usrwrk);

  cds->vFieldOffset = nrs->fieldOffset;
  cds->fieldOffset[0]  = nrs->fieldOffset;
  cds->fieldOffsetScan[0] = 0;
  dlong sum = cds->fieldOffset[0];
  for(int s = 1; s < cds->NSfields; ++s){
    cds->fieldOffset[s] = cds->fieldOffset[0];
    cds->fieldOffsetScan[s] = sum;
    sum += cds->fieldOffset[s];
    cds->mesh[s] = cds->mesh[0];
  }
  cds->fieldOffsetSum = sum;

  cds->gsh = nrs->gsh;
  
  if(nrs->cht) {
    meshParallelGatherScatterSetup(mesh, mesh->Nlocal, mesh->globalIds, platform->comm.mpiComm, 0);
    oogs_mode oogsMode = OOGS_AUTO; 
    //if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP") oogsMode = OOGS_DEFAULT;
    cds->gshT = oogs::setup(mesh->ogs, 1, cds->fieldOffset[0], ogsDfloat, NULL, oogsMode);
  } else {
    cds->gshT = cds->gsh;
  }

  // Solution storage at interpolation nodes
  cds->U     = nrs->U; // Point to INS side Velocity
  cds->S     =
    (dfloat*) calloc(std::max(cds->nBDF, cds->nEXT) * cds->fieldOffsetSum,sizeof(dfloat));
  cds->BF    = (dfloat*) calloc(cds->fieldOffsetSum,sizeof(dfloat));
  cds->FS    =
    (dfloat*) calloc(cds->nBDF * cds->fieldOffsetSum,sizeof(dfloat));

  cds->Nsubsteps = nrs->Nsubsteps;
  if(cds->Nsubsteps) {
    cds->nRK   = nrs->nRK;
    cds->coeffsfRK   = nrs->coeffsfRK;
    cds->weightsRK   = nrs->weightsRK;
    cds->nodesRK   = nrs->nodesRK;
    cds->o_coeffsfRK = nrs->o_coeffsfRK;
    cds->o_weightsRK = nrs->o_weightsRK;
  }

  cds->dt  = nrs->dt;
  cds->sdt = nrs->sdt;

  cds->prop = (dfloat*) calloc(2 * cds->fieldOffsetSum,sizeof(dfloat));



  for(int is = 0; is < cds->NSfields; is++) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << is;
    string sid = ss.str();

    if(options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) continue;

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
  cds->o_prop =
    device.malloc(2 * cds->fieldOffsetSum * sizeof(dfloat), cds->prop);
  cds->o_diff = cds->o_prop.slice(0 * cds->fieldOffsetSum * sizeof(dfloat));
  cds->o_rho  = cds->o_prop.slice(1 * cds->fieldOffsetSum * sizeof(dfloat));

  cds->var_coeff = 1; // use always var coeff elliptic
  cds->ellipticCoeff   = nrs->ellipticCoeff;
  cds->o_ellipticCoeff = nrs->o_ellipticCoeff;

  cds->o_U  = nrs->o_U;
  cds->o_Ue = nrs->o_Ue;
  cds->o_S  =
    platform->device.malloc(std::max(cds->nBDF, cds->nEXT) * cds->fieldOffsetSum * sizeof(dfloat), cds->S);
  cds->o_Se =
    platform->device.malloc(cds->fieldOffsetSum ,  sizeof(dfloat));
  cds->o_BF = platform->device.malloc(cds->fieldOffsetSum * sizeof(dfloat), cds->BF);
  cds->o_FS =
    platform->device.malloc(cds->nEXT * cds->fieldOffsetSum * sizeof(dfloat),
                        cds->FS);

  cds->o_relUrst = nrs->o_relUrst;
  cds->o_Urst = nrs->o_Urst;

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
    (is) ? mesh = cds->meshV : mesh = cds->mesh[0]; // only first scalar can be a CHT mesh
 
    cds->options[is] = options;

    cds->options[is].setArgs("RAMP CONSTANT", options.getArgs("SCALAR" + sid + " RAMP CONSTANT"));
    cds->options[is].setArgs("AVM C0", options.getArgs("SCALAR" + sid + " AVM C0"));
    cds->options[is].setArgs("FILTER STABILIZATION", options.getArgs("SCALAR" + sid + " FILTER STABILIZATION"));
    cds->options[is].setArgs("VISMAX COEFF", options.getArgs("SCALAR" + sid + " VISMAX COEFF"));
    cds->options[is].setArgs("HPFRT STRENGTH", options.getArgs("SCALAR" + sid + " HPFRT STRENGTH"));
    cds->options[is].setArgs("HPFRT MODES", options.getArgs("SCALAR" + sid + " HPFRT MODES"));
    cds->options[is].setArgs("KRYLOV SOLVER", options.getArgs("SCALAR" + sid + " KRYLOV SOLVER"));
    cds->options[is].setArgs("PGMRES RESTART", options.getArgs("SCALAR" + sid + " PGMRES RESTART"));
    cds->options[is].setArgs("DISCRETIZATION", options.getArgs("SCALAR DISCRETIZATION"));
    cds->options[is].setArgs("BASIS", options.getArgs("SCALAR BASIS"));
    cds->options[is].setArgs("PRECONDITIONER", options.getArgs("SCALAR" + sid + " PRECONDITIONER"));
    cds->options[is].setArgs("SOLVER TOLERANCE",
                         options.getArgs("SCALAR" + sid +  " SOLVER TOLERANCE"));
    cds->options[is].setArgs("RESIDUAL PROJECTION",  options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION"));
    cds->options[is].setArgs("RESIDUAL PROJECTION VECTORS",  options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION VECTORS"));
    cds->options[is].setArgs("RESIDUAL PROJECTION START",  options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION START"));
    cds->options[is].setArgs("RESIDUAL PROJECTION METHOD",  options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION METHOD"));
    cds->options[is].setArgs("MAXIMUM ITERATIONS", options.getArgs("SCALAR MAXIMUM ITERATIONS"));

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

  bool scalarFilteringEnabled = false;
  bool avmEnabled = false;
  {
    for(int is = 0; is < cds->NSfields; is++) {
      if(!cds->options[is].compareArgs("FILTER STABILIZATION", "NONE")) scalarFilteringEnabled = true;
      if(cds->options[is].compareArgs("FILTER STABILIZATION", "AVM")) avmEnabled = true;
    }
  }

  if(scalarFilteringEnabled)
  {
    const dlong Nmodes = cds->mesh[0]->N + 1;
    cds->o_filterMT = platform->device.malloc(cds->NSfields * Nmodes * Nmodes, sizeof(dfloat));
    for(int is = 0; is < cds->NSfields; is++)
    {
      if(cds->options[is].compareArgs("FILTER STABILIZATION", "NONE")) continue;
      int filterNc = -1;
      cds->options[is].getArgs("HPFRT MODES", filterNc);
      dfloat filterS;
      cds->options[is].getArgs("HPFRT STRENGTH", filterS);
      filterS = -1.0 * fabs(filterS);
      cds->filterS[is] = filterS;

      dfloat* A = filterSetup(cds->mesh[is], filterNc);

      const dlong Nmodes = cds->mesh[is]->N + 1;
      cds->o_filterMT.copyFrom(A, Nmodes * Nmodes * sizeof(dfloat), is * Nmodes * Nmodes * sizeof(dfloat));

      free(A);
    }
  }

  if(avmEnabled) avm::setup(cds);

  // build kernels
  occa::properties kernelInfo = *nrs->kernelInfo;
  //kernelInfo["defines/" "p_NSfields"]  = cds->NSfields;

  string fileName, kernelName;
  const string suffix = "Hex3D";
  const string oklpath = install_dir + "/okl/";

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("loading cds kernels ... "); fflush(stdout);

   {
      {
        occa::properties prop = kernelInfo;
        prop["defines/" "p_cubNq"] = cds->mesh[0]->cubNq;
        prop["defines/" "p_cubNp"] = cds->mesh[0]->cubNp;
        fileName = oklpath + "cds/advection" + suffix + ".okl";

	kernelName = "strongAdvectionVolume" + suffix;
        cds->advectionStrongVolumeKernel =
          device.buildKernel(fileName, kernelName, prop);

	kernelName = "strongAdvectionCubatureVolume" + suffix;
        cds->advectionStrongCubatureVolumeKernel =  
          device.buildKernel(fileName, kernelName, prop);
      }

      fileName = oklpath + "cds/advectMeshVelocityHex3D.okl";
  	kernelName = "advectMeshVelocityHex3D";
      cds->advectMeshVelocityKernel =
        platform->device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "core/mask.okl";
      kernelName = "maskCopy";
      cds->maskCopyKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      {
        occa::properties prop = kernelInfo;
        const int movingMesh = options.compareArgs("MOVING MESH", "TRUE");
        prop["defines/" "p_MovingMesh"] = movingMesh;
        prop["defines/" "p_nEXT"] =  cds->nEXT;
        prop["defines/" "p_nBDF"] =  cds->nBDF;
        if(cds->Nsubsteps)
          prop["defines/" "p_SUBCYCLING"] = 1;
        else
          prop["defines/" "p_SUBCYCLING"] = 0;
          
        fileName   = oklpath + "cds/sumMakef.okl";
        kernelName = "sumMakef";
        cds->sumMakefKernel =  device.buildKernel(fileName, kernelName, prop);

      }

      fileName = oklpath + "cds/helmholtzBC" + suffix + ".okl";
      kernelName = "helmholtzBC" + suffix;
      cds->helmholtzRhsBCKernel =  device.buildKernel(fileName, kernelName, kernelInfoBC);
      kernelName = "dirichletBC";
      cds->dirichletBCKernel =  device.buildKernel(fileName, kernelName, kernelInfoBC);

      fileName = oklpath + "core/setEllipticCoeff.okl";
      kernelName = "setEllipticCoeff";
      cds->setEllipticCoeffKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "cds/regularization/filterRT" + suffix + ".okl";
      kernelName = "filterRT" + suffix;
      cds->filterRTKernel =
        device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "core/nStagesSum.okl";
      kernelName = "nStagesSum3";
      cds->nStagesSum3Kernel =
        device.buildKernel(fileName, kernelName, platform->kernelInfo);

      if(cds->Nsubsteps) {
        occa::properties prop = kernelInfo;
        const int movingMesh = options.compareArgs("MOVING MESH", "TRUE");
        prop["defines/" "p_MovingMesh"] = movingMesh;
        prop["defines/" "p_nEXT"] =  cds->nEXT;
        prop["defines/" "p_nBDF"] =  cds->nBDF;
        prop["defines/" "p_cubNq"] =  cds->mesh[0]->cubNq;
        prop["defines/" "p_cubNp"] =  cds->mesh[0]->cubNp;
 

        fileName = oklpath + "cds/subCycle" + suffix + ".okl";
        occa::properties subCycleStrongCubatureProps = prop;
        if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP"){
          fileName = oklpath + "cds/subCycle" + suffix + ".c";
          subCycleStrongCubatureProps["okl/enabled"] = false;
        }
        kernelName = "subCycleStrongCubatureVolume" + suffix;
        cds->subCycleStrongCubatureVolumeKernel =
          device.buildKernel(fileName, kernelName, subCycleStrongCubatureProps);
        fileName = oklpath + "cds/subCycle" + suffix + ".okl";
        kernelName = "subCycleStrongVolume" + suffix;
        cds->subCycleStrongVolumeKernel =
          device.buildKernel(fileName, kernelName, prop);


        fileName = oklpath + "cds/subCycleRKUpdate.okl";
        kernelName = "subCycleLSERKUpdate";
        if(cds->nRK == 4) kernelName = "subCycleERKUpdate";
        cds->subCycleRKUpdateKernel =  platform->device.buildKernel(fileName, kernelName, prop);
        kernelName = "subCycleRK";
        cds->subCycleRKKernel =  platform->device.buildKernel(fileName, kernelName, prop);

        kernelName = "subCycleInitU0";
        cds->subCycleInitU0Kernel =  platform->device.buildKernel(fileName, kernelName, prop);
      }
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel); fflush(stdout);

  return cds;
}
}
