#include "nrs.hpp"
#include "meshSetup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "bcMap.hpp"
#include <vector>
#include <map>
#include "filter.hpp"
#include "avm.hpp"
#include <cctype>

namespace{
cds_t* cdsSetup(nrs_t* nrs, setupAide options);
}

std::vector<int>
determineMGLevels(std::string section)
{
  const std::string optionsPrefix = [section](){
    std::string prefix = section + std::string(" ");
    if(section.find("temperature") != std::string::npos){
      prefix = std::string("scalar00 ");
    }
    std::transform(prefix.begin(), prefix.end(), prefix.begin(), 
               [](unsigned char c){ return std::toupper(c); });
    return prefix;
  }();

  std::vector<int> levels;
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);

  std::string p_mglevels;
  if(platform->options.getArgs(optionsPrefix + "MULTIGRID COARSENING", p_mglevels)) {
    const std::vector<std::string> mgLevelList = serializeString(p_mglevels,',');
    for(auto && s : mgLevelList){
      levels.push_back(std::stoi(s));
    }

    
    bool invalid = false;
    invalid |= (levels[0] != N); // top level order must match
    for(unsigned i = 0U; i < levels.size(); ++i){
      invalid |= (levels[i] < 0); // each level must be positive
      if(i > 0)
        invalid |= (levels[i] >= levels[i-1]); // each successive level must be smaller
    }

    if(invalid){
      if(platform->comm.mpiRank == 0) printf("ERROR: Invalid multigrid coarsening!\n");
      ABORT(EXIT_FAILURE);;
    }
    if(levels.back() > 1)
    {
      if(platform->options.compareArgs(optionsPrefix + "MULTIGRID COARSE SOLVE", "TRUE")){
        // if the coarse level has p > 1 and requires solving the coarsest level,
        // rather than just smoothing, then use the SEMFEM discretization
        platform->options.setArgs(optionsPrefix + "MULTIGRID COARSE SEMFEM", "TRUE");
        platform->options.setArgs(optionsPrefix + "MULTIGRID COARSE SEMFEM", "TRUE");

        // However, if the user explicitly asked for the FEM discretization, bail
        if(platform->options.compareArgs(optionsPrefix + "USER SPECIFIED FEM COARSE SOLVER", "TRUE"))
        {
          if(platform->comm.mpiRank == 0){
            printf("Error! FEM coarse discretization only supports p=1 for the coarsest level!\n");
          }
          ABORT(1);
        }
      }
    }

    return levels;

  } else if(platform->options.compareArgs(optionsPrefix + "MULTIGRID DOWNWARD SMOOTHER","ASM") ||
            platform->options.compareArgs(optionsPrefix + "MULTIGRID DOWNWARD SMOOTHER","RAS")) {
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

    return mg_level_lookup.at(N);
  }
  else {
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

    return mg_level_lookup.at(N);
  }
}

void nrsSetup(MPI_Comm comm, setupAide &options, nrs_t *nrs)
{
  {
    int N;
    platform->options.getArgs("POLYNOMIAL DEGREE", N);
    const int Nq = N+1;
    if( BLOCKSIZE < Nq * Nq ){
      if(platform->comm.mpiRank == 0)
        printf("ERROR: several kernels require BLOCKSIZE >= Nq * Nq."
          "BLOCKSIZE = %d, Nq*Nq = %d\n", BLOCKSIZE, Nq * Nq);
      ABORT(EXIT_FAILURE);
    }
  }

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
  platform->options.getArgs("POLYNOMIAL DEGREE", N);
  platform->options.getArgs("CUBATURE POLYNOMIAL DEGREE", cubN);
  platform->options.getArgs("NUMBER OF SCALARS", nrs->Nscalar);
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

  if(platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE"))
  {
    platform->options.getArgs("FLOW RATE", nrs->flowRate);
    nrs->fromBID = -1;
    nrs->toBID = -1;
    platform->options.getArgs("CONSTANT FLOW FROM BID", nrs->fromBID);
    platform->options.getArgs("CONSTANT FLOW TO BID", nrs->toBID);
    if(platform->options.compareArgs("CONSTANT FLOW DIRECTION", "X"))
    {
      nrs->flowDirection[0] = 1.0;
      nrs->flowDirection[1] = 0.0;
      nrs->flowDirection[2] = 0.0;
    }
    if(platform->options.compareArgs("CONSTANT FLOW DIRECTION", "Y"))
    {
      nrs->flowDirection[0] = 0.0;
      nrs->flowDirection[1] = 1.0;
      nrs->flowDirection[2] = 0.0;
    }
    if(platform->options.compareArgs("CONSTANT FLOW DIRECTION", "Z"))
    {
      nrs->flowDirection[0] = 0.0;
      nrs->flowDirection[1] = 0.0;
      nrs->flowDirection[2] = 1.0;
    }
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
    if (platform->comm.mpiRank == 0) std::cout << "\n";
  }

  nrs->cht = 0;
  if (nekData.nelv != nekData.nelt && nrs->Nscalar) nrs->cht = 1;
  if (nrs->cht && !platform->options.compareArgs("SCALAR00 IS TEMPERATURE", "TRUE")) {
    if (platform->comm.mpiRank == 0) std::cout << "Conjugate heat transfer requires solving for temperature!\n"; 
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

    const int pageW = ALIGN_SIZE / sizeof(dfloat);
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
      if(platform->comm.mpiRank == 0) std::cout << "Unsupported subcycling scheme!\n";
      ABORT(1);
    }
    nrs->o_coeffsfRK = device.malloc(nrs->nRK * sizeof(dfloat), nrs->coeffsfRK);
    nrs->o_weightsRK = device.malloc(nrs->nRK * sizeof(dfloat), nrs->weightsRK);
  }

  // setup mempool 
  int ellipticMaxFields = 1;
  if(platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE"))
    ellipticMaxFields = nrs->NVfields;
  const int ellipticWrkFields = elliptic_t::NScratchFields*ellipticMaxFields;

  int wrkFields = 9;
  if(nrs->Nsubsteps) wrkFields += 3*nrs->NVfields;
  if(options.compareArgs("MOVING MESH", "TRUE")) wrkFields += nrs->NVfields;

  const int mempoolNflds = std::max(wrkFields, 2*nrs->NVfields + ellipticWrkFields);
  platform->create_mempool(nrs->fieldOffset, mempoolNflds);

  // offset mempool available for elliptic because also used it for ellipticSolve input/output  
  auto const o_mempoolElliptic = 
    platform->o_mempool.o_ptr.slice(2*nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));

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

  nrs->o_ellipticCoeff = device.malloc(2 * nrs->fieldOffset * sizeof(dfloat));

  int nProperties = 2;
  if(options.compareArgs("MESH SOLVER", "ELASTICITY")) nProperties = 4;
  nrs->prop =  (dfloat*) calloc(nProperties * nrs->fieldOffset,sizeof(dfloat));
  for (int e = 0; e < mesh->Nelements; e++)
    for (int n = 0; n < mesh->Np; n++) {
      nrs->prop[0 * nrs->fieldOffset + e * mesh->Np + n] = mue;
      nrs->prop[1 * nrs->fieldOffset + e * mesh->Np + n] = rho;
    }

  nrs->o_prop = device.malloc(nProperties * nrs->fieldOffset * sizeof(dfloat), nrs->prop);
  nrs->o_mue = nrs->o_prop.slice(0 * nrs->fieldOffset * sizeof(dfloat));
  nrs->o_rho = nrs->o_prop.slice(1 * nrs->fieldOffset * sizeof(dfloat));
  if(options.compareArgs("MESH SOLVER", "ELASTICITY")){
    nrs->o_meshMue = nrs->o_prop.slice(2 * nrs->fieldOffset * sizeof(dfloat));
    nrs->o_meshRho = nrs->o_prop.slice(3 * nrs->fieldOffset * sizeof(dfloat));
  }

  if(platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")){
    nrs->o_Uc  = platform->device.malloc(nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
    nrs->o_Pc  = platform->device.malloc(nrs->fieldOffset * sizeof(dfloat));
    nrs->o_prevProp = device.malloc(2 * nrs->fieldOffset * sizeof(dfloat), nrs->prop);
  }

  nrs->div   = (dfloat*) calloc(nrs->fieldOffset,sizeof(dfloat));
  nrs->o_div = device.malloc(nrs->fieldOffset * sizeof(dfloat), nrs->div);

  nrs->o_coeffEXT = platform->device.malloc(nrs->nEXT * sizeof(dfloat), nrs->coeffEXT);
  nrs->o_coeffBDF = platform->device.malloc(nrs->nBDF * sizeof(dfloat), nrs->coeffBDF);
  nrs->o_coeffSubEXT = platform->device.malloc(nrs->nEXT * sizeof(dfloat), nrs->coeffEXT);

  meshParallelGatherScatterSetup(mesh, mesh->Nlocal, mesh->globalIds, platform->comm.mpiComm, 0);
  oogs_mode oogsMode = OOGS_AUTO; 
  //if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP") oogsMode = OOGS_DEFAULT;
  nrs->gsh = oogs::setup(mesh->ogs, nrs->NVfields, nrs->fieldOffset, ogsDfloat, NULL, oogsMode);

  linAlg_t * linAlg = platform->linAlg;

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

  nrs->EToB = (int*) calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
  int cnt = 0;
  for (int e = 0; e < mesh->Nelements; e++) {
    for (int f = 0; f < mesh->Nfaces; f++) {
      int bc = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], "velocity");
      nrs->EToB[cnt] = bc;
      cnt++;
    }
  }
  nrs->o_EToB = device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int),nrs->EToB);

  if(platform->options.compareArgs("MESH SOLVER", "ELASTICITY")) {
    nrs->EToBMesh = (int*) calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
    int cnt = 0;
    for (int e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        int bc = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], "mesh");
        nrs->EToBMesh[cnt] = bc;
        cnt++;
      }
    }
    nrs->o_EToBMesh = device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int),nrs->EToBMesh);
  }

  if(platform->options.compareArgs("VELOCITY REGULARIZATION METHOD", "RELAXATION")){

    nrs->filterNc = -1;
    dfloat filterS;
    platform->options.getArgs("VELOCITY HPFRT STRENGTH", filterS);
    platform->options.getArgs("VELOCITY HPFRT MODES", nrs->filterNc);
    filterS = -1.0 * fabs(filterS);
    nrs->filterS = filterS;

    dfloat* A = filterSetup(nrs->meshV, nrs->filterNc);

    const dlong Nmodes = nrs->meshV->N + 1;

    nrs->o_filterMT = platform->device.malloc(Nmodes * Nmodes * sizeof(dfloat), A);

    free(A);
  }

  // build kernels
  std::string kernelName;
  const std::string suffix = "Hex3D";

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("loading ns kernels ... "); fflush(stdout);

  {
      const std::string section = "nrs-";
      kernelName = "nStagesSum3";
      nrs->nStagesSum3Kernel =
        platform->kernels.get( section + kernelName);

      kernelName = "computeFieldDotNormal";
      nrs->computeFieldDotNormalKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "computeFaceCentroid";
      nrs->computeFaceCentroidKernel =
        platform->kernels.get( section + kernelName);

      {
        kernelName = "strongAdvectionVolume" + suffix;
        nrs->strongAdvectionVolumeKernel = platform->kernels.get(section + kernelName);
        kernelName = "strongAdvectionCubatureVolume" + suffix;
        nrs->strongAdvectionCubatureVolumeKernel = platform->kernels.get(section + kernelName);
      }

      kernelName = "curl" + suffix;
      nrs->curlKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "gradientVolume" + suffix;
      nrs->gradientVolumeKernel =  platform->kernels.get( section + kernelName);

      kernelName = "wGradientVolume" + suffix;
      nrs->wgradientVolumeKernel =
        platform->kernels.get( section + kernelName);

      {
        kernelName = "sumMakef";
        nrs->sumMakefKernel =  platform->kernels.get( section + kernelName);
      }

      kernelName = "wDivergenceVolume" + suffix;
      nrs->wDivergenceVolumeKernel =
        platform->kernels.get( section + kernelName);
      kernelName = "divergenceVolume" + suffix;
      nrs->divergenceVolumeKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "divergenceSurface" + suffix;
      nrs->divergenceSurfaceKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "advectMeshVelocity" + suffix;
      nrs->advectMeshVelocityKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "pressureRhs" + suffix;
      nrs->pressureRhsKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "pressureStress" + suffix;
      nrs->pressureStressKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "pressureDirichletBC" + suffix;
      nrs->pressureDirichletBCKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "velocityRhs" + suffix;
      nrs->velocityRhsKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "velocityDirichletBC" + suffix;
      nrs->velocityDirichletBCKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "velocityNeumannBC" + suffix;
      nrs->velocityNeumannBCKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "UrstCubature" + suffix;
      nrs->UrstCubatureKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "Urst" + suffix;
      nrs->UrstKernel =
        platform->kernels.get( section + kernelName);


      if(nrs->Nsubsteps){
        kernelName = "subCycleStrongCubatureVolume" + suffix;
        nrs->subCycleStrongCubatureVolumeKernel =
          platform->kernels.get( section + kernelName);
        kernelName = "subCycleStrongVolume" + suffix;
        nrs->subCycleStrongVolumeKernel =
          platform->kernels.get( section + kernelName);

        kernelName = "subCycleRKUpdate";
        nrs->subCycleRKUpdateKernel =
          platform->kernels.get( section + kernelName);
        kernelName = "subCycleRK";
        nrs->subCycleRKKernel =
          platform->kernels.get( section + kernelName);

        kernelName = "subCycleInitU0";
        nrs->subCycleInitU0Kernel =  platform->kernels.get( section + kernelName);
      }

      kernelName = "extrapolate";
      nrs->extrapolateKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "maskCopy";
      nrs->maskCopyKernel =
        platform->kernels.get( section + kernelName);
      kernelName = "mask";
      nrs->maskKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "filterRT" + suffix;
      nrs->filterRTKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "cfl" + suffix;
      nrs->cflKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "pressureAddQtl";
      nrs->pressureAddQtlKernel =
        platform->kernels.get( section + kernelName);

      kernelName = "setEllipticCoeff";
      nrs->setEllipticCoeffKernel =
        platform->kernels.get( section + kernelName);
      kernelName = "setEllipticCoeffPressure";
      nrs->setEllipticCoeffPressureKernel =
        platform->kernels.get( section + kernelName);
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel); fflush(stdout);

  if(nrs->Nscalar) {
    nrs->cds = cdsSetup(nrs, platform->options);
  }

  // get IC + t0 from nek
  double startTime;
  nek::copyFromNek(startTime);
  platform->options.setArgs("START TIME", to_string_f(startTime));

  if(platform->comm.mpiRank == 0)  printf("calling udf_setup ... "); fflush(stdout);
  udf.setup(nrs);
  if(platform->comm.mpiRank == 0)  printf("done\n"); fflush(stdout);

  nrs->o_U.copyFrom(nrs->U);
  nrs->o_P.copyFrom(nrs->P);
  nrs->o_prop.copyFrom(nrs->prop);
  if(nrs->Nscalar) {
    nrs->cds->o_S.copyFrom(nrs->cds->S);
    nrs->cds->o_prop.copyFrom(nrs->cds->prop);
  }

  evaluateProperties(nrs, startTime);
  nrs->o_prop.copyTo(nrs->prop);
  if(nrs->Nscalar) nrs->cds->o_prop.copyTo(nrs->cds->prop);

  nek::ocopyToNek(startTime, 0);

  // setup elliptic solvers

  const int nbrBIDs = bcMap::size(0);
  int NBCType = nbrBIDs + 1;

  if(nrs->Nscalar) {
    cds_t* cds = nrs->cds;

    for (int is = 0; is < cds->NSfields; is++) {
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(2) << is;
      std::string sid = ss.str();
 
      if(!cds->compute[is]) continue;
 
      mesh_t* mesh;
      (is) ? mesh = cds->meshV : mesh = cds->mesh[0]; // only first scalar can be a CHT mesh

      if (platform->comm.mpiRank == 0)
        std::cout << "================= ELLIPTIC SETUP SCALAR" << sid << " ===============\n";

      int nbrBIDs = bcMap::size(0);
      if(nrs->cht && is == 0) nbrBIDs = bcMap::size(1);
      int* sBCType = (int*) calloc(nbrBIDs + 1, sizeof(int));
 
      for (int bID = 1; bID <= nbrBIDs; bID++) {
        std::string bcTypeText(bcMap::text(bID, "scalar" + sid));
        if(platform->comm.mpiRank == 0) printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str());
        sBCType[bID] = bcMap::type(bID, "scalar" + sid);
      }
 
      cds->solver[is] = new elliptic_t();
      cds->solver[is]->name = "scalar" + sid;
      cds->solver[is]->blockSolver = 0;
      cds->solver[is]->Nfields = 1;
      cds->solver[is]->Ntotal = nrs->fieldOffset;
      cds->solver[is]->o_wrk = o_mempoolElliptic;
      cds->solver[is]->mesh = mesh;
      cds->solver[is]->dim = cds->dim;
      cds->solver[is]->elementType = cds->elementType;
      cds->solver[is]->BCType = (int*) calloc(nbrBIDs + 1,sizeof(int));
      memcpy(cds->solver[is]->BCType,sBCType,(nbrBIDs + 1) * sizeof(int));
      free(sBCType);

      const int coeffField = platform->options.compareArgs("SCALAR" + sid + " COEFF FIELD", "TRUE");
      cds->solver[is]->coeffField = coeffField;
      cds->solver[is]->coeffFieldPreco = coeffField;
      cds->solver[is]->poisson = 0;

      platform->linAlg->fill(2*nrs->fieldOffset, 1.0, nrs->o_ellipticCoeff);
      cds->solver[is]->o_lambda = cds->o_ellipticCoeff;
      cds->solver[is]->loffset = 0;
 
      cds->solver[is]->options = cds->options[is];
      ellipticSolveSetup(cds->solver[is]);
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
      std::string bcTypeText(bcMap::text(bID, "velocity"));
      if(platform->comm.mpiRank == 0) printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str());

      uBCType[bID] = bcMap::type(bID, "x-velocity");
      vBCType[bID] = bcMap::type(bID, "y-velocity");
      wBCType[bID] = bcMap::type(bID, "z-velocity");
    }

    nrs->vOptions = options;
    nrs->vOptions.setArgs("PGMRES RESTART",        options.getArgs("VELOCITY PGMRES RESTART"));
    nrs->vOptions.setArgs("KRYLOV SOLVER",        options.getArgs("VELOCITY KRYLOV SOLVER"));
    nrs->vOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("VELOCITY SOLVER TOLERANCE"));
    nrs->vOptions.setArgs("LINEAR SOLVER STOPPING CRITERION",     options.getArgs("VELOCITY LINEAR SOLVER STOPPING CRITERION"));
    nrs->vOptions.setArgs("DISCRETIZATION",       options.getArgs("VELOCITY DISCRETIZATION"));
    nrs->vOptions.setArgs("BASIS",                options.getArgs("VELOCITY BASIS"));
    nrs->vOptions.setArgs("PRECONDITIONER",       options.getArgs("VELOCITY PRECONDITIONER"));
    nrs->vOptions.setArgs("INITIAL GUESS",       options.getArgs("VELOCITY INITIAL GUESS"));
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
    nrs->vOptions.setArgs("MAXIMUM ITERATIONS", options.getArgs("VELOCITY MAXIMUM ITERATIONS"));
    nrs->vOptions.setArgs("STABILIZATION METHOD", options.getArgs("VELOCITY STABILIZATION METHOD"));
    nrs->vOptions.setArgs("HPFRT STRENGTH", options.getArgs("VELOCITY HPFRT STRENGTH"));
    nrs->vOptions.setArgs("HPFRT MODES", options.getArgs("VELOCITY HPFRT MODES"));

    nrs->mOptions = options;
    nrs->mOptions.setArgs("PGMRES RESTART",        options.getArgs("MESH PGMRES RESTART"));
    nrs->mOptions.setArgs("KRYLOV SOLVER",        options.getArgs("MESH KRYLOV SOLVER"));
    nrs->mOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("MESH SOLVER TOLERANCE"));
    nrs->mOptions.setArgs("DISCRETIZATION",       options.getArgs("MESH DISCRETIZATION"));
    nrs->mOptions.setArgs("BASIS",                options.getArgs("MESH BASIS"));
    nrs->mOptions.setArgs("PRECONDITIONER",       options.getArgs("MESH PRECONDITIONER"));
    nrs->mOptions.setArgs("INITIAL GUESS",       options.getArgs("MESH INITIAL GUESS"));
    nrs->mOptions.setArgs("RESIDUAL PROJECTION VECTORS",       options.getArgs("MESH RESIDUAL PROJECTION VECTORS"));
    nrs->mOptions.setArgs("RESIDUAL PROJECTION START",       options.getArgs("MESH RESIDUAL PROJECTION START"));
    nrs->mOptions.setArgs("MULTIGRID COARSENING", options.getArgs("MESH MULTIGRID COARSENING"));
    nrs->mOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("MESH MULTIGRID SMOOTHER"));
    nrs->mOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE",
                          options.getArgs("MESH MULTIGRID CHEBYSHEV DEGREE"));
    nrs->mOptions.setArgs("PARALMOND CYCLE",      options.getArgs("MESH PARALMOND CYCLE"));
    nrs->mOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("MESH PARALMOND SMOOTHER"));
    nrs->mOptions.setArgs("PARALMOND PARTITION",  options.getArgs("MESH PARALMOND PARTITION"));
    nrs->mOptions.setArgs("PARALMOND CHEBYSHEV DEGREE",
                          options.getArgs("MESH PARALMOND CHEBYSHEV DEGREE"));
    nrs->mOptions.setArgs("PARALMOND AGGREGATION STRATEGY",
                          options.getArgs("MESH PARALMOND AGGREGATION STRATEGY"));
    nrs->mOptions.setArgs("MAXIMUM ITERATIONS", options.getArgs("MESH MAXIMUM ITERATIONS"));

    // coeff used by ellipticSetup to detect allNeumann
    platform->linAlg->fill(2*nrs->fieldOffset, 1.0, nrs->o_ellipticCoeff);
    
    const int velCoeffField = platform->options.compareArgs("VELOCITY COEFF FIELD", "TRUE");

    if(nrs->uvwSolver) {
      nrs->uvwSolver->blockSolver = 1;
      nrs->uvwSolver->stressForm = 0;
      if(options.compareArgs("STRESSFORMULATION", "TRUE"))
        nrs->uvwSolver->stressForm = 1;
      nrs->uvwSolver->Nfields = nrs->NVfields;
      nrs->uvwSolver->Ntotal = nrs->fieldOffset;
      nrs->uvwSolver->o_wrk = o_mempoolElliptic;
      nrs->uvwSolver->mesh = mesh;
      nrs->uvwSolver->options = nrs->vOptions;
      nrs->uvwSolver->dim = nrs->dim;
      nrs->uvwSolver->elementType = nrs->elementType;
      nrs->uvwSolver->NBCType = NBCType;
      nrs->uvwSolver->BCType = (int*) calloc(nrs->NVfields * NBCType,sizeof(int));
      memcpy(nrs->uvwSolver->BCType,uvwBCType,nrs->NVfields * NBCType * sizeof(int));
      nrs->uvwSolver->coeffField = velCoeffField;
      nrs->uvwSolver->coeffFieldPreco = velCoeffField;
      nrs->uvwSolver->o_lambda = nrs->o_ellipticCoeff;
      nrs->uvwSolver->loffset = 0; // use same ellipticCoeff for u,v and w
      nrs->uvwSolver->poisson = 0;

      ellipticSolveSetup(nrs->uvwSolver);
    } else {
      nrs->uSolver = new elliptic_t();
      nrs->uSolver->blockSolver = 0;
      nrs->uSolver->Nfields = 1;
      nrs->uSolver->Ntotal = nrs->fieldOffset;
      nrs->uSolver->o_wrk = o_mempoolElliptic;
      nrs->uSolver->mesh = mesh;
      nrs->uSolver->options = nrs->vOptions;
      nrs->uSolver->dim = nrs->dim;
      nrs->uSolver->elementType = nrs->elementType;
      nrs->uSolver->NBCType = NBCType;
      nrs->uSolver->BCType = (int*) calloc(NBCType,sizeof(int));
      memcpy(nrs->uSolver->BCType,uBCType,NBCType * sizeof(int));
      nrs->uSolver->coeffField = velCoeffField;
      nrs->uSolver->coeffFieldPreco = velCoeffField;
      nrs->uSolver->o_lambda = nrs->o_ellipticCoeff;
      nrs->uSolver->loffset = 0;
      nrs->uSolver->poisson = 0;

      ellipticSolveSetup(nrs->uSolver);

      nrs->vSolver = new elliptic_t();
      nrs->vSolver->blockSolver = 0;
      nrs->vSolver->Nfields = 1;
      nrs->vSolver->Ntotal = nrs->fieldOffset;
      nrs->vSolver->o_wrk = o_mempoolElliptic;
      nrs->vSolver->mesh = mesh;
      nrs->vSolver->options = nrs->vOptions;
      nrs->vSolver->dim = nrs->dim;
      nrs->vSolver->elementType = nrs->elementType;
      nrs->vSolver->NBCType = NBCType;
      nrs->vSolver->BCType = (int*) calloc(NBCType,sizeof(int));
      memcpy(nrs->vSolver->BCType,vBCType,NBCType * sizeof(int));
      nrs->vSolver->coeffField = velCoeffField;
      nrs->vSolver->coeffFieldPreco = velCoeffField;
      nrs->vSolver->o_lambda = nrs->o_ellipticCoeff;
      nrs->vSolver->loffset = 0;
      nrs->vSolver->poisson = 0;

      ellipticSolveSetup(nrs->vSolver);

      if (nrs->dim == 3) {
        nrs->wSolver = new elliptic_t();
        nrs->wSolver->blockSolver = 0;
        nrs->wSolver->Nfields = 1;
        nrs->wSolver->Ntotal = nrs->fieldOffset;
        nrs->wSolver->o_wrk = o_mempoolElliptic;
        nrs->wSolver->mesh = mesh;
        nrs->wSolver->options = nrs->vOptions;
        nrs->wSolver->dim = nrs->dim;
        nrs->wSolver->elementType = nrs->elementType;
        nrs->wSolver->NBCType = NBCType;
        nrs->wSolver->BCType = (int*) calloc(NBCType,sizeof(int));
        memcpy(nrs->wSolver->BCType,wBCType,NBCType * sizeof(int));
        nrs->wSolver->coeffField = velCoeffField;
        nrs->wSolver->coeffFieldPreco = velCoeffField;
        nrs->wSolver->o_lambda = nrs->o_ellipticCoeff;
        nrs->wSolver->loffset = 0;
        nrs->wSolver->poisson = 0;

        ellipticSolveSetup(nrs->wSolver);
      }
    }

    if(platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE")) {
      nrs->uvwSolver->name = "velocity";
    } else {
      nrs->uSolver->name = "x-velocity";
      nrs->vSolver->name = "y-velocity";
      nrs->wSolver->name = "z-velocity";
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
    nrs->pOptions.setArgs("LINEAR SOLVER STOPPING CRITERION",     options.getArgs("PRESSURE LINEAR SOLVER STOPPING CRITERION"));
    nrs->pOptions.setArgs("DISCRETIZATION",       options.getArgs("PRESSURE DISCRETIZATION"));
    nrs->pOptions.setArgs("BASIS",                options.getArgs("PRESSURE BASIS"));
    nrs->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRESSURE PRECONDITIONER"));
    nrs->pOptions.setArgs("SEMFEM SOLVER", options.getArgs("PRESSURE SEMFEM SOLVER"));
    nrs->pOptions.setArgs("SEMFEM SOLVER PRECISION", options.getArgs("PRESSURE SEMFEM SOLVER PRECISION"));
    nrs->pOptions.setArgs("MULTIGRID COARSENING", options.getArgs("PRESSURE MULTIGRID COARSENING"));
    nrs->pOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
    nrs->pOptions.setArgs("MULTIGRID COARSE SOLVE",   options.getArgs("PRESSURE MULTIGRID COARSE SOLVE"));
    nrs->pOptions.setArgs("MULTIGRID COARSE SEMFEM",   options.getArgs("PRESSURE MULTIGRID COARSE SEMFEM"));
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
    nrs->pOptions.setArgs("INITIAL GUESS", options.getArgs("PRESSURE INITIAL GUESS"));
    nrs->pOptions.setArgs("RESIDUAL PROJECTION VECTORS",
                          options.getArgs("PRESSURE RESIDUAL PROJECTION VECTORS"));
    nrs->pOptions.setArgs("RESIDUAL PROJECTION START",
                          options.getArgs("PRESSURE RESIDUAL PROJECTION START"));
    nrs->pOptions.setArgs("MULTIGRID VARIABLE COEFFICIENT", "FALSE");
    nrs->pOptions.setArgs("MAXIMUM ITERATIONS", options.getArgs("PRESSURE MAXIMUM ITERATIONS"));
    nrs->pOptions.setArgs("MULTIGRID CHEBYSHEV MAX EIGENVALUE BOUND FACTOR", options.getArgs("PRESSURE MULTIGRID CHEBYSHEV MAX EIGENVALUE BOUND FACTOR"));
    nrs->pOptions.setArgs("MULTIGRID CHEBYSHEV MIN EIGENVALUE BOUND FACTOR", options.getArgs("PRESSURE MULTIGRID CHEBYSHEV MIN EIGENVALUE BOUND FACTOR"));

    nrs->pSolver = new elliptic_t();
    nrs->pSolver->name = "pressure";
    nrs->pSolver->blockSolver = 0;
    nrs->pSolver->Nfields = 1;
    nrs->pSolver->Ntotal = nrs->fieldOffset;
    nrs->pSolver->o_wrk = o_mempoolElliptic;
    nrs->pSolver->mesh = mesh;
    nrs->pSolver->dim = nrs->dim;
    nrs->pSolver->elementType = nrs->elementType;
    nrs->pSolver->BCType = (int*) calloc(nbrBIDs + 1,sizeof(int));
    memcpy(nrs->pSolver->BCType,pBCType,(nbrBIDs + 1) * sizeof(int));

    int pCoeffField = 0;

    if(platform->options.compareArgs("LOWMACH", "TRUE"))
      pCoeffField = 1;

    nrs->pSolver->coeffField = pCoeffField;
    nrs->pSolver->coeffFieldPreco = pCoeffField;
    nrs->pSolver->poisson = 1;

    // lambda0 = 1/rho
    // lambda1 = 0
    platform->linAlg->fill(2*nrs->fieldOffset, 0.0, nrs->o_ellipticCoeff);
    nrs->o_ellipticCoeff.copyFrom(nrs->o_rho, nrs->fieldOffset * sizeof(dfloat));
    platform->linAlg->ady(mesh->Nlocal, 1.0, nrs->o_ellipticCoeff);
    nrs->pSolver->o_lambda = nrs->o_ellipticCoeff;
    nrs->pSolver->loffset = 0; // Poisson

    {
      const std::vector<int> levels = determineMGLevels("pressure");
      nrs->pSolver->nLevels = levels.size();
      nrs->pSolver->levels = (int*) calloc(nrs->pSolver->nLevels,sizeof(int));
      for(int i = 0; i < nrs->pSolver->nLevels; ++i)
        nrs->pSolver->levels[i] = levels.at(i);
    }

    nrs->pSolver->options = nrs->pOptions;
    ellipticSolveSetup(nrs->pSolver);

  } // flow
  if(nrs->flow){
    if(options.compareArgs("MESH SOLVER", "ELASTICITY")){
      if (platform->comm.mpiRank == 0) printf("================ ELLIPTIC SETUP MESH ================\n");
      int* uvwMeshBCType = (int*) calloc(3 * NBCType, sizeof(int));
      int* uMeshBCType = uvwMeshBCType + 0 * NBCType;
      int* vMeshBCType = uvwMeshBCType + 1 * NBCType;
      int* wMeshBCType = uvwMeshBCType + 2 * NBCType;
      for (int bID = 1; bID <= nbrBIDs; bID++) {
        std::string bcTypeText(bcMap::text(bID, "mesh"));
        if(platform->comm.mpiRank == 0) printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str());

        uMeshBCType[bID] = bcMap::type(bID, "x-mesh");
        vMeshBCType[bID] = bcMap::type(bID, "y-mesh");
        wMeshBCType[bID] = bcMap::type(bID, "z-mesh");
      }

      const int meshCoeffField = platform->options.compareArgs("MESH COEFF FIELD", "TRUE");
      platform->linAlg->fill(2*nrs->fieldOffset, 1.0, nrs->o_ellipticCoeff);

      nrs->meshSolver = new elliptic_t();
      nrs->meshSolver->name = "mesh";
      nrs->meshSolver->blockSolver = 1;
      nrs->meshSolver->stressForm = 1;
      nrs->meshSolver->Nfields = nrs->NVfields;
      nrs->meshSolver->Ntotal = nrs->fieldOffset;
      nrs->meshSolver->o_wrk = o_mempoolElliptic;
      nrs->meshSolver->mesh = mesh;
      nrs->meshSolver->options = nrs->mOptions;
      nrs->meshSolver->dim = nrs->dim;
      nrs->meshSolver->elementType = nrs->elementType;
      nrs->meshSolver->NBCType = NBCType;
      nrs->meshSolver->BCType = (int*) calloc(nrs->NVfields * NBCType,sizeof(int));
      memcpy(nrs->meshSolver->BCType,uvwMeshBCType,nrs->NVfields * NBCType * sizeof(int));
      nrs->meshSolver->coeffField = meshCoeffField;
      nrs->meshSolver->coeffFieldPreco = meshCoeffField;
      nrs->meshSolver->o_lambda = nrs->o_ellipticCoeff;
      nrs->meshSolver->loffset = 0; // use same ellipticCoeff for u,v and w
      nrs->meshSolver->poisson = 0;

      ellipticSolveSetup(nrs->meshSolver);
    }
  }
  // set I.C. for U, W
  if(platform->options.compareArgs("MESH SOLVER", "ELASTICITY"))
  {
    double startTime;
    platform->options.getArgs("START TIME", startTime);
    platform->linAlg->fill(nrs->NVfields*nrs->fieldOffset, -1.0*std::numeric_limits<dfloat>::max(), platform->o_mempool.slice0);
    for (int sweep = 0; sweep < 2; sweep++) {
      nrs->velocityDirichletBCKernel(mesh->Nelements,
                                     nrs->fieldOffset,
                                     startTime,
                                     mesh->o_sgeo,
                                     mesh->o_x,
                                     mesh->o_y,
                                     mesh->o_z,
                                     mesh->o_vmapM,
                                     mesh->o_EToB,
                                     nrs->o_EToB,
                                     nrs->o_usrwrk,
                                     nrs->o_U,
                                     platform->o_mempool.slice0);
      if (sweep == 0) oogs::startFinish(platform->o_mempool.slice0, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsMax, nrs->gsh);
      if (sweep == 1) oogs::startFinish(platform->o_mempool.slice0, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsMin, nrs->gsh);
    }
    platform->o_mempool.slice3.copyFrom(platform->o_mempool.slice0, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));

    platform->linAlg->fill(nrs->NVfields*nrs->fieldOffset, 0.0, platform->o_mempool.slice0);
    for (int sweep = 0; sweep < 2; sweep++) {
    nrs->meshV->velocityDirichletKernel(mesh->Nelements,
                                   nrs->fieldOffset,
                                   mesh->o_vmapM,
                                   nrs->o_EToBMesh,
                                   platform->o_mempool.slice3,
                                   platform->o_mempool.slice0);
      //take care of Neumann-Dirichlet shared edges across elements
      if(sweep == 0) oogs::startFinish(platform->o_mempool.slice0, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsMax, nrs->gsh);
      if(sweep == 1) oogs::startFinish(platform->o_mempool.slice0, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsMin, nrs->gsh);
    }
    oogs::startFinish(platform->o_mempool.slice0, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);
    platform->linAlg->axmyMany(
      mesh->Nlocal,
      nrs->NVfields,
      nrs->fieldOffset,
      0,
      1.0,
      nrs->meshSolver->o_invDegree,
      platform->o_mempool.slice0
    );
    mesh->o_U.copyFrom(platform->o_mempool.slice0, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
  }


}

namespace{
cds_t* cdsSetup(nrs_t* nrs, setupAide options)
{
  const std::string section = "cds-";
  cds_t* cds = new cds_t();
  platform_t* platform = platform_t::getInstance();
  device_t& device = platform->device;

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
    (dfloat*) calloc(cds->nEXT * cds->fieldOffsetSum,sizeof(dfloat));

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
    std::string sid = ss.str();

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
    std::string sid = ss.str();

    cds->compute[is] = 1;
    if (options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
      cds->compute[is] = 0;
      continue;
    }

    mesh_t* mesh; 
    (is) ? mesh = cds->meshV : mesh = cds->mesh[0]; // only first scalar can be a CHT mesh
 
    cds->options[is] = options;

    cds->options[is].setArgs("REGULARIZATION RAMP CONSTANT", options.getArgs("SCALAR" + sid + " REGULARIZATION RAMP CONSTANT"));
    cds->options[is].setArgs("REGULARIZATION AVM C0", options.getArgs("SCALAR" + sid + " REGULARIZATION AVM C0"));
    cds->options[is].setArgs("REGULARIZATION METHOD", options.getArgs("SCALAR" + sid + " REGULARIZATION METHOD"));
    cds->options[is].setArgs("REGULARIZATION VISMAX COEFF", options.getArgs("SCALAR" + sid + " REGULARIZATION VISMAX COEFF"));
    cds->options[is].setArgs("REGULARIZATION SCALING COEFF", options.getArgs("SCALAR" + sid + " REGULARIZATION SCALING COEFF"));
    cds->options[is].setArgs("HPFRT STRENGTH", options.getArgs("SCALAR" + sid + " HPFRT STRENGTH"));
    cds->options[is].setArgs("HPFRT MODES", options.getArgs("SCALAR" + sid + " HPFRT MODES"));
    cds->options[is].setArgs("KRYLOV SOLVER", options.getArgs("SCALAR" + sid + " KRYLOV SOLVER"));
    cds->options[is].setArgs("PGMRES RESTART", options.getArgs("SCALAR" + sid + " PGMRES RESTART"));
    cds->options[is].setArgs("DISCRETIZATION", options.getArgs("SCALAR DISCRETIZATION"));
    cds->options[is].setArgs("BASIS", options.getArgs("SCALAR BASIS"));
    cds->options[is].setArgs("PRECONDITIONER", options.getArgs("SCALAR" + sid + " PRECONDITIONER"));
    cds->options[is].setArgs("SOLVER TOLERANCE",
                         options.getArgs("SCALAR" + sid +  " SOLVER TOLERANCE"));
    cds->options[is].setArgs("LINEAR SOLVER STOPPING CRITERION",     options.getArgs("SCALAR" + sid + " LINEAR SOLVER STOPPING CRITERION"));
    cds->options[is].setArgs("INITIAL GUESS",  options.getArgs("SCALAR" + sid + " INITIAL GUESS"));
    cds->options[is].setArgs("RESIDUAL PROJECTION VECTORS",  options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION VECTORS"));
    cds->options[is].setArgs("RESIDUAL PROJECTION START",  options.getArgs("SCALAR" + sid + " RESIDUAL PROJECTION START"));
    cds->options[is].setArgs("MAXIMUM ITERATIONS", options.getArgs("SCALAR" + sid + " MAXIMUM ITERATIONS"));

    dfloat largeNumber = 1 << 20;
    cds->EToB[is] = (int*) calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
    int* EToB = cds->EToB[is];
    int cnt = 0;
    for (int e = 0; e < mesh->Nelements; e++) {
      for (int f = 0; f < mesh->Nfaces; f++) {
        int bc = bcMap::id(mesh->EToB[f + e * mesh->Nfaces], "scalar" + sid);
        EToB[cnt] = bc;
        cnt++;
      }
    }
    cds->o_EToB[is] = device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int), EToB);
  }

  bool scalarFilteringEnabled = false;
  bool avmEnabled = false;
  {
    for(int is = 0; is < cds->NSfields; is++) {
      if(!cds->options[is].compareArgs("REGULARIZATION METHOD", "NONE")) scalarFilteringEnabled = true;
      if(cds->options[is].compareArgs("REGULARIZATION METHOD", "HPF_RESIDUAL")) avmEnabled = true;
      if(cds->options[is].compareArgs("REGULARIZATION METHOD", "HIGHEST_MODAL_DECAY")) avmEnabled = true;
    }
  }

  if(scalarFilteringEnabled)
  {
    const dlong Nmodes = cds->mesh[0]->N + 1;
    cds->o_filterMT = platform->device.malloc(cds->NSfields * Nmodes * Nmodes, sizeof(dfloat));
    for(int is = 0; is < cds->NSfields; is++)
    {
      if(cds->options[is].compareArgs("REGULARIZATION METHOD", "NONE")) continue;
      if(!cds->compute[is]) continue;
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

  std::string kernelName;
  const std::string suffix = "Hex3D";

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("loading cds kernels ... "); fflush(stdout);

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
          kernelName = "subCycleStrongCubatureVolume" + suffix;
          cds->subCycleStrongCubatureVolumeKernel = platform->kernels.get(section + kernelName);
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

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel); fflush(stdout);

  return cds;
}
}
