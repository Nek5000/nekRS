#include "nrs.hpp"
#include "meshSetup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "filter.hpp"
#include "bcMap.hpp"

static dfloat *scratch;
static occa::memory o_scratch;

cds_t *cdsSetup(ins_t *ins, mesh_t *mesh, setupAide options, occa::properties &kernelInfoH);
              
 
ins_t *insSetup(MPI_Comm comm, setupAide &options, int buildOnly)
{

  ins_t *ins = new ins_t();
  ins->options = options;
  ins->kernelInfo = new occa::properties();
  occa::properties& kernelInfo = *ins->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  kernelInfo["include_paths"].asArray();

  int N;
  string install_dir;
  options.getArgs("POLYNOMIAL DEGREE", N);
  options.getArgs("NUMBER OF SCALARS", ins->Nscalar);
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  options.getArgs("MESH DIMENSION", ins->dim);
  options.getArgs("ELEMENT TYPE", ins->elementType);
 
  ins->flow = 1;
  if(options.compareArgs("VELOCITY", "FALSE")) ins->flow = 0;
  if(options.compareArgs("VELOCITY SOLVER", "NONE")) ins->flow = 0;

  ins->cht = 0;
  if (nekData.nelv != nekData.nelt && ins->Nscalar) ins->cht = 1;

  if (buildOnly) {
    ins->meshT = createMeshDummy(comm, N, options, kernelInfo);
    ins->mesh = ins->meshT;
  } else {
    ins->meshT = createMeshT(comm, N, ins->cht, options, kernelInfo);
    ins->mesh = ins->meshT;
    if (ins->cht) ins->mesh = createMeshV(comm, N, ins->meshT, options, kernelInfo);
  }
  mesh_t *mesh = ins->mesh;

  occa::properties kernelInfoV  = kernelInfo;
  occa::properties kernelInfoP  = kernelInfo;
  occa::properties kernelInfoS  = kernelInfo;

  ins->NVfields = (ins->dim==3) ? 3:2; // Total Number of Velocity Fields
  ins->NTfields = ins->NVfields + 1;   // Total Velocity + Pressure

  ins->SNrk = 0;
  options.getArgs("SUBCYCLING TIME STAGE NUMBER", ins->SNrk);

  mesh->Nfields = 1; 
  ins->g0 =  1.0;

  ins->extbdfA = (dfloat*) calloc(3, sizeof(dfloat));
  ins->extbdfB = (dfloat*) calloc(3, sizeof(dfloat));
  ins->extbdfC = (dfloat*) calloc(3, sizeof(dfloat));
 
  ins->extC = (dfloat*) calloc(3, sizeof(dfloat));

  if (options.compareArgs("TIME INTEGRATOR", "TOMBO1")) {
    ins->Nstages = 1;
    ins->temporalOrder = 1;
    ins->g0 = 1.0;
  } else if (options.compareArgs("TIME INTEGRATOR", "TOMBO2")) { 
    ins->Nstages = 2;
    ins->temporalOrder = 2;
    ins->g0 = 1.5;
  } else if (options.compareArgs("TIME INTEGRATOR", "TOMBO3")) { 
    ins->Nstages = 3;
    ins->temporalOrder = 3;
    ins->g0 = 11.f/6.f;
  }

  ins->readRestartFile = 0; 
  options.getArgs("RESTART FROM FILE", ins->readRestartFile);
  
  ins->writeRestartFile = 0; 
  options.getArgs("WRITE RESTART FILE", ins->writeRestartFile);

  dfloat mue = 1;
  dfloat rho = 1;
  options.getArgs("VISCOSITY", mue);
  options.getArgs("DENSITY", rho);

  options.getArgs("SUBCYCLING STEPS",ins->Nsubsteps);

  dfloat dt; 
  options.getArgs("DT", dt);
  ins->dt = dt; 

  options.getArgs("FINAL TIME", ins->finalTime);
  options.getArgs("START TIME", ins->startTime);
  if(ins->startTime > 0.0) {
    int numSteps;
    if(options.getArgs("NUMBER TIMESTEPS", numSteps)) 
      ins->finalTime += ins->startTime; 
  }
 
  ins->NtimeSteps = ceil((ins->finalTime-ins->startTime)/ins->dt);
  options.setArgs("NUMBER TIMESTEPS", std::to_string(ins->NtimeSteps)); 
  if(ins->Nsubsteps) ins->sdt = ins->dt/ins->Nsubsteps;

  // Hold some inverses for kernels
  ins->idt = 1.0/ins->dt;
  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", ins->outputStep);

  const dlong Nlocal = mesh->Np*mesh->Nelements;
  const dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

  ins->Nlocal = Nlocal;
  ins->Ntotal = Ntotal;

  // ensure that offset is large enough for v and t mesh and is properly aligned 
  {
    const dlong NtotalT = ins->meshT->Np*(ins->meshT->Nelements+ins->meshT->totalHaloPairs);
    ins->fieldOffset = mymax(Ntotal, NtotalT);

    int PAGESIZE = 4096; // default is 4kB 
    char *tmp;
    tmp = getenv("NEKRS_PAGE_SIZE");
    if (tmp != NULL) PAGESIZE = std::stoi(tmp);
    const int pageW = PAGESIZE/sizeof(dfloat);
    if (ins->fieldOffset%pageW) ins->fieldOffset = (ins->fieldOffset/pageW + 1)*pageW;
  }

  ins->Nblock = (Nlocal+blockSize-1)/blockSize;

  ins->U  = (dfloat*) calloc(ins->NVfields*ins->Nstages*ins->fieldOffset,sizeof(dfloat));
  ins->Ue = (dfloat*) calloc(ins->NVfields*ins->fieldOffset,sizeof(dfloat));

  ins->P  = (dfloat*) calloc(ins->fieldOffset,sizeof(dfloat));
  ins->PI = (dfloat*) calloc(ins->fieldOffset,sizeof(dfloat));

  ins->BF = (dfloat*) calloc(ins->NVfields*ins->fieldOffset,sizeof(dfloat));
  ins->FU = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*ins->fieldOffset,sizeof(dfloat));

  if(ins->Nsubsteps){
    int Sorder; 
    options.getArgs("SUBCYCLING TIME ORDER", Sorder);
    if(Sorder==4 && ins->SNrk==4){ // ERK(4,4)
      dfloat rka[4] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
      dfloat rkb[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
      dfloat rkc[4] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, ins->SNrk*sizeof(dfloat));
    }else{
      if(mesh->rank==0) cout << "Unsupported subcycling scheme!\n"; 
      ABORT(1);
    }
    ins->o_Srka = mesh->device.malloc(ins->SNrk*sizeof(dfloat), ins->Srka);
    ins->o_Srkb = mesh->device.malloc(ins->SNrk*sizeof(dfloat), ins->Srkb);
  }

  // setup scratch space
  const int wrkNflds = 6; 
  const int ellipticWrkNflds = 15;
  ins->ellipticWrkOffset = wrkNflds*ins->fieldOffset;

  const int scratchNflds = wrkNflds+ellipticWrkNflds;
  scratch   = (dfloat*) calloc(scratchNflds*ins->fieldOffset,sizeof(dfloat));
  o_scratch = mesh->device.malloc(scratchNflds*ins->fieldOffset*sizeof(dfloat), scratch);

  ins->o_wrk0  = o_scratch.slice( 0*ins->fieldOffset*sizeof(dfloat));
  ins->o_wrk1  = o_scratch.slice( 1*ins->fieldOffset*sizeof(dfloat));
  ins->o_wrk2  = o_scratch.slice( 2*ins->fieldOffset*sizeof(dfloat));
  ins->o_wrk3  = o_scratch.slice( 3*ins->fieldOffset*sizeof(dfloat));
  ins->o_wrk4  = o_scratch.slice( 4*ins->fieldOffset*sizeof(dfloat));
  ins->o_wrk5  = o_scratch.slice( 5*ins->fieldOffset*sizeof(dfloat));
  ins->o_wrk6  = o_scratch.slice( 6*ins->fieldOffset*sizeof(dfloat));
  ins->o_wrk9  = o_scratch.slice( 9*ins->fieldOffset*sizeof(dfloat));
  ins->o_wrk12 = o_scratch.slice(12*ins->fieldOffset*sizeof(dfloat));
  ins->o_wrk15 = o_scratch.slice(15*ins->fieldOffset*sizeof(dfloat));

  // dummy decleration for user work space 
  //ins->usrwrk   = (dfloat*) calloc(1, sizeof(dfloat));
  //ins->o_usrwrk = mesh->device.malloc(1*sizeof(dfloat), ins->usrwrk);

  ins->o_U  = mesh->device.malloc(ins->NVfields*ins->Nstages*ins->fieldOffset*sizeof(dfloat), ins->U);
  ins->o_Ue = mesh->device.malloc(ins->NVfields*ins->fieldOffset*sizeof(dfloat), ins->Ue);
  ins->o_P  = mesh->device.malloc(ins->fieldOffset*sizeof(dfloat), ins->P);
  ins->o_PI = mesh->device.malloc(ins->fieldOffset*sizeof(dfloat), ins->PI);

  ins->o_FU = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*ins->fieldOffset*sizeof(dfloat), ins->FU);
  ins->o_BF = mesh->device.malloc(ins->NVfields*ins->fieldOffset*sizeof(dfloat), ins->BF);

  ins->var_coeff = 1; // use always var coeff elliptic
  ins->ellipticCoeff = (dfloat*) calloc(2*ins->fieldOffset,sizeof(dfloat));
  ins->o_ellipticCoeff = mesh->device.malloc(2*ins->fieldOffset*sizeof(dfloat), ins->ellipticCoeff);  

  ins->prop =  (dfloat*) calloc(2*ins->fieldOffset,sizeof(dfloat));
  for (int e=0;e<mesh->Nelements;e++) { 
    for (int n=0;n<mesh->Np;n++) {
      ins->prop[0*ins->fieldOffset + e*mesh->Np + n] = mue;
      ins->prop[1*ins->fieldOffset + e*mesh->Np + n] = rho;
    }
  }
  ins->o_prop = mesh->device.malloc(2*ins->fieldOffset*sizeof(dfloat), ins->prop);  
  ins->o_mue = ins->o_prop.slice(0*ins->fieldOffset*sizeof(dfloat));
  ins->o_rho = ins->o_prop.slice(1*ins->fieldOffset*sizeof(dfloat));

  ins->lowMach = 0;
  if(options.compareArgs("LOWMACH", "TRUE")) ins->lowMach = 1;
  ins->div   = (dfloat*) calloc(ins->fieldOffset,sizeof(dfloat));
  ins->o_div = mesh->device.malloc(ins->fieldOffset*sizeof(dfloat), ins->div);  

  ins->elementInfo = (dlong*) calloc(ins->meshT->Nelements,sizeof(dlong));
  for (int e=0;e<ins->meshT->Nelements;e++) ins->elementInfo[e] = mesh->elementInfo[e]; 
  ins->o_elementInfo = mesh->device.malloc(ins->meshT->Nelements*sizeof(dlong), ins->elementInfo);  
  dfloat rkC[4]  = {1.0, 0.0, -1.0, -2.0};
  ins->o_rkC     = mesh->device.malloc(4*sizeof(dfloat),rkC);
  ins->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
  ins->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
  ins->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 
  ins->o_extC    = mesh->device.malloc(3*sizeof(dfloat)); 
  ins->o_prkA    = ins->o_extbdfC;
  ins->o_prkB    = ins->o_extbdfC;

  dfloat *lumpedMassMatrix  = (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  ins->o_InvM = 
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), lumpedMassMatrix);

  // define aux kernel constants
  kernelInfo["defines/" "p_eNfields"] = ins->NVfields;
  kernelInfo["defines/" "p_NTfields"]= ins->NTfields;
  kernelInfo["defines/" "p_NVfields"]= ins->NVfields;
  kernelInfo["defines/" "p_NfacesNfp"]=  mesh->Nfaces*mesh->Nfp;
  kernelInfo["defines/" "p_Nstages"]=  ins->Nstages;
  if(ins->Nsubsteps)
    kernelInfo["defines/" "p_SUBCYCLING"]=  1;
  else
    kernelInfo["defines/" "p_SUBCYCLING"]=  0;

  kernelInfo["defines/" "p_blockSize"]= blockSize;
  //kernelInfo["parser/" "automate-add-barriers"] =  "disabled";
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = mymax(1,256/mesh->Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1,256/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo["defines/" "p_maxNodesVolumeCub"]= maxNodesVolumeCub;
  int cubNblockV = mymax(1,256/maxNodesVolumeCub);

  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, 
                           mesh->Nfaces*mesh->intNfp));
  kernelInfo["defines/" "p_maxNodesSurfaceCub"]=maxNodesSurfaceCub;
  int cubNblockS = mymax(256/maxNodesSurfaceCub,1);
 
  kernelInfo["defines/" "p_cubNblockV"]=cubNblockV;
  kernelInfo["defines/" "p_cubNblockS"]=cubNblockS;

  // jit compile udf kernels
  if (udf.loadKernels) {
    if (mesh->rank == 0) cout << "building udf kernels ...";
    udf.loadKernels(ins);
    if (mesh->rank == 0) cout << " done" << endl;
  }
 
  occa::properties kernelInfoBC = kernelInfo;
  const string bcDataFile = install_dir + "/include/insBcData.h";
  kernelInfoBC["includes"] += bcDataFile.c_str();
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfoBC["includes"] += realpath(boundaryHeaderFileName.c_str(), NULL);

  if(options.compareArgs("FILTER STABILIZATION", "RELAXATION")) 
    filterSetup(ins); 

  const int nbrBIDs = bcMap::size(0);
  int NBCType = nbrBIDs+1;

  if (ins->flow) {

  if (mesh->rank==0) printf("==================VELOCITY SETUP=========================\n");

  ins->velTOL  = 1E-6;
  ins->uvwSolver = NULL;
  if(options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE"))
    ins->uvwSolver = new elliptic_t();

  int *uvwBCType = (int*) calloc(3*NBCType, sizeof(int));
  int *uBCType = uvwBCType + 0*NBCType;
  int *vBCType = uvwBCType + 1*NBCType;
  int *wBCType = uvwBCType + 2*NBCType;
  for (int bID=1; bID <= nbrBIDs; bID++) {
    string bcTypeText(bcMap::text(bID, "velocity"));
    if(mesh->rank == 0) printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str()); 

    uBCType[bID] = bcMap::type(bID, "x-velocity");
    vBCType[bID] = bcMap::type(bID, "y-velocity");
    wBCType[bID] = bcMap::type(bID, "z-velocity");
  }
  
  ins->vOptions = options;
  ins->vOptions.setArgs("KRYLOV SOLVER",        options.getArgs("VELOCITY KRYLOV SOLVER"));
  ins->vOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("VELOCITY SOLVER TOLERANCE"));
  ins->vOptions.setArgs("DISCRETIZATION",       options.getArgs("VELOCITY DISCRETIZATION"));
  ins->vOptions.setArgs("BASIS",                options.getArgs("VELOCITY BASIS"));
  ins->vOptions.setArgs("PRECONDITIONER",       options.getArgs("VELOCITY PRECONDITIONER"));
  ins->vOptions.setArgs("MULTIGRID COARSENING", options.getArgs("VELOCITY MULTIGRID COARSENING"));
  ins->vOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("VELOCITY MULTIGRID SMOOTHER"));
  ins->vOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE",  options.getArgs("VELOCITY MULTIGRID CHEBYSHEV DEGREE"));
  ins->vOptions.setArgs("PARALMOND CYCLE",      options.getArgs("VELOCITY PARALMOND CYCLE"));
  ins->vOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("VELOCITY PARALMOND SMOOTHER"));
  ins->vOptions.setArgs("PARALMOND PARTITION",  options.getArgs("VELOCITY PARALMOND PARTITION"));
  ins->vOptions.setArgs("PARALMOND CHEBYSHEV DEGREE",  options.getArgs("VELOCITY PARALMOND CHEBYSHEV DEGREE"));
  ins->vOptions.setArgs("PARALMOND AGGREGATION STRATEGY", options.getArgs("VELOCITY PARALMOND AGGREGATION STRATEGY"));
  ins->vOptions.setArgs("DEBUG ENABLE OGS", "1");
  ins->vOptions.setArgs("DEBUG ENABLE REDUCTIONS", "1");

  // coeff used by ellipticSetup to detect allNeumann 
  for (int i=0;i<2*ins->fieldOffset;i++) ins->ellipticCoeff[i] = 1; 

  if(ins->uvwSolver){
    ins->uvwSolver->blockSolver = 1;
    ins->uvwSolver->Nfields = ins->NVfields; 
    ins->uvwSolver->Ntotal = ins->fieldOffset;
    ins->uvwSolver->wrk = scratch + ins->ellipticWrkOffset; 
    ins->uvwSolver->o_wrk = o_scratch.slice(ins->ellipticWrkOffset*sizeof(dfloat)); 
    ins->uvwSolver->mesh = mesh;
    ins->uvwSolver->options = ins->vOptions;
    ins->uvwSolver->dim = ins->dim;
    ins->uvwSolver->elementType = ins->elementType;
    ins->uvwSolver->NBCType = NBCType;
    ins->uvwSolver->BCType = (int*) calloc(ins->NVfields*NBCType,sizeof(int));
    memcpy(ins->uvwSolver->BCType,uvwBCType,ins->NVfields*NBCType*sizeof(int));
    ins->uvwSolver->var_coeff = ins->var_coeff;
    ins->uvwSolver->lambda = ins->ellipticCoeff; 
    ins->uvwSolver->o_lambda = ins->o_ellipticCoeff; 
    ins->uvwSolver->loffset = 0; // use same ellipticCoeff for u,v and w
  
    ellipticSolveSetup(ins->uvwSolver, kernelInfoV); 
  } else {
    ins->uSolver = new elliptic_t();
    ins->uSolver->blockSolver = 0;
    ins->uSolver->Nfields = 1; 
    ins->uSolver->Ntotal = ins->fieldOffset;
    ins->uSolver->wrk = scratch + ins->ellipticWrkOffset; 
    ins->uSolver->o_wrk = o_scratch.slice(ins->ellipticWrkOffset*sizeof(dfloat)); 
    ins->uSolver->mesh = mesh;
    ins->uSolver->options = ins->vOptions;
    ins->uSolver->dim = ins->dim;
    ins->uSolver->elementType = ins->elementType;
    ins->uSolver->NBCType = NBCType;
    ins->uSolver->BCType = (int*) calloc(NBCType,sizeof(int));
    memcpy(ins->uSolver->BCType,uBCType,NBCType*sizeof(int));
    ins->uSolver->var_coeff = ins->var_coeff;
    ins->uSolver->lambda = ins->ellipticCoeff; 
    ins->uSolver->o_lambda = ins->o_ellipticCoeff; 
    ins->uSolver->loffset = 0;
 
    ellipticSolveSetup(ins->uSolver, kernelInfoV); 
  
    ins->vSolver = new elliptic_t();
    ins->vSolver->blockSolver = 0;
    ins->vSolver->Nfields = 1; 
    ins->vSolver->Ntotal = ins->fieldOffset;
    ins->vSolver->wrk = scratch + ins->ellipticWrkOffset; 
    ins->vSolver->o_wrk = o_scratch.slice(ins->ellipticWrkOffset*sizeof(dfloat)); 
    ins->vSolver->mesh = mesh;
    ins->vSolver->options = ins->vOptions;
    ins->vSolver->dim = ins->dim;
    ins->vSolver->elementType = ins->elementType;
    ins->vSolver->NBCType = NBCType;
    ins->vSolver->BCType = (int*) calloc(NBCType,sizeof(int));
    memcpy(ins->vSolver->BCType,vBCType,NBCType*sizeof(int));
    ins->vSolver->var_coeff = ins->var_coeff;
    ins->vSolver->lambda = ins->ellipticCoeff; 
    ins->vSolver->o_lambda = ins->o_ellipticCoeff; 
    ins->vSolver->loffset = 0;
    
    ellipticSolveSetup(ins->vSolver, kernelInfoV);
  
    if (ins->dim==3) {
      ins->wSolver = new elliptic_t();
      ins->wSolver->blockSolver = 0;
      ins->wSolver->Nfields = 1; 
      ins->wSolver->Ntotal = ins->fieldOffset;
      ins->wSolver->wrk = scratch + ins->ellipticWrkOffset; 
      ins->wSolver->o_wrk = o_scratch.slice(ins->ellipticWrkOffset*sizeof(dfloat)); 
      ins->wSolver->mesh = mesh;
      ins->wSolver->options = ins->vOptions;
      ins->wSolver->dim = ins->dim;
      ins->wSolver->elementType = ins->elementType;
      ins->wSolver->NBCType = NBCType;
      ins->wSolver->BCType = (int*) calloc(NBCType,sizeof(int));
      memcpy(ins->wSolver->BCType,wBCType,NBCType*sizeof(int));
      ins->wSolver->var_coeff = ins->var_coeff;
      ins->wSolver->lambda = ins->ellipticCoeff; 
      ins->wSolver->o_lambda = ins->o_ellipticCoeff; 
      ins->wSolver->loffset = 0;
 
      ellipticSolveSetup(ins->wSolver, kernelInfoV);
    }
  }

  } // flow

  // setup scalar solver
  if(ins->Nscalar) {
   mesh_t *msh;
   (ins->cht) ? msh = ins->meshT : msh = ins->mesh;
   ins->cds = cdsSetup(ins, msh, options, kernelInfoS); 
  }

  if (ins->flow) {

  if (mesh->rank==0) printf("==================PRESSURE SETUP=========================\n");

  ins->presTOL = 1E-4;

  int *pBCType = (int*) calloc(NBCType, sizeof(int));
  for (int bID=1; bID <= nbrBIDs; bID++) {
    pBCType[bID] = bcMap::type(bID, "pressure");
  }

  ins->pOptions = options;
  ins->pOptions.setArgs("KRYLOV SOLVER",        options.getArgs("PRESSURE KRYLOV SOLVER"));
  ins->pOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("PRESSURE SOLVER TOLERANCE"));
  ins->pOptions.setArgs("DISCRETIZATION",       options.getArgs("PRESSURE DISCRETIZATION"));
  ins->pOptions.setArgs("BASIS",                options.getArgs("PRESSURE BASIS"));
  ins->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRESSURE PRECONDITIONER"));
  ins->pOptions.setArgs("MULTIGRID COARSENING", options.getArgs("PRESSURE MULTIGRID COARSENING"));
  ins->pOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
  ins->pOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE",  options.getArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE"));
  ins->pOptions.setArgs("PARALMOND CYCLE",      options.getArgs("PRESSURE PARALMOND CYCLE"));
  ins->pOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("PRESSURE PARALMOND SMOOTHER"));
  ins->pOptions.setArgs("PARALMOND PARTITION",  options.getArgs("PRESSURE PARALMOND PARTITION"));
  ins->pOptions.setArgs("PARALMOND CHEBYSHEV DEGREE",  options.getArgs("PRESSURE PARALMOND CHEBYSHEV DEGREE"));
  ins->pOptions.setArgs("PARALMOND AGGREGATION STRATEGY", options.getArgs("PRESSURE PARALMOND AGGREGATION STRATEGY"));
  ins->pOptions.setArgs("DEBUG ENABLE OGS", "1");
  ins->pOptions.setArgs("DEBUG ENABLE REDUCTIONS", "1");
  ins->pOptions.setArgs("MULTIGRID VARIABLE COEFFICIENT", "FALSE");

  ins->pSolver = new elliptic_t();
  ins->pSolver->blockSolver = 0;
  ins->pSolver->Nfields = 1; 
  ins->pSolver->Ntotal = ins->fieldOffset;
  ins->pSolver->wrk = scratch + ins->ellipticWrkOffset; 
  ins->pSolver->o_wrk = o_scratch.slice(ins->ellipticWrkOffset*sizeof(dfloat)); 
  ins->pSolver->mesh = mesh;
  ins->pSolver->options = ins->pOptions;
  ins->pSolver->dim = ins->dim;
  ins->pSolver->elementType = ins->elementType;
  ins->pSolver->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
  memcpy(ins->pSolver->BCType,pBCType,(nbrBIDs+1)*sizeof(int));
  ins->pSolver->var_coeff = 1;
  // coeff used by ellipticSetup to detect allNeumann 
  // and coeff[0] to setup MG levels
  for (int i=0;i<2*ins->fieldOffset;i++) ins->ellipticCoeff[i] = 0; 
  ins->pSolver->lambda = ins->ellipticCoeff;
  ins->pSolver->o_lambda = ins->o_ellipticCoeff;
  ins->pSolver->loffset = 0;

  ellipticSolveSetup(ins->pSolver, kernelInfoP);

  // setup boundary mapping
  dfloat largeNumber = 1<<20;
  ins->VmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (int e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) ins->VmapB[n+e*mesh->Np] = largeNumber;
  }

  ins->EToB = (int*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(int));

  int cnt = 0;
  for (int e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = bcMap::id(mesh->EToB[f+e*mesh->Nfaces], "velocity");
      ins->EToB[cnt] = bc;
      if (bc>0) {
	for (int n=0;n<mesh->Nfp;n++) {
	  int fid = mesh->faceNodes[n+f*mesh->Nfp];
          ins->VmapB[fid+e*mesh->Np] = mymin(bc,ins->VmapB[fid+e*mesh->Np]); // Dirichlet wins
	}
      }
      cnt++;
    }
  }

  ogsGatherScatter(ins->VmapB, ogsInt, ogsMin, mesh->ogs); 
  for (int n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (ins->VmapB[n] == largeNumber) ins->VmapB[n] = 0;
  }
 
  ins->o_EToB = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int),ins->EToB);
  ins->o_VmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), ins->VmapB);

  } // flow

  // build inverse mass matrix  
  for(hlong e=0;e<mesh->Nelements;++e)
    for(int n=0;n<mesh->Np;++n)
      lumpedMassMatrix[e*mesh->Np+n] = mesh->vgeo[e*mesh->Np*mesh->Nvgeo+JWID*mesh->Np+n];
  ogsGatherScatter(lumpedMassMatrix, ogsDfloat, ogsAdd, mesh->ogs);
  for(int n=0;n<mesh->Np*mesh->Nelements;++n)
    lumpedMassMatrix[n] = 1./lumpedMassMatrix[n];
  ins->o_InvM.copyFrom(lumpedMassMatrix);
  free(lumpedMassMatrix);

  // build kernels
  string fileName, kernelName ;
  const string suffix = "Hex3D";
  const string oklpath = install_dir + "/okl/core/";

  for (int r=0;r<2;r++){
    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {
      
      fileName = oklpath + "insAdvection" + suffix + ".okl";

      kernelName = "insStrongAdvectionVolume" + suffix;
      ins->advectionStrongVolumeKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
      kernelName = "insStrongAdvectionCubatureVolume" + suffix;
      ins->advectionStrongCubatureVolumeKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
      
      fileName = oklpath + "insPressureAx" + suffix + ".okl";
      kernelName = "insPressureAx" + suffix;
      ins->pressureAxKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  
     
      fileName = oklpath + "insCurl" + suffix + ".okl";
      kernelName = "insCurl" + suffix;
      ins->curlKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);  

      fileName = oklpath + "insMassMatrix" + ".okl";
      kernelName = "insMassMatrix" + suffix;
      ins->massMatrixKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);  

      kernelName = "insInvMassMatrix" + suffix;
      ins->invMassMatrixKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);  

      fileName = oklpath + "insGradient" + suffix + ".okl";
      kernelName = "insGradientVolume" + suffix;
      ins->gradientVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "insSumMakef" + suffix + ".okl";
      kernelName = "insSumMakef" + suffix;
      ins->sumMakefKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insDivergence" + suffix + ".okl";
      kernelName = "insDivergenceVolumeTOMBO" + suffix;
      ins->divergenceVolumeKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      kernelName = "insDivergenceSurfaceTOMBO" + suffix;
      ins->divergenceSurfaceKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      fileName = oklpath + "insPressureRhs" + suffix + ".okl";
      kernelName = "insPressureRhsTOMBO" + suffix;
      ins->pressureRhsKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insPressureStress" + suffix + ".okl";
      kernelName = "insPressureStress" + suffix;
      ins->pressureStressKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insPressureBC" + suffix + ".okl";
      kernelName = "insPressureAddBCTOMBO" + suffix;
      ins->pressureAddBCKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      fileName = oklpath + "insPressureUpdate" + ".okl";
      kernelName = "insPressureUpdateTOMBO";
      ins->pressureUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
  
      fileName = oklpath + "insVelocityRhs" + suffix + ".okl";
      kernelName = "insVelocityRhsTOMBO" + suffix;
      ins->velocityRhsKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
     
      fileName = oklpath + "insVelocityBC" + suffix + ".okl";
      kernelName = "insVelocityBC" + suffix;
      ins->velocityRhsBCKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      kernelName = "insVelocityAddBC" + suffix;
      ins->velocityAddBCKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfoBC);

      fileName = oklpath + "insSubCycle" + suffix + ".okl";
      kernelName = "insSubCycleStrongCubatureVolume" + suffix;
      ins->subCycleStrongCubatureVolumeKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      kernelName = "insSubCycleStrongVolume" + suffix; 
      ins->subCycleStrongVolumeKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insSubCycleRKUpdate" + ".okl";
      kernelName = "insSubCycleLSERKUpdate";
      if(ins->SNrk==4) kernelName = "insSubCycleERKUpdate";
      ins->subCycleRKUpdateKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insVelocityExt" + ".okl";
      kernelName = "insVelocityExt";
      ins->velocityExtKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      // ===========================================================================

      fileName = install_dir + "/libparanumal/okl/scaledAdd.okl";
      kernelName = "scaledAddwOffset";
      ins->scaledAddKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = install_dir + "/libparanumal/okl/addScalar.okl";
      kernelName = "setScalar";
      ins->setScalarKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = install_dir + "/libparanumal/okl/dotMultiply.okl";
      kernelName = "dotMultiply";
      ins->dotMultiplyKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "math" + ".okl";
      kernelName = "max";
      ins->maxKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      kernelName = "scalarScaledAdd";
      ins->scalarScaledAddKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      // ===========================================================================

      fileName = oklpath + "insFilterRT" + suffix + ".okl";
      kernelName = "insFilterRT" + suffix;
      ins->filterRTKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insCfl" + suffix + ".okl";
      kernelName = "insCfl" + suffix;
      ins->cflKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insQtl" + suffix + ".okl";
      kernelName = "insQtl" + suffix;
      ins->qtlKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insPressureAddQtl" + ".okl";
      kernelName = "insPressureAddQtl";
      ins->pressureAddQtlKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "setEllipticCoeff.okl"; 
      kernelName = "setEllipticCoeff";
      ins->setEllipticCoeffKernel =  
        mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "setEllipticCoeffPressure";
      ins->setEllipticCoeffPressureKernel =  
        mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "insPQ.okl"; 
      kernelName = "insPQ";
      ins->pqKernel =  
        mesh->device.buildKernel(fileName, kernelName, kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }

  if(!buildOnly) {
    int err = 0;
    dlong gNelements = mesh->Nelements;
    MPI_Allreduce(MPI_IN_PLACE, &gNelements, 1, MPI_DLONG, MPI_SUM, mesh->comm);
    const dfloat sum2 = (dfloat)gNelements * mesh->Np;
    ins->setScalarKernel(ins->fieldOffset, 1.0, ins->o_wrk0);
    ogsGatherScatter(ins->o_wrk0, ogsDfloat, ogsAdd, mesh->ogs);
    ins->dotMultiplyKernel(
         Nlocal,
         mesh->ogs->o_invDegree,
         ins->o_wrk0,
         ins->o_wrk1);
    dfloat *tmp = (dfloat *) calloc(Nlocal, sizeof(dfloat));
    ins->o_wrk1.copyTo(tmp, Nlocal*sizeof(dfloat));
    dfloat sum1 = 0;
    for(int i=0; i<Nlocal; i++) sum1 += tmp[i];
    MPI_Allreduce(MPI_IN_PLACE, &sum1, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
    sum1 = abs(sum1-sum2)/sum2;
    if(sum1 > 1e-15) {
      if(mesh->rank==0) printf("ogsGatherScatter test err=%g!\n", sum1); 
      fflush(stdout);
      err++;
    }

    mesh->ogs->o_invDegree.copyTo(tmp, Nlocal*sizeof(dfloat));
    double *vmult = (double *) nek_ptr("vmult");
    sum1 = 0;
    for(int i=0; i<Nlocal; i++) sum1 += abs(tmp[i] - vmult[i]);
    MPI_Allreduce(MPI_IN_PLACE, &sum1, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
    if(sum1 > 1e-15) {
      if(mesh->rank==0) printf("multiplicity test err=%g!\n", sum1); 
      fflush(stdout);
      err++;
    }

    if(err) ABORT(1);
    free(tmp);
  }

  return ins;
}

cds_t *cdsSetup(ins_t *ins, mesh_t *mesh, setupAide options, occa::properties &kernelInfoH)
{
  cds_t *cds = new cds_t(); 
  cds->mesh = mesh;
 
  if (mesh->rank==0) 
    cout << "==================SCALAR SETUP===========================\n";
                          
  string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  
  // set mesh, options
  cds->meshV       = ins->mesh; 
  cds->elementType = ins->elementType; 
  cds->dim         = ins->dim; 
  cds->NVfields    = ins->NVfields;
  cds->NSfields    = ins->Nscalar;

  cds->extbdfA = ins->extbdfA;
  cds->extbdfB = ins->extbdfB;
  cds->extbdfC = ins->extbdfC;
  cds->extC    = ins->extC;

  cds->Nstages       = ins->Nstages; 
  cds->temporalOrder = ins->temporalOrder; 
  cds->g0            = ins->g0; 

  cds->o_usrwrk = &(ins->o_usrwrk);

  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  cds->Nlocal  = Nlocal;
  cds->Ntotal  = Ntotal;

  cds->vFieldOffset = ins->fieldOffset;
  cds->fieldOffset  = ins->fieldOffset;
  cds->Nblock       = (Nlocal+blockSize-1)/blockSize;

  cds->o_wrk0 = ins->o_wrk0;
  cds->o_wrk1 = ins->o_wrk1;
  cds->o_wrk2 = ins->o_wrk2;
  cds->o_wrk3 = ins->o_wrk3;
  cds->o_wrk4 = ins->o_wrk4;
  cds->o_wrk5 = ins->o_wrk5;
  cds->o_wrk6 = ins->o_wrk6;

  // Solution storage at interpolation nodes
  cds->U     = ins->U; // Point to INS side Velocity
  cds->S     = (dfloat*) calloc(cds->NSfields*(cds->Nstages+0)*cds->fieldOffset,sizeof(dfloat));
  cds->BF    = (dfloat*) calloc(cds->NSfields*cds->fieldOffset,sizeof(dfloat));
  cds->FS    = (dfloat*) calloc(cds->NSfields*(cds->Nstages+1)*cds->fieldOffset,sizeof(dfloat));

  cds->Nsubsteps = ins->Nsubsteps; 
  if(cds->Nsubsteps){
    cds->SNrk   = ins->SNrk;  
    cds->Srka   = ins->Srka; 
    cds->Srkb   = ins->Srkb; 
    cds->Srkc   = ins->Srkc;
    cds->o_Srka = ins->o_Srka;
    cds->o_Srkb = ins->o_Srkb; 
  }

  cds->startTime =ins->startTime;
  cds->dt  = ins->dt; 
  cds->idt = 1.0/cds->dt;
  cds->sdt = ins->sdt; 
  cds->NtimeSteps = ins->NtimeSteps; 

  cds->prop = (dfloat*) calloc(cds->NSfields*2*cds->fieldOffset,sizeof(dfloat));
  for(int is=0; is<cds->NSfields; is++) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << is;
    string sid = ss.str(); 

    if(options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) continue;

    dfloat diff = 1;
    dfloat rho = 1;
    options.getArgs("SCALAR" + sid + " DIFFUSIVITY", diff);
    options.getArgs("SCALAR" + sid + " DENSITY", rho);

    const dlong off = cds->NSfields*cds->fieldOffset;
    for (int e=0;e<mesh->Nelements;e++) { 
      for (int n=0;n<mesh->Np;n++) { 
        cds->prop[0*off + is*cds->fieldOffset + e*mesh->Np + n] = diff;
        cds->prop[1*off + is*cds->fieldOffset + e*mesh->Np + n] = rho;
      }
    }
  }
  cds->o_prop = mesh->device.malloc(cds->NSfields*2*cds->fieldOffset*sizeof(dfloat), cds->prop);  
  cds->o_diff = cds->o_prop.slice(0*cds->NSfields*cds->fieldOffset*sizeof(dfloat));
  cds->o_rho  = cds->o_prop.slice(1*cds->NSfields*cds->fieldOffset*sizeof(dfloat));

  cds->var_coeff = 1; // use always var coeff elliptic
  cds->ellipticCoeff   = ins->ellipticCoeff;
  cds->o_ellipticCoeff = ins->o_ellipticCoeff;  

  cds->o_U  = ins->o_U;
  cds->o_Ue = ins->o_Ue;
  cds->o_S  = mesh->device.malloc(cds->NSfields*(cds->Nstages+0)*cds->fieldOffset*sizeof(dfloat), cds->S);
  cds->o_BF = mesh->device.malloc(cds->NSfields*cds->fieldOffset*sizeof(dfloat), cds->BF);
  cds->o_FS = mesh->device.malloc(cds->NSfields*(cds->Nstages+1)*cds->fieldOffset*sizeof(dfloat), cds->FS);

  cds->options = options;
  cds->options.setArgs("KRYLOV SOLVER",        options.getArgs("SCALAR SOLVER"));
  cds->options.setArgs("DISCRETIZATION",       options.getArgs("SCALAR DISCRETIZATION"));
  cds->options.setArgs("BASIS",                options.getArgs("SCALAR BASIS"));
  /*
  cds->options.setArgs("MULTIGRID COARSENING", options.getArgs("SCALAR MULTIGRID COARSENING"));
  cds->options.setArgs("MULTIGRID SMOOTHER",   options.getArgs("SCALAR MULTIGRID SMOOTHER"));
  cds->options.setArgs("MULTIGRID CHEBYSHEV DEGREE",  options.getArgs("SCALAR MULTIGRID CHEBYSHEV DEGREE")); 
  cds->options.setArgs("PARALMOND CYCLE",      options.getArgs("SCALAR PARALMOND CYCLE"));
  cds->options.setArgs("PARALMOND SMOOTHER",   options.getArgs("SCALAR PARALMOND SMOOTHER"));
  cds->options.setArgs("PARALMOND PARTITION",  options.getArgs("SCALAR PARALMOND PARTITION"));
  cds->options.setArgs("PARALMOND CHEBYSHEV DEGREE",  options.getArgs("SCALAR PARALMOND CHEBYSHEV DEGREE"));
  cds->options.setArgs("PARALMOND AGGREGATION STRATEGY", options.getArgs("SCALAR PARALMOND AGGREGATION STRATEGY"));
  */
  cds->options.setArgs("DEBUG ENABLE OGS", "1");
  cds->options.setArgs("DEBUG ENABLE REDUCTIONS", "1");

  cds->TOL = 1e-6;

  for (int is=0; is<cds->NSfields; is++) {
    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh; // only first scalar can be a CHT mesh

    int nbrBIDs = bcMap::size(0);
    if(ins->cht && is==0) nbrBIDs = bcMap::size(1);
    int *sBCType = (int*) calloc(nbrBIDs+1, sizeof(int));

    std::stringstream ss;
    ss  << std::setfill('0') << std::setw(2) << is;
    string sid = ss.str(); 

    cds->compute[is] = 1;
    if (options.compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
      cds->compute[is] = 0;
      continue;
    }

    for (int bID=1; bID <= nbrBIDs; bID++) {
      string bcTypeText(bcMap::text(bID, "scalar" + sid));
      if(mesh->rank == 0) printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str()); 
      sBCType[bID] = bcMap::type(bID, "scalar" + sid);
    }

    cds->options.setArgs("PRECONDITIONER", options.getArgs("SCALAR" + sid + " PRECONDITIONER"));
    cds->options.setArgs("SOLVER TOLERANCE", options.getArgs("SCALAR" + sid +  " SOLVER TOLERANCE"));

    cds->solver[is] = new elliptic_t();
    cds->solver[is]->blockSolver = 0;
    cds->solver[is]->Nfields = 1;
    cds->solver[is]->Ntotal = ins->fieldOffset;
    cds->solver[is]->wrk = scratch + ins->ellipticWrkOffset;
    cds->solver[is]->o_wrk = o_scratch.slice(ins->ellipticWrkOffset*sizeof(dfloat));
    cds->solver[is]->mesh = mesh;
    cds->solver[is]->options = cds->options;
    cds->solver[is]->dim = cds->dim;
    cds->solver[is]->elementType = cds->elementType;
    cds->solver[is]->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
    memcpy(cds->solver[is]->BCType,sBCType,(nbrBIDs+1)*sizeof(int));

    cds->solver[is]->var_coeff = cds->var_coeff;
    for (int i=0;i<2*ins->fieldOffset;i++) ins->ellipticCoeff[i] = 1; 
    cds->solver[is]->lambda = cds->ellipticCoeff; 
    cds->solver[is]->o_lambda = cds->o_ellipticCoeff; 
    cds->solver[is]->loffset = 0; 
    ellipticSolveSetup(cds->solver[is], kernelInfoH); 

    // setup boundary mapping
    dfloat largeNumber = 1<<20;
    cds->mapB[is] = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
    int *mapB = cds->mapB[is];
    for (int e=0;e<mesh->Nelements;e++) {
      for (int n=0;n<mesh->Np;n++) mapB[n+e*mesh->Np] = largeNumber;
    }
  
    cds->EToB[is] = (int*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(int));
    int *EToB = cds->EToB[is]; 
 
    int cnt = 0;
    for (int e=0;e<mesh->Nelements;e++) {
      for (int f=0;f<mesh->Nfaces;f++) {
        int bc = bcMap::id(mesh->EToB[f+e*mesh->Nfaces], "scalar" + sid);
        EToB[cnt] = bc;
        if (bc>0) {
          for (int n=0;n<mesh->Nfp;n++) {
            int fid = mesh->faceNodes[n+f*mesh->Nfp];
            mapB[fid+e*mesh->Np] = mymin(bc,mapB[fid+e*mesh->Np]); // Dirichlet wins
	  }
        }
        cnt++;
      }
    }
  
    ogsGatherScatter(mapB, ogsInt, ogsMin, mesh->ogs);
    
    for (int n=0;n<mesh->Nelements*mesh->Np;n++) {
      if (mapB[n] == largeNumber) mapB[n] = 0;
    }
  
    cds->o_EToB[is] = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int), EToB);
    cds->o_mapB[is] = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), mapB);

    free(sBCType);
  }

  // build inverse mass matrix
  dfloat *lumpedMassMatrix = (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  for(hlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      lumpedMassMatrix[e*mesh->Np+n] = mesh->vgeo[e*mesh->Np*mesh->Nvgeo+JWID*mesh->Np+n];
    }
  }
  ogsGatherScatter(lumpedMassMatrix, ogsDfloat, ogsAdd, mesh->ogs);
  for(int n=0;n<mesh->Np*mesh->Nelements;++n)
    lumpedMassMatrix[n] = 1./lumpedMassMatrix[n];
  cds->o_InvM = 
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), lumpedMassMatrix);
  cds->o_InvMV = ins->o_InvM;
  free(lumpedMassMatrix);

  // time stepper
  dfloat rkC[4]  = {1.0, 0.0, -1.0, -2.0};
  cds->o_rkC     = ins->o_rkC;
  cds->o_extbdfA = ins->o_extbdfA;
  cds->o_extbdfB = ins->o_extbdfB;
  cds->o_extbdfC = ins->o_extbdfC; 
  cds->o_extC    = ins->o_extC;
  cds->o_prkA    = ins->o_extbdfC;
  cds->o_prkB    = ins->o_extbdfC;

  // build kernels
  occa::properties kernelInfo = *ins->kernelInfo;
  occa::properties kernelInfoBC = kernelInfo;
  //kernelInfo["defines/" "p_NSfields"]  = cds->NSfields;

  const string bcDataFile = install_dir + "/include/insBcData.h";
  kernelInfoBC["includes"] += bcDataFile.c_str();
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfoBC["includes"] += realpath(boundaryHeaderFileName.c_str(), NULL);

  string fileName, kernelName;
  const string suffix = "Hex3D";
  const string oklpath = install_dir + "/okl/core/";

  for (int r=0;r<2;r++){
    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {

      fileName = oklpath + "cdsAdvection" + suffix + ".okl";

      kernelName = "cdsStrongAdvectionVolume" + suffix;
      cds->advectionStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "cdsStrongAdvectionCubatureVolume" + suffix;
      cds->advectionStrongCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================

      fileName = install_dir + "/libparanumal/okl/addScalar.okl";
      kernelName = "setScalar";
      cds->setScalarKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName   = oklpath + "cdsSumMakef" + suffix + ".okl"; 
      kernelName = "cdsSumMakef" + suffix;
      cds->sumMakefKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      fileName = oklpath + "cdsHelmholtzBC" + suffix + ".okl"; 
      kernelName = "cdsHelmholtzBC" + suffix; 
      cds->helmholtzRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfoBC);

      kernelName = "cdsHelmholtzAddBC" + suffix;
      cds->helmholtzAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfoBC);

      fileName = oklpath + "setEllipticCoeff.okl"; 
      kernelName = "setEllipticCoeff";
      cds->setEllipticCoeffKernel =  
        mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = oklpath + "cdsMassMatrix.okl"; 
      kernelName = "cdsMassMatrix" + suffix;
      cds->massMatrixKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  

      kernelName = "cdsInvMassMatrix" + suffix;
      cds->invMassMatrixKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  

      fileName = oklpath + "cdsFilterRT" + suffix + ".okl"; 
      kernelName = "cdsFilterRT" + suffix;
      cds->filterRTKernel =
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = install_dir + "/libparanumal/okl/scaledAdd.okl";
      kernelName = "scaledAddwOffset";
      cds->scaledAddKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      if(cds->Nsubsteps){
        fileName = oklpath + "cdsSubCycle" + suffix + ".okl"; 
        kernelName = "cdsSubCycleStrongCubatureVolume" + suffix;
        cds->subCycleStrongCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        kernelName = "cdsSubCycleStrongVolume" + suffix;
        cds->subCycleStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        fileName = oklpath + "cdsSubCycleRKUpdate.okl";
        kernelName = "cdsSubCycleLSERKUpdate";
        if(cds->SNrk==4) kernelName = "cdsSubCycleERKUpdate";
        cds->subCycleRKUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }
      fileName = oklpath + "insVelocityExt" + ".okl";
      kernelName = "insVelocityExt";
      cds->velocityExtKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

    }
    MPI_Barrier(mesh->comm);
  }

  return cds;
}
