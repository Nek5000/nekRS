#include "nrs.hpp"
#include "meshSetup.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "filter.hpp"
#include "bcMap.hpp"

cds_t *cdsSetup(ins_t *ins, mesh_t *mesh, setupAide &options, occa::properties &kernelInfoH);
              
 
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

  int cht = 0;
  if (nekData.nelv != nekData.nelt) cht = 1;

  if (buildOnly) {
    ins->meshT = createMeshDummy(comm, N, options, kernelInfo);
    ins->mesh = ins->meshT;
  } else {
    ins->meshT = createMeshT(comm, N, cht, options, kernelInfo);
    ins->mesh = ins->meshT;
    if (cht) ins->mesh = createMeshV(comm, N, ins->meshT, options, kernelInfo);
  }
  mesh_t *mesh = ins->mesh;

  ins->NVfields = (ins->dim==3) ? 3:2; //  Total Number of Velocity Fields
  ins->NTfields = (ins->dim==3) ? 4:3; // Total Velocity + Pressure

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

  { // ensure that offset is large enough for v and t mesh and is properly aligned 
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

    if(Sorder==2 && ins->SNrk==2){
      dfloat rka[2] = {0.0,     1.0 };
      dfloat rkb[2] = {0.5,     0.5 };
      dfloat rkc[2] = {0.0,     1.0 };
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, ins->SNrk*sizeof(dfloat));
    }else if(Sorder ==3 && ins->SNrk==3){
      // Using Williamson 3rd order scheme converted to low storage since the better truncation 
      dfloat rka[3] = {0.0,     -5.0/9.0,  -153.0/128.0};
      dfloat rkb[3] = {1.0/3.0, 15.0/16.0,    8.0/15.0 };
      dfloat rkc[3] = {0.0,      1.0/3.0,     3.0/4.0  };
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, ins->SNrk*sizeof(dfloat));
    }else if(Sorder==4 && ins->SNrk==4){ // ERK(4,4)
      dfloat rka[4] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
      dfloat rkb[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
      dfloat rkc[4] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, ins->SNrk*sizeof(dfloat));
    }else if(Sorder==4 && ins->SNrk==5){ // LSERK(4,5)
      dfloat rka[5] = {0.0,
                      -567301805773.0/1357537059087.0,
                      -2404267990393.0/2016746695238.0,
                      -3550918686646.0/2091501179385.0,
                      -1275806237668.0/842570457699.0};
      dfloat rkb[5] = {1432997174477.0/9575080441755.0,
                      5161836677717.0/13612068292357.0,
                      1720146321549.0/2090206949498.0,
                      3134564353537.0/4481467310338.0,
                      2277821191437.0/14882151754819.0};
      dfloat rkc[6] = {0.0,
                      1432997174477.0/9575080441755.0,
                      2526269341429.0/6820363962896.0,
                      2006345519317.0/3224310063776.0,
                      2802321613138.0/2924317926251.0,
                      1.};
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk+1, sizeof(dfloat));
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, (ins->SNrk+1)*sizeof(dfloat));
    }else{
      if(mesh->rank==0) cout << "Unsupported subcycling scheme!\n"; 
      MPI_Finalize(); 
      exit(1);
    }
    ins->o_Srka = mesh->device.malloc(ins->SNrk*sizeof(dfloat), ins->Srka);
    ins->o_Srkb = mesh->device.malloc(ins->SNrk*sizeof(dfloat), ins->Srkb);
  }

  occa::properties kernelInfoV  = kernelInfo;
  occa::properties kernelInfoP  = kernelInfo;
  occa::properties kernelInfoS  = kernelInfo;

  // ADD-DEFINES
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

  int NblockV = mymax(1,256/mesh->Np); // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1,256/maxNodes); // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo["defines/" "p_maxNodesVolumeCub"]= maxNodesVolumeCub;
  int cubNblockV = mymax(1,256/maxNodesVolumeCub);
  //
  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo["defines/" "p_maxNodesSurfaceCub"]=maxNodesSurfaceCub;
  int cubNblockS = mymax(256/maxNodesSurfaceCub,1); // works for CUDA
  //
  kernelInfo["defines/" "p_cubNblockV"]=cubNblockV;
  kernelInfo["defines/" "p_cubNblockS"]=cubNblockS;

  // Struct for BC implementation
  string bcDataFile;
  bcDataFile = install_dir + "/include/insBcData.h";
  kernelInfo["includes"] += bcDataFile.c_str();

  //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += realpath(boundaryHeaderFileName.c_str(), NULL);

  // jit compile udf kernels
  if (udf.loadKernels) {
    if (mesh->rank == 0) cout << "building udf kernels ...";
    udf.loadKernels(ins);
    if (mesh->rank == 0) cout << " done" << endl;
  }  

  // setup scratch space
  const int ellipticWrkNflds = 9; 
  ins->ellipticWrkSize = ellipticWrkNflds*ins->fieldOffset;
  const int scratchNflds = 9+ellipticWrkNflds;
  ins->scratch   = (dfloat*) calloc(scratchNflds*ins->fieldOffset,sizeof(dfloat));
  ins->o_scratch = mesh->device.malloc(scratchNflds*ins->fieldOffset*sizeof(dfloat), ins->scratch);

  // dummy decleration for user work space 
  ins->usrwrk   = (dfloat*) calloc(1, sizeof(dfloat));
  ins->o_usrwrk = mesh->device.malloc(1*sizeof(dfloat), ins->usrwrk);

  ins->o_U  = mesh->device.malloc(ins->NVfields*ins->Nstages*ins->fieldOffset*sizeof(dfloat), ins->U);
  ins->o_Ue = mesh->device.malloc(ins->NVfields*ins->fieldOffset*sizeof(dfloat), ins->Ue);
  ins->o_P  = mesh->device.malloc(ins->fieldOffset*sizeof(dfloat), ins->P);
  ins->o_PI = mesh->device.malloc(ins->fieldOffset*sizeof(dfloat), ins->PI);

  ins->o_FU = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*ins->fieldOffset*sizeof(dfloat), ins->FU);
  ins->o_BF = mesh->device.malloc(ins->NVfields*ins->fieldOffset*sizeof(dfloat), ins->BF);

  ins->var_coeff = 1; // use always var coeff elliptic
  ins->ellipticCoeff = (dfloat*) calloc(2*ins->fieldOffset,sizeof(dfloat));
  for (int i=0;i<2*ins->fieldOffset;i++) // just to avoid devision by 0 in Jacobi setup 
      ins->ellipticCoeff[i] = 1;
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
  if(ins->options.compareArgs("LOWMACH", "TRUE")) ins->lowMach = 1;
  ins->qtl   = (dfloat*) calloc(ins->fieldOffset,sizeof(dfloat));
  ins->o_qtl = mesh->device.malloc(ins->fieldOffset*sizeof(dfloat), ins->qtl);  

  dfloat rkC[4]  = {1.0, 0.0, -1.0, -2.0};
  ins->o_rkC     = mesh->device.malloc(4*sizeof(dfloat),rkC);
  ins->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
  ins->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
  ins->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 
  ins->o_extC    = mesh->device.malloc(3*sizeof(dfloat)); 
  ins->o_prkA    = ins->o_extbdfC;
  ins->o_prkB    = ins->o_extbdfC;
 
  if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION")) filterSetup(ins); 
    
  if (mesh->rank==0) printf("==================VELOCITY SETUP=========================\n");

  //make option objects for elliptc solvers
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

  const int nbrBIDs = bcMap::size();
  int *uBCType = (int*) calloc(nbrBIDs+1, sizeof(int));
  int *vBCType = (int*) calloc(nbrBIDs+1, sizeof(int));
  int *wBCType = (int*) calloc(nbrBIDs+1, sizeof(int));
  int *pBCType = (int*) calloc(nbrBIDs+1, sizeof(int));

  for (int bID=1; bID <= nbrBIDs; bID++) {
    string bcTypeText(bcMap::text(bID, "velocity"));
    if(mesh->rank == 0) printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str()); 

    uBCType[bID] = bcMap::type(bID, "x-velocity");
    vBCType[bID] = bcMap::type(bID, "y-velocity");
    wBCType[bID] = bcMap::type(bID, "z-velocity");
    pBCType[bID] = bcMap::type(bID, "pressure");
  }

  //default solver tolerances
  ins->presTOL = 1E-4;
  ins->velTOL  = 1E-6;
 
  ins->uSolver = new elliptic_t();
  ins->uSolver->wrkOffset = ins->fieldOffset;
  ins->uSolver->wrk = ins->scratch; 
  ins->uSolver->o_wrk = ins->o_scratch; 
  ins->uSolver->mesh = mesh;
  ins->uSolver->options = ins->vOptions;
  ins->uSolver->dim = ins->dim;
  ins->uSolver->elementType = ins->elementType;
  ins->uSolver->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
  memcpy(ins->uSolver->BCType,uBCType,(nbrBIDs+1)*sizeof(int));

  ins->uSolver->var_coeff = ins->var_coeff;
  ins->uSolver->coeff = ins->ellipticCoeff; 
  ins->uSolver->o_coeff = ins->o_ellipticCoeff; 
  const dfloat lambda = 1; // not used if var_coeff

  ellipticSolveSetup(ins->uSolver, lambda, kernelInfoV); 

  ins->vSolver = new elliptic_t();
  ins->vSolver->wrkOffset = ins->fieldOffset;
  ins->vSolver->wrk = ins->scratch;
  ins->vSolver->o_wrk = ins->o_scratch;
  ins->vSolver->mesh = mesh;
  ins->vSolver->options = ins->vOptions;
  ins->vSolver->dim = ins->dim;
  ins->vSolver->elementType = ins->elementType;
  ins->vSolver->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
  memcpy(ins->vSolver->BCType,vBCType,(nbrBIDs+1)*sizeof(int));

  ins->vSolver->var_coeff = ins->var_coeff;
  ins->vSolver->coeff = ins->ellipticCoeff; 
  ins->vSolver->o_coeff = ins->o_ellipticCoeff; 
  
  ellipticSolveSetup(ins->vSolver, lambda, kernelInfoV); //!!!!!

  if (ins->dim==3) {
    ins->wSolver = new elliptic_t();
    ins->wSolver->wrkOffset = ins->fieldOffset;
    ins->wSolver->wrk = ins->scratch;
    ins->wSolver->o_wrk = ins->o_scratch;
    ins->wSolver->mesh = mesh;
    ins->wSolver->options = ins->vOptions;
    ins->wSolver->dim = ins->dim;
    ins->wSolver->elementType = ins->elementType;
    ins->wSolver->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
    memcpy(ins->wSolver->BCType,wBCType,(nbrBIDs+1)*sizeof(int));

    ins->wSolver->var_coeff = ins->var_coeff;
    ins->wSolver->coeff = ins->ellipticCoeff; 
    ins->wSolver->o_coeff = ins->o_ellipticCoeff; 
    
    ellipticSolveSetup(ins->wSolver, lambda, kernelInfoV);  //!!!!!
  }
  
  if (mesh->rank==0) printf("==================PRESSURE SETUP=========================\n");
  ins->pSolver = new elliptic_t();
  ins->pSolver->wrkOffset = ins->fieldOffset;
  ins->pSolver->wrk = ins->scratch;
  ins->pSolver->o_wrk = ins->o_scratch;
  ins->pSolver->mesh = mesh;

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

  ins->pSolver->options = ins->pOptions;
  ins->pSolver->dim = ins->dim;
  ins->pSolver->elementType = ins->elementType;
  ins->pSolver->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
  memcpy(ins->pSolver->BCType,pBCType,(nbrBIDs+1)*sizeof(int));

  ins->pSolver->var_coeff = 1;
  ins->pSolver->coeff = ins->ellipticCoeff; 
  ins->pSolver->o_coeff = ins->o_ellipticCoeff;
 
  ellipticSolveSetup(ins->pSolver, 0.0, kernelInfoP); //!!!!

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

  // build inverse mass matrix  
  dfloat *lumpedMassMatrix  = (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  for(hlong e=0;e<mesh->Nelements;++e)
    for(int n=0;n<mesh->Np;++n)
      lumpedMassMatrix[e*mesh->Np+n] = mesh->vgeo[e*mesh->Np*mesh->Nvgeo+JWID*mesh->Np+n];
  ogsGatherScatter(lumpedMassMatrix, ogsDfloat, ogsAdd, mesh->ogs);
  for(int n=0;n<mesh->Np*mesh->Nelements;++n)
    lumpedMassMatrix[n] = 1./lumpedMassMatrix[n];
  ins->o_InvM = 
    mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), lumpedMassMatrix);
 
  // halo setup
  if(mesh->totalHaloPairs){
    dlong vHaloBytes = mesh->totalHaloPairs*mesh->Np*ins->NVfields*sizeof(dfloat);
    dlong pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    dlong vGatherBytes = ins->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat);

    ins->o_vHaloBuffer = mesh->device.malloc(vHaloBytes);
    ins->o_pHaloBuffer = mesh->device.malloc(pHaloBytes);

    ins->vSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, 
                       NULL, ins->o_vSendBuffer, ins->h_vSendBuffer);
    ins->vRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, 
                       NULL, ins->o_vRecvBuffer, ins->h_vRecvBuffer);

    ins->pSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, 
                       NULL, ins->o_pSendBuffer, ins->h_pSendBuffer);
    ins->pRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, 
                       NULL, ins->o_pRecvBuffer, ins->h_pRecvBuffer);

    ins->velocityHaloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, vGatherBytes, NULL, 
                                           ins->o_gatherTmpPinned, ins->h_gatherTmpPinned);
    ins->o_velocityHaloGatherTmp = mesh->device.malloc(vGatherBytes,  ins->velocityHaloGatherTmp);
  }

  // build kernels
  string fileName, kernelName ;
  string suffix;

  if(ins->dim==3)
    suffix = "Hex3D";
  else
    suffix = "Quad2D";

  string oklpath;
  oklpath += install_dir + "/okl/";

  for (int r=0;r<2;r++){
    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {
      
      fileName = oklpath + "insHaloExchange.okl";
 
      kernelName = "insVelocityHaloExtract";
      ins->velocityHaloExtractKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
      kernelName = "insVelocityHaloScatter";
      ins->velocityHaloScatterKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
      kernelName = "insPressureHaloExtract";
      ins->pressureHaloExtractKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
      kernelName = "insPressureHaloScatter";
      ins->pressureHaloScatterKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      // ===========================================================================

      fileName = oklpath + "insAdvection" + suffix + ".okl";

      kernelName = "insStrongAdvectionVolume" + suffix;
      ins->advectionStrongVolumeKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
      kernelName = "insStrongAdvectionCubatureVolume" + suffix;
      ins->advectionStrongCubatureVolumeKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);
      
      fileName = oklpath + "insAx" + suffix + ".okl";
      kernelName = "insAx" + suffix;
      ins->AxKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  
     
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
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      kernelName = "insDivergenceSurfaceTOMBO" + suffix;
      ins->divergenceSurfaceKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insPressureRhs" + suffix + ".okl";
      kernelName = "insPressureRhsTOMBO" + suffix;
      ins->pressureRhsKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insPressureBC" + suffix + ".okl";
      kernelName = "insPressureAddBCTOMBO" + suffix;
      ins->pressureAddBCKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

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
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      kernelName = "insVelocityAddBC" + suffix;
      ins->velocityAddBCKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

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

      fileName = oklpath + "math" + ".okl";
      kernelName = "max";
      ins->maxKernel =  
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

      fileName = oklpath + "insNC.okl"; 
      kernelName = "insNC";
      ins->ncKernel =  
        mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================
      fileName = oklpath + "insHalo.okl";
      kernelName = "insHaloGet";
      ins->haloGetKernel =
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      kernelName = "insHaloPut";
      ins->haloPutKernel =
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

    }
    MPI_Barrier(mesh->comm);
  }

  // setup scalar solve
  if(ins->Nscalar) {
   mesh_t *msh;
   (cht) ? msh = ins->meshT : msh = ins->mesh;
   ins->cds = cdsSetup(ins, msh, options, kernelInfoS); 
  }


  return ins;
}

cds_t *cdsSetup(ins_t *ins, mesh_t *mesh, setupAide &options, occa::properties &kernelInfoH)
{
  cds_t *cds = new cds_t(); 
  cds->mesh = mesh;
 
  if (mesh->rank==0) 
    cout << "==================SCALARS SETUP==========================\n";
                          
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

  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  cds->Nlocal  = Nlocal;
  cds->Ntotal  = Ntotal;

  cds->vFieldOffset = ins->fieldOffset;
  cds->fieldOffset  = ins->fieldOffset;
  cds->Nblock       = (Nlocal+blockSize-1)/blockSize;

  // Solution storage at interpolation nodes
  cds->U     = ins->U; // Point to INS side Velocity
  cds->S     = (dfloat*) calloc(cds->NSfields*(cds->Nstages+0)*cds->fieldOffset,sizeof(dfloat));
  cds->BF    = (dfloat*) calloc(cds->NSfields*cds->fieldOffset,sizeof(dfloat));
  cds->FS    = (dfloat*) calloc(cds->NSfields*(cds->Nstages+1)*cds->fieldOffset,sizeof(dfloat));

  // Use Nsubsteps if INS does to prevent stability issues
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

  occa::properties& kernelInfo  = *ins->kernelInfo; 
  // ADD-DEFINES
  kernelInfo["defines/" "p_NSfields"]  = cds->NSfields;
  kernelInfo["defines/" "p_NTSfields"] = (cds->NVfields+cds->NSfields + 1);
 
  cds->o_U  = ins->o_U;
  cds->o_Ue = ins->o_Ue;
  cds->o_S  = mesh->device.malloc(cds->NSfields*(cds->Nstages+0)*cds->fieldOffset*sizeof(dfloat), cds->S);
  cds->o_BF = mesh->device.malloc(cds->NSfields*cds->fieldOffset*sizeof(dfloat), cds->BF);
  cds->o_FS = mesh->device.malloc(cds->NSfields*(cds->Nstages+1)*cds->fieldOffset*sizeof(dfloat), cds->FS);

  //make option objects for elliptc solvers
  cds->options = options;
  cds->options.setArgs("KRYLOV SOLVER",        options.getArgs("SCALAR KRYLOV SOLVER"));
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

  const int nbrBIDs = bcMap::size();
  int *sBCType = (int*) calloc(nbrBIDs+1, sizeof(int));

  cds->TOL = 1e-6;

  for (int is=0; is<cds->NSfields; is++) {
    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh; // only first scalar can be a CHT mesh

    std::stringstream ss;
    ss  << std::setfill('0') << std::setw(2) << is;
    string sid = ss.str(); 

    for (int bID=1; bID <= nbrBIDs; bID++) {
      string bcTypeText(bcMap::text(bID, "scalar" + sid));
      if(mesh->rank == 0) printf("bID %d -> bcType %s\n", bID, bcTypeText.c_str()); 
      sBCType[bID] = bcMap::type(bID, "scalar" + sid);
    }

    cds->options.setArgs("PRECONDITIONER", options.getArgs("SCALAR" + sid + " PRECONDITIONER"));
    cds->options.setArgs("SOLVER TOLERANCE", options.getArgs("SCALAR" + sid +  " SOLVER TOLERANCE"));

    cds->solver[is] = new elliptic_t();
    cds->solver[is]->wrkOffset = ins->fieldOffset;
    cds->solver[is]->wrk = ins->scratch;
    cds->solver[is]->o_wrk = ins->o_scratch;
    cds->solver[is]->mesh = mesh;
    cds->solver[is]->options = cds->options;
    cds->solver[is]->dim = cds->dim;
    cds->solver[is]->elementType = cds->elementType;
    cds->solver[is]->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
    memcpy(cds->solver[is]->BCType,sBCType,(nbrBIDs+1)*sizeof(int));

    cds->solver[is]->var_coeff = cds->var_coeff;
    cds->solver[is]->coeff = cds->ellipticCoeff; 
    cds->solver[is]->o_coeff = cds->o_ellipticCoeff; 
    const dfloat lambda = 1; // not used if var_coeff
    ellipticSolveSetup(cds->solver[is], lambda, kernelInfoH); 

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
            if(bc != 1 && mapB[fid+e*mesh->Np] != 1)
              mapB[fid+e*mesh->Np] = mapB[fid+e*mesh->Np]; // for Neumann BCs do nothing      
            else
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

  if(mesh->totalHaloPairs){//halo setup
    int npe = mesh->Nfp; 
    dlong haloBytes   = mesh->totalHaloPairs*npe*(cds->NSfields + cds->NVfields)*sizeof(dfloat);
    dlong gatherBytes = (cds->NSfields+cds->NVfields)*mesh->ogs->NhaloGather*sizeof(dfloat);
    cds->o_haloBuffer = mesh->device.malloc(haloBytes);

    cds->sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, haloBytes, NULL, 
                      cds->o_sendBuffer, cds->h_sendBuffer);
    cds->recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, haloBytes, NULL, 
                      cds->o_recvBuffer, cds->h_recvBuffer);
    cds->haloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, gatherBytes, NULL, 
                      cds->o_gatherTmpPinned, cds->h_gatherTmpPinned); 
    cds->o_haloGatherTmp = mesh->device.malloc(gatherBytes,  cds->haloGatherTmp);

    if(cds->Nsubsteps){
      dlong shaloBytes   = mesh->totalHaloPairs*npe*(cds->NSfields)*sizeof(dfloat);
      dlong sgatherBytes = (cds->NSfields)*mesh->ogs->NhaloGather*sizeof(dfloat);
      cds->o_shaloBuffer = mesh->device.malloc(shaloBytes);

      cds->ssendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, shaloBytes, NULL, 
                         cds->o_ssendBuffer, cds->h_ssendBuffer);
      cds->srecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, shaloBytes, NULL,
                         cds->o_srecvBuffer, cds->h_srecvBuffer);
      cds->shaloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, sgatherBytes, NULL, 
                         cds->o_sgatherTmpPinned, cds->h_sgatherTmpPinned);
      cds->o_shaloGatherTmp = mesh->device.malloc(sgatherBytes,  cds->shaloGatherTmp);
    }
  }

  // build kernels
  string suffix, fileName, kernelName;
  if(cds->elementType==QUADRILATERALS)
     suffix = "Quad2D";
  if(cds->elementType==HEXAHEDRA)
     suffix = "Hex3D";

  for (int r=0;r<2;r++){
    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {

      fileName = install_dir + "/okl/cdsHaloExchange.okl";
      
      kernelName = "cdsHaloGet";
      cds->haloGetKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo); 
      
      kernelName = "cdsHaloPut";
      cds->haloPutKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      if(cds->Nsubsteps){
        kernelName =  "cdsScalarHaloGet";
        cds->scalarHaloGetKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo); 
        
        kernelName = "cdsScalarHaloPut";
        cds->scalarHaloPutKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);    
      } 

      fileName = install_dir + "/okl/cdsAdvection" + suffix + ".okl";

      kernelName = "cdsStrongAdvectionVolume" + suffix;
      cds->advectionStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "cdsStrongAdvectionCubatureVolume" + suffix;
      cds->advectionStrongCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================

      fileName = install_dir + "/libparanumal/okl/addScalar.okl";
      kernelName = "setScalar";
      cds->setScalarKernel =  
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName   = install_dir + "/okl/cdsSumMakef" + suffix + ".okl"; 
      kernelName = "cdsSumMakef" + suffix;
      cds->sumMakefKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      fileName = install_dir + "/okl/cdsHelmholtzBC" + suffix + ".okl"; 
      kernelName = "cdsHelmholtzBC" + suffix; 
      cds->helmholtzRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "cdsHelmholtzAddBC" + suffix;
      cds->helmholtzAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = install_dir + "/okl/setEllipticCoeff.okl"; 
      kernelName = "setEllipticCoeff";
      cds->setEllipticCoeffKernel =  
        mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = install_dir + "/okl/cdsMassMatrix.okl"; 
      kernelName = "cdsMassMatrix" + suffix;
      cds->massMatrixKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  

      kernelName = "cdsInvMassMatrix" + suffix;
      cds->invMassMatrixKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  

      fileName = install_dir + "/okl/cdsFilterRT" + suffix + ".okl"; 
      kernelName = "cdsFilterRT" + suffix;
      cds->filterRTKernel =
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      if(cds->Nsubsteps){
        fileName = install_dir + "/libparanumal/okl/scaledAdd.okl";
        kernelName = "scaledAddwOffset";
        cds->scaledAddKernel = 
          mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

        fileName = install_dir + "/okl/cdsSubCycle" + suffix + ".okl"; 
        kernelName = "cdsSubCycleStrongCubatureVolume" + suffix;
        cds->subCycleStrongCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        kernelName = "cdsSubCycleStrongVolume" + suffix;
        cds->subCycleStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        fileName = install_dir + "/okl/cdsSubCycleRKUpdate.okl";
        kernelName = "cdsSubCycleLSERKUpdate";
        if(cds->SNrk==4) kernelName = "cdsSubCycleERKUpdate";
        cds->subCycleRKUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }
      fileName = install_dir + "/okl/insVelocityExt" + ".okl";
      kernelName = "insVelocityExt";
      cds->velocityExtKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

    }
    MPI_Barrier(mesh->comm);
  }

  return cds;
}
