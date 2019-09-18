#include "nekrs.hpp"
#include "nekInterfaceAdapter.hpp"

#define DIRICHLET 1
#define NEUMANN 2

ins_t *setup(mesh_t *mesh, setupAide &options){
  ins_t *ins = new ins_t();
  
  ins->mesh = mesh;
  ins->options = options;
  ins->kernelInfo = new occa::properties();

  string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));

  options.getArgs("MESH DIMENSION", ins->dim);
  options.getArgs("ELEMENT TYPE", ins->elementType);
  
  ins->NVfields = (ins->dim==3) ? 3:2; //  Total Number of Velocity Fields
  ins->NTfields = (ins->dim==3) ? 4:3; // Total Velocity + Pressure

  ins->SNrk = 0;
  
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

  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  
  ins->Ntotal = Ntotal;
  ins->fieldOffset = Ntotal;
  ins->Nblock = (Nlocal+blockSize-1)/blockSize;

  // compute samples of q at interpolation nodes
  ins->U     = (dfloat*) calloc(ins->NVfields*ins->Nstages*Ntotal,sizeof(dfloat));
  ins->P     = (dfloat*) calloc(              ins->Nstages*Ntotal,sizeof(dfloat));

  //rhs storage
  ins->rhsU  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsV  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsW  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsP  = (dfloat*) calloc(Ntotal,sizeof(dfloat));

  //additional field storage
  ins->NU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  ins->LU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  ins->GP   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));

  ins->GU   = (dfloat*) calloc(ins->NVfields*Ntotal*4,sizeof(dfloat));
  
  ins->rkU  = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkP  = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  ins->PI   = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  
  ins->rkNU = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkLU = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkGP = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));

  ins->FU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));

  //extra storage for interpolated fields
  ins->cU = (dfloat *) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));

  options.getArgs("SUBCYCLING STEPS",ins->Nsubsteps);

  if(ins->Nsubsteps){
    ins->Ud    = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->Ue    = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->resU  = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->rhsUd = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->cUd   = (dfloat*) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));

    // Prepare RK stages for Subcycling Part
   
    int Sorder; 
    options.getArgs("SUBCYCLING TIME ORDER", Sorder);
    if(Sorder==2){
      ins->SNrk     = 2; 
      dfloat rka[2] = {0.0,     1.0 };
      dfloat rkb[2] = {0.5,     0.5 };
      dfloat rkc[2] = {0.0,     1.0 };
      //
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      //
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, ins->SNrk*sizeof(dfloat));
    }else if(Sorder ==3){
      // Using Williamson 3rd order scheme converted to low storage since the better truncation 
      ins->SNrk     = 3; 
      dfloat rka[3] = {0.0,     -5.0/9.0,  -153.0/128.0};
      dfloat rkb[3] = {1.0/3.0, 15.0/16.0,    8.0/15.0 };
      dfloat rkc[3] = {0.0,      1.0/3.0,     3.0/4.0  };
      //
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      //
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, ins->SNrk*sizeof(dfloat));
    }else{
      ins->SNrk     = 5; 
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      // Asumes initialized in mesh, can be moved here
      for(int rk=0; rk<ins->SNrk; rk++){
        ins->Srka[rk] = mesh->rka[rk]; 
        ins->Srkb[rk] = mesh->rkb[rk]; 
        ins->Srkc[rk] = mesh->rkc[rk]; 
      }
    }
  }

  options.getArgs("VISCOSITY", ins->nu);

  ins->Re = 1/ins->nu;

  occa::properties& kernelInfo = *ins->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();
  kernelInfo["include_paths"].asArray();

  if(ins->dim==3)
      meshOccaSetup3D(mesh, options, kernelInfo);
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  occa::properties kernelInfoV  = kernelInfo;
  occa::properties kernelInfoP  = kernelInfo;

  // ADD-DEFINES
  kernelInfo["defines/" "p_nu"]= ins->nu;

  kernelInfo["defines/" "p_NTfields"]= ins->NTfields;
  kernelInfo["defines/" "p_NVfields"]= ins->NVfields;
  kernelInfo["defines/" "p_NfacesNfp"]=  mesh->Nfaces*mesh->Nfp;
  kernelInfo["defines/" "p_Nstages"]=  ins->Nstages;
  if(ins->Nsubsteps)
    kernelInfo["defines/" "p_SUBCYCLING"]=  1;
  else
    kernelInfo["defines/" "p_SUBCYCLING"]=  0;

  ins->ARKswitch = 0;

  // Struct for BC implementation
  string bcDataFile;
  bcDataFile = install_dir + "/include/insBcData.h";
  kernelInfo["includes"] += bcDataFile.c_str();

  //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += realpath(boundaryHeaderFileName.c_str(), NULL);

  ins->o_U = mesh->device.malloc(ins->NVfields*ins->Nstages*Ntotal*sizeof(dfloat), ins->U);
  ins->o_P = mesh->device.malloc(              ins->Nstages*Ntotal*sizeof(dfloat), ins->P);

  if(ins->Nsubsteps){
    // Note that resU and resV can be replaced with already introduced buffer
    ins->o_Ue    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ue);
    ins->o_Ud    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ud);
    ins->o_resU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->resU);
    ins->o_rhsUd = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rhsUd);
    ins->o_cUd = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cUd);
  }

  dfloat dt; 
  options.getArgs("DT", dt);
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(ins->dti), 1, MPI_DFLOAT, MPI_MIN, mesh->comm);
  ins->dt = ins->dti; 

  options.getArgs("FINAL TIME", ins->finalTime);
  options.getArgs("START TIME", ins->startTime);
 
  ins->NtimeSteps = ceil((ins->finalTime-ins->startTime)/ins->dt);

  if(ins->Nsubsteps){
    ins->dt         = ins->Nsubsteps*ins->dt;
    ins->NtimeSteps = ceil((ins->finalTime-ins->startTime)/ins->dt);
    ins->dt         = (ins->finalTime-ins->startTime)/ins->NtimeSteps;
    ins->sdt        = ins->dt/ins->Nsubsteps;
  } else{
    ins->NtimeSteps = ceil((ins->finalTime-ins->startTime)/ins->dt);
    ins->dt         = (ins->finalTime-ins->startTime)/ins->NtimeSteps;
  }

  // Hold some inverses for kernels
  ins->inu = 1.0/ins->nu; 
  ins->idt = 1.0/ins->dt;
  
  ins->lambda = ins->g0 / (ins->dt * ins->nu);

  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", ins->outputStep);

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


  const int velMap[7]  = {0,DIRICHLET,DIRICHLET,NEUMANN,NEUMANN,NEUMANN,NEUMANN};
  const int prMap[7]   = {0,NEUMANN,NEUMANN,DIRICHLET,NEUMANN,NEUMANN,NEUMANN};
  const int scalMap[4] = {0,DIRICHLET,NEUMANN,NEUMANN};
  const int nbrBIDs = nekData.NboundaryID;

  int *uBCType = (int*) calloc(nbrBIDs+1, sizeof(int));
  int *vBCType = (int*) calloc(nbrBIDs+1, sizeof(int));
  int *wBCType = (int*) calloc(nbrBIDs+1, sizeof(int));
  int *pBCType = (int*) calloc(nbrBIDs+1, sizeof(int));
  //int *sBCType = (int*) calloc(nbrBIDs+1, sizeof(int));

  // boundary IDs are contiguous and start from 1 
  for (int bID=1; bID <= nbrBIDs; bID++) {
    int bcID = nek_bcmap(bID, 1);
    uBCType[bID] = vBCType[bID] = wBCType[bID] = velMap[bcID];
    if (bcID == 4) uBCType[bID] = DIRICHLET; 
    if (bcID == 5) vBCType[bID] = DIRICHLET; 
    if (bcID == 6) wBCType[bID] = DIRICHLET;
    pBCType[bID] = prMap[bcID];;

    //bcID = nek_bcmap(bID, 2);
    //sBCType[bID] = scalMap[bcID];
  }
 
  //default solver tolerances 
  ins->presTOL = 1E-4;
  ins->velTOL  = 1E-6;

  if (mesh->rank==0) printf("==================VELOCITY SOLVE SETUP=========================\n");

  ins->uSolver = new elliptic_t();
  ins->uSolver->mesh = mesh;
  ins->uSolver->options = ins->vOptions;
  ins->uSolver->dim = ins->dim;
  ins->uSolver->elementType = ins->elementType;
  ins->uSolver->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
  memcpy(ins->uSolver->BCType,uBCType,(nbrBIDs+1)*sizeof(int));

  ellipticSolveSetup(ins->uSolver, ins->lambda, kernelInfoV); 

  ins->vSolver = new elliptic_t();
  ins->vSolver->mesh = mesh;
  ins->vSolver->options = ins->vOptions;
  ins->vSolver->dim = ins->dim;
  ins->vSolver->elementType = ins->elementType;
  ins->vSolver->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
  memcpy(ins->vSolver->BCType,vBCType,(nbrBIDs+1)*sizeof(int));
  
  ellipticSolveSetup(ins->vSolver, ins->lambda, kernelInfoV); //!!!!!

  if (ins->dim==3) {
    ins->wSolver = new elliptic_t();
    ins->wSolver->mesh = mesh;
    ins->wSolver->options = ins->vOptions;
    ins->wSolver->dim = ins->dim;
    ins->wSolver->elementType = ins->elementType;
    ins->wSolver->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
    memcpy(ins->wSolver->BCType,wBCType,(nbrBIDs+1)*sizeof(int));
    
    ellipticSolveSetup(ins->wSolver, ins->lambda, kernelInfoV);  //!!!!!
  }
  
  if (mesh->rank==0) printf("==================PRESSURE SOLVE SETUP=========================\n");
  ins->pSolver = new elliptic_t();
  ins->pSolver->mesh = mesh;
  ins->pSolver->options = ins->pOptions;
  ins->pSolver->dim = ins->dim;
  ins->pSolver->elementType = ins->elementType;
  ins->pSolver->BCType = (int*) calloc(nbrBIDs+1,sizeof(int));
  memcpy(ins->pSolver->BCType,pBCType,(nbrBIDs+1)*sizeof(int));
  
  ellipticSolveSetup(ins->pSolver, 0.0, kernelInfoP); //!!!!

  // TW: this code needs to be re-evaluated from here ====>
  //make node-wise boundary flags
  dfloat largeNumber = 1<<20;
  ins->VmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  ins->PmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (int e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) ins->VmapB[n+e*mesh->Np] = largeNumber;
  }

  ins->EToB = (int*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(int));

  int cnt = 0;
  for (int e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      ins->EToB[cnt] = nek_bcmap(mesh->EToB[f+e*mesh->Nfaces], 1);
      int bc = ins->EToB[cnt];
      if (bc>0) {
	for (int n=0;n<mesh->Nfp;n++) {
	  int fid = mesh->faceNodes[n+f*mesh->Nfp];
	  ins->VmapB[fid+e*mesh->Np] = mymin(bc,ins->VmapB[fid+e*mesh->Np]);
	  ins->PmapB[fid+e*mesh->Np] = mymax(bc,ins->PmapB[fid+e*mesh->Np]);
	}
      }
      cnt++;
    }
  }

  // this is potentially not good for Neumann (since it will )
  ogsGatherScatter(ins->VmapB, ogsInt, ogsMin, mesh->ogs); 
  ogsGatherScatter(ins->PmapB, ogsInt, ogsMax, mesh->ogs); 
  for (int n=0;n<mesh->Nelements*mesh->Np;n++) {

    if (ins->VmapB[n] == largeNumber) {
      ins->VmapB[n] = 0.;
    }
  }
 
  ins->o_VmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), ins->VmapB);
  ins->o_PmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), ins->PmapB);

  // TW: this code needs to be re-evaluated to here <======
  // and the kernels that use VmapB, PmapB

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


  if (mesh->rank==0) {
    printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);
  }
  
  dfloat rkC[4] = {1.0, 0.0, -1.0, -2.0};

  ins->o_rkC  = mesh->device.malloc(4*sizeof(dfloat),rkC);
  ins->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
  ins->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
  ins->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 

  ins->o_extC = mesh->device.malloc(3*sizeof(dfloat)); 

  ins->o_prkA = ins->o_extbdfC;
  ins->o_prkB = ins->o_extbdfC;

  // dummy decleration for scratch space 
  ins->Wrk     = (dfloat*) calloc(1, sizeof(dfloat));
  ins->o_Wrk   = mesh->device.malloc(1*sizeof(dfloat), ins->Wrk);
 
  // MEMORY ALLOCATION
  ins->o_rhsU  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsU);
  ins->o_rhsV  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsV);
  ins->o_rhsW  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsW);
  ins->o_rhsP  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsP);
  //storage for helmholtz solves
  ins->o_UH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  ins->o_VH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  ins->o_WH = mesh->device.malloc(Ntotal*sizeof(dfloat));

  ins->o_FU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->FU);
  ins->o_NU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->NU);
  ins->o_GP    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->GP);
  ins->o_rkGP  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkGP);

  ins->o_NC = ins->o_GP; // Use GP storage to store curl(curl(u)) history

  // build lumped mass matrix for NEK
  dfloat *lumpedMassMatrix     = (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  dfloat *copyLumpedMassMatrix = (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));

  for(hlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      lumpedMassMatrix[e*mesh->Np+n]     = mesh->vgeo[e*mesh->Np*mesh->Nvgeo+JWID*mesh->Np+n];
      copyLumpedMassMatrix[e*mesh->Np+n] = mesh->vgeo[e*mesh->Np*mesh->Nvgeo+JWID*mesh->Np+n];
    }
  }

  ogsGatherScatter(lumpedMassMatrix, ogsDfloat, ogsAdd, mesh->ogs);

  for(int n=0;n<mesh->Np*mesh->Nelements;++n)
    lumpedMassMatrix[n] = 1./lumpedMassMatrix[n];

  ins->o_invLumpedMassMatrix = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), lumpedMassMatrix);
  ins->o_InvM = ins->o_invLumpedMassMatrix; 
 
  ins->o_rkU   = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkU);
  ins->o_rkP   = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->rkP);
  ins->o_PI    = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->PI);

  ins->o_cU = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cU);

  // halo setup
  if(mesh->totalHaloPairs){
    dlong vHaloBytes = mesh->totalHaloPairs*mesh->Np*ins->NVfields*sizeof(dfloat);
    dlong pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    dlong vGatherBytes = ins->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat);

    ins->o_vHaloBuffer = mesh->device.malloc(vHaloBytes);
    ins->o_pHaloBuffer = mesh->device.malloc(pHaloBytes);

    ins->vSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, NULL, ins->o_vSendBuffer, ins->h_vSendBuffer);
    ins->vRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, NULL, ins->o_vRecvBuffer, ins->h_vRecvBuffer);

    ins->pSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, NULL, ins->o_pSendBuffer, ins->h_pSendBuffer);
    ins->pRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, NULL, ins->o_pRecvBuffer, ins->h_pRecvBuffer);

    ins->velocityHaloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, vGatherBytes, NULL, 
                                           ins->o_gatherTmpPinned, ins->h_gatherTmpPinned);
    ins->o_velocityHaloGatherTmp = mesh->device.malloc(vGatherBytes,  ins->velocityHaloGatherTmp);
  }

  // load kernels
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
      
      // Create lap(U) kernel for Ins for (lap(u) - lambda*u) operation
      fileName = oklpath + "insAx" + suffix + ".okl";
      kernelName = "insPressureAx" + suffix;
      ins->pressureAxKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  
     
      fileName = oklpath + "insCurl" + suffix + ".okl";
      kernelName = "insCurl" + suffix;
      ins->curlKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);  
      if(ins->dim==2){
        kernelName = "insCurlB" + suffix;
        ins->curlBKernel = 
          mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo); 
      }

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

      kernelName = "insGradientSurface" + suffix;
      ins->gradientSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

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
      kernelName = "insPressureBC" + suffix;
      ins->pressureRhsBCKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

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

      fileName = oklpath + "insVelocityUpdate" + ".okl";
      kernelName = "insVelocityUpdate";
      ins->velocityUpdateKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);   


      fileName = oklpath + "insSubCycle" + suffix + ".okl";
      kernelName = "insSubCycleStrongCubatureVolume" + suffix;
      ins->subCycleStrongCubatureVolumeKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      kernelName = "insSubCycleStrongVolume" + suffix; 
      ins->subCycleStrongVolumeKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);


      fileName = oklpath + "insSubCycle" + ".okl";
      kernelName = "insSubCycleRKUpdate";
      ins->subCycleRKUpdateKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      kernelName = "insSubCycleExt";
      ins->subCycleExtKernel = 
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

      fileName = oklpath + "insFilter" + ".okl";
      kernelName = "insFilterRT" + suffix;
      ins->filterKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

      fileName = oklpath + "insCfl" + suffix + ".okl";
      kernelName = "insCfl" + suffix;
      ins->cflKernel = 
        mesh->device.buildKernel(fileName.c_str(), kernelName.c_str(), kernelInfo);

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

  ins->o_EToB = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int),
                                    ins->EToB);

  // Initialize Cfl Stuff
  ins->cflComputed = 0; 
   
  if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
    insFilterSetup(ins); 

  return ins;
}
