/*
The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nekrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"

void extbdfCoefficents(ins_t *ins, int order) {

  if(order==1) {
     //advection, first order in time, increment
    ins->g0 =  1.0f;
    dfloat extbdfB[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfA[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};

    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder = 1;

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0;
  } else if(order==2) {
    //advection, second order in time, increment
    ins->g0 =  1.5f;
    dfloat extbdfB[3] = {2.0f,-0.5f, 0.0f};
    dfloat extbdfA[3] = {2.0f,-1.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};

    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder=2;

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0;
  } else if(order==3) {
    //advection, third order in time, increment
    ins->g0 =  11.f/6.f;
    dfloat extbdfB[3] = {3.0f,-1.5f, 1.0f/3.0f};
    dfloat extbdfA[3] = {3.0f,-3.0f, 1.0f};
    dfloat extbdfC[3] = {2.0f,-1.0f, 0.0f};

    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder=3;

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0;
  }
}

void runPlan4(ins_t *ins){

  mesh_t *mesh = ins->mesh;

  double time0 = MPI_Wtime();
  for(int tstep=0;tstep<ins->NtimeSteps;++tstep){

    if(tstep<1) 
      extbdfCoefficents(ins,tstep+1);
    else if(tstep<2 && ins->temporalOrder>=2) 
      extbdfCoefficents(ins,tstep+1);
    else if(tstep<3 && ins->temporalOrder>=3) 
      extbdfCoefficents(ins,tstep+1);

    dfloat time = ins->startTime + tstep*ins->dt;

    dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

    if(ins->Nsubsteps) {
      libParanumal::insStrongSubCycle(ins, time, ins->Nstages, ins->o_U, ins->o_NU);
    } else {
      libParanumal::insAdvection(ins, time, ins->o_U, ins->o_NU);
    }

    insCurlCurl(ins, time, ins->o_U, ins->o_NC);

    ins->setScalarKernel(ins->Ntotal*ins->NVfields, 0.0, ins->o_FU);
    if(udf.velocityForce) udf.velocityForce(ins, time+ins->dt, ins->o_U, ins->o_FU);

    if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
      ins->filterKernel(mesh->Nelements,
                        ins->o_filterMT,
                        ins->filterS,
                        ins->fieldOffset,
                        ins->o_U,
                        ins->o_FU);

    insPressureRhs  (ins, time+ins->dt, ins->Nstages);
    insPressureSolve(ins, time+ins->dt, ins->Nstages); 
    insPressureUpdate(ins, time+ins->dt, ins->Nstages, ins->o_rkP);
   
    // copy updated pressure
    ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat)); 

    insVelocityRhs  (ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
    insVelocitySolve(ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);

    //cycle velocity history after velocityUpdate
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_U.copyFrom(ins->o_U, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    }

    //copy updated pressure
    ins->o_U.copyFrom(ins->o_rkU, ins->NVfields*ins->Ntotal*sizeof(dfloat)); 

    //cycle rhs history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_NU.copyFrom(ins->o_NU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
       (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
       (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      
      ins->o_NC.copyFrom(ins->o_NC, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
       (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
       (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      
      ins->o_FU.copyFrom(ins->o_FU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
       (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
       (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    }

    dfloat cfl = insComputeCfl(ins, time+ins->dt, tstep+1);
    if (mesh->rank==0)
      printf("tstep = %d, time = %.5e, cfl = %.2f, iter: U - %3d, V - %3d, W - %3d, P - %3d, elapsed = %.5e s\n",
        tstep+1, time+ins->dt, cfl, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP, MPI_Wtime()-time0);

    if (ins->outputStep > 0)
      if (((tstep+1)%(ins->outputStep))==0  ||  tstep+1 == ins->NtimeSteps)
        report(ins, time+ins->dt, tstep+1);

    if (udf.executeStep) udf.executeStep(ins, time+ins->dt, tstep+1);

    if (mesh->rank==0 && (tstep+1)%5==0) fflush(stdout);
  }
}
