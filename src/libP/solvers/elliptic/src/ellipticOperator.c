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

#include "elliptic.h"
#include "ogsInterface.h"
#include <iostream>

#include "omp.h"
void ellipticAx(elliptic_t* elliptic,
                dlong NelementsList,
                occa::memory &o_elementsList,
                occa::memory &o_q,
                occa::memory &o_Aq,
                const char* precision)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide &options = elliptic->options;

  const int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  const int serial = options.compareArgs("THREAD MODEL", "SERIAL");
  const int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");
  const int mapType = (elliptic->elementType == HEXAHEDRA &&
                       options.compareArgs("ELEMENT MAP", "TRILINEAR")) ? 1:0;
  const int integrationType = (elliptic->elementType == HEXAHEDRA &&
                               options.compareArgs("ELLIPTIC INTEGRATION", "CUBATURE")) ? 1:0;

  {
    bool valid = true;
    valid &= continuous;
    if(!strstr(precision, dfloatString)){
      valid &= !elliptic->var_coeff;
      valid &= !elliptic->blockSolver;
      if(!serial){
        valid &= mapType == 0;
        valid &= integrationType == 0;
      }
    }
    if(!valid){
      printf("Encountered invalid configuration inside ellipticAx!\n");
      if(elliptic->var_coeff)
        printf("Precision level (%s) does not support variable coefficient\n", precision);
      if(elliptic->blockSolver)
        printf("Precision level (%s) does not support block solver\n", precision);
      if(!serial){
        if(mapType != 0)
          printf("Precision level (%s) does not support mapType %d\n", precision, mapType);
        if(integrationType != 0)
          printf("Precision level (%s) does not support integrationType %d\n", precision, integrationType);
      }
      exit(1);
    }
  }

  if(serial) {
    if(continuous) {
      if(elliptic->var_coeff) {
        if(elliptic->blockSolver){
          occa::memory & o_geom_factors = elliptic->stressForm ? mesh->o_vgeo : mesh->o_ggeo;
          if(!elliptic->stressForm){
            elliptic->AxKernel(mesh->Nelements, elliptic->Ntotal, elliptic->loffset, o_geom_factors,
                               mesh->o_Dmatrices, mesh->o_Smatrices, elliptic->o_lambda,
                               o_q, o_Aq);
          } else {
            elliptic->AxStressKernel(mesh->Nelements, elliptic->Ntotal, elliptic->loffset, o_geom_factors,
                               mesh->o_Dmatrices, mesh->o_Smatrices, elliptic->o_lambda,
                               o_q, o_Aq);
          }
        }
        else
          elliptic->AxKernel(mesh->Nelements, elliptic->Ntotal, mesh->o_ggeo, mesh->o_Dmatrices,
                             mesh->o_Smatrices, elliptic->o_lambda, o_q, o_Aq);
      }else{
        const dfloat lambda = elliptic->lambda[0];
        if(elliptic->blockSolver){
          occa::memory & o_geom_factors = elliptic->stressForm ? mesh->o_vgeo : mesh->o_ggeo;
          if(!elliptic->stressForm){
            elliptic->AxKernel(mesh->Nelements, elliptic->Ntotal, elliptic->loffset, o_geom_factors,
                               mesh->o_Dmatrices, mesh->o_Smatrices, elliptic->o_lambda,
                               o_q, o_Aq);
          } else {
            elliptic->AxStressKernel(mesh->Nelements, elliptic->Ntotal, elliptic->loffset, o_geom_factors,
                               mesh->o_Dmatrices, mesh->o_Smatrices, elliptic->o_lambda,
                               o_q, o_Aq);
          }
        }
        else{
          occa::memory &o_ggeo = (!strstr(precision,dfloatString)) ? mesh->o_ggeoPfloat : mesh->o_ggeo;
          occa::memory &o_Dmatrices = (!strstr(precision,dfloatString)) ? mesh->o_DmatricesPfloat : mesh->o_Dmatrices;
          occa::memory &o_Smatrices = (!strstr(precision,dfloatString)) ? mesh->o_SmatricesPfloat : mesh->o_Smatrices;
          occa::kernel &AxKernel = (!strstr(precision,dfloatString)) ? elliptic->AxPfloatKernel : elliptic->AxKernel;
          AxKernel(mesh->Nelements, o_ggeo, o_Dmatrices, o_Smatrices, elliptic->lambda[0],
                   o_q, o_Aq);
        }
      }
    } else {
      exit(1);
    }
    return;
  }

  if(continuous) {
    occa::kernel &partialAxKernel =
      (!strstr(precision, dfloatString)) ? elliptic->partialAxPfloatKernel : elliptic->partialAxKernel;

    if(NelementsList) {
      if(integrationType == 0) { // GLL or non-hex
        if(mapType == 0) {
          if(elliptic->var_coeff) {
            if(elliptic->blockSolver){
              occa::memory & o_geom_factors = elliptic->stressForm ? mesh->o_vgeo : mesh->o_ggeo;
              partialAxKernel(NelementsList,
                              elliptic->Ntotal,
                              elliptic->loffset,
                              o_elementsList,
                              o_geom_factors,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
            }
            else
              partialAxKernel(NelementsList,
                              elliptic->Ntotal,
                              o_elementsList,
                              mesh->o_ggeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
          }else{
            if(elliptic->blockSolver){
              occa::memory & o_geom_factors = elliptic->stressForm ? mesh->o_vgeo : mesh->o_ggeo;
              partialAxKernel(NelementsList,
                              elliptic->Ntotal,
                              elliptic->loffset,
                              o_elementsList,
                              o_geom_factors,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
            }
            else{
              occa::memory &o_ggeo = (!strstr(precision,dfloatString)) ? mesh->o_ggeoPfloat : mesh->o_ggeo;
              occa::memory &o_Dmatrices = (!strstr(precision,dfloatString)) ? mesh->o_DmatricesPfloat : mesh->o_Dmatrices;
              occa::memory &o_Smatrices = (!strstr(precision,dfloatString)) ? mesh->o_SmatricesPfloat : mesh->o_Smatrices;
              partialAxKernel(NelementsList,
                              o_elementsList,
                              o_ggeo,
                              o_Dmatrices,
                              o_Smatrices,
                              elliptic->lambda[0],
                              o_q,
                              o_Aq);
            }
          }
        }else{
          if(elliptic->var_coeff) {
            if(elliptic->blockSolver)
              printf("Trilinear version for block solver is not avalibale yet\n");
            else
              partialAxKernel(NelementsList,
                              elliptic->Ntotal,
                              o_elementsList,
                              elliptic->o_EXYZ,
                              elliptic->o_gllzw,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
          }else{
            if(elliptic->blockSolver)
              printf("Trilinear version for block solver is not avalibale yet\n");
            else
              partialAxKernel(NelementsList,
                              o_elementsList,
                              elliptic->o_EXYZ,
                              elliptic->o_gllzw,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              elliptic->lambda[0],
                              o_q,
                              o_Aq);
          }
        }
      }
    }
  }
}

void ellipticOperator(elliptic_t* elliptic,
                      occa::memory &o_q,
                      occa::memory &o_Aq,
                      const char* precision)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide &options = elliptic->options;
  oogs_t* oogsAx = elliptic->oogsAx;
  const char* ogsDataTypeString = (!strstr(precision, dfloatString)) ? 
    options.compareArgs("ENABLE FLOATCOMMHALF GS SUPPORT","TRUE") ? ogsFloatCommHalf : ogsPfloat
    :
    ogsDfloat;
  int serial = options.compareArgs("THREAD MODEL", "SERIAL");
  if(serial) {
    occa::memory o_dummy;
    ellipticAx(elliptic, mesh->Nelements, o_dummy, o_q, o_Aq, precision);
    oogs::startFinish(o_Aq, elliptic->Nfields, elliptic->Ntotal, ogsDataTypeString, ogsAdd, oogsAx);
  } else {
    ellipticAx(elliptic, mesh->NglobalGatherElements, mesh->o_globalGatherElementList, o_q, o_Aq, precision);
    oogs::start(o_Aq, elliptic->Nfields, elliptic->Ntotal, ogsDataTypeString, ogsAdd, oogsAx);
    ellipticAx(elliptic, mesh->NlocalGatherElements, mesh->o_localGatherElementList, o_q, o_Aq, precision);
    oogs::finish(o_Aq, elliptic->Nfields, elliptic->Ntotal, ogsDataTypeString, ogsAdd, oogsAx);
  }
  occa::kernel &maskKernel = (!strstr(precision, dfloatString)) ? mesh->maskPfloatKernel : mesh->maskKernel;
  if (elliptic->Nmasked) maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Aq);
}
