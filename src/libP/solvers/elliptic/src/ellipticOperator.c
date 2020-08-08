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

#include "omp.h"

//host gs
void ellipticSerialOperator(elliptic_t* elliptic,
                            occa::memory &o_q,
                            occa::memory &o_Aq,
                            const char* precision)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide &options = elliptic->options;

  int enableGatherScatters = 1;
  int enableReductions = 1;
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");

  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);
  options.getArgs("DEBUG ENABLE OGS", enableGatherScatters);

  dfloat* sendBuffer = elliptic->sendBuffer;
  dfloat* recvBuffer = elliptic->recvBuffer;
  dfloat* gradSendBuffer = elliptic->gradSendBuffer;
  dfloat* gradRecvBuffer = elliptic->gradRecvBuffer;

  dfloat alpha = 0., alphaG = 0.;

  if(continuous) {
    if(elliptic->var_coeff) {
      if(elliptic->blockSolver)
        elliptic->AxKernel(mesh->Nelements, elliptic->Ntotal, elliptic->loffset, mesh->o_ggeo,
                           mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, elliptic->o_lambda,
                           o_q, o_Aq);
      else
        elliptic->AxKernel(mesh->Nelements, elliptic->Ntotal, mesh->o_ggeo, mesh->o_Dmatrices,
                           mesh->o_Smatrices, mesh->o_MM, elliptic->o_lambda, o_q, o_Aq);
    }else{
      const dfloat lambda = elliptic->lambda[0];
      if(elliptic->blockSolver)
        elliptic->AxKernel(mesh->Nelements, elliptic->Ntotal, elliptic->loffset, mesh->o_ggeo,
                           mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM, elliptic->o_lambda,
                           o_q, o_Aq);
      else
        elliptic->AxKernel(mesh->Nelements,
                           mesh->o_ggeo,
                           mesh->o_Dmatrices,
                           mesh->o_Smatrices,
                           mesh->o_MM,
                           lambda,
                           o_q,
                           o_Aq);
    }

    oogs::startFinish(o_Aq, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);
    if (elliptic->Nmasked)
      mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Aq);

  } else if(ipdg) {

    printf("WARNING: DEBUGGING C0\n");
    MPI_Finalize();

    exit(-1);
  }
}

void ellipticOperator(elliptic_t* elliptic,
                      occa::memory &o_q,
                      occa::memory &o_Aq,
                      const char* precision)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide &options = elliptic->options;

  int enableGatherScatters = 1;
  int enableReductions = 1;
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int serial = options.compareArgs("THREAD MODEL", "SERIAL");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");

  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);
  options.getArgs("DEBUG ENABLE OGS", enableGatherScatters);

  dfloat* sendBuffer = elliptic->sendBuffer;
  dfloat* recvBuffer = elliptic->recvBuffer;
  dfloat* gradSendBuffer = elliptic->gradSendBuffer;
  dfloat* gradRecvBuffer = elliptic->gradRecvBuffer;

  dfloat alpha = 0., alphaG = 0.;
  dlong Nblock = elliptic->Nblock;
  dfloat* tmp = elliptic->tmp;
  occa::memory &o_tmp = elliptic->o_tmp;

  if(serial) {
    ellipticSerialOperator(elliptic, o_q, o_Aq, precision);
    return;
  }

  if(continuous) {
    oogs_t* oogsAx = elliptic->oogsAx;

    int mapType = (elliptic->elementType == HEXAHEDRA &&
                   options.compareArgs("ELEMENT MAP", "TRILINEAR")) ? 1:0;

    int integrationType = (elliptic->elementType == HEXAHEDRA &&
                           options.compareArgs("ELLIPTIC INTEGRATION", "CUBATURE")) ? 1:0;

    occa::kernel &partialAxKernel =
      (strstr(precision, "float")) ? elliptic->partialFloatAxKernel : elliptic->partialAxKernel;

    if(mesh->NglobalGatherElements) {
      if(integrationType == 0) { // GLL or non-hex
        if(mapType == 0) {
          if(elliptic->var_coeff) {
            if(elliptic->blockSolver)
              partialAxKernel(mesh->NglobalGatherElements,
                              elliptic->Ntotal,
                              elliptic->loffset,
                              mesh->o_globalGatherElementList,
                              mesh->o_ggeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
            else
              partialAxKernel(mesh->NglobalGatherElements,
                              elliptic->Ntotal,
                              mesh->o_globalGatherElementList,
                              mesh->o_ggeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
          }else{
            if(elliptic->blockSolver)
              partialAxKernel(mesh->NglobalGatherElements,
                              elliptic->Ntotal,
                              elliptic->loffset,
                              mesh->o_globalGatherElementList,
                              mesh->o_ggeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
            else
              partialAxKernel(mesh->NglobalGatherElements,
                              mesh->o_globalGatherElementList,
                              mesh->o_ggeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->lambda[0],
                              o_q,
                              o_Aq);
          }
        }else{
          if(elliptic->var_coeff) {
            if(elliptic->blockSolver)
              printf("Trilinear version for block solver is not avalibale yet\n");
            else
              partialAxKernel(mesh->NglobalGatherElements,
                              elliptic->Ntotal,
                              mesh->o_globalGatherElementList,
                              elliptic->o_EXYZ,
                              elliptic->o_gllzw,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
          }else{
            if(elliptic->blockSolver)
              printf("Trilinear version for block solver is not avalibale yet\n");
            else
              partialAxKernel(mesh->NglobalGatherElements,
                              mesh->o_globalGatherElementList,
                              elliptic->o_EXYZ,
                              elliptic->o_gllzw,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->lambda[0],
                              o_q,
                              o_Aq);
          }
        }
      }
    }

    if(enableGatherScatters)
      oogs::start(o_Aq, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, oogsAx);

    if(mesh->NlocalGatherElements) {
      if(integrationType == 0) { // GLL or non-hex
        if(mapType == 0) {
          if(elliptic->var_coeff) {
            if(elliptic->blockSolver)
              partialAxKernel(mesh->NlocalGatherElements,
                              elliptic->Ntotal,
                              elliptic->loffset,
                              mesh->o_localGatherElementList,
                              mesh->o_ggeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
            else
              partialAxKernel(mesh->NlocalGatherElements,
                              elliptic->Ntotal,
                              mesh->o_localGatherElementList,
                              mesh->o_ggeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
          }else{
            if(elliptic->blockSolver)
              partialAxKernel(mesh->NlocalGatherElements,
                              elliptic->Ntotal,
                              elliptic->loffset,
                              mesh->o_localGatherElementList,
                              mesh->o_ggeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
            else
              partialAxKernel(mesh->NlocalGatherElements,
                              mesh->o_localGatherElementList,
                              mesh->o_ggeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->lambda[0],
                              o_q,
                              o_Aq);
          }
        }else {
          if(elliptic->var_coeff) {
            if(elliptic->blockSolver)
              printf("Trilinear version for block solver is not avalibale yet\n");
            else
              partialAxKernel(mesh->NlocalGatherElements,
                              elliptic->Ntotal,
                              mesh->o_localGatherElementList,
                              elliptic->o_EXYZ,
                              elliptic->o_gllzw,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->o_lambda,
                              o_q,
                              o_Aq);
          }else{
            if(elliptic->blockSolver)
              printf("Trilinear version for block solver is not avalibale yet\n");
            else
              partialAxKernel(mesh->NlocalGatherElements,
                              mesh->o_localGatherElementList,
                              elliptic->o_EXYZ,
                              elliptic->o_gllzw,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              elliptic->lambda[0],
                              o_q,
                              o_Aq);
          }
        }
      }
    }

    if(enableGatherScatters)
      oogs::finish(o_Aq, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, oogsAx);

    if (elliptic->Nmasked)
      mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Aq);

  }
}
