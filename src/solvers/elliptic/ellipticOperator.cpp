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

#include <type_traits>
#include "elliptic.h"
#include <iostream>

void ellipticAx(elliptic_t* elliptic,
                dlong NelementsList,
                occa::memory &o_elementsList,
                occa::memory &o_q,
                occa::memory &o_Aq,
                const char* precision)
{

  if(NelementsList == 0) return;

  mesh_t* mesh = elliptic->mesh;
  setupAide &options = elliptic->options;

  const bool coeffField = options.compareArgs("ELLIPTIC COEFF FIELD", "TRUE");
  const bool continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  const int mapType = (elliptic->elementType == HEXAHEDRA &&
                       options.compareArgs("ELEMENT MAP", "TRILINEAR")) ? 1:0;
  const std::string precisionStr(precision);
  const std::string dFloatStr(dfloatString);


  bool valid = true;
  valid &= continuous;
  if(precisionStr != dFloatStr) {
    valid &= !elliptic->blockSolver;
    valid &= mapType == 0;
  }

  auto errTxt = [&]()
  {
    std::stringstream txt;
    txt << "Encountered invalid configuration inside ellipticAx!\n";
    if(elliptic->blockSolver)
      txt << "Precision level (" << precision << ") does not support block solver\n";
    if(mapType != 0)
      txt << "Precision level (" << precision << ") does not support mapType " << mapType << "\n";

    return txt.str();
  };
  nrsCheck(!valid, MPI_COMM_SELF, EXIT_FAILURE, errTxt().c_str(), "");

  occa::memory & o_geom_factors =
    (precisionStr != dFloatStr) ? mesh->o_ggeoPfloat :
      elliptic->stressForm ? mesh->o_vgeo : mesh->o_ggeo;
  occa::memory & o_D = (precisionStr != dFloatStr) ? mesh->o_DPfloat : mesh->o_D;
  occa::memory & o_DT = (precisionStr != dFloatStr) ? mesh->o_DTPfloat : mesh->o_DT;
  occa::memory & o_lambda0 = elliptic->o_lambda0;
  occa::memory & o_lambda1 = elliptic->o_lambda1;

  occa::kernel &AxKernel =
      (precisionStr != dFloatStr) ? elliptic->AxPfloatKernel : elliptic->AxKernel;

  AxKernel(NelementsList,
           elliptic->fieldOffset,
           elliptic->loffset,
           o_elementsList,
           o_geom_factors,
           o_D,
           o_DT,
           o_lambda0,
           o_lambda1,
           o_q,
           o_Aq);

  double flopCount = mesh->Np * 12 * mesh->Nq + 15 * mesh->Np;
  if(coeffField)
    flopCount += 3 * mesh->Np;
  else 
    flopCount += 1 * mesh->Np;

  if(!elliptic->poisson)
    flopCount += (2 + 1) * mesh->Np;

  if (elliptic->stressForm)
    flopCount += (15 + 6) * mesh->Np;

  flopCount *= elliptic->Nfields * static_cast<double>(NelementsList);


  const double factor = std::is_same<pfloat, float>::value && (precisionStr != dFloatStr) ? 0.5 : 1.0;

  platform->flopCounter->add(elliptic->name + " Ax, N=" + std::to_string(mesh->N) + ", " +
                             std::string(precision),
                             factor * flopCount);
}

void ellipticOperator(elliptic_t* elliptic,
                      occa::memory &o_q,
                      occa::memory &o_Aq,
                      const char* precision,
                      bool masked)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide &options = elliptic->options;
  oogs_t* oogs = elliptic->oogsAx;
  const bool overlap = (oogs != elliptic->oogs);
  const char* ogsDataTypeString = (!strstr(precision, dfloatString)) ? ogsPfloat : ogsDfloat;

  if(overlap) {

    ellipticAx(elliptic, mesh->NglobalGatherElements, mesh->o_globalGatherElementList, o_q, o_Aq, precision);
    if (masked) {
      ellipticApplyMask(elliptic,
                        mesh->NglobalGatherElements,
                        elliptic->NmaskedGlobal,
                        mesh->o_globalGatherElementList,
                        elliptic->o_maskIdsGlobal,
                        o_Aq,
                        precision);
    }
 
    oogs::start(o_Aq, elliptic->Nfields, elliptic->fieldOffset, ogsDataTypeString, ogsAdd, oogs);
    ellipticAx(elliptic, mesh->NlocalGatherElements, mesh->o_localGatherElementList, o_q, o_Aq, precision);
 
    if (masked) {
      ellipticApplyMask(elliptic,
                        mesh->NlocalGatherElements,
                        elliptic->NmaskedLocal,
                        mesh->o_localGatherElementList,
                        elliptic->o_maskIdsLocal,
                        o_Aq,
                        precision);
    }

    oogs::finish(o_Aq, elliptic->Nfields, elliptic->fieldOffset, ogsDataTypeString, ogsAdd, oogs);

  } else {

    ellipticAx(elliptic, mesh->Nelements, mesh->o_elementList, o_q, o_Aq, precision);
    if (masked) {
      ellipticApplyMask(elliptic,
                        mesh->Nelements,
                        elliptic->Nmasked,
                        mesh->o_elementList,
                        elliptic->o_maskIds,
                        o_Aq,
                        precision);
    }
    oogs::startFinish(o_Aq, elliptic->Nfields, elliptic->fieldOffset, ogsDataTypeString, ogsAdd, oogs);

  }
}
