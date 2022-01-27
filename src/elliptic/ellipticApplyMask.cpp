#include <elliptic.h>
#include <ellipticApplyMask.hpp>
void applyMask(elliptic_t *solver, occa::memory &o_x, std::string precision)
{
  mesh_t *mesh = solver->mesh;
  occa::kernel &maskKernel = (precision != dfloatString) ? mesh->maskPfloatKernel : mesh->maskKernel;
  occa::kernel &enforceUnKernel =
      (precision != dfloatString) ? solver->enforceUnPfloatKernel : solver->enforceUnKernel;

  if (solver->UNormalZero) {
    dlong Nelems = mesh->NlocalGatherElements;

    if (Nelems > 0) {
      occa::memory &o_elemList = solver->mesh->o_localGatherElementList;
      enforceUnKernel(Nelems,
                      solver->Ntotal,
                      o_elemList,
                      mesh->o_VT1,
                      mesh->o_VT2,
                      mesh->o_vmapM,
                      mesh->o_EToB,
                      solver->o_BCType,
                      o_x);
    }

    Nelems = mesh->NglobalGatherElements;
    if (Nelems > 0) {
      occa::memory &o_elemList = solver->mesh->o_globalGatherElementList;
      enforceUnKernel(Nelems,
                      solver->Ntotal,
                      o_elemList,
                      mesh->o_VT1,
                      mesh->o_VT2,
                      mesh->o_vmapM,
                      mesh->o_EToB,
                      solver->o_BCType,
                      o_x);
    }
  }

  const dlong Nmasked = solver->Nmasked;
  occa::memory &o_maskIds = solver->o_maskIds;
  if (Nmasked)
    maskKernel(Nmasked, o_maskIds, o_x);
}
void applyMaskInterior(elliptic_t *solver, occa::memory &o_x, std::string precision)
{
  mesh_t *mesh = solver->mesh;
  occa::kernel &maskKernel = (precision != dfloatString) ? mesh->maskPfloatKernel : mesh->maskKernel;
  occa::kernel &enforceUnKernel =
      (precision != dfloatString) ? solver->enforceUnPfloatKernel : solver->enforceUnKernel;

  if (solver->UNormalZero) {
    dlong Nelems = mesh->NlocalGatherElements;

    if (Nelems > 0) {
      occa::memory &o_elemList = solver->mesh->o_localGatherElementList;
      enforceUnKernel(Nelems,
                      solver->Ntotal,
                      o_elemList,
                      mesh->o_VT1,
                      mesh->o_VT2,
                      mesh->o_vmapM,
                      mesh->o_EToB,
                      solver->o_BCType,
                      o_x);
    }
  }

  const dlong Nmasked = solver->NmaskedLocal;
  occa::memory &o_maskIds = solver->o_maskIdsLocal;
  if (Nmasked)
    maskKernel(Nmasked, o_maskIds, o_x);
}

void applyMaskExterior(elliptic_t *solver, occa::memory &o_x, std::string precision)
{
  mesh_t *mesh = solver->mesh;
  occa::kernel &maskKernel = (precision != dfloatString) ? mesh->maskPfloatKernel : mesh->maskKernel;
  occa::kernel &enforceUnKernel =
      (precision != dfloatString) ? solver->enforceUnPfloatKernel : solver->enforceUnKernel;

  const dlong Nmasked = solver->NmaskedGlobal;
  occa::memory &o_maskIds = solver->o_maskIdsGlobal;

  if (solver->UNormalZero) {
    const dlong Nelems = mesh->NglobalGatherElements;

    if (Nelems > 0) {
      occa::memory &o_elemList = solver->mesh->o_globalGatherElementList;
      enforceUnKernel(Nelems,
                      solver->Ntotal,
                      o_elemList,
                      mesh->o_VT1,
                      mesh->o_VT2,
                      mesh->o_vmapM,
                      mesh->o_EToB,
                      solver->o_BCType,
                      o_x);
    }
  }
  if (Nmasked)
    maskKernel(Nmasked, o_maskIds, o_x);
}