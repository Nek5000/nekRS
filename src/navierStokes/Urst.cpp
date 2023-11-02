#include "Urst.hpp"
#include "nrs.hpp"

void computeUrst(nrs_t *nrs, bool relative, bool cubature)
{
  occa::memory &o_Urst = relative ? nrs->o_relUrst : nrs->o_Urst;
  auto mesh = nrs->meshV;
  double flopCount = 0.0;

  if (cubature) {
    nrs->UrstCubatureKernel(mesh->Nelements,
                            static_cast<int>(relative),
                            mesh->o_cubvgeo,
                            mesh->o_cubInterpT,
                            nrs->fieldOffset,
                            nrs->cubatureOffset,
                            nrs->o_U,
                            mesh->o_U,
                            o_Urst);
    flopCount += 6 * mesh->Np * mesh->cubNq;
    flopCount += 6 * mesh->Nq * mesh->Nq * mesh->cubNq * mesh->cubNq;
    flopCount += 6 * mesh->Nq * mesh->cubNp;
    flopCount += 24 * mesh->cubNp;
    flopCount *= mesh->Nelements;
  }
  else {
    nrs->UrstKernel(mesh->Nelements,
                    static_cast<int>(relative),
                    mesh->o_vgeo,
                    nrs->fieldOffset,
                    nrs->o_U,
                    mesh->o_U,
                    o_Urst);
    flopCount += 24 * static_cast<double>(mesh->Nlocal);
  }
  platform->flopCounter->add("Urst", flopCount);

}
