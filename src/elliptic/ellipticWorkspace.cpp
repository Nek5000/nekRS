#include "elliptic.h"
#include "ellipticPrecon.h"

void ellipticAllocateWorkspace(elliptic_t *elliptic)
{
  const auto Nlocal = elliptic->Nfields * static_cast<size_t>(elliptic->fieldOffset);

  elliptic->o_rPfloat = platform->deviceMemoryPool.reserve<pfloat>(Nlocal);
  elliptic->o_zPfloat = platform->deviceMemoryPool.reserve<pfloat>(Nlocal);

  if (elliptic->precon) {
    if (elliptic->precon->MGSolver) {
      elliptic->precon->MGSolver->allocateWorkStorage();
    }
  }
}

void ellipticFreeWorkspace(elliptic_t *elliptic)
{
  if (elliptic->precon) {
    if (elliptic->precon->MGSolver) {
      elliptic->precon->MGSolver->freeWorkStorage();
    }
  }

  elliptic->o_rPfloat.free();
  elliptic->o_zPfloat.free();
}
