#if !defined(nekrs_meshsetup_hpp_)
#define nekrs_meshsetup_hpp_

#include "nekrs.hpp"
mesh_t *createMeshDummy(MPI_Comm comm, int N);
mesh_t *createMeshT(MPI_Comm comm, int N, int isMeshT);
mesh_t *createMeshV(MPI_Comm comm, int N, mesh_t *meshT);
void meshVOccaSetup3D(mesh_t *mesh, setupAide &options, occa::properties &kernelInfo);

#endif
