#if !defined(nekrs_meshsetup_hpp_)
#define nekrs_meshsetup_hpp_

#include "nrs.hpp"
void createMeshDummy(mesh_t* mesh,
                     MPI_Comm comm,
                     int N,
                     int cubN,
                     setupAide &options,
                     occa::device device,
                     occa::properties &kernelInfo);

void createMesh(mesh_t* mesh,
                MPI_Comm comm,
                int N,
                int cubN,
                int isMeshT,
                setupAide &options,
                occa::device device,
                occa::properties &kernelInfo);

void createMeshV(mesh_t* mesh,
                 MPI_Comm comm,
                 int N,
                 int cubN,
                 mesh_t* meshT,
                 setupAide &options,
                 occa::properties &kernelInfo);

#endif
