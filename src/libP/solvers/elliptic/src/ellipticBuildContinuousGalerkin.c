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

// compare on global indices
static int parallelCompareRowColumn(const void* a, const void* b)
{
  nonZero_t* fa = (nonZero_t*) a;
  nonZero_t* fb = (nonZero_t*) b;

  if(fa->row < fb->row) return -1;
  if(fa->row > fb->row) return +1;

  if(fa->col < fb->col) return -1;
  if(fa->col > fb->col) return +1;

  return 0;
}

void ellipticBuildContinuousGalerkinHex3D (elliptic_t* elliptic,
                                           elliptic_t* ellipticFine,
                                           dfloat lambda,
                                           nonZero_t** A,
                                           dlong* nnz,
                                           ogs_t** ogs,
                                           hlong* globalStarts);

void ellipticBuildContinuousGalerkin(elliptic_t* elliptic,
                                     elliptic_t* ellipticFine,
                                     dfloat lambda,
                                     nonZero_t** A,
                                     dlong* nnz,
                                     ogs_t** ogs,
                                     hlong* globalStarts)
{
  if(elliptic->mesh->Nq != 2) {
    fprintf(stderr,"Coarse level order must equal 1 to use the Galerkin coarse system.");
    exit(1);
  }

  mesh_t* mesh = elliptic->mesh;
  switch(elliptic->elementType) {
  case TRIANGLES:
  case TETRAHEDRA:
  case QUADRILATERALS:
    printf("ellipticBuildContinuousGalerkin is only supported for HEX.\n");
    exit(1);
    break;
  case HEXAHEDRA:
    ellipticBuildContinuousGalerkinHex3D(elliptic,ellipticFine,
                                         lambda,A,nnz,ogs,globalStarts);
    break;
  default:
    break;
  }
}

void ellipticGenerateCoarseBasisHex3D(dfloat* b,int j_,elliptic_t* elliptic)
{
  // b = lx1**dim
  mesh_t* mesh = elliptic->mesh;
  dfloat* z = mesh->gllz;

  int jj = j_ + 1;
  dfloat* zr,* zs,* zt,* z0,* z1;
  zr = (dfloat*)calloc(mesh->Nq,sizeof(dfloat));
  zs = (dfloat*)calloc(mesh->Nq,sizeof(dfloat));
  zt = (dfloat*)calloc(mesh->Nq,sizeof(dfloat));
  z0 = (dfloat*)calloc(mesh->Nq,sizeof(dfloat));
  z1 = (dfloat*)calloc(mesh->Nq,sizeof(dfloat));

  for(int i = 0; i < mesh->Nq; i++) {
    z0[i] = 0.5 * (1 - z[i]);
    z1[i] = 0.5 * (1 + z[i]);
  }

  memcpy(zr,z0,mesh->Nq * sizeof(dfloat));
  memcpy(zs,z0,mesh->Nq * sizeof(dfloat));
  memcpy(zt,z0,mesh->Nq * sizeof(dfloat));

  if(jj % 2 == 0) memcpy(zr,z1,mesh->Nq * sizeof(dfloat));
  if(jj == 3 || jj == 4 || jj == 7 || jj == 8) memcpy(zs,z1,mesh->Nq * sizeof(dfloat));
  if(jj > 4) memcpy(zt,z1,mesh->Nq * sizeof(dfloat));

  int dim = mesh->dim;
  for(int k = 0; k < mesh->Nq; k++)
    for(int j = 0; j < mesh->Nq; j++)
      for(int i = 0; i < mesh->Nq; i++) {
        int n = i + mesh->Nq * j + mesh->Nq * mesh->Nq * k + j_ * mesh->Np;
        b[n] = zr[i] * zs[j] * zt[k];
      }
  free(zr);
  free(zs);
  free(zt);
  free(z0);
  free(z1);
}

void ellipticBuildContinuousGalerkinHex3D(elliptic_t* elliptic,
                                          elliptic_t* ellipticFine,
                                          dfloat lambda,
                                          nonZero_t** A,
                                          dlong* nnz,
                                          ogs_t** ogs,
                                          hlong* globalStarts)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int rank = mesh->rank;

  //use the masked gs handle to define a global ordering

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = elliptic->ogs->Ngather;
  dlong Ntotal  = mesh->Np * mesh->Nelements;

  // create a global numbering system
  hlong* globalIds = (hlong*) calloc(Ngather,sizeof(hlong));
  int* owner     = (int*) calloc(Ngather,sizeof(int));

  // every gathered degree of freedom has its own global id
  MPI_Allgather(&Ngather, 1, MPI_HLONG, globalStarts + 1, 1, MPI_HLONG, mesh->comm);
  for(int r = 0; r < mesh->size; ++r)
    globalStarts[r + 1] = globalStarts[r] + globalStarts[r + 1];

  //use the offsets to set a consecutive global numbering
  for (dlong n = 0; n < elliptic->ogs->Ngather; n++) {
    globalIds[n] = n + globalStarts[rank];
    owner[n] = rank;
  }

  //scatter this numbering to the original nodes
  hlong* globalNumbering = (hlong*) calloc(Ntotal,sizeof(hlong));
  int* globalOwners = (int*) calloc(Ntotal,sizeof(int));
  for (dlong n = 0; n < Ntotal; n++) globalNumbering[n] = -1;
  ogsScatter(globalNumbering, globalIds, ogsHlong, ogsAdd, elliptic->ogs);
  ogsScatter(globalOwners, owner, ogsInt, ogsAdd, elliptic->ogs);

  free(globalIds);
  free(owner);

  // 2. Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = mesh->Np * mesh->Np * mesh->Nelements;
  nonZero_t* sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  int* AsendCounts  = (int*) calloc(mesh->size, sizeof(int));
  int* ArecvCounts  = (int*) calloc(mesh->size, sizeof(int));
  int* AsendOffsets = (int*) calloc(mesh->size + 1, sizeof(int));
  int* ArecvOffsets = (int*) calloc(mesh->size + 1, sizeof(int));

  int* mask = (int*) calloc(mesh->Np * mesh->Nelements,sizeof(int));
  for (dlong n = 0; n < elliptic->Nmasked; n++) mask[elliptic->maskIds[n]] = 1;

  if(mesh->rank == 0) printf("Building full FEM matrix via Galerkin projection...");
  fflush(stdout);

  mesh_t* meshf = ellipticFine->mesh;

  dfloat* b,* q,* Aq;
  b = (dfloat*)calloc(meshf->Np * mesh->Np,sizeof(dfloat));
  q = (dfloat*)calloc(meshf->Np * mesh->Nelements,sizeof(dfloat));
  Aq = (dfloat*)calloc(meshf->Np * mesh->Nelements,sizeof(dfloat));

  occa::memory o_q = mesh->device.malloc(meshf->Np *
                                         mesh->Nelements * sizeof(dfloat), q);
  occa::memory o_Aq = mesh->device.malloc(meshf->Np *
                                          mesh->Nelements * sizeof(dfloat),Aq);

  for(int jj = 0; jj < mesh->Np; jj++)
    ellipticGenerateCoarseBasisHex3D(b,jj,ellipticFine);

  int enableGatherScatters = 1;
  int enableReductions = 1;

  setupAide &optionsFine = ellipticFine->options;
  optionsFine.getArgs("DEBUG ENABLE REDUCTIONS",
                      enableReductions);
  optionsFine.getArgs("DEBUG ENABLE OGS",
                      enableGatherScatters);

  dlong cnt = 0;
  for (int nz = 0; nz < mesh->Nq; nz++)
    for (int ny = 0; ny < mesh->Nq; ny++)
      for (int nx = 0; nx < mesh->Nq; nx++) {
        int idn = nx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
        // ok
        for (dlong e = 0; e < mesh->Nelements; e++)
          memcpy(&q[e * meshf->Np],&b[idn * meshf->Np],
                 meshf->Np * sizeof(dfloat));

        o_q.copyFrom(q);
        ellipticFine->AxKernel(mesh->Nelements,meshf->o_ggeo,
                               meshf->o_Dmatrices,meshf->o_Smatrices,
                               lambda,o_q,o_Aq);
        o_Aq.copyTo(Aq);

        for(dlong e = 0; e < mesh->Nelements; e++)
          for (int mz = 0; mz < mesh->Nq; mz++)
            for (int my = 0; my < mesh->Nq; my++)
              for (int mx = 0; mx < mesh->Nq; mx++) {
                int idm = mx + my * mesh->Nq + mz * mesh->Nq * mesh->Nq;
                if (mask[e * mesh->Np + idm])
                  continue; //skip masked nodes
                if (mask[e * mesh->Np + idn])
                  continue; //skip masked nodes

                dfloat val = 0.f;
                for(int mmm = 0; mmm < meshf->Np; mmm++)
                  val += Aq[e * meshf->Np + mmm] * b[mmm + idm * meshf->Np];

                // pack non-zero
                dfloat nonZeroThreshold = 1e-7;
                if (fabs(val) >= nonZeroThreshold) {
                  sendNonZeros[cnt].val = val;
                  sendNonZeros[cnt].row = globalNumbering[e * mesh->Np + idm];
                  sendNonZeros[cnt].col = globalNumbering[e * mesh->Np + idn];
                  sendNonZeros[cnt].ownerRank = globalOwners[e * mesh->Np + idm];
                  cnt++;
                }
              }
      }

  free(b);
  free(q);
  free(Aq);
  o_q.free();
  o_Aq.free();

  // Make the MPI_NONZERO_T data type
  MPI_Datatype MPI_NONZERO_T;
  MPI_Datatype dtype[4] = {MPI_HLONG, MPI_HLONG, MPI_INT, MPI_DFLOAT};
  int blength[4] = {1, 1, 1, 1};
  MPI_Aint addr[4], displ[4];
  MPI_Get_address ( &(sendNonZeros[0]          ), addr + 0);
  MPI_Get_address ( &(sendNonZeros[0].col      ), addr + 1);
  MPI_Get_address ( &(sendNonZeros[0].ownerRank), addr + 2);
  MPI_Get_address ( &(sendNonZeros[0].val      ), addr + 3);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  MPI_Type_create_struct (4, blength, displ, dtype, &MPI_NONZERO_T);
  MPI_Type_commit (&MPI_NONZERO_T);

  // count how many non-zeros to send to each process
  for(dlong n = 0; n < cnt; ++n)
    AsendCounts[sendNonZeros[n].ownerRank]++;

  // sort by row ordering
  qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_INT, ArecvCounts, 1, MPI_INT, mesh->comm);

  // find send and recv offsets for gather
  *nnz = 0;
  for(int r = 0; r < mesh->size; ++r) {
    AsendOffsets[r + 1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r + 1] = ArecvOffsets[r] + ArecvCounts[r];
    *nnz += ArecvCounts[r];
  }

  *A = (nonZero_t*) calloc(*nnz, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_NONZERO_T,
                (*A), ArecvCounts, ArecvOffsets, MPI_NONZERO_T,
                mesh->comm);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), *nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(dlong n = 1; n < *nnz; ++n) {
    if((*A)[n].row == (*A)[cnt].row &&
       (*A)[n].col == (*A)[cnt].col) {
      (*A)[cnt].val += (*A)[n].val;
    }else {
      ++cnt;
      (*A)[cnt] = (*A)[n];
    }
  }
  if (*nnz) cnt++;
  *nnz = cnt;

  if(mesh->rank == 0) printf("done.\n");

  MPI_Barrier(mesh->comm);
  MPI_Type_free(&MPI_NONZERO_T);

  free(sendNonZeros);
  free(globalNumbering);
  free(globalOwners);

  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);
}
