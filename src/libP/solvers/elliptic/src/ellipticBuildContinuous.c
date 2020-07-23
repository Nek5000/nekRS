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
int parallelCompareRowColumn(const void* a, const void* b)
{
  nonZero_t* fa = (nonZero_t*) a;
  nonZero_t* fb = (nonZero_t*) b;

  if(fa->row < fb->row) return -1;
  if(fa->row > fb->row) return +1;

  if(fa->col < fb->col) return -1;
  if(fa->col > fb->col) return +1;

  return 0;
}

// void ellipticBuildContinuousTri2D (elliptic_t *elliptic, dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts);
void ellipticBuildContinuousQuad2D(elliptic_t* elliptic,
                                   nonZero_t** A,
                                   dlong* nnz,
                                   ogs_t** ogs,
                                   hlong* globalStarts);
// void ellipticBuildContinuousQuad3D(elliptic_t *elliptic, dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts);
// void ellipticBuildContinuousTet3D (elliptic_t *elliptic, dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts);
void ellipticBuildContinuousHex3D (elliptic_t* elliptic,
                                   nonZero_t** A,
                                   dlong* nnz,
                                   ogs_t** ogs,
                                   hlong* globalStarts);

void ellipticBuildContinuous(elliptic_t* elliptic,
                             nonZero_t** A,
                             dlong* nnz,
                             ogs_t** ogs,
                             hlong* globalStarts)
{
  switch(elliptic->elementType) {
  // case TRIANGLES:
  //   ellipticBuildContinuousTri2D(elliptic, lambda, A, nnz, ogs, globalStarts); break;
  case QUADRILATERALS: {
    // if(elliptic->dim==2)
    ellipticBuildContinuousQuad2D(elliptic, A, nnz, ogs, globalStarts);
    // else
    // ellipticBuildContinuousQuad3D(elliptic, lambda, A, nnz, ogs, globalStarts);
    break;
  }
  // case TETRAHEDRA:
  //   ellipticBuildContinuousTet3D(elliptic, lambda, A, nnz, ogs, globalStarts); break;
  case HEXAHEDRA:
    ellipticBuildContinuousHex3D(elliptic, A, nnz, ogs, globalStarts);
    break;
  }
}

void ellipticBuildContinuousQuad2D(elliptic_t* elliptic,
                                   nonZero_t** A,
                                   dlong* nnz,
                                   ogs_t** ogs,
                                   hlong* globalStarts)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;
  // currently constant coefficient case only
  const dfloat lambda = elliptic->lambda[0];

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

  if(mesh->rank == 0) printf("Building full FEM matrix...");
  fflush(stdout);

  //Build unassembed non-zeros
  dlong cnt = 0;
  for (dlong e = 0; e < mesh->Nelements; e++)
    for (int ny = 0; ny < mesh->Nq; ny++)
      for (int nx = 0; nx < mesh->Nq; nx++) {
        if (mask[e * mesh->Np + nx + ny * mesh->Nq]) continue; //skip masked nodes
        for (int my = 0; my < mesh->Nq; my++)
          for (int mx = 0; mx < mesh->Nq; mx++) {
            if (mask[e * mesh->Np + mx + my * mesh->Nq]) continue; //skip masked nodes

            int id;
            dfloat val = 0.;

            if (ny == my) {
              for (int k = 0; k < mesh->Nq; k++) {
                id = k + ny * mesh->Nq;
                dfloat Grr = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G00ID * mesh->Np];

                val += Grr * mesh->D[nx + k * mesh->Nq] * mesh->D[mx + k * mesh->Nq];
              }
            }

            id = mx + ny * mesh->Nq;
            dfloat Grs = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G01ID * mesh->Np];
            val += Grs * mesh->D[nx + mx * mesh->Nq] * mesh->D[my + ny * mesh->Nq];

            id = nx + my * mesh->Nq;
            dfloat Gsr = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G01ID * mesh->Np];
            val += Gsr * mesh->D[mx + nx * mesh->Nq] * mesh->D[ny + my * mesh->Nq];

            if (nx == mx) {
              for (int k = 0; k < mesh->Nq; k++) {
                id = nx + k * mesh->Nq;
                dfloat Gss = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G11ID * mesh->Np];

                val += Gss * mesh->D[ny + k * mesh->Nq] * mesh->D[my + k * mesh->Nq];
              }
            }

            if ((nx == mx) && (ny == my)) {
              id = nx + ny * mesh->Nq;
              dfloat JW = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + GWJID * mesh->Np];
              val += JW * lambda;
            }

            dfloat nonZeroThreshold = 1e-7;
            if (fabs(val) > nonZeroThreshold) {
              // pack non-zero
              sendNonZeros[cnt].val = val;
              sendNonZeros[cnt].row = globalNumbering[e * mesh->Np + nx + ny * mesh->Nq];
              sendNonZeros[cnt].col = globalNumbering[e * mesh->Np + mx + my * mesh->Nq];
              sendNonZeros[cnt].ownerRank = globalOwners[e * mesh->Np + nx + ny * mesh->Nq];
              cnt++;
            }
          }
      }

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

#if 1
  // Write matlab dat for postprocess
  char fname[BUFSIZ];
  sprintf(fname, "Ax.dat");
  FILE* fp;
  fp = fopen(fname, "w");

  for(dlong n = 1; n < *nnz; ++n)
    fprintf(fp, hlongFormat " " hlongFormat " %.8e\n", (*A)[n].row + 1, (*A)[n].col + 1,
            (*A)[n].val);

  fclose(fp);
#endif

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

void ellipticBuildContinuousHex3D(elliptic_t* elliptic,
                                  nonZero_t** A,
                                  dlong* nnz,
                                  ogs_t** ogs,
                                  hlong* globalStarts)
{
  mesh2D* mesh = elliptic->mesh;
  setupAide options = elliptic->options;
  // currently constant coefficient case only
  const dfloat lambda = elliptic->lambda[0];

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

  if(mesh->rank == 0) printf("Building full FEM matrix...");
  fflush(stdout);

  dlong cnt = 0;
  for (dlong e = 0; e < mesh->Nelements; e++)
    for (int nz = 0; nz < mesh->Nq; nz++)
      for (int ny = 0; ny < mesh->Nq; ny++)
        for (int nx = 0; nx < mesh->Nq; nx++) {
          int idn = nx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
          if (mask[e * mesh->Np + idn]) continue; //skip masked nodes

          for (int mz = 0; mz < mesh->Nq; mz++)
            for (int my = 0; my < mesh->Nq; my++)
              for (int mx = 0; mx < mesh->Nq; mx++) {
                int idm = mx + my * mesh->Nq + mz * mesh->Nq * mesh->Nq;
                if (mask[e * mesh->Np + idm]) continue; //skip masked nodes

                int id;
                dfloat val = 0.;

                if ((ny == my) && (nz == mz)) {
                  for (int k = 0; k < mesh->Nq; k++) {
                    id = k + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                    dfloat Grr = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G00ID * mesh->Np];

                    val += Grr * mesh->D[nx + k * mesh->Nq] * mesh->D[mx + k * mesh->Nq];
                  }
                }

                if (nz == mz) {
                  id = mx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                  dfloat Grs = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G01ID * mesh->Np];
                  val += Grs * mesh->D[nx + mx * mesh->Nq] * mesh->D[my + ny * mesh->Nq];

                  id = nx + my * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                  dfloat Gsr = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G01ID * mesh->Np];
                  val += Gsr * mesh->D[mx + nx * mesh->Nq] * mesh->D[ny + my * mesh->Nq];
                }

                if (ny == my) {
                  id = mx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                  dfloat Grt = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G02ID * mesh->Np];
                  val += Grt * mesh->D[nx + mx * mesh->Nq] * mesh->D[mz + nz * mesh->Nq];

                  id = nx + ny * mesh->Nq + mz * mesh->Nq * mesh->Nq;
                  dfloat Gst = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G02ID * mesh->Np];
                  val += Gst * mesh->D[mx + nx * mesh->Nq] * mesh->D[nz + mz * mesh->Nq];
                }

                if ((nx == mx) && (nz == mz)) {
                  for (int k = 0; k < mesh->Nq; k++) {
                    id = nx + k * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                    dfloat Gss = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G11ID * mesh->Np];

                    val += Gss * mesh->D[ny + k * mesh->Nq] * mesh->D[my + k * mesh->Nq];
                  }
                }

                if (nx == mx) {
                  id = nx + my * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                  dfloat Gst = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G12ID * mesh->Np];
                  val += Gst * mesh->D[ny + my * mesh->Nq] * mesh->D[mz + nz * mesh->Nq];

                  id = nx + ny * mesh->Nq + mz * mesh->Nq * mesh->Nq;
                  dfloat Gts = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G12ID * mesh->Np];
                  val += Gts * mesh->D[my + ny * mesh->Nq] * mesh->D[nz + mz * mesh->Nq];
                }

                if ((nx == mx) && (ny == my)) {
                  for (int k = 0; k < mesh->Nq; k++) {
                    id = nx + ny * mesh->Nq + k * mesh->Nq * mesh->Nq;
                    dfloat Gtt = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + G22ID * mesh->Np];

                    val += Gtt * mesh->D[nz + k * mesh->Nq] * mesh->D[mz + k * mesh->Nq];
                  }
                }

                if ((nx == mx) && (ny == my) && (nz == mz)) {
                  id = nx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                  dfloat JW = mesh->ggeo[e * mesh->Np * mesh->Nggeo + id + GWJID * mesh->Np];
                  val += JW * lambda;
                }

                // pack non-zero
                dfloat nonZeroThreshold = 1e-7;
                if (fabs(val) >= nonZeroThreshold) {
                  sendNonZeros[cnt].val = val;
                  sendNonZeros[cnt].row = globalNumbering[e * mesh->Np + idn];
                  sendNonZeros[cnt].col = globalNumbering[e * mesh->Np + idm];
                  sendNonZeros[cnt].ownerRank = globalOwners[e * mesh->Np + idn];
                  cnt++;
                }
              }
        }

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

// void ellipticBuildContinuousTri2D(elliptic_t *elliptic, dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts) {

//   mesh2D *mesh = elliptic->mesh;
//   setupAide options = elliptic->options;

//   int rank = mesh->rank;

//   //use the masked gs handle to define a global ordering

//   // number of degrees of freedom on this rank (after gathering)
//   hlong Ngather = elliptic->ogs->Ngather;
//   dlong Ntotal  = mesh->Np*mesh->Nelements;

//   // create a global numbering system
//   hlong *globalIds = (hlong *) calloc(Ngather,sizeof(hlong));
//   int   *owner     = (int *) calloc(Ngather,sizeof(int));

//   // every gathered degree of freedom has its own global id
//   MPI_Allgather(&Ngather, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
//   for(int r=0;r<mesh->size;++r)
//     globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

//   //use the offsets to set a consecutive global numbering
//   for (dlong n =0;n<elliptic->ogs->Ngather;n++) {
//     globalIds[n] = n + globalStarts[rank];
//     owner[n] = rank;
//   }

//   //scatter this numbering to the original nodes
//   hlong *globalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
//   int *globalOwners = (int *) calloc(Ntotal,sizeof(int));
//   for (dlong n=0;n<Ntotal;n++) globalNumbering[n] = -1;
//   ogsScatter(globalNumbering, globalIds, ogsHlong, ogsAdd, elliptic->ogs);
//   ogsScatter(globalOwners, owner, ogsInt, ogsAdd, elliptic->ogs);

//   free(globalIds); free(owner);

//   // Build non-zeros of stiffness matrix (unassembled)
//   dlong nnzLocal = mesh->Np*mesh->Np*mesh->Nelements;

//   nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
//   int *AsendCounts  = (int*) calloc(mesh->size, sizeof(int));
//   int *ArecvCounts  = (int*) calloc(mesh->size, sizeof(int));
//   int *AsendOffsets = (int*) calloc(mesh->size+1, sizeof(int));
//   int *ArecvOffsets = (int*) calloc(mesh->size+1, sizeof(int));

//   dfloat *Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
//   dfloat *Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
//   dfloat *Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
//   dfloat *MM  = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

//   for (int n=0;n<mesh->Np;n++) {
//     for (int m=0;m<mesh->Np;m++) {
//       Srr[m+n*mesh->Np] = mesh->Srr[m+n*mesh->Np];
//       Srs[m+n*mesh->Np] = mesh->Srs[m+n*mesh->Np] + mesh->Ssr[m+n*mesh->Np];
//       Sss[m+n*mesh->Np] = mesh->Sss[m+n*mesh->Np];
//       MM[m+n*mesh->Np] = mesh->MM[m+n*mesh->Np];
//     }
//   }

//   if(mesh->rank==0) printf("Building full FEM matrix...");fflush(stdout);

//   //Build unassembed non-zeros
//   dlong cnt =0;
//   for (dlong e=0;e<mesh->Nelements;e++) {
//     dfloat Grr = mesh->ggeo[e*mesh->Nggeo + G00ID];
//     dfloat Grs = mesh->ggeo[e*mesh->Nggeo + G01ID];
//     dfloat Gss = mesh->ggeo[e*mesh->Nggeo + G11ID];
//     dfloat J   = mesh->ggeo[e*mesh->Nggeo + GWJID];

//     for (int n=0;n<mesh->Np;n++) {
//       if (globalNumbering[e*mesh->Np + n]<0) continue; //skip masked nodes
//       for (int m=0;m<mesh->Np;m++) {
//         if (globalNumbering[e*mesh->Np + m]<0) continue; //skip masked nodes

//         dfloat val = 0.;

//         val += Grr*Srr[m+n*mesh->Np];
//         val += Grs*Srs[m+n*mesh->Np];
//         val += Gss*Sss[m+n*mesh->Np];
//         val += J*lambda*MM[m+n*mesh->Np];

//         dfloat nonZeroThreshold = 1e-7;
//         if (fabs(val)>nonZeroThreshold) {
//           // pack non-zero
//           sendNonZeros[cnt].val = val;
//           sendNonZeros[cnt].row = globalNumbering[e*mesh->Np + n];
//           sendNonZeros[cnt].col = globalNumbering[e*mesh->Np + m];
//           sendNonZeros[cnt].ownerRank = globalOwners[e*mesh->Np + n];
//           cnt++;
//         }
//       }
//     }
//   }

//   // Make the MPI_NONZERO_T data type
//   MPI_Datatype MPI_NONZERO_T;
//   MPI_Datatype dtype[4] = {MPI_HLONG, MPI_HLONG, MPI_INT, MPI_DFLOAT};
//   int blength[4] = {1, 1, 1, 1};
//   MPI_Aint addr[4], displ[4];
//   MPI_Get_address ( &(sendNonZeros[0]          ), addr+0);
//   MPI_Get_address ( &(sendNonZeros[0].col      ), addr+1);
//   MPI_Get_address ( &(sendNonZeros[0].ownerRank), addr+2);
//   MPI_Get_address ( &(sendNonZeros[0].val      ), addr+3);
//   displ[0] = 0;
//   displ[1] = addr[1] - addr[0];
//   displ[2] = addr[2] - addr[0];
//   displ[3] = addr[3] - addr[0];
//   MPI_Type_create_struct (4, blength, displ, dtype, &MPI_NONZERO_T);
//   MPI_Type_commit (&MPI_NONZERO_T);

//   // count how many non-zeros to send to each process
//   for(dlong n=0;n<cnt;++n)
//     AsendCounts[sendNonZeros[n].ownerRank]++;

//   // sort by row ordering
//   qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);

//   // find how many nodes to expect (should use sparse version)
//   MPI_Alltoall(AsendCounts, 1, MPI_INT, ArecvCounts, 1, MPI_INT, mesh->comm);

//   // find send and recv offsets for gather
//   *nnz = 0;
//   for(int r=0;r<mesh->size;++r){
//     AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
//     ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
//     *nnz += ArecvCounts[r];
//   }

//   *A = (nonZero_t*) calloc(*nnz, sizeof(nonZero_t));

//   // determine number to receive
//   MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_NONZERO_T,
//                         (*A), ArecvCounts, ArecvOffsets, MPI_NONZERO_T,
//                         mesh->comm);

//   // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
//   qsort((*A), *nnz, sizeof(nonZero_t), parallelCompareRowColumn);

//   // compress duplicates
//   cnt = 0;
//   for(dlong n=1;n<*nnz;++n){
//     if((*A)[n].row == (*A)[cnt].row &&
//        (*A)[n].col == (*A)[cnt].col){
//       (*A)[cnt].val += (*A)[n].val;
//     }
//     else{
//       ++cnt;
//       (*A)[cnt] = (*A)[n];
//     }
//   }
//   if (*nnz) cnt++;
//   *nnz = cnt;

//   if(mesh->rank==0) printf("done.\n");

//   MPI_Barrier(mesh->comm);
//   MPI_Type_free(&MPI_NONZERO_T);

//   free(sendNonZeros);
//   free(globalNumbering); free(globalOwners);

//   free(AsendCounts);
//   free(ArecvCounts);
//   free(AsendOffsets);
//   free(ArecvOffsets);

//   free(Srr);
//   free(Srs);
//   free(Sss);
//   free(MM );
// }

// void ellipticBuildContinuousQuad3D(elliptic_t *elliptic, dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts) {

//   mesh2D *mesh = elliptic->mesh;
//   setupAide options = elliptic->options;

//   int rank = mesh->rank;

//   //use the masked gs handle to define a global ordering

//   // number of degrees of freedom on this rank (after gathering)
//   hlong Ngather = elliptic->ogs->Ngather;
//   dlong Ntotal  = mesh->Np*mesh->Nelements;

//   // create a global numbering system
//   hlong *globalIds = (hlong *) calloc(Ngather,sizeof(hlong));
//   int   *owner     = (int *) calloc(Ngather,sizeof(int));

//   // every gathered degree of freedom has its own global id
//   MPI_Allgather(&Ngather, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
//   for(int r=0;r<mesh->size;++r)
//     globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

//   //use the offsets to set a consecutive global numbering
//   for (dlong n =0;n<elliptic->ogs->Ngather;n++) {
//     globalIds[n] = n + globalStarts[rank];
//     owner[n] = rank;
//   }

//   //scatter this numbering to the original nodes
//   hlong *globalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
//   int *globalOwners = (int *) calloc(Ntotal,sizeof(int));
//   for (dlong n=0;n<Ntotal;n++) globalNumbering[n] = -1;
//   ogsScatter(globalNumbering, globalIds, ogsHlong, ogsAdd, elliptic->ogs);
//   ogsScatter(globalOwners, owner, ogsInt, ogsAdd, elliptic->ogs);

//   free(globalIds); free(owner);

//   // 2. Build non-zeros of stiffness matrix (unassembled)
//   dlong nnzLocal = mesh->Np*mesh->Np*mesh->Nelements;
//   nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
//   int *AsendCounts  = (int*) calloc(mesh->size, sizeof(int));
//   int *ArecvCounts  = (int*) calloc(mesh->size, sizeof(int));
//   int *AsendOffsets = (int*) calloc(mesh->size+1, sizeof(int));
//   int *ArecvOffsets = (int*) calloc(mesh->size+1, sizeof(int));

//   int *mask = (int *) calloc(mesh->Np*mesh->Nelements,sizeof(int));
//   for (dlong n=0;n<elliptic->Nmasked;n++) mask[elliptic->maskIds[n]] = 1;

//   if(mesh->rank==0) printf("Building full FEM matrix...");fflush(stdout);

// #if 0
//   hlong NTf = mesh->Nelements*mesh->Np * mesh->Nelements*mesh->Np ;
//   dfloat *Af = (dfloat *)calloc(NTf, sizeof(dfloat));
// #endif

//   //Build unassembed non-zeros
//   dlong cnt =0;
//   for (dlong e=0;e<mesh->Nelements;e++) {
//     for (int ny=0;ny<mesh->Nq;ny++) {
//       for (int nx=0;nx<mesh->Nq;nx++) {
//         if (mask[e*mesh->Np + nx+ny*mesh->Nq]) continue; //skip masked nodes
//         for (int my=0;my<mesh->Nq;my++) {
//           for (int mx=0;mx<mesh->Nq;mx++) {
//             if (mask[e*mesh->Np + mx+my*mesh->Nq]) continue; //skip masked nodes

//             int id;
//             dfloat val = 0.;

//              if (ny==my) {
//               for (int k=0;k<mesh->Nq;k++) {
//                 id = k+ny*mesh->Nq;
//                 dfloat Grr = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G00ID*mesh->Np];

//                 val += Grr*mesh->D[nx+k*mesh->Nq]*mesh->D[mx+k*mesh->Nq];
//               }
//             }

//             id = mx+ny*mesh->Nq;
//             dfloat Grs = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
//             val += Grs*mesh->D[nx+mx*mesh->Nq]*mesh->D[my+ny*mesh->Nq];

//             id = nx+my*mesh->Nq;
//             dfloat Gsr = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
//             val += Gsr*mesh->D[mx+nx*mesh->Nq]*mesh->D[ny+my*mesh->Nq];

//             // id = mx+ny*mesh->Nq;
//             // dfloat Grt = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G02ID*mesh->Np];
//             // val += Grt*mesh->D[nx+mx*mesh->Nq];

//             // id = nx+my*mesh->Nq;
//             // dfloat Gtr = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G02ID*mesh->Np];
//             // val += Gtr*mesh->D[mx+nx*mesh->Nq];

//             if (nx==mx) {
//               for (int k=0;k<mesh->Nq;k++) {
//                 id = nx+k*mesh->Nq;
//                 dfloat Gss = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G11ID*mesh->Np];

//                 val += Gss*mesh->D[ny+k*mesh->Nq]*mesh->D[my+k*mesh->Nq];
//               }
//             }

//             // double check following two: AK
//             // id = nx+my*mesh->Nq;
//             // dfloat Gst = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G12ID*mesh->Np];
//             // val += Gst*mesh->D[ny+my*mesh->Nq];

//             // id = mx+ny*mesh->Nq;
//             // dfloat Gts = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G12ID*mesh->Np];
//             // val += Gts*mesh->D[my+ny*mesh->Nq];

//             if ((nx==mx)&&(ny==my)) {
//               id = nx + ny*mesh->Nq;

//               // dfloat Gtt = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G22ID*mesh->Np];
//               // val += Gtt;

//               dfloat JW = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + GWJID*mesh->Np];
//               val += JW*lambda;
//             }

// #if 0
//             const hlong rowid = e*mesh->Np + nx + ny*mesh->Nq;
//             const hlong colid = e*mesh->Np + mx + my*mesh->Nq;

//             Af[rowid*mesh->Nelements*mesh->Np + colid] = val;
// #endif

//             dfloat nonZeroThreshold = 1e-7;
//             if (fabs(val)>nonZeroThreshold) {
//               // pack non-zero
//               sendNonZeros[cnt].val = val;
//               sendNonZeros[cnt].row = globalNumbering[e*mesh->Np + nx+ny*mesh->Nq];
//               sendNonZeros[cnt].col = globalNumbering[e*mesh->Np + mx+my*mesh->Nq];
//               sendNonZeros[cnt].ownerRank = globalOwners[e*mesh->Np + nx+ny*mesh->Nq];
//               cnt++;
//             }
//           }
//         }
//       }
//     }
//   }

// #if 0
//  // Write matlab dat for postprocess
//   char fname[BUFSIZ];
//   sprintf(fname, "Ax.dat");
//   FILE *fp;
//   fp = fopen(fname, "w");

//   for(hlong row = 0; row<(mesh->Nelements*mesh->Np); row++){
//     for(hlong col = 0; col<(mesh->Nelements*mesh->Np); col++){
//       dfloat val = Af[row*mesh->Nelements*mesh->Np + col];
//       fprintf(fp,"%.8e ", val);
//     }
//     fprintf(fp,"\n");
//   }

//  fclose(fp);

// #endif

//   // Make the MPI_NONZERO_T data type
//   MPI_Datatype MPI_NONZERO_T;
//   MPI_Datatype dtype[4] = {MPI_HLONG, MPI_HLONG, MPI_INT, MPI_DFLOAT};
//   int blength[4] = {1, 1, 1, 1};
//   MPI_Aint addr[4], displ[4];
//   MPI_Get_address ( &(sendNonZeros[0]          ), addr+0);
//   MPI_Get_address ( &(sendNonZeros[0].col      ), addr+1);
//   MPI_Get_address ( &(sendNonZeros[0].ownerRank), addr+2);
//   MPI_Get_address ( &(sendNonZeros[0].val      ), addr+3);
//   displ[0] = 0;
//   displ[1] = addr[1] - addr[0];
//   displ[2] = addr[2] - addr[0];
//   displ[3] = addr[3] - addr[0];
//   MPI_Type_create_struct (4, blength, displ, dtype, &MPI_NONZERO_T);
//   MPI_Type_commit (&MPI_NONZERO_T);

//   // count how many non-zeros to send to each process
//   for(dlong n=0;n<cnt;++n)
//     AsendCounts[sendNonZeros[n].ownerRank]++;

//   // sort by row ordering
//   qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);

//   // find how many nodes to expect (should use sparse version)
//   MPI_Alltoall(AsendCounts, 1, MPI_INT, ArecvCounts, 1, MPI_INT, mesh->comm);

//   // find send and recv offsets for gather
//   *nnz = 0;
//   for(int r=0;r<mesh->size;++r){
//     AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
//     ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
//     *nnz += ArecvCounts[r];
//   }

//   *A = (nonZero_t*) calloc(*nnz, sizeof(nonZero_t));

//   // determine number to receive
//   MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_NONZERO_T,
//                         (*A), ArecvCounts, ArecvOffsets, MPI_NONZERO_T,
//                         mesh->comm);

//   // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
//   qsort((*A), *nnz, sizeof(nonZero_t), parallelCompareRowColumn);

//   // compress duplicates
//   cnt = 0;
//   for(dlong n=1;n<*nnz;++n){
//     if((*A)[n].row == (*A)[cnt].row &&
//        (*A)[n].col == (*A)[cnt].col){
//       (*A)[cnt].val += (*A)[n].val;
//     }
//     else{
//       ++cnt;
//       (*A)[cnt] = (*A)[n];
//     }
//   }
//   if (*nnz) cnt++;
//   *nnz = cnt;

// #if 0
//   // Write matlab dat for postprocess
//   char fname[BUFSIZ];
//   sprintf(fname, "Ax.dat");
//   FILE *fp;
//   fp = fopen(fname, "w");

//   for(dlong n=1;n<*nnz;++n){
//       fprintf(fp,"%d %d %.8e\n", (*A)[n].row+1, (*A)[n].col+1, (*A)[n].val);
//   }

//  fclose(fp);
// #endif

//   if(mesh->rank==0) printf("done.\n");

//   MPI_Barrier(mesh->comm);
//   MPI_Type_free(&MPI_NONZERO_T);

//   free(sendNonZeros);
//   free(globalNumbering); free(globalOwners);

//   free(AsendCounts);
//   free(ArecvCounts);
//   free(AsendOffsets);
//   free(ArecvOffsets);
// }

// void ellipticBuildContinuousTet3D(elliptic_t *elliptic, dfloat lambda, nonZero_t **A, dlong *nnz, ogs_t **ogs, hlong *globalStarts) {

//   mesh2D *mesh = elliptic->mesh;
//   setupAide options = elliptic->options;

//   int rank = mesh->rank;

//   //use the masked gs handle to define a global ordering

//   // number of degrees of freedom on this rank (after gathering)
//   hlong Ngather = elliptic->ogs->Ngather;
//   dlong Ntotal  = mesh->Np*mesh->Nelements;

//   // create a global numbering system
//   hlong *globalIds = (hlong *) calloc(Ngather,sizeof(hlong));
//   int   *owner     = (int *) calloc(Ngather,sizeof(int));

//   // every gathered degree of freedom has its own global id
//   MPI_Allgather(&Ngather, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
//   for(int r=0;r<mesh->size;++r)
//     globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

//   //use the offsets to set a consecutive global numbering
//   for (dlong n =0;n<elliptic->ogs->Ngather;n++) {
//     globalIds[n] = n + globalStarts[rank];
//     owner[n] = rank;
//   }

//   //scatter this numbering to the original nodes
//   hlong *globalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
//   int *globalOwners = (int *) calloc(Ntotal,sizeof(int));
//   for (dlong n=0;n<Ntotal;n++) globalNumbering[n] = -1;
//   ogsScatter(globalNumbering, globalIds, ogsHlong, ogsAdd, elliptic->ogs);
//   ogsScatter(globalOwners, owner, ogsInt, ogsAdd, elliptic->ogs);

//   free(globalIds); free(owner);

//   // Build non-zeros of stiffness matrix (unassembled)
//   dlong nnzLocal = mesh->Np*mesh->Np*mesh->Nelements;

//   nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
//   int *AsendCounts  = (int*) calloc(mesh->size, sizeof(int));
//   int *ArecvCounts  = (int*) calloc(mesh->size, sizeof(int));
//   int *AsendOffsets = (int*) calloc(mesh->size+1, sizeof(int));
//   int *ArecvOffsets = (int*) calloc(mesh->size+1, sizeof(int));

//   int *mask = (int *) calloc(mesh->Np*mesh->Nelements,sizeof(int));
//   for (dlong n=0;n<elliptic->Nmasked;n++) mask[elliptic->maskIds[n]] = 1;

//   //Build unassembed non-zeros
//   if(mesh->rank==0) printf("Building full FEM matrix...");fflush(stdout);

//   dlong cnt =0;
//   #pragma omp parallel for
//   for (dlong e=0;e<mesh->Nelements;e++) {

//     dfloat Grr = mesh->ggeo[e*mesh->Nggeo + G00ID];
//     dfloat Grs = mesh->ggeo[e*mesh->Nggeo + G01ID];
//     dfloat Grt = mesh->ggeo[e*mesh->Nggeo + G02ID];
//     dfloat Gss = mesh->ggeo[e*mesh->Nggeo + G11ID];
//     dfloat Gst = mesh->ggeo[e*mesh->Nggeo + G12ID];
//     dfloat Gtt = mesh->ggeo[e*mesh->Nggeo + G22ID];
//     dfloat J   = mesh->ggeo[e*mesh->Nggeo + GWJID];

//     for (int n=0;n<mesh->Np;n++) {
//       if (mask[e*mesh->Np + n]) continue; //skip masked nodes
//       for (int m=0;m<mesh->Np;m++) {
//         if (mask[e*mesh->Np + m]) continue; //skip masked nodes
//         dfloat val = 0.;

//         val += Grr*mesh->Srr[m+n*mesh->Np];
//         val += Grs*mesh->Srs[m+n*mesh->Np];
//         val += Grt*mesh->Srt[m+n*mesh->Np];
//         val += Grs*mesh->Ssr[m+n*mesh->Np];
//         val += Gss*mesh->Sss[m+n*mesh->Np];
//         val += Gst*mesh->Sst[m+n*mesh->Np];
//         val += Grt*mesh->Str[m+n*mesh->Np];
//         val += Gst*mesh->Sts[m+n*mesh->Np];
//         val += Gtt*mesh->Stt[m+n*mesh->Np];
//         val += J*lambda*mesh->MM[m+n*mesh->Np];

//         dfloat nonZeroThreshold = 1e-7;
//         if (fabs(val)>nonZeroThreshold) {
//           #pragma omp critical
//           {
//             // pack non-zero
//             sendNonZeros[cnt].val = val;
//             sendNonZeros[cnt].row = globalNumbering[e*mesh->Np + n];
//             sendNonZeros[cnt].col = globalNumbering[e*mesh->Np + m];
//             sendNonZeros[cnt].ownerRank = globalOwners[e*mesh->Np + n];
//             cnt++;
//           }
//         }
//       }
//     }
//   }

//   // Make the MPI_NONZERO_T data type
//   MPI_Datatype MPI_NONZERO_T;
//   MPI_Datatype dtype[4] = {MPI_HLONG, MPI_HLONG, MPI_INT, MPI_DFLOAT};
//   int blength[4] = {1, 1, 1, 1};
//   MPI_Aint addr[4], displ[4];
//   MPI_Get_address ( &(sendNonZeros[0]          ), addr+0);
//   MPI_Get_address ( &(sendNonZeros[0].col      ), addr+1);
//   MPI_Get_address ( &(sendNonZeros[0].ownerRank), addr+2);
//   MPI_Get_address ( &(sendNonZeros[0].val      ), addr+3);
//   displ[0] = 0;
//   displ[1] = addr[1] - addr[0];
//   displ[2] = addr[2] - addr[0];
//   displ[3] = addr[3] - addr[0];
//   MPI_Type_create_struct (4, blength, displ, dtype, &MPI_NONZERO_T);
//   MPI_Type_commit (&MPI_NONZERO_T);

//   // count how many non-zeros to send to each process
//   for(dlong n=0;n<cnt;++n)
//     AsendCounts[sendNonZeros[n].ownerRank] += 1;

//   // sort by row ordering
//   qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);

//   // find how many nodes to expect (should use sparse version)
//   MPI_Alltoall(AsendCounts, 1, MPI_INT, ArecvCounts, 1, MPI_INT, mesh->comm);

//   // find send and recv offsets for gather
//   *nnz = 0;
//   for(int r=0;r<mesh->size;++r){
//     AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
//     ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
//     *nnz += ArecvCounts[r];
//   }

//   *A = (nonZero_t*) calloc(*nnz, sizeof(nonZero_t));

//   // determine number to receive
//   MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_NONZERO_T,
//                         (*A), ArecvCounts, ArecvOffsets, MPI_NONZERO_T,
//                         mesh->comm);

//   // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
//   qsort((*A), *nnz, sizeof(nonZero_t), parallelCompareRowColumn);

//   // compress duplicates
//   cnt = 0;
//   for(dlong n=1;n<*nnz;++n){
//     if((*A)[n].row == (*A)[cnt].row &&
//        (*A)[n].col == (*A)[cnt].col){
//       (*A)[cnt].val += (*A)[n].val;
//     }
//     else{
//       ++cnt;
//       (*A)[cnt] = (*A)[n];
//     }
//   }
//   if (*nnz) cnt++;
//   *nnz = cnt;

//   if(mesh->rank==0) printf("done.\n");

//   MPI_Barrier(mesh->comm);
//   MPI_Type_free(&MPI_NONZERO_T);

//   free(sendNonZeros);
//   free(globalNumbering); free(globalOwners);

//   free(AsendCounts);
//   free(ArecvCounts);
//   free(AsendOffsets);
//   free(ArecvOffsets);

//   free(mask);
// }
