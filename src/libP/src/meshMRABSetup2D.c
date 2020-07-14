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


#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"


dfloat meshMRABSetup2D(mesh2D *mesh, dfloat *EToDT, int maxLevels, dfloat finalTime) {

  int rank, size;
  rank = mesh->rank;
  size = mesh->size;

  //find global min and max dt
  dfloat dtmin, dtmax;
  dtmin = EToDT[0];
  dtmax = EToDT[0];
  for (dlong e=1;e<mesh->Nelements;e++) {
    dtmin = mymin(dtmin,EToDT[e]);
    dtmax = mymax(dtmax,EToDT[e]);
  }
  dfloat dtGmin, dtGmax;
  MPI_Allreduce(&dtmin, &dtGmin, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);    
  MPI_Allreduce(&dtmax, &dtGmax, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);    


  if (rank==0) {
    printf("----------- MRAB Setup ------------------------------\n");
    printf("dtmin = %g, dtmax = %g\n", dtGmin, dtGmax);
  }

  //number of levels
  mesh->MRABNlevels = mymin(floor(log2(dtGmax/dtGmin))+1,maxLevels);

  //shift dtGmin so that we have an integer number of steps
  int NtimeSteps = finalTime/(pow(2,mesh->MRABNlevels-1)*dtGmin);
  dtGmin = finalTime/(pow(2,mesh->MRABNlevels-1)*NtimeSteps);

  

  //compute the level of each element
  mesh->MRABlevel = (dlong *) calloc(mesh->Nelements+mesh->totalHaloPairs,sizeof(int));
  int *MRABsendBuffer;
  for(int lev=0; lev<mesh->MRABNlevels; lev++){             
    dfloat dtlev = dtGmin*pow(2,lev);   
    for(dlong e=0;e<mesh->Nelements;++e){
      if(EToDT[e] >=dtlev) 
        mesh->MRABlevel[e] = lev;
    }
  }

  //enforce one level difference between neighbours
  if (mesh->totalHaloPairs) 
    MRABsendBuffer = (int *) calloc(mesh->totalHaloPairs,sizeof(int));
  
  for (int lev=0; lev < mesh->MRABNlevels; lev++){
    if (mesh->totalHaloPairs) 
      meshHaloExchange(mesh, sizeof(int), mesh->MRABlevel, MRABsendBuffer, mesh->MRABlevel+mesh->Nelements);
    for (dlong e =0; e<mesh->Nelements;e++) {
      if (mesh->MRABlevel[e] > lev+1) { //find elements at least 2 levels higher than lev
        for (int f=0;f<mesh->Nfaces;f++) { //check for a level lev neighbour
          int eP = mesh->EToE[mesh->Nfaces*e+f];
          if (eP > -1) 
            if (mesh->MRABlevel[eP] == lev)
              mesh->MRABlevel[e] = lev + 1;  //if one exists, lower the level of this element to lev-1
        }
      }
    }
  }

  

  if (mesh->totalHaloPairs) free(MRABsendBuffer);

  //this could change the number of MRAB levels there are, so find the new max level
  mesh->MRABNlevels = 0;
  for (dlong e=0;e<mesh->Nelements;e++)
    mesh->MRABNlevels = (mesh->MRABlevel[e]>mesh->MRABNlevels) ? mesh->MRABlevel[e] : mesh->MRABNlevels;
  mesh->MRABNlevels++;
  int localNlevels = mesh->MRABNlevels;
  MPI_Allreduce(&localNlevels, &(mesh->MRABNlevels), 1, MPI_INT, MPI_MAX, mesh->comm);    
  // mesh->NtimeSteps = mesh->finalTime/(pow(2,mesh->MRABNlevels-1)*dtGmin);

  printf("MRABNlevels %d \n", mesh->MRABNlevels);

  //now we need to perform a weighted repartitioning of the mesh to optimize MRAB
  if (size>1) {
    //for the moment, just weigth the elements by the number or RHS evals per MRAB step
    // TODO: We should make this an input parameter later to handle other problems. 
    dfloat *weights = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));
    for (dlong e=0; e<mesh->Nelements;e++) {
      weights[e] = pow(2,mesh->MRABNlevels-mesh->MRABlevel[e]);
    }
    
    if (rank==0) printf("Repartitioning for MRAB...\n");
    meshMRABWeightedPartition2D(mesh,weights,mesh->MRABNlevels, mesh->MRABlevel);
  }

  //construct element and halo lists
  mesh->MRABelementIds = (int **) calloc(mesh->MRABNlevels,sizeof(int*));
  mesh->MRABhaloIds = (int **) calloc(mesh->MRABNlevels,sizeof(int*));
  
  mesh->MRABNelements = (int *) calloc(mesh->MRABNlevels,sizeof(int));
  mesh->MRABNhaloElements = (int *) calloc(mesh->MRABNlevels,sizeof(int));

  for (dlong e=0;e<mesh->Nelements;e++) {
    mesh->MRABNelements[mesh->MRABlevel[e]]++;
    for (int f=0;f<mesh->Nfaces;f++) { 
      int eP = mesh->EToE[mesh->Nfaces*e+f];
      if (eP > -1) {
        if (mesh->MRABlevel[eP] == mesh->MRABlevel[e]-1) {//check for a level lev-1 neighbour
          mesh->MRABNhaloElements[mesh->MRABlevel[e]]++;
          break;
        }
      }
    }
  }

  for (int lev =0;lev<mesh->MRABNlevels;lev++){
    mesh->MRABelementIds[lev] = (int *) calloc(mesh->MRABNelements[lev],sizeof(int));
    mesh->MRABhaloIds[lev] = (int *) calloc(mesh->MRABNhaloElements[lev],sizeof(int));
    int cnt  =0;
    int cnt2 =0;
    for (dlong e=0;e<mesh->Nelements;e++){
      if (mesh->MRABlevel[e] == lev) {
        mesh->MRABelementIds[lev][cnt++] = e;
      
        for (int f=0;f<mesh->Nfaces;f++) { 
          dlong eP = mesh->EToE[mesh->Nfaces*e+f];
          if (eP > -1) {
            if (mesh->MRABlevel[eP] == lev-1) {//check for a level lev-1 neighbour
              mesh->MRABhaloIds[lev][cnt2++] = e;
              break;
            }
          }
        }
      }
    }
  }


  
  //offset index
  mesh->MRABshiftIndex = (int *) calloc(mesh->MRABNlevels,sizeof(int));

  if (rank==0){
    printf("| Rank | Level | Nelements | Level/Level Boundary Elements | \n");
    printf("------------------------------------------------------------\n");
  }
  MPI_Barrier(mesh->comm);
  for (int r =0;r<size;r++) {
    if (r==rank) {
      for (int lev =0; lev<mesh->MRABNlevels; lev++) 
        printf("|  %d,    %d,      %d,        %d     \n", rank, lev, mesh->MRABNelements[lev], mesh->MRABNhaloElements[lev]);
      printf("------------------------------------------------------------\n");
    }
    MPI_Barrier(mesh->comm);
  }
  MPI_Barrier(mesh->comm);


  return dtGmin;
}



#if 0

void meshMRABSetup2D(mesh2D *mesh, dfloat *EToDT, int maxLevels) {

  int rank, size;
  MPI_Comm_rank(mesh->comm, &rank);
  MPI_Comm_size(mesh->comm, &size);

  //find global min and max dt
  dfloat dtmin, dtmax;
  dtmin = EToDT[0];
  dtmax = EToDT[0];
  for (int e=1;e<mesh->Nelements;e++) {
    dtmin = mymin(dtmin,EToDT[e]);
    dtmax = mymax(dtmax,EToDT[e]);
  }
  dfloat dtGmin, dtGmax;
  MPI_Allreduce(&dtmin, &dtGmin, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);    
  MPI_Allreduce(&dtmax, &dtGmax, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);    


  if (rank==0) {
    printf("----------- MRAB Setup ------------------------------\n");
    printf("dtmin = %g, dtmax = %g\n", dtGmin, dtGmax);
  }

  //number of levels
  mesh->MRABNlevels = mymin(floor(log2(dtGmax/dtGmin))+1,maxLevels);

  //shift dtGmin so that we have an integer number of steps
  mesh->NtimeSteps = mesh->finalTime/(pow(2,mesh->MRABNlevels-1)*dtGmin);
  dtGmin = mesh->finalTime/(pow(2,mesh->MRABNlevels-1)*mesh->NtimeSteps);

  mesh->dt = dtGmin; 

  //compute the level of each element
  mesh->MRABlevel = (int *) calloc(mesh->Nelements+mesh->totalHaloPairs,sizeof(int));
  int *MRABsendBuffer;
  for(int lev=0; lev<mesh->MRABNlevels; lev++){             
    dfloat dtlev = dtGmin*pow(2,lev);   
    for(int e=0;e<mesh->Nelements;++e){
      if(EToDT[e] >=dtlev) 
        mesh->MRABlevel[e] = lev;
    }
  }

  //enforce one level difference between neighbours
  if (mesh->totalHaloPairs) 
    MRABsendBuffer = (int *) calloc(mesh->totalHaloPairs,sizeof(int));
  
  for (int lev=0; lev < mesh->MRABNlevels; lev++){
    if (mesh->totalHaloPairs) 
      meshHaloExchange(mesh, sizeof(int), mesh->MRABlevel, MRABsendBuffer, mesh->MRABlevel+mesh->Nelements);
    for (int e =0; e<mesh->Nelements;e++) {
      if (mesh->MRABlevel[e] > lev+1) { //find elements at least 2 levels higher than lev
        for (int f=0;f<mesh->Nfaces;f++) { //check for a level lev neighbour
          int eP = mesh->EToE[mesh->Nfaces*e+f];
          if (eP > -1) 
            if (mesh->MRABlevel[eP] == lev)
              mesh->MRABlevel[e] = lev + 1;  //if one exists, lower the level of this element to lev-1
        }
      }
    }
  }

  

  if (mesh->totalHaloPairs) free(MRABsendBuffer);

  //this could change the number of MRAB levels there are, so find the new max level
  mesh->MRABNlevels = 0;
  for (int e=0;e<mesh->Nelements;e++)
    mesh->MRABNlevels = (mesh->MRABlevel[e]>mesh->MRABNlevels) ? mesh->MRABlevel[e] : mesh->MRABNlevels;
  mesh->MRABNlevels++;
  int localNlevels = mesh->MRABNlevels;
  MPI_Allreduce(&localNlevels, &(mesh->MRABNlevels), 1, MPI_INT, MPI_MAX, mesh->comm);    
  mesh->NtimeSteps = mesh->finalTime/(pow(2,mesh->MRABNlevels-1)*dtGmin);

  printf("MRABNlevels %d \n", mesh->MRABNlevels);

  //now we need to perform a weighted repartitioning of the mesh to optimize MRAB
  if (size>1) {
    //for the moment, just weigth the elements by the number or RHS evals per MRAB step
    // TODO: We should make this an input parameter later to handle other problems. 
    dfloat *weights = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));
    for (int e=0; e<mesh->Nelements;e++) {
      weights[e] = pow(2,mesh->MRABNlevels-mesh->MRABlevel[e]);
    }
    
    if (rank==0) printf("Repartitioning for MRAB...\n");
    meshMRABWeightedPartitionTri2D(mesh,weights,mesh->MRABNlevels, mesh->MRABlevel);
  }

  //construct element and halo lists
  mesh->MRABelementIds = (int **) calloc(mesh->MRABNlevels,sizeof(int*));
  mesh->MRABhaloIds = (int **) calloc(mesh->MRABNlevels,sizeof(int*));
  
  mesh->MRABNelements = (int *) calloc(mesh->MRABNlevels,sizeof(int));
  mesh->MRABNhaloElements = (int *) calloc(mesh->MRABNlevels,sizeof(int));

  for (int e=0;e<mesh->Nelements;e++) {
    mesh->MRABNelements[mesh->MRABlevel[e]]++;
    for (int f=0;f<mesh->Nfaces;f++) { 
      int eP = mesh->EToE[mesh->Nfaces*e+f];
      if (eP > -1) {
        if (mesh->MRABlevel[eP] == mesh->MRABlevel[e]-1) {//check for a level lev-1 neighbour
          mesh->MRABNhaloElements[mesh->MRABlevel[e]]++;
          break;
        }
      }
    }
  }

  for (int lev =0;lev<mesh->MRABNlevels;lev++){
    mesh->MRABelementIds[lev] = (int *) calloc(mesh->MRABNelements[lev],sizeof(int));
    mesh->MRABhaloIds[lev] = (int *) calloc(mesh->MRABNhaloElements[lev],sizeof(int));
    int cnt  =0;
    int cnt2 =0;
    for (int e=0;e<mesh->Nelements;e++){
      if (mesh->MRABlevel[e] == lev) {
        mesh->MRABelementIds[lev][cnt++] = e;
      
        for (int f=0;f<mesh->Nfaces;f++) { 
          int eP = mesh->EToE[mesh->Nfaces*e+f];
          if (eP > -1) {
            if (mesh->MRABlevel[eP] == lev-1) {//check for a level lev-1 neighbour
              mesh->MRABhaloIds[lev][cnt2++] = e;
              break;
            }
          }
        }
      }
    }
  }


  
  //offset index
  mesh->MRABshiftIndex = (int *) calloc(mesh->MRABNlevels,sizeof(int));

  if (rank==0){
    printf("| Rank | Level | Nelements | Level/Level Boundary Elements | \n");
    printf("------------------------------------------------------------\n");
  }
  MPI_Barrier(mesh->comm);
  for (int r =0;r<size;r++) {
    if (r==rank) {
      for (int lev =0; lev<mesh->MRABNlevels; lev++) 
        printf("|  %d,    %d,      %d,        %d     \n", rank, lev, mesh->MRABNelements[lev], mesh->MRABNhaloElements[lev]);
      printf("------------------------------------------------------------\n");
    }
    MPI_Barrier(mesh->comm);
  }
  MPI_Barrier(mesh->comm);
}
#endif
