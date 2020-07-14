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
#include "mesh3D.h"

typedef struct {
  int id;
  int level;
  dfloat weight;

  // 8 for maximum number of vertices per element in 3D
  int v[8];
  dfloat EX[8], EY[8], EZ[8];

  int cRank;
  int cId;
  int type;
} cElement_t;

typedef struct {
  int Nelements;
  int offSet;
} cluster_t;

void meshBuildMRABClusters3D(mesh3D *mesh, int lev, dfloat *weights, int *levels,
            int *Nclusters, cluster_t **clusters, int *Nelements, cElement_t **newElements);

// geometric partition of clusters of elements in 3D mesh using Morton ordering + parallelSort
dfloat meshClusteredGeometricPartition3D(mesh3D *mesh, int Nclusters, cluster_t *clusters, 
                              int *Nelements, cElement_t **elements);


/* ---------------------------------------------------------

This function is a bit spaghetti, but the general idea is 
we cluster low-MRAB-level elements together along with a 
halo and partition the mesh of clusters. This reduces the MPI 
costs of communicating on the low levels.

The algorithm performs the following steps
  - cluster elements of level lev or lower
  - put clusters together a single 'owning' process
  - sort the list of clusters using a space-filling curve
  - partition the SFC between the processors, exchange the
    elements along the processor boundaries to improve the
    partitioning.
  - If the resulting partition is acceptable, save it.
  - If not, return to the last acceptable partition, and rerun 
    the mesh setup. 

------------------------------------------------------------ */
void meshMRABWeightedPartition3D(mesh3D *mesh, dfloat *weights,
                                      int numLevels, int *levels) {

  const dfloat TOL = 0.8; //tolerance on what partitions are ruled 'acceptable'
                          // min_{ranks} totalWeight > TOL*max_{ranks} totalWeight => accepted

  int rank, size;
  rank = mesh->rank;
  size = mesh->size;
  
  int Nelements, Nclusters;

  cElement_t *elements, *acceptedPartition;
  cluster_t *clusters;

  if (!levels) numLevels = 1;

  //perform the first weigthed partitioning with no clustering
  meshBuildMRABClusters3D(mesh, -1, weights, levels, &Nclusters, &clusters, &Nelements, &elements);
  meshClusteredGeometricPartition3D(mesh, Nclusters, clusters, &Nelements, &elements);

  //initialize the accepted partition
  int acceptedNelements = Nelements;
  acceptedPartition = elements;

  for (int lev = 0; lev<mesh->MRABNlevels; lev++) {
    if (rank==0) printf("Clustering level %d...", lev);
    meshBuildMRABClusters3D(mesh, lev, weights, levels, &Nclusters, &clusters, &Nelements, &elements);
    if (rank==0) printf("done.\n");
    dfloat partQuality = meshClusteredGeometricPartition3D(mesh, Nclusters, clusters, &Nelements, &elements);

    if (partQuality > TOL) {
      if (rank ==0) printf("Accepting level %d clustered partition...(quality = %g)\n", lev, partQuality);
      free(acceptedPartition); //discard the old partition
      acceptedNelements = Nelements;
      acceptedPartition = elements; //good partition
    } else {
      if (rank ==0) printf("Regecting level %d clustered partition...(quality = %g)\n", lev, partQuality);
      free(elements); //discard this partition
      break;  
    }
  }

  //save this partition, and perform the mesh setup again.  
  mesh->Nelements = acceptedNelements;

  mesh->EToV = (hlong*) realloc(mesh->EToV,mesh->Nelements*mesh->Nverts*sizeof(hlong));
  mesh->EX = (dfloat*) realloc(mesh->EX,mesh->Nelements*mesh->Nverts*sizeof(dfloat));
  mesh->EY = (dfloat*) realloc(mesh->EY,mesh->Nelements*mesh->Nverts*sizeof(dfloat));
  mesh->EZ = (dfloat*) realloc(mesh->EZ,mesh->Nelements*mesh->Nverts*sizeof(dfloat));
  mesh->elementInfo = (int *) realloc(mesh->elementInfo,mesh->Nelements*sizeof(int));
  mesh->MRABlevel = (int *) realloc(mesh->MRABlevel,mesh->Nelements*sizeof(int));

  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Nverts;++n){
      mesh->EToV[e*mesh->Nverts + n] = acceptedPartition[e].v[n];
      mesh->EX  [e*mesh->Nverts + n] = acceptedPartition[e].EX[n];
      mesh->EY  [e*mesh->Nverts + n] = acceptedPartition[e].EY[n];
      mesh->EZ  [e*mesh->Nverts + n] = acceptedPartition[e].EZ[n];
    }
    mesh->elementInfo[e] = acceptedPartition[e].type;
    mesh->MRABlevel[e] = acceptedPartition[e].level;
  }

  // connect elements using parallel sort
  meshParallelConnect(mesh);

  // print out connectivity statistics
  meshPartitionStatistics(mesh);
  
  // connect elements to boundary faces
  meshConnectBoundary(mesh);


  if(mesh->NfaceVertices==2){ // Quad 3D
    meshLoadReferenceNodesQuad2D(mesh, mesh->N);
    meshPhysicalNodesQuad3D(mesh);
    meshGeometricFactorsQuad3D(mesh);
  }
  else if(mesh->NfaceVertices==4){ // tet
    // compute physical (x,y) locations of the element nodes
    meshPhysicalNodesTet3D(mesh);
    // compute geometric factors
    meshGeometricFactorsTet3D(mesh);
  }else{                         // Hex
    // compute physical (x,y) locations of the element nodes
    meshPhysicalNodesHex3D(mesh);
    // compute geometric factors
    meshGeometricFactorsHex3D(mesh);
  }
  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(mesh);
  


  // connect face nodes (find trace indices)
  meshConnectFaceNodes3D(mesh);
  
  // compute surface geofacs
  if(mesh->NfaceVertices==2)
    meshSurfaceGeometricFactorsQuad3D(mesh);
  else if(mesh->NfaceVertices==3)
    meshSurfaceGeometricFactorsTet3D(mesh);
  else
    meshSurfaceGeometricFactorsHex3D(mesh);
  
  // global nodes
  meshParallelConnectNodes(mesh);

  if (mesh->totalHaloPairs) {
    mesh->MRABlevel = (int *) realloc(mesh->MRABlevel,(mesh->Nelements+mesh->totalHaloPairs)*sizeof(int));
    int *MRABsendBuffer = (int *) calloc(mesh->totalHaloPairs,sizeof(int));
    meshHaloExchange(mesh, sizeof(int), mesh->MRABlevel, MRABsendBuffer, mesh->MRABlevel+mesh->Nelements);
    free(MRABsendBuffer);
  }
  
  free(acceptedPartition);
}

