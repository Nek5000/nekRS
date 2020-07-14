#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include "mesh.h"

#if 0
#include "parmetis.h"
#include "defs.h"
#endif

void meshRecursiveSpectralBisectionPartition(mesh_t *mesh){

  printf("meshRecursiveSpectralBisectionPartition is disabled, exiting\n");
  exit(-1);

  #if 0
  
  int rank = mesh->rank;
  int size = mesh->size;
  int Nverts = mesh->Nverts;
  
  int *allNelements = (int*) calloc(size, sizeof(int));
  
  /* local number of elements */
  int Nelements = mesh->Nelements;
  
  /* find number of elements on all processors */
  MPI_Allgather(&Nelements, 1, MPI_INT, allNelements, 1, MPI_INT, mesh->comm);

  /* element distribution -- cumulative element count on processes */
  idx_t *elmdist = (idx_t*) calloc(size+1, sizeof(idx_t)); // element starts
  
  int e,r,n;
  
  elmdist[0] = 0;
  for(r=0;r<size;++r)
    elmdist[r+1] = elmdist[r] + allNelements[r];
  
  /* list of element starts */
  idx_t *eptr = (idx_t*) calloc(Nelements+1, sizeof(idx_t)); // element starts
  
  eptr[0] = 0;
  for(e=0;e<Nelements;++e)
    eptr[e+1] = eptr[e] + Nverts;
  
  /* local element to vertex */
  idx_t *eind = (idx_t*) calloc(Nverts*Nelements, sizeof(idx_t)); // element starts
  
  for(e=0;e<Nelements;++e)
    for(n=0;n<Nverts;++n)
      eind[e*Nverts+n] = mesh->EToV[e*Nverts+n];
  
  /* weight per element */
  idx_t *elmwgt = (idx_t*) calloc(Nelements, sizeof(idx_t)); // element starts

 for(e=0;e<Nelements;++e)
   elmwgt[e] = 1.;

 /* weight flag */
 int wgtflag = 0;

 /* number flag (1=fortran, 0=c) */
 int numflag = 0;
 
 /* ncon = 1 */
 int ncon = 1;
 
 /* nodes on element face */
 int ncommonnodes = mesh->NfaceVertices;

 /* number of partitions */
 int nparts = size;
 
 /* tpwgts */
 float *tpwgts = (float*) calloc(nparts, sizeof(float));
 
 for(e=0;e<nparts;++e)
   tpwgts[e] = 1./(float)nparts;
 
 #define MAXNCON 32
 float ubvec[MAXNCON];
 
 for (n=0; n<ncon; ++n)
   ubvec[n] = UNBALANCE_FRACTION;
 
 int options[10];
 
 options[0] = 1;
 options[PMV3_OPTION_DBGLVL] = 7;
 
 options[PMV3_OPTION_SEED] = 0;
 
 int edgecut;
 
 idx_t *part = (idx_t*) calloc(Nelements, sizeof(idx_t)); // element starts
 
 MPI_Comm comm;
 MPI_Comm_dup(MPI_COMM_WORLD, &comm);

 ParMETIS_V3_PartMeshKway
   (elmdist,
    eptr,
    eind,
    elmwgt,
    &wgtflag,
    &numflag,
    &ncon,
    &ncommonnodes,
    &nparts,
    tpwgts,
    ubvec,
    options,
    &edgecut,
    part,
    &comm);

 /* now repartition EToV */

 /* add up how many ints need to be sent to each process (element count in each partition) */
 int *outNdata = (int*) calloc(size, sizeof(int));
 for(e=0;e<Nelements;++e)
   outNdata[part[e]] += Nverts;
 
 /* get count of incoming elements from each process */
 int *inNdata = (int*) calloc(size, sizeof(int));
 MPI_Alltoall(outNdata, 1, MPI_INT,
	      inNdata,  1, MPI_INT,
	      MPI_COMM_WORLD);

 /* find offsets into outgoing array for each rank data */
 int *outStarts = (int*) calloc(size, sizeof(int));
 for(r=1;r<size;++r)
   outStarts[r] = outStarts[r-1]+outNdata[r-1];

 /* find offsets into incoming array for each rank's data */
 int *inStarts = (int*) calloc(size, sizeof(int));
 for(r=1;r<size;++r)
   inStarts[r] = inStarts[r-1]+inNdata[r-1];
 
 /* create array for outgoing data */
 int *outEToV = (int*) calloc(Nelements*Nverts, sizeof(int));
 int *outElementInfo = (int*) calloc(Nelements*Nverts, sizeof(int));
 
 int *outCnt  = (int*) calloc(size, sizeof(int));
 dfloat *outEX = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
 dfloat *outEY = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
 dfloat *outEZ = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
 
 for(r=0;r<size;++r)
   outCnt[r] = outStarts[r];
 
 for(e=0;e<Nelements;++e){
   for(n=0;n<Nverts;++n){
     outEToV[outCnt[part[e]]] = mesh->EToV[e*Nverts+n];
     outEX[outCnt[part[e]]] = mesh->EX[e*Nverts+n];
     outEY[outCnt[part[e]]] = mesh->EY[e*Nverts+n];
     outEZ[outCnt[part[e]]] = mesh->EZ[e*Nverts+n];
     outElementInfo[outCnt[part[e]]] = mesh->elementInfo[e]; // yes this is lazy
     ++outCnt[part[e]];
   }
 }
 
 // reset number of elements
 Nelements = 0;
 for(r=0;r<size;++r){
   //   printf("rank %d gets %d new elements from rank %d \n", rank, inNdata[r]/Nverts, r);
   Nelements += inNdata[r]/Nverts;
 }

 /* send elements to their new rank */
 hlong *inEToV = (hlong*) calloc(Nelements*Nverts, sizeof(hlong));
 int *inElementInfo = (int*) calloc(Nelements*Nverts, sizeof(int));
 dfloat *inEX = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
 dfloat *inEY = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
 dfloat *inEZ = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
 
 MPI_Alltoallv(outEToV, outNdata, outStarts, MPI_HLONG,
	       inEToV,   inNdata,  inStarts, MPI_HLONG, 
	       mesh->comm);

 MPI_Alltoallv(outElementInfo, outNdata, outStarts, MPI_INT,
	       inElementInfo,   inNdata,  inStarts, MPI_INT,
	       mesh->comm);

 MPI_Alltoallv(outEX, outNdata, outStarts, MPI_DFLOAT,
	       inEX,   inNdata,  inStarts, MPI_DFLOAT, 
	       mesh->comm);

 MPI_Alltoallv(outEY, outNdata, outStarts, MPI_DFLOAT,
	       inEY,   inNdata,  inStarts, MPI_DFLOAT, 
	       mesh->comm);


 MPI_Alltoallv(outEZ, outNdata, outStarts, MPI_DFLOAT,
	       inEZ,   inNdata,  inStarts, MPI_DFLOAT, 
	       mesh->comm);

 free(mesh->EToV);
 free(mesh->EX);
 free(mesh->EY);
 free(mesh->EZ);
 free(mesh->elementInfo);
 
 // scrape EToV from inEToV (may be different type hlong to EToV)
 mesh->EToV = (hlong*) calloc(Nelements*Nverts, sizeof(hlong));
 mesh->elementInfo = (hlong*) calloc(Nelements, sizeof(hlong));
 mesh->EX   = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
 mesh->EY   = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
 mesh->EZ   = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
 mesh->elementInfo = (hlong*) calloc(Nelements, sizeof(hlong));
 for(e=0;e<Nelements;++e)
   for(n=0;n<Nverts;++n){
     mesh->EToV[e*Nverts+n] = inEToV[e*Nverts+n];
     mesh->elementInfo[e] = inElementInfo[e*Nverts]; // lazy
     mesh->EX[e*Nverts+n] = inEX[e*Nverts+n];
     mesh->EY[e*Nverts+n] = inEY[e*Nverts+n];
     mesh->EZ[e*Nverts+n] = inEZ[e*Nverts+n];

   }
 
 // reset element count
 mesh->Nelements = Nelements;
 
 // free temporaries
 free(allNelements);
 free(tpwgts);
 free(outNdata);
 free(outStarts);
 free(outEToV);
 free(outCnt);
 free(inNdata);
 free(inStarts);
 free(inEToV);

 free(elmdist);
 free(eptr);
 free(eind);
 free(elmwgt);
 free(part);
#endif
 
}
