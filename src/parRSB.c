#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "genmap.h"
#include "genmap-impl.h"
#include "parRSB.h"

void fparRSB_partMesh(int *part, long long *vtx, int *nel, int *nve,
                      int *options, int *comm, int *err) {
  *err = 1;

  GenmapCommExternal c;
#if defined(GENMAP_MPI)
  c = MPI_Comm_f2c(*comm);
#else
  c = 0;
#endif

  *err = parRSB_partMesh(part, vtx, *nel, *nve, options, c);
}

int parRSB_partMesh(int *part, long long *vtx, int nel, int nve, int *options,
                    MPI_Comm comm) {
  int bin = nel > 0;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  MPI_Comm commRSB;
  MPI_Comm_split(comm, bin, rank, &commRSB);

  if(bin > 0) {
    double time0 = comm_time();
    GenmapHandle h;
    GenmapInit(&h, commRSB);

    if(options[0] != 0) {
      h->dbgLevel = options[1];
      h->printStat = options[2];
    }

    GenmapSetNLocalElements(h, (GenmapInt)nel);
    GenmapScan(h, GenmapGetGlobalComm(h));
    GenmapSetNVertices(h, nve);

    GenmapLong nelg = GenmapGetNGlobalElements(h);
    GenmapInt id = GenmapCommRank(GenmapGetGlobalComm(h));
    GenmapInt size_ = GenmapCommSize(GenmapGetGlobalComm(h));
    if((GenmapLong)size_ > nelg) {
      if(id == 0)
        printf("Total number of elements is smaller than the number of processors.\n"
               "Run with smaller number of processors.\n");
      return 1;
    }
    GenmapElements e = GenmapGetElements(h);
    GenmapLong start = GenmapGetLocalStartIndex(h);

    GenmapInt i, j;
    for(i = 0; i < nel; i++) {
      e[i].globalId = start + i;
      e[i].origin = id;
      for(j = 0; j < nve; j++) {
        e[i].vertices[j] = vtx[i * nve + j];
      }
    }

    GenmapRSB(h);

    e = GenmapGetElements(h);
    for(j = 0; j < GenmapGetNLocalElements(h); j++) {
      e[j].proc = GenmapCommRank(GenmapGetGlobalComm(h));
    }

    GenmapCrystalInit(h, GenmapGetGlobalComm(h));
    GenmapCrystalTransfer(h, GENMAP_ORIGIN);
    GenmapCrystalFinalize(h);

    // This should hold true
    assert(GenmapGetNLocalElements(h) == nel);

    e = GenmapGetElements(h);
    buffer buf; buffer_init(&buf, 1024);
    sarray_sort(struct GenmapElement_private, e, (unsigned int)nel, globalId,
                TYPE_LONG, &buf);
    buffer_free(&buf);

    for(i = 0; i < nel; i++) {
      part[i] = e[i].proc;
    }

    if(id == 0 && h->dbgLevel > 0)
      printf("\nfinished in %lfs\n", comm_time() - time0);

    if(h->printStat > 0)
      GenmapPartitionQuality(h);

    GenmapFinalize(h);
    fflush(stdout);
  }

  MPI_Comm_free(&commRSB);

  return 0;
}
