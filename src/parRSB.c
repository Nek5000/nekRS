#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "genmap.h"
#include "genmap-impl.h"
#include "genmap-io.h"
#include "parRSB.h"

void fparRSB_partMesh(int *part, long long *vtx, int *nel,
                      int *nve, int *options, int *comm, int *err) {
  *err = 1;

  GenmapCommExternal c;
#if defined(GENMAP_MPI)
  c = MPI_Comm_f2c(*comm);
#else
  c = 0;
#endif

  *err = parRSB_partMesh(part, vtx, *nel, *nve, options, c);
}

int parRSB_partMesh(int *part, long long *vtx, int nel, int nve,
                     int *options, MPI_Comm comm) {
  GenmapHandle h;
  GenmapInit(&h, comm, "interface");

  if(options[0] != 0) {
    h->dbgLevel = options[1];
    h->printStat = options[2];
  }

  // Assert that nel is greater then zero. Will be remove in future.
  assert(nel > 0);

  GenmapSetNLocalElements(h, (GenmapInt)nel);
  GenmapScan(h, GenmapGetGlobalComm(h));
  GenmapSetNVertices(h, nve);

  GenmapInt id = GenmapCommRank(GenmapGetGlobalComm(h));
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

  GenmapCrystalInit(h, GenmapGetGlobalComm(h));
  GenmapCrystalTransfer(h, GENMAP_ORIGIN);
  GenmapCrystalFinalize(h);

  // This should hold true
  assert(GenmapGetNLocalElements(h) == nel);

  for(i = 0; i < nel; i++) {
    part[i] = e[i].proc;
  }

  GenmapFinalize(h);

  return 0;
}
