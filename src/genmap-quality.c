#include <genmap-impl.h>

#include <math.h>
#include <stdio.h>

GenmapInt GenmapPartitionQuality(GenmapHandle h) {
  GenmapInt id = GenmapCommRank(GenmapGetGlobalComm(h));
  GenmapInt np = GenmapCommSize(GenmapGetGlobalComm(h));
  GenmapInt nel = GenmapGetNLocalElements(h);
  GenmapInt nv = h->nv;

  GenmapLong *data;
  GenmapInt numPoints = nel * nv;
  GenmapMalloc(numPoints, &data);

  GenmapElements elements = GenmapGetElements(h);
  GenmapInt i, j;
  for(i = 0; i < nel; i++) {
    for(j = 0; j < nv; j++) {
      data[i * nv + j] = elements[i].vertices[j];
    }
  }

  GenmapComm c;
  GenmapCreateComm(&c, GenmapGetGlobalComm(h)->gsComm.c);
  c->verticesHandle = gs_setup(data, numPoints, &c->gsComm, 0,
                               gs_pairwise,
                               0);

  GenmapInt neighborsCount = 0;
  for(i = 0; i < np; i++) {
    if(i != id) {
      for(j = 0; j < numPoints; j++) {
        data[j] = -1;
      }
    } else {
      for(j = 0; j < numPoints; j++) {
        data[j] = id + 1;
      }
    }

    gs(data, genmap_gs_long, gs_max, 0, c->verticesHandle, &c->buf);

    for(j = 0; j < numPoints; j++) {
      if(data[j] > 0) {
        neighborsCount++;
        break;
      }
    }
  }

  GenmapInt ncMax = neighborsCount;
  GenmapInt ncMin = neighborsCount;
  GenmapInt ncSum = neighborsCount;
  GenmapGop(c, &ncMax, 1, GENMAP_INT, GENMAP_MAX);
  GenmapGop(c, &ncMin, 1, GENMAP_INT, GENMAP_MIN);
  GenmapGop(c, &ncSum, 1, GENMAP_INT, GENMAP_SUM);

  GenmapInt nelMax = nel;
  GenmapInt nelMin = nel;
  GenmapGop(c, &nelMin, 1, GENMAP_INT, GENMAP_MIN);
  GenmapGop(c, &nelMax, 1, GENMAP_INT, GENMAP_MAX);

  if (GenmapCommRank(GenmapGetGlobalComm(h)) == 0) {
    printf(
      " Max neighbors: %d | Min neighbors: %d | Avg neighbors: %lf\n",
      ncMax, ncMin, (double)ncSum/np);
    printf(
      " Max elements: %d | Min elements: %d | Balance: %lf\n",
      nelMax, nelMin, (double)nelMax/nelMin);
    fflush(stdout);
  }


  GenmapFree(data);
  GenmapDestroyComm(c);

  return neighborsCount;
}
