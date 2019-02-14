#include <genmap-impl.h>

#include <math.h>
#include <stdio.h>

void GenmapPartitionQuality(GenmapHandle h) {
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
  
  c->verticesHandle = gs_setup(data, numPoints, &c->gsComm, 0, gs_pairwise, 0);

  GenmapInt Nmsg;
  pw_data_nmsg(c->verticesHandle, &Nmsg);

  GenmapInt ncMax = Nmsg;
  GenmapInt ncMin = Nmsg;
  GenmapInt ncSum = Nmsg;

  GenmapGop(c, &ncMax, 1, GENMAP_INT, GENMAP_MAX);
  GenmapGop(c, &ncMin, 1, GENMAP_INT, GENMAP_MIN);
  GenmapGop(c, &ncSum, 1, GENMAP_INT, GENMAP_SUM);

  GenmapInt *Ncomm;
  GenmapMalloc(Nmsg, &Ncomm);
  pw_data_size(c->verticesHandle, Ncomm);

  GenmapInt nsMax = Ncomm[0];
  GenmapInt nsMin = Ncomm[0];
  GenmapInt nsSum = Ncomm[0];
  for (i=1; i<Nmsg; ++i){
    nsMax = Ncomm[i] > Ncomm[i-1] ? Ncomm[i] : Ncomm[i-1];
    nsMin = Ncomm[i] < Ncomm[i-1] ? Ncomm[i] : Ncomm[i-1];
    nsSum += Ncomm[i];
  }
  GenmapGop(c, &nsMax, 1, GENMAP_INT, GENMAP_MAX);
  GenmapGop(c, &nsMin, 1, GENMAP_INT, GENMAP_MIN);


  GenmapInt nssMin = nsSum;
  GenmapInt nssMax = nsSum;
  GenmapInt nssSum = nsSum;
  GenmapGop(c, &nssMax, 1, GENMAP_INT, GENMAP_MAX);
  GenmapGop(c, &nssMin, 1, GENMAP_INT, GENMAP_MIN);
  GenmapGop(c, &nssSum, 1, GENMAP_INT, GENMAP_SUM);

  nsSum = nsSum/Nmsg;
  GenmapGop(c, &nsSum, 1, GENMAP_INT, GENMAP_SUM);

  GenmapInt nelMax = nel;
  GenmapInt nelMin = nel;
  GenmapGop(c, &nelMax, 1, GENMAP_INT, GENMAP_MAX);
  GenmapGop(c, &nelMin, 1, GENMAP_INT, GENMAP_MIN);

  if (GenmapCommRank(GenmapGetGlobalComm(h)) == 0) {
    printf(
      " Max neighbors: %d | Min neighbors: %d | Avg neighbors: %lf\n",
      ncMax, ncMin, (double)ncSum/np);
   printf(
      " Max nvolume: %d | Min nvolume: %d | Avg nvolume: %lf\n",
      nsMax, nsMin, (double)nsSum/np);
    printf(
      " Max volume: %d | Min volume: %d | Avg volume: %lf\n",
      nssMax, nssMin, (double)nssSum/np);
    printf(
      " Max elements: %d | Min elements: %d | Balance: %lf\n",
      nelMax, nelMin, (double)nelMax/nelMin);
    fflush(stdout);
  }

  GenmapFree(data);
  GenmapDestroyComm(c);
}
