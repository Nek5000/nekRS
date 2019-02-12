#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

void GenmapFiedlerMinMax(GenmapHandle h, GenmapScalar *min, GenmapScalar *max) {
  *min = 1; *max = -1;

  GenmapElements e = GenmapGetElements(h);
  GenmapInt i;
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    if(e[i].fiedler < *min) {
      *min = e[i].fiedler;
    }
    if(e[i].fiedler > *max) {
      *max = e[i].fiedler;
    }
  }

  GenmapGop(GenmapGetLocalComm(h), min, 1, GENMAP_SCALAR, GENMAP_MIN);
  GenmapGop(GenmapGetLocalComm(h), max, 1, GENMAP_SCALAR, GENMAP_MAX);
}

void GenmapGlobalIdMinMax(GenmapHandle h, GenmapLong *min, GenmapLong *max) {
  *min = LONG_MAX; *max = LONG_MIN;

  GenmapElements e = GenmapGetElements(h);
  GenmapInt i;
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    if(e[i].globalId < *min) {
      *min = e[i].globalId;
    }
    if(e[i].globalId > *max) {
      *max = e[i].globalId;
    }
  }

  GenmapGop(GenmapGetLocalComm(h), min, 1, GENMAP_SCALAR, GENMAP_MIN);
  GenmapGop(GenmapGetLocalComm(h), max, 1, GENMAP_SCALAR, GENMAP_MAX);
}

GenmapInt GenmapSetFiedlerBin(GenmapHandle h) {
  GenmapScalar min, max;
  GenmapFiedlerMinMax(h, &min, &max);
  GenmapScalar range = max - min;

  GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));
  GenmapInt nbins = np;
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapElements elements = GenmapGetElements(h);

  GenmapElements p, e;
  for(p = elements, e = p + lelt; p != e; p++) {
    GenmapInt id;
    for(id = 0; id < np; id++) {
      GenmapScalar start = min + (range * id) / nbins;
      GenmapScalar end = min + (range * (id + 1)) / nbins;
      if(start <= p->fiedler && p->fiedler < end) {
        p->proc = id;
        p->procGlobal = GenmapCommRank(GenmapGetGlobalComm(h));
        break;
      }
    }
    if(id == np) p->proc = np - 1;
  }

  return 0;
}

GenmapInt GenmapSetGlobalIdBin(GenmapHandle h) {
  GenmapLong min, max;
  GenmapGlobalIdMinMax(h, &min, &max);
  GenmapLong range = max - min;

  GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));
  GenmapInt nbins = np;
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapElements elements = GenmapGetElements(h);

  GenmapElements p, e;
  for(p = elements, e = p + lelt; p != e; p++) {
    GenmapInt id;
    for(id = 0; id < np; id++) {
      GenmapScalar start = min + (range * id) / nbins;
      GenmapScalar end = min + (range * (id + 1)) / nbins;
      if(start <= p->globalId && p->globalId < end) {
        p->proc = id;
        p->procGlobal = GenmapCommRank(GenmapGetGlobalComm(h));
        break;
      }
    }
    if(id == np) p->proc = np - 1;
  }

  return 0;
}

int GenmapFiedler(GenmapHandle h, GenmapComm c, int maxIter, int global) {
  // 1. Do lanczos in local communicator.
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapVector initVec, alphaVec, betaVec;

  GenmapCreateVector(&initVec, GenmapGetNLocalElements(h));
  GenmapElements elements = GenmapGetElements(h);

  GenmapInt i;
#if defined(GENMAP_PAUL)
  if(global > 0) {
    for(i = 0;  i < lelt; i++) {
      if(GenmapGetLocalStartIndex(h) + i + 1  < GenmapGetNGlobalElements(h) / 2)
        initVec->data[i] = GenmapGetLocalStartIndex(h) + i + 1 + 1000. *
                           GenmapGetNGlobalElements(h);
      else
        initVec->data[i] = GenmapGetLocalStartIndex(h) + i + 1;
    }
  } else {
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }
#else
  if(global > 0) {
    initVec->data[i] = (GenmapScalar) elements[i].globalId;
  } else {
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }
#endif

  GenmapCreateVector(&alphaVec, maxIter);
  GenmapCreateVector(&betaVec, maxIter - 1);
  GenmapVector *q = NULL;

#if defined(GENMAP_PAUL)
  GenmapOrthogonalizebyOneVector(h, c, initVec, GenmapGetNGlobalElements(h));
  GenmapScalar rtr = GenmapDotVector(initVec, initVec);
  GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScalar rni = 1.0 / sqrt(rtr);
  GenmapScaleVector(initVec, initVec, rni);
  int iter = GenmapLanczosLegendary(h, c, initVec, maxIter, &q, alphaVec,
                                    betaVec);
#else
  int iter = GenmapLanczos(h, c, initVec, maxIter, &q, alphaVec, betaVec);
#endif

  GenmapVector evLanczos, evTriDiag;
  GenmapCreateVector(&evTriDiag, iter);
#if defined(GENMAP_PAUL)
  // 2. Use TQLI and find the minimum eigenvalue and associated vector
  GenmapVector *eVectors, eValues;
  GenmapTQLI(h, alphaVec, betaVec, &eVectors, &eValues);

  GenmapScalar eValMin = fabs(eValues->data[0]);
  GenmapInt eValMinI = 0;
  for(GenmapInt i = 1; i < iter; i++) {
    if(fabs(eValues->data[i]) < eValMin) {
      eValMin = fabs(eValues->data[i]);
      eValMinI = i;
    }
  }

  GenmapCopyVector(evTriDiag, eVectors[eValMinI]);
#else
  // 2. Do inverse power iteration on local communicator and find
  // local Fiedler vector.
  GenmapVector evInit;
  GenmapCreateVector(&evInit, iter);

  // Setup initial vector and orthogonalize in 1-norm to (1,1,1...)
  for(i = 0; i < iter; i++) {
    evInit->data[i] = i + 1;
  }
  GenmapOrthogonalizebyOneVector(h, c, evInit, (GenmapLong)iter);

  GenmapInvPowerIter(evTriDiag, alphaVec, betaVec, evInit, 100);
#endif

  // Multiply tri-diagonal matrix by [q1, q2, ...q_{iter}]
  GenmapInt j;
  GenmapCreateZerosVector(&evLanczos, lelt);
  for(i = 0; i < lelt; i++) {
    for(j = 0; j < iter; j++) {
      evLanczos->data[i] += q[j]->data[i] * evTriDiag->data[j];
    }
  }

  GenmapScalar lNorm = 0;
  for(i = 0; i < lelt; i++) {
    lNorm += evLanczos->data[i] * evLanczos->data[i];
  }

  GenmapGop(c, &lNorm, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScaleVector(evLanczos, evLanczos, 1. / sqrt(lNorm));
  for(i = 0; i < lelt; i++) {
    elements[i].fiedler = evLanczos->data[i];
  }

  // n. Destory the data structures
  GenmapDestroyVector(initVec);
  GenmapDestroyVector(alphaVec);
  GenmapDestroyVector(betaVec);
  GenmapDestroyVector(evLanczos);
  GenmapDestroyVector(evTriDiag);
#if defined(GENMAP_PAUL)
  GenmapDestroyVector(eValues);
  for(i = 0; i < iter; i++) {
    GenmapDestroyVector(eVectors[i]);
  }
  GenmapFree(eVectors);
#else
  GenmapDestroyVector(evInit);

  for(i = 0; i < iter; i++) {
    GenmapDestroyVector(q[i]);
  }
  GenmapFree(q);
#endif

  return iter;
}

void GenmapSplitByGlobalId(GenmapHandle h) {
  GenmapLong start = GenmapGetLocalStartIndex(h);
  GenmapLong nel = GenmapGetNGlobalElements(h);
  GenmapInt id = GenmapCommRank(GenmapGetLocalComm(h));
  GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapElements elements = GenmapGetElements(h);

  GenmapInt pNel = (GenmapInt)(nel / np);
  GenmapInt nrem = (GenmapInt)(nel - pNel * np);
  GenmapInt idCount = 0;
  while(idCount * pNel + ((idCount < nrem) ? idCount : nrem) < start)
    idCount++;

  GenmapLong upLimit = idCount * pNel + ((idCount < nrem) ? idCount : nrem);
  GenmapLong downLimit = start;
  do {
    GenmapInt end = upLimit - start < lelt ? (GenmapInt)(upLimit - start) : lelt;
    GenmapInt i;
    for(i = (GenmapInt)(downLimit - start); i < end; i++)
      elements[i].proc = idCount - 1;
    downLimit = upLimit;
    idCount++;
    upLimit = idCount * pNel + ((idCount < nrem) ? idCount : nrem);
  } while(downLimit - start < lelt);
}

void GenmapSplitByMedian(GenmapHandle h) {
  GenmapLong start = GenmapGetLocalStartIndex(h);
  GenmapLong nel = GenmapGetNGlobalElements(h);
  GenmapInt id = GenmapCommRank(GenmapGetLocalComm(h));
  GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapElements elements = GenmapGetElements(h);

  GenmapInt pNel = (GenmapInt)(nel / np);
  GenmapInt nrem = (GenmapInt)(nel - pNel * np);
  GenmapInt idCount = 0;
  while(idCount * pNel + ((idCount < nrem) ? idCount : nrem) < start)
    idCount++;

  GenmapLong upLimit = idCount * pNel + ((idCount < nrem) ? idCount : nrem);
  GenmapLong downLimit = start;
  do {
    GenmapInt end = upLimit - start < lelt ? (GenmapInt)(upLimit - start) : lelt;
    GenmapInt i;
    for(i = (GenmapInt)(downLimit - start); i < end; i++)
      elements[i].proc = idCount - 1;
    downLimit = upLimit;
    idCount++;
    upLimit = idCount * pNel + ((idCount < nrem) ? idCount : nrem);
  } while(downLimit - start < lelt);
}

void GenmapAssignBins(GenmapHandle h, int field, buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  if(field == GENMAP_FIEDLER) { // Fiedler
    // sort locally according to Fiedler vector
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
                  TYPE_DOUBLE, globalId, TYPE_LONG, buf0);
    // set the bin based on Fiedler vector
    GenmapSetFiedlerBin(h);
  } else if(GENMAP_GLOBALID) {
    // sort locally according to globalId
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                  globalId, TYPE_LONG, globalId, TYPE_LONG, buf0);
    // set the bin based on globalId
    GenmapSetGlobalIdBin(h);
  }
}

void GenmapTransferToBins(GenmapHandle h, int field, buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  if(field == GENMAP_FIEDLER) { // Fiedler
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0,
                    &(h->cr));
    elements = GenmapGetElements(h);
    lelt = GenmapGetNLocalElements(h);
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
                  TYPE_DOUBLE, globalId, TYPE_LONG, buf0);
  } else if(field == GENMAP_GLOBALID) {
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0,
                    &(h->cr));
    elements = GenmapGetElements(h);
    lelt = GenmapGetNLocalElements(h);
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                  globalId, TYPE_LONG, globalId, TYPE_LONG, buf0);
  }
}

void GenmapBinSort(GenmapHandle h, int field, buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  GenmapAssignBins(h, field, buf0);
  GenmapTransferToBins(h, field, buf0);
  GenmapScan(h, GenmapGetLocalComm(h));
  if(field == GENMAP_FIEDLER) {
    GenmapSplitByMedian(h);
  } else if(field == GENMAP_GLOBALID) {
    GenmapSplitByGlobalId(h);
  }
  GenmapTransferToBins(h, field, buf0);
}

void GenmapRSB(GenmapHandle h) {
  int maxIter = 50;
  int npass = 50;

  if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0 && h->dbgLevel > 0)
    printf("running RSB "), fflush(stdout);

  crystal_init(&(h->cr), &(h->local->gsComm));
  buffer buf0 = null_buffer;

  while(GenmapCommSize(GenmapGetLocalComm(h)) > 1) {
    if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0
        && h->dbgLevel > 1) printf("."), fflush(stdout);

#if defined(GENMAP_PAUL)
    int global = 1;
#else
    int global = (GenmapCommSize(GenmapGetLocalComm(h)) == GenmapCommSize(
                    GenmapGetGlobalComm(h)));
#endif

    int ipass = 0;
    int iter;
    do {
      iter = GenmapFiedler(h, GenmapGetLocalComm(h), maxIter, global);
      ipass++;
      global = 0;
    } while(ipass < npass && iter == maxIter);

    GenmapBinSort(h, 0, &buf0);
    int bin;
    GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));
    GenmapInt id = GenmapCommRank(GenmapGetLocalComm(h));
    if(id < (np + 1) / 2) bin = 0;
    else bin = 1;

    GenmapComm c = GenmapGetLocalComm(h);
    GenmapSplitComm(h, &c, bin);
    GenmapSetLocalComm(h, c);

#if defined(GENMAP_PAUL)
    GenmapBinSort(h, 1, &buf0);
#endif
  }

  crystal_free(&(h->cr));
  buffer_free(&buf0);
}
