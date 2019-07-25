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
    if(e[i].globalId0 < *min) {
      *min = e[i].globalId0;
    }
    if(e[i].globalId0 > *max) {
      *max = e[i].globalId0;
    }
  }

  GenmapGop(GenmapGetLocalComm(h), min, 1, GENMAP_SCALAR, GENMAP_MIN);
  GenmapGop(GenmapGetLocalComm(h), max, 1, GENMAP_SCALAR, GENMAP_MAX);
}

GenmapInt GenmapSetFiedlerBin(GenmapHandle h) {
  GenmapScalar min, max;
  GenmapFiedlerMinMax(h, &min, &max);
  GenmapScalar range = max - min;

  int np = GenmapCommSize(GenmapGetLocalComm(h));
  GenmapInt nbins = np;
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapElements elements = GenmapGetElements(h);

  GenmapElements p, e;
  for(p = elements, e = p + lelt; p != e; p++) {
    int id;
    for(id = 0; id < np; id++) {
      GenmapScalar start = min + (range * id) / nbins;
      GenmapScalar end = min + (range * (id + 1)) / nbins;
      if(start <= p->fiedler && p->fiedler < end) {
        p->proc = id;
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
      GenmapLong start = min + (range * id) / nbins;
      GenmapLong end = min + (range * (id + 1)) / nbins;
      if(start <= p->globalId0 && p->globalId0 < end) {
        p->proc = id;
        break;
      }
    }
    if(id == np) p->proc = np - 1;
  }

  return 0;
}

void GenmapAssignBins(GenmapHandle h, int field, buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  if(field == GENMAP_FIEDLER) { // Fiedler
    // sort locally according to Fiedler vector
    sarray_sort(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
                TYPE_DOUBLE, buf0);
    //sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
    //            TYPE_DOUBLE, globalId, TYPE_LONG, buf0);
    // set the bin based on Fiedler vector
    GenmapSetFiedlerBin(h);
  } else if(GENMAP_GLOBALID) {
    // sort locally according to globalId
    sarray_sort(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                globalId0, TYPE_LONG, buf0);
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
    GenmapScan(h, GenmapGetLocalComm(h));
    elements = GenmapGetElements(h);
    lelt = GenmapGetNLocalElements(h);
    sarray_sort(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
                TYPE_DOUBLE, buf0);
    //sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
    //            TYPE_DOUBLE, globalId0, TYPE_LONG, buf0);
  } else if(field == GENMAP_GLOBALID) {
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0,
                    &(h->cr));
    GenmapScan(h, GenmapGetLocalComm(h));
    elements = GenmapGetElements(h);
    lelt = GenmapGetNLocalElements(h);
    sarray_sort(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                globalId0, TYPE_LONG, buf0);
  }
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
