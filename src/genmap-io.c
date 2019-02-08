#include "genmap-impl.h"

#include <stdio.h>

int GenmapCreateHandle_interface(GenmapHandle h) {
  h->dbgLevel = 0;
  h->printStat = 0;

  return 0;
}
