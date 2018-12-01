#ifndef _GENMAP_READERS_H_
#define _GENMAP_READERS_H_

#include "genmap-gmsh.h"

int GenmapCreateHandle_interface(GenmapHandle h);
int GenmapRead_interface(GenmapHandle h, void *data);

#endif
