#include <genmap-impl.h>
#include <genmap-io.h>

#define GENMAP_GMSH_HEX 5
#define GENMAP_GMSH_QUAD 3
#define GENMAP_ROOT 0
//
// GenmapHandle: Create
//
int GenmapCreateHandle_gmsh(GenmapHandle h) {
  h->Read = GenmapRead_gmsh;

  return 0;
}
//
// Do File I/O in Parallel
//
int GenmapRead_gmsh(GenmapHandle h, void *data) {
  char *fileName = (char *) data;

  GenmapLong NglobalElements, NglobalNodes, NlocalElementsLong;
  GenmapInt NlocalElements;
  GenmapInt rank, size;
  rank = GenmapId(h->global);
  size = GenmapNp(h->global);

  struct comm *c = &(h->global->gsComm);
  struct crystal cr;
  crystal_init(&cr, c);

  GenmapScalar *VX, *VY, *VZ;

  char *status;
  char buf[BUFSIZ];
  FILE *fp = NULL;

  if(rank == GENMAP_ROOT) {
    fp = fopen(fileName, "r");
    if(fp == NULL) {
      printf("GenmapRead_gmsh:%s:%d can't open file.\n", __FILE__, __LINE__);
    }

    do {
      status = fgets(buf, BUFSIZ, fp);
    } while(!strstr(buf, "$Nodes"));

    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, GenmapLongFormat, &NglobalNodes);

    GenmapCalloc(NglobalNodes, &VX);
    GenmapCalloc(NglobalNodes, &VY);
    GenmapCalloc(NglobalNodes, &VZ);
    GenmapInt i;
    for(i = 0; i < NglobalNodes; i++) {
      status = fgets(buf, BUFSIZ, fp);
      if(status == NULL) {
        printf("Error reading data in %s:%d.\n", __FILE__, __LINE__);
        exit(1);
      }
      sscanf(buf, "%*d "GenmapScalarFormat GenmapScalarFormat
             GenmapScalarFormat, VX + i, VY + i, VZ + i);
    }

    do {
      status = fgets(buf, BUFSIZ, fp);
    } while(!strstr(buf, "$Elements"));

    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, GenmapLongFormat, &NglobalElements);

    array_init(struct GenmapElement_private, &h->elementArray,
               NglobalElements);
    h->elementArray.n = NglobalElements;

    GenmapElements elements = GenmapGetElements(h);
    GenmapLong count = 0;

    GenmapLong j;
    for(j = 0; j < NglobalElements; j++) {
      status = fgets(buf, BUFSIZ, fp);
      GenmapLong elementId;
      GenmapInt elementType;
      GenmapLong v1, v2, v3, v4, v5, v6, v7, v8;

      sscanf(buf, GenmapLongFormat GenmapIntFormat, &elementId, &elementType);
      if(elementType == GENMAP_GMSH_HEX) {
        sscanf(buf, "%*d %*d %*d %*d %*d "GenmapLongFormat" "GenmapLongFormat" "
               GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat
               " "GenmapLongFormat" "GenmapLongFormat, &v1, &v2, &v3, &v4, &v5, &v6,
               &v7,
               &v8);
//      printf(GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "
//             GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "
//             GenmapLongFormat" "GenmapLongFormat"\n", v1, v2, v3, v4, v5, v6, v7,
//             v8);
        elements[count].globalId = elementId;
        elements[count].vertices[0] = v1;
        elements[count].vertices[1] = v2;
        elements[count].vertices[2] = v3;
        elements[count].vertices[3] = v4;
        elements[count].vertices[4] = v5;
        elements[count].vertices[5] = v6;
        elements[count].vertices[6] = v7;
        elements[count].vertices[7] = v8;
//      printf(GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "
//             GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "
//             GenmapLongFormat" "GenmapLongFormat"\n",
//             elements[count].vertices[0],
//             elements[count].vertices[1],
//             elements[count].vertices[2],
//             elements[count].vertices[3],
//             elements[count].vertices[4],
//             elements[count].vertices[5],
//             elements[count].vertices[6],
//             elements[count].vertices[7]
//             );
        count++;
      }
    }
    NglobalElements = count;
    h->elementArray.n = count;
  }

  comm_bcast(c, &NglobalElements, sizeof(GenmapLong), GENMAP_ROOT);
  comm_bcast(c, &NglobalNodes, sizeof(GenmapLong), GENMAP_ROOT);

  NlocalElements = NglobalElements / size;
  GenmapInt nrem = NglobalElements - NlocalElements * size;
  if(rank < nrem) NlocalElements++;

  GenmapLong out[2][1], bufs[2][1], *in = NULL;
  NlocalElementsLong = (GenmapLong)NlocalElements;
  comm_scan(out, &(h->global->gsComm), genmap_gs_long, gs_add,
            &NlocalElementsLong, 1, bufs);
  h->header->start = out[0][0];
  GenmapLong end = out[0][0] + NlocalElements;

#if defined(GENMAP_DEBUG)
  printf("rank="GenmapIntFormat" out="GenmapLongFormat" end="GenmapLongFormat"\n",
         rank, out[0][0], end);
#endif

  if(rank == GENMAP_ROOT) {
    GenmapMalloc(size, &in);
  } else {
    array_init(struct GenmapElement_private, &(h->elementArray),
               NlocalElements);
    h->elementArray.n = 0;
  }

  comm_gather(c, &end, sizeof(GenmapLong) * 1, in, sizeof(GenmapLong) * 1,
              \
              GENMAP_ROOT);

  if(rank == GENMAP_ROOT) {
#if defined(GENMAP_DEBUG)
    for(int i = 0; i < size; i++) {
      printf("in[%d]="GenmapLongFormat"\n", i, in[i]);
    }
#endif
    GenmapElements elements = GenmapGetElements(h);
    GenmapInt proc = 0;
    GenmapLong j;
    for(j = 0; j < NglobalElements; j++) {
      if(j >= in[proc]) proc++;
      elements[j].proc = proc;
    }
  }

  sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc,
                  0, &cr);

  h->elementArray.n = NlocalElements;
  h->header->Nnodes = NglobalNodes;
  h->header->nel = NglobalElements;
  h->header->ndim = 3;
  h->header->nv = 8;
  h->header->lelt = NlocalElements;

  if(rank == GENMAP_ROOT) {
    fclose(fp);
    GenmapFree(in);
    GenmapFree(VX);
    GenmapFree(VY);
    GenmapFree(VZ);
  }

  crystal_free(&cr);

  return 0;
}


int GenmapWriteGmsh(GenmapScalar *VX, GenmapScalar *VY,
                    GenmapScalar *VZ,
                    GenmapInt *EToV, GenmapLong Nnodes, GenmapLong Nelements,
                    char *fileName) {

  FILE *fp;
  fp = fopen(fileName, "w");

  fprintf(fp, "$MeshForman\n2.2 0 8\n$EndMeshFormat\n");
  fprintf(fp, "$PhysicalNames\n$EndPhysicalNames\n");

  fprintf(fp, "$Nodes\n");
  fprintf(fp, GenmapLongFormat"\n", Nnodes);
  GenmapLong i;
  for(i = 0; i < Nnodes; i++) {
    fprintf(fp,
            GenmapLongFormat" "GenmapScalarFormat" "GenmapScalarFormat" "GenmapScalarFormat
            "\n",
            i, VX[i], VY[i], VZ[i]);
  }
  fprintf(fp, "$EndNodes\n");

  fprintf(fp, "$Elements\n");
  fprintf(fp, GenmapLongFormat"\n", Nelements);
  GenmapLong j1, j2, j3, hexid;
  j1 = j2 = j3 = 9; hexid = GENMAP_GMSH_HEX;
  for(i = 0; i < Nelements; i++) {
    GenmapLong v1 = EToV[8 * i + 0];
    GenmapLong v2 = EToV[8 * i + 1];
    GenmapLong v3 = EToV[8 * i + 2];
    GenmapLong v4 = EToV[8 * i + 3];
    GenmapLong v5 = EToV[8 * i + 4];
    GenmapLong v6 = EToV[8 * i + 5];
    GenmapLong v7 = EToV[8 * i + 6];
    GenmapLong v8 = EToV[8 * i + 7];
    fprintf(fp,
            GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "\
            GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "\
            GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "GenmapLongFormat" "\
            GenmapLongFormat"\n",
            i, hexid, j1, j2, j3, v1, v2, v3, v4, v5, v6, v7, v8);
  }
  fprintf(fp, "$EndElements\n");

  fclose(fp);

  return 0;
}

