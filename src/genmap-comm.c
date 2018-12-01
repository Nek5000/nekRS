#include <genmap-impl.h>

int GenmapCreateComm(GenmapComm *c, GenmapCommExternal ce) {
  GenmapMalloc(1, c);
  comm_init(&(*c)->gsComm, ce);
  (*c)->verticesHandle = NULL;
  (*c)->laplacianWeights = NULL;

  return 0;
}

int GenmapDestroyComm(GenmapComm c) {
  if(c->verticesHandle)
    gs_free(c->verticesHandle);
  if(c->laplacianWeights)
    GenmapFree(c->laplacianWeights);
  comm_free(&c->gsComm);
  GenmapFree(c);

  return 0;
}

int GenmapNp(GenmapComm c) {
  return (int) c->gsComm.np;
}

int GenmapId(GenmapComm c) {
  return (int) c->gsComm.id;
}

int GenmapAx(GenmapHandle h, GenmapComm c, GenmapVector u,
             GenmapVector weights, GenmapVector v) {
  assert(u->size == v->size);

  GenmapInt lelt = u->size;
  GenmapInt nv = h->header->nv;

  GenmapScalar *ucv;
  GenmapMalloc((size_t)(nv * lelt), &ucv);

  for(GenmapInt i = 0; i < lelt; i++)
    for(GenmapInt j = 0; j < nv; j++)
      ucv[nv * i + j] = u->data[i];

  gs(ucv, genmap_gs_scalar, gs_add, 0, c->verticesHandle, NULL);

  for(GenmapInt i = 0; i < lelt; i++) {
    v->data[i] = weights->data[i] * u->data[i];
    for(GenmapInt j = 0; j < nv; j ++) {
      v->data[i] += ucv[nv * i + j];
    }
  }

  GenmapFree(ucv);

  return 0;
}

int GenmapAxInit(GenmapHandle h, GenmapComm c,
                 GenmapVector weights) {
  GenmapInt lelt = h->header->lelt;
  GenmapInt nv = h->header->nv;
  GenmapUInt numPoints = (GenmapUInt) nv * lelt;

  GenmapLong *vertices;
  GenmapMalloc(numPoints, &vertices);

  GenmapElements elements = GenmapGetElements(h);
  for(GenmapInt i = 0; i < lelt; i++) {
    for(int j = 0; j < nv; j++) {
      vertices[i * nv + j] = elements[i].vertices[j];
    }
  }

  if(c->verticesHandle)
    gs_free(c->verticesHandle);

#if defined(GENMAP_DEBUG)
  double t1 = GenmapGetMaxRss();
  if(h->Id(h->local) == 0) printf("RSS before gs_setup: %lf\n", t1);
#endif

  c->verticesHandle = gs_setup(vertices, numPoints, &c->gsComm, 0,
                               gs_crystal_router, 0);
#if defined(GENMAP_DEBUG)
  t1 = GenmapGetMaxRss();
  if(h->Id(h->local) == 0) printf("RSS after gs_setup: %lf\n", t1);
#endif

  GenmapScalar *u;
  GenmapMalloc(numPoints, &u);

  for(GenmapInt i = 0; i < lelt; i++)
    for(GenmapInt j = 0; j < nv; j++)
      u[nv * i + j] = 1.;

  gs(u, genmap_gs_scalar, gs_add, 0, c->verticesHandle, NULL);

  assert(weights->size == lelt);

  for(GenmapInt i = 0; i < lelt; i++) {
    weights->data[i] = 0.;
    for(GenmapInt j = 0; j < nv; j++) {
      weights->data[i] += u[nv * i + j];
    }
  }

  for(GenmapInt i = 0; i < lelt; i++) {
    weights->data[i] *= -1;
  }

  GenmapFree(u);
  GenmapFree(vertices);

  return 0;
}

int GenmapGop(GenmapComm c, void *v, GenmapInt size,
              GenmapDataType type, GenmapInt op) {
#ifdef GENMAP_MPI
  if(op == GENMAP_SUM) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_SUM,
                  c->gsComm.c);
  } else if(op == GENMAP_MAX) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_MAX,
                  c->gsComm.c);
  } else if(op == GENMAP_MIN) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_MIN,
                  c->gsComm.c);
  }
#endif

  return 0;
}
