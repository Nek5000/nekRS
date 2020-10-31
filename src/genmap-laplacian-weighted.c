#include <genmap-impl.h>

int GenmapInitLaplacianWeighted(GenmapHandle h, GenmapComm c, GenmapVector weights) {
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapInt nv = GenmapGetNVertices(h);
  GenmapUInt numPoints = (GenmapUInt) nv * lelt;

  GenmapLong *vertices;
  GenmapMalloc(numPoints, &vertices);

  GenmapElements elements = GenmapGetElements(h);
  GenmapInt i, j;
  for(i = 0; i < lelt; i++) {
    for(j = 0; j < nv; j++) {
      vertices[i * nv + j] = elements[i].vertices[j];
    }
  }

  if(c->gsh)
    gs_free(c->gsh);

#if defined(GENMAP_DEBUG)
  double t1 = GenmapGetMaxRss();
  if(GenmapCommRank(GenmapGetLocalComm(h)) == 0)
    printf("RSS before gs_setup: %lf\n", t1);
#endif

  c->gsh = gs_setup(vertices, numPoints, &c->gsc, 0,
                               gs_crystal_router, 0);
#if defined(GENMAP_DEBUG)
  t1 = GenmapGetMaxRss();
  if(GenmapCommRank(GenmapGetLocalComm(h)) == 0)
    printf("RSS after gs_setup: %lf\n", t1);
#endif

  GenmapScalar *u;
  GenmapMalloc(numPoints, &u);

  for(i = 0; i < lelt; i++)
    for(j = 0; j < nv; j++)
      u[nv * i + j] = 1.;

  gs(u, genmap_gs_scalar, gs_add, 0, c->gsh, &c->buf);

  assert(weights->size == lelt);

  for(i = 0; i < lelt; i++) {
    weights->data[i] = 0.;
    for(j = 0; j < nv; j++) {
      weights->data[i] += u[nv * i + j];
    }
  }

  for(i = 0; i < lelt; i++) {
    weights->data[i] *= -1;
  }

  GenmapFree(u);
  GenmapFree(vertices);

  return 0;
}

int GenmapLaplacianWeighted(GenmapHandle h, GenmapComm c, GenmapVector u,
                    GenmapVector weights, GenmapVector v) {
  assert(u->size == v->size);
  assert(u->size == GenmapGetNLocalElements(h));

  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapInt nv = GenmapGetNVertices(h);

  GenmapScalar *ucv;
  GenmapMalloc((size_t)(nv * lelt), &ucv);

  GenmapInt i, j;
  for(i = 0; i < lelt; i++)
    for(j = 0; j < nv; j++)
      ucv[nv * i + j] = u->data[i];

  gs(ucv, genmap_gs_scalar, gs_add, 0, c->gsh, &c->buf);

  for(i = 0; i < lelt; i++) {
    v->data[i] = weights->data[i] * u->data[i];
    for(j = 0; j < nv; j ++) {
      v->data[i] += ucv[nv * i + j];
    }
  }

  GenmapFree(ucv);

  return 0;
}
