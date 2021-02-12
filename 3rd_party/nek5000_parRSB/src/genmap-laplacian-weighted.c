#include <genmap-impl.h>

int GenmapInitLaplacianWeighted(genmap_handle h, GenmapComm c){
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapInt nv = GenmapGetNVertices(h);

  GenmapRealloc(lelt, &h->weights);
  GenmapUInt numPoints = (GenmapUInt) nv * lelt;

  GenmapLong *vertices;
  GenmapMalloc(numPoints, &vertices);

  GenmapElements elements = GenmapGetElements(h);
  GenmapInt i, j;
  for(i = 0; i < lelt; i++) {
    for(j = 0; j < nv; j++)
      vertices[i * nv + j] = elements[i].vertices[j];
  }

  if(c->gsw)
    gs_free(c->gsw);

#if defined(GENMAP_DEBUG)
  double t1 = GenmapGetMaxRss();
  if(GenmapCommRank(GenmapGetLocalComm(h)) == 0)
    printf("RSS before gs_setup: %lf\n", t1);
#endif

  c->gsw=gs_setup(vertices,numPoints,&c->gsc,0,gs_crystal_router,0);

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

  gs(u, gs_double, gs_add, 0, c->gsw, &c->buf);

  for(i = 0; i < lelt; i++) {
    h->weights[i] = 0.0;
    for(j = 0; j < nv; j++)
      h->weights[i] += u[nv * i + j];
  }

  for(i = 0; i < lelt; i++)
    h->weights[i] *= -1;

  GenmapFree(u);
  GenmapFree(vertices);

  return 0;
}

int GenmapLaplacianWeighted(genmap_handle h, GenmapComm c, GenmapScalar *u, GenmapScalar *v) {
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapInt nv = GenmapGetNVertices(h);

  GenmapScalar *ucv;
  GenmapMalloc((size_t)(nv * lelt), &ucv);

  GenmapInt i, j;
  for(i = 0; i < lelt; i++)
    for(j = 0; j < nv; j++)
      ucv[nv * i + j] = u[i];

  gs(ucv, gs_double, gs_add, 0, c->gsw, &c->buf);

  for(i = 0; i < lelt; i++) {
    v[i] = h->weights[i] * u[i];
    for(j = 0; j < nv; j ++)
      v[i] += ucv[nv * i + j];
  }

  GenmapFree(ucv);

  return 0;
}
