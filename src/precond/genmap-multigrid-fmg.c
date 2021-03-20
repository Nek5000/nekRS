#include <math.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#if 0
int inverse(double *x, int level, int iter, genmap_handle h, genmap_comm c,
            mgData d) {
  assert(level < d->nlevels - 2);

  int lelt = d->level_off[1] - d->level_off[0];
  GenmapScalar *y;
  GenmapMalloc(lelt, &y);
  lelt = d->level_off[level + 1] - d->level_off[level];

  uint i, j;
  for (i = 0; i < 10; i++) {
    // TODO: 1-orthogonalize

    // solve Ay=x
    metric_tic(&c->gsc, PROJECT);
    j = project_lvl(h, c, d, x, iter, level, y);
    metric_toc(&c->gsc, PROJECT);
    metric_acc(NPROJECT, j);

    // TODO: 1-orthogonalize

    // x=y/||y||
    GenmapScalar norm = 0.0;
    for (j = 0; j < lelt; j++)
      norm += y[j] * y[j];
    GenmapGop(c, &norm, 1, GENMAP_SCALAR, GENMAP_SUM);
    GenmapScalar normi = 1.0 / sqrt(norm);
    for (j = 0; j < lelt; j++)
      x[j] = y[j] * normi;
  }

  GenmapFree(y);

  return 0;
}

int fmg(genmap_handle h, genmap_comm c, mgData d, GenmapScalar *z, int iter,
        GenmapScalar *y) {
  int nlevels = d->nlevels;
  mgLevel *lvls = d->levels;
  uint *lvl_off = d->level_off;

  int lvl = nlevels - 2;
  uint off = lvl_off[lvl];
  int n = lvl_off[lvl + 1] - off;
  uint size = lvl_off[nlevels] - lvl_off[0];

  GenmapScalar *r;
  GenmapCalloc(size, &r);
  uint i;
  for (i = 0; i < n; i++)
    r[off + i] = z[i];

  buffer buf;
  buffer_init(&buf, 1024);
  for (lvl = nlevels - 3; lvl >= 0; lvl--) {
    // interpolate x to lvl+1
    mgLevel l = lvls[lvl];
    off = lvl_off[lvl];
    gs(r + off, gs_double, gs_add, 0, l->J, &buf);

    // solve by inverse iteration
    inverse(r + off, lvl, iter, h, c, d);
  }

  for (i = lvl_off[0]; i < lvl_off[1]; i++)
    y[i] = r[i];

  GenmapFree(r);
  buffer_free(&buf);

  return 0;
}
#endif
