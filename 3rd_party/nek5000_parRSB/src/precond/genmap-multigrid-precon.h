#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

#include <genmap-impl.h>

typedef struct csr_mat_  *csr_mat;
typedef struct mgData_  *mgData;
typedef struct mgLevel_ *mgLevel;

struct csr_mat_{
  uint rn;
  ulong row_start;

  uint *row_off;
  ulong *col;
  GenmapScalar *v,*diag;

  struct gs_data *gsh;
};

// for the coarse level
void csr_mat_setup(struct array *entries,struct comm *c,csr_mat *M);
void csr_mat_apply(GenmapScalar *y,csr_mat M,GenmapScalar *x);
void csr_mat_print(csr_mat M,struct comm *c);
int  csr_mat_free(csr_mat M);
void csr_mat_gather(csr_mat M,struct gs_data *gsh,GenmapScalar *x,
  GenmapScalar *buf,buffer *bfr);
struct gs_data *get_csr_top(csr_mat M,struct comm *c);

struct mgLevel_{
  mgData data;
  int nsmooth;
  GenmapScalar sigma;
  struct gs_data *J; // interpolation from level i to i+1
  struct gs_data *Q; // global to local conversion of a vector
  csr_mat M;
};

struct mgData_{
  struct comm c;
  genmap_handle h;
  struct gs_data *top;
  buffer bfr;
  int nlevels;
  mgLevel *levels;
  uint *level_off;
  GenmapScalar *y,*x,*b,*u,*rhs,*buf;
};

void mgSetup(GenmapComm c,csr_mat M,mgData *d);
void mgLevelSetup(mgData data,uint level);
void mgFree(mgData d);

int log2i(sint i);

typedef struct{
  ulong r,c;
  uint proc;
} csr_entry;

typedef struct{
  ulong r,c,rn,cn;
  uint p;
  GenmapScalar v;
} entry;

#define GETLNG(p,i,off) (*((ulong*)((char*)(p)+(off)+(i)*sizeof(entry))))

#define GETPTR(p,i,off)            ((char*)(p)+(off)+(i)*sizeof(entry))

void setOwner(char *ptr,sint n,size_t inOffset,size_t outOffset,
  slong lelg,sint np);

void mg_vcycle(GenmapScalar *u,GenmapScalar *rhs,mgData d);
void mg_vcycle_lvl(GenmapScalar *u1,GenmapScalar *rhs,mgData d,
  int lvl_start);

int flex_cg(genmap_handle h,GenmapComm c,mgData d,GenmapVector r,
  int maxIter,int verbose,GenmapVector x);

int project_pf(genmap_handle h,GenmapComm c,mgData d,GenmapVector r,
  int maxIter,int verbose,GenmapVector x);
int project_pf_lvl(genmap_handle h,GenmapComm c,mgData d,GenmapScalar *ri,
  int maxIter,int verbose,int lvl_start,GenmapScalar *xo);

int rqi(genmap_handle h,GenmapComm c,mgData d,GenmapVector z,
  int maxIter,int verbose,GenmapVector fiedler);
int fmg(genmap_handle h,GenmapComm c,mgData d,GenmapScalar *z,
  int maxIter,int verbose,GenmapScalar *fiedler);

#endif
