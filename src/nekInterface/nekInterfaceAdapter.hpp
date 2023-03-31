#if !defined(nek_interface_)
#define nek_interface_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dlfcn.h>
#include <mpi.h>

#include "nrs.hpp"

#define DECLARE_USER_FUNC(a) void nek_##a(void);
#define DEFINE_USER_FUNC(a)                                                                                  \
void nek_##a(void)                                                                                           \
{                                                                                                            \
(*a##_ptr)();                                                                                                \
}

struct setupAide;
struct session_data_t;

struct nekdata_private {
  double *param;

  int *istep;
  int *ifield;

  /* x,y and z co-ordinates */
  double *xm1, *ym1, *zm1;
  double *xc, *yc, *zc;

  double *unx, *uny, *unz;

  double *time;

  /* solution */
  double *vx, *vy, *vz;
  double *pr;
  double *t;

  double *qtl;

  int *ifgetu, *ifgetp, *ifgett, *ifgetps;

  /* global vertex ids */
  long long *glo_num;

  /* Boundary data */
  char *cbc;
  int *boundaryID;
  int *boundaryIDt;

  int NboundaryID;
  int NboundaryIDt;

  /* id to face mapping */
  int *eface1, *eface, *icface;

  /* dimension of the problem */
  int ndim;
  /* local problem size */
  int nelv, nelt;
  int lelt;

  int ldimt;

  /* polynomial order + 1*/
  int nx1;

  /* MPI communicator */
  MPI_Comm comm;

  /* multigrid levels */
  int *mg_nx;
  int mg_lmax;

  /* thermodynamic pressure */
  double *p0th;
  double *dp0thdt;

  /* mesh velocities */
  double *wx, *wy, *wz;
};

extern nekdata_private nekData;

#ifdef __cplusplus
extern "C" {
#endif

DECLARE_USER_FUNC(usrdat)
DECLARE_USER_FUNC(usrdat2)
DECLARE_USER_FUNC(usrdat3)
DECLARE_USER_FUNC(uservp)
DECLARE_USER_FUNC(userf)
DECLARE_USER_FUNC(userq)
DECLARE_USER_FUNC(userbc)
DECLARE_USER_FUNC(useric)
DECLARE_USER_FUNC(usrsetvert)
DECLARE_USER_FUNC(userqtl)

#ifdef __cplusplus
}
#endif

void buildNekInterface(const char *casename, int nFields, int N, int np, setupAide &options);
namespace nek {
void *ptr(const char *id);
void *scPtr(int id);
void outSolutionFld(double time, double outputTime);
void outfld(const char *filename,
            dfloat t,
            int step,
            int coords,
            int FP64,
            void *o_u,
            void *o_p,
            void *o_s,
            int NSfields);
void uic(int ifield);
void finalize(void);
void map_m_to_n(double *a, int na, double *b, int nb);
void outpost(double *v1, double *v2, double *v3, double *vp, double *vt, char *name);
int lglel(int e);
void uf(double *u, double *v, double *w);
int setup(nrs_t *nrs);
void bootstrap();
void ifoutfld(int i);
void setic(void);
void userchk(void);
int bcmap(int bid, int ifld);

void copyToNek(dfloat time, int tstep);
void ocopyToNek(void);
void ocopyToNek(dfloat time, int tstep);
void copyToNek(dfloat time);
void copyFromNek(dfloat &time);
void ocopyFromNek(dfloat &time);
long long set_glo_num(int npts, int isTMesh);

void bdfCoeff(double *g0, double *coeff, double *dt, int order);
void extCoeff(double *coeff, double *dt, int nAB, int nBDF);
void coeffAB(double *coeff, double *dt, int order);
void recomputeGeometry();
void printMeshMetrics();
} // namespace nek
#endif
