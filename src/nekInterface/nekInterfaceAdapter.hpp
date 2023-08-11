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
void outfld(const char *filename,
            double t,
            int step,
            int coords,
            int FP64,
            const occa::memory& o_u,
            const occa::memory& o_p,
            const occa::memory& o_s,
            int NSfields);
void uic(int ifield);
void finalize(void);
void xm1N(dfloat *x, dfloat *y, dfloat *z, int Nq, dlong Nelements);
int lglel(int e);
int setup(nrs_t *nrs);
void bootstrap();
void ifoutfld(int i);
void getIC(void);
void restartFromFile(const std::string& fileName);
void userchk(void);
int bcmap(int bid, int ifld);

void copyToNek(double time, int tstep);
void ocopyToNek(void);
void ocopyToNek(double time, int tstep);
void copyToNek(double time);
void copyFromNek(double &time);
void ocopyFromNek(double &time);
long long set_glo_num(int npts, int isTMesh);

void bdfCoeff(dfloat *g0, dfloat *coeff, dfloat *dt, int order);
void extCoeff(dfloat *coeff, dfloat *dt, int nAB, int nBDF);
void coeffAB(dfloat *coeff, dfloat *dt, int order);
void recomputeGeometry();
void printMeshMetrics();
} // namespace nek
#endif
