#if !defined(nek_interface_)
#define nek_interface_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dlfcn.h>
#include <mpi.h>

#include "setupAide.hpp"
#include "nekrs.hpp"

#define DECLARE_USER_FUNC(a) void nek_ ## a(void);
#define DEFINE_USER_FUNC(a) void nek_ ## a(void) { (* a ## _ptr)(); }

typedef struct {
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

  int *ifgetu, *ifgetp;

  double *cbscnrs;

  /* global vertex ids */
  long long *glo_num, ngv;

  /* Boundary data */
  char *cbc;
  int *boundaryID;

  int NboundaryIDs;

  /* id to face mapping */
  int *eface1, *eface, *icface; 

  /* dimension of the problem */
  int ndim;
  /* local problem size */
  int nelv, nelt;
  int lelt;
  /* polynomial order + 1*/
  int nx1;

  /* MPI communicator */
  MPI_Comm comm;
} nekdata_private;

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

void*  nek_ptr(const char *id);
void   nek_outfld(void);
void   nek_uic(int ifield);
void   nek_end(void);
void   nek_map_m_to_n(double *a, int na, double *b, int nb);
void   nek_outpost(double *v1, double *v2, double *v3, double *vp, double *vt, char *name);
int    nek_lglel(int e);
void   nek_uf(double *u, double *v, double *w);
int    nek_setup(MPI_Comm c, setupAide &options);
void   nek_ifoutfld(int i);
void   nek_setic(void);
void   nek_userchk(void);
int    nek_bcmap(int bid, int ifld);

#ifdef __cplusplus
}
#endif

int buildNekInterface(const char *casename, int nFields, int N, int np);
void nek_copyFrom(ins_t *ins, dfloat time, int tstep);
void nek_ocopyFrom(ins_t *ins, dfloat time, int tstep);
void nek_copyFrom(ins_t *ins, dfloat time);
void nek_copyTo(ins_t *ins, dfloat &time);
void nek_ocopyTo(ins_t *ins, dfloat &time);
void nek_copyRestart(ins_t *ins);

#endif
