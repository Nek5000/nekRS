#if !defined(nek_interface_)
#define nek_interface_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dlfcn.h>
#include <mpi.h>

#include "nrs.hpp"

#define DECLARE_USER_FUNC(a) void nek_ ## a(void);
#define DEFINE_USER_FUNC(a) void nek_ ## a(void) { (*a ## _ptr)(); }

typedef struct
{
  double* param;

  int* istep;
  int* ifield;

  /* x,y and z co-ordinates */
  double* xm1, * ym1, * zm1;
  double* xc, * yc, * zc;

  double* unx, * uny, * unz;

  double* time;

  /* solution */
  double* vx, * vy, * vz;
  double* pr;
  double* t;

  double* qtl;

  int* ifgetu, * ifgetp, * ifgett, * ifgetps;

  double* cbscnrs;

  /* global vertex ids */
  long long* glo_num;

  /* Boundary data */
  char* cbc;
  int* boundaryID;
  int* boundaryIDt;

  int NboundaryID;
  int NboundaryIDt;

  /* id to face mapping */
  int* eface1, * eface, * icface;

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
  int* mg_nx;
  int mg_lmax;
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

#ifdef __cplusplus
}
#endif

void*  nek_ptr(const char* id);
void*  nek_scPtr(int id);
void   nek_outSolutionFld(double time, double outputTime);
void   nek_outfld(const char* suffix, dfloat t, int coords, int FP64,
                  void* o_u, void* o_p, void* o_s,
                  int NSfields);
void   nek_uic(int ifield);
void   nek_end(void);
void   nek_map_m_to_n(double* a, int na, double* b, int nb);
void   nek_outpost(double* v1, double* v2, double* v3, double* vp, double* vt, char* name);
int    nek_lglel(int e);
void   nek_uf(double* u, double* v, double* w);
int    nek_setup(MPI_Comm c, setupAide &options, nrs_t* nrs);
void   nek_ifoutfld(int i);
void   nek_setic(void);
void   nek_userchk(void);
int    nek_bcmap(int bid, int ifld);

int buildNekInterface(const char* casename, int nFields, int N, int np);
void nek_copyFrom(dfloat time, int tstep);
void nek_ocopyFrom(void);
void nek_ocopyFrom(dfloat time, int tstep);
void nek_copyFrom(dfloat time);
void nek_copyTo(dfloat &time);
void nek_ocopyTo(dfloat &time);
long long nek_set_glo_num(int npts, int isTMesh);

void nek_bdfCoeff(double *g0, double *coeff, double *dt, int order);
void nek_extCoeff(double *coeff, double *dt, int order);
#endif
