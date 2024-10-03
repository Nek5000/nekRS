#if !defined(nek_interface_)
#define nek_interface_

#include "nekrsSys.hpp"
#include <dlfcn.h>

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

void registerNekDoublePtr(char*, double*);
void registerNekIntPtr(char*, int*);

#ifdef __cplusplus
}
#endif

void buildNekInterface(const char *casename, int nFields, int N, int np, setupAide &options);

namespace nek {

struct fldData {
  double time = 0;
  double p0th = 0;
  std::vector<occa::memory> o_x; 
  std::vector<occa::memory> o_u;
  std::vector<occa::memory> o_p; 
  std::vector<occa::memory> o_t; 
  std::vector<std::vector<occa::memory>> o_s; 
};

fldData openFld(const std::string& filename, std::vector<std::string>& _availableVariables);
void readFld(fldData& data);
void writeFld(const std::string& filename,
              const fldData& data,
              bool FP64 = false,
              const std::vector<int>& elementMask = {},
              int Nout = 0,
              bool uniform = false);

void finalize(void);
void xm1N(dfloat *x, dfloat *y, dfloat *z, int Nq, dlong Nelements);
long long int localElementIdToGlobal(int _id);
int lglel(int e);
int setup(int numberActiveFields);
void bootstrap();
void ifoutfld(int i);
void getIC(void);
void getIC(int ifield);

void restartFromFile(const std::string& fileName);
void userchk(void);
int bcmap(int bid, int ifld);

int globalElementIdToRank(long long id);
int globalElementIdToLocal(long long id);

long long set_glo_num(int npts, int isTMesh);

void bdfCoeff(dfloat *g0, dfloat *coeff, dfloat *dt, int order);
void extCoeff(dfloat *coeff, dfloat *dt, int nAB, int nBDF);
void coeffAB(dfloat *coeff, dfloat *dt, int order);
void recomputeGeometry();
void printMeshMetrics();

const std::map<std::string, void*>& ptrList();

template <typename T>
T *ptr(const std::string& name)
{
  auto entry = ptrList().find(name);
  nekrsCheck(entry == ptrList().end(), MPI_COMM_SELF, EXIT_FAILURE, "Cannot find %s\n", name.c_str());
  return static_cast<T*>(entry->second);
}

} // namespace nek

#endif
