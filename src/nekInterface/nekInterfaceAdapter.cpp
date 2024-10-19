#include "platform.hpp"
#include "nekInterfaceAdapter.hpp"
#include "bcMap.hpp"
#include "fileUtils.hpp"
#include "re2Reader.hpp"

nekdata_private nekData;
static int rank;
static setupAide *options;

static void (*usrdat_ptr)(void);
static void (*usrdat2_ptr)(void);
static void (*usrdat3_ptr)(void);
static void (*userchk_ptr)(void);
static void (*uservp_ptr)(void);
static void (*userf_ptr)(void);
static void (*userq_ptr)(void);
static void (*userbc_ptr)(void);
static void (*useric_ptr)(void);
static void (*userqtl_ptr)(void);
static void (*usrsetvert_ptr)(void);
static void (*nek_outfld_ptr)(char *,
                              double *,
                              int *,
                              int *,
                              int *,
                              double *,
                              double *,
                              double *,
                              double *,
                              double *,
                              double *,
                              double *,
                              double *,
                              double *,
                              int *,
                              int);
static void (*nek_openfld_ptr)(char *, double *, double *, int);
static void (*nek_readfld_ptr)(double *,
                               double *,
                               double *,
                               double *,
                               double *,
                               double *,
                               double *,
                               double *,
                               double *);
static void (*nek_uic_ptr)(int *);
static void (*nek_end_ptr)(void);
static void (*nek_restart_ptr)(char *, int *);
static void (*nek_map_m_to_n_ptr)(double *a, int *na, double *b, int *nb, int *if3d, double *w, int *nw);
static int (*nek_lglel_ptr)(int *);
static void (*nek_bootstrap_ptr)(int *, char *, char *, char *, int, int, int);
static void (*nek_setup_ptr)(int *,
                             int *,
                             int *,
                             int *,
                             int *,
                             int *,
                             int *,
                             int *,
                             int *,
                             double *,
                             double *,
                             double *,
                             double *,
                             double *,
                             int *);
static void (*nek_ifoutfld_ptr)(int *);
static void (*nek_setics_ptr)(void);
static int (*nek_bcmap_ptr)(int *, int *, int *);
static void (*nek_gen_bcmap_ptr)(void);
static int (*nek_nbid_ptr)(int *);
static long long (*nek_set_vert_ptr)(int *, int *);

static void (*nek_setbd_ptr)(double *, double *, int *);
static void (*nek_setabbd_ptr)(double *, double *, int *, int *);

static void (*nek_updggeom_ptr)(void);
static void (*nek_meshmetrics_ptr)(void);

static int (*nek_gllnid_ptr)(int *);
static int (*nek_gllel_ptr)(int *);

static std::map<std::string, void *> ptrListData;

void noop_func(void) {}

void check_error(const char *error)
{
  nekrsCheck(error != NULL, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", error);
}

static void map_m_to_n(double *a, int na, double *b, int nb)
{
  // output a, intput b
  int if3d = (nekData.ndim == 3);

  int maxn = na > nb ? na : nb;
  int N = 2 * maxn * maxn * maxn;
  double *w = (double *)malloc(sizeof(double) * N);

  (*nek_map_m_to_n_ptr)(a, &na, b, &nb, &if3d, w, &N);
  free(w);
}

namespace nek
{

int globalElementIdToRank(long long _id)
{
  auto id = static_cast<int>(_id + 1);
  return (*nek_gllnid_ptr)(&id);
}

int globalElementIdToLocal(long long _id)
{
  auto id = static_cast<int>(_id + 1);
  return (*nek_gllel_ptr)(&id) - 1;
}

long long int localElementIdToGlobal(int _id)
{
  int id = _id + 1;
  int gid = (*nek_lglel_ptr)(&id);
  return static_cast<long long int>(gid - 1);
}

fldData openFld(const std::string &filename, std::vector<std::string> &_availableVariables)
{
  auto fname = const_cast<char *>(filename.c_str());

  double time_;
  double p0th_;

  (*nek_openfld_ptr)(fname, &time_, &p0th_, static_cast<int>(filename.size()));

  if (*ptr<int>("getxr")) {
    _availableVariables.push_back("mesh");
  }
  if (*ptr<int>("getur")) {
    _availableVariables.push_back("velocity");
  }
  if (*ptr<int>("getpr")) {
    _availableVariables.push_back("pressure");
  }
  if (*ptr<int>("gettr")) {
    _availableVariables.push_back("temperature");
  }

  const auto nsr = *ptr<int>("npsr");
  for (int i = 0; i < nsr; i++) {
    _availableVariables.push_back("scalar" + scalarDigitStr(i));
  }

  fldData data;
  data.time = time_;
  data.p0th = p0th_;

  return data;
}

void readFld(fldData &data)
{
  const auto nxyz = nekData.nx1 * nekData.nx1 * nekData.nx1;
  const auto Nlocal = nekData.nelt * nxyz;
  const size_t nekFieldOffset = nekData.lelt * nxyz;

  occa::memory xm, ym, zm;
  if (*ptr<int>("getxr")) {
    xm = platform->memoryPool.reserve<double>(Nlocal);
    ym = platform->memoryPool.reserve<double>(Nlocal);
    zm = platform->memoryPool.reserve<double>(Nlocal);
  }

  occa::memory vx, vy, vz;
  if (*ptr<int>("getur")) {
    vx = platform->memoryPool.reserve<double>(Nlocal);
    vy = platform->memoryPool.reserve<double>(Nlocal);
    vz = platform->memoryPool.reserve<double>(Nlocal);
  }

  occa::memory pr;
  if (*ptr<int>("getpr")) {
    pr = platform->memoryPool.reserve<double>(Nlocal);
  }

  occa::memory t;
  if (*ptr<int>("gettr")) {
    t = platform->memoryPool.reserve<double>(Nlocal);
  }

  occa::memory s;
  const auto nsr = *ptr<int>("npsr");
  if (nsr) {
    s = platform->memoryPool.reserve<double>(nekFieldOffset * nsr);
  }

  (*nek_readfld_ptr)(static_cast<double *>(xm.ptr()),
                     static_cast<double *>(ym.ptr()),
                     static_cast<double *>(zm.ptr()),
                     static_cast<double *>(vx.ptr()),
                     static_cast<double *>(vy.ptr()),
                     static_cast<double *>(vz.ptr()),
                     static_cast<double *>(pr.ptr()),
                     static_cast<double *>(t.ptr()),
                     static_cast<double *>(s.ptr()));

  auto populate = [&](const std::vector<occa::memory> &fields, std::vector<occa::memory> &o_u) {
    auto o_tmpDouble = platform->device.malloc<double>(Nlocal);
    o_u.resize(fields.size());
    for (int i = 0; i < fields.size(); i++) {
      o_tmpDouble.copyFrom(fields[i]);
      o_u[i] = platform->deviceMemoryPool.reserve<dfloat>(Nlocal);
      platform->copyDoubleToDfloatKernel(Nlocal, o_tmpDouble, o_u[i]);
    }
  };

  if (xm.isInitialized()) {
    std::vector<occa::memory> fields{xm, ym, zm};
    populate(fields, data.o_x);
  }

  if (vx.isInitialized()) {
    std::vector<occa::memory> fields{vx, vy, vz};
    populate(fields, data.o_u);
  }

  if (pr.isInitialized()) {
    std::vector<occa::memory> fields{pr};
    populate(fields, data.o_p);
  }

  if (t.isInitialized()) {
    std::vector<occa::memory> fields{t};
    populate(fields, data.o_t);
  }

  for (int i = 0; i < nsr; i++) {
    std::vector<occa::memory> fields = {s.slice(i * nekFieldOffset, Nlocal)};
    std::vector<occa::memory> o_Si;
    data.o_s.push_back(o_Si);
    populate(fields, data.o_s.at(i));
  }
}

void writeFld(const std::string &filename,
              const fldData &data,
              bool FP64,
              const std::vector<int> &elementMask,
              int Nout,
              bool uniform)
{
  int step = 0;
  const auto nxyz = nekData.nx1 * nekData.nx1 * nekData.nx1;
  const auto Nlocal = nekData.nelt * nxyz;

  const auto &time = data.time;
  const auto &p0th = data.p0th;

  const auto &o_x = data.o_x;
  const auto &o_u = data.o_u;
  const auto &o_p = data.o_p;
  const auto &o_t = data.o_t;
  const auto &o_s = data.o_s;

  auto copyField = [&](occa::memory o_fldIn, double *fldOut, std::string tag) {
    nekrsCheck(o_fldIn.size() < Nlocal,
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s%s%s\n",
               "outfld: ",
               tag.c_str(),
               " is too short on T-mesh!");
    auto o_tmpDouble = platform->deviceMemoryPool.reserve<double>(Nlocal);
    if (o_fldIn.dtype() == occa::dtype::get<dfloat>()) {
      platform->copyDfloatToDoubleKernel(Nlocal, o_fldIn, o_tmpDouble);
      o_tmpDouble.copyTo(fldOut, Nlocal);
    } else {
      o_fldIn.copyTo(fldOut, Nlocal);
    }
  };

  std::vector<double> xm;
  std::vector<double> ym;
  std::vector<double> zm;
  if (o_x.size()) {
    nekrsCheck(o_x.size() < 3,
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s\n",
               "Mesh coordinates must have 3 components");

    xm.resize(Nlocal, 0);
    copyField(o_x[0], xm.data(), "o_x[0]");

    ym.resize(Nlocal, 0);
    copyField(o_x[1], ym.data(), "o_x[1]");

    zm.resize(Nlocal, 0);
    copyField(o_x[2], zm.data(), "o_x[2]");
  }

  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;
  if (o_u.size()) {
    nekrsCheck(o_u.size() < 3,
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s\n",
               "Velocity must have 3 components");

    vx.resize(Nlocal, 0);
    copyField(o_u[0], vx.data(), "o_u[0]");

    vy.resize(Nlocal, 0);
    copyField(o_u[1], vy.data(), "o_u[1]");

    vz.resize(Nlocal, 0);
    copyField(o_u[2], vz.data(), "o_u[2]");
  }

  std::vector<double> pr;
  if (o_p.size()) {
    pr.resize(Nlocal, 0);
    copyField(o_p[0], pr.data(), "o_p");
  }

  std::vector<double> temp;
  if (o_t.size()) {
    temp.resize(Nlocal, 0);
    copyField(o_t[0], temp.data(), "o_t");
  }

  int nps = 0;
  std::vector<double> ps;
  if (o_s.size()) {
    const auto nekFieldOffset = static_cast<size_t>(nekData.lelt) * nxyz;
    ps.resize(o_s.size() * nekFieldOffset, 0);
    for (int is = 0; is < o_s.size(); is++) {
      auto &o_Si = o_s[is][0];
      copyField(o_Si, ps.data() + nps * nekFieldOffset, "o_S[" + scalarDigitStr(is) + "]");
      nps++;
    }
  }

  {
    const int step_s = *(nekData.istep);
    *(nekData.istep) = step;

    const auto p0th_s = *(nekData.p0th);
    *(nekData.p0th) = p0th;

    std::vector<int> fldWriteFlag;
    fldWriteFlag.push_back(xm.size() ? 1 : 0);
    fldWriteFlag.push_back(vx.size() ? 1 : 0);
    fldWriteFlag.push_back(pr.size() ? 1 : 0);
    fldWriteFlag.push_back(temp.size() ? 1 : 0);
    for (int is = 0; is < nps; is++) {
      fldWriteFlag.push_back(1);
    }

    auto &p63 = nekData.param[62];
    const auto p63_s = p63;
    p63 = (FP64) ? 1 : 0;

    int nxo = Nout + 1;
    int ifreg = uniform;
    auto nek_out_mask = ptr<int>("out_mask");
    for (int i = 0; i < nekData.lelt; i++) {
      nek_out_mask[i] = 1;
    }

    // filter elements
    int filterEnabled = elementMask.size() ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &filterEnabled, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
    if (filterEnabled) {
      for (int i = 0; i < nekData.lelt; i++) {
        nek_out_mask[i] = 0;
      }
      for (auto &entry : elementMask) {
        nek_out_mask[entry] = 1;
      }
    }

    (*nek_outfld_ptr)(const_cast<char *>(filename.c_str()),
                      const_cast<double *>(&time),
                      fldWriteFlag.data(),
                      &nxo,
                      &ifreg,
                      xm.data(),
                      ym.data(),
                      zm.data(),
                      vx.data(),
                      vy.data(),
                      vz.data(),
                      pr.data(),
                      temp.data(),
                      ps.data(),
                      &nps,
                      filename.size());

    // filter reset
    for (int i = 0; i < nekData.lelt; i++) {
      nek_out_mask[i] = 1;
    }

    *(nekData.p0th) = p0th_s;
    *(nekData.istep) = step_s;
    p63 = p63_s;
  }
}

void getIC(int ifield)
{
  if (!useric_ptr) {
    return;
  }
  (*nek_uic_ptr)(&ifield);
}

void finalize()
{
  (*nek_end_ptr)();
}

void restartFromFile(const std::string &str_in)
{
  std::string str = str_in;
  str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());

  auto pos = str.find('+');
  if (pos == std::string::npos) {
    pos = str.length();
  }

  auto fileName = str.substr(0, pos);

  std::string options;
  if (pos != std::string::npos) {
    options = str.substr(pos);
  }
  upperCase(options);
  std::replace_copy(options.begin(), options.end(), options.begin(), '+', ' ');

  auto str_nek = fileName + " " + options;
  int len = str_nek.length();

  (*nek_restart_ptr)((char *)str_nek.c_str(), &len);
}

void getIC(void)
{
  int Nscalar;
  options->getArgs("NUMBER OF SCALARS", Nscalar);

  for (int ifield = 1; ifield <= 1 + Nscalar; ifield++) {
    getIC(ifield);
  }
}

void xm1N(dfloat *_x, dfloat *_y, dfloat *_z, int N, dlong Nelements)
{
  const int Np = (N + 1) * (N + 1) * (N + 1);
  const int nxyz = nekData.nx1 * nekData.nx1 * nekData.nx1;

  if ((N + 1) == nekData.nx1) {
    for (int i = 0; i < Nelements * Np; i++) {
      _x[i] = nekData.xm1[i];
      _y[i] = nekData.ym1[i];
      _z[i] = nekData.zm1[i];
    }
    return;
  }

  std::vector<double> x(Np);
  std::vector<double> y(Np);
  std::vector<double> z(Np);

  for (dlong e = 0; e < Nelements; ++e) {
    map_m_to_n(x.data(), N + 1, &nekData.xm1[e * nxyz], nekData.nx1);
    map_m_to_n(y.data(), N + 1, &nekData.ym1[e * nxyz], nekData.nx1);
    map_m_to_n(z.data(), N + 1, &nekData.zm1[e * nxyz], nekData.nx1);
    for (int i = 0; i < Np; i++) {
      _x[i + e * Np] = x[i];
      _y[i + e * Np] = y[i];
      _z[i + e * Np] = z[i];
    }
  }
}

void ifoutfld(int i)
{
  (*nek_ifoutfld_ptr)(&i);
}

void setics(void)
{
  (*nek_setics_ptr)();
}

void userchk(void)
{
  if (!userchk_ptr) {
    return;
  }

  if (rank == 0) {
    printf("calling nek_userchk ...\n");
  }
  (*userchk_ptr)();
}

} // namespace nek

DEFINE_USER_FUNC(usrdat)
DEFINE_USER_FUNC(usrdat2)
DEFINE_USER_FUNC(usrdat3)
DEFINE_USER_FUNC(uservp)
DEFINE_USER_FUNC(userf)
DEFINE_USER_FUNC(userq)
DEFINE_USER_FUNC(userbc)
DEFINE_USER_FUNC(useric)
DEFINE_USER_FUNC(usrsetvert)
DEFINE_USER_FUNC(userqtl)

void set_usr_handles(const char *session_in, int verbose)
{
  std::string cache_dir(getenv("NEKRS_CACHE_DIR"));
  if (platform->cacheBcast) {
    cache_dir = platform->tmpDir;
  }

  const std::string lib = cache_dir + "/nek5000/lib" + session_in + ".so";

  if (platform->comm.mpiRank == 0 && platform->verbose) {
    std::cout << "\nloading " << lib << std::endl;
  }
  void *handle = dlopen(lib.c_str(), RTLD_NOW | RTLD_LOCAL);

  nekrsCheck(!handle, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", dlerror());

  // check if we need to append an underscore
  auto us = [handle] {
    auto fptr = (void (*)(void))dlsym(handle, "usrdat_");
    if (handle) {
      return "_";
    } else {
      return "";
    }
  }();
  dlerror(); /* Clear any existing error */

  char func[100];
#define fname(s) (strcpy(func, (s)), strcat(func, us), func)

  usrdat_ptr = (void (*)(void))dlsym(handle, fname("usrdat"));
  check_error(dlerror());
  usrdat2_ptr = (void (*)(void))dlsym(handle, fname("usrdat2"));
  check_error(dlerror());
  usrdat3_ptr = (void (*)(void))dlsym(handle, fname("usrdat3"));
  check_error(dlerror());
  userchk_ptr = (void (*)(void))dlsym(handle, fname("userchk"));
  check_error(dlerror());

  nek_bootstrap_ptr =
      (void (*)(int *, char *, char *, char *, int, int, int))dlsym(handle, fname("nekf_bootstrap"));
  check_error(dlerror());
  nek_setup_ptr = (void (*)(int *,
                            int *,
                            int *,
                            int *,
                            int *,
                            int *,
                            int *,
                            int *,
                            int *,
                            double *,
                            double *,
                            double *,
                            double *,
                            double *,
                            int *))dlsym(handle, fname("nekf_setup"));
  check_error(dlerror());
  nek_uic_ptr = (void (*)(int *))dlsym(handle, fname("nekf_uic"));
  check_error(dlerror());
  nek_end_ptr = (void (*)(void))dlsym(handle, fname("nekf_end"));
  check_error(dlerror());
  nek_outfld_ptr = (void (*)(char *,
                             double *,
                             int *,
                             int *,
                             int *,
                             double *,
                             double *,
                             double *,
                             double *,
                             double *,
                             double *,
                             double *,
                             double *,
                             double *,
                             int *,
                             int))dlsym(handle, fname("nekf_outfld"));
  check_error(dlerror());
  nek_openfld_ptr = (void (*)(char *, double *, double *, int))dlsym(handle, fname("nekf_openfld"));
  check_error(dlerror());
  nek_readfld_ptr =
      (void (*)(double *, double *, double *, double *, double *, double *, double *, double *, double *))
          dlsym(handle, fname("nekf_readfld"));
  check_error(dlerror());

  nek_restart_ptr = (void (*)(char *, int *))dlsym(handle, fname("nekf_restart"));
  check_error(dlerror());
  nek_lglel_ptr = (int (*)(int *))dlsym(handle, fname("nekf_lglel"));
  check_error(dlerror());
  nek_ifoutfld_ptr = (void (*)(int *))dlsym(handle, fname("nekf_ifoutfld"));
  check_error(dlerror());
  nek_setics_ptr = (void (*)(void))dlsym(handle, fname("nekf_setics"));
  check_error(dlerror());
  nek_bcmap_ptr = (int (*)(int *, int *, int *))dlsym(handle, fname("nekf_bcmap"));
  check_error(dlerror());
  nek_gen_bcmap_ptr = (void (*)(void))dlsym(handle, fname("nekf_gen_bcmap"));
  check_error(dlerror());
  nek_map_m_to_n_ptr =
      (void (*)(double *, int *, double *, int *, int *, double *, int *))dlsym(handle, fname("map_m_to_n"));
  check_error(dlerror());
  nek_nbid_ptr = (int (*)(int *))dlsym(handle, fname("nekf_nbid"));
  check_error(dlerror());
  nek_set_vert_ptr = (long long (*)(int *, int *))dlsym(handle, fname("nekf_set_vert"));
  check_error(dlerror());

  nek_setbd_ptr = (void (*)(double *, double *, int *))dlsym(handle, fname("setbd"));
  check_error(dlerror());
  nek_setabbd_ptr = (void (*)(double *, double *, int *, int *))dlsym(handle, fname("setabbd"));
  check_error(dlerror());

  nek_updggeom_ptr = (void (*)(void))dlsym(handle, fname("nekf_updggeom"));
  check_error(dlerror());
  nek_meshmetrics_ptr = (void (*)(void))dlsym(handle, fname("mesh_metrics"));
  check_error(dlerror());

  nek_gllnid_ptr = (int (*)(int *))dlsym(handle, fname("gllnid"));
  check_error(dlerror());
  nek_gllel_ptr = (int (*)(int *))dlsym(handle, fname("gllel"));
  check_error(dlerror());

#define postfix(x) x##_ptr
#define load_or_noop(s)                                                                                      \
  do {                                                                                                       \
    postfix(s) = (void (*)(void))dlsym(handle, fname(#s));                                                   \
    if (!(postfix(s))) {                                                                                     \
      postfix(s) = noop_func;                                                                                \
    } else if (verbose && rank == 0) {                                                                       \
      printf("Loading " #s "\n");                                                                            \
    }                                                                                                        \
  } while (0)

  load_or_noop(uservp);
  load_or_noop(userf);
  load_or_noop(userq);
  load_or_noop(userbc);
  load_or_noop(useric);
  load_or_noop(userqtl);
  load_or_noop(usrsetvert);

#undef fname
#undef postfix
#undef load_or_noop
}

void mkSIZE(int lx1,
            int lxd,
            int lelt,
            int lelg,
            int ldim,
            int lpmin,
            int ldimt,
            setupAide &options,
            const char *SIZE)
{
  const int lineSize = 1024;
  char line[lineSize];
  const char *cache_dir = getenv("NEKRS_CACHE_DIR");

  const std::string installDir(getenv("NEKRS_HOME"));
  const std::string nek5000_dir = installDir + "/nek5000";

  const int verbose = options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  // Read and generate the new size file.
  snprintf(line, lineSize, "%s/core/SIZE.template", nek5000_dir.c_str());
  FILE *fp = fopen(line, "r");
  nekrsCheck(!fp, MPI_COMM_SELF, EXIT_FAILURE, "Cannot open %s!\n", line);

  char *sizeFile;
  {
    fseek(fp, 0, SEEK_END);
    const long length = ftell(fp);
    rewind(fp);
    sizeFile = (char *)calloc(4096 + length, sizeof(char));
  }

  int lx1m = lx1;
  constexpr int nsessmax = 1;
  constexpr int nMaxObj = 20;

  int count = 0;
  while (fgets(line, lineSize, fp) != NULL) {
    if (strstr(line, "parameter (lx1=") != NULL) {
      snprintf(line, lineSize, "      parameter (lx1=%d)\n", lx1);
    } else if (strstr(line, "parameter (lxd=") != NULL) {
      snprintf(line, lineSize, "      parameter (lxd=%d)\n", lxd);
    } else if (strstr(line, "parameter (lelt=") != NULL) {
      snprintf(line, lineSize, "      parameter (lelt=%d)\n", lelt);
    } else if (strstr(line, "parameter (lelg=") != NULL) {
      snprintf(line, lineSize, "      parameter (lelg=%d)\n", lelg);
    } else if (strstr(line, "parameter (ldim=") != NULL) {
      snprintf(line, lineSize, "      parameter (ldim=%d)\n", ldim);
    } else if (strstr(line, "parameter (lpmin=") != NULL) {
      snprintf(line, lineSize, "      parameter (lpmin=%d)\n", lpmin);
    } else if (strstr(line, "parameter (ldimt=") != NULL) {
      snprintf(line, lineSize, "      parameter (ldimt=%d)\n", ldimt);
    } else if (strstr(line, "parameter (mxprev=") != NULL) {
      snprintf(line, lineSize, "      parameter (mxprev=%d)\n", 1);
    } else if (strstr(line, "parameter (lgmres=") != NULL) {
      snprintf(line, lineSize, "      parameter (lgmres=%d)\n", 1);
    } else if (strstr(line, "parameter (lxo=") != NULL) {
      snprintf(line, lineSize, "      parameter (lxo=%d)\n", lx1 + 4);
    } else if (strstr(line, "parameter (lorder=") != NULL) {
      snprintf(line, lineSize, "      parameter (lorder=%d)\n", 1);
    } else if (strstr(line, "parameter (lhis=") != NULL) {
      snprintf(line, lineSize, "      parameter (lhis=%d)\n", 1);
    } else if (strstr(line, "parameter (lelr=") != NULL) {
      snprintf(line, lineSize, "      parameter (lelr=%d)\n", std::min(128 * lelt, lelg));
    } else if (strstr(line, "parameter (lx1m=") != NULL) {
      snprintf(line, lineSize, "      parameter (lx1m=%d)\n", lx1m);
    } else if (strstr(line, "parameter (nsessmax=") != NULL) {
      snprintf(line, lineSize, "      parameter (nsessmax=%d)\n", nsessmax);
    } else if (strstr(line, "parameter (maxobj=") != NULL) {
      snprintf(line, lineSize, "      parameter (maxobj=%d)\n", nMaxObj);
    }

    strcpy(sizeFile + count, line);
    count += strlen(line);
  }
  fclose(fp);

  int writeSize = 1;

  // read size if exists
  std::ifstream osize;
  osize.open(SIZE, std::ifstream::in);
  if (osize.is_open()) {
    writeSize = 0;
    std::string line;
    while (getline(osize, line)) {
      if (line.find("lelg=") != std::string::npos) {
        hlong oldval;
        sscanf(line.c_str(), "%*[^=]=%lld", &oldval);
        if (oldval < lelg) {
          writeSize = 1;
        }
      }
      if (line.find("lelt=") != std::string::npos) {
        int oldval;
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if (oldval < lelt) {
          writeSize = 1;
        }
      }
      if (line.find("lx1m=") != std::string::npos) {
        int oldval;
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if (oldval < lx1m) {
          writeSize = 1;
        }
      }
      if (line.find("lx1=") != std::string::npos) {
        int oldval;
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if (oldval != lx1) {
          writeSize = 1;
        }
      }
      if (line.find("ldimt=") != std::string::npos) {
        int oldval;
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if (oldval < ldimt) {
          writeSize = 1;
        }
      }
    }
  }
  osize.close();

  if (writeSize) {
    fp = fopen(SIZE, "w");
    fputs(sizeFile, fp);
    fclose(fp);
    free(sizeFile);
    if (verbose) {
      printf("using new SIZE\n");
    }
  } else {
    if (verbose) {
      printf("using existing SIZE file %s/SIZE\n", cache_dir);
    }
  }

  fflush(stdout);
}

void buildNekInterface(int ldimt, int N, int np, setupAide &options)
{
  int buildRank = rank;
  const bool buildNodeLocal = platform->cacheLocal;
  if (buildNodeLocal) {
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);
  }

  const int verbose = platform->verbose;

  const std::string installDir(getenv("NEKRS_HOME"));
  const std::string nek5000_dir = installDir + "/nek5000";
  const std::string cache_dir = fs::path(getenv("NEKRS_CACHE_DIR")) / "nek5000";

  const std::string casename = options.getArgs("CASENAME");
  std::string usrFile = fs::absolute(fs::path(options.getArgs("NEK USR FILE")));
  const std::string usrCaseName = casename + ".usr";
  const std::string usrFileCache = cache_dir / fs::path(usrCaseName).filename();
  const std::string libFile = cache_dir + "/lib" + casename + ".so";

  int usrFileExists;
  if (platform->comm.mpiRank == 0) {
    usrFileExists = fs::exists(usrFile) && fs::file_size(usrFile) > 0;
  }
  MPI_Bcast(&usrFileExists, 1, MPI_INT, 0, platform->comm.mpiComm);
  if (!usrFileExists) {
    usrFile = nek5000_dir + "/core/zero.usr";
  }

  int buildRequired = 0;
  if (platform->comm.mpiRank == 0) {
    if (platform->options.compareArgs("BUILD ONLY", "TRUE")) {
      buildRequired = 1;
    } else if (isFileNewer(usrFile.c_str(), libFile.c_str())) {
      buildRequired = 1;
    }
  }
  MPI_Bcast(&buildRequired, 1, MPI_INT, 0, platform->comm.mpiComm);

  int err = [&]() {
    if (buildRank == 0) {
      mkdir(cache_dir.c_str(), S_IRWXU);

      int nelgt, nelgv;
      const int ndim = 3;
      const std::string meshFile = options.getArgs("MESH FILE");
      re2::nelg(meshFile, nelgt, nelgv, MPI_COMM_SELF);

      int lelt = (nelgt / np) + 3;
      if (lelt > nelgt) {
        lelt = nelgt;
      }

      const std::string sizeFile = cache_dir + "/SIZE";
      mkSIZE(N + 1, 1, lelt, nelgt, ndim, np, ldimt, options, sizeFile.c_str());

      if (buildRequired) {
        copyFile(usrFile.c_str(), usrFileCache.c_str());
      }

      if (isFileNewer(sizeFile.c_str(), libFile.c_str())) {
        buildRequired = 1;
      }

      if (buildRequired) {
        const double tStart = MPI_Wtime();

        if (rank == 0) {
          printf("building nekInterface for lx1=%d, lelt=%d and lelg=%d ... ", N + 1, lelt, nelgt);
        }
        fflush(stdout);

        const std::string makeOutput = (rank == 0) ? cache_dir + "/make.out" : "/dev/null";
        const std::string case_dir(fs::current_path());
        const std::string include_dirs = "./ " + case_dir + " " + installDir + "/include/core/bdry";
        const std::string nekInterface_dir = installDir + "/nekInterface";
        std::string out_args = std::string(">" + makeOutput + " 2>&1");
        std::string make_args = "-j4 ";
        if (!verbose) {
          make_args += "-s ";
        } else {
          out_args = "";
        }

        const int bufSize = 4096;
        char buf[bufSize];
        snprintf(buf,
                 bufSize,
                 "cd %s"
                 " && cp -f %s/makefile.template makefile"
                 " && make %s"
                 "S=%s "
                 "OPT_INCDIR=\"%s\" "
                 "CASENAME=%s "
                 "CASEDIR=%s "
                 "-f %s/Makefile lib usr libnekInterface "
                 "%s",
                 cache_dir.c_str(),
                 nek5000_dir.c_str(),
                 make_args.c_str(),
                 nek5000_dir.c_str(),
                 include_dirs.c_str(),
                 casename.c_str(),
                 cache_dir.c_str(),
                 nekInterface_dir.c_str(),
                 out_args.c_str());

        if (verbose && rank == 0) {
          printf("\n%s\n", buf);
        }
        if (system(buf)) {
          if (rank == 0) {
            std::ifstream f(makeOutput);
            std::cerr << f.rdbuf();
          }

          return EXIT_FAILURE;
        }
        fileSync(libFile.c_str());

        if (rank == 0) {
          printf("done (%gs)\n", MPI_Wtime() - tStart);
        }
        fflush(stdout);
      } else {
        if (rank == 0) {
          printf("skip building nekInterface (SIZE requires no update)\n");
        }
        fflush(stdout);
      }
    } // buildRank

    if (platform->cacheBcast) {
      fileBcast(libFile,
                fs::path(platform->tmpDir) / fs::path("nek5000"),
                platform->comm.mpiComm,
                platform->verbose);
    }

    return 0;
  }();

  nekrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "");
}

namespace nek
{
int bcmap(int bid, int ifld, int isMesh)
{
  return (*nek_bcmap_ptr)(&bid, &ifld, &isMesh);
}

void gen_bcmap()
{
  (*nek_gen_bcmap_ptr)();
}

void bootstrap()
{
  options = &platform->options;

  int size;
  MPI_Comm_rank(platform->comm.mpiComm, &rank);
  MPI_Comm_size(platform->comm.mpiComm, &size);

  int buildRank = rank;
  const bool buildNodeLocal = platform->cacheLocal;
  if (buildNodeLocal) {
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);
  }

  int N;
  options->getArgs("POLYNOMIAL DEGREE", N);

  int Nscalar;
  options->getArgs("NUMBER OF SCALARS", Nscalar);

  int npTarget = size;
  options->getArgs("NP TARGET", npTarget);

  buildNekInterface(std::max(3, Nscalar), N, npTarget, *options);

  if (platform->options.compareArgs("BUILD ONLY", "FALSE")) {
    if (rank == 0) {
      printf("loading nek ... ");
      fflush(stdout);
    }

    std::string usrname;
    options->getArgs("CASENAME", usrname);
    std::string meshFile;
    options->getArgs("MESH FILE", meshFile);

    std::string cwd(fs::current_path());

    MPI_Fint nek_comm = MPI_Comm_c2f(platform->comm.mpiComm);

    set_usr_handles(usrname.c_str(), platform->verbose);

    (*nek_bootstrap_ptr)(&nek_comm,
                         (char *)cwd.c_str(),
                         /* basename */ (char *)usrname.c_str(),
                         (char *)meshFile.c_str(),
                         cwd.length(),
                         usrname.length(),
                         meshFile.length());
    if (rank == 0) {
      printf("done\n");
      fflush(stdout);
    }
  }
}

int setup(int numberActiveFields)
{
  MPI_Comm_rank(platform->comm.mpiComm, &rank);

  bool meshSolver = !options->compareArgs("MESH SOLVER", "NONE");

  std::string casename;
  options->getArgs("CASENAME", casename);

  int nscal;
  options->getArgs("NUMBER OF SCALARS", nscal);
  int nscalSolve = nscal;
  auto scalarCompute = (int *)calloc(nscal, sizeof(int));
  {
    int cnt = 0;
    for (int is = 0; is < nscal; is++) {
      std::string sid = scalarDigitStr(is);
      if (options->compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
        scalarCompute[is] = -1;
        nscalSolve--;
      }
    }
  }

  std::string velocitySolver;

  int meshPartType = 0; // RCB+RSB
  if (options->compareArgs("MESH PARTITIONER", "rcb")) {
    meshPartType = 1;
  }

  double meshConTol = 0.2;
  options->getArgs("MESH CONNECTIVITY TOL", meshConTol);

  int nBcRead = 1; // at the very least we have to read the boundaryIDs
  if (bcMap::useNekBCs()) {
    nBcRead += nscalSolve;
  }

  double rho;
  options->getArgs("DENSITY", rho);

  double mue;
  options->getArgs("VISCOSITY", mue);

  double rhoCp;
  options->getArgs("SCALAR00 DENSITY", rhoCp);

  double lambda;
  options->getArgs("SCALAR00 DIFFUSIVITY", lambda);

  int stressForm = 1; // avoid recompilation + bypass unligned SYM/SHL check

  // for now velocityExists is always true
  int velocityExists = (options->getArgs("VELOCITY SOLVER").empty()) ? 0 : 1;

  int nelgt, nelgv;
  re2::nelg(options->getArgs("MESH FILE"), nelgt, nelgv, platform->comm.mpiComm);
  const int cht = (nelgt > nelgv) && nscal;

  auto boundaryIDMap = [&](bool vMesh = false) {
    const std::string prefix = (cht && vMesh) ? "MESHV " : "MESH ";

    std::vector<std::string> list;
    options->getArgs(prefix + "BOUNDARY ID MAP", list, ",");

    std::vector<int> map;
    for (auto &entry : list) {
      map.push_back(std::stoi(entry));
    }
    return map;
  };

  auto bMapV = boundaryIDMap(true);
  int bMapVSize = bMapV.size();
  auto bMapT = boundaryIDMap();
  int bMapTSize = (cht) ? bMapT.size() : 0;

  (*nek_setup_ptr)(&velocityExists,
                   bMapV.data(),
                   &bMapVSize,
                   bMapT.data(),
                   &bMapTSize,
                   &nscal,
                   scalarCompute,
                   &nBcRead,
                   &meshPartType,
                   &meshConTol,
                   &rho,
                   &mue,
                   &rhoCp,
                   &lambda,
                   &stressForm);

  nekData.comm = MPI_Comm_f2c(*ptr<int>("nekcomm"));
  nekData.param = ptr<double>("param");

  nekData.ifield = ptr<int>("ifield");
  nekData.istep = ptr<int>("istep");
  nekData.time = ptr<double>("time");

  nekData.ndim = *ptr<int>("ndim");
  nekData.lelt = *ptr<int>("lelt");
  nekData.nelt = *ptr<int>("nelt");
  nekData.nelv = *ptr<int>("nelv");

  nekData.nx1 = *ptr<int>("nx1");
  nekData.ldimt = *ptr<int>("ldimt");

  nekData.glo_num = ptr<long long>("glo_num");

  nekData.p0th = ptr<double>("p0th");
  nekData.vx = ptr<double>("vx");
  nekData.vy = ptr<double>("vy");
  nekData.vz = ptr<double>("vz");
  nekData.pr = ptr<double>("pr");
  nekData.t = ptr<double>("t");
  nekData.wx = ptr<double>("wx");
  nekData.wy = ptr<double>("wy");
  nekData.wz = ptr<double>("wz");

  nekData.xm1 = ptr<double>("xm1");
  nekData.ym1 = ptr<double>("ym1");
  nekData.zm1 = ptr<double>("zm1");
  nekData.xc = ptr<double>("xc");
  nekData.yc = ptr<double>("yc");
  nekData.zc = ptr<double>("zc");

  nekData.unx = ptr<double>("unx");
  nekData.uny = ptr<double>("uny");
  nekData.unz = ptr<double>("unz");

  nekData.cbc = ptr<char>("cbc");
  nekData.boundaryID = ptr<int>("boundaryID");
  nekData.boundaryIDt = ptr<int>("boundaryIDt");

  if (bcMap::useNekBCs() && numberActiveFields > 0) {
    if (rank == 0) {
      printf("importing BCs from nek\n");
    }

    gen_bcmap();

    auto flow = true;
    if (platform->options.compareArgs("VELOCITY SOLVER", "NONE")) {
      flow = false;
    }

    if (flow) {
      if (rank == 0) {
        printf(" velocity\n");
        fflush(stdout);
      }

      int isTMesh = 0;
      int nIDs = (*nek_nbid_ptr)(&isTMesh);
      int *map = (int *)calloc(nIDs, sizeof(int));
      for (int id = 0; id < nIDs; id++) {
        map[id] = bcmap(id + 1, 1, 0);
      }
      bcMap::setBcMap("velocity", map, nIDs);

      if (meshSolver) {
        if (rank == 0) {
          printf(" mesh\n");
        }
        for (int id = 0; id < nIDs; id++) {
          map[id] = bcmap(id + 1, 1, 1);
        }
        bcMap::setBcMap("mesh", map, nIDs);
      }

      free(map);
    }
    for (int is = 0; is < nscal; is++) {
      std::string sid = scalarDigitStr(is);

      if (options->compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
        continue;
      }

      if (rank == 0) {
        printf(" scalar%02d\n", is);
        fflush(stdout);
      }

      int isTMesh = (cht && is == 0) ? 1 : 0;
      int nIDs = (*nek_nbid_ptr)(&isTMesh);

      int *map = (int *)calloc(nIDs, sizeof(int));

      for (int id = 0; id < nIDs; id++) {
        map[id] = bcmap(id + 1, is + 2, 0);
      }
      bcMap::setBcMap("scalar" + sid, map, nIDs);
      free(map);
    }
  }

  int isTMesh = 0;
  nekData.NboundaryID = (*nek_nbid_ptr)(&isTMesh);
  isTMesh = 1;
  nekData.NboundaryIDt = (*nek_nbid_ptr)(&isTMesh);

  dfloat startTime;
  options->getArgs("START TIME", startTime);
  *(nekData.time) = startTime;

  {
    hlong NelementsV = nekData.nelv;
    MPI_Allreduce(MPI_IN_PLACE, &NelementsV, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);
    nekrsCheck(NelementsV != nelgv, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "Invalid element partitioning");

    if (cht) {
      hlong NelementsT = nekData.nelt;
      MPI_Allreduce(MPI_IN_PLACE, &NelementsT, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);
      nekrsCheck(NelementsT <= NelementsV || NelementsT != nelgt,
                 MPI_COMM_SELF,
                 EXIT_FAILURE,
                 "%s\n",
                 "Invalid solid element partitioning");
    }
  }

  return 0;
}

long long set_glo_num(int nx, int isTMesh)
{
  return (*nek_set_vert_ptr)(&nx, &isTMesh);
}

void bdfCoeff(dfloat *g0, dfloat *coeff, dfloat *_dt, int order)
{
  double dt[3] = {_dt[0], _dt[1], _dt[2]};
  double nekCoeff[4];
  (*nek_setbd_ptr)(nekCoeff, dt, &order);
  *g0 = nekCoeff[0];
  for (int i = 0; i < order; i++) {
    coeff[i] = nekCoeff[i + 1];
  }
}

void extCoeff(dfloat *coeff, dfloat *_dt, int nAB, int nBDF)
{
  double dt[3] = {_dt[0], _dt[1], _dt[2]};
  double nekCoeff[3];
  (*nek_setabbd_ptr)(nekCoeff, dt, &nAB, &nBDF);
  for (int i = 0; i < nAB; i++) {
    coeff[i] = nekCoeff[i];
  }
}

void coeffAB(dfloat *coeff, dfloat *_dt, int order)
{
  double dt[3] = {_dt[0], _dt[1], _dt[2]};
  double nekCoeff[3];
  int one = 1;
  (*nek_setabbd_ptr)(nekCoeff, dt, &order, &one);
  for (int i = 0; i < order; i++) {
    coeff[i] = nekCoeff[i];
  }
}

void recomputeGeometry()
{
  (*nek_updggeom_ptr)();
}

void printMeshMetrics()
{
  (*nek_meshmetrics_ptr)();
}

const std::map<std::string, void *> &ptrList()
{
  return ptrListData;
}

} // namespace nek

// called from nek
extern "C" {

void nekf_registerptr(char *id, void *val, int *nameLen)
{
  auto name = std::string(id, *nameLen); // id comes from Fortran and is not null terminated
  auto entry = ptrListData.find(name);
  auto entryFound = (entry != ptrListData.end());
  nekrsCheck(entryFound && entry->second != val,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s exists already but is pointing to a different memory address\n",
             name.c_str());
  if (!entryFound) {
    const auto [it, success] = ptrListData.insert(std::make_pair(name, val));
    nekrsCheck(!success, MPI_COMM_SELF, EXIT_FAILURE, "Adding %s failed\n", name.c_str());
  }
}

void nekf_registerptr_(char *id, void *val, int *nameLen)
{
  nekf_registerptr(id, val, nameLen);
}
}
