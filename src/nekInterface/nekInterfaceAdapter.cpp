#include <unistd.h>
#include <fstream>
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "bcMap.hpp"
#include "neknek.hpp"
#include "fileUtils.hpp"
#include "re2Reader.hpp"
#include "fileUtils.hpp"

nekdata_private nekData;
static int rank;
static setupAide *options;
static nrs_t *nrs;

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

static void (*nek_ptr_ptr)(void **, char *, int *);
static void (*nek_scptr_ptr)(int *, void *);
static void (*nek_outfld_ptr)(char *, int);
static void (*nek_resetio_ptr)(void);
static void (*nek_setio_ptr)(double *, int *, int *, int *, int *, int *, int *);
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

static void (*nek_storesol_ptr)(void);
static void (*nek_restoresol_ptr)(void);
static void (*nek_updggeom_ptr)(void);
static void (*nek_meshmetrics_ptr)(void);

void noop_func(void) {}

void check_error(const char *error)
{
  nrsCheck(error != NULL, MPI_COMM_SELF, EXIT_FAILURE,
           "%s\n", error);
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

namespace nek {
void *ptr(const char *id)
{
  void *ptr;
  int len = strlen(id);
  (*nek_ptr_ptr)(&ptr, (char *)id, &len);
  return ptr;
}

void *scPtr(int id)
{
  nrsCheck(!id, MPI_COMM_SELF, EXIT_FAILURE,
           "Accessing invalid index %d!\n", id); 

  void *ptr;
  (*nek_scptr_ptr)(&id, &ptr);
  return ptr;
}

void outfld(const char *filename,
            double t,
            int step,
            int coords,
            int FP64,
            const occa::memory& o_uu,
            const occa::memory& o_pp,
            const occa::memory& o_ss,
            int NSfields)
{
  nrsCheck(NSfields > nekData.ldimt, platform->comm.mpiComm, EXIT_FAILURE,
           "Adjust ldimt in SIZE to %d and larger\n", NSfields);

  occa::memory o_u, o_p, o_s;
  if (o_uu.isInitialized()) {
    o_u = o_uu.cast(occa::dtype::byte);
  }
  if (o_pp.isInitialized()) {
    o_p = o_pp.cast(occa::dtype::byte);
  }
  if (o_ss.isInitialized() && NSfields) {
    o_s = o_ss.cast(occa::dtype::byte);
  }

  int xo = 0;
  int vo = 0;
  int po = 0;
  int so = 0;

  (*nek_storesol_ptr)();

  const int stepSave = *(nekData.istep);
  *(nekData.istep) = step;

  platform->timer.tic("checkpointing", 1);

  // nek5000 writes all fields using nelt
  // note, nrs->fieldOffset >= nelt*nxyz
  auto mesh = nrs->_mesh; 
  dlong Nlocal = mesh->Nelements * mesh->Np;

  // nek5000 uses lelv = lelt
  const dlong nekFieldOffset = nekData.lelt * mesh->Np;

  const auto p0thSave = *(nekData.p0th);
  *(nekData.p0th) = nrs->p0th[0];

  std::vector<dfloat> tmp1(Nlocal);
  std::vector<dfloat> tmp2(Nlocal);
  std::vector<dfloat> tmp3(Nlocal);

  if (coords) {
    mesh->o_x.copyTo(tmp1.data(), Nlocal);
    mesh->o_y.copyTo(tmp2.data(), Nlocal);
    mesh->o_z.copyTo(tmp3.data(), Nlocal);
    for(int i=0; i < Nlocal; i++) {
      nekData.xm1[i] = tmp1[i];   
      nekData.ym1[i] = tmp2[i];    
      nekData.zm1[i] = tmp3[i];    
    }
    xo = 1;
  }

  if (o_u.isInitialized()) {
    occa::memory o_vx = o_u + (0 * sizeof(dfloat)) * nrs->fieldOffset;
    occa::memory o_vy = o_u + (1 * sizeof(dfloat)) * nrs->fieldOffset;
    occa::memory o_vz = o_u + (2 * sizeof(dfloat)) * nrs->fieldOffset;
    o_vx.copyTo(tmp1.data(), Nlocal * sizeof(dfloat));
    o_vy.copyTo(tmp2.data(), Nlocal * sizeof(dfloat));
    o_vz.copyTo(tmp3.data(), Nlocal * sizeof(dfloat));
    for(int i=0; i < Nlocal; i++) {
      nekData.vx[i] = tmp1[i];   
      nekData.vy[i] = tmp2[i];    
      nekData.vz[i] = tmp3[i];    
    }
    vo = 1;
  }
  if (o_p.isInitialized()) {
    o_p.copyTo(tmp1.data(), Nlocal * sizeof(dfloat));
    for(int i=0; i < Nlocal; i++) {
      nekData.pr[i] = tmp1[i];   
    }
    po = 1;
  }
  if (o_s.isInitialized()) {
    for (int is = 0; is < NSfields; is++) {
      occa::memory o_Si = o_s + (is * sizeof(dfloat)) * nrs->fieldOffset;
      o_Si.copyTo(tmp1.data(), Nlocal * sizeof(dfloat));
      auto Ti = nekData.t + is * nekFieldOffset;
      for(int i=0; i < Nlocal; i++) {
        Ti[i] = tmp1[i];   
      }
    }
    so = 1;
  }

  (*nek_setio_ptr)(&t, &xo, &vo, &po, &so, &NSfields, &FP64);
  (*nek_outfld_ptr)((char *)filename, strlen(filename));
  (*nek_resetio_ptr)();

  platform->timer.toc("checkpointing");

  (*nek_restoresol_ptr)();
  *(nekData.p0th) = p0thSave; 

  *(nekData.istep) = stepSave;
}

void uic(int ifield) { (*nek_uic_ptr)(&ifield); }

void finalize() { (*nek_end_ptr)(); }

void restartFromFile(const std::string& fileName)
{
  std::string str(fileName.size(), '\0');
  std::replace_copy(fileName.begin(), fileName.end(), str.begin(), '+', ' ');
  int len = str.length();
  (*nek_restart_ptr)((char *)str.c_str(), &len);
}

void getIC(void)
{
  int Nscalar;
  options->getArgs("NUMBER OF SCALARS", Nscalar);

  for (int ifield = 1; ifield <= 1+Nscalar; ifield++)
    uic(ifield);

  double startTime;
  copyFromNek(startTime);
  platform->options.setArgs("START TIME", to_string_f(startTime));  
}

void xm1N(dfloat *_x, dfloat *_y, dfloat *_z, int N, dlong Nelements)
{
  const int Np = (N+1) * (N+1) * (N+1);
  const int nxyz = nekData.nx1 * nekData.nx1 * nekData.nx1;

  auto x = (double*) calloc(Np, sizeof(double));
  auto y = (double*) calloc(Np, sizeof(double));
  auto z = (double*) calloc(Np, sizeof(double));

  for(dlong e = 0; e < Nelements; ++e) {
    map_m_to_n(x, N+1, &nekData.xm1[e * nxyz], nekData.nx1); 
    map_m_to_n(y, N+1, &nekData.ym1[e * nxyz], nekData.nx1);
    map_m_to_n(z, N+1, &nekData.zm1[e * nxyz], nekData.nx1);
    for(int i = 0; i < Np; i++) {
      _x[i + e*Np] = x[i];
      _y[i + e*Np] = y[i];
      _z[i + e*Np] = z[i];
    }
  }
  free(x);
  free(y);
  free(z);
}

int lglel(int e)
{
  int ee = e + 1;
  return (*nek_lglel_ptr)(&ee) - 1;
}

void ifoutfld(int i) { (*nek_ifoutfld_ptr)(&i); }

void setics(void) { (*nek_setics_ptr)(); }

void userchk(void)
{
  if (rank == 0)
    printf("calling nek_userchk ...\n");
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
  if (platform->cacheBcast)
    cache_dir = platform->tmpDir;

  const std::string lib = cache_dir + "/nek5000/lib" + session_in + ".so";

  if(platform->comm.mpiRank == 0 && platform->verbose)
    std::cout << "\nloading " << lib << std::endl;
  void *handle = dlopen(lib.c_str(), RTLD_NOW | RTLD_LOCAL);

  nrsCheck(!handle, MPI_COMM_SELF, EXIT_FAILURE, 
           "%s\n", dlerror());

  // check if we need to append an underscore
  auto us = [handle]
  { 
    (void (*)(void))dlsym(handle, "usrdat_");
    if(handle)
      return "_";
    else
      return "";
  } ();
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

  nek_ptr_ptr = (void (*)(void **, char *, int *))dlsym(handle, fname("nekf_ptr"));
  check_error(dlerror());
  nek_scptr_ptr = (void (*)(int *, void *))dlsym(handle, fname("nekf_scptr"));
  check_error(dlerror());
  nek_bootstrap_ptr =
      (void (*)(int *, char *, char *, char *, int, int, int))dlsym(handle, fname("nekf_bootstrap"));
  check_error(dlerror());
  nek_setup_ptr =
      (void (*)(int *, int *, int *, int *, int *, double *, double *, double *, double *, double *, int *))
          dlsym(handle, fname("nekf_setup"));
  check_error(dlerror());
  nek_uic_ptr = (void (*)(int *))dlsym(handle, fname("nekf_uic"));
  check_error(dlerror());
  nek_end_ptr = (void (*)(void))dlsym(handle, fname("nekf_end"));
  check_error(dlerror());
  nek_outfld_ptr = (void (*)(char *, int))dlsym(handle, fname("nekf_outfld"));
  check_error(dlerror());
  nek_resetio_ptr = (void (*)(void))dlsym(handle, fname("nekf_resetio"));
  check_error(dlerror());
  nek_setio_ptr =
      (void (*)(double *, int *, int *, int *, int *, int *, int *))dlsym(handle, fname("nekf_setio"));
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

  nek_storesol_ptr = (void (*)(void))dlsym(handle, fname("nekf_storesol"));
  check_error(dlerror());
  nek_restoresol_ptr = (void (*)(void))dlsym(handle, fname("nekf_restoresol"));
  check_error(dlerror());
  nek_updggeom_ptr = (void (*)(void))dlsym(handle, fname("nekf_updggeom"));
  check_error(dlerror());
  nek_meshmetrics_ptr = (void (*)(void))dlsym(handle, fname("mesh_metrics"));
  check_error(dlerror());

#define postfix(x) x##_ptr
#define load_or_noop(s)                                                                                      \
do {                                                                                                         \
postfix(s) = (void (*)(void))dlsym(handle, fname(#s));                                                       \
if (!(postfix(s))) {                                                                                         \
postfix(s) = noop_func;                                                                                      \
if (verbose)                                                                                                 \
printf("Setting function " #s " to noop_func.\n");                                                           \
}                                                                                                            \
else if (verbose && rank == 0) {                                                                             \
printf("Loading " #s"\n");                                                                                   \
}                                                                                                            \
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
  sprintf(line, "%s/core/SIZE.template", nek5000_dir.c_str());
  FILE *fp = fopen(line, "r");
  nrsCheck(!fp, MPI_COMM_SELF, EXIT_FAILURE,
           "Cannot open %s!\n", line); 

  char *sizeFile;
  {
    fseek(fp, 0, SEEK_END);
    const long length = ftell(fp);
    rewind(fp);
    sizeFile = (char *)calloc(4096 + length, sizeof(char));
  }

  int lx1m = lx1;
  constexpr int nsessmax = 2;
  constexpr int nMaxObj = 20;

  int count = 0;
  while (fgets(line, lineSize, fp) != NULL) {
    if (strstr(line, "parameter (lx1=") != NULL)
      sprintf(line, "      parameter (lx1=%d)\n", lx1);
    else if (strstr(line, "parameter (lxd=") != NULL)
      sprintf(line, "      parameter (lxd=%d)\n", lxd);
    else if (strstr(line, "parameter (lelt=") != NULL)
      sprintf(line, "      parameter (lelt=%d)\n", lelt);
    else if (strstr(line, "parameter (lelg=") != NULL)
      sprintf(line, "      parameter (lelg=%d)\n", lelg);
    else if (strstr(line, "parameter (ldim=") != NULL)
      sprintf(line, "      parameter (ldim=%d)\n", ldim);
    else if (strstr(line, "parameter (lpmin=") != NULL)
      sprintf(line, "      parameter (lpmin=%d)\n", lpmin);
    else if (strstr(line, "parameter (ldimt=") != NULL)
      sprintf(line, "      parameter (ldimt=%d)\n", ldimt);
    else if (strstr(line, "parameter (mxprev=") != NULL)
      sprintf(line, "      parameter (mxprev=%d)\n", 1);
    else if (strstr(line, "parameter (lgmres=") != NULL)
      sprintf(line, "      parameter (lgmres=%d)\n", 1);
    else if (strstr(line, "parameter (lorder=") != NULL)
      sprintf(line, "      parameter (lorder=%d)\n", 1);
    else if (strstr(line, "parameter (lhis=") != NULL)
      sprintf(line, "      parameter (lhis=%d)\n", 100000);
    else if (strstr(line, "parameter (lelr=") != NULL)
      sprintf(line, "      parameter (lelr=%d)\n", 128 * lelt);
    else if (strstr(line, "parameter (lx1m=") != NULL)
      sprintf(line, "      parameter (lx1m=%d)\n", lx1m);
    else if(strstr(line, "parameter (nsessmax=") != NULL)
      sprintf(line, "      parameter (nsessmax=%d)\n", nsessmax);
    else if(strstr(line, "parameter (maxobj=") != NULL)
      sprintf(line, "      parameter (maxobj=%d)\n", nMaxObj);

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
        if (oldval < lelg)
          writeSize = 1;
      }
      if (line.find("lelt=") != std::string::npos) {
        int oldval;
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if (oldval < lelt)
          writeSize = 1;
      }
      if (line.find("lx1m=") != std::string::npos) {
        int oldval;
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if (oldval < lx1m)
          writeSize = 1;
      }
      if (line.find("lx1=") != std::string::npos) {
        int oldval;
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if (oldval != lx1)
          writeSize = 1;
      }
      if (line.find("ldimt=") != std::string::npos) {
        int oldval;
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if (oldval < ldimt)
          writeSize = 1;
      }
    }
  }
  osize.close();

  if (writeSize) {
    fp = fopen(SIZE, "w");
    fputs(sizeFile, fp);
    fclose(fp);
    free(sizeFile);
    if (verbose)
      printf("using new SIZE\n");
  }
  else {
    if (verbose)
      printf("using existing SIZE file %s/SIZE\n", cache_dir);
  }

  fflush(stdout);
}

void buildNekInterface(int ldimt, int N, int np, setupAide &options)
{
  int buildRank = rank;
  const bool buildNodeLocal = platform->cacheLocal;
  if (buildNodeLocal)
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);

  const int verbose = options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

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
  if (!usrFileExists)
    usrFile = nek5000_dir + "/core/zero.usr"; 

  int buildRequired = 0;
  if(platform->comm.mpiRank == 0) {
    if (platform->options.compareArgs("BUILD ONLY", "TRUE")) {
      buildRequired = 1;
    } else if (isFileNewer(usrFile.c_str(), libFile.c_str())) {
      buildRequired = 1;
    }
  }
  MPI_Bcast(&buildRequired, 1, MPI_INT, 0, platform->comm.mpiComm);

  if (platform->cacheBcast || platform->cacheLocal) {
    fileBcast(usrFile, platform->tmpDir, platform->comm.mpiComm, platform->verbose);
    usrFile = platform->tmpDir / fs::path(usrFile);
  }

  int err = [&]() {
    if (buildRank == 0) {
      mkdir(cache_dir.c_str(), S_IRWXU);

      int nelgt, nelgv;
      const int ndim = 3;
      const std::string meshFile = options.getArgs("MESH FILE");
      re2::nelg(meshFile, nelgt, nelgv, MPI_COMM_SELF);

      int lelt = (nelgt / np) + 3;
      if (lelt > nelgt)
        lelt = nelgt;

      const std::string sizeFile = cache_dir + "/SIZE";
      mkSIZE(N + 1, 1, lelt, nelgt, ndim, np, ldimt, options, sizeFile.c_str());

      if (buildRequired)
        copyFile(usrFile.c_str(), usrFileCache.c_str());

      if (isFileNewer(sizeFile.c_str(), libFile.c_str()))
        buildRequired = 1;

      if (buildRequired) {
        const double tStart = MPI_Wtime();

        const std::string makeOutput = (rank == 0) ? cache_dir + "/make.out" : "/dev/null"; 
        const std::string case_dir(fs::current_path());
        const std::string include_dirs = "./ " + case_dir + " " + installDir + "/include/bdry";
        const std::string nekInterface_dir = installDir + "/nekInterface";
        std::string make_args = "-j8 ";
        if (!verbose)
          make_args += "-s ";

        if (rank == 0)
          printf("building nekInterface for lx1=%d, lelt=%d and lelg=%d ... ", N + 1, lelt, nelgt);

        fflush(stdout);

        char buf[4096];
        sprintf(buf,
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
                std::string(">" + makeOutput + " 2>&1").c_str());

        if (verbose && rank == 0)
          printf("\n%s\n", buf);
        if (system(buf)) {
          if(rank == 0) {
            std::ifstream f(makeOutput);
            std::cerr << f.rdbuf();
          }

          return EXIT_FAILURE;
        }
        fileSync(libFile.c_str());

        if (rank == 0)
          printf("done (%gs)\n", MPI_Wtime() - tStart);
        fflush(stdout);
      }
      else {
        if (rank == 0)
          printf("skip building nekInterface (SIZE requires no update)\n");
        fflush(stdout);
      }
    } // buildRank

    if (platform->cacheBcast)
      fileBcast(cache_dir, platform->tmpDir, platform->comm.mpiComm, platform->verbose);

    return 0;
  }();

  nrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "");
}

namespace nek {
int bcmap(int bid, int ifld, int isMesh) { return (*nek_bcmap_ptr)(&bid, &ifld, &isMesh); }

void gen_bcmap() { (*nek_gen_bcmap_ptr)(); }

void bootstrap()
{
  options = &platform->options;

  int size;
  MPI_Comm_rank(platform->comm.mpiComm, &rank);
  MPI_Comm_size(platform->comm.mpiComm, &size);

  int buildRank = rank;
  const bool buildNodeLocal = platform->cacheLocal;
  if (buildNodeLocal)
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);

  int N;
  options->getArgs("POLYNOMIAL DEGREE", N);

  int Nscalar;
  options->getArgs("NUMBER OF SCALARS", Nscalar);

  int npTarget = size;
  options->getArgs("NP TARGET", npTarget);

  buildNekInterface(std::max(5, Nscalar), N, npTarget, *options);

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

int setup(nrs_t *nrs_in)
{
  nrs = nrs_in;
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
    for(int is = 0; is < nscal; is++) {
      std::string sid = scalarDigitStr(is);
      if (options->compareArgs("SCALAR" + sid + " SOLVER", "NONE")) {
        scalarCompute[is] = -1;
        nscalSolve--;
      }
    }
  }

  std::string velocitySolver;

  int meshPartType = 3; // RCB+RSB
  if (options->compareArgs("MESH PARTITIONER", "rcb"))
    meshPartType = 2;
  if (options->compareArgs("MESH PARTITIONER", "rcb+rsb"))
    meshPartType = 3;

  double meshConTol = 0.2;
  options->getArgs("MESH CONNECTIVITY TOL", meshConTol);

  int nBcRead = 1; // at the very least we have to read the boundaryIDs
  if (bcMap::useNekBCs())
    nBcRead += nscalSolve;

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

  (*nek_setup_ptr)(&velocityExists,
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

  nekData.param = (double *)ptr("param");
  nekData.ifield = (int *)ptr("ifield");
  nekData.istep = (int *)ptr("istep");
  nekData.time = (double *)ptr("time");

  nekData.ndim = *(int *)ptr("ndim");
  nekData.nelt = *(int *)ptr("nelt");
  nekData.nelv = *(int *)ptr("nelv");
  nekData.lelt = *(int *)ptr("lelt");
  nekData.ldimt = *(int *)ptr("ldimt");
  nekData.nx1 = *(int *)ptr("nx1");

  nekData.vx = (double *)ptr("vx");
  nekData.vy = (double *)ptr("vy");
  nekData.vz = (double *)ptr("vz");
  nekData.pr = (double *)ptr("pr");
  nekData.t = (double *)ptr("t");

  nekData.qtl = (double *)ptr("qtl");

  nekData.ifgetu = (int *)ptr("ifgetu");
  nekData.ifgetp = (int *)ptr("ifgetp");
  nekData.ifgett = (int *)ptr("ifgett");
  nekData.ifgetps = (int *)ptr("ifgetps");

  nekData.unx = (double *)ptr("unx");
  nekData.uny = (double *)ptr("uny");
  nekData.unz = (double *)ptr("unz");

  nekData.xm1 = (double *)ptr("xm1");
  nekData.ym1 = (double *)ptr("ym1");
  nekData.zm1 = (double *)ptr("zm1");
  nekData.xc = (double *)ptr("xc");
  nekData.yc = (double *)ptr("yc");
  nekData.zc = (double *)ptr("zc");

  nekData.glo_num = (long long *)ptr("glo_num");
  nekData.cbc = (char *)ptr("cbc");

  nekData.boundaryID = (int *)ptr("boundaryID");
  nekData.boundaryIDt = (int *)ptr("boundaryIDt");

  nekData.eface1 = (int *)ptr("eface1");
  nekData.eface = (int *)ptr("eface");
  nekData.icface = (int *)ptr("icface");
  nekData.comm = MPI_Comm_f2c(*(int *)ptr("nekcomm"));

  nekData.p0th = (double *)ptr("p0th");
  nekData.dp0thdt = (double *)ptr("dp0thdt");

  nekData.wx = (double *)ptr("wx");
  nekData.wy = (double *)ptr("wy");
  nekData.wz = (double *)ptr("wz");

  int cht = 0;
  if (nekData.nelv != nekData.nelt && nscal)
    cht = 1;

  if (bcMap::useNekBCs() && numberActiveFields(nrs) > 0) {
    if (rank == 0)
      printf("importing BCs from nek\n");

    gen_bcmap();

    if (nrs->flow) {
      if (rank == 0) printf(" velocity\n");

      int isTMesh = 0;
      int nIDs = (*nek_nbid_ptr)(&isTMesh);
      int *map = (int *)calloc(nIDs, sizeof(int));
      for (int id = 0; id < nIDs; id++)
        map[id] = bcmap(id + 1, 1, 0);
      bcMap::setBcMap("velocity", map, nIDs);

      if (meshSolver) {
        if (rank == 0) printf(" mesh\n");
        for (int id = 0; id < nIDs; id++)
          map[id] = bcmap(id + 1, 1, 1);
        bcMap::setBcMap("mesh", map, nIDs);
      }

      free(map);
    }
    for(int is = 0; is < nscal; is++) {
      std::string sid = scalarDigitStr(is);

      if (options->compareArgs("SCALAR" + sid + " SOLVER", "NONE"))
        continue;

      if (rank == 0) printf(" scalar%02d\n", is);

      int isTMesh = 0;
      if (cht && is == 0)
        isTMesh = 1;
      int nIDs = (*nek_nbid_ptr)(&isTMesh);

      int *map = (int *)calloc(nIDs, sizeof(int));

      for (int id = 0; id < nIDs; id++)
        map[id] = bcmap(id + 1, is + 2, 0);
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

  return 0;
}

static void updateMesh()
{
  auto mesh = nrs->_mesh;

  for(int i=0; i < mesh->Nlocal; i++) {
    nekData.xm1[i] = mesh->x[i];
    nekData.ym1[i] = mesh->y[i];
    nekData.zm1[i] = mesh->z[i];
  }
  recomputeGeometry();
}

void copyToNek(double time)
{
  if (rank == 0) {
    printf("copying solution to nek\n");
    fflush(stdout);
  }

  if (*(nekData.istep) == 0) {
    updateMesh();
  }

  mesh_t *mesh = nrs->meshV;
  const dlong Nlocal = mesh->Nelements * mesh->Np;

  auto vx = nrs->U + 0 * nrs->fieldOffset;
  auto vy = nrs->U + 1 * nrs->fieldOffset;
  auto vz = nrs->U + 2 * nrs->fieldOffset;

  *(nekData.time) = time;
  *(nekData.p0th) = nrs->p0th[0];
  *(nekData.dp0thdt) = nrs->dp0thdt;

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    mesh_t *mesh = nrs->meshV;
    if (nrs->cht)
      mesh = nrs->cds->mesh[0];
    const dlong Nlocal = mesh->Nelements * mesh->Np;
    auto wx = mesh->U + 0 * nrs->fieldOffset;
    auto wy = mesh->U + 1 * nrs->fieldOffset;
    auto wz = mesh->U + 2 * nrs->fieldOffset;
    for(int i = 0; i < Nlocal; i++) {
      nekData.wx[i] = wx[i];
      nekData.wy[i] = wy[i];
      nekData.wz[i] = wz[i];
    }
    updateMesh();
  }

  for(int i = 0; i < Nlocal; i++) {
    nekData.vx[i] = vx[i];
    nekData.vy[i] = vy[i];
    nekData.vz[i] = vz[i];
    nekData.pr[i] = nrs->P[i];
  }
  if (nrs->Nscalar) {
    const dlong nekFieldOffset = nekData.lelt * mesh->Np;
    for (int is = 0; is < nrs->Nscalar; is++) {
      mesh_t *mesh;
      (is) ? mesh = nrs->cds->meshV : mesh = nrs->cds->mesh[0];
      const dlong Nlocal = mesh->Nelements * mesh->Np;
      auto Ti = nekData.t + is * nekFieldOffset;
      auto Si = nrs->cds->S + nrs->cds->fieldOffsetScan[is];
      for(int i = 0; i < Nlocal; i++) {
        Ti[i] = Si[i];
      }
    }
  }
}

void ocopyToNek(void)
{
  nrs->o_U.copyTo(nrs->U);
  nrs->o_P.copyTo(nrs->P);
  if (nrs->Nscalar) {
    nrs->cds->o_S.copyTo(nrs->cds->S);
  }
  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    mesh_t *mesh = nrs->meshV;
    if (nrs->cht)
      mesh = nrs->cds->mesh[0];
    mesh->o_U.copyTo(mesh->U);
    mesh->o_x.copyTo(mesh->x);
    mesh->o_y.copyTo(mesh->y);
    mesh->o_z.copyTo(mesh->z);
  }
  copyToNek(0.0);
}

void ocopyToNek(double time, int tstep)
{
  if (tstep == 0) {
    auto mesh = nrs->_mesh;
    mesh->o_x.copyTo(mesh->x);
    mesh->o_y.copyTo(mesh->y);
    mesh->o_z.copyTo(mesh->z);
  }

  nrs->o_U.copyTo(nrs->U);
  nrs->o_P.copyTo(nrs->P);
  if (nrs->Nscalar) {
    nrs->cds->o_S.copyTo(nrs->cds->S);
  }
  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    mesh_t *mesh = nrs->meshV;
    if (nrs->cht)
      mesh = nrs->cds->mesh[0];
    mesh->o_U.copyTo(mesh->U);
    mesh->o_x.copyTo(mesh->x);
    mesh->o_y.copyTo(mesh->y);
    mesh->o_z.copyTo(mesh->z);
  }
  copyToNek(time, tstep);
}

void copyToNek(double time, int tstep)
{
  *(nekData.istep) = tstep;
  copyToNek(time);
}

void ocopyFromNek(double &time)
{
  copyFromNek(time);
  nrs->o_P.copyFrom(nrs->P);
  nrs->o_U.copyFrom(nrs->U);
  if (nrs->Nscalar) {
    nrs->cds->o_S.copyFrom(nrs->cds->S);
  }
  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    mesh_t *mesh = nrs->meshV;
    if (nrs->cht)
      mesh = nrs->cds->mesh[0];
    mesh->o_x.copyFrom(mesh->x);
    mesh->o_y.copyFrom(mesh->y);
    mesh->o_z.copyFrom(mesh->z);
    mesh->o_U.copyFrom(mesh->U);
  }
}

void copyFromNek(double &time)
{
  if (rank == 0) {
    printf("copying solution from nek\n");
    fflush(stdout);
  }

  mesh_t *mesh = nrs->meshV;
  dlong Nlocal = mesh->Nelements * mesh->Np;

  time = *(nekData.time);
  nrs->p0th[0] = *(nekData.p0th);
  nrs->dp0thdt = *(nekData.dp0thdt);

  auto vx = nrs->U + 0 * nrs->fieldOffset;
  auto vy = nrs->U + 1 * nrs->fieldOffset;
  auto vz = nrs->U + 2 * nrs->fieldOffset;

  for(int i=0; i < Nlocal; i++) {
    vx[i] = nekData.vx[i];
    vy[i] = nekData.vy[i];
    vz[i] = nekData.vz[i];
  }

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    mesh_t *mesh = nrs->meshV;
    if (nrs->cht)
      mesh = nrs->cds->mesh[0];

    const dlong Nlocal = mesh->Nelements * mesh->Np;

    auto wx = mesh->U + 0 * nrs->fieldOffset;
    auto wy = mesh->U + 1 * nrs->fieldOffset;
    auto wz = mesh->U + 2 * nrs->fieldOffset;

    for(int i=0; i < Nlocal; i++) {
      wx[i] = nekData.wx[i];
      wy[i] = nekData.wy[i];
      wz[i] = nekData.wz[i];

      nekData.xm1[i] = mesh->x[i];
      nekData.ym1[i] = mesh->y[i];
      nekData.zm1[i] = mesh->z[i];
    }
    recomputeGeometry();
  }

  for(int i=0; i < Nlocal; i++) nrs->P[i] = nekData.pr[i]; 

  if (nrs->Nscalar) {
    const dlong nekFieldOffset = nekData.lelt * mesh->Np;
    for (int is = 0; is < nrs->Nscalar; is++) {
      mesh_t *mesh;
      (is) ? mesh = nrs->cds->meshV : mesh = nrs->cds->mesh[0];
      const dlong Nlocal = mesh->Nelements * mesh->Np;
      auto Ti = nekData.t + is * nekFieldOffset;
      auto Si = nrs->cds->S + nrs->cds->fieldOffsetScan[is];
      for(int i=0; i < Nlocal; i++) Si[i] = Ti[i]; 
    }
  }
}

long long set_glo_num(int nx, int isTMesh) { return (*nek_set_vert_ptr)(&nx, &isTMesh); }

void bdfCoeff(dfloat *g0, dfloat *coeff, dfloat *_dt, int order)
{
  double dt[3] = {_dt[0], _dt[1], _dt[2]};
  double nekCoeff[4];
  (*nek_setbd_ptr)(nekCoeff, dt, &order);
  *g0 = nekCoeff[0];
  for(int i=0; i < order; i++)
    coeff[i] = nekCoeff[i+1];
}

void extCoeff(dfloat *coeff, dfloat *_dt, int nAB, int nBDF) 
{ 
  double dt[3] = {_dt[0], _dt[1], _dt[2]};
  double nekCoeff[3];
  (*nek_setabbd_ptr)(nekCoeff, dt, &nAB, &nBDF); 
  for(int i=0; i < nAB; i++)
    coeff[i] = nekCoeff[i] ;
}

void coeffAB(dfloat *coeff, dfloat *_dt, int order)
{
  double dt[3] = {_dt[0], _dt[1], _dt[2]};
  double nekCoeff[3];
  int one = 1;
  (*nek_setabbd_ptr)(nekCoeff, dt, &order, &one);
  for(int i=0; i < order; i++)
    coeff[i] = nekCoeff[i]; 
}

void recomputeGeometry() { (*nek_updggeom_ptr)(); }

void printMeshMetrics() { (*nek_meshmetrics_ptr)(); }

} // namespace nek
