#include <unistd.h>
#include <fstream>
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "bcMap.hpp"

nekdata_private nekData;
static int rank;
static setupAide* options;
static ins_t** ins;
static dfloat timeLast = -1;

static void (* usrdat_ptr)(void);
static void (* usrdat2_ptr)(void);
static void (* usrdat3_ptr)(void);
static void (* userchk_ptr)(void);
static void (* uservp_ptr)(void);
static void (* userf_ptr)(void);
static void (* userq_ptr)(void);
static void (* userbc_ptr)(void);
static void (* useric_ptr)(void);
static void (* userqtl_ptr)(void);
static void (* usrsetvert_ptr)(void);

static void (* nek_ptr_ptr)(void**, char*, int*);
static void (* nek_outfld_ptr)(char*);
static void (* nek_resetio_ptr)(void);
static void (* nek_setio_ptr)(double*, int*, int*, int*, int*, int*, int*);
static void (* nek_uic_ptr)(int*);
static void (* nek_end_ptr)(void);
static void (* nek_restart_ptr)(char*, int*);
static void (* nek_map_m_to_n_ptr)(double* a, int* na, double* b, int* nb, int* if3d,
                                   double* w, int* nw);
static void (* nek_outpost_ptr)(double* v1, double* v2, double* v3, double* vp,
                                double* vt, char* name, int);
static void (* nek_uf_ptr)(double*, double*, double*);
static int (* nek_lglel_ptr)(int*);
static void (* nek_setup_ptr)(int*, char*, char*, int*, int*, int*, int*, int, int);
static void (* nek_ifoutfld_ptr)(int*);
static void (* nek_setics_ptr)(void);
static int (* nek_bcmap_ptr)(int*, int*);
static void (* nek_gen_bcmap_ptr)(void);
static int (* nek_nbid_ptr)(int*);
static long long (* nek_set_vert_ptr)(int*, int*);

void noop_func(void) {}

void* nek_ptr(const char* id)
{
  void* ptr;
  int len = strlen(id);
  (*nek_ptr_ptr)(&ptr, (char*)id, &len);
  return ptr;
}

void nek_outfld()
{
  const char suffix[] = "   ";
  (*nek_outfld_ptr)((char*)suffix);
}

void nek_outfld(const char* suffix)
{
  (*nek_outfld_ptr)((char*)suffix);
}

void nek_outfld(const char* suffix, dfloat t, int coords,
                occa::memory o_u, occa::memory o_p, occa::memory o_s,
                int NSfields, int FP64)
{
  mesh_t* mesh = (*ins)->mesh;
  cds_t* cds = (*ins)->cds;
  dlong Nlocal = mesh->Nelements * mesh->Np;

  double time = t;

  int xo = 0;
  int vo = 0;
  int po = 0;
  int so = 0;

  if(coords)
    xo = 1;
  if(o_u.ptr()) {
    occa::memory o_vx = o_u + 0 * (*ins)->fieldOffset * sizeof(dfloat);
    occa::memory o_vy = o_u + 1 * (*ins)->fieldOffset * sizeof(dfloat);
    occa::memory o_vz = o_u + 2 * (*ins)->fieldOffset * sizeof(dfloat);
    o_vx.copyTo(nekData.vx, Nlocal * sizeof(dfloat));
    o_vy.copyTo(nekData.vy, Nlocal * sizeof(dfloat));
    o_vz.copyTo(nekData.vz, Nlocal * sizeof(dfloat));
    vo = 1;
  }
  if(o_p.ptr()) {
    o_p.copyTo(nekData.pr, Nlocal * sizeof(dfloat));
    po = 1;
  }
  if(o_s.ptr() && NSfields) {
    const dlong nekFieldOffset = nekData.lelt * mesh->Np;
    for(int is = 0; is < NSfields; is++) {
      mesh_t* mesh;
      (is) ? mesh = (*ins)->cds->meshV : mesh = (*ins)->cds->mesh;
      const dlong Nlocal = mesh->Nelements * mesh->Np;
      dfloat* Ti = nekData.t + is * nekFieldOffset;
      occa::memory o_Si = o_s + is * cds->fieldOffset * sizeof(dfloat);
      o_Si.copyTo(Ti, Nlocal * sizeof(dfloat));
    }
    so = 1;
  }

  (*nek_setio_ptr)(&t, &xo, &vo, &po, &so, &NSfields, &FP64);
  (*nek_outfld_ptr)((char*)suffix);
  (*nek_resetio_ptr)();
}

void nek_uic(int ifield)
{
  (*nek_uic_ptr)(&ifield);
}

void nek_end()
{
  (*nek_end_ptr)();
}

void nek_setic(void)
{
  int readRestartFile;
  options->getArgs("RESTART FROM FILE", readRestartFile);

  if (readRestartFile) {
    std::string str1;
    options->getArgs("RESTART FILE NAME", str1);
    std::string str2(str1.size(), '\0');
    std::replace_copy(str1.begin(), str1.end(), str2.begin(), '+', ' ');
    int len = str2.length();
    (*nek_restart_ptr)((char*)str2.c_str(), &len);
  }

  (*nek_setics_ptr)();

  if (readRestartFile) {
    double startTime = *(nekData.time);
    options->setArgs("START TIME", to_string_f(startTime));
  }
}

void nek_map_m_to_n(double* a, int na, double* b, int nb)
{
  // output a, intput b
  int if3d = (nekData.ndim == 3);

  int maxn = na > nb ? na : nb;
  int N = 2 * maxn * maxn * maxn;
  double* w = (double*) malloc(sizeof(double) * N);

  (*nek_map_m_to_n_ptr)(a, &na, b, &nb, &if3d, w, &N);
  free(w);
}

void nek_uf(double* u, double* v, double* w)
{
  (*nek_uf_ptr)(u, v, w);
}

int nek_lglel(int e)
{
  int ee = e + 1;
  return (*nek_lglel_ptr)(&ee) - 1;
}

void nek_ifoutfld(int i)
{
  (*nek_ifoutfld_ptr)(&i);
}

void nek_setics(void)
{
  (*nek_setics_ptr)();
}

void nek_userchk(void)
{
  if(rank == 0) printf("calling nek_userchk\n");
  (*userchk_ptr)();
}

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

void check_error(char* error)
{
  if(error != NULL) {
    fprintf(stderr, "Error: %s\n", error);
    ABORT(EXIT_FAILURE);
  }
}

void set_function_handles(const char* session_in,int verbose)
{
  // load lib{session_in}.so
  char lib_session[BUFSIZ], * error;

  const char* cache_dir = getenv("NEKRS_CACHE_DIR");
  sprintf(lib_session, "%s/lib%s.so", cache_dir, session_in);

  void* handle = dlopen(lib_session,RTLD_NOW | RTLD_GLOBAL);
  if(!handle) {
    fprintf(stderr, "%s\n", dlerror());
    ABORT(EXIT_FAILURE);
  }

  // check if we need to append an underscore
  char us[2] = "";
  char func[20];
  usrdat_ptr = (void (*)(void))dlsym(handle, "usrdat_");
  if (usrdat_ptr) strcpy(us,"_");
  dlerror();

#define fname(s) (strcpy(func,(s)), strcat(func, us), func)

  usrdat_ptr = (void (*)(void))dlsym(handle, fname("usrdat"));
  check_error(dlerror());
  usrdat2_ptr = (void (*)(void))dlsym(handle, fname("usrdat2"));
  check_error(dlerror());
  usrdat3_ptr = (void (*)(void))dlsym(handle, fname("usrdat3"));
  check_error(dlerror());
  userchk_ptr = (void (*)(void))dlsym(handle, fname("userchk"));
  check_error(dlerror());

  nek_ptr_ptr = (void (*)(void**, char*, int*))dlsym(handle, fname("nekf_ptr"));
  check_error(dlerror());
  nek_setup_ptr =
    (void (*)(int*, char*, char*, int*, int*, int*, int*, int, int))dlsym(handle, fname("nekf_setup"));
  check_error(dlerror());
  nek_uic_ptr = (void (*)(int*))dlsym(handle, fname("nekf_uic"));
  check_error(dlerror());
  nek_end_ptr = (void (*)(void))dlsym(handle, fname("nekf_end"));
  check_error(dlerror());
  nek_outfld_ptr = (void (*)(char*))dlsym(handle, fname("nekf_outfld"));
  check_error(dlerror());
  nek_resetio_ptr = (void (*)(void))dlsym(handle, fname("nekf_resetio"));
  check_error(dlerror());
  nek_setio_ptr =
    (void (*)(double*, int*, int*, int*, int*, int*, int*))dlsym(handle, fname("nekf_setio"));
  check_error(dlerror());
  nek_restart_ptr = (void (*)(char*, int*))dlsym(handle, fname("nekf_restart"));
  check_error(dlerror());
  check_error(dlerror());
  nek_uf_ptr = (void (*)(double*, double*, double*))dlsym(handle,
                                                          fname("nekf_uf"));
  check_error(dlerror());
  nek_lglel_ptr = (int (*)(int*))dlsym(handle,fname("nekf_lglel"));
  check_error(dlerror());
  nek_ifoutfld_ptr = (void (*)(int*))dlsym(handle,fname("nekf_ifoutfld"));
  check_error(dlerror());
  nek_setics_ptr = (void (*)(void))dlsym(handle,fname("nekf_setics"));
  check_error(dlerror());
  nek_bcmap_ptr = (int (*)(int*, int*))dlsym(handle,fname("nekf_bcmap"));
  check_error(dlerror());
  nek_gen_bcmap_ptr = (void (*)(void))dlsym(handle,fname("nekf_gen_bcmap"));
  check_error(dlerror());
  nek_map_m_to_n_ptr = (void (*)(double*, int*, double*, int*, int*, double*, int*)) \
                       dlsym(handle, fname("map_m_to_n"));
  check_error(dlerror());
  nek_nbid_ptr = (int (*)(int*))dlsym(handle,fname("nekf_nbid"));
  check_error(dlerror());
  nek_set_vert_ptr = (long long (*)(int*, int*))dlsym(handle,fname("nekf_set_vert"));
  check_error(dlerror());

#define postfix(x) x ## _ptr
#define load_or_noop(s) \
  do { \
    postfix(s) = (void (*)(void))dlsym(handle,fname(#s)); \
    if(!(postfix(s))) { \
      postfix(s) = noop_func; \
      if(verbose) printf("Setting function " #s " to noop_func.\n"); \
    } else if(verbose) { \
      printf("Loading " #s " from lib%s.so\n",session_in); \
    } \
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

void mkSIZE(int lx1, int lxd, int lelt, int lelg, int ldim, int lpmin, int ldimt)
{
  //printf("generating SIZE file ... "); fflush(stdout);

  char line[BUFSIZ];
  char cmd[BUFSIZ];

  const char* cache_dir = getenv("NEKRS_CACHE_DIR");
  const char* nekrs_nek5000_dir = getenv("NEKRS_NEK5000_DIR");

  // Read and generate the new size file.
  sprintf(line,"%s/core/SIZE.template", nekrs_nek5000_dir);
  FILE* fp = fopen(line, "r");
  char* sizeFile, * curSizeFile;
  size_t result;

  if (fp) {
    fseek(fp, 0, SEEK_END);
    long length = ftell(fp);
    rewind(fp);
    // allocate actual length + some buffer
    sizeFile = (char*) calloc(length + 500,sizeof(char));
    if(!sizeFile) {
      fprintf(stderr, "Error allocating space for SIZE file.\n");
      ABORT(EXIT_FAILURE);
    }
  } else {
    fprintf(stderr, "Error opening %s/core/SIZE.template!\n", nekrs_nek5000_dir);
    ABORT(EXIT_FAILURE);
  }

  int count = 0;
  while(fgets(line, BUFSIZ, fp) != NULL) {
    if(strstr(line, "parameter (lx1=") != NULL)
      sprintf(line, "      parameter (lx1=%d)\n", lx1);
    else if(strstr(line, "parameter (lxd=") != NULL)
      sprintf(line, "      parameter (lxd=%d)\n", lxd);
    else if(strstr(line, "parameter (lelt=") != NULL)
      sprintf(line, "      parameter (lelt=%d)\n", lelt);
    else if(strstr(line, "parameter (lelg=") != NULL)
      sprintf(line, "      parameter (lelg=%ld)\n", lelg);
    else if(strstr(line, "parameter (ldim=") != NULL)
      sprintf(line, "      parameter (ldim=%d)\n", ldim);
    else if(strstr(line, "parameter (lpmin=") != NULL)
      sprintf(line, "      parameter (lpmin=%d)\n", lpmin);
    else if(strstr(line, "parameter (ldimt=") != NULL)
      sprintf(line, "      parameter (ldimt=%d)\n", ldimt);
    else if(strstr(line, "parameter (mxprev=") != NULL)
      sprintf(line, "      parameter (mxprev=%d)\n", 1);
    else if(strstr(line, "parameter (lgmres=") != NULL)
      sprintf(line, "      parameter (lgmres=%d)\n", 1);
    else if(strstr(line, "parameter (lorder=") != NULL)
      sprintf(line, "      parameter (lorder=%d)\n", 1);
    else if(strstr(line, "parameter (lelr=") != NULL)
      sprintf(line, "      parameter (lelr=%d)\n", 128 * lelt);

    strcpy(sizeFile + count, line);
    count += strlen(line);
  }
  fclose(fp);

  int writeSize = 1;

  // read size if exists
  std::ifstream osize;
  sprintf(line,"%s/SIZE", cache_dir);

  osize.open(line, std::ifstream::in);
  if(osize.is_open()) {
    writeSize = 0;
    string line;
    int oldval;
    while(getline( osize, line )) {
      if(line.find( "lelg=") != string::npos ) {
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if(oldval < lelg) writeSize = 1;
      }
      if(line.find( "lelt=") != string::npos ) {
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if(oldval < lelt) writeSize = 1;
      }
      if(line.find( "lx1=") != string::npos ) {
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if(oldval != lx1) writeSize = 1;
      }
      if(line.find( "ldimt=") != string::npos ) {
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if(oldval < ldimt) writeSize = 1;
      }
    }
  }
  osize.close();

  if(writeSize) {
    fp = fopen(line, "w");
    fputs(sizeFile, fp);
    fclose(fp);
    free(sizeFile);
    //printf("done\n");
  } else {
    //printf("using existing SIZE file %s/SIZE\n", cache_dir);
  }

  fflush(stdout);
}

int buildNekInterface(const char* casename, int ldimt, int N, int np)
{
  printf("building nek ... ");
  fflush(stdout);

  char buf[BUFSIZ];
  char fflags[BUFSIZ];
  char cflags[BUFSIZ];

  const char* cache_dir = getenv("NEKRS_CACHE_DIR");
  const char* nekInterface_dir = getenv("NEKRS_NEKINTERFACE_DIR");
  const char* nek5000_dir = getenv("NEKRS_NEK5000_DIR");

  FILE* fp;
  int retval;

  sprintf(buf, "%s.re2", casename);
  fp = fopen(buf, "r");
  if (!fp) {
    printf("\nERROR: Cannot find %s!\n", buf);
    ABORT(EXIT_FAILURE);
  }
  fgets(buf, 80, fp);
  fclose(fp);

  char ver[10];
  int nelgv, nelgt, ndim;
  sscanf(buf, "%5s %9d %1d %9d", ver, &nelgt, &ndim, &nelgv);
  int lelt = nelgt / np + 3;
  if(lelt > nelgt) lelt = nelgt;
  mkSIZE(N + 1, 1, lelt, nelgt, ndim, np, ldimt);

  // Copy case.usr file to cache_dir
  sprintf(buf,"%s.usr",casename);
  if(access(buf,F_OK) != -1)
    sprintf(buf, "cp -pf %s.usr %s",casename,cache_dir);
  else
    sprintf(buf, "cp -pf %s/core/zero.usr %s/%s.usr",nek5000_dir,cache_dir,casename);
  retval = system(buf);
  if (retval) goto err;

  // Copy Nek5000/core from install_dir to cache_dir
  sprintf(buf, "cp -pr %s %s", nek5000_dir, cache_dir); 
  retval = system(buf);
  if (retval) goto err;

  //TODO: Fix hardwired compiler flags
  sprintf(fflags, "\"${NEKRS_FFLAGS} -mcmodel=medium -fPIC -fcray-pointer -I../ \"");
  sprintf(cflags, "\"${NEKRS_CXXFLAGS} -fPIC -I${NEKRS_NEKINTERFACE_DIR}\"");

  sprintf(buf, "cd %s && yes n | FC=\"${NEKRS_FC}\" CC=\"${NEKRS_CC}\" FFLAGS=%s "
          "CFLAGS=%s PPLIST=\"${NEKRS_NEK5000_PPLIST}\" NEK_SOURCE_ROOT=%s/nek5000 "
          "%s/nek5000/bin/nekconfig %s >build.log 2>&1", cache_dir, fflags,
          cflags, cache_dir, cache_dir, casename);
  retval = system(buf);
  //if (retval) goto err;
  sprintf(buf, "cd %s && NEKRS_WORKING_DIR=%s make -j4 -f %s/Makefile nekInterface "
          ">>build.log 2>&1", cache_dir, cache_dir, nekInterface_dir);
  retval = system(buf);
  if (retval) goto err;

  printf("done\n\n");
  fflush(stdout);
  return 0;

err:
  printf("\nAn ERROR occured, see %s/build.log for details!\n", cache_dir);
  ABORT(EXIT_FAILURE);
}

int nek_bcmap(int bid, int ifld)
{
  return (*nek_bcmap_ptr)(&bid, &ifld);
}

void nek_gen_bcmap()
{
  (*nek_gen_bcmap_ptr)();
}

int nek_setup(MPI_Comm c, setupAide &options_in, ins_t** ins_in)
{
  options = &options_in;
  ins = ins_in;
  MPI_Comm_rank(c,&rank);
  MPI_Fint nek_comm = MPI_Comm_c2f(c);

  string casename;
  options->getArgs("CASENAME", casename);

  char buf[FILENAME_MAX];
  getcwd(buf, sizeof(buf));
  string cwd;
  cwd.assign(buf);

  set_function_handles(casename.c_str(), 0);

  int nscal = 0;
  options->getArgs("NUMBER OF SCALARS", nscal);

  string velocitySolver;
  int flow = 1;
  if(options->compareArgs("VELOCITY", "FALSE")) flow = 0;

  int meshPartType = 1; // RSB
  if(options->compareArgs("MESH PARTITIONER", "RCB")) meshPartType = 2;
  if(options->compareArgs("MESH PARTITIONER", "RCB+RSB")) meshPartType = 3;

  int nBcRead = 1;
  int bcInPar = 1;
  if(bcMap::size(0) == 0 && bcMap::size(1) == 0) {
    bcInPar = 0;
    nBcRead = flow + nscal;
  }

  (*nek_setup_ptr)(&nek_comm, (char*)cwd.c_str(), (char*)casename.c_str(),
                   &flow, &nscal, &nBcRead, &meshPartType, 
                   cwd.length(), casename.length()); 

  nekData.param = (double*) nek_ptr("param");
  nekData.ifield = (int*) nek_ptr("ifield");
  nekData.istep = (int*) nek_ptr("istep");
  nekData.time = (double*) nek_ptr("time");

  nekData.ndim = *(int*) nek_ptr("ndim");
  nekData.nelt = *(int*) nek_ptr("nelt");
  nekData.nelv = *(int*) nek_ptr("nelv");
  nekData.lelt = *(int*) nek_ptr("lelt");
  nekData.nx1 =  *(int*) nek_ptr("nx1");

  nekData.vx = (double*) nek_ptr("vx");
  nekData.vy = (double*) nek_ptr("vy");
  nekData.vz = (double*) nek_ptr("vz");
  nekData.pr = (double*) nek_ptr("pr");
  nekData.t  = (double*) nek_ptr("t");

  nekData.qtl = (double*) nek_ptr("qtl");

  nekData.ifgetu  = (int*) nek_ptr("ifgetu");
  nekData.ifgetp  = (int*) nek_ptr("ifgetp");
  nekData.ifgett  = (int*) nek_ptr("ifgett");
  nekData.ifgetps = (int*) nek_ptr("ifgetps");

  nekData.unx = (double*) nek_ptr("unx");
  nekData.uny = (double*) nek_ptr("uny");
  nekData.unz = (double*) nek_ptr("unz");

  nekData.xm1 = (double*) nek_ptr("xm1");
  nekData.ym1 = (double*) nek_ptr("ym1");
  nekData.zm1 = (double*) nek_ptr("zm1");
  nekData.xc = (double*) nek_ptr("xc");
  nekData.yc = (double*) nek_ptr("yc");
  nekData.zc = (double*) nek_ptr("zc");

  nekData.glo_num = (long long*) nek_ptr("glo_num");
  nekData.cbscnrs = (double*) nek_ptr("cb_scnrs");
  nekData.cbc = (char*) nek_ptr("cbc");

  nekData.boundaryID  = (int*) nek_ptr("boundaryID");
  nekData.boundaryIDt = (int*) nek_ptr("boundaryIDt");

  nekData.eface1 = (int*) nek_ptr("eface1");
  nekData.eface = (int*) nek_ptr("eface");
  nekData.icface = (int*) nek_ptr("icface");
  nekData.comm = MPI_Comm_f2c(*(int*) nek_ptr("nekcomm"));

  int cht = 0;
  if (nekData.nelv != nekData.nelt && nscal) cht = 1;

  // import BCs from nek if not specified in par
  if(!bcInPar) {
    nek_gen_bcmap();
    if(flow) {
      int isTMesh = 0;
      int nIDs = (*nek_nbid_ptr)(&isTMesh);
      int* map = (int*) calloc(nIDs, sizeof(int));
      for(int id = 0; id < nIDs; id++) map[id] = nek_bcmap(id + 1, 1);
      bcMap::setBcMap("velocity", map, nIDs);
      free(map);
    }
    for(int is = 0; is < nscal; is++) {
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(2) << is;
      string sid = ss.str();

      int isTMesh = 0;
      if (cht && is == 0) isTMesh = 1;
      int nIDs = (*nek_nbid_ptr)(&isTMesh);

      int* map = (int*) calloc(nIDs, sizeof(int));
      for(int id = 0; id < nIDs; id++) map[id] = nek_bcmap(id + 1, is + 2);
      bcMap::setBcMap("scalar" + sid, map, nIDs);
      free(map);
    }
  }

  int isTMesh = 0;
  nekData.NboundaryID  = (*nek_nbid_ptr)(&isTMesh);
  isTMesh = 1;
  nekData.NboundaryIDt = (*nek_nbid_ptr)(&isTMesh);

  dfloat rho;
  options->getArgs("DENSITY", rho);
  nekData.param[0] = rho;

  dfloat mue;
  options->getArgs("VISCOSITY", mue);
  nekData.param[1] = mue;

  dfloat rhoCp;
  options->getArgs("SCALAR00 DENSITY", rhoCp);
  nekData.param[6] = rhoCp;

  dfloat diff;
  options->getArgs("SCALAR00 DIFFUSIVITY", diff);
  nekData.param[7] = diff;

  return 0;
}

void nek_copyFrom(dfloat time)
{
  if(rank == 0) {
    printf("copying solution to nek\n");
    fflush(stdout);
  }

  timeLast = time;

  mesh_t* mesh = (*ins)->mesh;
  dlong Nlocal = mesh->Nelements * mesh->Np;

  dfloat* vx = (*ins)->U + 0 * (*ins)->fieldOffset;
  dfloat* vy = (*ins)->U + 1 * (*ins)->fieldOffset;
  dfloat* vz = (*ins)->U + 2 * (*ins)->fieldOffset;

  *(nekData.time) = time;

  memcpy(nekData.vx, vx, sizeof(dfloat) * Nlocal);
  memcpy(nekData.vy, vy, sizeof(dfloat) * Nlocal);
  memcpy(nekData.vz, vz, sizeof(dfloat) * Nlocal);
  memcpy(nekData.pr, (*ins)->P, sizeof(dfloat) * Nlocal);
  if((*ins)->Nscalar) {
    const dlong nekFieldOffset = nekData.lelt * mesh->Np;
    for(int is = 0; is < (*ins)->Nscalar; is++) {
      mesh_t* mesh;
      (is) ? mesh = (*ins)->cds->meshV : mesh = (*ins)->cds->mesh;
      const dlong Nlocal = mesh->Nelements * mesh->Np;
      dfloat* Ti = nekData.t   + is * nekFieldOffset;
      dfloat* Si = (*ins)->cds->S + is * (*ins)->cds->fieldOffset;
      memcpy(Ti, Si, Nlocal * sizeof(dfloat));
    }
  }
}

void nek_ocopyFrom(dfloat time, int tstep)
{
  if(time != timeLast) {
    (*ins)->o_U.copyTo((*ins)->U);
    (*ins)->o_P.copyTo((*ins)->P);
    if((*ins)->Nscalar) (*ins)->cds->o_S.copyTo((*ins)->cds->S);
    nek_copyFrom(time, tstep);
  }
}

void nek_copyFrom(dfloat time, int tstep)
{
  nek_copyFrom(time);
  *(nekData.istep) = tstep;
}

void nek_ocopyTo(dfloat &time)
{
  nek_copyTo(time);
  (*ins)->o_P.copyFrom((*ins)->P);
  (*ins)->o_U.copyFrom((*ins)->U);
  if((*ins)->Nscalar) (*ins)->cds->o_S.copyFrom((*ins)->cds->S);
}

void nek_copyTo(dfloat &time)
{
  if(rank == 0) {
    printf("copying solution from nek\n");
    fflush(stdout);
  }

  mesh_t* mesh = (*ins)->mesh;
  dlong Nlocal = mesh->Nelements * mesh->Np;

  time = *(nekData.time);

  dfloat* vx = (*ins)->U + 0 * (*ins)->fieldOffset;
  dfloat* vy = (*ins)->U + 1 * (*ins)->fieldOffset;
  dfloat* vz = (*ins)->U + 2 * (*ins)->fieldOffset;

  memcpy(vx, nekData.vx, sizeof(dfloat) * Nlocal);
  memcpy(vy, nekData.vy, sizeof(dfloat) * Nlocal);
  memcpy(vz, nekData.vz, sizeof(dfloat) * Nlocal);
  memcpy((*ins)->P, nekData.pr, sizeof(dfloat) * Nlocal);
  if((*ins)->Nscalar) {
    const dlong nekFieldOffset = nekData.lelt * mesh->Np;
    for(int is = 0; is < (*ins)->Nscalar; is++) {
      mesh_t* mesh;
      (is) ? mesh = (*ins)->cds->meshV : mesh = (*ins)->cds->mesh;
      const dlong Nlocal = mesh->Nelements * mesh->Np;
      dfloat* Ti = nekData.t   + is * nekFieldOffset;
      dfloat* Si = (*ins)->cds->S + is * (*ins)->cds->fieldOffset;
      memcpy(Si, Ti, Nlocal * sizeof(dfloat));
    }
  }
}

void nek_copyRestart()
{
  mesh_t* mesh = (*ins)->mesh;
  dlong Nlocal = mesh->Nelements * mesh->Np;
  if (*(nekData.ifgetu)) {
    dfloat* vx = (*ins)->U + 0 * (*ins)->fieldOffset;
    dfloat* vy = (*ins)->U + 1 * (*ins)->fieldOffset;
    dfloat* vz = (*ins)->U + 2 * (*ins)->fieldOffset;
    memcpy(vx, nekData.vx, sizeof(dfloat) * Nlocal);
    memcpy(vy, nekData.vy, sizeof(dfloat) * Nlocal);
    memcpy(vz, nekData.vz, sizeof(dfloat) * Nlocal);
  }
  if (*(nekData.ifgetp)) memcpy((*ins)->P, nekData.pr, sizeof(dfloat) * Nlocal);
  if((*ins)->Nscalar) {
    const dlong nekFieldOffset = nekData.lelt * mesh->Np;
    for(int is = 0; is < (*ins)->Nscalar; is++) {
      mesh_t* mesh;
      (is) ? mesh = (*ins)->cds->meshV : mesh = (*ins)->cds->mesh;
      const dlong Nlocal = mesh->Nelements * mesh->Np;
      dfloat* Ti = nekData.t   + is * nekFieldOffset;
      dfloat* Si = (*ins)->cds->S + is * (*ins)->cds->fieldOffset;
      memcpy(Si, Ti, Nlocal * sizeof(dfloat));
    }
  }
}

long long nek_set_glo_num(int nx, int isTMesh)
{
  return (*nek_set_vert_ptr)(&nx, &isTMesh);
}
