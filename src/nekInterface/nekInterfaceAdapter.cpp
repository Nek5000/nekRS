#include <unistd.h>
#include <fstream>
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "bcMap.hpp"
#include "fileUtils.hpp"
#include "re2Reader.hpp"

nekdata_private nekData;
static int rank;
static setupAide* options;
static nrs_t *nrs;

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
static void (* nek_scptr_ptr)(int*, void*);
static void (* nek_outfld_ptr)(char*, int);
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
static void (* nek_bootstrap_ptr)(int*, char*, char*, char*, int, int, int);
static void (
    *nek_setup_ptr)(int *, int *, int *, int *, double *, double *, double *, double *, double *, int *);
static void (* nek_ifoutfld_ptr)(int*);
static void (* nek_setics_ptr)(void);
static int (* nek_bcmap_ptr)(int*, int*,int*);
static void (* nek_gen_bcmap_ptr)(void);
static int (* nek_nbid_ptr)(int*);
static long long (* nek_set_vert_ptr)(int*, int*);

static void (* nek_setbd_ptr)(double *, double*, int*);
static void (* nek_setabbd_ptr)(double *, double*, int*, int*);

static void (* nek_storesol_ptr)(void);
static void (* nek_restoresol_ptr)(void);
static void (* nek_updggeom_ptr)(void);
static void (* nek_meshmetrics_ptr)(void);

void noop_func(void) {}

void check_error(const char* error)
{
  if(error != NULL) {
    if(rank == 0) fprintf(stderr, "Error: %s!\n", error);
    ABORT(EXIT_FAILURE);
  }
}

namespace nek{
void* ptr(const char* id)
{
  void* ptr;
  int len = strlen(id);
  (*nek_ptr_ptr)(&ptr, (char*)id, &len);
  return ptr;
}

void* scPtr(int id)
{
  if(!id) {
    if(rank == 0) fprintf(stderr, "Error: nek::scPtr index <1!\n");
    ABORT(EXIT_FAILURE);
  }

  void* ptr;
  (*nek_scptr_ptr)(&id, &ptr);
  return ptr;
}

void outfld(const char *filename, dfloat t, int step, int coords, int FP64,
            void *o_uu, void *o_pp, void *o_ss,
            int NSfields)
{

  mesh_t* mesh = nrs->meshV;
  dlong Nlocal = mesh->Nelements * mesh->Np;

  double time = t;

  if(NSfields > nekData.ldimt) {
    const char *errTxt = "NSfields > ldimt in nek_outfld";
    check_error(errTxt);
  }

  occa::memory o_u, o_p, o_s;
  if(o_uu) o_u = *((occa::memory *) o_uu);
  if(o_pp) o_p = *((occa::memory *) o_pp);
  if(o_ss && NSfields) o_s = *((occa::memory *) o_ss);

  int xo = 0;
  int vo = 0;
  int po = 0;
  int so = 0;

  (*nek_storesol_ptr)();

  const int stepSave = *(nekData.istep);
  *(nekData.istep) = step;

  platform->timer.tic("checkpointing", 1);

  if(coords){
    mesh_t *mesh = nrs->meshV;
    if(nrs->cht) mesh = nrs->cds->mesh[0];	  
    mesh->o_x.copyTo(nekData.xm1, Nlocal * sizeof(dfloat));
    mesh->o_y.copyTo(nekData.ym1, Nlocal * sizeof(dfloat));
    mesh->o_z.copyTo(nekData.zm1, Nlocal * sizeof(dfloat));
    xo = 1;
  }
  if(o_u.ptr()) {
    occa::memory o_vx = o_u + (0 * sizeof(dfloat)) * nrs->fieldOffset;
    occa::memory o_vy = o_u + (1 * sizeof(dfloat)) * nrs->fieldOffset;
    occa::memory o_vz = o_u + (2 * sizeof(dfloat)) * nrs->fieldOffset;
    o_vx.copyTo(nekData.vx, Nlocal * sizeof(dfloat));
    o_vy.copyTo(nekData.vy, Nlocal * sizeof(dfloat));
    o_vz.copyTo(nekData.vz, Nlocal * sizeof(dfloat));
    vo = 1;
  }
  if(o_p.ptr()) {
    o_p.copyTo(nekData.pr, Nlocal * sizeof(dfloat));
    po = 1;
  }
  if(o_s.ptr()) {
    const dlong nekFieldOffset = nekData.lelt * mesh->Np;
    for(int is = 0; is < NSfields; is++) {
      mesh_t* mesh = nrs->meshV;
      if(nrs->cds)
        (is) ? mesh = nrs->meshV: mesh = nrs->cds->mesh[0];
      const dlong Nlocal = mesh->Nelements * mesh->Np;
      dfloat* Ti = nekData.t + is * nekFieldOffset;
      occa::memory o_Si = o_s + (is * sizeof(dfloat)) * nrs->fieldOffset;
      o_Si.copyTo(Ti, Nlocal * sizeof(dfloat));
    }
    so = 1;
  }

  (*nek_setio_ptr)(&t, &xo, &vo, &po, &so, &NSfields, &FP64);
  (*nek_outfld_ptr)((char*) filename, strlen(filename));
  (*nek_resetio_ptr)();

  platform->timer.toc("checkpointing");

  (*nek_restoresol_ptr)();

  *(nekData.istep) = stepSave;
}

void uic(int ifield)
{
  (*nek_uic_ptr)(&ifield);
}

void end()
{
  (*nek_end_ptr)();
}

void setic(void)
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
}

void map_m_to_n(double* a, int na, double* b, int nb)
{
  // output a, intput b
  int if3d = (nekData.ndim == 3);

  int maxn = na > nb ? na : nb;
  int N = 2 * maxn * maxn * maxn;
  double* w = (double*) malloc(sizeof(double) * N);

  (*nek_map_m_to_n_ptr)(a, &na, b, &nb, &if3d, w, &N);
  free(w);
}

void uf(double* u, double* v, double* w)
{
  (*nek_uf_ptr)(u, v, w);
}

int lglel(int e)
{
  int ee = e + 1;
  return (*nek_lglel_ptr)(&ee) - 1;
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
  if(rank == 0) printf("calling nek_userchk ...\n");
  (*userchk_ptr)();
}
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

void set_usr_handles(const char* session_in,int verbose)
{
  // load lib{session_in}.so
  char lib_session[BUFSIZ], * error;

  const char* cache_dir = getenv("NEKRS_CACHE_DIR");
  sprintf(lib_session, "%s/nek5000/lib%s.so", cache_dir, session_in);

  void* handle = dlopen(lib_session,RTLD_NOW | RTLD_LOCAL);
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
  nek_scptr_ptr = (void (*)(int*, void*))dlsym(handle, fname("nekf_scptr"));
  check_error(dlerror());
  nek_bootstrap_ptr =
    (void (*)(int*, char*, char*, char*, int, int, int))dlsym(handle, fname("nekf_bootstrap"));
  check_error(dlerror());
  nek_setup_ptr =
      (void (*)(int *, int *, int *, int *, double *, double *, double *, double *, double *, int *))
          dlsym(handle, fname("nekf_setup"));
  check_error(dlerror());
  nek_uic_ptr = (void (*)(int*))dlsym(handle, fname("nekf_uic"));
  check_error(dlerror());
  nek_end_ptr = (void (*)(void))dlsym(handle, fname("nekf_end"));
  check_error(dlerror());
  nek_outfld_ptr = (void (*)(char*, int))dlsym(handle, fname("nekf_outfld"));
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
  nek_bcmap_ptr = (int (*)(int*, int*,int*))dlsym(handle,fname("nekf_bcmap"));
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

  nek_setbd_ptr = (void (*)(double *, double*, int*))dlsym(handle, fname("setbd"));
  check_error(dlerror());
  nek_setabbd_ptr = (void (*)(double *, double*, int*, int*))dlsym(handle, fname("setabbd"));
  check_error(dlerror());

  nek_storesol_ptr = (void (*)(void))dlsym(handle, fname("nekf_storesol"));
  check_error(dlerror());
  nek_restoresol_ptr = (void (*)(void))dlsym(handle, fname("nekf_restoresol"));
  check_error(dlerror());
  nek_updggeom_ptr = (void (*)(void))dlsym(handle, fname("nekf_updggeom"));
  check_error(dlerror());
  nek_meshmetrics_ptr = (void (*)(void))dlsym(handle, fname("mesh_metrics"));
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

void mkSIZE(int lx1, int lxd, int lelt, int lelg, int ldim, int lpmin, int ldimt, setupAide& options, char* SIZE)
{
  char line[BUFSIZ];
  const char *cache_dir = getenv("NEKRS_CACHE_DIR");

  const std::string installDir(getenv("NEKRS_HOME"));
  const std::string nek5000_dir = installDir + "/nek5000";

  const int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;

  // Read and generate the new size file.
  sprintf(line,"%s/core/SIZE.template", nek5000_dir.c_str());
  FILE* fp = fopen(line, "r");
  if (!fp) {
    fprintf(stderr, "Error opening %s!\n", line);
    ABORT(EXIT_FAILURE);
  }

  char *sizeFile;
  {
    fseek(fp, 0, SEEK_END);
    const long length = ftell(fp);
    rewind(fp);
    sizeFile = (char*) calloc(length + 500,sizeof(char));
    if(!sizeFile) {
      fprintf(stderr, "Error allocating space for SIZE file.\n");
      ABORT(EXIT_FAILURE);
    }
  }

  int lx1m = lx1;

  constexpr int nMaxObj = 20;

  int count = 0;
  while(fgets(line, BUFSIZ, fp) != NULL) {
    if(strstr(line, "parameter (lx1=") != NULL)
      sprintf(line, "      parameter (lx1=%d)\n", lx1);
    else if(strstr(line, "parameter (lxd=") != NULL)
      sprintf(line, "      parameter (lxd=%d)\n", lxd);
    else if(strstr(line, "parameter (lelt=") != NULL)
      sprintf(line, "      parameter (lelt=%d)\n", lelt);
    else if(strstr(line, "parameter (lelg=") != NULL)
      sprintf(line, "      parameter (lelg=%d)\n", lelg);
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
    else if(strstr(line, "parameter (lhis=") != NULL)
      sprintf(line, "      parameter (lhis=%d)\n", 100000);
    else if(strstr(line, "parameter (lelr=") != NULL)
      sprintf(line, "      parameter (lelr=%d)\n", 128 * lelt);
    else if(strstr(line, "parameter (lx1m=") != NULL)
      sprintf(line, "      parameter (lx1m=%d)\n", lx1m);
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
  if(osize.is_open()) {
    writeSize = 0;
    std::string line;
    while(getline( osize, line )) {
      if(line.find( "lelg=") != std::string::npos ) {
	hlong oldval;      
        sscanf(line.c_str(), "%*[^=]=%lld", &oldval);
        if(oldval < lelg) writeSize = 1;
      }
      if(line.find( "lelt=") != std::string::npos ) {
	int oldval;      
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if(oldval < lelt) writeSize = 1;
      }
      if(line.find( "lx1m=") != std::string::npos ) {
	int oldval;      
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if(oldval < lx1m) writeSize = 1;
      }
      if(line.find( "lx1=") != std::string::npos ) {
	int oldval;      
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if(oldval != lx1) writeSize = 1;
      }
      if(line.find( "ldimt=") != std::string::npos ) {
	int oldval;      
        sscanf(line.c_str(), "%*[^=]=%d", &oldval);
        if(oldval < ldimt) writeSize = 1;
      }
    }
  }
  osize.close();

  if(writeSize) {
    fp = fopen(SIZE, "w");
    fputs(sizeFile, fp);
    fclose(fp);
    free(sizeFile);
    if(verbose) printf("using new SIZE\n");
  } else {
    if(verbose) printf("using existing SIZE file %s/SIZE\n", cache_dir);
  }

  fflush(stdout);
}


void buildNekInterface(int ldimt, int N, int np, setupAide& options)
{
  int buildRank = rank;
  const bool buildNodeLocal = useNodeLocalCache();
  if(buildNodeLocal)
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);    
  
  int err = [&](){
  
    if(buildRank == 0){
      const std::string cache_dir = std::string(getenv("NEKRS_CACHE_DIR")) + "/nek5000";
      mkdir(cache_dir.c_str(), S_IRWXU); 

      const std::string installDir(getenv("NEKRS_HOME"));
      const std::string nekInterface_dir = installDir + "/nekInterface";
      const std::string nek5000_dir = installDir + "/nek5000";

      char buf[10*BUFSIZ];
      char *ret = getcwd(buf, sizeof(buf));
      const std::string case_dir(buf);

      const int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;

      const std::string usrname = options.getArgs("CASENAME");
      const std::string meshFile = options.getArgs("MESH FILE");

      int nelgt, nelgv;
      const int ndim = 3;
      re2::nelg(meshFile, nelgt, nelgv, MPI_COMM_NULL); 

      int lelt = (nelgt/np) + 3;
      if(lelt > nelgt) lelt = nelgt;
      sprintf(buf,"%s/SIZE",cache_dir.c_str());
      mkSIZE(N + 1, 1, lelt, nelgt, ndim, np, ldimt, options, buf);

      // generate usr
      char usrFileCache[BUFSIZ];
      {
        char usrFile[BUFSIZ];
        const std::string usrFileStr = options.getArgs("NEK USR FILE");
        char usrFileCaseName[BUFSIZ];
        sprintf(usrFileCaseName,"%s.usr",usrname.c_str());
        strcpy(usrFile, usrFileStr.c_str());
        if(!fileExists(usrFile) || isFileEmpty(usrFile))
          sprintf(usrFile, "%s/core/zero.usr", nek5000_dir.c_str());

        sprintf(usrFileCache,"%s/%s",cache_dir.c_str(),usrFileCaseName);
        if(isFileNewer(usrFile, usrFileCache))
          copyFile(usrFile, usrFileCache);
      }

      // build
      char libFile[BUFSIZ];
      sprintf(libFile,"%s/lib%s.so",cache_dir.c_str(),usrname.c_str());
      int recompile = 0;
      if(isFileNewer(usrFileCache, libFile)) recompile = 1;  
      sprintf(buf,"%s/SIZE",cache_dir.c_str());
      if(isFileNewer(buf, libFile)) recompile = 1;  
      if(recompile) {
        const double tStart = MPI_Wtime();
        const std::string pipeToNull = (rank == 0) ? std::string("") :  std::string(">/dev/null 2>&1");
	    const std::string include_dirs = "./ " + case_dir;
        std::string make_args = "-j8 "; 
        if(!verbose) make_args += "-s "; 
        if(rank == 0) 
	      printf("building nekInterface for lx1=%d, lelt=%d and lelg=%d ...", N+1, lelt, nelgt); fflush(stdout);

        sprintf(buf, 
             "cd %s"
             " && cp -f %s/makefile.template makefile"
             "  && make %s"
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
		     usrname.c_str(), 
		     cache_dir.c_str(), 
             nekInterface_dir.c_str(), 
		     pipeToNull.c_str());
        if(verbose && rank == 0) printf("\n%s\n", buf);
        if(system(buf)) return EXIT_FAILURE;
        fileSync(libFile);

        if(rank == 0) printf("done (%gs)\n", MPI_Wtime() - tStart);
        fflush(stdout);
      } else {
        if(rank == 0) printf("skip building nekInterface (SIZE requires no update)\n");
        fflush(stdout);
      }
    } 
    return 0;
  }();

  MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, platform->comm.mpiComm);
  if(err) ABORT(EXIT_FAILURE);
}

namespace nek{
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
  MPI_Comm_rank(platform->comm.mpiComm,&rank);
  MPI_Comm_size(platform->comm.mpiComm,&size);

  int buildRank = rank;
  const bool buildNodeLocal = useNodeLocalCache();
  if(buildNodeLocal)
    MPI_Comm_rank(platform->comm.mpiCommLocal, &buildRank);    

  int N;
  options->getArgs("POLYNOMIAL DEGREE", N);
 
  int Nscalar;
  options->getArgs("NUMBER OF SCALARS", Nscalar);

  int npTarget = size;
  options->getArgs("NP TARGET", npTarget);

  buildNekInterface(mymax(5, Nscalar), N, npTarget, *options);

  if (platform->options.compareArgs("BUILD ONLY", "FALSE")) {
    if (rank == 0) { 
     printf("loading nek ... \n"); 
     fflush(stdout);
    }

    std::string usrname;
    options->getArgs("CASENAME", usrname);
    std::string meshFile;
    options->getArgs("MESH FILE", meshFile);

    char buf[FILENAME_MAX];
    getcwd(buf, sizeof(buf));
    std::string cwd;
    cwd.assign(buf);

    MPI_Fint nek_comm = MPI_Comm_c2f(platform->comm.mpiComm);

    set_usr_handles(usrname.c_str(), 0);

    (*nek_bootstrap_ptr)(&nek_comm, (char*)cwd.c_str(), 
                         /* basename */ (char*)usrname.c_str(),
                         (char*)meshFile.c_str(),
                         cwd.length(), usrname.length(),
                         meshFile.length());
    if (rank == 0) { 
     printf("done\n"); 
     fflush(stdout);
    }


  }
}

int setup(nrs_t* nrs_in)
{
  nrs = nrs_in;
  MPI_Comm_rank(platform->comm.mpiComm, &rank);

  bool meshSolver = options->compareArgs("MESH SOLVER", "ELASTICITY");

  std::string casename;
  options->getArgs("CASENAME", casename);

  int nscal = 0;
  options->getArgs("NUMBER OF SCALARS", nscal);

  std::string velocitySolver;
  int flow = 1;
  if(options->compareArgs("VELOCITY", "FALSE")) flow = 0;

  int meshPartType = 3; // RCB+RSB
  if(options->compareArgs("MESH PARTITIONER", "rcb")) meshPartType = 2;
  if(options->compareArgs("MESH PARTITIONER", "rcb+rsb")) meshPartType = 3;

  double meshConTol = 0.2;
  options->getArgs("MESH CONNECTIVITY TOL", meshConTol);

  int nBcRead = 1;
  if (bcMap::useNekBCs())
    nBcRead = flow + nscal;

  dfloat rho;
  options->getArgs("DENSITY", rho);

  dfloat mue;
  options->getArgs("VISCOSITY", mue);

  dfloat rhoCp;
  options->getArgs("SCALAR00 DENSITY", rhoCp);

  dfloat lambda;
  options->getArgs("SCALAR00 DIFFUSIVITY", lambda);

  int stressForm = 1; // avoid recompilation + bypass unligned SYM/SHL check 

  (*nek_setup_ptr)(&flow,
                   &nscal,
                   &nBcRead,
                   &meshPartType,
                   &meshConTol,
                   &rho,
                   &mue,
                   &rhoCp,
                   &lambda,
                   &stressForm);

  nekData.param = (double*) ptr("param");
  nekData.ifield = (int*) ptr("ifield");
  nekData.istep = (int*) ptr("istep");
  nekData.time = (double*) ptr("time");

  nekData.ndim = *(int*) ptr("ndim");
  nekData.nelt = *(int*) ptr("nelt");
  nekData.nelv = *(int*) ptr("nelv");
  nekData.lelt = *(int*) ptr("lelt");
  nekData.ldimt = *(int*) ptr("ldimt");
  nekData.nx1 =  *(int*) ptr("nx1");

  nekData.vx = (double*) ptr("vx");
  nekData.vy = (double*) ptr("vy");
  nekData.vz = (double*) ptr("vz");
  nekData.pr = (double*) ptr("pr");
  nekData.t  = (double*) ptr("t");

  nekData.qtl = (double*) ptr("qtl");

  nekData.ifgetu  = (int*) ptr("ifgetu");
  nekData.ifgetp  = (int*) ptr("ifgetp");
  nekData.ifgett  = (int*) ptr("ifgett");
  nekData.ifgetps = (int*) ptr("ifgetps");

  nekData.unx = (double*) ptr("unx");
  nekData.uny = (double*) ptr("uny");
  nekData.unz = (double*) ptr("unz");

  nekData.xm1 = (double*) ptr("xm1");
  nekData.ym1 = (double*) ptr("ym1");
  nekData.zm1 = (double*) ptr("zm1");
  nekData.xc = (double*) ptr("xc");
  nekData.yc = (double*) ptr("yc");
  nekData.zc = (double*) ptr("zc");

  nekData.glo_num = (long long*) ptr("glo_num");
  nekData.cbc = (char*) ptr("cbc");

  nekData.boundaryID  = (int*) ptr("boundaryID");
  nekData.boundaryIDt = (int*) ptr("boundaryIDt");

  nekData.eface1 = (int*) ptr("eface1");
  nekData.eface = (int*) ptr("eface");
  nekData.icface = (int*) ptr("icface");
  nekData.comm = MPI_Comm_f2c(*(int*) ptr("nekcomm"));
  
  nekData.p0th = (double*) ptr("p0th");
  nekData.dp0thdt = (double*) ptr("dp0thdt");

  nekData.wx = (double*) ptr("wx");
  nekData.wy = (double*) ptr("wy");
  nekData.wz = (double*) ptr("wz");

  int cht = 0;
  if (nekData.nelv != nekData.nelt && nscal) cht = 1;

  if (bcMap::useNekBCs()) {
    if (rank == 0)
      printf("import BCs from nek\n");
    gen_bcmap();
    if(flow) {
      int isTMesh = 0;
      int nIDs = (*nek_nbid_ptr)(&isTMesh);
      int* map = (int*) calloc(nIDs, sizeof(int));
      for(int id = 0; id < nIDs; id++) map[id] = bcmap(id + 1, 1, 0);
      bcMap::setBcMap("velocity", map, nIDs);

      if(meshSolver){
        for(int id = 0; id < nIDs; id++) map[id] = bcmap(id + 1, 1, 1);
        bcMap::setBcMap("mesh", map, nIDs);
      }

      free(map);
    }
    for(int is = 0; is < nscal; is++) {
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(getDigitsRepresentation(NSCALAR_MAX - 1)) << is;
      std::string sid = ss.str();

      int isTMesh = 0;
      if (cht && is == 0) isTMesh = 1;
      int nIDs = (*nek_nbid_ptr)(&isTMesh);

      int* map = (int*) calloc(nIDs, sizeof(int));
      for(int id = 0; id < nIDs; id++) map[id] = bcmap(id + 1, is + 2, 0);
      bcMap::setBcMap("scalar" + sid, map, nIDs);
      free(map);
    }
  }

  int isTMesh = 0;
  nekData.NboundaryID  = (*nek_nbid_ptr)(&isTMesh);
  isTMesh = 1;
  nekData.NboundaryIDt = (*nek_nbid_ptr)(&isTMesh);

  dfloat startTime;
  options->getArgs("START TIME", startTime);
  *(nekData.time) = startTime; 

  return 0;
}

void copyToNek(dfloat time)
{
  if(rank == 0) {
    printf("copying solution to nek\n");
    fflush(stdout);
  }

  mesh_t* mesh = nrs->meshV;
  dlong Nlocal = mesh->Nelements * mesh->Np;

  dfloat* vx = nrs->U + 0 * nrs->fieldOffset;
  dfloat* vy = nrs->U + 1 * nrs->fieldOffset;
  dfloat* vz = nrs->U + 2 * nrs->fieldOffset;

  *(nekData.time) = time;
  *(nekData.p0th) = nrs->p0th[0];
  *(nekData.dp0thdt) = nrs->dp0thdt;

  if(platform->options.compareArgs("MOVING MESH", "TRUE")){
    mesh_t *mesh = nrs->meshV;
    if(nrs->cht) mesh = nrs->cds->mesh[0];
    const dlong Nlocal = mesh->Nelements * mesh->Np;
    dfloat* wx = mesh->U + 0 * nrs->fieldOffset;
    dfloat* wy = mesh->U + 1 * nrs->fieldOffset;
    dfloat* wz = mesh->U + 2 * nrs->fieldOffset;
    memcpy(nekData.wx, wx, sizeof(dfloat) * Nlocal);
    memcpy(nekData.wy, wy, sizeof(dfloat) * Nlocal);
    memcpy(nekData.wz, wz, sizeof(dfloat) * Nlocal);
    memcpy(nekData.xm1, mesh->x, sizeof(dfloat) * Nlocal);
    memcpy(nekData.ym1, mesh->y, sizeof(dfloat) * Nlocal);
    memcpy(nekData.zm1, mesh->z, sizeof(dfloat) * Nlocal);
    recomputeGeometry();
  }

  memcpy(nekData.vx, vx, sizeof(dfloat) * Nlocal);
  memcpy(nekData.vy, vy, sizeof(dfloat) * Nlocal);
  memcpy(nekData.vz, vz, sizeof(dfloat) * Nlocal);
  memcpy(nekData.pr, nrs->P, sizeof(dfloat) * Nlocal);
  if(nrs->Nscalar) {
    if(platform->options.compareArgs("LOWMACH", "TRUE")) memcpy(nekData.qtl, nrs->div, sizeof(dfloat) * Nlocal);
    const dlong nekFieldOffset = nekData.lelt * mesh->Np;
    for(int is = 0; is < nrs->Nscalar; is++) {
      mesh_t* mesh;
      (is) ? mesh = nrs->cds->meshV : mesh = nrs->cds->mesh[0];
      const dlong Nlocal = mesh->Nelements * mesh->Np;
      dfloat* Ti = nekData.t   + is * nekFieldOffset;
      dfloat* Si = nrs->cds->S + nrs->cds->fieldOffsetScan[is];
      memcpy(Ti, Si, Nlocal * sizeof(dfloat));
    }
  }
}

void ocopyToNek(void)
{
  nrs->o_U.copyTo(nrs->U);
  nrs->o_P.copyTo(nrs->P);
  if(nrs->Nscalar){
    nrs->o_div.copyTo(nrs->div);
    nrs->cds->o_S.copyTo(nrs->cds->S);
  }
  if(platform->options.compareArgs("MOVING MESH", "TRUE")){
    mesh_t *mesh = nrs->meshV;
    if(nrs->cht) mesh = nrs->cds->mesh[0];
    mesh->o_U.copyTo(mesh->U);
    mesh->o_x.copyTo(mesh->x);
    mesh->o_y.copyTo(mesh->y);
    mesh->o_z.copyTo(mesh->z);
  }
  copyToNek(0.0);
}

void ocopyToNek(dfloat time, int tstep)
{
  nrs->o_U.copyTo(nrs->U);
  nrs->o_P.copyTo(nrs->P);
  if(nrs->Nscalar){
    nrs->o_div.copyTo(nrs->div);
    nrs->cds->o_S.copyTo(nrs->cds->S);
  }
  if(platform->options.compareArgs("MOVING MESH", "TRUE")){
    mesh_t *mesh = nrs->meshV;
    if(nrs->cht) mesh = nrs->cds->mesh[0];
    mesh->o_U.copyTo(mesh->U);
    mesh->o_x.copyTo(mesh->x);
    mesh->o_y.copyTo(mesh->y);
    mesh->o_z.copyTo(mesh->z);
  }
  copyToNek(time, tstep);
}

void copyToNek(dfloat time, int tstep)
{
  *(nekData.istep) = tstep;
  copyToNek(time);
}

void ocopyFromNek(dfloat &time)
{
  copyFromNek(time);
  nrs->o_P.copyFrom(nrs->P);
  nrs->o_U.copyFrom(nrs->U);
  if(nrs->Nscalar){
    nrs->cds->o_S.copyFrom(nrs->cds->S);
  }
  if(platform->options.compareArgs("MOVING MESH", "TRUE")){
    mesh_t *mesh = nrs->meshV;
    if(nrs->cht) mesh = nrs->cds->mesh[0];
    mesh->o_x.copyFrom(mesh->x);
    mesh->o_y.copyFrom(mesh->y);
    mesh->o_z.copyFrom(mesh->z);
    mesh->o_U.copyFrom(mesh->U);
  }
}

void copyFromNek(dfloat &time)
{
  if(rank == 0) {
    printf("copying solution from nek\n");
    fflush(stdout);
  }

  mesh_t* mesh = nrs->meshV;
  dlong Nlocal = mesh->Nelements * mesh->Np;

  time = *(nekData.time);
  nrs->p0th[0] = *(nekData.p0th);
  nrs->dp0thdt = *(nekData.dp0thdt);

  dfloat* vx = nrs->U + 0 * nrs->fieldOffset;
  dfloat* vy = nrs->U + 1 * nrs->fieldOffset;
  dfloat* vz = nrs->U + 2 * nrs->fieldOffset;

  memcpy(vx, nekData.vx, sizeof(dfloat) * Nlocal);
  memcpy(vy, nekData.vy, sizeof(dfloat) * Nlocal);
  memcpy(vz, nekData.vz, sizeof(dfloat) * Nlocal);
  if(platform->options.compareArgs("MOVING MESH", "TRUE")){
    mesh_t *mesh = nrs->meshV;
    if(nrs->cht) mesh = nrs->cds->mesh[0];
    const dlong Nlocal = mesh->Nelements * mesh->Np;
    dfloat* wx = mesh->U + 0 * nrs->fieldOffset;
    dfloat* wy = mesh->U + 1 * nrs->fieldOffset;
    dfloat* wz = mesh->U + 2 * nrs->fieldOffset;
    memcpy(wx, nekData.wx, sizeof(dfloat) * Nlocal);
    memcpy(wy, nekData.wy, sizeof(dfloat) * Nlocal);
    memcpy(wz, nekData.wz, sizeof(dfloat) * Nlocal);
    memcpy(nekData.xm1, mesh->x, sizeof(dfloat) * Nlocal);
    memcpy(nekData.ym1, mesh->y, sizeof(dfloat) * Nlocal);
    memcpy(nekData.zm1, mesh->z, sizeof(dfloat) * Nlocal);
    recomputeGeometry();
  }
  memcpy(nrs->P, nekData.pr, sizeof(dfloat) * Nlocal);
  if(nrs->Nscalar) {
    const dlong nekFieldOffset = nekData.lelt * mesh->Np;
    for(int is = 0; is < nrs->Nscalar; is++) {
      mesh_t* mesh;
      (is) ? mesh = nrs->cds->meshV : mesh = nrs->cds->mesh[0];
      const dlong Nlocal = mesh->Nelements * mesh->Np;
      dfloat* Ti = nekData.t   + is * nekFieldOffset;
      dfloat* Si = nrs->cds->S + nrs->cds->fieldOffsetScan[is];
      memcpy(Si, Ti, Nlocal * sizeof(dfloat));
    }
  }
}

long long set_glo_num(int nx, int isTMesh)
{
  return (*nek_set_vert_ptr)(&nx, &isTMesh);
}

void bdfCoeff(double *g0, double *coeff, double *dt, int order) 
{
  double nekCoeff[4];
  (*nek_setbd_ptr)(nekCoeff, dt, &order);
  *g0 = nekCoeff[0];
  memcpy(coeff, &nekCoeff[1], order*sizeof(double));
}

void extCoeff(double *coeff, double *dt, int nAB, int nBDF)
{
  (*nek_setabbd_ptr)(coeff, dt, &nAB, &nBDF);
}
void coeffAB(double *coeff, double *dt, int order)
{
  int one = 1;
  (*nek_setabbd_ptr)(coeff, dt, &order, &one);
}

void recomputeGeometry() { (*nek_updggeom_ptr)(); }

void printMeshMetrics() { (*nek_meshmetrics_ptr)(); }

}
