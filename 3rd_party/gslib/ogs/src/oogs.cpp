#include <limits>
#include <list>
#include <occa.hpp>

#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "ogsInterface.h"

//#define DISABLE_OOGS
//#define OGS_ENABLE_TIMER

#ifdef OGS_ENABLE_TIMER
#include "timer.hpp"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "gslib.h"

// hardwired for now
static const unsigned transpose = 0;
static const unsigned recv = 0^transpose, send = 1^transpose;

typedef enum { mode_plain, mode_vec, mode_many,
               mode_dry_run } gs_mode;

struct pw_comm_data {
  uint n;      /* number of messages */
  uint *p;     /* message source/dest proc */
  uint *size;  /* size of message */
  uint total;  /* sum of message sizes */
};

struct pw_data {
  struct pw_comm_data comm[2];
  const uint *map[2];
  comm_req *req;
  uint buffer_size;
};

typedef void exec_fun(
  void *data, gs_mode mode, unsigned vn, gs_dom dom, gs_op op,
  unsigned transpose, const void *execdata, const struct comm *comm, char *buf);
typedef void fin_fun(void *data);

struct gs_remote {
  uint buffer_size, mem_size;
  void *data;
  exec_fun *exec;
  exec_fun *exec_irecv;
  exec_fun *exec_isend;
  exec_fun *exec_wait;
  fin_fun *fin;
};

struct gs_data {
  struct comm comm;
  const uint *map_local[2]; /* 0=unflagged, 1=all */
  const uint *flagged_primaries;
  struct gs_remote r;
  uint handle_size;
};

#ifdef __cplusplus
}
#endif

static void convertPwMap(const uint *restrict map,
                         int *restrict starts,
                         int *restrict ids)
{
  uint i,j; 
  int n=0, s=0;
  while((i=*map++)!=UINT_MAX) { // end of map
    starts[s] = n;
    j=*map++; 
    do {
      ids[n] = j;
      n++;
    } while((j=*map++)!=UINT_MAX); // new message
    starts[s+1] = n;
    s++;
  }
}

static void pairwiseExchange(int unit_size, oogs_t *gs)
{
  ogs_t *ogs = gs->ogs;
  struct gs_data *hgs = (gs_data*) ogs->haloGshSym;
  const void* execdata = hgs->r.data; 
  const struct pw_data *pwd = (pw_data*) execdata; 
  const struct comm *comm = &hgs->comm;

  if(!gs->earlyPrepostRecv) { 
    comm_req *req = pwd->req; 
    const struct pw_comm_data *c = &pwd->comm[recv];
    const uint *p, *pe, *size=c->size;
    uint bufOffset = 0;
    for(p=c->p,pe=p+c->n;p!=pe;++p) {
      const size_t len = *(size++);
      unsigned char *recvbuf = (unsigned char *)gs->bufRecv + bufOffset;
      if(gs->mode == OOGS_DEVICEMPI) recvbuf = (unsigned char*)gs->o_bufRecv.ptr() + bufOffset;
      MPI_Irecv((void*)recvbuf,len*unit_size,MPI_UNSIGNED_CHAR,*p,*p,comm->c,req++);
      bufOffset += len*unit_size;
    }
  }

  { // pw exchange
    if(gs->mode != OOGS_DEVICEMPI) ogs->device.finish(); // waiting for send buffers to be ready

    comm_req *req = &pwd->req[pwd->comm[recv].n]; 
    const struct pw_comm_data *c = &pwd->comm[send];
    const uint *p, *pe, *size=c->size;
    uint bufOffset = 0;
    for(p=c->p,pe=p+c->n;p!=pe;++p) {
      const size_t len = *(size++);
      unsigned char *sendbuf = (unsigned char*)gs->bufSend + bufOffset;
      if(gs->mode == OOGS_DEVICEMPI) sendbuf = (unsigned char*)gs->o_bufSend.ptr() + bufOffset;
      MPI_Isend((void*)sendbuf,len*unit_size,MPI_UNSIGNED_CHAR,*p,comm->id,comm->c,req++);
      bufOffset += len*unit_size;
    }
    MPI_Waitall(pwd->comm[send].n + pwd->comm[recv].n,pwd->req,MPI_STATUSES_IGNORE);
  }
}

oogs_t* oogs::setup(ogs_t *ogs, int nVec, dlong stride, const char *type, std::function<void()> callback, oogs_mode gsMode)
{
  oogs_t *gs = new oogs_t[1];
  gs->ogs = ogs;

  occa::device device = gs->ogs->device;

  struct gs_data *hgs = (gs_data*) ogs->haloGshSym;
  const void* execdata = hgs->r.data;
  const struct pw_data *pwd = (pw_data*) execdata;
  const unsigned unit_size = nVec*sizeof(double); // just need to be big enough to run callcack

  gs->comm = hgs->comm.c;

  int rank;
  MPI_Comm_rank(gs->comm, &rank);
  gs->rank = rank; 
  gs->mode = gsMode;

  if(gsMode == OOGS_DEFAULT) return gs; 

  for (int r=0;r<2;r++) {
    if ((r==0 && gs->rank==0) || (r==1 && gs->rank>0)) {
      gs->packBufFloatAddKernel = device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_floatAdd", ogs::kernelInfo);
      gs->unpackBufFloatAddKernel = device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_floatAdd", ogs::kernelInfo);
      gs->packBufDoubleAddKernel = device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_doubleAdd", ogs::kernelInfo);
      gs->unpackBufDoubleAddKernel = device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_doubleAdd", ogs::kernelInfo);
      gs->packBufDoubleMinKernel = device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_doubleMin", ogs::kernelInfo);
      gs->unpackBufDoubleMinKernel = device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_doubleMin", ogs::kernelInfo);
      gs->packBufDoubleMaxKernel = device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_doubleMax", ogs::kernelInfo);
      gs->unpackBufDoubleMaxKernel = device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_doubleMax", ogs::kernelInfo);

      if(device.mode() == "HIP" || device.mode() == "CUDA") {
        occa::properties halfKernelInfo = ogs::kernelInfo;
        std::string fileName = DOGS;
        if(device.mode() == "CUDA") fileName += "/okl/oogs-half.cu";
        if(device.mode() == "HIP") fileName += "/okl/oogs-half.hip";
        halfKernelInfo["okl/enabled"] = false;
        gs->packBufFloatToHalfAddKernel = device.buildKernel(fileName.c_str(), "packBuf_halfAdd", halfKernelInfo);
        gs->unpackBufHalfToFloatAddKernel = device.buildKernel(fileName.c_str(), "unpackBuf_halfAdd", halfKernelInfo);
      }
    }
    MPI_Barrier(gs->comm);
  }

  if(ogs->NhaloGather == 0) return gs;

  occa::properties props;
  props["mapped"] = true;

  gs->h_buffSend = ogs->device.malloc(pwd->comm[send].total*unit_size, props);
  gs->bufSend = (unsigned char*)gs->h_buffSend.ptr(props); 
  int *scatterOffsets = (int*) calloc(ogs->NhaloGather+1,sizeof(int));
  int *scatterIds = (int*) calloc(pwd->comm[send].total,sizeof(int));
  convertPwMap(pwd->map[send], scatterOffsets, scatterIds);

  gs->o_bufSend = ogs->device.malloc(pwd->comm[send].total*unit_size);
  gs->o_scatterOffsets = ogs->device.malloc((ogs->NhaloGather+1)*sizeof(int), scatterOffsets);
  gs->o_scatterIds = ogs->device.malloc(pwd->comm[send].total*sizeof(int), scatterIds);
  free(scatterOffsets);
  free(scatterIds);

  gs->h_buffRecv = ogs->device.malloc(pwd->comm[recv].total*unit_size, props);
  gs->bufRecv = (unsigned char*)gs->h_buffRecv.ptr(props);
  int* gatherOffsets  = (int*) calloc(ogs->NhaloGather+1,sizeof(int));
  int *gatherIds  = (int*) calloc(pwd->comm[recv].total,sizeof(int));
  convertPwMap(pwd->map[recv], gatherOffsets, gatherIds);

  gs->o_bufRecv = ogs->device.malloc(pwd->comm[recv].total*unit_size);
  gs->o_gatherOffsets  = ogs->device.malloc((ogs->NhaloGather+1)*sizeof(int), gatherOffsets);
  gs->o_gatherIds  = ogs->device.malloc(pwd->comm[recv].total*sizeof(int), gatherIds);
  free(gatherOffsets);
  free(gatherIds);
 
  std::list<oogs_mode> oogs_mode_list;
  oogs_mode_list.push_back(OOGS_DEFAULT);
  oogs_mode_list.push_back(OOGS_HOSTMPI);
  const char* env_val = std::getenv ("OGS_MPI_SUPPORT");
  if(env_val != NULL) { 
    if(std::stoi(env_val)) oogs_mode_list.push_back(OOGS_DEVICEMPI);; 
  }

  if(gsMode == OOGS_AUTO) {
    if(gs->rank == 0) printf("timing oogs modes: ");
    const int Ntests = 10;
    double elapsedMin = std::numeric_limits<double>::max();
    oogs_mode fastestMode;
    occa::memory o_q;
    if(!stride)
      o_q = device.malloc(ogs->N*unit_size);
    else
      o_q = device.malloc(stride*unit_size);

    int* prepostRecv = (int*) calloc(oogs_mode_list.size(), sizeof(int));
    for (auto const& mode : oogs_mode_list)
    {
      gs->mode = mode;
      prepostRecv[gs->mode] = 0; 
      // warum-up
      oogs::start (o_q, nVec, stride, type, ogsAdd, gs);
      if(callback) callback();
      oogs::finish(o_q, nVec, stride, type, ogsAdd, gs);

      for(int pass = 0; pass < 2; pass++) {
        gs->earlyPrepostRecv = pass;
        device.finish();
        MPI_Barrier(gs->comm);
        const double tStart = MPI_Wtime();
        for(int test=0;test<Ntests;++test) {
          oogs::start (o_q, nVec, stride, type, ogsAdd, gs);
          if(callback) callback();
          oogs::finish(o_q, nVec, stride, type, ogsAdd, gs);
        }
        device.finish();
        double elapsed = (MPI_Wtime() - tStart)/Ntests;
        MPI_Allreduce(MPI_IN_PLACE, &elapsed, 1, MPI_DOUBLE, MPI_MAX, gs->comm);
        if(gs->rank == 0) printf("%gs ", elapsed);
        if(elapsed < elapsedMin){
          fastestMode = gs->mode;
          elapsedMin = elapsed;
          if(pass) prepostRecv[fastestMode] = 1;
        }
      }
    }
    MPI_Bcast(&fastestMode, 1, MPI_INT, 0, gs->comm);
    MPI_Bcast(prepostRecv, oogs_mode_list.size(), MPI_INT, 0, gs->comm);
    gs->mode = fastestMode;
    gs->earlyPrepostRecv = prepostRecv[gs->mode];
    o_q.free();
    free(prepostRecv);
  } else {
    gs->mode = gsMode;
    gs->earlyPrepostRecv = 0;
  }

#ifdef DISABLE_OOGS
  gs->mode = OOGS_DEFAULT;
#endif
  if(gs->rank == 0) printf("used mode: %d.%d\n", gs->mode, gs->earlyPrepostRecv);

  return gs; 
}

oogs_t* oogs::setup(dlong N, hlong *ids, int nVec, dlong stride, const char *type, MPI_Comm &comm,
                    int verbose, occa::device device, std::function<void()> callback, oogs_mode gsMode)
{
  ogs_t *ogs  = ogsSetup(N, ids, comm, verbose, device);
  return setup(ogs, nVec, stride, type, callback, gsMode);
}

static void packBuf(oogs_t *gs,
                    const  dlong Ngather,
                    const int k,
                    const dlong stride,
                    occa::memory o_gstarts,
                    occa::memory o_gids,
                    occa::memory o_sstarts,
                    occa::memory o_sids,
                    const char* type,
                    const char* op,
                    occa::memory  o_v,
                    occa::memory  o_gv)
{
  if ((!strcmp(type, "floatCommHalf"))&&(!strcmp(op, ogsAdd))) {
    const int threadBlockSize = 256;
    occa::dim outer, inner;
    outer.dims = 1;
    inner.dims = 1;
    outer[0] = (Ngather * k + threadBlockSize - 1) / threadBlockSize;
    inner[0] = threadBlockSize;
    gs->packBufFloatToHalfAddKernel.setRunDims(outer, inner);
    gs->packBufFloatToHalfAddKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  } else if (!strcmp(type, "float") && !strcmp(op, ogsAdd)) {
    gs->packBufFloatAddKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  } else if (!strcmp(type, "double") && !strcmp(op, ogsAdd)) {
    gs->packBufDoubleAddKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  } else if (!strcmp(type, "double") && !strcmp(op, ogsMin)) {
    gs->packBufDoubleMinKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  } else if (!strcmp(type, "double") && !strcmp(op, ogsMax)) {
    gs->packBufDoubleMaxKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  } else {
    printf("oogs: unsupported operation or datatype!\n");
    exit(1);
  }
}

static void unpackBuf(oogs_t *gs,
                      const  dlong Ngather,
                      const int k,
                      const dlong stride,
                      occa::memory o_gstarts,
                      occa::memory o_gids,
                      occa::memory o_sstarts,
                      occa::memory o_sids,
                      const char* type,
                      const char* op,
                      occa::memory  o_v,
                      occa::memory  o_gv)
{
  if ((!strcmp(type, "floatCommHalf"))&&(!strcmp(op, ogsAdd))) {
    const int threadBlockSize = 256;
    occa::dim outer, inner;
    outer.dims = 1;
    inner.dims = 1;
    outer[0] = (Ngather * k + threadBlockSize - 1) / threadBlockSize;
    inner[0] = threadBlockSize;
    gs->unpackBufHalfToFloatAddKernel.setRunDims(outer, inner);
    gs->unpackBufHalfToFloatAddKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  } else if (!strcmp(type, "float") && !strcmp(op, ogsAdd)) {
    gs->unpackBufFloatAddKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  } else if (!strcmp(type, "double") && !strcmp(op, ogsAdd)) {
    gs->unpackBufDoubleAddKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  } else if (!strcmp(type, "double") && !strcmp(op, ogsMin)) {
    gs->unpackBufDoubleMinKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  } else if (!strcmp(type, "double") && !strcmp(op, ogsMax)) {
    gs->unpackBufDoubleMaxKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  } else {
    printf("oogs: unsupported operation or datatype!\n");
    exit(1);
  }
}

void reallocBuffers(int unit_size, oogs_t *gs)
{
  ogs_t *ogs = gs->ogs; 
  struct gs_data *hgs = (gs_data*) ogs->haloGshSym;
  const void* execdata = hgs->r.data;
  const struct pw_data *pwd = (pw_data*) execdata;

  if (ogs::o_haloBuf.size() < ogs->NhaloGather*unit_size) {
    if (ogs::o_haloBuf.size()) ogs::o_haloBuf.free();
    ogs::haloBuf = ogsHostMallocPinned(ogs->device, ogs->NhaloGather*unit_size, NULL, ogs::o_haloBuf, ogs::h_haloBuf);
  }
  if (gs->o_bufSend.size() < pwd->comm[send].total*unit_size) {
    occa::properties props;
    props["mapped"] = true;
    if(gs->o_bufSend.size()) gs->o_bufSend.free();
    gs->o_bufSend = ogs->device.malloc(pwd->comm[send].total*unit_size);
    if(gs->h_buffSend.size()) gs->h_buffSend.free();
    gs->h_buffSend = ogs->device.malloc(pwd->comm[send].total*unit_size, props);
    gs->bufSend = (unsigned char*)gs->h_buffSend.ptr(props);
  }
  if (gs->o_bufRecv.size() < pwd->comm[recv].total*unit_size) {
    occa::properties props;
    props["mapped"] = true;
    if(gs->o_bufRecv.size()) gs->o_bufRecv.free();
    gs->o_bufRecv = ogs->device.malloc(pwd->comm[recv].total*unit_size);
    if(gs->h_buffRecv.size()) gs->h_buffRecv.free();
    gs->h_buffRecv = ogs->device.malloc(pwd->comm[recv].total*unit_size, props);
    gs->bufRecv = (unsigned char*)gs->h_buffRecv.ptr(props);
  }
}

void oogs::start(occa::memory o_v, const int k, const dlong stride, const char *_type, const char *op, oogs_t *gs) 
{
  size_t Nbytes;
  ogs_t *ogs = gs->ogs; 
  const char* type = (!strcmp(_type,"floatCommHalf")) ? "float" : _type;
  if (!strcmp(_type, "floatCommHalf"))
    Nbytes = sizeof(float)/2;
  else if (!strcmp(type, "float"))
    Nbytes = sizeof(float);
  else if (!strcmp(type, "double"))
    Nbytes = sizeof(double);
  else if (!strcmp(type, "int"))
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int"))
    Nbytes = sizeof(long long int);

  if (!strcmp(_type, "floatCommHalf") && ogs->device.mode() == "SERIAL"){
    if(gs->rank == 0)
      std::cout << "ERROR: Current backend does not support floatCommHalf!\n";
    exit(-1);
  }

  if(gs->mode == OOGS_DEFAULT) { 
    if(k>1)
      ogsGatherScatterManyStart(o_v, k, stride, type, op, ogs);
    else
      ogsGatherScatterStart(o_v, type, op, ogs);
    
    return;
  }

  if (ogs->NhaloGather) {
    reallocBuffers(Nbytes*k, gs);

    packBuf(gs, ogs->NhaloGather, k, stride, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, 
            gs->o_scatterOffsets, gs->o_scatterIds, _type, op, o_v, gs->o_bufSend);

    if(gs->earlyPrepostRecv) { 
      struct gs_data *hgs = (gs_data*) ogs->haloGshSym;
      const void* execdata = hgs->r.data;
      const struct pw_data *pwd = (pw_data*) execdata;

      comm_req *req = pwd->req;
      int unit_size = Nbytes*k;
      const struct pw_comm_data *c = &pwd->comm[recv];
      const uint *p, *pe, *size=c->size;
      uint bufOffset = 0;
      for(p=c->p,pe=p+c->n;p!=pe;++p) {
        const size_t len = *(size++);
        unsigned char *recvbuf = (unsigned char *)gs->bufRecv + bufOffset;
        if(gs->mode == OOGS_DEVICEMPI) recvbuf = (unsigned char*)gs->o_bufRecv.ptr() + bufOffset;
        MPI_Irecv((void*)recvbuf,len*unit_size,MPI_UNSIGNED_CHAR,*p,*p,gs->comm,req++);
        bufOffset += len*unit_size;
      }
    }

    ogs->device.finish();
  }
}

void oogs::finish(occa::memory o_v, const int k, const dlong stride, const char *_type, const char *op, oogs_t *gs) 
{
  size_t Nbytes;
  ogs_t *ogs = gs->ogs; 
  const char* type = (!strcmp(_type,"floatCommHalf")) ? "float" : _type;
  if (!strcmp(_type, "floatCommHalf"))
    Nbytes = sizeof(float)/2;
  else if (!strcmp(_type, "float"))
    Nbytes = sizeof(float);
  else if (!strcmp(_type, "double"))
    Nbytes = sizeof(double);

  if (!strcmp(_type, "floatCommHalf") && ogs->device.mode() == "SERIAL"){
    if(gs->rank == 0)
      std::cout << "ERROR: Current backend does not support floatCommHalf!\n";
    exit(-1);
  }

  if(gs->mode == OOGS_DEFAULT) { 
    if(k>1)
      ogsGatherScatterManyFinish(o_v, k, stride, type, op, ogs);
    else
      ogsGatherScatterFinish(o_v, type, op, ogs);
    
    return;
  }

  if(ogs->NlocalGather) 
    occaGatherScatterMany(ogs->NlocalGather, k, stride, ogs->o_localGatherOffsets, 
		          ogs->o_localGatherIds, type, op, o_v);

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);

    struct gs_data *hgs = (gs_data*) ogs->haloGshSym;
    const void* execdata = hgs->r.data;
    const struct pw_data *pwd = (pw_data*) execdata;
    if(gs->mode == OOGS_HOSTMPI)
      gs->o_bufSend.copyTo(gs->bufSend, pwd->comm[send].total*Nbytes*k, 0, "async: true");

#ifdef OGS_ENABLE_TIMER
    timer::hostTic("oogsMPI",1);
#endif
    pairwiseExchange(Nbytes*k, gs);
#ifdef OGS_ENABLE_TIMER
    timer::hostToc("oogsMPI");
#endif

    if(gs->mode == OOGS_HOSTMPI)
      gs->o_bufRecv.copyFrom(gs->bufRecv,pwd->comm[recv].total*Nbytes*k, 0, "async: true");

    unpackBuf(gs, ogs->NhaloGather, k, stride, gs->o_gatherOffsets, gs->o_gatherIds, 
              ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, _type, op, gs->o_bufRecv, o_v);
    ogs->device.finish();

    ogs->device.setStream(ogs::defaultStream);
  }
}

void oogs::startFinish(void *v, const int k, const dlong stride, const char *type, const char *op, oogs_t *h)
{
  ogsGatherScatterMany(v, k, stride, type, op, h->ogs);
}
void oogs::startFinish(occa::memory o_v, const int k, const dlong stride, const char *type, const char *op, oogs_t *h)
{
   start(o_v, k, stride, type, op, h);
   finish(o_v, k, stride, type, op, h);
}

void oogs::destroy(oogs_t *gs)
{
  //ogsFree(gs->ogs);

  gs->h_buffSend.free();
  gs->h_buffRecv.free();

  gs->o_scatterIds.free();
  gs->o_gatherIds.free();

  gs->o_scatterOffsets.free();
  gs->o_gatherOffsets.free();

  gs->o_bufRecv.free();
  gs->o_bufSend.free();

  gs->packBufFloatAddKernel.free();
  gs->unpackBufFloatAddKernel.free();
  gs->packBufDoubleAddKernel.free() ;
  gs->unpackBufDoubleAddKernel.free();
  gs->packBufDoubleMinKernel.free();
  gs->unpackBufDoubleMinKernel.free();
  gs->packBufDoubleMaxKernel.free();
  gs->unpackBufDoubleMaxKernel.free();

  gs->packBufFloatToHalfAddKernel.free();
  gs->unpackBufHalfToFloatAddKernel.free();
  
  free(gs);
}
