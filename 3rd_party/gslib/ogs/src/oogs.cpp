#include <limits>
#include <list>
#include <occa.hpp>

#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "ogsInterface.h"

//#define ENABLE_EARLY_PREPOST

#ifdef __cplusplus
extern "C" {
#endif

#include "gslib.h"

// hardwired for now
static const unsigned transpose = 0;
static const unsigned recv = 0 ^ transpose, send = 1 ^ transpose;

static int OGS_MPI_SUPPORT = 0;
static int OGS_OVERLAP = 1;
static int compiled = 0;

typedef enum { mode_plain, mode_vec, mode_many, mode_dry_run } gs_mode;

struct pw_comm_data {
  uint n;     /* number of messages */
  uint *p;    /* message source/dest proc */
  uint *size; /* size of message */
  uint total; /* sum of message sizes */
};

struct pw_data {
  struct pw_comm_data comm[2];
  const uint *map[2];
  comm_req *req;
  uint buffer_size;
};

typedef void exec_fun(void *data,
                      gs_mode mode,
                      unsigned vn,
                      gs_dom dom,
                      gs_op op,
                      unsigned transpose,
                      const void *execdata,
                      const struct comm *comm,
                      char *buf);
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

static void convertPwMap(const uint *restrict map, int *restrict starts, int *restrict ids)
{
  uint i, j;
  int n = 0, s = 0;
  while ((i = *map++) != UINT_MAX) { // end of map
    starts[s] = n;
    j = *map++;
    do {
      ids[n] = j;
      n++;
    } while ((j = *map++) != UINT_MAX); // new message
    starts[s + 1] = n;
    s++;
  }
}

static void neighborAllToAll(int unit_size, oogs_t *gs)
{
  ogs_t *ogs = gs->ogs;
  struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
  const void *execdata = hgs->r.data;
  const struct pw_data *pwd = (pw_data *)execdata;

  {
    uint bufOffset = 0;
    const struct pw_comm_data *c = &pwd->comm[send];
    for (int i = 0; i < c->n; i++) {
      const int len = c->size[i] * unit_size;
      gs->nbc.sendcounts[i] = len;
      gs->nbc.senddispls[i] = bufOffset;
      bufOffset += len;
    }
  }
  {
    uint bufOffset = 0;
    const struct pw_comm_data *c = &pwd->comm[recv];
    for (int i = 0; i < c->n; i++) {
      const int len = c->size[i] * unit_size;
      gs->nbc.recvcounts[i] = len;
      gs->nbc.recvdispls[i] = bufOffset;
      bufOffset += len;
    }
  }

  unsigned char *bufRecv = (unsigned char *)gs->o_bufRecv.ptr();
  unsigned char *bufSend = (unsigned char *)gs->o_bufSend.ptr();
  if (gs->mode != OOGS_DEVICEMPI) {
    ogs->device.finish(); // waiting for send buffers to be ready
    bufRecv = (unsigned char *)gs->bufRecv;
    bufSend = (unsigned char *)gs->bufSend;
  }
  MPI_Neighbor_alltoallv(bufSend,
                         gs->nbc.sendcounts,
                         gs->nbc.senddispls,
                         MPI_UNSIGNED_CHAR,
                         bufRecv,
                         gs->nbc.recvcounts,
                         gs->nbc.recvdispls,
                         MPI_UNSIGNED_CHAR,
                         gs->nbc.comm);
}

static void pairwiseExchange(int unit_size, oogs_t *gs)
{
  ogs_t *ogs = gs->ogs;
  struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
  const void *execdata = hgs->r.data;
  const struct pw_data *pwd = (pw_data *)execdata;

  if (!gs->earlyPrepostRecv) {
    unsigned char *buf = (unsigned char *)gs->o_bufRecv.ptr();
    if (gs->mode != OOGS_DEVICEMPI)
      buf = (unsigned char *)gs->bufRecv;

    comm_req *req = pwd->req;
    const struct pw_comm_data *c = &pwd->comm[recv];
    const uint *p, *pe, *size = c->size;
    for (p = c->p, pe = p + c->n; p != pe; ++p) {
      const int len = *(size++) * unit_size;
      MPI_Irecv((void *)buf, len, MPI_UNSIGNED_CHAR, *p, *p, gs->comm, req++);
      buf += len;
    }
  }

  {
    unsigned char *buf = (unsigned char *)gs->o_bufSend.ptr();
    if (gs->mode != OOGS_DEVICEMPI) {
      ogs->device.finish(); // waiting for send buffers to be ready
      buf = (unsigned char *)gs->bufSend;
    }

    comm_req *req = &pwd->req[pwd->comm[recv].n];
    const struct pw_comm_data *c = &pwd->comm[send];
    const uint *p, *pe, *size = c->size;
    for (p = c->p, pe = p + c->n; p != pe; ++p) {
      const int len = *(size++) * unit_size;
      MPI_Isend((void *)buf, len, MPI_UNSIGNED_CHAR, *p, gs->rank, gs->comm, req++);
      buf += len;
    }
    MPI_Waitall(pwd->comm[send].n + pwd->comm[recv].n, pwd->req, MPI_STATUSES_IGNORE);
  }
}
void occaGatherScatterLocal(const dlong NlocalGather,
                            const dlong NrowBlocks,
                            occa::memory &o_bstart,
                            occa::memory &o_gstart,
                            occa::memory &o_gids,
                            const int Nvectors,
                            const dlong stride,
                            const char *type,
                            const char *op,
                            occa::memory &o_v)
{
#if 1
  occaGatherScatterMany(NlocalGather, Nvectors, stride, o_gstart, o_gids, type, op, o_v);
#else
  const int Nentries = 1;
  if (!strcmp(type, "float") && !strcmp(op, "add")) {
    ogs::gatherScatterNewKernel_floatAdd(NrowBlocks,
                                         Nentries,
                                         Nvectors,
                                         stride,
                                         o_bstart,
                                         o_gstart,
                                         o_gids,
                                         o_v);
  }
  else if (!strcmp(type, "double") && !strcmp(op, "add")) {
    ogs::gatherScatterNewKernel_doubleAdd(NrowBlocks,
                                          Nentries,
                                          Nvectors,
                                          stride,
                                          o_bstart,
                                          o_gstart,
                                          o_gids,
                                          o_v);
  }
  else if (!strcmp(type, "double") && !strcmp(op, "min")) {
    ogs::gatherScatterNewKernel_doubleMin(NrowBlocks,
                                          Nentries,
                                          Nvectors,
                                          stride,
                                          o_bstart,
                                          o_gstart,
                                          o_gids,
                                          o_v);
  }
  else if (!strcmp(type, "double") && !strcmp(op, "max")) {
    ogs::gatherScatterNewKernel_doubleMax(NrowBlocks,
                                          Nentries,
                                          Nvectors,
                                          stride,
                                          o_bstart,
                                          o_gstart,
                                          o_gids,
                                          o_v);
  }
  else {
    printf("occaGatherScatterNewKernel: unsupported operation or datatype!\n");
    exit(1);
  }
#endif
}

void oogs::gpu_mpi(int val) { OGS_MPI_SUPPORT = val; }

void oogs::overlap(int val) { OGS_OVERLAP = val; }

int oogs::gpu_mpi() { return OGS_MPI_SUPPORT; }

void oogs::compile(const occa::device &device, std::string mode, MPI_Comm comm, bool verbose)
{
  ogs::initKernels(comm, device, verbose);
  occa::properties props = ogs::kernelInfo;
  if (verbose) {
    props["verbose"] = true;
  }

  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_floatAdd", props);
    device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_floatAdd", props);
    device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_doubleAdd", props);
    device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_doubleAdd", props);
    device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_doubleMin", props);
    device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_doubleMin", props);
    device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_doubleMax", props);
    device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_doubleMax", props);
  }
  compiled++;
}

void reallocBuffers(int unit_size, oogs_t *gs)
{
  ogs_t *ogs = gs->ogs;
  struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
  const void *execdata = hgs->r.data;
  const struct pw_data *pwd = (pw_data *)execdata;

  if (gs->o_bufSend.size() < pwd->comm[send].total * unit_size) {
    if (gs->o_bufSend.size())
      gs->o_bufSend.free();
    if (gs->h_buffSend.size())
      gs->h_buffSend.free();
    gs->bufSend = (unsigned char *)ogsHostMallocPinned(ogs->device,
                                                       pwd->comm[send].total * unit_size,
                                                       NULL,
                                                       gs->o_bufSend,
                                                       gs->h_buffSend);
  }
  if (gs->o_bufRecv.size() < pwd->comm[recv].total * unit_size) {
    if (gs->o_bufRecv.size())
      gs->o_bufRecv.free();
    if (gs->h_buffRecv.size())
      gs->h_buffRecv.free();
    gs->bufRecv = (unsigned char *)ogsHostMallocPinned(ogs->device,
                                                       pwd->comm[recv].total * unit_size,
                                                       NULL,
                                                       gs->o_bufRecv,
                                                       gs->h_buffRecv);
  }
}

size_t typeBytes(const char *type)
{
  if (!strcmp(type, "float"))
    return sizeof(float);
  else if (!strcmp(type, "double"))
    return sizeof(double);
  else if (!strcmp(type, "int"))
    return sizeof(int);
  else if (!strcmp(type, "long long int"))
    return sizeof(long long int);
  else
    return -1;
}

oogs_t *oogs::setup(ogs_t *ogs,
                    int nVec,
                    dlong stride,
                    const char *type,
                    std::function<void()> callback,
                    oogs_mode gsMode)
{
  const size_t Nbytes = typeBytes(type);

  oogs_t *gs = new oogs_t[1];
  gs->ogs = ogs;

  occa::device device = gs->ogs->device;

  struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
  const void *execdata = hgs->r.data;
  const struct pw_data *pwd = (pw_data *)execdata;

  gs->comm = hgs->comm.c;

  int rank;
  MPI_Comm_rank(gs->comm, &rank);
  gs->rank = rank;
  gs->mode = gsMode;

  if (!compiled)
    oogs::compile(device, device.mode(), gs->comm);

  if (gsMode == OOGS_DEFAULT)
    return gs;

  gs->packBufFloatAddKernel = device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_floatAdd", ogs::kernelInfo);
  gs->unpackBufFloatAddKernel =
      device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_floatAdd", ogs::kernelInfo);
  gs->packBufDoubleAddKernel = device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_doubleAdd", ogs::kernelInfo);
  gs->unpackBufDoubleAddKernel =
      device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_doubleAdd", ogs::kernelInfo);
  gs->packBufDoubleMinKernel = device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_doubleMin", ogs::kernelInfo);
  gs->unpackBufDoubleMinKernel =
      device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_doubleMin", ogs::kernelInfo);
  gs->packBufDoubleMaxKernel = device.buildKernel(DOGS "/okl/oogs.okl", "packBuf_doubleMax", ogs::kernelInfo);
  gs->unpackBufDoubleMaxKernel =
      device.buildKernel(DOGS "/okl/oogs.okl", "unpackBuf_doubleMax", ogs::kernelInfo);

  std::list<oogs_mode> oogs_mode_list;
  oogs_mode_list.push_back(OOGS_LOCAL);
  std::list<oogs_modeExchange> oogs_modeExchange_list;
  oogs_modeExchange_list.push_back(OOGS_EX_PW);

  if (ogs->NhaloGather > 0) {
    gs->bufSend = (unsigned char *)ogsHostMallocPinned(ogs->device,
                                                       pwd->comm[send].total * nVec * sizeof(double),
                                                       NULL,
                                                       gs->o_bufSend,
                                                       gs->h_buffSend);
    int *scatterOffsets = (int *)calloc(ogs->NhaloGather + 1, sizeof(int));
    int *scatterIds = (int *)calloc(pwd->comm[send].total, sizeof(int));
    convertPwMap(pwd->map[send], scatterOffsets, scatterIds);
    gs->o_scatterOffsets = ogs->device.malloc((ogs->NhaloGather + 1) * sizeof(int), scatterOffsets);
    gs->o_scatterIds = ogs->device.malloc(pwd->comm[send].total * sizeof(int), scatterIds);
    free(scatterOffsets);
    free(scatterIds);

    gs->bufRecv = (unsigned char *)ogsHostMallocPinned(ogs->device,
                                                       pwd->comm[recv].total * nVec * sizeof(double),
                                                       NULL,
                                                       gs->o_bufRecv,
                                                       gs->h_buffRecv);
    int *gatherOffsets = (int *)calloc(ogs->NhaloGather + 1, sizeof(int));
    int *gatherIds = (int *)calloc(pwd->comm[recv].total, sizeof(int));
    convertPwMap(pwd->map[recv], gatherOffsets, gatherIds);
    gs->o_gatherOffsets = ogs->device.malloc((ogs->NhaloGather + 1) * sizeof(int), gatherOffsets);
    gs->o_gatherIds = ogs->device.malloc(pwd->comm[recv].total * sizeof(int), gatherIds);
    free(gatherOffsets);
    free(gatherIds);

    const int reorder = 0;
    int *src = (int *)calloc(pwd->comm[recv].n, sizeof(int));
    int *dst = (int *)calloc(pwd->comm[send].n, sizeof(int));
    for (int i = 0; i < pwd->comm[recv].n; ++i) {
      src[i] = pwd->comm[recv].p[i];
    }
    for (int i = 0; i < pwd->comm[send].n; ++i) {
      dst[i] = pwd->comm[send].p[i];
    }
    MPI_Dist_graph_create_adjacent(gs->comm,
                                   pwd->comm[recv].n,
                                   src,
                                   MPI_UNWEIGHTED,
                                   pwd->comm[send].n,
                                   dst,
                                   MPI_UNWEIGHTED,
                                   MPI_INFO_NULL,
                                   reorder,
                                   &gs->nbc.comm);
    free(src);
    free(dst);
    gs->nbc.sendcounts = (int *)calloc(pwd->comm[send].n, sizeof(int));
    gs->nbc.senddispls = (int *)calloc(pwd->comm[send].n, sizeof(int));
    gs->nbc.recvcounts = (int *)calloc(pwd->comm[recv].n, sizeof(int));
    gs->nbc.recvdispls = (int *)calloc(pwd->comm[recv].n, sizeof(int));

    oogs_mode_list.push_back(OOGS_DEFAULT);
    oogs_mode_list.push_back(OOGS_HOSTMPI);
    if (ogs->device.mode() != "Serial") {
      if(OGS_MPI_SUPPORT) oogs_mode_list.push_back(OOGS_DEVICEMPI);
    }
    oogs_modeExchange_list.push_back(OOGS_EX_NBC);
  }

  if (gsMode == OOGS_AUTO) {
    if (gs->rank == 0)
      printf("timing gs: ");
    const int Ntests = 10;
    double elapsedMin = std::numeric_limits<double>::max();
    oogs_mode fastestMode = OOGS_DEFAULT;
    oogs_modeExchange fastestModeExchange = OOGS_EX_PW;
    int fastestPrepostRecv = 0;

    void * q = calloc(std::max(stride, ogs->N) * nVec * Nbytes, 1);
    occa::memory o_q = device.malloc(std::max(stride, ogs->N) * nVec * Nbytes, q);
    free(q);

    for (auto const &mode : oogs_mode_list) {
      gs->mode = mode;

      for (auto const &modeExchange : oogs_modeExchange_list) {
        gs->modeExchange = modeExchange;

        if (gs->modeExchange == OOGS_EX_NBC && gs->mode == OOGS_DEVICEMPI)
          continue; // not supported yet by all MPI implementations

#ifdef ENABLE_EARLY_PREPOST
        const int nPass = 2;
#else
        const int nPass = 1;
#endif

        MPI_Barrier(gs->comm);
        for (int pass = 0; pass < nPass; pass++) {
          gs->earlyPrepostRecv = pass;

          // skip invalid combinations
          if (gs->modeExchange != OOGS_EX_PW && gs->earlyPrepostRecv)
            continue;
          if (gs->mode == OOGS_DEFAULT || gs->mode == OOGS_LOCAL) {
            if (gs->modeExchange != OOGS_EX_PW)
              continue;
            if (gs->earlyPrepostRecv)
              continue;
          }

          // run Ntests measurements to eliminate runtime variations
          double elapsedTest = std::numeric_limits<double>::max();
          for (int test = 0; test < Ntests; ++test) {
            device.finish();
            MPI_Barrier(gs->comm);
            const double tStart = MPI_Wtime();

            oogs::start(o_q, nVec, stride, type, ogsAdd, gs);
            if (callback)
              callback();
            oogs::finish(o_q, nVec, stride, type, ogsAdd, gs);

            device.finish();
            elapsedTest = std::min(elapsedTest, MPI_Wtime() - tStart);
          }
          MPI_Allreduce(MPI_IN_PLACE, &elapsedTest, 1, MPI_DOUBLE, MPI_MAX, gs->comm);

          if (gs->rank == 0)
            printf("%.2es ", elapsedTest);
          fflush(stdout);

          if (elapsedTest < elapsedMin) {
            if (gs->mode != OOGS_LOCAL) {
              elapsedMin = elapsedTest;
              fastestMode = gs->mode;
              fastestModeExchange = gs->modeExchange;
              fastestPrepostRecv = gs->earlyPrepostRecv;
            }
          }
        }
      }
    }
    MPI_Bcast(&fastestMode, 1, MPI_INT, 0, gs->comm);
    MPI_Bcast(&fastestModeExchange, 1, MPI_INT, 0, gs->comm);
    MPI_Bcast(&fastestPrepostRecv, 1, MPI_INT, 0, gs->comm);
    gs->mode = fastestMode;
    gs->modeExchange = fastestModeExchange;
    gs->earlyPrepostRecv = fastestPrepostRecv;
    o_q.free();
  }
  else {
    gs->mode = gsMode;
    gs->modeExchange = OOGS_EX_PW;
    gs->earlyPrepostRecv = 0;
  }

#ifdef DISABLE_OOGS
  gs->mode = OOGS_DEFAULT;
#endif


  double elapsedMinMPI = std::numeric_limits<double>::max();
  {
    const int earlyPrepostRecv = gs->earlyPrepostRecv;
    gs->earlyPrepostRecv = 0;
    const int Ntests = 10;

    const size_t unit_size = nVec * Nbytes;
    reallocBuffers(unit_size, gs);

    for (int test = 0; test < Ntests; ++test) {
      device.finish();
      MPI_Barrier(gs->comm);
      const double tStart = MPI_Wtime();
      if (gs->modeExchange == OOGS_EX_NBC)
        neighborAllToAll(unit_size, gs);
      else
        pairwiseExchange(unit_size, gs);
      elapsedMinMPI = std::min(elapsedMinMPI, MPI_Wtime() - tStart);
    }
    gs->earlyPrepostRecv = earlyPrepostRecv;

    int size;
    MPI_Comm_size(gs->comm, &size);
    double nBytesExchange = unit_size * (pwd->comm[send].total + pwd->comm[recv].total);
    MPI_Allreduce(MPI_IN_PLACE, &nBytesExchange, 1, MPI_DOUBLE, MPI_SUM, gs->comm);

    double tmin, tmax, tavg;
    MPI_Allreduce(&elapsedMinMPI, &tmin, 1, MPI_DOUBLE, MPI_MIN, gs->comm);
    MPI_Allreduce(&elapsedMinMPI, &tmax, 1, MPI_DOUBLE, MPI_MAX, gs->comm);
    MPI_Allreduce(&elapsedMinMPI, &tavg, 1, MPI_DOUBLE, MPI_SUM, gs->comm);

    std::string configStr = (gs->modeExchange == OOGS_EX_NBC) ? "nbc" : "pw";
    configStr += (gs->earlyPrepostRecv) ? "+early" : "";
    if (gs->rank == 0) {
      if (ogs->NhaloGather > 0) {
        switch (gs->mode) {
        case OOGS_DEFAULT:
          if (ogs->device.mode() != "Serial") configStr += "+host";
          break;
        case OOGS_HOSTMPI:
          if (ogs->device.mode() != "Serial")  configStr += "+hybrid"; 
          break;
        case OOGS_DEVICEMPI:
           configStr += "+device";
          break;
        }
        printf("\nused config: %s ", configStr.c_str());
        if (tavg/size > MPI_Wtick())
          printf("(MPI min/max/avg: %.2es %.2es %.2es / avg bi-bw: %.1fGB/s/rank)\n",
                 tmin, tmax, tavg/size,
                 nBytesExchange / tavg / 1e9);
        else
          printf("\n");
      }
      else {
        printf("\nused config: local\n");
      }
    }
    fflush(stdout);
  }
  return gs;
}

oogs_t *oogs::setup(dlong N,
                    hlong *ids,
                    int nVec,
                    dlong stride,
                    const char *type,
                    MPI_Comm &comm,
                    int verbose,
                    occa::device device,
                    std::function<void()> callback,
                    oogs_mode gsMode)
{
  ogs_t *ogs = ogsSetup(N, ids, comm, verbose, device);
  return setup(ogs, nVec, stride, type, callback, gsMode);
}

static void packBuf(oogs_t *gs,
                    const dlong Ngather,
                    const int k,
                    const dlong stride,
                    occa::memory &o_gstarts,
                    occa::memory &o_gids,
                    occa::memory &o_sstarts,
                    occa::memory &o_sids,
                    const char *type,
                    const char *op,
                    occa::memory &o_v,
                    occa::memory &o_gv)
{
  if (!strcmp(type, "float") && !strcmp(op, ogsAdd)) {
    gs->packBufFloatAddKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  }
  else if (!strcmp(type, "double") && !strcmp(op, ogsAdd)) {
    gs->packBufDoubleAddKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  }
  else if (!strcmp(type, "double") && !strcmp(op, ogsMin)) {
    gs->packBufDoubleMinKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  }
  else if (!strcmp(type, "double") && !strcmp(op, ogsMax)) {
    gs->packBufDoubleMaxKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  }
  else {
    printf("oogs: unsupported operation or datatype!\n");
    exit(1);
  }
}

static void unpackBuf(oogs_t *gs,
                      const dlong Ngather,
                      const int k,
                      const dlong stride,
                      occa::memory &o_gstarts,
                      occa::memory &o_gids,
                      occa::memory &o_sstarts,
                      occa::memory &o_sids,
                      const char *type,
                      const char *op,
                      occa::memory &o_v,
                      occa::memory &o_gv)
{
  if (!strcmp(type, "float") && !strcmp(op, ogsAdd)) {
    gs->unpackBufFloatAddKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  }
  else if (!strcmp(type, "double") && !strcmp(op, ogsAdd)) {
    gs->unpackBufDoubleAddKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  }
  else if (!strcmp(type, "double") && !strcmp(op, ogsMin)) {
    gs->unpackBufDoubleMinKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  }
  else if (!strcmp(type, "double") && !strcmp(op, ogsMax)) {
    gs->unpackBufDoubleMaxKernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
  }
  else {
    printf("oogs: unsupported operation or datatype!\n");
    exit(1);
  }
}

void oogs::start(occa::memory &o_v,
                 const int k,
                 const dlong stride,
                 const char *type,
                 const char *op,
                 oogs_t *gs)
{
  ogs_t *ogs = gs->ogs;
  const int unit_size = (int) typeBytes(type) * k;

  if (gs->mode == OOGS_DEFAULT) {
    if (k > 1)
      ogsGatherScatterManyStart(o_v, k, stride, type, op, ogs);
    else
      ogsGatherScatterStart(o_v, type, op, ogs);

    return;
  }

  if (ogs->NhaloGather && gs->mode != OOGS_LOCAL) {
    reallocBuffers(unit_size, gs);

    packBuf(gs,
            ogs->NhaloGather,
            k,
            stride,
            ogs->o_haloGatherOffsets,
            ogs->o_haloGatherIds,
            gs->o_scatterOffsets,
            gs->o_scatterIds,
            type,
            op,
            o_v,
            gs->o_bufSend);

    if (gs->earlyPrepostRecv) {
      unsigned char *buf = (unsigned char *)gs->o_bufRecv.ptr();
      if (gs->mode != OOGS_DEVICEMPI)
        buf = (unsigned char *)gs->bufRecv;

      struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
      const void *execdata = hgs->r.data;
      const struct pw_data *pwd = (pw_data *)execdata;

      comm_req *req = pwd->req;
      const struct pw_comm_data *c = &pwd->comm[recv];
      const uint *p, *pe, *size = c->size;
      for (p = c->p, pe = p + c->n; p != pe; ++p) {
        const size_t len = *(size++) * unit_size;
        MPI_Irecv((void *)buf, len * unit_size, MPI_UNSIGNED_CHAR, *p, *p, gs->comm, req++);
        buf += len;
      }
    }

    ogs->device.finish();
  }
}

void oogs::finish(occa::memory &o_v,
                  const int k,
                  const dlong stride,
                  const char *type,
                  const char *op,
                  oogs_t *gs)
{
  ogs_t *ogs = gs->ogs;
  const int unit_size = (int) typeBytes(type) * k;

  if (gs->mode == OOGS_DEFAULT) {
    if (k > 1)
      ogsGatherScatterManyFinish(o_v, k, stride, type, op, ogs);
    else
      ogsGatherScatterFinish(o_v, type, op, ogs);

    return;
  }

  if (ogs->NlocalGather)
    occaGatherScatterLocal(ogs->NlocalGather,
                           ogs->NrowBlocks,
                           ogs->o_blockRowStarts,
                           ogs->o_localGatherOffsets,
                           ogs->o_localGatherIds,
                           k,
                           stride,
                           type,
                           op,
                           o_v);


  if (ogs->NhaloGather && !OGS_OVERLAP)
    ogs->device.finish();

  if (ogs->NhaloGather && gs->mode == OOGS_HOSTMPI) {
    ogs->device.setStream(ogs::dataStream);

    struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
    const void *execdata = hgs->r.data;
    const struct pw_data *pwd = (pw_data *)execdata;

    gs->o_bufSend.copyTo(gs->bufSend, pwd->comm[send].total * unit_size, 0, "async: true");

    // syncs device on exit (MPI is not stream-aware)
    ogsHostTic(gs->comm, 1);
    if (gs->modeExchange == OOGS_EX_NBC)
      neighborAllToAll(unit_size, gs);
    else
      pairwiseExchange(unit_size, gs);
    ogsHostToc();

    gs->o_bufRecv.copyFrom(gs->bufRecv, pwd->comm[recv].total * unit_size, 0, "async: true");

    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);
  }


  if (ogs->NhaloGather && gs->mode == OOGS_DEVICEMPI) {
    // syncs device on exit (MPI is not stream-aware)
    ogsHostTic(gs->comm, 1);
    if (gs->modeExchange == OOGS_EX_NBC)
      neighborAllToAll(unit_size, gs);
    else
      pairwiseExchange(unit_size, gs);
    ogsHostToc();
  }

  if (ogs->NhaloGather)
     unpackBuf(gs,
              ogs->NhaloGather,
              k,
              stride,
              gs->o_gatherOffsets,
              gs->o_gatherIds,
              ogs->o_haloGatherOffsets,
              ogs->o_haloGatherIds,
              type,
              op,
              gs->o_bufRecv,
              o_v);

}

void oogs::startFinish(void *v, const int k, const dlong stride, const char *type, const char *op, oogs_t *h)
{
  ogsGatherScatterMany(v, k, stride, type, op, h->ogs);
}

void oogs::startFinish(occa::memory &o_v,
                       const int k,
                       const dlong stride,
                       const char *type,
                       const char *op,
                       oogs_t *h)
{
  start(o_v, k, stride, type, op, h);
  finish(o_v, k, stride, type, op, h);
}

void oogs::destroy(oogs_t *gs)
{
  // ogsFree(gs->ogs);

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
  gs->packBufDoubleAddKernel.free();
  gs->unpackBufDoubleAddKernel.free();
  gs->packBufDoubleMinKernel.free();
  gs->unpackBufDoubleMinKernel.free();
  gs->packBufDoubleMaxKernel.free();
  gs->unpackBufDoubleMaxKernel.free();

  free(gs);
}
