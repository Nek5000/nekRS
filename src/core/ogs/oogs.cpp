#include <limits>
#include <list>

#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "ogsInterface.h"

#include "platform.hpp"

#ifdef __cplusplus
extern "C" {
#endif

#include "gslib.h"

#define MPI_CHECK(x)                                                                                         \
do {                                                                                                         \
int __ret = (x);                                                                                             \
if (MPI_SUCCESS != __ret) {                                                                                  \
fprintf(stderr, "(%s:%d) ERROR: MPI call returned error code %d", __FILE__, __LINE__, __ret);                \
MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);                                                                     \
}                                                                                                            \
} while (0)

// hardwired for now
static constexpr unsigned transpose = 0;
static constexpr unsigned recv = 0 ^ transpose;
static constexpr unsigned send = 1 ^ transpose;

static int OGS_MPI_SUPPORT = 0;
static int OGS_OVERLAP = 1;
static int OGS_SYNC_RECV = 0;

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

namespace oogs
{
occa::kernel packBufFloatAddKernel, unpackBufFloatAddKernel;
occa::kernel packBufFloatMinKernel, unpackBufFloatMinKernel;
occa::kernel packBufFloatMaxKernel, unpackBufFloatMaxKernel;

occa::kernel packBufDoubleAddKernel, unpackBufDoubleAddKernel;
occa::kernel packBufDoubleMinKernel, unpackBufDoubleMinKernel;
occa::kernel packBufDoubleMaxKernel, unpackBufDoubleMaxKernel;
} // namespace oogs

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
  MPI_CHECK(MPI_Neighbor_alltoallv(bufSend,
                                   gs->nbc.sendcounts,
                                   gs->nbc.senddispls,
                                   MPI_UNSIGNED_CHAR,
                                   bufRecv,
                                   gs->nbc.recvcounts,
                                   gs->nbc.recvdispls,
                                   MPI_UNSIGNED_CHAR,
                                   gs->nbc.comm));
}

static void pairwiseExchange(int unit_size, oogs_t *gs)
{
  ogs_t *ogs = gs->ogs;
  struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
  const void *execdata = hgs->r.data;
  const struct pw_data *pwd = (pw_data *)execdata;

  if (!gs->earlyPrepostRecv) {
    unsigned char *buf = (unsigned char *)gs->o_bufRecv.ptr();
    if (gs->mode != OOGS_DEVICEMPI) {
      buf = (unsigned char *)gs->bufRecv;
    }

    comm_req *req = pwd->req;
    const struct pw_comm_data *c = &pwd->comm[recv];
    const uint *p, *pe, *size = c->size;
    for (p = c->p, pe = p + c->n; p != pe; ++p) {
      const int len = *(size++) * unit_size;
      MPI_CHECK(MPI_Irecv((void *)buf, len, MPI_UNSIGNED_CHAR, *p, *p, gs->comm, req++));
      buf += len;
    }
  }

  {
    unsigned char *buf = (unsigned char *)gs->o_bufSend.ptr();
    if (gs->mode != OOGS_DEVICEMPI) {
      ogs->device.finish(); // waiting for send buffers to be ready
      buf = (unsigned char *)gs->bufSend;
    }

    if (OGS_SYNC_RECV) {
      MPI_Barrier(gs->comm);
    }

    comm_req *req = &pwd->req[pwd->comm[recv].n];
    const struct pw_comm_data *c = &pwd->comm[send];
    const uint *p, *pe, *size = c->size;
    for (p = c->p, pe = p + c->n; p != pe; ++p) {
      const int len = *(size++) * unit_size;
      MPI_CHECK(MPI_Isend((void *)buf, len, MPI_UNSIGNED_CHAR, *p, gs->rank, gs->comm, req++));
      buf += len;
    }
  }

  MPI_CHECK(MPI_Waitall(pwd->comm[send].n + pwd->comm[recv].n, pwd->req, MPI_STATUSES_IGNORE));
}

void occaGatherScatterLocal(const dlong NlocalGather,
                            const dlong NrowBlocks,
                            const occa::memory &o_bstart,
                            const occa::memory &o_gstart,
                            const occa::memory &o_gids,
                            const int Nvectors,
                            const dlong stride,
                            const char *type,
                            const char *op,
                            const occa::memory &o_v)
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
  } else if (!strcmp(type, "double") && !strcmp(op, "add")) {
    ogs::gatherScatterNewKernel_doubleAdd(NrowBlocks,
                                          Nentries,
                                          Nvectors,
                                          stride,
                                          o_bstart,
                                          o_gstart,
                                          o_gids,
                                          o_v);
  } else if (!strcmp(type, "double") && !strcmp(op, "min")) {
    ogs::gatherScatterNewKernel_doubleMin(NrowBlocks,
                                          Nentries,
                                          Nvectors,
                                          stride,
                                          o_bstart,
                                          o_gstart,
                                          o_gids,
                                          o_v);
  } else if (!strcmp(type, "double") && !strcmp(op, "max")) {
    ogs::gatherScatterNewKernel_doubleMax(NrowBlocks,
                                          Nentries,
                                          Nvectors,
                                          stride,
                                          o_bstart,
                                          o_gstart,
                                          o_gids,
                                          o_v);
  } else {
    printf("occaGatherScatterNewKernel: unsupported operation or datatype!\n");
    exit(1);
  }
#endif
}

void oogs::gpu_mpi(int val)
{
  OGS_MPI_SUPPORT = val;
}

void oogs::overlap(int val)
{
  OGS_OVERLAP = val;
}

int oogs::gpu_mpi()
{
  return OGS_MPI_SUPPORT;
}

void oogs::sync_recv(int val)
{
  OGS_SYNC_RECV = val;  
}

void oogs::registerKernels()
{
  ogs::defaultStream = platform->device.occaDevice().getStream();
  ogs::dataStream = platform->device.occaDevice().createStream();

  occa::properties& props = ogs::kernelInfo;

  static bool firstTime = true;
  if (firstTime) {
    props["defines"].asObject();
    props["includes"].asArray();
    props["header"].asArray();
    props["flags"].asObject();
 
    props["defines/ "
          "p_blockSize"] = BLOCKSIZE;
    props["defines/ "
          "dlong"] = dlongString;
    props["defines/ "
          "hlong"] = hlongString;
 
    if ("OpenCL" == platform->device.mode())
      props["defines/"
            "hlong"] = "long";

    props["includes"] += std::string(getenv("NEKRS_KERNEL_DIR")) + "/core/ogs/ogsDefs.h";

    props["defines/init_"
          "float"
          "_add"] = (float)0;
    props["defines/init_"
          "float"
          "_mul"] = (float)1;
    props["defines/init_"
          "float"
          "_min"] = (float)std::numeric_limits<float>::max() / 100;
    props["defines/init_"
          "float"
          "_max"] = (float)std::numeric_limits<float>::lowest() / 100;
 
    props["defines/init_"
          "double"
          "_add"] = (double)0;
    props["defines/init_"
          "double"
          "_mul"] = (double)1;
    props["defines/init_"
          "double"
          "_min"] = (double)std::numeric_limits<double>::max() / 100;
    props["defines/init_"
           "double"
           "_max"] = (double)std::numeric_limits<double>::lowest() / 100;
 
    props["defines/init_"
          "int"
          "_add"] = (int)0;
    props["defines/init_"
          "int"
          "_mul"] = (int)1;
    props["defines/init_"
          "int"
          "_min"] = (int)std::numeric_limits<int>::max() / 100;
    props["defines/init_"
           "int"
           "_max"] = (int)std::numeric_limits<int>::lowest() / 100;
 
    props["defines/init_"
          "long_long"
          "_add"] = static_cast<int64_t>(0);
    props["defines/init_"
          "long_long"
          "_mul"] = static_cast<int64_t>(1);
    props["defines/init_"
          "long_long"
          "_min"] = std::numeric_limits<int64_t>::max() / 100;
    props["defines/init_"
           "long_long"
           "_max"] = std::numeric_limits<int64_t>::lowest() / 100;
 
    props["defines/p_gatherNodesPerBlock"] = ogs::gatherNodesPerBlock;
    firstTime = false;
  }

  ogs::initKernels();

  auto buildKernel = [&props](const std::string& kernelName)
  {
    const auto oklpath = std::string(getenv("NEKRS_KERNEL_DIR")) + "/core/ogs/";
    const auto fileName = oklpath + "oogs.okl";
    const auto reqName = fileName;

    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      platform->kernelRequests.add(reqName, fileName, props); 
      return occa::kernel();
    } else {
      return platform->kernelRequests.load(reqName, kernelName);
    }
  };

  packBufFloatAddKernel = buildKernel("packBuf_float_add");
  unpackBufFloatAddKernel = buildKernel("unpackBuf_float_add");
  packBufDoubleAddKernel = buildKernel("packBuf_double_add");
  unpackBufDoubleAddKernel = buildKernel("unpackBuf_double_add");

  packBufFloatMinKernel = buildKernel("packBuf_float_min");
  unpackBufFloatMinKernel = buildKernel("unpackBuf_float_min");
  packBufDoubleMinKernel = buildKernel("packBuf_double_min");
  unpackBufDoubleMinKernel = buildKernel("unpackBuf_double_min");

  packBufFloatMaxKernel = buildKernel("packBuf_float_max");
  unpackBufFloatMaxKernel = buildKernel("unpackBuf_float_max");
  packBufDoubleMaxKernel = buildKernel("packBuf_double_max");
  unpackBufDoubleMaxKernel = buildKernel("unpackBuf_double_max");
}

void reallocBuffers(int unit_size, oogs_t *gs)
{
  ogs_t *ogs = gs->ogs;
  struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
  const void *execdata = hgs->r.data;
  const struct pw_data *pwd = (pw_data *)execdata;

  if (gs->o_bufSend.size() < pwd->comm[send].total * unit_size) {
    if (gs->o_bufSend.size()) {
      gs->o_bufSend.free();
    }
    if (gs->h_buffSend.size()) {
      gs->h_buffSend.free();
    }
    gs->bufSend = (unsigned char *)ogsHostMallocPinned(ogs->device,
                                                       pwd->comm[send].total * unit_size,
                                                       NULL,
                                                       gs->o_bufSend,
                                                       gs->h_buffSend);
  }
  if (gs->o_bufRecv.size() < pwd->comm[recv].total * unit_size) {
    if (gs->o_bufRecv.size()) {
      gs->o_bufRecv.free();
    }
    if (gs->h_buffRecv.size()) {
      gs->h_buffRecv.free();
    }
    gs->bufRecv = (unsigned char *)ogsHostMallocPinned(ogs->device,
                                                       pwd->comm[recv].total * unit_size,
                                                       NULL,
                                                       gs->o_bufRecv,
                                                       gs->h_buffRecv);
  }
}

size_t typeBytes(const char *type)
{
  if (!strcmp(type, "float")) {
    return sizeof(float);
  } else if (!strcmp(type, "double")) {
    return sizeof(double);
  } else if (!strcmp(type, "int")) {
    return sizeof(int);
  } else if (!strcmp(type, "long long int")) {
    return sizeof(long long int);
  } else {
    return -1;
  }
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
  const std::string oklpath = std::string(getenv("NEKRS_KERNEL_DIR")) + "/core/ogs/";

  struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
  const void *execdata = hgs->r.data;
  const struct pw_data *pwd = (pw_data *)execdata;

  gs->comm = hgs->comm.c;

  int rank;
  MPI_Comm_rank(gs->comm, &rank);
  gs->rank = rank;
  gs->mode = gsMode;

  if (gsMode == OOGS_DEFAULT) {
    return gs;
  }

  int nbcDeviceEnabled = 0;
  if (getenv("OOGS_ENABLE_NBC_DEVICE")) {
    if (std::stoi(getenv("OOGS_ENABLE_NBC_DEVICE")) > 0) {
      nbcDeviceEnabled = 1;
    }
  }

  std::list<oogs_mode> oogs_mode_list;
  oogs_mode_list.push_back(OOGS_LOCAL);
  std::list<oogs_modeExchange> oogs_modeExchange_list;
  oogs_modeExchange_list.push_back(OOGS_EX_PW);

  {
    const hlong NhaloGather = ogs->NhaloGather;
    MPI_Allreduce(&NhaloGather, &(ogs->NhaloGatherGlobal), 1, MPI_HLONG, MPI_SUM, gs->comm);
  }

  if (ogs->NhaloGatherGlobal > 0) {
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
      if (OGS_MPI_SUPPORT) {
        oogs_mode_list.push_back(OOGS_DEVICEMPI);
      }
    }
    oogs_modeExchange_list.push_back(OOGS_EX_NBC);
  }

  const auto ogsModeEnv = (getenv("OOGS_MODE")) ? std::string(getenv("OOGS_MODE")) : "";
  if (!ogsModeEnv.empty() && ogsModeEnv != "OOGS_AUTO") {
    oogs_mode_list.clear();
    oogs_mode_list.push_back(OOGS_LOCAL);
    oogs_modeExchange_list.push_back(OOGS_EX_PW);

    int err = 0;
    if (ogsModeEnv == "OOGS_DEFAULT") {
      oogs_mode_list.push_back(OOGS_DEFAULT);

    } else if (ogsModeEnv.find("OOGS_HOSTMPI") != std::string::npos) {
      oogs_mode_list.push_back(OOGS_HOSTMPI);

      oogs_modeExchange_list.clear();
      if (ogsModeEnv == "OOGS_HOSTMPI+OOGS_EX_PW") {
        oogs_modeExchange_list.push_back(OOGS_EX_PW);
      } else if (ogsModeEnv == "OOGS_HOSTMPI+OOGS_EX_NBC") {
        oogs_modeExchange_list.push_back(OOGS_EX_NBC);
      } else {
        err++;
      }
    } else if (ogsModeEnv.find("OOGS_DEVICEMPI") != std::string::npos) {
      oogs_mode_list.push_back(OOGS_DEVICEMPI);

      oogs_modeExchange_list.clear();
      if (ogsModeEnv == "OOGS_DEVICEMPI+OOGS_EX_PW") {
        oogs_modeExchange_list.push_back(OOGS_EX_PW);
      } else if (ogsModeEnv == "OOGS_DEVICEMPI+OOGS_EX_NBC") {
        oogs_modeExchange_list.push_back(OOGS_EX_NBC);
      } else {
        err++;
      }
    } else {
      err++;
    }

    if (err) {
      printf("OOGS_MODE set to invalid value %s!\n", ogsModeEnv.c_str());
      exit(1);
    }
  }

  if (gsMode == OOGS_AUTO) {
    auto knlOverlapStr = (callback) ? "userKnlOverlap" : ""; 
    if (gs->rank == 0) {
      printf("autotuning gs for wordSize=%d nFields=%d %s\n", (int)Nbytes, nVec, knlOverlapStr);
    }

    double elapsedMin = std::numeric_limits<double>::max();
    oogs_mode fastestMode = OOGS_DEFAULT;
    oogs_modeExchange fastestModeExchange = OOGS_EX_PW;
    int fastestPrepostRecv = 0;

    void *q = calloc(std::max(stride, ogs->N) * nVec * Nbytes, 1);
    occa::memory o_q = device.malloc(std::max(stride, ogs->N) * nVec * Nbytes, q);
    free(q);

    for (auto const &mode : oogs_mode_list) {
      gs->mode = mode;

      for (auto const &modeExchange : oogs_modeExchange_list) {
        gs->modeExchange = modeExchange;

        MPI_Barrier(gs->comm);

        // skip invalid combinations
        if (gs->mode == OOGS_LOCAL && gs->modeExchange != OOGS_EX_PW) {
          continue;
        }
        if (gs->mode == OOGS_DEFAULT && gs->modeExchange != OOGS_EX_PW) {
          continue;
        }
        if (gs->modeExchange == OOGS_EX_NBC && gs->mode == OOGS_DEVICEMPI) {
          if (!nbcDeviceEnabled) {
            continue; // not supported yet by all MPI implementations
          }
        }

        if ((gs->mode == OOGS_DEVICEMPI || gs->mode == OOGS_HOSTMPI) && 
            gs->modeExchange == OOGS_EX_PW) {
          gs->earlyPrepostRecv = 1;
        } else {
          gs->earlyPrepostRecv = 0;
        }

        if (gs->rank == 0) {
          if (gs->mode == OOGS_LOCAL) {
            printf("local:");
          }
          if (gs->mode == OOGS_DEFAULT) {
            printf("pack/unpack host + hostBuffer MPI using pw:");
          }

          const auto exchangeMethod = (gs->modeExchange == OOGS_EX_NBC) ? "nbc" : "pw";
          if (gs->mode == OOGS_HOSTMPI) {
            printf("pack/unpack device + hostBuffer MPI using %s:", exchangeMethod);
          }
          if (gs->mode == OOGS_DEVICEMPI) {
            printf("pack/unpack device + deviceBuffer MPI using %s:", exchangeMethod);
          }
          fflush(stdout);
        }

        // run Ntests measurements and take min to eliminate runtime variations
        constexpr int Ntests = 100;
        double elapsedTest = std::numeric_limits<double>::max();
        for (int test = 0; test < Ntests; ++test) {
          device.finish();
          MPI_Barrier(gs->comm);
          const double tStart = MPI_Wtime();

          oogs::start(o_q, nVec, stride, type, ogsAdd, gs);
          if (callback) {
            callback();
          }
          oogs::finish(o_q, nVec, stride, type, ogsAdd, gs);

          device.finish();
          elapsedTest = std::min(elapsedTest, MPI_Wtime() - tStart);
        }
        MPI_Allreduce(MPI_IN_PLACE, &elapsedTest, 1, MPI_DOUBLE, MPI_MAX, gs->comm);

        if (gs->rank == 0) {
          printf(" %.4es ", elapsedTest);
        }

        if (gs->mode == OOGS_LOCAL && !callback) {
          double rowSizeSum = 0;
          for (dlong i = 0; i < ogs->NlocalGather; i++) {
            rowSizeSum += ogs->localGatherOffsets[i + 1] - ogs->localGatherOffsets[i];
          }
          double localGsBw = (2 * nVec * Nbytes) * rowSizeSum;
          localGsBw += 2 * rowSizeSum * sizeof(int); // index
          MPI_Allreduce(MPI_IN_PLACE, &localGsBw, 1, MPI_DOUBLE, MPI_MIN, gs->comm);
          localGsBw /= elapsedTest;

          int commSize;
          MPI_Comm_size(gs->comm, &commSize);
          if (gs->rank == 0) {
            printf("(%.1fGB/s)", localGsBw / 1e9);
            fflush(stdout);

          }
        }

        if (gs->rank == 0) {
          printf("\n");
          fflush(stdout);
        }

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
    MPI_Bcast(&fastestMode, 1, MPI_INT, 0, gs->comm);
    MPI_Bcast(&fastestModeExchange, 1, MPI_INT, 0, gs->comm);
    MPI_Bcast(&fastestPrepostRecv, 1, MPI_INT, 0, gs->comm);

    gs->mode = fastestMode;
    gs->modeExchange = fastestModeExchange;
    gs->earlyPrepostRecv = fastestPrepostRecv;
    o_q.free();
  } else {
    gs->mode = gsMode;
    gs->modeExchange = OOGS_EX_PW;
    if ((gs->mode == OOGS_DEVICEMPI || gs->mode == OOGS_HOSTMPI)) { 
       gs->earlyPrepostRecv = 1;
    } else {
      gs->earlyPrepostRecv = 0;
    }
  }

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

      if (gs->modeExchange == OOGS_EX_NBC) {
        neighborAllToAll(unit_size, gs);
      } else {
        pairwiseExchange(unit_size, gs);
      }
      elapsedMinMPI = std::min(elapsedMinMPI, MPI_Wtime() - tStart);
    }
    gs->earlyPrepostRecv = earlyPrepostRecv;

    int size;
    MPI_Comm_size(gs->comm, &size);
    double nBytesExchange = unit_size * (pwd->comm[send].total + pwd->comm[recv].total);
    MPI_Allreduce(MPI_IN_PLACE, &nBytesExchange, 1, MPI_DOUBLE, MPI_SUM, gs->comm);
    nBytesExchange /= size;

    double tmin, tmax, tavg;
    MPI_Allreduce(&elapsedMinMPI, &tmin, 1, MPI_DOUBLE, MPI_MIN, gs->comm);
    MPI_Allreduce(&elapsedMinMPI, &tmax, 1, MPI_DOUBLE, MPI_MAX, gs->comm);
    MPI_Allreduce(&elapsedMinMPI, &tavg, 1, MPI_DOUBLE, MPI_SUM, gs->comm);
    tavg /= size;

    if (tmin > MPI_Wtick() && ogs->NhaloGatherGlobal > 0) {
      if (gs->rank == 0 ) {
        printf("MPI min/max/avg: %.2es %.2es %.2es / avg bi-bw: %.1fGB/s/rank\n",
               tmin,
               tmax,
               tavg,
               nBytesExchange / tmax / 1e9);
      }
    }
  }

  fflush(stdout);
  return gs;
}

oogs_t *oogs::setup(dlong N,
                    hlong *ids,
                    int nVec,
                    dlong stride,
                    const char *type,
                    const MPI_Comm &comm,
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
                    const occa::memory &o_gstarts,
                    const occa::memory &o_gids,
                    const occa::memory &o_sstarts,
                    const occa::memory &o_sids,
                    const char *type,
                    const char *op,
                    const occa::memory &o_v,
                    const occa::memory &o_gv)
{
  if (Ngather == 0) {
    return;
  }

  occa::kernel kernel;

  if (!strcmp(type, "float") && !strcmp(op, ogsAdd)) {
    kernel = oogs::packBufFloatAddKernel;
  } else if (!strcmp(type, "float") && !strcmp(op, ogsMin)) {
    kernel = oogs::packBufFloatMinKernel;
  } else if (!strcmp(type, "float") && !strcmp(op, ogsMax)) {
    kernel = oogs::packBufFloatMaxKernel;
  } else if (!strcmp(type, "double") && !strcmp(op, ogsAdd)) {
    kernel = oogs::packBufDoubleAddKernel;
  } else if (!strcmp(type, "double") && !strcmp(op, ogsMin)) {
    kernel = oogs::packBufDoubleMinKernel;
  } else if (!strcmp(type, "double") && !strcmp(op, ogsMax)) {
    kernel = oogs::packBufDoubleMaxKernel;
  } else {
    printf("oogs: unsupported operation %s or datatype %s!\n", op, type);
    exit(1);
  }

  kernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
}

static void unpackBuf(oogs_t *gs,
                      const dlong Ngather,
                      const int k,
                      const dlong stride,
                      const occa::memory &o_gstarts,
                      const occa::memory &o_gids,
                      const occa::memory &o_sstarts,
                      const occa::memory &o_sids,
                      const char *type,
                      const char *op,
                      const occa::memory &o_v,
                      const occa::memory &o_gv)
{
  if (Ngather == 0) {
    return;
  }

  occa::kernel kernel;

  if (!strcmp(type, "float") && !strcmp(op, ogsAdd)) {
    kernel = oogs::unpackBufFloatAddKernel;
  } else if (!strcmp(type, "float") && !strcmp(op, ogsMin)) {
    kernel = oogs::unpackBufFloatMinKernel;
  } else if (!strcmp(type, "float") && !strcmp(op, ogsMax)) {
    kernel = oogs::unpackBufFloatMaxKernel;
  } else if (!strcmp(type, "double") && !strcmp(op, ogsAdd)) {
    kernel = oogs::unpackBufDoubleAddKernel;
  } else if (!strcmp(type, "double") && !strcmp(op, ogsMin)) {
    kernel = oogs::unpackBufDoubleMinKernel;
  } else if (!strcmp(type, "double") && !strcmp(op, ogsMax)) {
    kernel = oogs::unpackBufDoubleMaxKernel;
  } else {
    printf("oogs: unsupported operation %s or datatype %s!\n", op, type);
    exit(1);
  }

  kernel(Ngather, k, stride, o_gstarts, o_gids, o_sstarts, o_sids, o_v, o_gv);
}

void oogs::start(occa::memory &o_v,
                 const int k,
                 const dlong stride,
                 const char *type,
                 const char *op,
                 oogs_t *gs)
{
  ogs_t *ogs = gs->ogs;
  const int unit_size = (int)typeBytes(type) * k;

  if (gs->mode == OOGS_DEFAULT) {
    if (k > 1) {
      ogsGatherScatterManyStart(o_v, k, stride, type, op, ogs);
    } else {
      ogsGatherScatterStart(o_v, type, op, ogs);
    }

    return;
  }

  if (gs->mode != OOGS_LOCAL) {
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

    ogs->device.finish(); // buffers (send/recv) ready for MPI

    if (gs->earlyPrepostRecv) {
      unsigned char *buf = (unsigned char *)gs->o_bufRecv.ptr();
      if (gs->mode != OOGS_DEVICEMPI) {
        buf = (unsigned char *)gs->bufRecv;
      }

      struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
      const void *execdata = hgs->r.data;
      const struct pw_data *pwd = (pw_data *)execdata;

      comm_req *req = pwd->req;
      const struct pw_comm_data *c = &pwd->comm[recv];
      const uint *p, *pe, *size = c->size;
      for (p = c->p, pe = p + c->n; p != pe; ++p) {
        const int len = *(size++) * unit_size;
        MPI_CHECK(MPI_Irecv((void *)buf, len, MPI_UNSIGNED_CHAR, *p, *p, gs->comm, req++));
        buf += len;
      }
    }
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
  const int unit_size = (int)typeBytes(type) * k;

  if (gs->mode == OOGS_DEFAULT) {
    if (k > 1) {
      ogsGatherScatterManyFinish(o_v, k, stride, type, op, ogs);
    } else {
      ogsGatherScatterFinish(o_v, type, op, ogs);
    }

    return;
  }

  if (ogs->NlocalGather) {
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
  }

  if (ogs->NhaloGatherGlobal && !OGS_OVERLAP) {
    ogs->device.finish();
  }

  if (gs->mode == OOGS_HOSTMPI) {
    ogs->device.setStream(ogs::dataStream);

    struct gs_data *hgs = (gs_data *)ogs->haloGshSym;
    const void *execdata = hgs->r.data;
    const struct pw_data *pwd = (pw_data *)execdata;

    if (pwd->comm[send].total) {
      gs->o_bufSend.copyTo(gs->bufSend, pwd->comm[send].total * unit_size, 0, "async: true");
    }

    ogsHostTic(gs->comm, 1);
    if (gs->modeExchange == OOGS_EX_NBC) {
      neighborAllToAll(unit_size, gs);
    } else {
      pairwiseExchange(unit_size, gs);
    }
    ogsHostToc();

    if (pwd->comm[recv].total) {
      gs->o_bufRecv.copyFrom(gs->bufRecv, pwd->comm[recv].total * unit_size, 0, "async: true");
    }

    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);
  }

  if (gs->mode == OOGS_DEVICEMPI) {
    ogsHostTic(gs->comm, 1);
    if (gs->modeExchange == OOGS_EX_NBC) {
      neighborAllToAll(unit_size, gs);
    } else {
      pairwiseExchange(unit_size, gs);
    }
    ogsHostToc();
  }

  if (gs->mode != OOGS_LOCAL)
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

  free(gs);
}
