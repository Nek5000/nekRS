/*
 * TODO:
 *  - Support other operations than just add 
 */

#include <limits>
#include <occa.hpp>
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "ogsInterface.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "gslib.h"

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
                         int *restrict ids){

  uint i,j; 
  int n=0, s=0;
  while((i=*map++)!=UINT_MAX) {
    starts[s] = n;
    j=*map++; 
    do {
      ids[n] = j;
      n++;
    } while((j=*map++)!=UINT_MAX);
    starts[s+1] = n;
    s++;
  }
}

static void _ogsHostGatherScatter(occa::memory o_u,
                                  const char *type, const char *op, 
                                  oogs_t *gs)
{
  ogs_t *ogs = gs->ogs;
  struct gs_data *hgs = (gs_data*) ogs->haloGshSym;
  const void* execdata = hgs->r.data; 
  const struct pw_data *pwd = (pw_data*) execdata; 
  const struct comm *comm = &hgs->comm;
  const unsigned Nhalo = ogs->NhaloGather;

  // hardwired for now
  const unsigned transpose = 0;
  const unsigned recv = 0^transpose, send = 1^transpose;

  size_t unit_size;
  if (!strcmp(type, "float"))
    unit_size = sizeof(float);
  else if (!strcmp(type, "double"))
    unit_size  = sizeof(double);
  else if (!strcmp(type, "int"))
    unit_size  = sizeof(int);
  else if (!strcmp(type, "long long int"))
    unit_size  = sizeof(long long int);

  { // prepost recv
    comm_req *req = pwd->req; 
    const struct pw_comm_data *c = &pwd->comm[recv];
    const uint *p, *pe, *size=c->size;
    uint bufOffset = 0;
    for(p=c->p,pe=p+c->n;p!=pe;++p) {
      size_t len = *(size++)*unit_size;
      unsigned char *recvbuf = (unsigned char *)gs->bufRecv + bufOffset;
      if(gs->mode == OOGS_DEVICEMPI) recvbuf = (unsigned char*)gs->o_bufRecv.ptr() + bufOffset;
      MPI_Irecv((void*)recvbuf,len,MPI_UNSIGNED_CHAR,*p,*p,comm->c,req++);
      bufOffset += len;
    }
  }

  { // scatter
    occaScatter(Nhalo, gs->o_scatterOffsets, gs->o_scatterIds, type, op, o_u, gs->o_bufSend);
    if(gs->mode == OOGS_HOSTMPI) {
      gs->o_bufSend.copyTo(gs->bufSend, pwd->comm[send].total*unit_size, 0, "async: true");
    }
  }

  { // pw exchange
    ogs->device.finish(); // waiting for buffers to be ready
    MPI_Barrier(comm->c);

    comm_req *req = &pwd->req[pwd->comm[recv].n]; 
    const struct pw_comm_data *c = &pwd->comm[send];
    const uint *p, *pe, *size=c->size;
    uint bufOffset = 0;
    for(p=c->p,pe=p+c->n;p!=pe;++p) {
      size_t len = *(size++)*unit_size;
      unsigned char *sendbuf = (unsigned char*)gs->bufSend + bufOffset;
      if(gs->mode == OOGS_DEVICEMPI) sendbuf = (unsigned char*)gs->o_bufSend.ptr() + bufOffset;
      MPI_Isend((void*)sendbuf,len,MPI_UNSIGNED_CHAR,*p,comm->id,comm->c,req++);
      bufOffset += len;
    }
    MPI_Waitall(pwd->comm[send].n + pwd->comm[recv].n,pwd->req,MPI_STATUSES_IGNORE);
  }

  { // gather
    if(gs->mode == OOGS_HOSTMPI){
      gs->o_bufRecv.copyFrom(gs->bufRecv,pwd->comm[recv].total*unit_size, 0, "async: true");
    }

    // op hardwired for now!!!
    occaGather(Nhalo, gs->o_gatherOffsets, gs->o_gatherIds, type, "add+self", gs->o_bufRecv, o_u);
  }

}

oogs_t* oogs::setup(dlong N, hlong *ids, const char *type, MPI_Comm &comm,
                    int verbose, occa::device device, oogs_mode gsMode)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  oogs_t *gs = new oogs_t[1];
  gs->ogs = ogsSetup(N, ids, comm, verbose, device); 
  ogs_t *ogs = gs->ogs;

  const unsigned transpose = 0;
  struct gs_data *hgs = (gs_data*) ogs->haloGshSym;
  const unsigned recv = 0^transpose, send = 1^transpose;
  const void* execdata = hgs->r.data;
  const struct pw_data *pwd = (pw_data*) execdata;
  const unsigned Nhalo = ogs->NhaloGather;
  const unsigned unit_size = sizeof(double); // hardwire just need to be big enough

  if(Nhalo == 0) return gs;

  occa::properties props;
  props["mapped"] = true;


  gs->h_buffSend = ogs->device.malloc(pwd->comm[send].total*unit_size, props);
  gs->bufSend = (unsigned char*)gs->h_buffSend.ptr(props); 
  int *scatterOffsets = (int*) calloc(2*Nhalo,sizeof(int));
  int *scatterIds = (int*) calloc(pwd->comm[send].total,sizeof(int));
  convertPwMap(pwd->map[send], scatterOffsets, scatterIds);

  gs->o_bufSend = ogs->device.malloc(pwd->comm[send].total*unit_size);
  gs->o_scatterOffsets = ogs->device.malloc(2*Nhalo*sizeof(int), scatterOffsets);
  gs->o_scatterIds = ogs->device.malloc(pwd->comm[send].total*sizeof(int), scatterIds);
  free(scatterOffsets);
  free(scatterIds);

  gs->h_buffRecv = ogs->device.malloc(pwd->comm[recv].total*unit_size, props);
  gs->bufRecv = (unsigned char*)gs->h_buffRecv.ptr(props);
  int* gatherOffsets  = (int*) calloc(2*Nhalo,sizeof(int));
  int *gatherIds  = (int*) calloc(pwd->comm[recv].total,sizeof(int));
  convertPwMap(pwd->map[recv], gatherOffsets, gatherIds);

  gs->o_bufRecv = ogs->device.malloc(pwd->comm[recv].total*unit_size);
  gs->o_gatherOffsets  = ogs->device.malloc(2*Nhalo*sizeof(int), gatherOffsets);
  gs->o_gatherIds  = ogs->device.malloc(pwd->comm[recv].total*sizeof(int), gatherIds);
  free(gatherOffsets);
  free(gatherIds);

  if(gsMode == OOGS_AUTO) {
    if(rank == 0) printf("timing gs modes: ");
    const int Ntests = 10;
    double elapsedLast = std::numeric_limits<double>::max();
    oogs_mode fastestMode; 
    occa::memory o_q = device.malloc(N*unit_size);
    for (auto const& mode : {OOGS_DEFAULT, OOGS_HOSTMPI})
    {
      gs->mode = mode;
      device.finish();
      MPI_Barrier(comm);
      const double tStart = MPI_Wtime();
      for(int test=0;test<Ntests;++test) {
        oogs::start (o_q, type, ogsAdd, gs);
        oogs::finish(o_q, type, ogsAdd, gs);
      }
      device.finish();
      MPI_Barrier(comm);
      const double elapsed = (MPI_Wtime() - tStart)/Ntests;
      if(rank == 0) printf("%gs ", elapsed);
      if(elapsed < elapsedLast) fastestMode = gs->mode;
      elapsedLast = elapsed;
    }
    MPI_Bcast(&fastestMode, 1, MPI_INT, 0, comm);
    gs->mode = fastestMode;
    o_q.free();
  } else {
    gs->mode = gsMode;
  }
  if(rank == 0) printf("\nused mode: %d\n", gs->mode);

  return gs; 
}

void oogs::start(occa::memory o_v, const char *type, const char *op, oogs_t *gs) 
{
  size_t Nbytes;
  if (!strcmp(type, "float")) 
    Nbytes = sizeof(float);
  else if (!strcmp(type, "double")) 
    Nbytes = sizeof(double);
  else if (!strcmp(type, "int")) 
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int")) 
    Nbytes = sizeof(long long int);

  ogs_t *ogs = gs->ogs; 

  if (ogs->NhaloGather) {
    if (ogs::o_haloBuf.size() < ogs->NhaloGather*Nbytes) {
      if (ogs::o_haloBuf.size()) ogs::o_haloBuf.free();
      ogs::haloBuf = ogsHostMallocPinned(ogs->device, ogs->NhaloGather*Nbytes, NULL, ogs::o_haloBuf, ogs::h_haloBuf);
    }
  }

  if (ogs->NhaloGather) {
    occaGather(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, o_v, ogs::o_haloBuf);
    ogs->device.finish(); // just in case dataStream is non-blocking 

    if(gs->mode == OOGS_DEFAULT) {
      ogs->device.setStream(ogs::dataStream);
      ogs::o_haloBuf.copyTo(ogs::haloBuf, ogs->NhaloGather*Nbytes, 0, "async: true");
      ogs->device.setStream(ogs::defaultStream);
    }
  }
}

void oogs::finish(occa::memory o_v, const char *type, const char *op, oogs_t *gs) 
{
  size_t Nbytes;
  if (!strcmp(type, "float")) 
    Nbytes = sizeof(float);
  else if (!strcmp(type, "double")) 
    Nbytes = sizeof(double);
  else if (!strcmp(type, "int")) 
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int")) 
    Nbytes = sizeof(long long int);

  ogs_t *ogs = gs->ogs; 

  if(ogs->NlocalGather) {
    occaGatherScatter(ogs->NlocalGather, ogs->o_localGatherOffsets, ogs->o_localGatherIds, type, op, o_v);
  }

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);

#ifdef OGS_ENABLE_TIMER
  timer::tic("gsMPI",1);
#endif
    if(gs->mode == OOGS_DEFAULT) {
      ogs->device.finish(); // waiting for gs::haloBuf copy to finish 
      ogsHostGatherScatter(ogs::haloBuf, type, op, ogs->haloGshSym);
    } else {
      _ogsHostGatherScatter(ogs::o_haloBuf, type, op, gs);
    }   
#ifdef OGS_ENABLE_TIMER
  timer::toc("gsMPI");
#endif
 
    if(gs->mode == OOGS_DEFAULT) { 
      ogs::o_haloBuf.copyFrom(ogs::haloBuf, ogs->NhaloGather*Nbytes, 0, "async: true");
    }

    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);

    occaScatter(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, ogs::o_haloBuf, o_v);
  }
}


void oogs::gatherScatter(occa::memory o_v, const char *type, const char *op, oogs_t *gs){
  oogs::start(o_v, type, op, gs);
  oogs::finish(o_v, type, op, gs);
}

void oogs::gatherScatter(void *v, const char *type, const char *op, oogs_t *gs){
  ogsHostGatherScatter(v, type, op, gs->ogs->hostGsh);
} 
