#if !defined(neknek_hpp_)
#define neknek_hpp_

#include <mpi.h>
#include "nrssys.hpp"
#include "ogs_FINDPTS.hpp"

struct neknek_t {
  dlong nsessions, sessionID;
  MPI_Comm global_comm;
  bool connected;

  dlong NcorrectorSteps;

  dlong Nscalar;
  bool globalMovingMesh;
  ogs_findpts_t *ogs_handle = nullptr;
  dlong npt;

  dlong *point_map;
  occa::memory o_point_map;

  dfloat *val_interp = nullptr;
  occa::memory o_val_interp;

  dlong  *code = nullptr;
  dlong  *proc = nullptr;
  dlong  *el = nullptr;
  dfloat *r = nullptr;
  dfloat *eval_buffer = nullptr;
};


struct nrs_t;
void neknek_setup(nrs_t *nrs);
void neknek_update_boundary(nrs_t *nrs);


#endif
