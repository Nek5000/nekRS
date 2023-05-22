#if !defined(neknek_hpp_)
#define neknek_hpp_

#include <mpi.h>
#include "nrssys.hpp"
#include "findpts.hpp"
#include "pointInterpolation.hpp"
#include <vector>
#include <memory>

struct nrs_t;

bool neknekCoupled();

class neknek_t {
public:
  neknek_t(nrs_t *nrs, dlong _nsessions, dlong _sessionID);
  dlong nsessions, sessionID;

  dlong Nscalar;
  dlong nEXT;
  std::vector<dfloat> coeffEXT;
  occa::memory o_coeffEXT;
  bool globalMovingMesh;
  dlong npt;

  // page-aligned length >= npt
  dlong fieldOffset;

  // maps from E-vector point to interpolated point
  std::vector<dlong> pointMap;
  occa::memory o_pointMap;

  occa::memory o_U;
  occa::memory o_S;

  occa::memory o_x;
  occa::memory o_y;
  occa::memory o_z;

  // metrics for printing
  dfloat tSync;
  dfloat tExch;
  dfloat ratio;

  std::shared_ptr<pointInterpolation_t> interpolator;
  void updateBoundary(nrs_t *nrs, int tstep, int stage);

  occa::kernel copyNekNekPointsKernel;
};

#endif
