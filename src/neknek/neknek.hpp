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
  neknek_t(nrs_t *_nrs, dlong _nsessions, dlong _sessionID);

  void updateBoundary(int tstep, int stage);

  // provide partition of unity function
  // this may be used, e.g., to do global integrations across the domain
  // when collocated with the mass matrix
  occa::memory partitionOfUnity();

  void fixCoupledSurfaceFlux(occa::memory o_U);

  // performance metrics for printing
  dfloat tSync() const
  {
    return tSync_;
  }

  dfloat tExch() const
  {
    return tExch_;
  }

  dfloat ratio() const
  {
    return ratio_;
  }

  // extrapolation order used in neknek
  dlong nEXT() const
  {
    return nEXT_;
  }

  // number of interpolation points
  dlong npt() const
  {
    return npt_;
  }

  // page-aligned length >= npt
  dlong fieldOffset() const
  {
    return fieldOffset_;
  }

  // number of coupled sessions
  dlong nsessions() const
  {
    return nsessions_;
  }

  // id of the current session
  dlong sessionID() const
  {
    return sessionID_;
  }

  // number of scalar fields
  dlong Nscalar() const
  {
    return Nscalar_;
  }

  const occa::memory &o_x() const
  {
    return o_x_;
  }

  const occa::memory &o_y() const
  {
    return o_y_;
  }

  const occa::memory &o_z() const
  {
    return o_z_;
  }

  const occa::memory &o_U() const
  {
    return o_U_;
  }

  const occa::memory &o_S() const
  {
    return o_S_;
  }

  const occa::memory &o_pointMap() const
  {
    return o_pointMap_;
  }

  const occa::memory &o_session() const
  {
    return o_session_;
  }

  const occa::memory &o_partition() const
  {
    return o_partition_;
  }

private:
  void reserveAllocation();
  void updateInterpPoints();
  void findIntPoints();
  void setup();

  nrs_t *nrs = nullptr;
  std::vector<dfloat> coeffEXT;
  occa::memory o_coeffEXT;
  bool globalMovingMesh;

  occa::memory o_x_;
  occa::memory o_y_;
  occa::memory o_z_;
  occa::memory o_U_;
  occa::memory o_S_;
  occa::memory o_pointMap_;
  occa::memory o_session_;
  occa::memory o_partition_;

  // performance metrics for printing
  dfloat tSync_;
  dfloat tExch_;
  dfloat ratio_;

  // extrapolation order used in neknek
  dlong nEXT_;

  // number of interpolation points
  dlong npt_;

  // page-aligned length >= npt
  dlong fieldOffset_;

  dlong nsessions_;
  dlong sessionID_;

  dlong Nscalar_;

  std::shared_ptr<pointInterpolation_t> interpolator;

  occa::kernel copyNekNekPointsKernel;
  occa::kernel computeFluxKernel;
  occa::kernel fixSurfaceFluxKernel;

  bool recomputePartition = true;
};

#endif
