
#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include "platform.hpp"
#include <vector>

#include "findpts.hpp"

#include "pointInterpolation.hpp"
#include <algorithm>
#include <inttypes.h>

#include "bcMap.hpp"

pointInterpolation_t::pointInterpolation_t(mesh_t *mesh,
                                           double bb_tol,
                                           double newton_tol_,
                                           bool mySession_,
                                           dlong sessionID_,
                                           bool multipleSessionSupport_)
    : pointInterpolation_t(mesh,
                           platform->comm.mpiCommParent,
                           mesh->Nlocal,
                           mesh->Nlocal,
                           bb_tol,
                           newton_tol_,
                           mySession_,
                           sessionID_,
                           multipleSessionSupport_)
{
}

pointInterpolation_t::pointInterpolation_t(mesh_t *mesh_,
                                           MPI_Comm comm,
                                           dlong localHashSize,
                                           dlong globalHashSize,
                                           double bb_tol,
                                           double newton_tol_,
                                           bool mySession_,
                                           dlong sessionID_,
                                           bool multipleSessionSupport_)
    : mesh(mesh_), newton_tol(newton_tol_), mySession(mySession_), sessionID(sessionID_),
      multipleSessionSupport(multipleSessionSupport_), nPoints(0)
{
  // communicator is implicitly required to be either platform->comm.mpiComm or platform->comm.mpiCommParent
  // due to other communicator synchronous calls, such as platform->timer.tic
  bool supported = false;
  for (auto &&supportedCommunicator : {platform->comm.mpiComm, platform->comm.mpiCommParent}) {
    int same = 0;
    MPI_Comm_compare(comm, supportedCommunicator, &same);
    supported |= (same != MPI_UNEQUAL);
  }
  nekrsCheck(!supported,
             comm,
             EXIT_FAILURE,
             "%s",
             "Communicator must be either platform->comm.mpiComm or platform->comm.mpiCommParent");

  newton_tol =
      (sizeof(dfloat) == sizeof(double)) ? std::max(5e-13, newton_tol_) : std::max(1e-6, newton_tol_);

  const int npt_max = 1;

  if (mySession) {
    mesh->o_x.copyTo(mesh->x, mesh->Nlocal);
    mesh->o_y.copyTo(mesh->y, mesh->Nlocal);
    mesh->o_z.copyTo(mesh->z, mesh->Nlocal);
  }

  std::vector<dfloat> dist;
  if (multipleSessionSupport_) {
    std::vector<int> bID;
    for (auto &[key, bcID] : bcMap::map()) {
      const auto field = key.first;
      if (field == "velocity") {
        if (bcID == bcMap::bcTypeINT) {
          bID.push_back(key.second + 1);
        }
      }
    }

    const auto nbID = bID.size();
    auto o_bID = platform->o_memPool.reserve<dlong>(nbID);
    o_bID.copyFrom(bID.data());

    _o_distance = mesh->minDistance(nbID, o_bID, "cheap_dist");
    dist.resize(mesh->Nlocal);
    _o_distance.copyTo(dist.data(), mesh->Nlocal);
    o_bID.free();
  }

  findpts_ = std::make_unique<findpts::findpts_t>(comm,
                                                  mySession ? mesh->x : nullptr,
                                                  mySession ? mesh->y : nullptr,
                                                  mySession ? mesh->z : nullptr,
                                                  mesh->Nq,
                                                  mySession ? mesh->Nelements : 0,
                                                  2 * mesh->Nq,
                                                  bb_tol,
                                                  localHashSize,
                                                  globalHashSize,
                                                  npt_max,
                                                  newton_tol,
                                                  sessionID,
                                                  dist.data());
}

occa::memory pointInterpolation_t::distance()
{
  nekrsCheck(!multipleSessionSupport,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "distance requires multipleSessionSupport to be enabled in the pointInterpolation_t object!");
  return _o_distance;
}

void pointInterpolation_t::find(pointInterpolation_t::VerbosityLevel verbosity, bool matchSession)
{
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic("pointInterpolation_t::find");
  }

  int iErr = 0;
  iErr += !pointsAdded;
  nekrsCheck(iErr, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "find called without any points added!");

  const auto n = nPoints;
  const dlong sessionIDMatch = matchSession;

  if (useHostPoints) {
    findpts_->find(&data_, _x, _y, _z, _session, sessionIDMatch, n);
  } else {
    findpts_->find(&data_, _o_x, _o_y, _o_z, _o_session, sessionIDMatch, n);
  }

  if (verbosity != VerbosityLevel::None) {

    auto *h_x = _x;
    auto *h_y = _y;
    auto *h_z = _z;
    if (useDevicePoints && verbosity == VerbosityLevel::Detailed) {
      h_x = h_x_vec.data();
      h_y = h_y_vec.data();
      h_z = h_z_vec.data();
      _o_x.copyTo(h_x, n);
      _o_y.copyTo(h_y, n);
      _o_z.copyTo(h_z, n);
    }

    dlong nOutside = 0;
    dlong nBoundary = 0;
    for (int in = 0; in < n; ++in) {
      if (data_.code_base[in] == findpts::CODE_BORDER) {
        if (data_.dist2_base[in] > 10 * newton_tol) {
          nBoundary += 1;
          if (nBoundary < 5 && verbosity == VerbosityLevel::Detailed) {
            std::cout << " WARNING: point on boundary or outside the mesh distNorm2: " << h_x[in] << ","
                      << h_y[in] << ", " << h_z[in] << ", " << data_.dist2_base[in] << std::endl;
          }
        }
      } else if (data_.code_base[in] == findpts::CODE_NOT_FOUND) {
        nOutside += 1;
        if (nOutside < 5 && verbosity == VerbosityLevel::Detailed) {
          std::cout << " WARNING: point not within mesh xy[z]: " << h_x[in] << "," << h_y[in] << ", "
                    << h_z[in] << std::endl;
        }
      }
    }
    std::array<hlong, 3> counts = {n, nBoundary, nOutside};
    MPI_Allreduce(MPI_IN_PLACE, counts.data(), counts.size(), MPI_HLONG, MPI_SUM, platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0 && counts[2] > 0) {
      std::cout << "WARNING interp::find - total = " << counts[0] << ", boundary = " << counts[1]
                << ", outside = " << counts[2] << "\n";
    }
  }

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc("pointInterpolation_t::find");
  }

  findCalled = true;
}

void pointInterpolation_t::eval(dlong nFields,
                                dlong inputFieldOffset,
                                const occa::memory& o_in,
                                dlong outputFieldOffset,
                                occa::memory &o_out)
{
  nekrsCheck(!findCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "find has not been called prior to eval!");

  if (timerLevel != TimerLevel::None) {
    platform->timer.tic("pointInterpolation_t::eval");
  }

  nekrsCheck(mesh->Nlocal > inputFieldOffset,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "pointInterpolation_t::eval inputFieldOffset (%d) is less than mesh->Nlocal (%d)\n",
             inputFieldOffset,
             mesh->Nlocal);

  nekrsCheck(o_in.byte_size() < nFields * inputFieldOffset * sizeof(dfloat),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "pointInterpolation_t::eval input size (%" PRId64 ") is smaller than expected (%ld)\n",
             o_in.byte_size(),
             nFields * inputFieldOffset * sizeof(dfloat));

  nekrsCheck(o_out.byte_size() < nFields * outputFieldOffset * sizeof(dfloat),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "pointInterpolation_t::eval output size (%" PRId64 ") is smaller than expected (%ld)\n",
             o_out.byte_size(),
             nFields * outputFieldOffset * sizeof(dfloat));

  findpts_->eval(nPoints, nFields, inputFieldOffset, outputFieldOffset, o_in, &data_, o_out);

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc("pointInterpolation_t::eval");
  }
}

void pointInterpolation_t::eval(dlong nFields,
                                dlong inputFieldOffset,
                                const std::vector<dfloat>& in,
                                dlong outputFieldOffset,
                                std::vector<dfloat>& out)
{
  nekrsCheck(!findCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "find has not been called prior to eval!");

  if (timerLevel != TimerLevel::None) {
    platform->timer.tic("pointInterpolation_t::eval");
  }

  nekrsCheck(mesh->Nlocal > inputFieldOffset,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "pointInterpolation_t::eval inputFieldOffset (%d) is less than mesh->Nlocal (%d)\n",
             inputFieldOffset,
             mesh->Nlocal);

  findpts_->eval(nPoints, nFields, inputFieldOffset, outputFieldOffset, const_cast<dfloat*>(in.data()), &data_, out.data());

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc("pointInterpolation_t::eval");
  }
}

void pointInterpolation_t::setPoints(const std::vector<dfloat>& x, const std::vector<dfloat>& y, const std::vector<dfloat>& z)
{
  std::vector<dlong> session;
  this->setPoints(x, y, z, session);
}

void pointInterpolation_t::setPoints(const std::vector<dfloat>& x, const std::vector<dfloat>& y, const std::vector<dfloat>& z, const std::vector<dlong>& session)
{
  auto o_x = platform->device.malloc<dfloat>(x.size());
  o_x.copyFrom(x.data());
  auto o_y = platform->device.malloc<dfloat>(y.size());
  o_y.copyFrom(y.data());
  auto o_z = platform->device.malloc<dfloat>(z.size());
  o_z.copyFrom(z.data());

  occa::memory o_session;
  if (session.size()) {
    o_session = platform->device.malloc<dlong>(session.size());
  } 
  this->setPoints(o_x, o_y, o_z, o_session);
}

void pointInterpolation_t::setPoints(const occa::memory &o_x,
                                     const occa::memory &o_y,
                                     const occa::memory &o_z)
{
  this->setPoints(o_x, o_y, o_z, o_NULL);
}

void pointInterpolation_t::setPoints(const occa::memory &o_x,
                                     const occa::memory &o_y,
                                     const occa::memory &o_z,
                                     const occa::memory &o_session)
{
  const int n = o_x.size();

  pointsAdded = true;
  useHostPoints = false;
  useDevicePoints = true;

  if (n > nPoints) {
    data_ = findpts::data_t(n);
  }

  if (n > 0) {
    _o_x = o_x;
    _o_y = o_y;
    _o_z = o_z;
    _o_session = o_session;

    h_x_vec.resize(n);
    h_y_vec.resize(n);
    h_z_vec.resize(n);
  }

  nPoints = n;
}

void pointInterpolation_t::setTimerLevel(TimerLevel level)
{
  timerLevel = level;
  findpts_->setTimerLevel(level);
}

TimerLevel pointInterpolation_t::getTimerLevel() const
{
  return timerLevel;
}

void pointInterpolation_t::setTimerName(std::string name)
{
  timerName = name;
  findpts_->setTimerName(name);
}

void pointInterpolation_t::o_update()
{
  findCalled = true;
  findpts_->o_update(data_);
}
