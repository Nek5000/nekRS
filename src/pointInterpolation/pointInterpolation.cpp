#include <inttypes.h>
#include "platform.hpp"
#include "findpts.hpp"
#include "pointInterpolation.hpp"

pointInterpolation_t::pointInterpolation_t(mesh_t *mesh_,
                                           MPI_Comm comm,
                                           bool mySession_,
                                           std::vector<int> bIntID,
                                           double bb_tol,
                                           double newton_tol_,
                                           dlong localHashSize,
                                           dlong globalHashSize)

    : mesh(mesh_), mySession(mySession_), nPoints(0)
{
  if (localHashSize == 0) localHashSize = mesh->Nlocal; 
  if (globalHashSize == 0) globalHashSize = mesh->Nlocal;

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
             "%s\n",
             "Communicator must be either platform->comm.mpiComm or platform->comm.mpiCommParent");

  newton_tol =
    (sizeof(dfloat) == sizeof(double)) ? std::max(5e-13, newton_tol_) : std::max(1e-6, newton_tol_);

  auto x = platform->memoryPool.reserve<dfloat>(mesh->Nlocal);
  auto y = platform->memoryPool.reserve<dfloat>(mesh->Nlocal);
  auto z = platform->memoryPool.reserve<dfloat>(mesh->Nlocal);

  if (mySession) {
    mesh->o_x.copyTo(x, mesh->Nlocal);
    mesh->o_y.copyTo(y, mesh->Nlocal);
    mesh->o_z.copyTo(z, mesh->Nlocal);
  }

  std::vector<dfloat> distanceINT;
  if (bIntID.size()) {
    auto o_bIntID = platform->deviceMemoryPool.reserve<int>(bIntID.size());
    o_bIntID.copyFrom(bIntID.data());
    _o_distanceINT = mesh->minDistance(bIntID.size(), o_bIntID, "cheap_dist");
    distanceINT.resize(mesh->Nlocal);
    _o_distanceINT.copyTo(distanceINT.data(), mesh->Nlocal);
  }

  // number of points to iterate on simultaneously
  const int npt_max = 1;

  int sessionID = 0;
  platform->options.getArgs("NEKNEK SESSION ID", sessionID);

  auto xPtr = x.ptr<dfloat>();
  auto yPtr = y.ptr<dfloat>();
  auto zPtr = z.ptr<dfloat>();

  findpts_ = std::make_unique<findpts::findpts_t>(comm,
                                                  mySession ? xPtr : nullptr,
                                                  mySession ? yPtr : nullptr,
                                                  mySession ? zPtr : nullptr,
                                                  mesh->Nq,
                                                  mySession ? mesh->Nelements : 0,
                                                  2 * mesh->Nq,
                                                  bb_tol,
                                                  localHashSize,
                                                  globalHashSize,
                                                  npt_max,
                                                  newton_tol,
                                                  sessionID,
                                                  distanceINT.data());
}

occa::memory pointInterpolation_t::distanceINT()
{
  nekrsCheck(!_o_distanceINT.isInitialized(),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "No INT boundary IDs provided on setup!");

  return _o_distanceINT;
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

    const auto maxVerbosePoints = 5;

    dlong nOutside = 0;
    dlong nBoundary = 0;
    dfloat maxDistNorm = 0;
    for (int in = 0; in < n; ++in) {
      if (data_.code_base[in] == findpts::CODE_BORDER) {
        if (data_.dist2_base[in] > 10 * newton_tol) {
          const auto distNorm = data_.dist2_base[in];
          maxDistNorm = std::max(maxDistNorm, distNorm);
          nBoundary++;
          if (nBoundary < maxVerbosePoints && verbosity == VerbosityLevel::Detailed) {
            std::cout << "pointInterpolation_t::find: WARNING point on boundary or outside the mesh"
                      << " xyz= " << h_x[in] << " " << h_y[in] << " " << h_z[in]
                      << " distNorm= " << std::scientific << std::setprecision(3) << distNorm << std::endl;
          }
        }
      } else if (data_.code_base[in] == findpts::CODE_NOT_FOUND) {
        nOutside++;
        if (nOutside < maxVerbosePoints && verbosity == VerbosityLevel::Detailed) {
          std::cout << "pointInterpolation_t::find: WARNING point outside the mesh"
                    << " xyz= " << h_x[in] << " " << h_y[in] << " " << h_z[in] << std::endl;
        }
      }
    }

    std::array<hlong, 3> counts = {n, nBoundary, nOutside};
    MPI_Allreduce(MPI_IN_PLACE, counts.data(), counts.size(), MPI_HLONG, MPI_SUM, platform->comm.mpiComm);
    MPI_Allreduce(MPI_IN_PLACE, &maxDistNorm, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);

    if (platform->comm.mpiRank == 0 && verbosity == VerbosityLevel::Detailed) {
      std::cout << "pointInterpolation_t::find:"
                << " total= " << counts[0] << " boundary= " << counts[1] << " (max distNorm=" << maxDistNorm
                << ")"
                << " outside= " << counts[2] << std::endl;
    }
  }

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc("pointInterpolation_t::find");
  }

  findCalled = true;
  data_.updateCache = true;
}

void pointInterpolation_t::eval(dlong nFields,
                                dlong inputFieldOffset,
                                const occa::memory &o_in,
                                dlong outputFieldOffset,
                                occa::memory &o_out,
                                dlong nPointsIn,
                                dlong offset)
{
  if (inputFieldOffset == 0) {
    inputFieldOffset = o_in.size();
  }
  if (outputFieldOffset == 0) {
    outputFieldOffset = o_out.size();
  }

  auto nPoints_ = (nPointsIn > -1) ? nPointsIn : nPoints;
  if (nPointsIn >= 0) {
    data_.updateCache = true; // enforce update as cache cannot be used
  }

  nekrsCheck(!findCalled, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "find has not been called prior to eval!");

  nekrsCheck(nFields > 1 && mesh->Nlocal > inputFieldOffset,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "pointInterpolation_t::eval inputFieldOffset (%d) is less than mesh->Nlocal (%d)\n",
             inputFieldOffset,
             mesh->Nlocal);

  nekrsCheck(nFields > 1 && nPoints_ > outputFieldOffset,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "pointInterpolation_t::eval outputFieldOffset (%d) is less than nPoints (%d)\n",
             inputFieldOffset,
             nPoints_);

  nekrsCheck(o_in.byte_size() < nFields * inputFieldOffset * sizeof(dfloat),
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "pointInterpolation_t::eval input size (%" PRId64 ") is smaller than expected (%ld)\n",
             o_in.byte_size(),
             nFields * inputFieldOffset * sizeof(dfloat));

  nekrsCheck(o_out.byte_size() < nFields * outputFieldOffset * sizeof(dfloat),
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "pointInterpolation_t::eval output size (%" PRId64 ") is smaller than expected (%ld)\n",
             o_out.byte_size(),
             nFields * outputFieldOffset * sizeof(dfloat));

  if (timerLevel != TimerLevel::None) {
    platform->timer.tic("pointInterpolation_t::eval");
  }

  findpts_->eval(nPoints_, offset, nFields, inputFieldOffset, outputFieldOffset, o_in, &data_, o_out);

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc("pointInterpolation_t::eval");
  }
}

void pointInterpolation_t::setPoints(const std::vector<dfloat> &x,
                                     const std::vector<dfloat> &y,
                                     const std::vector<dfloat> &z)
{
  std::vector<dlong> session;
  this->setPoints(x, y, z, session);
}

void pointInterpolation_t::setPoints(const std::vector<dfloat> &x,
                                     const std::vector<dfloat> &y,
                                     const std::vector<dfloat> &z,
                                     const std::vector<dlong> &session)
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
    o_session.copyFrom(session.data());
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
