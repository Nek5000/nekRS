
#include <cstdlib>
#include <mpi.h>
#include "nrs.hpp"
#include "platform.hpp"
#include <vector>

#include "findpts.hpp"

#include "pointInterpolation.hpp"
#include <algorithm>
#include <inttypes.h>

pointInterpolation_t::pointInterpolation_t(nrs_t *nrs_, double bb_tol, double newton_tol_, bool mySession_)
    : pointInterpolation_t(nrs_, platform->comm.mpiCommParent, nrs_->_mesh->Nlocal, nrs_->_mesh->Nlocal, bb_tol, newton_tol_, mySession_)
{
}

pointInterpolation_t::pointInterpolation_t(nrs_t *nrs_,
                                           MPI_Comm comm,
                                           dlong localHashSize,
                                           dlong globalHashSize,
                                           double bb_tol,
                                           double newton_tol_,
                                           bool mySession_)
    : nrs(nrs_), newton_tol(newton_tol_), mySession(mySession_), nPoints(0)
{

  // communicator is implicitly required to be either platform->comm.mpiComm or platform->comm.mpiCommParent
  // due to other communicator synchronous calls, such as platform->timer.tic
  bool supported = false;
  for(auto && supportedCommunicator : {platform->comm.mpiComm, platform->comm.mpiCommParent}){
    int same = 0;
    MPI_Comm_compare(comm, supportedCommunicator, &same);
    supported |= (same != MPI_UNEQUAL);
  }
  nrsCheck(!supported,
    comm,
    EXIT_FAILURE,
    "%s",
    "Communicator passed to pointInterpolation_t must be either platform->comm.mpiComm or platform->comm.mpiCommParent");

  newton_tol = std::max(5e-13, newton_tol_);

  const int npt_max = 128;

  mesh_t *mesh = nrs->_mesh;

  if (mySession) {
    mesh->o_x.copyTo(mesh->x, mesh->Nlocal * sizeof(dfloat));
    mesh->o_y.copyTo(mesh->y, mesh->Nlocal * sizeof(dfloat));
    mesh->o_z.copyTo(mesh->z, mesh->Nlocal * sizeof(dfloat));
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
                                                  newton_tol);
}

void pointInterpolation_t::find(pointInterpolation_t::VerbosityLevel verbosity)
{
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic("pointInterpolation_t::find", 1);
  }

  int iErr = 0;
  iErr += !pointsAdded;
  MPI_Allreduce(MPI_IN_PLACE, &iErr, 1, MPI_DLONG, MPI_MAX, platform->comm.mpiComm);
  if (iErr) {
    if (platform->comm.mpiRank == 0) {
      std::cout << "pointInterpolation_t::find called without any points added!\n";
    }
  }
  
  nrsCheck(iErr, platform->comm.mpiComm, EXIT_FAILURE, "%s", "");

  const auto n = nPoints;

  if (useHostPoints) {
    findpts_->find(&data_, _x, _y, _z, n);
  }
  else {
    findpts_->find(&data_, _o_x, _o_y, _o_z, n);
  }

  if (verbosity != VerbosityLevel::None) {

    auto *h_x = _x;
    auto *h_y = _y;
    auto *h_z = _z;
    if (useDevicePoints && verbosity == VerbosityLevel::Detailed) {
      h_x = h_x_vec.data();
      h_y = h_y_vec.data();
      h_z = h_z_vec.data();
      _o_x.copyTo(h_x, n * sizeof(dfloat));
      _o_y.copyTo(h_y, n * sizeof(dfloat));
      _o_z.copyTo(h_z, n * sizeof(dfloat));
    }

    dlong nOutside = 0;
    dlong nBoundary = 0;
    for (int in = 0; in < n; ++in) {
      if (data_.code_base[in] == findpts::CODE_BORDER) {
        if (data_.dist2_base[in] > 10 * newton_tol) {
          nBoundary += 1;
          if (nBoundary < 5 && verbosity == VerbosityLevel::Detailed) {
            std::cout << " WARNING: point on boundary or outside the mesh xy[z]d^2: " << h_x[in] << ","
                      << h_y[in] << ", " << h_z[in] << ", " << data_.dist2_base[in] << std::endl;
          }
        }
      }
      else if (data_.code_base[in] == findpts::CODE_NOT_FOUND) {
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
}

void pointInterpolation_t::eval(dlong nFields,
                                dlong inputFieldOffset,
                                occa::memory o_in,
                                dlong outputFieldOffset,
                                occa::memory o_out)
{
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic("pointInterpolation_t::eval", 1);
  }
  
  nrsCheck(nrs->fieldOffset > inputFieldOffset,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "pointInterpolation_t::eval inputFieldOffset (%d) is less than nrs->fieldOffset (%d)\n",
           inputFieldOffset,
           nrs->fieldOffset);
  
  nrsCheck(o_in.size() < nFields * inputFieldOffset * sizeof(dfloat),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "pointInterpolation_t::eval input size (%" PRId64 ") is smaller than expected (%ld)\n",
           o_in.size(),
           nFields * inputFieldOffset * sizeof(dfloat));
  
  nrsCheck(o_out.size() < nFields * outputFieldOffset * sizeof(dfloat),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "pointInterpolation_t::eval output size (%" PRId64 ") is smaller than expected (%ld)\n",
           o_out.size(),
           nFields * outputFieldOffset * sizeof(dfloat));

  findpts_->eval(nPoints, nFields, inputFieldOffset, outputFieldOffset, o_in, &data_, o_out);

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc("pointInterpolation_t::eval");
  }
}

void pointInterpolation_t::eval(dlong nFields,
                                dlong inputFieldOffset,
                                dfloat *in,
                                dlong outputFieldOffset,
                                dfloat *out)
{
  if (timerLevel != TimerLevel::None) {
    platform->timer.tic("pointInterpolation_t::eval", 1);
  }

  nrsCheck(nrs->fieldOffset > inputFieldOffset,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "pointInterpolation_t::eval inputFieldOffset (%d) is less than nrs->fieldOffset (%d)\n",
           inputFieldOffset,
           nrs->fieldOffset);

  findpts_->eval(nPoints, nFields, inputFieldOffset, outputFieldOffset, in, &data_, out);

  if (timerLevel != TimerLevel::None) {
    platform->timer.toc("pointInterpolation_t::eval");
  }
}

void pointInterpolation_t::setPoints(int n, dfloat *x, dfloat *y, dfloat *z)
{

  pointsAdded = true;
  useHostPoints = true;
  useDevicePoints = false;

  if(n > nPoints){
    data_ = findpts::data_t(n);
  }

  nPoints = n;

  _x = x;
  _y = y;
  _z = z;
}

void pointInterpolation_t::setPoints(int n, occa::memory o_x, occa::memory o_y, occa::memory o_z)
{

  pointsAdded = true;
  useHostPoints = false;
  useDevicePoints = true;

  if (n > nPoints) {
    data_ = findpts::data_t(n);
  }

  nPoints = n;

  _o_x = o_x;
  _o_y = o_y;
  _o_z = o_z;

  h_x_vec.resize(n);
  h_y_vec.resize(n);
  h_z_vec.resize(n);
}

void pointInterpolation_t::setTimerLevel(TimerLevel level)
{
  timerLevel = level;
  findpts_->setTimerLevel(level);
}

TimerLevel pointInterpolation_t::getTimerLevel() const { return timerLevel; }

void pointInterpolation_t::setTimerName(std::string name)
{
  timerName = name;
  findpts_->setTimerName(name);
}

void pointInterpolation_t::update() { findpts_->update(data_); }
