#include "platform.hpp"
#include "QQt.hpp"

QQt::QQt(std::vector<hlong> ids, const timingConfig& config, oogs_mode gsMode, const MPI_Comm& comm)
{
  int verbose = 0;
  gsH = oogs::setup(static_cast<int>(ids.size()),
                    ids.data(),
                    config.nVec,
                    config.stride,
                    config.type.c_str(),
                    comm,
                    verbose,
                    platform->device.occaDevice(),
                    config.callback,
                    gsMode);
}

QQt::QQt(oogs_t *h)

{
  gsH = h;
}

void QQt::startFinish(const std::string& op_,
                      occa::memory& o_v,
                      const dlong stride,
                      const int k)
{
  auto type = [o_v]()
  {
    std::string val;
    if (o_v.dtype().name() == ogsFloat)
      val = ogsFloat;
    else if (o_v.dtype().name() == ogsDouble)
      val = ogsDouble;

    nekrsCheck(val.empty(), MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "invalid datatype"); 
    return val;
  }();

  auto op = [&op_]()
  {
    std::string val;
    if (op_ == "+")
      val = ogsAdd;
    else if (op_ == "*")
      val = ogsMul;
    else if (op_ == "min")
      val = ogsMin;
    else if (op_ == "max")
      val = ogsMax;

    nekrsCheck(val.empty(), MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "invalid operation"); 
    return val;
  }();

  oogs::start(o_v, k, stride, type.c_str(), op.c_str(), gsH);
  oogs::finish(o_v, k, stride, type.c_str(), op.c_str(), gsH);
}

void QQt::startFinish(const std::string& op,
                      occa::memory& o_v,
                      const dlong stride)
{
  startFinish(op, o_v, stride, o_v.size()/stride);
}

void QQt::startFinish(const std::string& op,
                      occa::memory& o_v)
{
  startFinish(op, o_v, 0, 1);
}

