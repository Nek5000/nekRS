#if !defined(nekrs_interp_hpp_)
#define nekrs_interp_hpp_

#include "nekrsSys.hpp"
#include "platform.hpp"
#include "findpts.hpp"
#include "mesh.h"

using findpts::TimerLevel;

class pointInterpolation_t
{
public:
  enum class VerbosityLevel { None, Basic, Detailed };
  pointInterpolation_t(mesh_t *mesh,
                       MPI_Comm comm,
                       bool mySession_ = true,
                       std::vector<int> bIntID = {});
  pointInterpolation_t(mesh_t *mesh_,
                       MPI_Comm comm,
                       dlong localHashSize,
                       dlong globalHashSize,
                       double bb_tol = 0.01,
                       double newton_tol = 0,
                       bool mySession_ = true,
                       std::vector<int> bIntID = {});
  ~pointInterpolation_t() = default;

  // Finds the process, element, and reference coordinates of the given points
  void find(VerbosityLevel verbosity = VerbosityLevel::Basic, bool matchSession = false);

  void eval(dlong nFields,
            dlong inputFieldOffset,
            const occa::memory& o_in,
            dlong outputFieldOffset,
            occa::memory& o_out);

  void eval(dlong nFields, dlong inputFieldOffset, const std::vector<dfloat>& in, dlong outputFieldOffset, std::vector<dfloat> &out);

  auto *ptr()
  {
    return findpts_.get();
  }

  auto &data()
  {
    return data_;
  }

  int numPoints() const
  {
    return nPoints;
  }

  // calls underlying findpts_t::update
  void o_update();

  void setPoints(const occa::memory &o_x, const occa::memory &o_y, const occa::memory &o_z);
  void setPoints(const occa::memory &o_x, const occa::memory &o_y, const occa::memory &o_z, const occa::memory &o_session);

  void setPoints(const std::vector<dfloat>& x, const std::vector<dfloat>& y, const std::vector<dfloat>& z);
  void setPoints(const std::vector<dfloat>& x, const std::vector<dfloat>& y, const std::vector<dfloat>& z, const std::vector<dlong>& session);

  // set timer level
  void setTimerLevel(TimerLevel level);
  TimerLevel getTimerLevel() const;

  // set timer name
  // this is used to prefix the timer names
  void setTimerName(std::string name);

  // for multi-session cases, query distance function to "INT" bounds
  occa::memory distanceINT();

private:
  mesh_t *mesh;
  double newton_tol;
  std::string timerName = "";
  TimerLevel timerLevel = TimerLevel::None;
  std::unique_ptr<findpts::findpts_t> findpts_;
  findpts::data_t data_;
  bool mySession;

  bool findCalled = false;

  bool pointsAdded = false;

  // correponds  to which setPoints overload is called
  bool useHostPoints = false;
  bool useDevicePoints = false;

  int nPoints;

  dfloat *_x;
  dfloat *_y;
  dfloat *_z;
  dlong *_session;

  occa::memory _o_x;
  occa::memory _o_y;
  occa::memory _o_z;
  occa::memory _o_session;
  occa::memory _o_distanceINT;

  // for storing host points to output when a particle leaves the domain
  std::vector<dfloat> h_x_vec;
  std::vector<dfloat> h_y_vec;
  std::vector<dfloat> h_z_vec;
};

#endif
