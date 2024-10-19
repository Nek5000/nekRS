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
  static constexpr int CODE_INTERNAL = 0;
  static constexpr int CODE_BORDER = 1;
  static constexpr int CODE_NOT_FOUND = 2;

  enum class VerbosityLevel { None, Basic, Detailed };
  pointInterpolation_t(mesh_t *mesh_,
                       MPI_Comm comm,
                       bool mySession_ = true,
                       std::vector<int> bIntID = {},
                       double bb_tol = 0.05,
                       double newton_tol = 0,
                       dlong localHashSize = 0,
                       dlong globalHashSize = 0);

  ~pointInterpolation_t() = default;

  // Finds the process, element, and reference coordinates of the given points
  void find(VerbosityLevel verbosity = VerbosityLevel::Basic, bool matchSession = false);

  void eval(dlong nFields,
            dlong inputFieldOffset,
            const occa::memory& o_in,
            dlong outputFieldOffset,
            occa::memory& o_out,
            dlong nPoints = -1,
            dlong findPtsOffset = 0);

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
  mesh_t *mesh = nullptr;
  double newton_tol = 0;
  std::string timerName;
  TimerLevel timerLevel = TimerLevel::None;

  std::unique_ptr<findpts::findpts_t> findpts_;
  findpts::data_t data_;

  bool mySession = true;
  bool findCalled = false;
  bool pointsAdded = false;

  // correponds  to which setPoints overload is called
  bool useHostPoints = false;
  bool useDevicePoints = false;

  int nPoints;

  dfloat *_x = nullptr;
  dfloat *_y = nullptr;
  dfloat *_z = nullptr;
  dlong *_session = nullptr;

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
