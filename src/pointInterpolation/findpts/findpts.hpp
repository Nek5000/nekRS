#ifndef FINDPTS_HPP
#define FINDPTS_HPP

#include "occa.hpp"
#include <mpi.h>
#include <limits>
#include <tuple>
#include <vector>
#include "nekrsSys.hpp"

struct crystal;
struct hash_data_3;
using hashData_t = hash_data_3;
struct evalSrcPt_t;

namespace findpts
{

static constexpr int CODE_INTERNAL = 0;
static constexpr int CODE_BORDER = 1;
static constexpr int CODE_NOT_FOUND = 2;
static constexpr int dim = 3;

// src cache on target
struct cache_t { 
  occa::memory o_el;
  occa::memory o_r; 
  std::vector<dlong> proc;
  std::vector<dlong> index;
};

struct data_t {
  bool updateCache = true;
  cache_t cache;

  std::vector<dlong> code;
  std::vector<dlong> proc;
  std::vector<dlong> el;
  std::vector<dfloat> r;
  std::vector<dfloat> dist2; // distance squared from found to requested point (in xyz space)

  dlong *code_base;
  dlong *proc_base;
  dlong *el_base;
  dfloat *r_base;
  dfloat *dist2_base;

  data_t() {}

  data_t(int npt)
  {
    code = std::vector<dlong>(npt, 0);
    proc = std::vector<dlong>(npt, 0);
    el = std::vector<dlong>(npt, 0);
    r = std::vector<dfloat>(3 * npt, 0);
    dist2 = std::vector<dfloat>(npt, 0);

    code_base = code.data();
    proc_base = proc.data();
    el_base = el.data();
    r_base = r.data();
    dist2_base = dist2.data();

    for (dlong i = 0; i < npt; ++i) {
      dist2_base[i] = 1e30;
      code_base[i] = CODE_NOT_FOUND;
    }
  }
};

enum class TimerLevel { None, Basic, Detailed };

namespace impl
{
struct gslibFindptsData_t;
}

class findpts_t
{
public:
  findpts_t(MPI_Comm comm,
            const dfloat *const x,
            const dfloat *const y,
            const dfloat *const z,
            const dlong Nq,
            const dlong Nelements,
            const dlong m,
            const dfloat bbox_tol,
            const dlong local_hash_size,
            const dlong global_hash_size,
            const dlong npt_max,
            const dfloat newt_tol);

  // ctor with multi-session support
  findpts_t(MPI_Comm comm,
            const dfloat *const x,
            const dfloat *const y,
            const dfloat *const z,
            const dlong Nq,
            const dlong Nelements,
            const dlong m,
            const dfloat bbox_tol,
            const dlong local_hash_size,
            const dlong global_hash_size,
            const dlong npt_max,
            const dfloat newt_tol,
            const dlong sessionId_,
            const dfloat *const distfint);

  ~findpts_t();

  void find(data_t *findPtsData,
            const occa::memory &o_x,
            const occa::memory &o_y,
            const occa::memory &o_z,
            const dlong npt);
  void find(data_t *findPtsData,
            const occa::memory &o_x,
            const occa::memory &o_y,
            const occa::memory &o_z,
            const occa::memory &o_sess,
            const dlong sessionIdMatch,
            const dlong npt);

  void find(data_t *findPtsData,
            const dfloat *const x,
            const dfloat *const y,
            const dfloat *const z,
            const dlong npt);
  void find(data_t *findPtsData,
            const dfloat *const x,
            const dfloat *const y,
            const dfloat *const z,
            const dlong *const session,
            const dlong sessionIdMatch,
            const dlong npt);

  // Device versions
  void eval(const dlong npt, const occa::memory &o_in, data_t *findPtsData, occa::memory &o_out);

  void eval(const dlong npt,
            const dlong offset,
            const dlong nFields,
            const dlong inputOffset,
            const dlong outputOffset,
            const occa::memory &o_in,
            data_t *findPtsData,
            occa::memory &o_out);

  // set timer level
  void setTimerLevel(TimerLevel level)
  {
    timerLevel = level;
  }

  TimerLevel getTimerLevel() const
  {
    return timerLevel;
  }

  // set timer name
  // this is used to prefix the timer names
  void setTimerName(std::string name)
  {
    timerName = name;
  }

  crystal *crystalRouter();

private:
  static constexpr int maxFields = 30;

  bool useMultiSessionSupport = false;
  int sessionId = 0;

  MPI_Comm comm;
  int rank;
  TimerLevel timerLevel = TimerLevel::None;
  std::string timerName = "";
  dfloat tol;
  crystal *cr;
  hashData_t *hash;
  occa::kernel localEvalKernel;
  occa::kernel localEvalMaskKernel;
  occa::kernel localKernel;

  occa::stream defaultStream;
  occa::stream localEvalStream;

  // data for elx
  occa::memory o_x;
  occa::memory o_y;
  occa::memory o_z;

  // distance function -- used for multi-session support
  occa::memory o_distfint;

  // results after find call

  // 0 - inside an element
  // 1 - closest point on a border
  //     (perhaps exactly, or maybe just near --- check dist2)
  // 2 - not found (bbox_tol controls cut-off between code 1 and 2)
  occa::memory o_code;

  // remote processor on which the point was found
  occa::memory o_proc;
  // element on remote processor in which the point was found
  occa::memory o_el;
  // parametric coordinates for point
  occa::memory o_r;

  // data for wtend
  occa::memory o_wtend_x;
  occa::memory o_wtend_y;
  occa::memory o_wtend_z;

  // SoA variant of obbox
  occa::memory o_c;
  occa::memory o_A;
  occa::memory o_min;
  occa::memory o_max;

  // hash data
  occa::memory o_offset;
  dlong hash_n;

  occa::memory o_hashMin;
  occa::memory o_hashFac;

  // void* is required here to avoid having gslib header dependency
  // since this will ultimately hold a C-style object,
  // we cannot simply forward declare it
  void *_findptsData;

  void findptsLocal(int *const code,
                    int *const el,
                    int *const elsid,
                    dfloat *const r,
                    dfloat *const dist2,
                    dfloat *const disti,
                    const dfloat *const x,
                    const dfloat *const y,
                    const dfloat *const z,
                    const dlong *const sess,
                    const int sessionIdMatch,
                    const int pn);

  void findptsLocal(int *const code,
                    int *const el,
                    int *const elsid,
                    dfloat *const r,
                    dfloat *const dist2,
                    dfloat *const disti,
                    occa::memory o_xint,
                    occa::memory o_yint,
                    occa::memory o_zint,
                    occa::memory o_sess,
                    const int sessionIdMatch,
                    const int pn);

  template <typename OutputType>
  void findptsEvalImpl(occa::memory &o_out,
                       dlong offset,
                       data_t *findPtsData,
                       const dlong npt,
                       const int nFields,
                       const dlong inputOffset,
                       const dlong outputOffset,
                       const occa::memory &o_in,
                       hashData_t &hash,
                       crystal &cr);

};

} // namespace findpts

#endif
