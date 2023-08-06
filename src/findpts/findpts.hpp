#ifndef FINDPTS_HPP
#define FINDPTS_HPP

#include "occa.hpp"
#include <mpi.h>
#include <limits>
#include <tuple>
#include <vector>
#include "nrssys.hpp"

struct crystal;
struct hash_data_3;
using hashData_t = hash_data_3;
struct evalSrcPt_t;

namespace findpts {

static constexpr int CODE_INTERNAL = 0;
static constexpr int CODE_BORDER = 1;
static constexpr int CODE_NOT_FOUND = 2;
static constexpr int dim = 3;

struct data_t {
  std::vector<dlong> code;
  std::vector<dlong> proc;
  std::vector<dlong> el;
  std::vector<dfloat> r;
  std::vector<dfloat> dist2;

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
      dist2_base[i] = std::numeric_limits<dfloat>::max();
      code_base[i] = CODE_NOT_FOUND;
    }
  }
};

enum class TimerLevel { None, Basic, Detailed };

namespace impl {
struct gslibFindptsData_t;
}

class findpts_t {
public:
  findpts_t(MPI_Comm comm,
            const dfloat *const x,
            const dfloat *const y,
            const dfloat *const z,
            const dlong Nq,
            const dlong Nelements,
            const dlong m,
            const dfloat bbox_tol,
            const hlong local_hash_size,
            const hlong global_hash_size,
            const dlong npt_max,
            const dfloat newt_tol);

  ~findpts_t();

  void find(data_t *findPtsData, occa::memory o_x, occa::memory o_y, occa::memory o_z, const dlong npt);

  void find(data_t *findPtsData,
            const dfloat *const x,
            const dfloat *const y,
            const dfloat *const z,
            const dlong npt);

  // Device versions
  void eval(const dlong npt, occa::memory o_in, data_t *findPtsData, occa::memory o_out);

  void eval(const dlong npt,
            const dlong nFields,
            const dlong inputOffset,
            const dlong outputOffset,
            occa::memory o_in,
            data_t *findPtsData,
            occa::memory o_out);

  // Host versions (copies to device when needed)
  void eval(const dlong npt, dfloat *in, data_t *findPtsData, dfloat *out);

  void eval(const dlong npt,
            const dlong nFields,
            const dlong inputOffset,
            const dlong outputOffset,
            dfloat *in,
            data_t *findPtsData,
            dfloat *out);

  // set timer level
  void setTimerLevel(TimerLevel level) { timerLevel = level; }
  TimerLevel getTimerLevel() const { return timerLevel; }

  // set timer name
  // this is used to prefix the timer names
  void setTimerName(std::string name) { timerName = name; }

  crystal *crystalRouter();

  // For use in, e.g., nek-nek
  // If altering code, proc, el, r, or dist2 after a find call,
  // update device arrays with this function
  void update(data_t &data);

private:
  static constexpr int maxFields = 30;

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

  // data for elx
  occa::memory o_x;
  occa::memory o_y;
  occa::memory o_z;

  // results after find call
  occa::memory o_code;
  occa::memory o_proc;
  occa::memory o_el;
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
                    dfloat *const r,
                    dfloat *const dist2,
                    const dfloat *const x,
                    const dfloat *const y,
                    const dfloat *const z,
                    const int pn);

  void findptsLocal(int *const code,
                    int *const el,
                    dfloat *const r,
                    dfloat *const dist2,
                    occa::memory o_xint,
                    occa::memory o_yint,
                    occa::memory o_zint,
                    const int pn);

  template <typename OutputType>
  void findptsEvalImpl(dfloat *out,
                       const int *const code_base,
                       const int *const proc_base,
                       const int *const el_base,
                       const dfloat *const r_base,
                       const int npt,
                       const int nFields,
                       const int inputOffset,
                       const int outputOffset,
                       const dfloat *const in,
                       hashData_t &hash,
                       crystal &cr);

  template <typename OutputType>
  void findptsEvalImpl(occa::memory &o_out,
                       const int *const code_base,
                       const int *const proc_base,
                       const int *const el_base,
                       const dfloat *const r_base,
                       const int npt,
                       const int nFields,
                       const int inputOffset,
                       const int outputOffset,
                       occa::memory &o_in,
                       hashData_t &hash,
                       crystal &cr);

  template <typename OutputType>
  void findptsLocalEvalInternal(OutputType *opt,
                                const evalSrcPt_t *spt,
                                const int pn,
                                const int nFields,
                                const int inputOffset,
                                const int outputOffset,
                                occa::memory &o_in);
};

} // namespace findpts

#endif
