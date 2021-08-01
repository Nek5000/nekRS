#if !defined(nekrs_interp_hpp_)
#define nekrs_interp_hpp_

#include "nrs.hpp"
#include "ogs_FINDPTS.hpp"

// Contains data for doing interpolations on a particular mesh
class interp_t{
private:
  nrs_t* nrs;
  double newton_tol;
  ogs_findpts_t* findpts;

public:

  interp_t(nrs_t* nrs_, dfloat newton_tol_ = 0);
  ~interp_t();

  // Finds the process, element, and reference coordinates of the given points
  void interp_t::findPoints(dfloat* x[], dlong xStride[],
                            dlong*  code,  dlong  codeStride,
                            dlong*  proc,  dlong  procStride,
                            dlong*  el,    dlong    elStride,
                            dfloat* r,     dlong     rStride,
                            dfloat* dist2, dlong dist2Stride,
                            dlong n,
                            bool printWarnings=true);

  // Evaluates the points using the (code, proc, el, r) tuples computed by findPoints
  void interp_t::evalPoints(dfloat* fields, dlong nFields,
                            dlong*   code,  dlong  code_stride,
                            dlong*   proc,  dlong  proc_stride,
                            dlong*   el,    dlong    el_stride,
                            dfloat*  r,     dlong     r_stride,
                            dfloat** out,   dlong   out_stride,
                            dlong n);
  // Evaluates the points using the (code, proc, el, r) tuples computed by findPoints
  void interp_t::evalPoints(occa::memory fields, dlong nFields,
                            dlong*   code,  dlong  code_stride,
                            dlong*   proc,  dlong  proc_stride,
                            dlong*   el,    dlong    el_stride,
                            dfloat*  r,     dlong     r_stride,
                            dfloat** out,   dlong   out_stride,
                            dlong n);

  // Evalutes points located on this process
  // Given a (code, proc, el, r) tuple computed by findPoints, proc must be this
  // process's rank, code must be 0 or 1, and el and r are passed to this function
  void interp_t::evalLocalPoints(dfloat* fields, dlong nFields,
                                 dlong*   el,    dlong    el_stride,
                                 dfloat*  r,     dlong     r_stride,
                                 dfloat** out,   dlong   out_stride,
                                 dlong n);
  // Evalutes points located on this process
  // Given a (code, proc, el, r) tuple computed by findPoints, proc must be this
  // process's rank, code must be 0 or 1, and el and r are passed to this function
  void interp_t::evalLocalPoints(occa::memory fields, dlong nFields,
                                 dlong*   el,    dlong    el_stride,
                                 dfloat*  r,     dlong     r_stride,
                                 dfloat** out,   dlong   out_stride,
                                 dlong n);

  // Evaluates a field at the given points
  template<typename fld_ptr>
  void evalField(fld_ptr fld, dlong nfld,
                 const dfloat* x[],   const dlong x_stride[],
                       dfloat* out[], const dlong out_stride[],
                 dlong n)
  {
    dlong *iwork = new dlong[3*n];
    dfloat *rwork = new dfloat[4*n];

    dlong*  code  = iwork;
    dlong*  proc  = iwork+n;
    dlong*  el    = iwork+2*n;
    dfloat* r     = rwork;
    dfloat* dist2 = rwork+3*n;

    findPoints(    x,     x_stride,
                code,  code_stride,
                proc,  proc_stride,
                  el,    el_stride,
                   r,     r_stride,
               dist2, dist2_stride,
               n, handle);

    evalPoints(fields, nFields
                code,  code_stride,
                proc,  proc_stride,
                  el,    el_stride,
                   r,     r_stride,
               dist2, dist2_stride,
                 out,   out_stride,
               n, handle);

    delete[] iwork;
    delete[] rwork;
  }

  // Evaluates the velocity at the given points
  void evalVelocity(      dfloat *uvw_base[], const dlong uvw_stride[],
                    const dfloat *xyz_base[], const dlong xyz_stride[],
                    dlong n)
  {
    evalField(this->nrs->o_U, 3, xyz_base, xyz_stride, uvw_base, uvw_stride, n);
  }
}

#endif
