#if !defined(nekrs_interp_hpp_)
#define nekrs_interp_hpp_

#include "nrs.hpp"
#include "ogs_FINDPTS.hpp"

// Contains data for doing interpolations on a particular mesh
struct interp_t{
  nrs_t* nrs;
  double newton_tol;
  ogs_findpts_t* findpts;


  interp_t(nrs_t* nrs_, dfloat newton_tol_ = 0);
  ~interp_t();

  // Finds the process, element, and reference coordinates of the given points
  void findPoints(const dfloat*const * x,   const dlong xStride[],
                        dlong*  code,  const dlong  codeStride,
                        dlong*  proc,  const dlong  procStride,
                        dlong*  el,    const dlong    elStride,
                        dfloat* r,     const dlong     rStride,
                        dfloat* dist2, const dlong dist2Stride,
                  dlong n, bool printWarnings=true);

  // Evaluates the points using the (code, proc, el, r) tuples computed by findPoints
  void evalPoints(const dfloat* fields, const dlong nFields,
                  const dlong*   code,  const dlong  code_stride,
                  const dlong*   proc,  const dlong  proc_stride,
                  const dlong*   el,    const dlong    el_stride,
                  const dfloat*  r,     const dlong     r_stride,
                        dfloat** out,   const dlong   out_stride[],
                  dlong n);
  // Evaluates the points using the (code, proc, el, r) tuples computed by findPoints
  void evalPoints(occa::memory fields, dlong nFields,
                  const dlong*   code,  const dlong  code_stride,
                  const dlong*   proc,  const dlong  proc_stride,
                  const dlong*   el,    const dlong    el_stride,
                  const dfloat*  r,     const dlong     r_stride,
                        dfloat** out,   const dlong   out_stride[],
                  dlong n);

  // Evalutes points located on this process
  // Given a (code, proc, el, r) tuple computed by findPoints, proc must be this
  // process's rank, code must be 0 or 1, and el and r are passed to this function
  void evalLocalPoints(const dfloat* fields, const dlong nFields,
                       const dlong*   el,    const dlong    el_stride,
                       const dfloat*  r,     const dlong     r_stride,
                             dfloat** out,   const dlong   out_stride[],
                       dlong n);
  // Evalutes points located on this process
  // Given a (code, proc, el, r) tuple computed by findPoints, proc must be this
  // process's rank, code must be 0 or 1, and el and r are passed to this function
  void evalLocalPoints(occa::memory fields, const dlong nFields,
                       const dlong*   el,   const dlong    el_stride,
                       const dfloat*  r,    const dlong     r_stride,
                             dfloat** out,  const dlong   out_stride[],
                       dlong n);

  // Evaluates a field at the given points
  template<typename fieldPtr>
  void interpField(fieldPtr fields, dlong nFields,
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

    findPoints(    x, x_stride,
                code, 1*sizeof(dlong),
                proc, 1*sizeof(dlong),
                  el, 1*sizeof(dlong),
                   r, 3*sizeof(dfloat),
               dist2, 1*sizeof(dfloat),
               n);

    evalPoints(fields, nFields,
                code, 1*sizeof(dlong),
                proc, 1*sizeof(dlong),
                  el, 1*sizeof(dlong),
                   r, 3*sizeof(dfloat),
                 out, out_stride,
               n);

    delete[] iwork;
    delete[] rwork;
  }

  // Evaluates the velocity at the given points
  void interpVelocity(      dfloat *uvw_base[], const dlong uvw_stride[],
                      const dfloat *xyz_base[], const dlong xyz_stride[],
                      dlong n)
  {
    interpField(this->nrs->o_U, 3, xyz_base, xyz_stride, uvw_base, uvw_stride, n);
  }
};

#endif
