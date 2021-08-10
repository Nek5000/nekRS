#if !defined(nekrs_particle_hpp_)
#define nekrs_particle_hpp_

#include "nrs.hpp"
#include "gslib.h"
#include "ogs_FINDPTS.hpp"

// Contains a set of particles and the information needed to interpolate on the mesh
template<class Extra=char>
struct particle_set {
  // helper class to allow making a vector of arrays
  struct dfloat_array {
    dfloat val[3];

    dfloat_array()
    {
    }

    dfloat_array(dfloat val_[3])
    {
      for (int i = 0; i < 3; ++i) {
        val[i] = val_[i];
      }
    }

    dfloat& operator[](int i)
    {
      return val[i];
    }
    dfloat operator[](int i) const
    {
      return val[i];
    }
  };
  // helper class to pass components of a single particle
  struct particle_t {
    dfloat x[3];
    dfloat r[3];
    dlong  code;
    dlong  proc;
    dlong  el;
    Extra  extra;

    particle_t()
    {
    }

    particle_t(dfloat x_[3], dlong code_, dlong proc_, dlong el_, dfloat r_[3], Extra extra_)
     : code(code_), proc(proc_), el(el_), extra(extra_) {
       for (int i = 0; i < 3; ++i) {
         x[i] = x_[i];
         r[i] = r_[i];
       }
    }

    particle_t(dfloat x_[3], dlong code_, dlong proc_, dlong el_, dfloat_array r_, Extra extra_)
     : code(code_), proc(proc_), el(el_), extra(extra_) {
       for (int i = 0; i < 3; ++i) {
         x[i] = x_[i];
         r[i] = r_[i];
       }
    }
  };


  std::shared_ptr<interp::handle_t> handle;


  std::vector<dfloat>       x[3];
  std::vector<dlong>        code;
  std::vector<dlong>        proc;
  std::vector<dlong>        el;
  std::vector<Extra>        extra;
  std::vector<dfloat_array> r;

  particle_set(nrs_t *nrs_, double newton_tol_)
      : handle(interp::setup(nrs_, newton_tol_))
  {
  }

  particle_set(particle_set& set)
      : handle(set.handle),
        x(set.x), code(set.code), proc(set.proc), el(set.el), extra(set.extra), r(set.r)
  {
  }

  ~particle_set() {
  }

  //// Set management ////
  void reserve(int n) {
    for (int i = 0; i < D; ++i) {
      x[i].reserve(n);
    }
    code.reserve(n);
    proc.reserve(n);
    el.reserve(n);
    r.reserve(n);
    extra.reserve(n);
  }

  dlong size() const {
    return x[0].size();
  }

  dlong capacity() const {
    return x[0].capacity();
  }

  particle_t operator[](int i) {
    dfloat x_[3];
    for (int j = 0; j < 3; ++j) {
      x_[j] = x[j][i];
    }
    return particle_t(x_, code[i], proc[i], el[i], r[i], extra[i]);
  }

  void push(particle_t particle) {
    for (int j = 0; j < 3; ++j) {
      x[j].push_back(particle.x[j]);
    }
    code.push_back(particle.code);
    proc.push_back(particle.proc);
    el.push_back(particle.el);
    r.push_back(particle.r);
    extra.push_back(particle.extra);
  }

  particle_t remove(int i) {
    particle_t part;
    if (i == size() -1) {
      // just pop the last element
      for (int j = 0; j < 3; ++j) {
        part.x[j] = x[j].back(); x[j].pop_back();
        part.r[j] = r.back()[j];
      }
      r.pop_back();
      part.code   = code.back();  code.pop_back();
      part.proc   = proc.back();  proc.pop_back();
      part.el     = el.back();    el.pop_back();
      part.extra  = extra.back(); extra.pop_back();
    } else {
      // swap last element to i'th position
      for (int j = 0; j < 3; ++j) {
        part.x[j] = x[j][i];  x[j][i]  = x[j].back();  x[j].pop_back();
        part.r[j] = r[i][j];  r[i][j]  = r.back()[j];
      }
      r.pop_back();
      part.code   = code[i];  code[i]  = code.back();  code.pop_back();
      part.proc   = proc[i];  proc[i]  = proc.back();  proc.pop_back();
      part.el     = el[i];    el[i]    = el.back();    el.pop_back();
      part.extra  = extra[i]; extra[i] = extra.back(); extra.pop_back();
    }
    return part;
  }

  void swap(int i, int j) {
    if (i == j) return;

    for (int d = 0; d < 3; ++d) {
      std::swap(x[d][i], x[d][j]);
      std::swap(r[d][i], r[d][j]);
    }
    std::swap(code[i],  code[j]);
    std::swap(proc[i],  proc[j]);
    std::swap(el  [i],  el  [j]);
    std::swap(extra[i], extra[j]);
  }

  //// particle operations ////

  // Locates the element and process for each particle
  void find(bool printWarnings=true, dfloat *dist2In=nullptr, dlong dist2Stride = 1) {
    dlong n = size();
    dfloat *dist2;
    if(dist2In != nullptr){
      dist2 = dist2In;
    }else{
      dist2 = new dfloat[n];
      dist2Stride = 1;
    }
    dfloat *xBase[3]; dlong xStride[3];
    for (int i = 0; i < 3; ++i){
      xBase[i] = x[i].data();
      xStride[i] = 1;
    }
    interp.findPoints(xBase, xStride,
                      code.data(),       1,
                      proc.data(),       1,
                        el.data(),       1,
                      &(r.data()[0][0]), 3,
                      dist2,             dist2Stride,
                      size(), printWarnings);
    if(dist2In == nullptr) {
      delete[] dist2;
    }
  }

  // Moves each particle to the process that owns it's current element
  // this->find must have been called since the last change in position
  void migrate() {
    int mpi_rank = platform_t::getInstance()->comm.mpiRank;

    struct array transfer;
    array_init(particle_t, &transfer, 128);

    int index = 0;
    int unfound_count = 0;
    while (index < size()) {
      if (code[index] == 2) {
        swap(index, unfound_count);
        ++unfound_count;
        ++index;
      } else if (proc[index] != mpi_rank) {
        // remove index'th element and move the last point to index'th storage
        array_reserve(particle_t, &transfer, transfer.n+1);
        ((particle_t*)transfer.ptr)[transfer.n] = remove(index);
        ++transfer.n;
      } else {
        // keep point on this process
        ++index;
      }
    }

    sarray_transfer(particle_t, &transfer, proc, true, ogsCrystalRouter(interp->findpts));

    reserve(size() + transfer.n);
    particle_t *transfer_ptr = (particle_t*)transfer.ptr;
    for (int i = 0; i < transfer.n; ++i) {
      transfer_ptr[i].proc = mpi_rank; // sarray_transfer sets proc to be the sender
      push(transfer_ptr[i]);
    }

    array_free(&transfer);
  }

  // Interpolates the fields at each particle
  // this->find must have been called since the last change in position
  //   fld            ... source field(s), may be host pointer or occa::memory (dfloat[nrs->fieldOffset*nfld])
  //   out            ... array of pointers to the output arrays (dfloat[n][3])
  //   nfld           ... number of fields
  template<typename fieldPtr>
  void interp(fieldPtr field, dfloat *out[], dlong nField)
  {
    dlong offset = 0;
    while (offset < pn && code[offset] == 2) ++offset;
    pn -= offset;

    dfloat **outOffset = new dfloat*[nFields];
    for (dlong i = 0; i < nFields; ++i) {
      outOffset[i] = out[i] + offset;
    }

    interp.evalPoints(field, nField,
                      code.data()+offset,     1*sizeof(dlong),
                      proc.data()+offset,     1*sizeof(dlong),
                      el.data()+offset,       1*sizeof(dlong),
                      &(r.data()[offset][0]), D*sizeof(dfloat),
                      outOffset,              1*sizeof(dfloat),
                      size());
    delete[] outOffset;
  }

  // Interpolates the fields at each particle with the assumption that all particles belong to local elements
  // this->migrate must have been called since the last change in position
  //   fld            ... source field(s), may be host pointer or occa::memory (dfloat[nrs->fieldOffset*nfld])
  //   out            ... array of pointers to the output arrays (dfloat[n][3])
  //   nfld           ... number of fields
  template<class fieldPtr>
  void interpLocal(fieldPtr field, dfloat *out[], dlong nFields)
  {
    dlong pn = size();
    dlong out_stride = 1*sizeof(dfloat);
    dlong   r_stride = D*sizeof(dfloat);
    dlong  el_stride = 1*sizeof(dlong);

    dlong offset = 0;
    while (offset < pn && code[offset] == 2) ++offset;
    pn -= offset;


    dfloat **outOffset = new dfloat*[nFields];
    for (dlong i = 0; i < nFields; ++i) {
      outOffset[i] = out[i] + offset;
    }

    interp.evalLocalPoints(field, nField
                           el.data()+offset, 1,
                            r.data()+offset, 3,
                           outOffset,        1,
                           pn);
    delete[] outOffset
  }
};

#endif
