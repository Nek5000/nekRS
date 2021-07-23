#if !defined(nekrs_particle_hpp_)
#define nekrs_particle_hpp_

#include "nrs.hpp"
#include "gslib.h"
#include "ogs_FINDPTS.hpp"

extern "C" {
  struct hash_data_2 {
    ulong hash_n;
    struct dbl_range bnd[2];
    double fac[2];
    uint *offset;
  };
  struct findpts_data_2 {
    struct crystal cr;
    struct findpts_local_data_2 local;
    struct hash_data_2 hash;
  };
  struct hash_data_3 {
    ulong hash_n;
    struct dbl_range bnd[3];
    double fac[3];
    uint *offset;
  };
  struct findpts_data_3 {
    struct crystal cr;
    struct findpts_local_data_3 local;
    struct hash_data_3 hash;
  };
}

// Contains a set of particles and the information needed to interpolate on the mesh
template<unsigned int D, class Extra=char>
struct particle_set {
  // helper class to allow making a vector of arrays
  struct dfloat_array {
    dfloat val[D];

    dfloat_array()
    {
    }

    dfloat_array(dfloat val_[D])
    {
      for (int i = 0; i < D; ++i) {
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
    dfloat x[D];
    dfloat r[D];
    dlong  code;
    dlong  proc;
    dlong  el;
    Extra  extra;

    particle_t()
    {
    }

    particle_t(dfloat x_[D], dlong code_, dlong proc_, dlong el_, dfloat r_[D], Extra extra_)
     : code(code_), proc(proc_), el(el_), extra(extra_) {
       for (int i = 0; i < D; ++i) {
         x[i] = x_[i];
         r[i] = r_[i];
       }
    }

    particle_t(dfloat x_[D], dlong code_, dlong proc_, dlong el_, dfloat_array r_, Extra extra_)
     : code(code_), proc(proc_), el(el_), extra(extra_) {
       for (int i = 0; i < D; ++i) {
         x[i] = x_[i];
         r[i] = r_[i];
       }
    }
  };


  nrs_t *nrs;
  double newton_tol;
  ogs_findpts_t *findpts;


  std::vector<dfloat>       x[D];
  std::vector<dlong>        code;
  std::vector<dlong>        proc;
  std::vector<dlong>        el;
  std::vector<Extra>        extra;
  std::vector<dfloat_array> r;

private:
  void setup_findpts() {
    if (newton_tol < 5e-13) {
      newton_tol = 5e-13;
    }
    int npt_max = 128;
    int bb_tol = 0.01;

    mesh_t *mesh = nrs->meshV;

    dlong nmsh = mesh->N;
    dlong nelm = mesh->Nelements;

    // element geometry
    dfloat *elx[3] = {mesh->x, mesh->y, mesh->z};

    // element dimensions
    dlong n1[3] = {mesh->N+1, mesh->N+1, mesh->N+1};
    dlong m1[3] = {2*n1[0], 2*n1[1], 2*n1[2]};

    // used for # of cells in hash tables
    dlong hash_size = nelm*n1[0]*n1[1];
    if (D == 3) hash_size *= n1[2];

    MPI_Comm comm = platform_t::getInstance()->comm.mpiComm;

    findpts = ogsFindptsSetup(D, comm, elx, n1, nelm, m1, bb_tol,
                              hash_size, hash_size, npt_max, newton_tol,
                              (occa::device*)&platform_t::getInstance()->device);
  }

public:

  particle_set(nrs_t *nrs_, double newton_tol_)
      : nrs(nrs_), newton_tol(newton_tol_) {
    setup_findpts();
  }

  particle_set(particle_set& set)
      : nrs(set.nrs), newton_tol(newton_tol),
        x(set.x), code(set.code), proc(set.proc), el(set.el), extra(set.extra), r(set.r) {
    setup_findpts();
  }

  ~particle_set() {
    if (this->findpts) {
      ogsFindptsFree(this->findpts);
    }
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
    dfloat x_[D];
    for (int j = 0; j < D; ++j) {
      x_[j] = x[j][i];
    }
    return particle_t(x_, code[i], proc[i], el[i], r[i], extra[i]);
  }

  void push(particle_t particle) {
    for (int j = 0; j < D; ++j) {
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
      for (int j = 0; j < D; ++j) {
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
      for (int j = 0; j < D; ++j) {
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

    for (int d = 0; d < D; ++d) {
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
  void find(bool use_device=true, dfloat *dist2_in=nullptr) {
    dlong n = size();
    dfloat *dist2;
    if(dist2_in != nullptr){
      dist2 = dist2_in;
    }else{
      dist2 = new dfloat[n];
    }
    dfloat *x_base[D]; dlong x_stride[D];
    for (int i = 0; i < D; ++i){
      x_base[i] = x[i].data();
      x_stride[i] = sizeof(dfloat);
    }
    ogsFindpts(code.data(),       1*sizeof(dlong),
               proc.data(),       1*sizeof(dlong),
               el.data(),         1*sizeof(dlong),
               &(r.data()[0][0]), D*sizeof(dfloat),
               dist2,             1*sizeof(dfloat),
               x_base,            x_stride,
               n, findpts,
               use_device);

    dlong nfail = 0;
    for (int in = 0; in < n; ++in) {
      if (code[in] == 1) {
        if (dist2[in] > 10*newton_tol) {
          nfail += 1;
          //if (nfail < 5) write(6,'(a,1p4e15.7)')     ' WARNING: point on boundary or outside the mesh xy[z]d^2: ',     xp(in),yp(in),zp(in),rwk(in,1)
          if (nfail < 5){
            std::cerr << " WARNING: point on boundary or outside the mesh xy[z]d^2: "
                      << x[0][in] << "," << x[1][in] << ", " << x[2][in] << ", " << dist2[in] << std::endl;
          }
        }
      } else if (code[in] == 2) {
        nfail += 1;
        //if (nfail < 5) write(6,'(a,1p3e15.7)')        ' WARNING: point not within mesh xy[z]: !',        xp(in),yp(in),zp(in)
        if (nfail < 5){
          std::cerr << " WARNING: point not within mesh xy[z]: "
                    << x[0][in] << "," << x[1][in] << ", " << x[2][in] << std::endl;
        }
      }
    }
    hlong lcounts[2] = {n, nfail}, gcounts[2];
    MPI_Reduce(lcounts, gcounts, 2, MPI_HLONG, MPI_SUM, 0, platform_t::getInstance()->comm.mpiComm);
    if (platform_t::getInstance()->comm.mpiRank == 0) {
      if (gcounts[1] > 0) {
        std::cout << "Total number of points = " << gcounts[0] << ", failed = " << gcounts[1] << " done :: particle_set::find" << std::endl;
      }
    }
    if(dist2_in == nullptr) {
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

    typedef typename std::conditional<D==2, findpts_data_2, findpts_data_3>::type findpts_handle_D;
    sarray_transfer(particle_t, &transfer, proc, true, &((findpts_handle_D*)findpts->findpts_data)->cr);

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
  //   out            ... array of pointers to the output arrays (dfloat[n][D])
  //   nfld           ... number of fields
  template<typename fld_ptr>
  void interp(fld_ptr fld, dfloat *out[], dlong nfld)
  {
    for (int ifld = 0; ifld < nfld; ++ifld) {
      ogsFindptsEval(out[ifld],         1*sizeof(dfloat),
                     code.data(),       1*sizeof(dlong),
                     proc.data(),       1*sizeof(dlong),
                     el.data(),         1*sizeof(dlong),
                     &(r.data()[0][0]), D*sizeof(dfloat),
                     size(), fld, findpts);

      fld += nrs->fieldOffset;
    }
  }

  // Interpolates the fields at each particle with the assumption that all particles belong to local elements
  // this->migrate must have been called since the last change in position
  //   fld            ... source field(s), may be host pointer or occa::memory (dfloat[nrs->fieldOffset*nfld])
  //   out            ... array of pointers to the output arrays (dfloat[n][D])
  //   nfld           ... number of fields
  void interp_local(const dfloat* fld, dfloat *out[], dlong nfld)
  {
    dlong pn = size();

    if (pn == 0 || nfld == 0) {
       return;
    }

    for (int ifld = 0; ifld < nfld; ++ifld) {
      ogsFindptsLocalEval(out[ifld],         1*sizeof(dfloat),
                          el.data(),         1*sizeof(dlong),
                          &(r.data()[0][0]), D*sizeof(dfloat),
                          pn, fld, findpts);

      fld += nrs->fieldOffset;
    }
  }

  void interp_local(occa::memory fld, dfloat *out[], dlong nfld)
  {
    dlong pn = size();
    dlong out_stride = 1*sizeof(dfloat);
    dlong   r_stride = D*sizeof(dfloat);
    dlong  el_stride = 1*sizeof(dlong);

    dlong start = 0;
    while (code[start] == 2 && start < pn) ++start;
    pn -= start;

    if (pn == 0 || nfld == 0) {
       return;
    }

    occa::device device = *findpts->device;
    occa::memory workspace = device.malloc((nfld*out_stride+r_stride+el_stride)*pn,
                                           occa::dtype::byte);
    occa::memory d_out = workspace; workspace += nfld*out_stride*pn;
    occa::memory d_r   = workspace; workspace +=        r_stride*pn;
    occa::memory d_el  = workspace; workspace +=       el_stride*pn;
    d_r .copyFrom(&(r.data()[start][0]),  r_stride*pn);
    d_el.copyFrom(el.data()+start,        el_stride*pn);

    occa::memory d_out_i = d_out;
    for (int ifld = 0; ifld < nfld; ++ifld) {
      ogsFindptsLocalEval(d_out_i, out_stride,
                          d_el,     el_stride,
                          d_r,       r_stride,
                          pn, fld, findpts);

      d_out_i += out_stride*pn;
      fld += nrs->fieldOffset;
    }
    d_out_i = d_out;
    for (int ifld = 0; ifld < nfld; ++ifld) {
      d_out_i.copyTo(out[ifld]+start, out_stride*pn);
      d_out_i += out_stride*pn;
    }
  }
};

#endif
