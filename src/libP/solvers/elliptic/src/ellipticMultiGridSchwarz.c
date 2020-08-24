/*

   The MIT License (MIT)

   Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

 */

#include "elliptic.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <cctype>
#include <string>
#include <sstream>
#include <exception>

struct ElementLengths
{
  std::vector < double > length_left_x;
  std::vector < double > length_left_y;
  std::vector < double > length_left_z;
  std::vector < double > length_middle_x;
  std::vector < double > length_middle_y;
  std::vector < double > length_middle_z;
  std::vector < double > length_right_x;
  std::vector < double > length_right_y;
  std::vector < double > length_right_z;
};
struct FDMOperators
{
  std::vector < double > Sx;
  std::vector < double > Sy;
  std::vector < double > Sz;
  std::vector < double > D;
};
void harmonic_mean_element_length(
  ElementLengths* lengths,
  elliptic_t* elliptic)
{
  mesh_t* mesh = elliptic->mesh;
  const dlong Nelements = mesh->Nelements;
  const int Nq = mesh->Nq;
  dfloat const* const w = mesh->gllw;
  const int i1 = 0;
  const int i2 = Nq - 1;
  const int j1 = i1;
  const int j2 = i2;
  const int k1 = i1;
  const int k2 = i2;
  const int nx = Nq - 2;
  for(dlong e = 0; e < Nelements; ++e) {
    const dlong elem_offset = e * Nq * Nq * Nq;
    /** r **/
    double lr2 = 0.0;
    double wsum = 0.0;
    for(int k = 1; k <= nx; ++k)
      for(int j = 1; j <= nx; ++j) {
        const dlong i2jk = i2 + j * Nq + k * Nq * Nq + elem_offset;
        const dlong i1jk = i1 + j * Nq + k * Nq * Nq + elem_offset;
        const double weight = w[j - 1] * w[k - 1];
        const double dist_x = mesh->x[i2jk] - mesh->x[i1jk];
        const double dist_y = mesh->y[i2jk] - mesh->y[i1jk];
        const double dist_z = mesh->z[i2jk] - mesh->z[i1jk];
        const double denom =
          dist_x * dist_x +
          dist_y * dist_y +
          dist_z * dist_z;
        lr2 += weight / denom;
        wsum += weight;
      }
    lr2 /= wsum;
    lengths->length_middle_x[e] = 1.0 / sqrt(lr2);

    /** s **/
    double ls2 = 0.0;
    wsum = 0.0;
    for(int k = 1; k <= nx; ++k)
      for(int i = 1; i <= nx; ++i) {
        const dlong ij2k = i + j2 * Nq + k * Nq * Nq + elem_offset;
        const dlong ij1k = i + j1 * Nq + k * Nq * Nq + elem_offset;
        const double weight = w[i - 1] * w[k - 1];
        const double dist_x = mesh->x[ij2k] - mesh->x[ij1k];
        const double dist_y = mesh->y[ij2k] - mesh->y[ij1k];
        const double dist_z = mesh->z[ij2k] - mesh->z[ij1k];
        const double denom =
          dist_x * dist_x +
          dist_y * dist_y +
          dist_z * dist_z;
        ls2 += weight / denom;
        wsum += weight;
      }
    ls2 /= wsum;
    lengths->length_middle_y[e] = 1.0 / sqrt(ls2);

    /** t **/
    double lt2 = 0.0;
    wsum = 0.0;
    for(int j = 1; j <= nx; ++j)
      for(int i = 1; i <= nx; ++i) {
        const dlong ijk2 = i + j * Nq + k2 * Nq * Nq + elem_offset;
        const dlong ijk1 = i + j * Nq + k1 * Nq * Nq + elem_offset;
        const double weight = w[i - 1] * w[j - 1];
        const double dist_x = mesh->x[ijk2] - mesh->x[ijk1];
        const double dist_y = mesh->y[ijk2] - mesh->y[ijk1];
        const double dist_z = mesh->z[ijk2] - mesh->z[ijk1];
        const double denom =
          dist_x * dist_x +
          dist_y * dist_y +
          dist_z * dist_z;
        lt2 += weight / denom;
        wsum += weight;
      }
    lt2 /= wsum;
    lengths->length_middle_z[e] = 1.0 / sqrt(lt2);
  }
}
ElementLengths*
compute_element_lengths(elliptic_t* elliptic)
{
  ElementLengths* lengths = new ElementLengths();
  mesh_t* mesh = elliptic->mesh;
  const dlong Nelements = elliptic->mesh->Nelements;
  dfloat* gllw = mesh->gllw;
  dfloat* gllz = mesh->gllz;
  lengths->length_left_x.resize(Nelements);
  lengths->length_left_y.resize(Nelements);
  lengths->length_left_z.resize(Nelements);
  lengths->length_middle_x.resize(Nelements);
  lengths->length_middle_y.resize(Nelements);
  lengths->length_middle_z.resize(Nelements);
  lengths->length_right_x.resize(Nelements);
  lengths->length_right_y.resize(Nelements);
  lengths->length_right_z.resize(Nelements);

  const int N = mesh->N;
  const int Nq = mesh->Nq;

  harmonic_mean_element_length(lengths,elliptic);

  // add check for small values in middle elements
  const double tol = 1e-12;
  for(dlong e = 0; e < Nelements; ++e) {
    bool failed = false;
    failed |= std::abs(lengths->length_middle_x[e]) < tol;
    failed |= std::abs(lengths->length_middle_y[e]) < tol;
    failed |= std::abs(lengths->length_middle_z[e]) < tol;
    if(failed) {
      std::cout << "Encountered length of zero in middle for element e = " << e << "!\n";
      std::cout << "x,y,z = " << lengths->length_middle_x[e] << ", "
                << lengths->length_middle_y[e] << ", " << lengths->length_middle_z[e] << "\n";
      exit(-1);
    }
    bool negative = false;
    negative |= lengths->length_middle_x[e] < -tol;
    negative |= lengths->length_middle_y[e] < -tol;
    negative |= lengths->length_middle_z[e] < -tol;
    if(negative) {
      std::cout << "Encountered negative length in middle for element e = " << e << "!\n";
      std::cout << "x,y,z = " << lengths->length_middle_x[e] << ", "
                << lengths->length_middle_y[e] << ", " << lengths->length_middle_z[e] << "\n";
      exit(-1);
    }
  }

  std::vector < double > l(mesh->Np * Nelements);
  std::fill(l.begin(), l.end(), 1.0);

  for(dlong e = 0; e < Nelements; ++e) {
    const dlong elem_offset = Nq * Nq * Nq * e;
    for(int j = 1; j < N; ++j)
      for(int k = 1; k < N; ++k) {
        l[k * Nq + j * Nq * Nq + elem_offset] = lengths->length_middle_x[e];
        l[Nq - 1 + k * Nq + j * Nq * Nq + elem_offset] = lengths->length_middle_x[e];
        l[k + 0 * Nq + j * Nq * Nq + elem_offset] = lengths->length_middle_y[e];
        l[k + (Nq - 1) * Nq + j * Nq * Nq + elem_offset] = lengths->length_middle_y[e];
        l[k + j * Nq + elem_offset] = lengths->length_middle_z[e];
        l[k + j * Nq + (Nq - 1) * Nq * Nq + elem_offset] = lengths->length_middle_z[e];
      }
  }

  ogsGatherScatter(l.data(), "double", ogsAdd, mesh->ogs);
  for(dlong e = 0; e < Nelements; ++e) {
    const dlong elem_offset = e * Nq * Nq * Nq;
    lengths->length_left_x[e] = l[1 * Nq + 1 * Nq * Nq + elem_offset] - lengths->length_middle_x[e];
    lengths->length_right_x[e] = l[Nq - 1 + 1 * Nq + 1 * Nq * Nq + elem_offset] -
                                 lengths->length_middle_x[e];
    lengths->length_left_y[e] = l[1 + 1 * Nq * Nq + elem_offset] - lengths->length_middle_y[e];
    lengths->length_right_y[e] = l[1 + (Nq - 1) * Nq + 1 * Nq * Nq + elem_offset] -
                                 lengths->length_middle_y[e];
    lengths->length_left_z[e] = l[1 + Nq + elem_offset] - lengths->length_middle_z[e];
    lengths->length_right_z[e] = l[1 + Nq + (Nq - 1) * Nq * Nq + elem_offset] -
                                 lengths->length_middle_z[e];
  }
  for(dlong e = 0; e < Nelements; ++e) {
    double length = lengths->length_left_x[e];
    if(std::abs(length) < tol || length < -tol)
      lengths->length_left_x[e] = lengths->length_middle_x[e];
    length = lengths->length_left_y[e];
    if(std::abs(length) < tol || length < -tol)
      lengths->length_left_y[e] = lengths->length_middle_y[e];
    length = lengths->length_left_z[e];
    if(std::abs(length) < tol || length < -tol)
      lengths->length_left_z[e] = lengths->length_middle_z[e];
  }
  for(dlong e = 0; e < Nelements; ++e) {
    double length = lengths->length_right_x[e];
    if(std::abs(length) < tol || length < -tol)
      lengths->length_right_x[e] = lengths->length_middle_x[e];
    length = lengths->length_right_y[e];
    if(std::abs(length) < tol || length < -tol)
      lengths->length_right_y[e] = lengths->length_middle_y[e];
    length = lengths->length_right_z[e];
    if(std::abs(length) < tol || length < -tol)
      lengths->length_right_z[e] = lengths->length_middle_z[e];
  }
  return lengths;
}
void compute_element_boundary_conditions(int* lbr,
                                         int* rbr,
                                         int* lbs,
                                         int* rbs,
                                         int* lbt,
                                         int* rbt,
                                         const int e,
                                         elliptic_t* elliptic)
{
  int fbc[6];
  const int lookup[] = {4,2,1,3,0,5};
  for(int iface = 0; iface < 6; ++iface) {
    const int id = lookup[iface];
    int bc = elliptic->EToB[6 * e + id];
    assert(bc > -1 && bc < 3);
    fbc[iface] = bc;
  }
  *lbr = fbc[0];
  *rbr = fbc[1];
  *lbs = fbc[2];
  *rbs = fbc[3];
  *lbt = fbc[4];
  *rbt = fbc[5];
}
void row_zero(
  std::vector < double >& S,
  const int nl,
  const int offset
  )
{
  for(int i = 0; i < nl; ++i)
    S.at(offset + nl * i) = 0.0;
}
void compute_1d_stiffness_matrix(
  std::vector < double >& a,
  const int lbc,
  const int rbc,
  const double ll,
  const double lm,
  const double lr,
  const dlong e,
  elliptic_t* elliptic
  )
{
  /** build ahat matrix here **/
  // Ah = D^T B D

  const int n =  elliptic->mesh->N;
  const int nl = n + 3;
  std::vector < double > ah((n + 1) * (n + 1));
  std::vector < double > tmp((n + 1) * (n + 1));
  for(int i = 0; i < n + 1; ++i)
    for(int j = 0; j < n + 1; ++j)
      tmp[i * (n + 1) + j] = elliptic->mesh->D[i * (n + 1) + j];
  for(int i = 0; i < n + 1; ++i)
    for(int j = 0; j < n + 1; ++j)
      tmp[i * (n + 1) + j] *= elliptic->mesh->gllw[i];
  // A = D^T B D
  for(int i = 0; i < n + 1; ++i)
    for(int j = 0; j < n + 1; ++j) {
      double aij = 0.0;
      for(int k = 0; k < n + 1; ++k)
        aij += elliptic->mesh->D[k * (n + 1) + i] * tmp[k * (n + 1) + j];
      ah[i + j * (n + 1)] = aij;
    }

#define ah(i,j) ah.at((i) + (n + 1) * (j))
#define a(id1,id2) a.at((id1) + nl * (id2))
  int i0 = 0;
  if(lbc == 1) i0 = 1;
  int i1 = n;
  if(rbc == 1) i1 = n - 1;

  std::fill(a.begin(),a.end(),0.);
  double fac = 2.0 / lm;
  a(1,1) = 1.0;
  a(n + 1,n + 1) = 1.0;
  for(int j = i0; j <= i1; ++j)
    for(int i = i0; i <= i1; ++i)
      a(i + 1,j + 1) = fac * ah(i,j);
  if(lbc == 0) {
    fac = 2.0 / ll;
    a(0,0) = fac * ah(n - 1,n - 1);
    a(1,0) = fac * ah(n,n - 1);
    a(0,1) = fac * ah(n - 1,n  );
    a(1,1) = a(1,1) + fac * ah(n,n  );
  }else {
    a(0,0) = 1.0;
  }
  if(rbc == 0) {
    fac = 2.0 / lr;
    a(n + 1,n + 1) = a(n + 1,n + 1) + fac * ah(0,0);
    a(n + 2,n + 1) = fac * ah(1,0);
    a(n + 1,n + 2) = fac * ah(0,1);
    a(n + 2,n + 2) = fac * ah(1,1);
  }else {
    a(n + 2,n + 2) = 1.0;
  }
#undef a
#undef ah
}
void compute_1d_mass_matrix(
  std::vector < double >& b,
  const int lbc,
  const int rbc,
  const double ll,
  const double lm,
  const double lr,
  const dlong e,
  elliptic_t* elliptic
  )
{
  const int n = elliptic->mesh->Nq - 1;
  const int nl = n + 3;
#define bh(i) elliptic->mesh->gllw[i]
#define b(id1,id2) b.at((id1) + nl * (id2))
  int i0 = 0;
  if(lbc == 1) i0 = 1;
  int i1 = n;
  if(rbc == 1) i1 = n - 1;

  std::fill(b.begin(),b.end(),0.0);

  double fac = 0.5 * lm;
  b(1,1) = 1.0;
  b(n + 1,n + 1) = 1.0;
  for(int i = i0; i <= i1; ++i)
    b(i + 1,i + 1) = fac * bh(i);
  if(lbc == 0) {
    fac = 0.5 * ll;
    b(0,0) = fac * bh(n - 1);
    b(1,1) = b(1,1) + fac * bh(n  );
  }else {
    b(0,0) = 1.0;
  }
  if(rbc == 0) {
    fac = 0.5 * lr;
    b(n + 1,n + 1) = b(n + 1,n + 1) + fac * bh(0);
    b(n + 2,n + 2) = fac * bh(1);
  } else {
    b(n + 2,n + 2) = 1.0;
  }
#undef b
}
extern "C"
{
void dsygv_ (
  int* ITYPE,
  char* JOBZ,
  char* UPLO,
  int* N,
  double* A,
  int* LDA,
  double* B,
  int* LDB,
  double* W,
  double* WORK,
  int* LWORK,
  int* INFO
  );
}
void solve_generalized_ev(
  std::vector < double >& a,
  std::vector < double >& b,
  std::vector < double >& lam
  )
{
  int n = lam.size();
  int info = 0;
  int worksize = n * n;
  double* work_arr = (double*) calloc(worksize, sizeof(double));
  int itype = 1;
  char JOBZ = 'V';
  char UPLO = 'U';
  // copy of A, B in case anything goes wrong
  std::vector < double > a_copy;
  std::vector < double > b_copy;
  std::copy(a.begin(), a.end(), std::back_inserter(a_copy));
  std::copy(b.begin(), b.end(), std::back_inserter(b_copy));
  dsygv_(&itype,&JOBZ,&UPLO,&n,a.data(),&n,b.data(), &n, lam.data(), work_arr, &worksize, &info);
  if(info != 0) {
    std::ostringstream err_logger;
    err_logger << "Error encountered in solve_generalized_ev!\n";
    if(info < 0) {
      err_logger << "Argument " << -info << " had an illegal value!\n";
    } else {
      if(info <= n) {
        err_logger <<
        "DSYEV failed to converge, as i off-diagonal elements of an intermediate tridiagonal form did not converge to zero\n";
      } else {
        info -= n;
        err_logger << "The leading minor of order " << info << " of B is not positive definite.\n"
                   <<
        "The factorization of B could not be completed and no eigenvalues/eigenvectors were computed.\n";
      }
    }

    // dump the operators
    err_logger << "B:\n";
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j)
        err_logger << b_copy.at(i * n + j) << "\t";
      err_logger << "\n";
    }
    err_logger << "A:\n";
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j)
        err_logger << a_copy.at(i * n + j) << "\t";
      err_logger << "\n";
    }
    throw std::runtime_error(err_logger.str().c_str());
  }
  free(work_arr);
}
void compute_1d_matrices(
  std::vector < double >& S,
  std::vector < double >& lam,
  const int lbc,
  const int rbc,
  const double ll,
  const double lm,
  const double lr,
  const dlong e,
  elliptic_t* elliptic,
  std::string direction
  )
{
  const int nl = lam.size();
  std::vector < double > b(nl * nl);
  compute_1d_stiffness_matrix(S,lbc,rbc,ll,lm,lr,e,elliptic);
  compute_1d_mass_matrix(b,lbc,rbc,ll,lm,lr,e,elliptic);
  try {
    solve_generalized_ev(S,b,lam);
  } catch(std::exception& failure) {
    std::cout << "Encountered error:\n";
    std::cout << failure.what();
    std::cout << "Direction " << direction << "\n";
    std::cout << "e = " << e << "\n";
    std::cout << "lbc = " << lbc << ", rbc = " << rbc << "\n";
    for(int iface = 0; iface < 6; ++iface)
      std::cout << "EToB[iface] = " << elliptic->EToB[6 * e + iface] << "\n";
    exit(-1);
  }
  if(lbc > 0)
    row_zero(S,nl,0);
  if(lbc == 1)
    row_zero(S,nl,1);
  if(rbc > 0)
    row_zero(S,nl,nl - 1);
  if(rbc == 1)
    row_zero(S,nl,nl - 2);

}
FDMOperators* gen_operators(ElementLengths* lengths, elliptic_t* elliptic)
{
  FDMOperators* op = new FDMOperators();
  const int Nq_e = elliptic->mesh->Nq + 2;
  const int Np_e = Nq_e * Nq_e * Nq_e;
  const dlong Nelements = elliptic->mesh->Nelements;
  op->Sx.resize(Nq_e * Nq_e * Nelements);
  op->Sy.resize(Nq_e * Nq_e * Nelements);
  op->Sz.resize(Nq_e * Nq_e * Nelements);
  op->D.resize(Np_e * Nelements);

#define df(a,b) (op->D.at(((a) + (b) * Np_e)))
  const double eps = 1e-5;
  for(dlong e = 0; e < Nelements; ++e) {
    std::vector < double > lr(Nq_e);
    std::vector < double > ls(Nq_e);
    std::vector < double > lt(Nq_e);
    std::vector < double > Sx(Nq_e * Nq_e);
    std::vector < double > Sy(Nq_e * Nq_e);
    std::vector < double > Sz(Nq_e * Nq_e);
    int lbr = -1, rbr = -1, lbs = -1, rbs = -1, lbt = -1, rbt = -1;
    compute_element_boundary_conditions(&lbr,&rbr,&lbs,&rbs,&lbt,&rbt,e,elliptic);
    compute_1d_matrices(Sx,
                        lr,
                        lbr,
                        rbr,
                        lengths->length_left_x[e],
                        lengths->length_middle_x[e],
                        lengths->length_right_x[e],
                        e,
                        elliptic,
                        "r");
    compute_1d_matrices(Sy,
                        ls,
                        lbs,
                        rbs,
                        lengths->length_left_y[e],
                        lengths->length_middle_y[e],
                        lengths->length_right_y[e],
                        e,
                        elliptic,
                        "s");
    compute_1d_matrices(Sz,
                        lt,
                        lbt,
                        rbt,
                        lengths->length_left_z[e],
                        lengths->length_middle_z[e],
                        lengths->length_right_z[e],
                        e,
                        elliptic,
                        "t");
    // store the transposes
    for(int i = 0; i < Nq_e; ++i)
      for(int j = 0; j < Nq_e; ++j) {
        const int elem_offset = Nq_e * Nq_e * e;
        const int ij = i + j * Nq_e;
        const int ji = j + i * Nq_e;
        op->Sx.at(elem_offset + ij) = Sx.at(ji);
        op->Sy.at(elem_offset + ij) = Sy.at(ji);
        op->Sz.at(elem_offset + ij) = Sz.at(ji);
      }
    unsigned l = 0;
    for(int k = 0; k < Nq_e; ++k)
      for(int j = 0; j < Nq_e; ++j)
        for(int i = 0; i < Nq_e; ++i) {
          const double diag = lr.at(i) + ls.at(j) + lt.at(k);
          if(diag > eps)
            df(l,e) = 1.0 / diag;
          else
            df(l,e) = 0.0;
          l += 1;
        }
  }
#undef df

  return op;
}

mesh_t* create_extended_mesh(elliptic_t* elliptic)
{
  mesh_t* meshRoot = elliptic->mesh;

  mesh_t* mesh = new mesh_t();
  mesh->rank = meshRoot->rank;
  mesh->size = meshRoot->size;
  mesh->comm = meshRoot->comm;
  mesh->device = meshRoot->device;
  mesh->N = meshRoot->N + 2;
  mesh->Np = (mesh->N + 1) * (mesh->N + 1) * (mesh->N + 1);
  mesh->Nelements = meshRoot->Nelements;
  mesh->Nverts = meshRoot->Nverts;
  mesh->EX = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  memcpy(mesh->EX, meshRoot->EX, mesh->Nverts * mesh->Nelements * sizeof(dfloat));
  memcpy(mesh->EY, meshRoot->EY, mesh->Nverts * mesh->Nelements * sizeof(dfloat));
  memcpy(mesh->EZ, meshRoot->EZ, mesh->Nverts * mesh->Nelements * sizeof(dfloat));
  mesh->EToV = (hlong*) calloc(mesh->Nverts * mesh->Nelements, sizeof(hlong));
  memcpy(mesh->EToV, meshRoot->EToV, mesh->Nverts * mesh->Nelements * sizeof(hlong));

  meshLoadReferenceNodesHex3D(mesh, mesh->N);

  int buildOnly  = 0;
  if(elliptic->options.compareArgs("BUILD ONLY", "TRUE")) buildOnly = 1;

  meshPhysicalNodesHex3D(mesh, buildOnly);
  meshHaloSetup(mesh);
  meshConnectFaceNodes3D(mesh);
  meshParallelConnectNodes(mesh, 0, buildOnly);
  mesh->ogs = ogsSetup(mesh->Nelements * mesh->Np, mesh->globalIds, mesh->comm, 1, mesh->device);

  const int bigNum = 1E9;
  dlong Ntotal = mesh->Np * mesh->Nelements;
  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  int* mapB = (int*) calloc(mesh->Nelements * mesh->Np,sizeof(int));
  for (dlong e = 0; e < mesh->Nelements; e++) {
    for (int n = 0; n < mesh->Np; n++) mapB[n + e * mesh->Np] = bigNum;
    for (int f = 0; f < mesh->Nfaces; f++) {
      int bc = mesh->EToB[f + e * mesh->Nfaces];
      if (bc > 0) {
        for (int n = 0; n < mesh->Nfp; n++) {
          int BCFlag = elliptic->BCType[bc];
          int fid = mesh->faceNodes[n + f * mesh->Nfp];
          mapB[fid + e * mesh->Np] = mymin(BCFlag,mapB[fid + e * mesh->Np]);
        }
      }
    }
  }
  ogsGatherScatter(mapB, ogsInt, ogsMin, mesh->ogs);
  //use the bc flags to find masked ids
  dlong Nmasked = 0;
  for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++) {
    if (mapB[n] == bigNum)
      mapB[n] = 0.;
    else if (mapB[n] == 1)     //Dirichlet boundary
      Nmasked++;
  }
  dlong* maskIds = (dlong*) calloc(Nmasked, sizeof(dlong));
  Nmasked = 0;
  for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++)
    if (mapB[n] == 1) maskIds[Nmasked++] = n;
  //make a masked version of the global id numbering
  free(mesh->maskedGlobalIds);
  mesh->maskedGlobalIds = (hlong*) calloc(Ntotal,sizeof(hlong));
  memcpy(mesh->maskedGlobalIds, mesh->globalIds, Ntotal * sizeof(hlong));
  for (dlong n = 0; n < Nmasked; n++)
    mesh->maskedGlobalIds[maskIds[n]] = 0;

  free(mapB);
  free(maskIds);

  return mesh;
}
// convenience function
void to_reg(pfloat* arr1,
            pfloat* arr2,
            mesh_t* mesh)
{
  const dlong Nelements = mesh->Nelements;
  const int nx = mesh->Nq;
  const int nx_e = mesh->Nq + 2;
#define arr1(r,s,t,e) (arr1[((r) + nx * ((s) + nx * ((t) + nx * (e))))])
#define arr2(r,s,t,e) (arr2[((r + 1) + nx_e * ((s + 1) + nx_e * ((t + 1) + nx_e * (e))))])
  for(dlong ie = 0; ie < Nelements; ++ie)
    for(int k = 0; k < nx; ++k)
      for(int j = 0; j < nx; ++j)
        for(int i = 0; i < nx; ++i)
          arr1(i,j,k,ie) = arr2(i,j,k,ie);

#undef arr1
#undef arr2
}
void extrude(pfloat* arr1,
             const int l1,
             const pfloat f1,
             const pfloat* arr2,
             const int l2,
             const pfloat f2,
             mesh_t* mesh)
{
  const dlong Nelements = mesh->Nelements;
  const int nx = mesh->Nq + 2;
  const int i0 = 1;
  const int i1 = nx - 1;
#define arr1(r,s,t,e) (arr1[((r) + nx * ((s) + nx * ((t) + nx * (e))))])
#define arr2(r,s,t,e) (arr2[((r) + nx * ((s) + nx * ((t) + nx * (e))))])
  for(dlong ie = 0; ie < Nelements; ++ie) {
    for(int k = i0; k < i1; ++k)
      for(int j = i0; j < i1; ++j) {
        arr1(l1,j,k,ie) = f1 * arr1(l1,j,k,ie)
                          + f2 * arr2(l2,j,k,ie);
        arr1(nx - l1 - 1,j,k,ie) = f1 * arr1(nx - l1 - 1,j,k,ie)
                                   + f2 * arr2(nx - l2 - 1,j,k,ie);
      }
    for(int k = i0; k < i1; ++k)
      for(int i = i0; i < i1; ++i) {
        arr1(i,l1,k,ie) = f1 * arr1(i,l1,k,ie)
                          + f2 * arr2(i,l2,k,ie);
        arr1(i,nx - l1 - 1,k,ie) = f1 * arr1(i,nx - l1 - 1,k,ie)
                                   + f2 * arr2(i,nx - l2 - 1,k,ie);
      }
    for(int j = i0; j < i1; ++j)
      for(int i = i0; i < i1; ++i) {
        arr1(i,j,l1,ie) = f1 * arr1(i,j,l1,ie)
                          + f2 * arr2(i,j,l2,ie);
        arr1(i,j,nx - l1 - 1,ie) = f1 * arr1(i,j,nx - l1 - 1,ie)
                                   + f2 * arr2(i,j,nx - l2 - 1,ie);
      }
  }
#undef arr1
#undef arr2
}
void MGLevel::generate_weights()
{
  const pfloat one = 1.0;
  const pfloat zero = 0.0;
  const pfloat onem = -1.0;
  const dlong Nelements = elliptic->mesh->Nelements;
  const int Nq_e = elliptic->mesh->Nq + 2;
  const int Nq = elliptic->mesh->Nq;
  dlong weightSize = Nq * Nq * Nq * Nelements;
  dlong extendedSize = Nq_e * Nq_e * Nq_e * Nelements;
  pfloat* wts = (pfloat* ) calloc(weightSize, sizeof(pfloat));
  pfloat* work1 = (pfloat* ) calloc(extendedSize, sizeof(pfloat));
  pfloat* work2 = (pfloat* ) calloc(extendedSize, sizeof(pfloat));
  for(dlong i = 0; i < extendedSize; ++i) {
    work1[i] = 1.0;
    work2[i] = 1.0;
  }
  extrude(work2, 0, zero, work1, 0, one, elliptic->mesh);

  oogs::startFinish(work1, 1, 0, ogsPfloat, ogsAdd, (oogs_t*) extendedOgs);

  extrude(work1, 0, one, work2, 0, onem, elliptic->mesh);
  extrude(work1, 2, one, work1, 0, one, elliptic->mesh);
  to_reg(wts, work1, elliptic->mesh);

  oogs::startFinish(wts, 1, 0, ogsPfloat, ogsAdd, (oogs_t*) ogs);

  for(dlong i = 0; i < weightSize; ++i)
    wts[i] = 1.0 / wts[i];
  o_wts = mesh->device.malloc < pfloat > (weightSize);
  o_wts.copyFrom(wts, weightSize * sizeof(pfloat));
  free(work1);
  free(work2);
  free(wts);
}
void MGLevel::build(
  elliptic_t* pSolver)
{
  if(elliptic->elementType != HEXAHEDRA) {
    printf("ERROR: Unsupported elements type!");
    exit(-1);
  }

  overlap = false;
  //if(Nq > 5) overlap = true

  const dlong Nelements = elliptic->mesh->Nelements;
  const int N = elliptic->mesh->Nq;
  const int Nq = elliptic->mesh->Nq;
  const int Np = elliptic->mesh->Np;

  mesh_t* extendedMesh = create_extended_mesh(elliptic);

  const int Nq_e = extendedMesh->Nq;
  const int Np_e = extendedMesh->Np;
  const dlong Nlocal_e = Nelements * Np_e;

  oogs_mode oogsMode = OOGS_AUTO; 
  if(options.compareArgs("THREAD MODEL", "SERIAL")) oogsMode = OOGS_DEFAULT;

  extendedOgs = (void*) oogs::setup(Nelements * Np_e, extendedMesh->maskedGlobalIds, 1, 0, 
                                    ogsPfloat, extendedMesh->comm, 1, extendedMesh->device,
                                    NULL, oogsMode);
  meshFree(extendedMesh);

  //if(overlap) callback = ...
  ogs = (void*) oogs::setup(Nelements * Np, elliptic->mesh->maskedGlobalIds, 1, 0,
                            ogsPfloat, elliptic->mesh->comm, 1, elliptic->mesh->device,
                            NULL, oogsMode);

  /** create the element lengths, using the most refined level **/
  ElementLengths* lengths = compute_element_lengths(pSolver);

  std::vector < pfloat > casted_Sx(Nq_e * Nq_e * Nelements);
  std::vector < pfloat > casted_Sy(Nq_e * Nq_e * Nelements);
  std::vector < pfloat > casted_Sz(Nq_e * Nq_e * Nelements);
  std::vector < pfloat > casted_D(Np_e * Nelements);
  std::vector < pfloat > casted_wts(N * N * N * Nelements);

  FDMOperators* op = gen_operators(lengths, elliptic);
  for(dlong i = 0; i < Nq_e * Nq_e * Nelements; ++i) {
    casted_Sx[i] = static_cast < pfloat > (op->Sx[i]);
    casted_Sy[i] = static_cast < pfloat > (op->Sy[i]);
    casted_Sz[i] = static_cast < pfloat > (op->Sz[i]);
  }
  for(dlong i = 0; i < Np_e * Nelements; ++i)
    casted_D[i] = static_cast < pfloat > (op->D[i]);
  delete op;
  delete lengths;


  const dlong weightSize = Np * Nelements;
  o_Sx = mesh->device.malloc < pfloat > (Nq_e * Nq_e * Nelements);
  o_Sy = mesh->device.malloc < pfloat > (Nq_e * Nq_e * Nelements);
  o_Sz = mesh->device.malloc < pfloat > (Nq_e * Nq_e * Nelements);
  o_invL = mesh->device.malloc < pfloat > (Nlocal_e);
  o_work1 = mesh->device.malloc < pfloat > (Nlocal_e);
  if(!options.compareArgs("MULTIGRID SMOOTHER","RAS"))
    o_work2 = mesh->device.malloc < pfloat > (Nlocal_e);
  o_Sx.copyFrom(casted_Sx.data(), Nq_e * Nq_e * Nelements * sizeof(pfloat));
  o_Sy.copyFrom(casted_Sy.data(), Nq_e * Nq_e * Nelements * sizeof(pfloat));
  o_Sz.copyFrom(casted_Sz.data(), Nq_e * Nq_e * Nelements * sizeof(pfloat));
  o_invL.copyFrom(casted_D.data(), Nlocal_e * sizeof(pfloat));

  generate_weights();

  char fileName[BUFSIZ], kernelName[BUFSIZ];
  for (int r = 0; r < 2; r++) {
    if ((r == 0 && mesh->rank == 0) || (r == 1 && mesh->rank > 0)) {
      occa::properties properties;
      properties += mesh->device.properties();
      properties["defines/p_Nq_e"] = Nq_e;
      properties["defines/p_threadBlockSize"] = 256;
      properties["defines/p_Nq"] = Nq;
      properties["defines/pfloat"] = pfloatString;
      properties["defines/dfloat"] = dfloatString;
      properties["defines/dlong"] = dlongString;
      properties["defines/p_restrict"] = 0;
      properties["defines/p_overlap"] = (int) overlap;
      if(options.compareArgs("MULTIGRID SMOOTHER","RAS"))
        properties["defines/p_restrict"] = 1;

      sprintf(fileName, DELLIPTIC "/okl/ellipticSchwarzSolverHex3D.okl");
      preFDMKernel = mesh->device.buildKernel(fileName, "preFDM", properties);
      fusedFDMKernel = mesh->device.buildKernel(fileName, "fusedFDM", properties);
      postFDMKernel = mesh->device.buildKernel(fileName, "postFDM", properties);
      collocateKernel = mesh->device.buildKernel(fileName, "collocate", properties);
    }
    MPI_Barrier(mesh->comm);
  }
}
void MGLevel::smoothSchwarz(occa::memory& o_u, occa::memory& o_Su, bool xIsZero)
{
  if(xIsZero) {
    const dlong Nelements = elliptic->mesh->Nelements;
    preFDMKernel(Nelements, o_u, o_work1);

    oogs::startFinish(o_work1, 1, 0, ogsPfloat, ogsAdd, (oogs_t*) extendedOgs);

    if(options.compareArgs("MULTIGRID SMOOTHER","RAS")) {
      if(mesh->NglobalGatherElements || !overlap)
        fusedFDMKernel(Nelements,mesh->NglobalGatherElements,mesh->o_globalGatherElementList,
                       o_Su,o_Sx,o_Sy,o_Sz,o_invL,o_work1, elliptic->ogs->o_invDegree);

      oogs::start(o_Su, 1, 0, ogsPfloat, ogsAdd, (oogs_t*) ogs);

      if(overlap)
        fusedFDMKernel(Nelements,mesh->NlocalGatherElements,mesh->o_localGatherElementList,
                       o_Su,o_Sx,o_Sy,o_Sz,o_invL,o_work1, elliptic->ogs->o_invDegree);

      oogs::finish(o_Su, 1, 0, ogsPfloat, ogsAdd, (oogs_t*) ogs);
    } else {
      if(mesh->NglobalGatherElements || !overlap)
        fusedFDMKernel(Nelements,mesh->NglobalGatherElements,mesh->o_globalGatherElementList,
                       o_work2,o_Sx,o_Sy,o_Sz,o_invL,o_work1);

      oogs::start(o_work2, 1, 0, ogsPfloat, ogsAdd, (oogs_t*) extendedOgs);

      if(overlap)
        fusedFDMKernel(Nelements,mesh->NlocalGatherElements,mesh->o_localGatherElementList,
                       o_work2,o_Sx,o_Sy,o_Sz,o_invL,o_work1);

      oogs::finish(o_work2, 1, 0, ogsPfloat, ogsAdd, (oogs_t*) extendedOgs);

      postFDMKernel(Nelements,o_work1,o_work2,o_Su, o_wts);

      oogs::startFinish(o_Su, 1, 0, ogsPfloat, ogsAdd, (oogs_t*) ogs);
    }
  }
}
