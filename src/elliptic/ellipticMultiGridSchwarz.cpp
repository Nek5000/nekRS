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

#include <type_traits>
#include "elliptic.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <cctype>
#include <string>
#include <sstream>
#include <exception>
#include "platform.hpp"

struct ElementLengths
{
  dfloat* length_left_x;
  dfloat* length_left_y;
  dfloat* length_left_z;
  dfloat* length_middle_x;
  dfloat* length_middle_y;
  dfloat* length_middle_z;
  dfloat* length_right_x;
  dfloat* length_right_y;
  dfloat* length_right_z;
};
struct FDMOperators
{
  dfloat* Sx;
  dfloat* Sy;
  dfloat* Sz;
  dfloat* D;
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

void
compute_element_lengths(ElementLengths* lengths, elliptic_t* elliptic)
{
  mesh_t* mesh = elliptic->mesh;
  const dlong Nelements = elliptic->mesh->Nelements;
  dfloat* gllw = mesh->gllw;
  dfloat* gllz = mesh->gllz;
  lengths->length_left_x = (dfloat*) calloc(Nelements, sizeof(dfloat));
  lengths->length_left_y = (dfloat*) calloc(Nelements, sizeof(dfloat));
  lengths->length_left_z = (dfloat*) calloc(Nelements, sizeof(dfloat));
  lengths->length_middle_x = (dfloat*) calloc(Nelements, sizeof(dfloat));
  lengths->length_middle_y = (dfloat*) calloc(Nelements, sizeof(dfloat));
  lengths->length_middle_z = (dfloat*) calloc(Nelements, sizeof(dfloat));
  lengths->length_right_x = (dfloat*) calloc(Nelements, sizeof(dfloat));
  lengths->length_right_y = (dfloat*) calloc(Nelements, sizeof(dfloat));
  lengths->length_right_z = (dfloat*) calloc(Nelements, sizeof(dfloat));

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
      ABORT(EXIT_FAILURE);;
    }
    bool negative = false;
    negative |= lengths->length_middle_x[e] < -tol;
    negative |= lengths->length_middle_y[e] < -tol;
    negative |= lengths->length_middle_z[e] < -tol;
    if(negative) {
      std::cout << "Encountered negative length in middle for element e = " << e << "!\n";
      std::cout << "x,y,z = " << lengths->length_middle_x[e] << ", "
                << lengths->length_middle_y[e] << ", " << lengths->length_middle_z[e] << "\n";
      ABORT(EXIT_FAILURE);;
    }
  }

  dfloat* l = (dfloat*) calloc(mesh->Np * Nelements, sizeof(dfloat));
  for(unsigned i = 0; i < mesh->Np * Nelements; ++i)
    l[i] = 0.0;

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

  ogsGatherScatter(l, dfloatString, ogsAdd, mesh->ogs);
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
  free(l);
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
    assert(bc == NO_OP || DIRICHLET || NEUMANN);
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
  dfloat* S,
  const int nl,
  const int offset
  )
{
  for(int i = 0; i < nl; ++i)
    S[offset + nl * i] = 0.0;
}

void compute_1d_stiffness_matrix(
  dfloat* a,
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
  dfloat* ah = (dfloat*) calloc((n + 1) * (n + 1), sizeof(dfloat));
  dfloat* tmp = (dfloat*) calloc((n + 1) * (n + 1), sizeof(dfloat));
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

#define ah(i,j) ah[(i) + (n + 1) * (j)]
#define a(id1,id2) a[(id1) + nl * (id2)]
  int i0 = 0;
  if(lbc == 1) i0 = 1;
  int i1 = n;
  if(rbc == 1) i1 = n - 1;

  for(unsigned int i = 0; i < nl * nl; ++i)
    a[i] = 0.0;
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
  free(ah);
  free(tmp);
}

void compute_1d_mass_matrix(
  dfloat* b,
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
#define b(id1,id2) b[(id1) + nl * (id2)]
  int i0 = 0;
  if(lbc == 1) i0 = 1;
  int i1 = n;
  if(rbc == 1) i1 = n - 1;

  for(unsigned int i = 0; i < nl * nl; ++i)
    b[i] = 0.0;

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
void ssygv_ (
  int* ITYPE,
  char* JOBZ,
  char* UPLO,
  int* N,
  float* A,
  int* LDA,
  float* B,
  int* LDB,
  float* W,
  float* WORK,
  int* LWORK,
  int* INFO
  );
}
void solve_generalized_ev(
  dfloat* a,
  dfloat* b,
  dfloat* lam,
  int n
  )
{
  int info = 0;
  int worksize = n * n;
  dfloat* work_arr = (dfloat*) calloc(worksize, sizeof(dfloat));
  int itype = 1;
  char JOBZ = 'V';
  char UPLO = 'U';
  // copy of A, B in case anything goes wrong
  dfloat* a_copy = (dfloat*) calloc(n * n, sizeof(dfloat));
  dfloat* b_copy = (dfloat*) calloc(n * n, sizeof(dfloat));
  for(unsigned i = 0; i < n * n; ++i) {
    a_copy[i] = a[i];
    b_copy[i] = b[i];
  }
#ifdef DFLOAT_DOUBLE
  dsygv_(&itype,&JOBZ,&UPLO,&n,a,&n,b, &n, lam, work_arr, &worksize, &info);
#else
  ssygv_(&itype,&JOBZ,&UPLO,&n,a,&n,b, &n, lam, work_arr, &worksize, &info);
#endif
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
        err_logger << b_copy[i * n + j] << "\t";
      err_logger << "\n";
    }
    err_logger << "A:\n";
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j)
        err_logger << a_copy[i * n + j] << "\t";
      err_logger << "\n";
    }
    throw std::runtime_error(err_logger.str().c_str());
  }
  free(work_arr);
  free(a_copy);
  free(b_copy);
}

void compute_1d_matrices(
  dfloat* S,
  dfloat* lam,
  const int lbc,
  const int rbc,
  const double ll,
  const double lm,
  const double lr,
  const dlong e,
  elliptic_t* elliptic,
  std::string direction,
  const int nl
  )
{
  dfloat* b = (dfloat*) calloc(nl * nl, sizeof(dfloat));
  compute_1d_stiffness_matrix(S,lbc,rbc,ll,lm,lr,e,elliptic);
  compute_1d_mass_matrix(b,lbc,rbc,ll,lm,lr,e,elliptic);
  try {
    solve_generalized_ev(S,b,lam,nl);
  } catch(std::exception& failure) {
    std::cout << "Encountered error:\n";
    std::cout << failure.what();
    std::cout << "Direction " << direction << "\n";
    std::cout << "e = " << e << "\n";
    std::cout << "lbc = " << lbc << ", rbc = " << rbc << "\n";
    for(int iface = 0; iface < 6; ++iface)
      std::cout << "EToB[iface] = " << elliptic->EToB[6 * e + iface] << "\n";
    ABORT(EXIT_FAILURE);;
  }
  if(lbc > 0)
    row_zero(S,nl,0);
  if(lbc == 1)
    row_zero(S,nl,1);
  if(rbc > 0)
    row_zero(S,nl,nl - 1);
  if(rbc == 1)
    row_zero(S,nl,nl - 2);
  free(b);
}

void gen_operators(FDMOperators* op, ElementLengths* lengths, elliptic_t* elliptic)
{
  const int Nq_e = elliptic->mesh->Nq + 2;
  const int Np_e = Nq_e * Nq_e * Nq_e;
  const dlong Nelements = elliptic->mesh->Nelements;
  op->Sx = (dfloat*) calloc(Nq_e * Nq_e * Nelements,sizeof(dfloat));
  op->Sy = (dfloat*) calloc(Nq_e * Nq_e * Nelements,sizeof(dfloat));
  op->Sz = (dfloat*) calloc(Nq_e * Nq_e * Nelements,sizeof(dfloat));
  op->D = (dfloat*) calloc(Np_e * Nelements,sizeof(dfloat));

#define df(a,b) (op->D[((a) + (b) * Np_e)])
  const double eps = 1e-5;
  dfloat* lr = (dfloat*) calloc(Nq_e,sizeof(dfloat));
  dfloat* ls = (dfloat*) calloc(Nq_e,sizeof(dfloat));
  dfloat* lt = (dfloat*) calloc(Nq_e,sizeof(dfloat));
  dfloat* Sx = (dfloat*) calloc(Nq_e * Nq_e,sizeof(dfloat));
  dfloat* Sy = (dfloat*) calloc(Nq_e * Nq_e,sizeof(dfloat));
  dfloat* Sz = (dfloat*) calloc(Nq_e * Nq_e,sizeof(dfloat));
  for(dlong e = 0; e < Nelements; ++e) {
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
                        "r",
                        Nq_e);
    compute_1d_matrices(Sy,
                        ls,
                        lbs,
                        rbs,
                        lengths->length_left_y[e],
                        lengths->length_middle_y[e],
                        lengths->length_right_y[e],
                        e,
                        elliptic,
                        "s",
                        Nq_e);
    compute_1d_matrices(Sz,
                        lt,
                        lbt,
                        rbt,
                        lengths->length_left_z[e],
                        lengths->length_middle_z[e],
                        lengths->length_right_z[e],
                        e,
                        elliptic,
                        "t",
                        Nq_e);
    // store the transposes
    for(int i = 0; i < Nq_e; ++i)
      for(int j = 0; j < Nq_e; ++j) {
        const int elem_offset = Nq_e * Nq_e * e;
        const int ij = i + j * Nq_e;
        const int ji = j + i * Nq_e;
        op->Sx[elem_offset + ij] = Sx[ji];
        op->Sy[elem_offset + ij] = Sy[ji];
        op->Sz[elem_offset + ij] = Sz[ji];
      }
    unsigned l = 0;
    for(int k = 0; k < Nq_e; ++k)
      for(int j = 0; j < Nq_e; ++j)
        for(int i = 0; i < Nq_e; ++i) {
          const double diag = lr[i] + ls[j] + lt[k];
          if(diag > eps)
            df(l,e) = 1.0 / diag;
          else
            df(l,e) = 0.0;
          l += 1;
        }
  }
#undef df
  free(lr);
  free(ls);
  free(lt);
  free(Sx);
  free(Sy);
  free(Sz);
}

mesh_t* create_extended_mesh(elliptic_t* elliptic, hlong* maskedGlobalIds)
{
  mesh_t* meshRoot = elliptic->mesh;

  mesh_t* mesh = new mesh_t();
  mesh->N = meshRoot->N + 2;
  mesh->Np = (mesh->N + 1) * (mesh->N + 1) * (mesh->N + 1);
  mesh->Nelements = meshRoot->Nelements;
  mesh->Nverts = meshRoot->Nverts;
  mesh->Nfaces = meshRoot->Nfaces;
  mesh->NfaceVertices = meshRoot->NfaceVertices;
  mesh->Nnodes = meshRoot->Nnodes;

  mesh->EX = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  memcpy(mesh->EX, meshRoot->EX, mesh->Nverts * mesh->Nelements * sizeof(dfloat));
  memcpy(mesh->EY, meshRoot->EY, mesh->Nverts * mesh->Nelements * sizeof(dfloat));
  memcpy(mesh->EZ, meshRoot->EZ, mesh->Nverts * mesh->Nelements * sizeof(dfloat));
  
  mesh->faceVertices = (int*) calloc(mesh->NfaceVertices * mesh->Nfaces, sizeof(int));
  memcpy(mesh->faceVertices, meshRoot->faceVertices, mesh->NfaceVertices * mesh->Nfaces * sizeof(int));

  mesh->EToV = (hlong*) calloc(mesh->Nverts * mesh->Nelements, sizeof(hlong));
  memcpy(mesh->EToV, meshRoot->EToV, mesh->Nverts * mesh->Nelements * sizeof(hlong));

  meshParallelConnect(mesh);

  mesh->EToB = (int*) calloc(mesh->Nfaces * mesh->Nelements, sizeof(int));
  memcpy(mesh->EToB, meshRoot->EToB, mesh->Nfaces * mesh->Nelements * sizeof(int));

  meshLoadReferenceNodesHex3D(mesh, mesh->N, 1);
  meshHaloSetup(mesh);
  meshPhysicalNodesHex3D(mesh);
  meshHaloPhysicalNodes(mesh);
  meshConnectFaceNodes3D(mesh);
  meshGlobalIds(mesh);

  mesh->ogs = ogsSetup(mesh->Nelements * mesh->Np, mesh->globalIds, platform->comm.mpiComm, 1, platform->device.occaDevice());

  const int bigNum = 1E9;
  dlong Ntotal = mesh->Np * mesh->Nelements;

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  int* mapB = (int*) calloc(mesh->Nelements * mesh->Np,sizeof(int));
  for (dlong e = 0; e < mesh->Nelements; e++) {
    for (int n = 0; n < mesh->Np; n++) mapB[n + e * mesh->Np] = bigNum;
    for (int f = 0; f < mesh->Nfaces; f++) {
      const int bc = elliptic->EToB[f + e * mesh->Nfaces];
      if (bc > 0) {
        for (int n = 0; n < mesh->Nfp; n++) {
          const int fid = mesh->faceNodes[n + f * mesh->Nfp];
          mapB[fid + e * mesh->Np] = mymin(bc, mapB[fid + e * mesh->Np]);
        }
      }
    }
  }

  ogsGatherScatter(mapB, ogsInt, ogsMin, mesh->ogs);
  //use the bc flags to find masked ids
  dlong Nmasked = 0;
  for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++) {
    const dlong node = n % mesh->Np;
    bool isEdgeNode = false;
    for(int edge = 0 ; edge < mesh->NedgeNodes; ++edge)
      if(mesh->edgeNodes[edge] == node) isEdgeNode = true;
    if (isEdgeNode)
      Nmasked++;
    else if (mapB[n] == bigNum)
      mapB[n] = 0.;
    else if (mapB[n] == DIRICHLET)
      Nmasked++;
  }
  dlong* maskIds = (dlong*) calloc(Nmasked, sizeof(dlong));
  Nmasked = 0;
  for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++){
    const dlong node = n % mesh->Np;
    bool isEdgeNode = false;
    for(int edge = 0 ; edge < mesh->NedgeNodes; ++edge)
      if(mesh->edgeNodes[edge] == node) isEdgeNode = true;
    if (mapB[n] == 1) maskIds[Nmasked++] = n;
    else if (isEdgeNode) maskIds[Nmasked++] = n;
  }
  //make a masked version of the global id numbering
  memcpy(maskedGlobalIds, mesh->globalIds, mesh->Nlocal * sizeof(hlong));
  for (dlong n = 0; n < Nmasked; n++)
    maskedGlobalIds[maskIds[n]] = 0;

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
  //platform_t* platform = platform_t::getInstance();
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

  oogs::startFinish(work1, 1, 0, ogsPfloat, ogsAdd, (oogs_t*) ogsExt);

  extrude(work1, 0, one, work2, 0, onem, elliptic->mesh);
  extrude(work1, 2, one, work1, 0, one, elliptic->mesh);
  to_reg(wts, work1, elliptic->mesh);

  oogs::startFinish(wts, 1, 0, ogsPfloat, ogsAdd, (oogs_t*) ogs);

  for(dlong i = 0; i < weightSize; ++i)
    wts[i] = 1.0 / wts[i];
  o_wts = platform->device.malloc(weightSize * sizeof(pfloat), wts);
  free(work1);
  free(work2);
  free(wts);
}

void MGLevel::build(
  elliptic_t* pSolver)
{
  if(elliptic->elementType != HEXAHEDRA) {
    printf("ERROR: Unsupported element type!");
    ABORT(EXIT_FAILURE);
  }

  const dlong Nelements = elliptic->mesh->Nelements;
  const int Nq = elliptic->mesh->Nq;
  const int Np = elliptic->mesh->Np;

  if(Nq == 2 && elliptic->options.compareArgs("MULTIGRID COARSE SOLVE", "TRUE")){
    return;
  }

  hlong* maskedGlobalIdsExt;
  maskedGlobalIdsExt = (hlong*) calloc(Nelements*(Nq+2)*(Nq+2)*(Nq+2),sizeof(hlong));
  mesh_t* extendedMesh = create_extended_mesh(elliptic, maskedGlobalIdsExt);

  const int Nq_e = extendedMesh->Nq;
  const int Np_e = extendedMesh->Np;
  const dlong Nlocal_e = Nelements * Np_e;

  overlap = false;
  if (Nlocal_e * sizeof(pfloat) > elliptic_t::minFDMBytesOverlap && 
      platform->comm.mpiCommSize)
    overlap = true;

  /** create the element lengths, using the most refined level **/
  ElementLengths* lengths = (ElementLengths*) calloc(1,sizeof(ElementLengths));
  compute_element_lengths(lengths, pSolver);

  pfloat* casted_Sx = (pfloat*) calloc(Nq_e * Nq_e * Nelements, sizeof(pfloat));
  pfloat* casted_Sy = (pfloat*) calloc(Nq_e * Nq_e * Nelements, sizeof(pfloat));
  pfloat* casted_Sz = (pfloat*) calloc(Nq_e * Nq_e * Nelements, sizeof(pfloat));
  pfloat* casted_D = (pfloat*) calloc(Np_e * Nelements, sizeof(pfloat));

  FDMOperators* op = (FDMOperators*) calloc(1, sizeof(FDMOperators));
  gen_operators(op, lengths, elliptic);
  for(dlong i = 0; i < Nq_e * Nq_e * Nelements; ++i) {
    casted_Sx[i] = static_cast < pfloat > (op->Sx[i]);
    casted_Sy[i] = static_cast < pfloat > (op->Sy[i]);
    casted_Sz[i] = static_cast < pfloat > (op->Sz[i]);
  }
  for(dlong i = 0; i < Np_e * Nelements; ++i)
    casted_D[i] = static_cast < pfloat > (op->D[i]);
  free(op->Sx);
  free(op->Sy);
  free(op->Sz);
  free(op->D);
  free(op);
  free(lengths->length_left_x);
  free(lengths->length_left_y);
  free(lengths->length_left_z);
  free(lengths->length_middle_x);
  free(lengths->length_middle_y);
  free(lengths->length_middle_z);
  free(lengths->length_right_x);
  free(lengths->length_right_y);
  free(lengths->length_right_z);
  free(lengths);

  const dlong weightSize = Np * Nelements;
  o_Sx = platform->device.malloc  (Nq_e * Nq_e * Nelements * sizeof(pfloat));
  o_Sy = platform->device.malloc  (Nq_e * Nq_e * Nelements * sizeof(pfloat));
  o_Sz = platform->device.malloc  (Nq_e * Nq_e * Nelements * sizeof(pfloat));
  o_invL = platform->device.malloc  (Nlocal_e * sizeof(pfloat));
  o_work1 = platform->device.malloc  (Nlocal_e * sizeof(pfloat));
  if(!options.compareArgs("MULTIGRID SMOOTHER","RAS"))
    o_work2 = platform->device.malloc  (Nlocal_e * sizeof(pfloat));
  o_Sx.copyFrom(casted_Sx, Nq_e * Nq_e * Nelements * sizeof(pfloat));
  o_Sy.copyFrom(casted_Sy, Nq_e * Nq_e * Nelements * sizeof(pfloat));
  o_Sz.copyFrom(casted_Sz, Nq_e * Nq_e * Nelements * sizeof(pfloat));
  o_invL.copyFrom(casted_D, Nlocal_e * sizeof(pfloat));

  {
    const std::string suffix = std::string("_") + std::to_string(Nq_e-1) + std::string("pfloat");
    preFDMKernel = platform->kernels.get("preFDM" + suffix);
    fusedFDMKernel = platform->kernels.get("fusedFDM" + suffix);
    postFDMKernel = platform->kernels.get("postFDM" + suffix);
  }

  const oogs_mode oogsMode = OOGS_AUTO;
  ogsExt = (void*) oogs::setup(Nelements * Np_e, maskedGlobalIdsExt, 1, 0,
                               ogsPfloat, platform->comm.mpiComm, 1, platform->device.occaDevice(),
                               NULL, oogsMode);

  ogsExtOverlap = NULL;
  if(overlap) {
    occa::memory o_Su = platform->device.malloc(mesh->Nlocal * sizeof(pfloat));
    const dlong Nelements = elliptic->mesh->Nelements;
    auto callback = [&]() {
      if(options.compareArgs("MULTIGRID SMOOTHER","RAS")) {
        if(mesh->NlocalGatherElements)
          fusedFDMKernel(mesh->NlocalGatherElements,mesh->o_localGatherElementList,
                         o_Su,o_Sx,o_Sy,o_Sz,o_invL,elliptic->o_invDegree,o_work1);
      } else {
         if(mesh->NlocalGatherElements)
           fusedFDMKernel(mesh->NlocalGatherElements,mesh->o_localGatherElementList,
                          o_work2,o_Sx,o_Sy,o_Sz,o_invL,o_work1);
      }
    };
    ogsExtOverlap = (void*) oogs::setup(Nelements * Np_e, maskedGlobalIdsExt, 1, 0,
                                        ogsPfloat, platform->comm.mpiComm, 1, platform->device.occaDevice(),
                                        callback, oogsMode);
    o_Su.free();
  }

  free(maskedGlobalIdsExt); 
  meshFree(extendedMesh);

  ogs = (void*) elliptic->oogs;

  generate_weights();

  free(casted_Sx);
  free(casted_Sy);
  free(casted_Sz);
  free(casted_D);
  
}

void MGLevel::smoothSchwarz(occa::memory& o_u, occa::memory& o_Su, bool xIsZero)
{
  const char* ogsDataTypeString =  ogsPfloat;
  const dlong Nelements = elliptic->mesh->Nelements;
  preFDMKernel(Nelements, o_u, o_work1);

  oogs_t *hogsExt = (overlap) ? (oogs_t*)ogsExtOverlap : (oogs_t*)ogsExt;

  oogs::startFinish(o_work1, 1, 0, ogsDataTypeString, ogsAdd, hogsExt);

  if(options.compareArgs("MULTIGRID SMOOTHER","RAS")) {
    if(!overlap){
      fusedFDMKernel(Nelements, mesh->o_elementList,
                     o_Su,o_Sx,o_Sy,o_Sz,o_invL,elliptic->o_invDegree,o_work1);
    } else {
      if(mesh->NglobalGatherElements)
        fusedFDMKernel(mesh->NglobalGatherElements,mesh->o_globalGatherElementList,
                       o_Su,o_Sx,o_Sy,o_Sz,o_invL,elliptic->o_invDegree,o_work1);
    }

    oogs::start(o_Su, 1, 0, ogsDataTypeString, ogsAdd, (oogs_t*) ogs);

    if(overlap && mesh->NlocalGatherElements)
      fusedFDMKernel(mesh->NlocalGatherElements,mesh->o_localGatherElementList,
                     o_Su,o_Sx,o_Sy,o_Sz,o_invL,elliptic->o_invDegree,o_work1);

    oogs::finish(o_Su, 1, 0, ogsDataTypeString, ogsAdd, (oogs_t*) ogs);
  } else {
    if(!overlap){
      fusedFDMKernel(Nelements, mesh->o_elementList,
                     o_work2,o_Sx,o_Sy,o_Sz,o_invL,o_work1);
    } else {
      if(mesh->NglobalGatherElements)
        fusedFDMKernel(mesh->NglobalGatherElements,mesh->o_globalGatherElementList,
                       o_work2,o_Sx,o_Sy,o_Sz,o_invL,o_work1);
    }

    oogs::start(o_work2, 1, 0, ogsDataTypeString, ogsAdd, hogsExt);

    if(overlap) {
      if(mesh->NlocalGatherElements)
        fusedFDMKernel(mesh->NlocalGatherElements,mesh->o_localGatherElementList,
                       o_work2,o_Sx,o_Sy,o_Sz,o_invL,o_work1);
    }

    oogs::finish(o_work2, 1, 0, ogsDataTypeString, ogsAdd, hogsExt);

    postFDMKernel(Nelements,o_work1,o_work2,o_Su, o_wts);

    oogs::startFinish(o_Su, 1, 0, ogsDataTypeString, ogsAdd, (oogs_t*) ogs);
  }
  ellipticApplyMask(elliptic, o_Su, pfloatString);

  const auto Nqe = mesh->Nq + 2;
  const auto Npe = Nqe * Nqe * Nqe;
  const double flopsPerElem = 12 * Nqe * Npe + Npe;
  const double flops = static_cast<double>(mesh->Nelements) * flopsPerElem;

  const double factor = std::is_same<pfloat, float>::value ? 0.5 : 1.0;
  platform->flopCounter->add(elliptic->name + " Schwarz, N=" + std::to_string(mesh->N), factor * flops);
}
