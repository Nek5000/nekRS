/*
 * -----------------------------------------------------------------
 * Clone of CVODE's SPGMR solver with some modifications;
 * - mixed precision (basis in reduced precision pfloat)
 * - no restarts
 * - no preconditioning
 * - always scale and cache s2Inv 
 * -----------------------------------------------------------------
 */

// - make Kthre user parameter

#ifdef ENABLE_CVODE

#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_math.h>

#include "sundials_context_impl.h"
#include "sundials_logger_impl.h"

#include "nrssys.hpp"
#include "platform.hpp"
#include "linAlg.hpp"

#include "cbGMRES.hpp"
#include "getN_VectorMemory.hpp"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

#define SPGMR_CONTENT(S)  ( (SUNLinearSolverContent_SPGMR)(S->content) )
#define LASTFLAG(S)       ( SPGMR_CONTENT(S)->last_flag )

static const std::string timerName = "cvode_t::";

static occa::kernel axmyzKernel;
static occa::kernel axpbyKernel;
static occa::kernel linearCombinationKernel;
static occa::kernel innerProdMultiKernel;

static dlong N;


static N_Vector s2Inv;
static N_Vector vtemp2;

static occa::memory o_V;

static occa::memory o_b;
static occa::memory o_v;
static occa::memory o_x;
static occa::memory o_xcor;
static occa::memory o_vtemp;
static occa::memory o_vtemp2;
static occa::memory o_s1;
static occa::memory o_s2;
static occa::memory o_s2Inv;

static occa::memory o_yg;
static occa::memory o_stemp;

static occa::memory o_scratch;
realtype *scratch;

// o_z[n] = alpha*o_x[n]*o_y[n]
static void axmyz(const dfloat alpha,
                  const occa::memory &o_x,
                  const occa::memory &o_y,
                  occa::memory &o_z,
                  const dlong xOffset = 0)
{
  axmyzKernel(N, xOffset*N, alpha, o_x, o_y, o_z);
}

// o_y[n] = beta*o_y[n] + alpha*o_x[n]
static void axpby(const dfloat alpha,
                  const occa::memory &o_x,
                  const dfloat beta,
                  occa::memory &o_y,
                  const dlong xOffset = 0,
                  const dlong yOffset = 0)
{
  axpbyKernel(N, xOffset*N, yOffset*N, alpha, o_x, beta, o_y);
}

// o_z[n] = y_{Nfields} * c_{Nfields} + \sum_{i=0}^{Nfields-1} c_i * x_i 
static void linearCombination(const dlong Nfields,
                              const occa::memory &o_coeff,
                              const occa::memory &o_x,
                              const occa::memory &o_y,
                              occa::memory &o_z)
{
  linearCombinationKernel(N, Nfields, N, o_coeff, o_x, o_y, o_z);
}

static void innerProdMulti(const dlong NVec,
                           const occa::memory &o_x,
                           const occa::memory &o_y,
                           MPI_Comm _comm,
                           dfloat *result,
                           const dlong yOffset = 0)
{
  const int Nblock = (N + BLOCKSIZE - 1) / BLOCKSIZE;
  const size_t Nbytes = NVec * Nblock * sizeof(dfloat);

  {
    innerProdMultiKernel(Nblock,
                         N,
                         1, // Nfields
                         N, // fieldOffset
                         NVec,
                         yOffset,
                         o_x,
                         o_y,
                         o_scratch);

    o_scratch.copyTo(scratch, Nbytes);

    for (int field = 0; field < NVec; ++field) {
      dfloat dot = 0;
      for (dlong n = 0; n < Nblock; ++n) {
        dot += scratch[n + field * Nblock];
      }
      result[field] = dot;
    }
  }

  if (_comm != MPI_COMM_SELF) {
    MPI_Allreduce(MPI_IN_PLACE, result, NVec, MPI_DFLOAT, MPI_SUM, _comm);
  }
}

static void setup(N_Vector x, N_Vector b, N_Vector xcor, N_Vector vtemp,
                  N_Vector s1, N_Vector s2, N_Vector *V, int l_max)
{
  static_assert(sizeof(realtype) == sizeof(dfloat), "realtype has to match dfloat!");

  auto &kernels = platform->kernels;
  const std::string prefix = "cvode_t::";
  axpbyKernel = kernels.get(prefix + "axpby");
  axmyzKernel = kernels.get(prefix + "axmyz");
  linearCombinationKernel = kernels.get(prefix + "linearCombination");
  innerProdMultiKernel = kernels.get(prefix + "innerProdMulti");

  N = N_VGetLocalLength(x); 

  const int Nblock = (N + BLOCKSIZE - 1) / BLOCKSIZE;
  const size_t Nbytes = (l_max+1) * Nblock * sizeof(dfloat);
  o_scratch = platform->device.malloc(Nbytes); 
  scratch = (sunrealtype*) std::malloc(o_scratch.size());

  s2Inv = N_VClone(s2);
  vtemp2 = N_VClone(vtemp);

  o_b = getN_VectorMemory(sunrealtype, b);
  o_x = getN_VectorMemory(sunrealtype, x);
  o_xcor = getN_VectorMemory(sunrealtype, xcor);
  o_vtemp = getN_VectorMemory(sunrealtype, vtemp);
  o_vtemp2 = getN_VectorMemory(sunrealtype, vtemp2);
  o_s1 = getN_VectorMemory(sunrealtype, s1);
  o_s2 = getN_VectorMemory(sunrealtype, s2);
  o_s2Inv = getN_VectorMemory(sunrealtype, s2Inv);

  o_yg = platform->device.malloc<sunrealtype>(l_max+1);
  o_stemp = platform->device.malloc<sunrealtype>(l_max+1);

  N_VDestroyVectorArray(V, l_max+1);
  o_V = platform->device.malloc<pfloat>((l_max+1) * N_VGetLocalLength(x));
}


/* two pass (iterative) classical Gram-Schmidt */
static void CGS2(realtype **h, int k, int p, realtype *stemp, realtype *new_vk_norm, occa::memory& o_omega)  
{
  const int k_minus_1 = k - 1;

  const int nIterMax = 2;
  const realtype Kthres = 100;

  const realtype vk_norm0  = SUNRsqrt(platform->linAlg->innerProd(N, o_omega, o_omega, platform->comm.mpiComm)); 

  for (int iter = 0; iter < nIterMax; iter++) {
    if(((*new_vk_norm) * Kthres < vk_norm0) || iter == 0) {
      innerProdMulti(k, o_V, o_omega, platform->comm.mpiComm, stemp); 
  
      for (int i=0; i < k; i++) {
        h[i][k_minus_1] = (iter > 0) ? h[i][k_minus_1] + stemp[i] : stemp[i];
        stemp[i] = -stemp[i];
      }
  
      // o_omega = stemp_[k] * o_omega + sum (o_stemp_[i] * o_v[i]) for i = 1 ... k-1 
      stemp[k] = ONE;
      o_stemp.copyFrom(stemp, k+1);
      linearCombination(k+1, o_stemp, o_V, o_omega, o_omega);

      *new_vk_norm = SUNRsqrt(platform->linAlg->innerProd(N, o_omega, o_omega, platform->comm.mpiComm));
#if 0
      if(platform->comm.mpiRank == 0) {
        std::cout << "vk_norm: " << *new_vk_norm <<  ", " << vk_norm0 << std::endl;
      }
#endif
    }
  }
}

int cbGMRES(SUNLinearSolver S, N_Vector x, N_Vector b, realtype delta)
{
  /* local data and shortcut variables */
  N_Vector *V, xcor, vtemp, s1, s2;
  realtype **Hes, *givens, *yg, *res_norm;
  realtype beta, rotation_product, r_norm, s_product, rho;
  booleantype converged;
  booleantype *zeroguess;
  int i, j, k, l, l_plus_1, l_max, krydim, ier, ntries, max_restarts, gstype;
  int *nli;
  void *A_data, *P_data;
  SUNATimesFn atimes;
  SUNPSolveFn psolve;

  /* local shortcuts for fused vector operations */
  realtype* cv;
  N_Vector* Xv;

  /* Initialize some variables */
  l_plus_1 = 0;
  krydim = 0;


  /* Make local shorcuts to solver variables. */
  if (S == NULL) return(SUNLS_MEM_NULL);
  l_max        = SPGMR_CONTENT(S)->maxl;
  max_restarts = SPGMR_CONTENT(S)->max_restarts;
  gstype       = SPGMR_CONTENT(S)->gstype;
  V            = SPGMR_CONTENT(S)->V;
  Hes          = SPGMR_CONTENT(S)->Hes;
  givens       = SPGMR_CONTENT(S)->givens;
  xcor         = SPGMR_CONTENT(S)->xcor;
  yg           = SPGMR_CONTENT(S)->yg;
  vtemp        = SPGMR_CONTENT(S)->vtemp;
  s1           = SPGMR_CONTENT(S)->s1;
  s2           = SPGMR_CONTENT(S)->s2;


  A_data       = SPGMR_CONTENT(S)->ATData;
  P_data       = SPGMR_CONTENT(S)->PData;
  atimes       = SPGMR_CONTENT(S)->ATimes;
  psolve       = SPGMR_CONTENT(S)->Psolve;
  zeroguess    = &(SPGMR_CONTENT(S)->zeroguess);
  nli          = &(SPGMR_CONTENT(S)->numiters);
  res_norm     = &(SPGMR_CONTENT(S)->resnorm);
  cv           = SPGMR_CONTENT(S)->cv;
  Xv           = SPGMR_CONTENT(S)->Xv;

  static int firstTime = 1;
  if (firstTime) {
    setup(x, b, xcor, vtemp, s1, s2, V, l_max); 
    firstTime = 0; 
  } 

  /* Initialize counters and convergence flag */
  *nli = 0;
  converged = SUNFALSE;

  /* Check if Atimes function has been set */
  if (atimes == NULL) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_ATIMES_NULL;
    return(LASTFLAG(S));
  }

  /* cache s2Inv */
  platform->linAlg->adyz(N, ONE, o_s2, o_s2Inv);

  /* Set V[0] to initial (unscaled) residual r_0 = b - A*x_0 */
  if (*zeroguess) {
    o_vtemp.copyFrom(o_b, o_b.length());  
  } else {
    ier = atimes(A_data, x, vtemp);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
      return(LASTFLAG(S));
    }

    platform->linAlg->axpbyz(N, ONE, o_b, -ONE, o_vtemp, o_vtemp);
  }

  /* Set r_norm = beta to L2 norm of V[0] = s1 r_0, and return if small */
  platform->linAlg->axmy(N, ONE, o_s1, o_vtemp);

  realtype rdotr = platform->linAlg->innerProd(N, o_vtemp, o_vtemp, platform->comm.mpiComm); 
  *res_norm = r_norm = beta = SUNRsqrt(rdotr);

  if (r_norm <= delta) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_SUCCESS;
    return(LASTFLAG(S));
  }

  /* Initialize rho to avoid compiler warning message */
  rho = beta;

  /* Initialize the Hessenberg matrix Hes and Givens rotation product */
  for (i=0; i<=l_max; i++) {
    for (j=0; j<l_max; j++) {
      Hes[i][j] = ZERO;
    }
  }
  rotation_product = ONE;

  /* Set V[0] and normalize */
  axpby(ONE/r_norm, o_vtemp, ZERO, o_V);

  /* Inner loop: generate Krylov sequence and Arnoldi basis */
  for (l=0; l<l_max; l++) {
    (*nli)++;
    krydim = l_plus_1 = l + 1;

    /* Apply right scaling: vtemp = s2_inv V[l] */
    axmyz(ONE, o_V, o_s2Inv, o_vtemp, l);

    /* Apply A: V[l+1] = A s2_inv V[l] */
    ier = atimes(A_data, vtemp, vtemp2);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
      return(LASTFLAG(S));
    }

    /* Apply left scaling: V[l+1] = s1 A s2_inv V[l] */
    platform->linAlg->axmy(N, ONE, o_s1, o_vtemp2);

    /* Orthogonalize vtemp2 (V[l+1]) against previous V[i] */
    platform->timer.tic(timerName + "solve::cvode::linearSolve::cgs2", 0);
    CGS2(Hes, l_plus_1, l_max, cv, &(Hes[l_plus_1][l]), o_vtemp2);
    platform->timer.toc(timerName + "solve::cvode::linearSolve::cgs2");

    /* Update the QR factorization of Hes */
    if(SUNQRfact(krydim, Hes, givens, l) != 0 ) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = SUNLS_QRFACT_FAIL;
      return(LASTFLAG(S));
    }

    /* Update residual norm estimate; break if convergence test passes */
    rotation_product *= givens[2*l+1];
    *res_norm = rho = SUNRabs(rotation_product*r_norm);
    if (rho <= delta) { converged = SUNTRUE; break; }

    /* Set V[l+1] and normalize with norm value from the Gram-Schmidt routine */
    axpby(ONE/Hes[l_plus_1][l], o_vtemp2, ZERO, o_V, ZERO, l_plus_1);
  }
  /* Inner loop is done.  Compute the new correction vector xcor */

  /* Construct g, then solve for y */
  yg[0] = r_norm;
  for (i=1; i<=krydim; i++) yg[i]=ZERO;
  if (SUNQRsol(krydim, Hes, givens, yg) != 0) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_QRSOL_FAIL;
    return(LASTFLAG(S));
  }
  yg[krydim] = ZERO;

  /* Set xcor to  correction vector V_l y */
  o_yg.copyFrom(yg);
  linearCombination(krydim+1, o_yg, o_V, o_xcor, o_xcor);

  /* Add unscaled xcor to initial x to get final solution x, and return */
  if (converged || rho < beta) {
    if (*zeroguess) {
      platform->linAlg->axmyz(N, ONE, o_s2Inv, o_xcor, o_x);
    } else {
      platform->linAlg->axmy(N, ONE, o_s2Inv, o_xcor);
      platform->linAlg->axpby(N, ONE, o_xcor, ONE, o_x);
    }

    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = (converged) ? SUNLS_SUCCESS : SUNLS_RES_REDUCED;
    return(LASTFLAG(S));
  }

  *zeroguess  = SUNFALSE;
  LASTFLAG(S) = SUNLS_CONV_FAIL;
  return(LASTFLAG(S));
}

#endif
