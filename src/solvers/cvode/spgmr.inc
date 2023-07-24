/*
 * -----------------------------------------------------------------
 * Clone of CVODE's SPGMR solver with some modifications;
 * - no restarts
 * - no preconditioning
 * - always scale and cache s2Inv 
 * - use CGS2 instead of CGS
 * -----------------------------------------------------------------
 */


#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_math.h>

#include "sundials_context_impl.h"
#include "sundials_logger_impl.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

#define SPGMR_CONTENT(S)  ( (SUNLinearSolverContent_SPGMR)(S->content) )
#define LASTFLAG(S)       ( SPGMR_CONTENT(S)->last_flag )

/* two pass (iterative) classical Gram-Schmidt */
int CGS2(N_Vector *v, realtype **h, int k, int p, realtype *new_vk_norm,
         realtype *stemp, N_Vector *vtemp)
{
  int  i, i0, k_minus_1, retval;
  realtype vk_norm;

  k_minus_1 = k - 1;
  i0 = SUNMAX(k-p,0);

  retval = N_VDotProdMulti(k-i0+1, v[k], v+i0, stemp);
  if (retval != 0) return(-1);

  vk_norm = SUNRsqrt(stemp[k-i0]);
  for (i=k-i0-1; i >= 0; i--) {
    h[i][k_minus_1] = stemp[i];
    stemp[i+1] = -stemp[i];
    vtemp[i+1] = v[i];
  }
  stemp[0] = ONE;
  vtemp[0] = v[k];

  retval = N_VLinearCombination(k-i0+1, stemp, vtemp, v[k]);
  if (retval != 0) return(-1);

  /* re-orthogonalize */
  {
    retval = N_VDotProdMulti(k-i0, v[k], v+i0, stemp+1);
    if (retval != 0) return(-1);

    stemp[0] = ONE;
    vtemp[0] = v[k];
    for (i=i0; i < k; i++) {
      h[i][k_minus_1] += stemp[i-i0+1];
      stemp[i-i0+1] = -stemp[i-i0+1];
      vtemp[i-i0+1] = v[i-i0];
    }

    retval = N_VLinearCombination(k+1, stemp, vtemp, v[k]);
    if (retval != 0) return(-1);
  }

  *new_vk_norm = SUNRsqrt(N_VDotProd(v[k],v[k]));

  return(0);
}

int SPGMR(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype delta)
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

  static int firstTime = 1;

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

  static N_Vector s2Inv;
  if (firstTime) s2Inv = N_VClone(s2); 

  A_data       = SPGMR_CONTENT(S)->ATData;
  P_data       = SPGMR_CONTENT(S)->PData;
  atimes       = SPGMR_CONTENT(S)->ATimes;
  psolve       = SPGMR_CONTENT(S)->Psolve;
  zeroguess    = &(SPGMR_CONTENT(S)->zeroguess);
  nli          = &(SPGMR_CONTENT(S)->numiters);
  res_norm     = &(SPGMR_CONTENT(S)->resnorm);
  cv           = SPGMR_CONTENT(S)->cv;
  Xv           = SPGMR_CONTENT(S)->Xv;

#if 0
  auto o_V = getN_VectorMemory(sunrealtype, V);
  auto o_xcor = getN_VectorMemory(sunrealtype, xcor);
  auto o_vtemp = getN_VectorMemory(sunrealtype, vtemp);
  auto o_s1 = getN_VectorMemory(sunrealtype, s1);
  auto o_s2 = getN_VectorMemory(sunrealtype, s2);
#endif


  /* Initialize counters and convergence flag */
  firstTime = 0; 
  *nli = 0;
  converged = SUNFALSE;

  /* Check if Atimes function has been set */
  if (atimes == NULL) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_ATIMES_NULL;
    return(LASTFLAG(S));
  }

  /* Set vtemp to initial (unscaled) residual r_0 = b - A*x_0 */
  if (*zeroguess) {
    N_VScale(ONE, b, vtemp);
  } else {
    ier = atimes(A_data, x, vtemp);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
      return(LASTFLAG(S));
    }
    N_VLinearSum(ONE, b, -ONE, vtemp, vtemp);
  }

  /* Set r_norm = beta to L2 norm of V[0] = s1 r_0, and return if small */
  N_VProd(s1, vtemp, V[0]);
  *res_norm = r_norm = beta = SUNRsqrt(N_VDotProd(V[0], V[0]));
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

  /* Normalize the initial vector V[0] */
  N_VScale(ONE/r_norm, V[0], V[0]);

  /* cache s2Inv */
  N_VInv(s2, s2Inv);

  /* Inner loop: generate Krylov sequence and Arnoldi basis */
  for (l=0; l<l_max; l++) {
    (*nli)++;
    krydim = l_plus_1 = l + 1;

    /* Apply right scaling: vtemp = s2_inv V[l] */
    N_VProd(V[l], s2Inv, vtemp);

    /* Apply A: V[l+1] = A s2_inv V[l] */
    ier = atimes( A_data, vtemp, V[l_plus_1] );
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
      return(LASTFLAG(S));
    }

    /* Apply left scaling: V[l+1] = s1 A s2_inv V[l] */
    N_VProd(s1, V[l_plus_1], V[l_plus_1]);

    /* Orthogonalize V[l+1] against previous V[i]: V[l+1] = w_tilde */
    if (gstype == SUN_CLASSICAL_GS) {
      if (CGS2(V, Hes, l_plus_1, l_max, &(Hes[l_plus_1][l]), cv, Xv) != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = SUNLS_GS_FAIL;
        return(LASTFLAG(S));
      }
    } else {
      if (SUNModifiedGS(V, Hes, l_plus_1, l_max, &(Hes[l_plus_1][l])) != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = SUNLS_GS_FAIL;
        return(LASTFLAG(S));
      }
    }

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

    /* Normalize V[l+1] with norm value from the Gram-Schmidt routine */
    N_VScale(ONE/Hes[l_plus_1][l], V[l_plus_1], V[l_plus_1]);
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

  /* Add correction vector V_l y to xcor */
  for (k=0; k<krydim; k++) {
    cv[k] = yg[k];
    Xv[k] = V[k];
  }
  ier = N_VLinearCombination(krydim, cv, Xv, xcor);
  if (ier != SUNLS_SUCCESS) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_VECTOROP_ERR;
    return(SUNLS_VECTOROP_ERR);
  }

  /* Add unscaled xcor to initial x to get final solution x, and return */
  if (converged || rho < beta) {
    if (*zeroguess) {
      N_VProd(s2Inv, xcor, x);
    } else {
      N_VProd(xcor, s2Inv, xcor);
      N_VLinearSum(ONE, x, ONE, xcor, x);
    }

    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = (converged) ? SUNLS_SUCCESS : SUNLS_RES_REDUCED;
    return(LASTFLAG(S));
  }

  *zeroguess  = SUNFALSE;
  LASTFLAG(S) = SUNLS_CONV_FAIL;
  return(LASTFLAG(S));
}
