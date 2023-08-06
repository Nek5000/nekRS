/*
 *
 * cblas_dger.c
 * This program is a C interface to dger.
 * Written by Keita Teranishi
 * 4/6/1998
 *
 */

#include "cblas.h"
#include "cblas_f77.h"
void cblas_dger(const CBLAS_LAYOUT layout, const CBLAS_INT M, const CBLAS_INT N,
                const double alpha, const double  *X, const CBLAS_INT incX,
                const double  *Y, const CBLAS_INT incY, double  *A, const CBLAS_INT lda)
{
#ifdef F77_INT
   F77_INT F77_M=M, F77_N=N, F77_lda=lda, F77_incX=incX, F77_incY=incY;
#else
   #define F77_M M
   #define F77_N N
   #define F77_incX incX
   #define F77_incY incY
   #define F77_lda lda
#endif

   extern int CBLAS_CallFromC;
   extern int RowMajorStrg;
   RowMajorStrg = 0;

   CBLAS_CallFromC = 1;
   if (layout == CblasColMajor)
   {
      F77_dger( &F77_M, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A,
                      &F77_lda);
   }
   else if (layout == CblasRowMajor)
   {
      RowMajorStrg = 1;
      F77_dger( &F77_N, &F77_M ,&alpha, Y, &F77_incY, X, &F77_incX, A,
                      &F77_lda);

   }
   else cblas_xerbla(1, "cblas_dger", "Illegal layout setting, %d\n", layout);
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}
