/*
 *
 * cblas_sgbmv.c
 * This program is a C interface to sgbmv.
 * Written by Keita Teranishi
 * 4/6/1998
 *
 */

#include "cblas.h"
#include "cblas_f77.h"
void cblas_sgbmv(const CBLAS_LAYOUT layout,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_INT M, const CBLAS_INT N,
                 const CBLAS_INT KL, const CBLAS_INT KU,
                 const float alpha, const float *A, const CBLAS_INT lda,
                 const float  *X, const CBLAS_INT incX, const float beta,
                 float  *Y, const CBLAS_INT incY)
{
   char TA;
#ifdef F77_CHAR
   F77_CHAR F77_TA;
#else
   #define F77_TA &TA
#endif
#ifdef F77_INT
   F77_INT F77_M=M, F77_N=N, F77_lda=lda, F77_incX=incX, F77_incY=incY;
   F77_INT F77_KL=KL,F77_KU=KU;
#else
   #define F77_M M
   #define F77_N N
   #define F77_lda lda
   #define F77_KL KL
   #define F77_KU KU
   #define F77_incX incX
   #define F77_incY incY
#endif
   extern int CBLAS_CallFromC;
   extern int RowMajorStrg;
   RowMajorStrg = 0;

   CBLAS_CallFromC = 1;
   if (layout == CblasColMajor)
   {
      if (TransA == CblasNoTrans) TA = 'N';
      else if (TransA == CblasTrans) TA = 'T';
      else if (TransA == CblasConjTrans) TA = 'C';
      else
      {
         cblas_xerbla(2, "cblas_sgbmv","Illegal TransA setting, %d\n", TransA);
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      #ifdef F77_CHAR
         F77_TA = C2F_CHAR(&TA);
      #endif
      F77_sgbmv(F77_TA, &F77_M, &F77_N, &F77_KL, &F77_KU, &alpha,
                     A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY);
   }
   else if (layout == CblasRowMajor)
   {
      RowMajorStrg = 1;
      if (TransA == CblasNoTrans) TA = 'T';
      else if (TransA == CblasTrans) TA = 'N';
      else if (TransA == CblasConjTrans) TA = 'N';
      else
      {
         cblas_xerbla(2, "cblas_sgbmv","Illegal TransA setting, %d\n", TransA);
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      #ifdef F77_CHAR
         F77_TA = C2F_CHAR(&TA);
      #endif
      F77_sgbmv(F77_TA, &F77_N, &F77_M, &F77_KU, &F77_KL, &alpha,
                     A ,&F77_lda, X, &F77_incX, &beta, Y, &F77_incY);
   }
   else cblas_xerbla(1, "cblas_sgbmv", "Illegal layout setting, %d\n", layout);
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}
