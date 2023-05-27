#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#if defined(PARRSB_UNDERSCORE)
#define FNAME(x) TOKEN_PASTE(x, _)
#else
#define FNAME(x) x
#endif

#define FDGETRF FNAME(dgetrf)
#define FDGETRI FNAME(dgetri)

#if defined(PARRSB_BLAS)

void FDGETRF(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
void FDGETRI(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork,
             int *INFO);

void matrix_inverse(int N, double *A) {
  int size = N * N;
  int info;

  int *ipiv = (int *)calloc(N, sizeof(int));
  double *work = (double *)calloc(N * N, sizeof(double));

  FDGETRF(&N, &N, A, &N, ipiv, &info);
  if (info != 0)
    printf("dgetrf: %d\n", info);

  FDGETRI(&N, A, &N, ipiv, work, &size, &info);
  if (info != 0)
    printf("dgetri: %d\n", info);

  free(ipiv);
  free(work);
}

#else

void matrix_inverse(int N, double *A) {}

#endif // PARRSB_BLAS

#undef FDGETRF
#undef FDGETRI
#undef FNAME
