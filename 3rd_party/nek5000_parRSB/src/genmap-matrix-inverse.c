#include <genmap-impl.h>

#if defined(GENMAP_UNDERSCORE)
#  define FNAME(x) TOKEN_PASTE(x,_)
#else
#  define FNAME(x) x
#endif

#if defined(GENMAP_BLAS)

#define FDGETRF FNAME(dgetrf)
void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
#define FDGETRI FNAME(dgetri)
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

void matrix_inverse(int N,double *A){
  int size=N*N;
  int info;

  int *ipiv=(int*) calloc(N,sizeof(int));
  double *work=(double *) calloc(N*N,sizeof(double));

  FDGETRF(&N,&N,A,&N,ipiv,&info);
  if(info!=0)
    printf("dgetrf: %d\n",info);

  FDGETRI(&N,A,&N,ipiv,work,&size,&info);
  if(info!=0)
    printf("dgetri: %d\n",info);

  free(ipiv);
  free(work);
}

#undef FDGETRF
#undef FDGETRI

#endif // GENMAP_BLAS

#undef FNAME
