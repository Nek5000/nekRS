#ifndef BOOMERAMG_H
#define BOOMERAMG_H

#define BOOMERAMG_NPARAM 10

#ifdef __cplusplus
extern "C" {
#endif

int boomerAMGSetup(int nrows,
                   int nz, const long long int *Ai, const long long int *Aj, const double *Av,
                   const int null_space, const MPI_Comm ce, int Nthreads, int deviceID,
                   const int useFP32, const double *param, const int verbose);

int boomerAMGSolve(void *x, void *b);

void boomerAMGFree();

#ifdef __cplusplus
}
#endif

#endif
