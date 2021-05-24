#ifndef BOOMERAMG_H
#define BOOMERAMG_H

#define BOOMERAMG_NPARAM 9

#ifdef __cplusplus
extern "C" {
#endif

int boomerAMGSetup(int nrows,
                   int nz, const long long int *Ai, const long long int *Aj, const double *Av,
                   const int null_space, const MPI_Comm ce, int Nthreads, int deviceID,
                   const double *param);

int boomerAMGSolve(void *x, void *b);

void boomerAMGFree();

#ifdef __cplusplus
}
#endif

#endif
