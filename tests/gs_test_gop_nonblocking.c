#include <assert.h>
#include "gslib.h"

struct gs_data *gop_handle;
int np;

//------------------------------------------------------------------------------
void gop_init(struct comm *gop_comm, comm_ext world) {
  comm_init(gop_comm, world);

  const long long gop_id = 1;

  gop_handle = gs_setup(&gop_id, 1, gop_comm, 0, gs_pairwise, 0);
}
//------------------------------------------------------------------------------
void igop(void *u, gs_dom dom, gs_op op, unsigned transpose) {
  // In a real case, these calls will be split across other code
  int handle;
  igs(u, dom, op, transpose, gop_handle, NULL, &handle);
  gs_wait (handle);
}
//------------------------------------------------------------------------------
void gop_free(struct comm* gop_comm) {
  comm_free(gop_comm);

  gs_free(gop_handle);
}
//------------------------------------------------------------------------------
int test_imin(int rank) {
  int min = rank;
  igop(&min, gs_int, gs_min, 0);

  if (rank == 0) printf("\ngop min test: ");
  if (min == 0) {
    if (rank == 0) printf("[Passed]");
    return 0;
  } else {
    if (rank == 0) printf("[Failed]");
    return 1;
  }
}
//------------------------------------------------------------------------------
int test_imax(int rank) {
  int max = rank;
  igop(&max, gs_int, gs_max, 0);

  if (rank == 0) printf("\ngop max test: ");
  if (max == np-1) {
    if (rank == 0) printf("[Passed]");
    return 0;
  } else {
    if (rank == 0) printf("[Failed]");
    return 1;
  }
}
//------------------------------------------------------------------------------
int test_iadd(int rank) {
  int sum = rank;
  igop(&sum, gs_int, gs_add, 0);
  sum *= 2;

  if (rank == 0) printf("\ngop add test: ");
  if (sum == np*(np-1)) {
    if (rank == 0) printf("[Passed]");
    return 0;
  } else {
    if (rank == 0) printf("[Failed]");
    return 1;
  }
}
//------------------------------------------------------------------------------
int test_imul(int rank) {
  int mul = rank + 1;
  igop(&mul, gs_int, gs_mul, 0);

  int answer=1, i;
  for(i = 2; i <= np; i++) {
    answer*=i;
  }
  if (rank == 0) printf("\ngop mul test: ");
  if (mul == answer) {
    if (rank == 0) printf("[Passed]");
    return 0;
  } else {
    if (rank == 0) printf("[Failed]");
    return 1;
  }
}
//------------------------------------------------------------------------------
int main(int narg, char *arg[])
{
  comm_ext world; int rank, result;
  struct comm comm;

#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);
  MPI_Comm_rank(world,&rank);
#else
  world=0, np=1; rank = 0;
#endif

  gop_init(&comm,world);

  result  = test_imin(rank);
  result += test_imax(rank);
  result += test_iadd(rank);
  result += test_imul(rank);

  gop_free(&comm);

  if (rank == 0) printf("\n");

#ifdef MPI
  MPI_Finalize();
#endif

  return result;
}
