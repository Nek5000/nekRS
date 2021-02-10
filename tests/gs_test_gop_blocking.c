#include <assert.h>
#include "gslib.h"

struct gs_data *gop_handle;
int np;

//------------------------------------------------------------------------------
void gop_init(struct comm *gop_comm, comm_ext world) {
  comm_init(gop_comm, world);

  const long long gop_id = 1;

  gop_handle = gs_setup(&gop_id, 1, gop_comm, 0, gs_auto, 0);
}
//------------------------------------------------------------------------------
void gop(void *u, gs_dom dom, gs_op op, unsigned transpose) {
  gs(u, dom, op, transpose, gop_handle, NULL);
}
//------------------------------------------------------------------------------
void gop_free(struct comm* gop_comm) {
  comm_free(gop_comm);

  gs_free(gop_handle);
}
//------------------------------------------------------------------------------
int test_min(int rank) {
  int min = rank;
  gop(&min, gs_int, gs_min, 0);

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
int test_max(int rank) {
  int max = rank;
  gop(&max, gs_int, gs_max, 0);

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
int test_add(int rank) {
  int sum = rank;
  gop(&sum, gs_int, gs_add, 0);
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

  result  = test_min(rank);
  result += test_max(rank);
  result += test_add(rank);

  gop_free(&comm);

#ifdef MPI
  MPI_Finalize();
#endif

  return result;
}
