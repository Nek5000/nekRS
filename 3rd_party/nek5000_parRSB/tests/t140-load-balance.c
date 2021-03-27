#include <genmap-impl.h>
#include <sort-impl.h>
#include <time.h>

typedef struct {
  double ds;
  uint proc;
} Data;

#define N 10

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  struct comm c;
  comm_init(&c, MPI_COMM_WORLD);

  srand(time(0));

  struct array arr;
  array_init(Data, &arr, N);

  int i, cnt;
  Data d;
  Data *ptr = arr.ptr;
  for (i = 0; i < N; i++)
    d.ds = (rand() % 100) / 100.0, ptr[i] = d;
  arr.n = N;

  parallel_sort(Data, &arr, ds, gs_double, 0, 1, &c);
  assert(arr.n == N);

  array_free(&arr);

  comm_free(&c);
  MPI_Finalize();

  return 0;
}
