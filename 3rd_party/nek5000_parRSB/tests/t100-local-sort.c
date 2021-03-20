#include <genmap-impl.h>
#include <sort-impl.h>
#include <time.h>

typedef struct {
  double ds;
  ulong dl;
  sint proc;
} Data;

#define N 10

int main(int argc, char *argv[]) {
  struct array arr;
  array_init(Data, &arr, N);

  srand(time(0));

  int i, cnt;
  Data d;
  Data *ptr = arr.ptr;
  for (i = cnt = 0; i < N / 2; i++) {
    d.dl = rand() % 100, d.ds = (rand() % 100) / 100.0, ptr[cnt++] = d;
    d.dl = rand() % 100, ptr[cnt++] = d;
  }
  arr.n = N;

  struct sort sd;
  sd.align = ALIGNOF(Data), sd.unit_size = sizeof(Data);
  sd.nfields = 2;
  sd.offset[0] = offsetof(Data, ds), sd.t[0] = gs_double;
  sd.offset[1] = offsetof(Data, dl), sd.t[1] = gs_long;
  sd.a = &arr;
  sd.balance = 1;
  sd.algo = bin_sort;
  buffer_init(&sd.buf, 1024);

  sort_local(&sd);

  buffer_free(&sd.buf);

  ptr = arr.ptr;
  for (i = 0; i < N - 1; i++) {
    assert(ptr[i].ds <= ptr[i + 1].ds && "Field ds is not sorted");
    if (i % 2 == 0)
      assert(ptr[i].dl <= ptr[i + 1].dl && "Field dl is not sorted");
  }

  array_free(&arr);
  return 0;
}
