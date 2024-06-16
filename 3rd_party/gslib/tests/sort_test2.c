#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "gslib.h"

#if 1

#define N (1<<20)

#define TYPE_INT    0
#define TYPE_LONG   1
#define TYPE_FLOAT  2
#define TYPE_DOUBLE 3

#define PI acos(-1.0)
#define SIZE 20

typedef double T;
struct test_struct {
  T data;
  uint i;
};
typedef struct test_struct S;

ulong A[N], out[N];
uint P[N];

int test_real_sort(buffer *buf)
{
  S u[SIZE];

  u[0].data=1.0,u[0].i=1;
  sarray_sort(S,u,1,data,TYPE_DOUBLE,buf);
  assert(fabs(u[0].data - 1.0) < 1e-16);

  u[1].data=-1.0,u[1].i=2;
  sarray_sort(S,u,2,data,TYPE_DOUBLE,buf);
  assert(fabs(u[0].data + 1.0) < 1e-16 && u[0].i == 2);
  assert(fabs(u[1].data - 1.0) < 1e-16 && u[1].i == 1);

  uint i;
  for(i=0;i<SIZE;i++) {
   u[i].data=sin(2*PI*i/(double)SIZE);
   u[i].i=i;
  }
  sarray_sort(S,u,SIZE,data,TYPE_DOUBLE,buf);
  for(i=1;i<SIZE;i++)
    assert(u[i-1].data <= u[i].data);

  return 0;
}

int main()
{
  buffer buf = null_buffer;
  uint i;
  unsigned long long tic, toc;
  unsigned r;

  for(i=0;i!=N;++i) {
    A[i]=rand();
    A[i]<<=CHAR_BIT*sizeof(int)-1;
    A[i]^=rand();
    A[i]<<=CHAR_BIT*sizeof(int)-1;
    A[i]^=rand();
    if(0) A[i]&=0x000ff00;
  }

  for(i=N;i;i>>=1) {
    unsigned long long t;
    sortv_long(out, A,i,sizeof(ulong), &buf);
  }

  for(i=N;i;i>>=1) {
    unsigned long long t;
    sortp_long(&buf,0, A,i,sizeof(ulong));
  }

  // Test local real-sort implementation
  int result = 1;
  result = test_real_sort(&buf);

  buffer_free(&buf);

  return result;
}

#else

int main()
{
  return 0;
}

#endif
