#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "sort.h"

#if 1

#define N (1<<20)

ulong A[N], out[N];
uint P[N];

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

  buffer_free(&buf);
  return 0;
}

#else

int main()
{
  return 0;
}

#endif

