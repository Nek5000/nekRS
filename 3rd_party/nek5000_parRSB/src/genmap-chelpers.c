#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>

#include <genmap.h>

#if defined __GLIBC__

#include <execinfo.h>

/* Obtain a backtrace and print it to stdout. */
void GenmapPrintStack(void) {
  void *bt[50];
  int i;
  int bt_size = backtrace(bt, 50);
  char **symbols = backtrace_symbols(bt, bt_size);

  printf("backtrace(): obtained %d stack frames.\n", bt_size);
  for (i = 0; i < bt_size; i++)
    printf("%s\n", symbols[i]);
  free(symbols);
}
#else
void GenmapPrintStack(){};
#endif // defined __GLIBC__

double GenmapGetMaxRss() {
  struct rusage r_usage;

  getrusage(RUSAGE_SELF, &r_usage);
#if defined(__APPLE__) && defined(__MACH__)
  return (double)r_usage.ru_maxrss;
#else
  return (double)(r_usage.ru_maxrss * 1024L);
#endif
}

int log2ll(long long n) {
  int k = 0;
  while (n > 1)
    n /= 2, k++;

  return k;
}
