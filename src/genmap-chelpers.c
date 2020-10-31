#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <string.h>

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
  for(i = 0; i < bt_size; i++) printf("%s\n", symbols[i]);
  free(symbols);
}
#else
void GenmapPrintStack() {};
#endif

double GenmapGetMaxRss() {
  struct rusage r_usage;

  getrusage(RUSAGE_SELF, &r_usage);
#if defined(__APPLE__) && defined(__MACH__)
  return (double)r_usage.ru_maxrss;
#else
  return (double)(r_usage.ru_maxrss * 1024L);
#endif
}

void set_stdout(char *f, int *sid, int flen) {
  char *logfile = (char *) malloc((flen + 2 + 5 + 1) * sizeof(char));
  strncpy(logfile, f, flen);
  int i;
  for(i = flen - 1; i >= 0; i--) if(logfile[i] != ' ') break;
  logfile[i + 1] = '\0';

  int redirect = 0;
  char *envvar;

  if(logfile[0] != '\0') {
    redirect = 1;
  } else if((envvar = getenv("NEK_LOGFILE")) != NULL) {
    if(*sid >= 0) sprintf(logfile, "s%05d_", *sid);
    strcat(logfile + strlen(logfile), envvar);
    redirect = 1;
  }

  if(redirect) {
    printf("redirecting stdout to %s\n", logfile);
    FILE *fp = freopen(logfile, "w+", stdout);
    assert(fp != NULL);
  }
}
