#include <sys/types.h>
#if defined(FUNCPROTO) || defined(__STDC__) || defined(__cplusplus) || defined(c_plusplus)
#ifndef	NULL
#define NULL	((void *) 0)
#endif
#else
#ifndef	NULL
#define NULL	0
#endif
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

extern char *getlogin();
#ifdef _MSC_VER
    #define strdup _strdup
    #include <io.h>
#pragma warning(disable: 4996)
#endif
