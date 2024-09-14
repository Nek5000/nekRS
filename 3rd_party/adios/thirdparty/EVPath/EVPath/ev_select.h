#ifndef EVPATH_SELECT_H
#define EVPATH_SELECT_H

#ifdef HAVE_SYS_SELECT_H
#include <sys/select.h>
#endif

#if defined(__has_feature)
#if __has_feature(memory_sanitizer)
#define USE_MEMSET
#endif
#endif
#if defined(__NVCOMPILER)
#define USE_MEMSET 1
#endif
#ifdef USE_MEMSET
#include <string.h>
#define EVPATH_FD_ZERO(set) \
do { \
    memset((set), 0, sizeof(fd_set)); \
} while (0)
#endif

#ifndef EVPATH_FD_ZERO
#define EVPATH_FD_ZERO FD_ZERO
#endif

#endif
