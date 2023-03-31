#ifndef nrsbacktrace_hpp_
#define nrsbacktrace_hpp_

#include <execinfo.h>  // for backtrace
#include <dlfcn.h>     // for dladdr
#include <cxxabi.h>    // for __cxa_demangle

#include <string>
#include <sstream>

static std::string nrsbacktrace(int skip = 1)
{
  void *callstack[128];
  const int nMaxFrames = sizeof(callstack) / sizeof(callstack[0]);
  char buf[1024];
  int nFrames = backtrace(callstack, nMaxFrames);

  std::ostringstream trace_buf;
  for (int i = skip; i < nFrames; i++) {
    Dl_info info;
    if (dladdr(callstack[i], &info)) {
      char *demangled = NULL;
      int status = -1;
      demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);
      snprintf(buf, sizeof(buf), "%-3d %0*p %s + %zd\n",
           i, 2 + sizeof(void*) * 2, callstack[i],
           status == 0 ? demangled : info.dli_sname,
           (char *)callstack[i] - (char *)info.dli_saddr);
      free(demangled);
    } else {
      snprintf(buf, sizeof(buf), "%-3d %0*p\n",
           i, 2 + sizeof(void*) * 2, callstack[i]);
    }
    trace_buf << buf;
  }
  if (nFrames == nMaxFrames)
    trace_buf << "  [truncated]\n";
  return trace_buf.str();
}

#endif
