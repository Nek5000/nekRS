#ifndef _CM_CONFIG_H
#define _CM_CONFIG_H

/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 if you have the <enet/enet.h> header file. */
#cmakedefine HAVE_ENET_ENET_H

/* Define to desired port range Low:High */
#define EVPATH_DEFAULT_PORT_RANGE "@EVPATH_DEFAULT_PORT_RANGE@"

/* Define to 1 if you have the enet header and libraries */
#cmakedefine ENET_FOUND

/* Define to 1 to if ZPL ENET transport was built */
#cmakedefine ZPL_ENET_AVAILABLE

/* Define to 1 if you have the udt4 header and libraries */
#cmakedefine UDT4_FOUND

/* Define to 1 if you have the IB header and libraries */
#cmakedefine IB_FOUND 

/* Define to 1 if you have the LIBFABRIC header and libraries */
#cmakedefine LIBFABRIC_FOUND 

/* Define to 1 if you have the df_shm header and libraries */
#cmakedefine DF_SHM_FOUND

/* Place where evpath transport modules are installed */
#cmakedefine EVPATH_MODULE_INSTALL_DIR "@EVPATH_MODULE_INSTALL_DIR@"

/* Place where evpath tests are to be installed */
#cmakedefine EVPATH_TEST_INSTALL_DIR "@EVPATH_TEST_INSTALL_DIR@"

/* Place where the ssh binary can be found */
#define SSH_PATH "@SSH@"

/* Define to 1 if you have the <dlfcn.h> header file. */
#cmakedefine HAVE_DLFCN_H

/* Define if fd_set has element fds_bits */
#cmakedefine HAVE_FDS_BITS

/* Define to 1 if you have the `getdomainname' function. */
#cmakedefine HAVE_GETDOMAINNAME

/* Define to 1 if you have the `getloadavg' function. */
#cmakedefine HAVE_GETLOADAVG

/* Define to 1 if you have the `gettimeofday' function. */
#cmakedefine HAVE_GETTIMEOFDAY

/* Define to 1 if you have the `pthread' library (-lpthread). */
#cmakedefine HAVE_LIBPTHREAD

/* Uses the Mac OS X style conventions for sysctl */
#cmakedefine HAVE_MAC_SYSCTL

/* Define if system uses LINUX-style sysctl */
#cmakedefine HAVE_LINUX_SYSCTL 

/* Define to 1 if you have the <malloc.h> header file. */
#cmakedefine HAVE_MALLOC_H

/* Define to 1 if you have the <math.h> header file. */
#cmakedefine HAVE_MATH_H

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H

/* Define to 1 if you have the <netdb.h> header file. */
#cmakedefine HAVE_NETDB_H

/* Define to 1 if you have the <stdarg.h> header file. */
#cmakedefine HAVE_STDARG_H

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H

/* Define to 1 if you have the `sysconf' function. */
#cmakedefine HAVE_SYSCONF

/* Define to 1 if you have the `sysinfo' function. */
#cmakedefine HAVE_SYSINFO

/* Define to 1 if you have the <sys/epoll.h> header file. */
#cmakedefine HAVE_SYS_EPOLL_H

/* Define to 1 if you have the <cod.h> header file. */
#cmakedefine HAVE_COD_H

/* Define to 1 if you have the <sys/select.h> header file. */
#cmakedefine HAVE_SYS_SELECT_H

/* Define to 1 if you have the <sys/sockio.h> header file. */
#cmakedefine HAVE_SYS_SOCKIO_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/sysctl.h> header file. */
#cmakedefine HAVE_SYS_SYSCTL_H

/* Define to 1 if you have the <sys/times.h> header file. */
#cmakedefine HAVE_SYS_TIMES_H

/* Define to 1 if you have the <sys/time.h> header file. */
#cmakedefine HAVE_SYS_TIME_H

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H

/* Define to 1 if you have the <sys/uio.h> header file. */
#cmakedefine HAVE_SYS_UIO_H

/* Define to 1 if you have the <sys/un.h> header file. */
#cmakedefine HAVE_SYS_UN_H

/* Define to 1 if you have the `uname' function. */
#cmakedefine HAVE_UNAME

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H

/* Define to 1 if you have the <windows.h> header file. */
#cmakedefine HAVE_WINDOWS_H

/* Define to 1 if you have the <winsock2.h> header file. */
#cmakedefine HAVE_WINSOCK2_H

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT @SIZEOF_INT@

/* The size of `long', as computed by sizeof. */
#cmakedefine SIZEOF_LONG @SIZEOF_LONG@

/* Define to 1 if you have the `writev' function. */
#cmakedefine HAVE_WRITEV

/* Define to 1 if you have the `getifaddrs' function. */
#cmakedefine HAVE_GETIFADDRS

/* Define to 1 if you have the `clock_gettime' function. */
#cmakedefine HAVE_CLOCK_GETTIME

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS

/* Define this if Pthreads should be used for running tests */
#cmakedefine USE_PTHREADS

/* Version number of package */
#cmakedefine VERSION

/* Define so that glibc/gnulib argp.h does not typedef error_t. */
#cmakedefine __error_t_defined

/* Define to empty if `const' does not conform to ANSI C. */
#cmakedefine const

/* Define to a type to use for `error_t' if it is not otherwise available. */
#cmakedefine error_t

/* Define to `int' if <sys/types.h> does not define. */
#cmakedefine pid_t

/* Define to `int' if <sys/types.h> does not define. */
#cmakedefine pid_t

/* Define to the shared module suffix in use */
#cmakedefine CMAKE_SHARED_MODULE_SUFFIX "@CMAKE_SHARED_MODULE_SUFFIX@"

/* Set to 1 to build without dynamic linking  */
#cmakedefine NO_DYNAMIC_LINKING @NO_DYNAMIC_LINKING@

/* Set to 1 if NNTI libraries and include file are found */
#cmakedefine NNTI_FOUND

/* Set to 1 if NVML libraries and include file are found */
#cmakedefine NVML_FOUND

/* Set to 1 if CM should default to CMSelfFormats */
#define CM_SELF_FORMATS @CM_SELF_FORMATS@

/* Set to 1 if Cmake found MPI for C */
#cmakedefine MPI_C_FOUND

/* Define if byteorder is bigendian */
#cmakedefine WORDS_BIGENDIAN

#define CM_DEFAULT_TRANSPORT "@CM_DEFAULT_TRANSPORT@"

#define CM_LIBRARY_PREFIX "@EVPATH_LIBRARY_PREFIX@"

#define IPCONFIG_ENVVAR_PREFIX "@IPCONFIG_ENVVAR_PREFIX@"

#endif
