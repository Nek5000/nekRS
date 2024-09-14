#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#undef FD_SETSIZE
#define FD_SETSIZE 1024
#include <winsock.h>
#define getpid()	_getpid()
#define getcwd(a,b)	_getcwd(a,b)
#define __ANSI_CPP__
#else
#include <unistd.h>
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif
#include <sys/socket.h>
#ifdef HAVE_SYS_SELECT_H
#include <sys/select.h>
#endif
#ifdef HAVE_SYS_UIO_H
#include <sys/uio.h>
#define HAVE_IOVEC_DEFINE
#endif
#ifdef HAVE_SYS_UN_H
#include <sys/un.h>
#endif
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#include <signal.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <signal.h>
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif

#include "fm.h"
#include "fm_internal.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>
#include <assert.h>

static void
die_with_error(int sig)
{
    fprintf(stderr, "Format server not responding.  Timeout.\n");
    exit(1);
}

extern int serverAtomicWrite(void *fd, const void *buffer, int length);

static int slave = 0;

static void check_for_running_server_and_fork();
static int quiet = 0;
static int do_proxy = 0;
static int no_fork = 0;
static int no_log = 0;

extern void ffs_set_no_log(int in_no_log);
int
main(int argc, char **argv)
{
    int i;
    int fs_port = DEFAULT_FS_PORT;
    int do_restart = 0;
    int launchd_mode = 0;

/*    if (os_sockets_init_func != NULL) os_sockets_init_func();*/

    if (strcmp(argv[0] + strlen(argv[0]) - 5, "proxy") == 0) {
	/* if the command name ends in "proxy", assume we're should be a proxy */
	do_proxy++;
    }
    for (i = 1; i < argc; i++) {
	if (strcmp(argv[i], "-no_fork") == 0) {
	    no_fork++;
	} else if (strcmp(argv[i], "-no_log") == 0) {
	    no_log++;
	} else if (strcmp(argv[i], "-quiet") == 0) {
	    quiet++;
	} else if (strcmp(argv[i], "-slave") == 0) {
	    slave++;
	} else if (strcmp(argv[i], "-restart") == 0) {
	    do_restart++;
	} else if (strcmp(argv[i], "-proxy") == 0) {
	    do_proxy++;
	} else if (strcmp(argv[i], "-launchd") == 0) {
	    launchd_mode++;
	    do_proxy++;
	} else {
	    fprintf(stderr, "Unknown argument \"%s\"\n", argv[i]);
	    fprintf(stderr, "Usage:  format_server [-no_fork] [-quiet] [-restart]\n");
	    exit(1);
	}
    }

    if (do_proxy && !no_fork) quiet++;

    if (!launchd_mode) check_for_running_server_and_fork();

    if (!launchd_mode) {
#ifdef RLIMIT_NOFILE
	struct rlimit lim;
	if (getrlimit(RLIMIT_NOFILE, &lim) != 0) {
	    perror("Getrlimit");
	}
	lim.rlim_cur = lim.rlim_max;
	if (setrlimit(RLIMIT_NOFILE, &lim) != 0) {
#if defined(__APPLE__) && defined(RLIMIT_NOFILE) && defined(OPEN_MAX)
	  if (lim.rlim_cur > OPEN_MAX) {
	    lim.rlim_cur = OPEN_MAX;
	    if (setrlimit(RLIMIT_NOFILE, &lim) != 0) {
	      perror("Setrlimit");
	    }
	  }
#else
          perror("Setrlimit");
#endif
        }

	if (getrlimit(RLIMIT_CORE, &lim) != 0) {
	    perror("Setrlimit");
	}
	lim.rlim_cur = lim.rlim_max;
	if (setrlimit(RLIMIT_CORE, &lim) != 0) {
	    perror("Setrlimit");
	}
	chdir("/tmp");
#endif
    }
    if (getenv("FORMAT_SERVER_PORT") != NULL) {
	char *port_string = getenv("FORMAT_SERVER_PORT");
	int tmp_port;
	if (sscanf(port_string, "%d", &tmp_port) != 1) {
	    printf("FORMAT_SERVER_PORT spec \"%s\" not understood.\n", 
		   port_string);
	} else {
	    fs_port = tmp_port;
	}
    }
    ffs_set_no_log(no_log);
    general_format_server(fs_port, do_restart, no_fork, do_proxy);
    return 0;
}

void check_for_running_server_and_fork()
{
    char testing_char = 'T';
    /* test to see if format server is running */
    FMContext test = create_FMcontext();
#ifndef HAVE_WINDOWS_H
    signal(SIGALRM, die_with_error);
    alarm(30);
#endif
    if (!quiet) {
	putenv(strdup("FORMAT_SERVER_VERBOSE=1"));
    }
    if (do_proxy) {
	establish_server_connection(test, local_only);
    } else {
	establish_server_connection(test, host_only);
    }
    if (serverAtomicWrite(test->server_fd, &testing_char, 1)
	== 1) {
	/* already running */
	if (!quiet) {
	    printf("Use FORMAT_SERVER_HOST environment variable to change format server location\n");
	    printf("Format server already running.  Exiting\n");
	}
	exit(0);
    } else {
	char format_server_pid[] = "/tmp/format_server_pid";

	/* make really, really certain there's no format_server running */
	FILE *format_server_pid_file = fopen(format_server_pid, "r");
	alarm(0);
	if (format_server_pid_file != NULL) {
	    long server_pid;
	    if (fscanf(format_server_pid_file, "%lx", &server_pid) == 1) {
		if (kill(server_pid, 0) == 0) {
		    fprintf(stderr, "Format Server %lx not responding, but still running\n",
			    server_pid);
		    exit(0);
		}
	    }
	    fclose(format_server_pid_file);
	}
	if (!no_fork) {
	    if (fork() != 0) {
		/* I'm the parent, return now */
		exit(0);
	    }
	}
	signal(SIGPIPE, SIG_IGN);
	format_server_pid_file = fopen(format_server_pid, "w");
	if (format_server_pid_file) {
	    fprintf(format_server_pid_file, "%lx\n", (long) getpid());
	    fclose(format_server_pid_file);
	}
	if (!no_fork) {
	    if (!quiet) {
		printf("Forking server to background\n");
	    }
	    if (fork() != 0) {
		/* I'm the parent, return now */
		exit(0);
	    }
	} else {
	    if (!quiet) {
		printf("Running...\n");
	    }
	}
    }
#ifndef HAVE_WINDOWS_H
    alarm(0);
#endif
}
