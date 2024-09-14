#ifndef PACKAGE
#include "config.h"
#endif

#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#ifdef HAVE_SYS_SOCKIO_H
#include <sys/sockio.h>
#endif
#ifdef HAVE_SYS_SOCKET_H
#include <sys/socket.h>
#endif
#ifdef HAVE_WINDOWS_H
#include <winsock.h>
#define __ANSI_CPP__
#else
#include <netinet/in.h>
#include <arpa/inet.h>
#include <net/if.h>
#endif
/*
 *  self_ip_addr.c
 *
 *  This file lives in the CM project.  It is used by many transports 
 *  and servers to abstract the cumbersomeness in determining a 
 *  reasonable IP contact address.  
 *
 *  If you found this file in ECho, or in another project that has a 
 *  need for determining IP or fully qualified hostnames, 
 *                 -- DON'T CHANGE THE LOCAL COPY --
 *  Instead, change it in CM and reimport the changed version.
 */

static int
get_self_ip_addr()
{
    struct hostent *host;
    char buf[256];
    char **p;
#ifdef SIOCGIFCONF
    char *ifreqs;
    struct ifreq *ifr;
    struct sockaddr_in *sai;
    struct ifconf ifaces;
    int ifrn;
    int ss;
#endif
    int rv = 0;
    char *IP_string = getenv("CERCS_IP");
    if (IP_string != NULL) {
	in_addr_t ip = inet_addr(IP_string);
	if (ip != (in_addr_t) -1) return ntohl(ip);
    }
    gethostname(buf, sizeof(buf));
    host = gethostbyname(buf);
    if (host != NULL) {
	for (p = host->h_addr_list; *p != 0; p++) {
	    struct in_addr *in = *(struct in_addr **) p;
	    if (ntohl(in->s_addr) != INADDR_LOOPBACK) {
		return (ntohl(in->s_addr));
	    }
	}
    }
    /*
     *  Since we couldn't find an IP address in some logical way, we'll open
     *  a DGRAM socket and ask it first for the list of its interfaces, and
     *  then checking for an interface we can use, and then finally asking that
     *  interface what its address is.
     */
#ifdef SIOCGIFCONF
    ss = socket(AF_INET, SOCK_DGRAM, 0);
    ifaces.ifc_len = 64 * sizeof(struct ifreq);
    ifaces.ifc_buf = ifreqs = malloc(ifaces.ifc_len);
    /*
     *  if we can't SIOCGIFCONF we're kind of dead anyway, bail.
     */
    if (ioctl(ss, SIOCGIFCONF, &ifaces) >= 0) {
	ifr = ifaces.ifc_req;
	ifrn = ifaces.ifc_len / sizeof(struct ifreq);
	for (; ifrn--; ifr++) {
	    /*
	     * Basically we'll take the first interface satisfying 
	     * the following: 
	     *   up, running, not loopback, address family is INET.
	     */
	    ioctl(ss, SIOCGIFFLAGS, ifr);
	    sai = (struct sockaddr_in *) &(ifr->ifr_addr);
	    if (ifr->ifr_flags & IFF_LOOPBACK) {
		continue;
	    }
	    if (!(ifr->ifr_flags & IFF_UP)) {
		continue;
	    }
	    if (!(ifr->ifr_flags & IFF_RUNNING)) {
		continue;
	    }
	    /*
	     * sure would be nice to test for AF_INET here but it doesn't
	     * cooperate and I've done enough for now ...
	     * if (sai->sin_addr.s.addr != AF_INET) continue;
	    */
	    if (sai->sin_addr.s_addr == INADDR_ANY)
		continue;
	    if (sai->sin_addr.s_addr == INADDR_LOOPBACK)
		continue;
	    rv = ntohl(sai->sin_addr.s_addr);
	    break;
	}
    }
    close(ss);
    free(ifreqs);
#endif
    /*
     *  Absolute last resort.  If we can't figure out anything else, look
     *  for the CM_LAST_RESORT_IP_ADDR environment variable.
     */
    if (rv == 0) {
	char *c = getenv("CM_LAST_RESORT_IP_ADDR");
	if (c != NULL) {
	    rv = inet_addr(c);
	}
    }
    /*
     *	hopefully by now we've set rv to something useful.  If not,
     *  GET A BETTER NETWORK CONFIGURATION.
     */
    return rv;
}

