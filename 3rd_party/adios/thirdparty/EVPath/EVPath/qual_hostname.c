#include "config.h"

#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#ifdef HAVE_SYS_SOCKIO_H
#include <sys/sockio.h>
#endif
#ifdef HAVE_WINDOWS_H
#define FD_SETSIZE 1024
#include <winsock2.h>
#define __ANSI_CPP__
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <net/if.h>
#include <sys/ioctl.h>
#include <ifaddrs.h>
#endif

#if defined (__INTEL_COMPILER)
/*  Allow unordered operations */
#  pragma warning (disable: 981)
//  allow int conversions
#  pragma warning (disable: 2259)
#endif

static int ipv4_is_loopback(int addr)
{
  return (htonl(addr) & htonl(0xff000000)) == htonl(0x7f000000);
}

/*
 *  qual_hostname.c
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
get_self_ip_addr(void *cm, CMtrans_services svc)
{
    struct hostent *host;
    char hostname_buf[256];
    char **p;
#ifdef HAVE_GETIFADDRS
  struct ifaddrs *if_addrs = NULL;
  struct ifaddrs *if_addr = NULL;
  void *tmp = NULL;
  char buf[INET6_ADDRSTRLEN];
#endif
#ifdef SIOCGIFCONF
    char *ifreqs;
    struct ifreq *ifr;
    struct sockaddr_in *sai;
    struct ifconf ifaces;
    int ifrn;
    int ss;
    int ipv4_count = 0;
    int ipv6_count = 0;
#endif
    int rv = 0;
#ifdef HAVE_GETIFADDRS
    if (getifaddrs(&if_addrs) == 0) {    
        char *interface;
	// Print possible addresses
	for (if_addr = if_addrs; if_addr != NULL; if_addr = if_addr->ifa_next) {
	    int family;
	    if (!if_addr->ifa_addr) continue;
	    family = if_addr->ifa_addr->sa_family;
	    if ((family != AF_INET) && (family != AF_INET6)) continue;
	    if (if_addr->ifa_addr->sa_family == AF_INET) {
	        tmp = &((struct sockaddr_in *)if_addr->ifa_addr)->sin_addr;
		ipv4_count++;
	    } else {
	        tmp = &((struct sockaddr_in6 *)if_addr->ifa_addr)->sin6_addr;
		ipv6_count++;
	    }
	    if (svc) {
	        svc->trace_out(cm, "CM<transport> IP possibility -> %s : %s",
			       if_addr->ifa_name,
			       inet_ntop(family, tmp, buf, sizeof(buf)));
	    }
	}
	if ((interface = getenv("CM_INTERFACE")) != NULL) {
	    for (if_addr = if_addrs; if_addr != NULL; if_addr = if_addr->ifa_next) {
	        int family;
	        if (!if_addr->ifa_addr) continue;
		family = if_addr->ifa_addr->sa_family;
		if (family != AF_INET) continue;  /* currently not looking for ipv6 */
		if (strcmp(if_addr->ifa_name, interface) != 0) continue;
		tmp = &((struct sockaddr_in *)if_addr->ifa_addr)->sin_addr;
		if (svc) {
		    svc->trace_out(cm, "CM<transport> Interface specified, returning ->%s : %s",
				   if_addr->ifa_name,
				   inet_ntop(family, tmp, buf, sizeof(buf)));
		}
		freeifaddrs(if_addrs);
		return (ntohl(*(uint32_t*)tmp));
	    }
	    printf("Warning!  CM_INTERFACE specified as \"%s\", but no active interface by that name found\n", interface);
	}
	    
	gethostname(hostname_buf, sizeof(hostname_buf));
	host = gethostbyname(hostname_buf);
	if (host != NULL) {
	    for (p = host->h_addr_list; *p != 0; p++) {
		struct in_addr *in = *(struct in_addr **) p;
		if (!ipv4_is_loopback(ntohl(in->s_addr))) {
		    if (svc)
			svc->trace_out(cm, "CM<transport> Prefer IP associated with hostname net -> %d.%d.%d.%d",
				       *((unsigned char *) &in->s_addr),
				       *(((unsigned char *) &in->s_addr) + 1),
				       *(((unsigned char *) &in->s_addr) + 2),
				       *(((unsigned char *) &in->s_addr) + 3));
		    freeifaddrs(if_addrs);
		    return (ntohl(in->s_addr));
		}
	    }
	}
	/* choose the first thing that's not a loopback interface */
	for (if_addr = if_addrs; if_addr != NULL; if_addr = if_addr->ifa_next) {
	    int family;
	    uint32_t ret_ip;
	    if (!if_addr->ifa_addr) continue;
	    family = if_addr->ifa_addr->sa_family;
	    if (family != AF_INET) continue;  /* currently not looking for ipv6 */
	    if ((if_addr->ifa_flags & IFF_LOOPBACK) != 0)  continue;
	    tmp = &((struct sockaddr_in *)if_addr->ifa_addr)->sin_addr;
	    if (svc) {
		svc->trace_out(cm, "CM<transport> get_self_ip_addr returning first avail -> %s : %s",
			       if_addr->ifa_name,
			       inet_ntop(family, tmp, buf, sizeof(buf)));
	    }
	    ret_ip = (ntohl(*(uint32_t*)tmp));
	    freeifaddrs(if_addrs);
	    return ret_ip;
	}
    }
    if (if_addrs) freeifaddrs(if_addrs);
#endif	
    gethostname(hostname_buf, sizeof(hostname_buf));
    host = gethostbyname(hostname_buf);
    if (host != NULL) {
	for (p = host->h_addr_list; *p != 0; p++) {
	    struct in_addr *in = *(struct in_addr **) p;
	    if (!ipv4_is_loopback(ntohl(in->s_addr))) {
		if (svc)
		    svc->trace_out(cm, "CM<transport> - Get self IP addr %lx, net %d.%d.%d.%d",
				   ntohl(in->s_addr),
				   *((unsigned char *) &in->s_addr),
				   *(((unsigned char *) &in->s_addr) + 1),
				   *(((unsigned char *) &in->s_addr) + 2),
				   *(((unsigned char *) &in->s_addr) + 3));
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
		if (svc)
		    svc->trace_out(cm, "CM<transport> - Get self IP addr %lx, rejected, loopback",
				   ntohl(sai->sin_addr.s_addr));
		continue;
	    }
	    if (!(ifr->ifr_flags & IFF_UP)) {
		if (svc)
		    svc->trace_out(cm, "CM<transport> - Get self IP addr %lx, rejected, not UP",
				   ntohl(sai->sin_addr.s_addr));
		continue;
	    }
	    if (!(ifr->ifr_flags & IFF_RUNNING)) {
		if (svc)
		    svc->trace_out(cm, "CM<transport> - Get self IP addr %lx, rejected, not RUNNING",
				   ntohl(sai->sin_addr.s_addr));
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
	    if (svc)
		svc->trace_out(cm, "CM<transport> - Get self IP addr DHCP %lx, net %d.%d.%d.%d",
			       ntohl(sai->sin_addr.s_addr),
			       *((unsigned char *) &sai->sin_addr.s_addr),
			       *(((unsigned char *) &sai->sin_addr.s_addr) + 1),
			       *(((unsigned char *) &sai->sin_addr.s_addr) + 2),
			       *(((unsigned char *) &sai->sin_addr.s_addr) + 3));
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
	if (svc)
	    svc->trace_out(cm, "CM<transport> - Get self IP addr at last resort");
	if (c != NULL) {
	    if (svc)
		svc->trace_out(cm, "CM<transport> - Translating last resort %s", c);
	    rv = inet_addr(c);
	}
    }
    /*
     *	hopefully by now we've set rv to something useful.  If not,
     *  GET A BETTER NETWORK CONFIGURATION.
     */
    return rv;
}

static int
is_private_IP(int IP)
{
    if ((IP & 0xffff0000) == 0xC0A80000) return 1;	/* equal 192.168.x.x */
    if ((IP & 0xffff0000) == 0xB6100000) return 1;	/* equal 182.16.x.x */
    if ((IP & 0xff000000) == 0x0A000000) return 1;	/* equal 10.x.x.x */
    return 0;
}

static void
get_qual_hostname(void *cm, char *buf, int len, CMtrans_services svc, attr_list attrs,
		  int *network_p)
{
    struct hostent *host = NULL;

    char *network_string = getenv("CM_NETWORK");
    char *hostname_string = getenv("CERCS_HOSTNAME");
    if (hostname_string != NULL) {
	strncpy(buf, hostname_string, len);
	return;
    }
    (void)get_qual_hostname;
    gethostname(buf, len);
    buf[len - 1] = '\0';
    if (memchr(buf, '.', strlen(buf)) == NULL) {
	/* no dots, probably not fully qualified */
#ifdef HAVE_GETDOMAINNAME
	int end = strlen(buf);
	buf[end] = '.';
	if (getdomainname((&buf[end]) + 1, len - end - 1) == -1) {
	    buf[end+1]=0;
	}
	if (buf[end + 1] == 0) {
	    char *tmp_name;
	    struct hostent *host = gethostbyname(buf);
	    buf[end] = 0;
	    /* getdomainname was useless, hope that gethostbyname helps */
	    if (host) {
		tmp_name = (gethostbyname(buf))->h_name;
		strncpy(buf, tmp_name, len);
	    }		
	}
#else
	{
	    /* no getdomainname, hope that gethostbyname will help */
	    struct hostent *he = gethostbyname(buf);
	    char *tmp_name;
	    if (he) {
		tmp_name = (gethostbyname(buf))->h_name;
		strncpy(buf, tmp_name, len);
	    }
	}
#endif
	buf[len - 1] = '\0';
    }
    svc->trace_out(cm, "CM<transport> - Tentative Qualified hostname %s", buf);
    if (memchr(buf, '.', strlen(buf)) == NULL) {
	/* useless hostname if it's not fully qualified */
	buf[0] = 0;
    }
    if ((buf[0] != 0) && ((host = gethostbyname(buf)) == NULL)) {
	/* useless hostname if we can't translate it */
	buf[0] = 0;
    }
    if (host != NULL) {
	char **p;
	int good_addr = 0;
	for (p = host->h_addr_list; *p != 0; p++) {
	    struct in_addr *in = *(struct in_addr **) p;
	    if (!ipv4_is_loopback(ntohl(in->s_addr))) {
		good_addr++;
		svc->trace_out(cm,
			       "CM<transport> - Hostname gets good addr %lx, %d.%d.%d.%d",
			       ntohl(in->s_addr),
			       *((unsigned char *) &in->s_addr),
			       *(((unsigned char *) &in->s_addr) + 1),
			       *(((unsigned char *) &in->s_addr) + 2),
			       *(((unsigned char *) &in->s_addr) + 3));
	    }
	}
	if (good_addr == 0) {
	    /* 
	     * even a fully qualifiedhostname that doesn't get us a valid
	     * IP addr is useless
	     */
	    buf[0] = 0;
	}
    }
    if (buf[0] == 0) {
	/* bloody hell, what do you have to do? */
	struct in_addr IP;
	extern int h_errno;
	IP.s_addr = htonl(get_self_ip_addr(cm, svc));
	svc->trace_out(cm, "CM<transport> - No hostname yet, trying gethostbyaddr on IP %lx", IP);
	if (!is_private_IP(ntohl(IP.s_addr))) {
	    host = gethostbyaddr((char *) &IP, sizeof(IP), AF_INET);
	    if (host != NULL) {
		svc->trace_out(cm, "     result was %s", host->h_name);
		strncpy(buf, host->h_name, len);
	    } else {
		svc->trace_out(cm, "     FAILED, errno %d", h_errno);
	    }
	}
    }
    if (network_string == NULL) {
	static atom_t CM_NETWORK_POSTFIX = -1;
	if (CM_NETWORK_POSTFIX == -1) {
	    CM_NETWORK_POSTFIX = attr_atom_from_string("CM_NETWORK_POSTFIX");
	}
	if (!get_string_attr(attrs, CM_NETWORK_POSTFIX, &network_string)) {
	    svc->trace_out(cm, "TCP/IP transport found no NETWORK POSTFIX attribute");
	} else {
	    svc->trace_out(cm, "TCP/IP transport found NETWORK POSTFIX attribute %s", network_string);
	}
    }
    if (network_string != NULL) {
	size_t name_len = strlen(buf) + 2 + strlen(network_string);
	char *new_name_str = svc->malloc_func(name_len);
	char *first_dot = strchr(buf, '.');

	/* stick the CM_NETWORK value in there */
	memset(new_name_str, 0, name_len);
	*first_dot = 0;
	first_dot++;
	sprintf(new_name_str, "%s%s.%s", buf, network_string, first_dot);
	if (gethostbyname(new_name_str) != NULL) {
	    /* host has no NETWORK interface */
	    strcpy(buf, new_name_str);
	    if (network_p) (*network_p)++;
	}
	svc->free_func(new_name_str);
    }

    if ((buf[0] == 0) ||
	((host = gethostbyname(buf)) == NULL) ||
	(memchr(buf, '.', strlen(buf)) == NULL)) {
	/* just use the bloody IP address */
	struct in_addr IP;
	IP.s_addr = htonl(get_self_ip_addr(cm, svc));
	if (IP.s_addr != 0) {
	    struct in_addr ip;
	    char *tmp;
	    ip.s_addr = htonl(get_self_ip_addr(cm, svc));
	    tmp = inet_ntoa(ip);
	    strncpy(buf, tmp, len);
	} else {
	    static int warn_once = 0;
	    if (warn_once == 0) {
		warn_once++;
		svc->trace_out(cm, "Attempts to establish your fully qualified hostname, or indeed any\nuseful network name, have failed horribly.  using localhost.\n");
	    }
	    strncpy(buf, "localhost", len);
	}
    }
    svc->trace_out(cm, "CM<transport> - GetQualHostname returning %s", buf);
}

