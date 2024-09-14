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
#include <Ws2def.h>
#include <ws2tcpip.h>
#define __ANSI_CPP__
#ifndef INET_ADDRSTRLEN
#define INET_ADDRSTRLEN 50
#endif
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <net/if.h>
#include <sys/ioctl.h>
#include <ifaddrs.h>
#include <ctype.h>
#endif
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>
#include <stdio.h>
#ifdef STDC_HEADERS
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#include "evpath.h"
#include "cm_transport.h"

#if defined (__INTEL_COMPILER)
/*  Allow unordered operations */
#  pragma warning (disable: 981)
//  allow int conversions
#  pragma warning (disable: 2259)
#endif

static int atom_init = 0;
static int CM_IP_INTERFACE = 0;
static int CM_IP_PORT = 0;

static int ipv4_is_loopback(int addr)
{
  return (htonl(addr) & htonl(0xff000000)) == htonl(0x7f000000);
}

static void dump_output(int length_estimate, char *format, ...);

static int
get_self_ip_iface(CMTransport_trace trace_func, void* trace_data, char *iface)
{
    struct hostent *host = NULL;
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
    static int first_call = 1;
#endif
    int rv = 0;
#ifdef HAVE_GETIFADDRS
    if (getifaddrs(&if_addrs) == 0) {    
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
	    trace_func(trace_data, "CM<IP_CONFIG> IP possibility -> %s : %s",
		       if_addr->ifa_name,
		       inet_ntop(family, tmp, buf, sizeof(buf)));
	    if ((family == AF_INET) && first_call) {
		dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "IP_CONFIG Possible interface %s : IPV4 %s\n",
			    if_addr->ifa_name,
			    inet_ntop(family, tmp, buf, sizeof(buf)));
	    } else {
		// until we support IPV6, don't dump info
		//dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "IP_CONFIG Possible interface %s : IPV6 %s\n",
		//	    if_addr->ifa_name,
		//	    inet_ntop(family, tmp, buf, sizeof(buf)));
	    }
	}
	if (!iface) iface = getenv(IPCONFIG_ENVVAR_PREFIX "INTERFACE");
	if (iface != NULL) {
	    trace_func(trace_data, "CM<IP_CONFIG> searching for interface %s\n", iface);
	    if (first_call) dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "IP_CONFIG interface %s requested\n", iface);
	    for (if_addr = if_addrs; if_addr != NULL; if_addr = if_addr->ifa_next) {
	        int family;
		uint32_t IP;
	        if (!if_addr->ifa_addr) continue;
		family = if_addr->ifa_addr->sa_family;
		if (family != AF_INET) continue;  /* currently not looking for ipv6 */
		if (strncmp(if_addr->ifa_name, iface, strlen(iface)) != 0) continue;
		tmp = &((struct sockaddr_in *)if_addr->ifa_addr)->sin_addr;
		trace_func(trace_data, "CM<IP_CONFIG> Interface specified, returning ->%s : %s",
			   if_addr->ifa_name,
			   inet_ntop(family, tmp, buf, sizeof(buf)));
		if (first_call)
		dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "IP_CONFIG interface %s found, using IP %s\n", iface,
			   inet_ntop(family, tmp, buf, sizeof(buf)));
		IP = ntohl(*(uint32_t*)tmp);
		free(if_addrs);
		first_call = 0;
		return IP;
	    }
	    printf("Warning!  " IPCONFIG_ENVVAR_PREFIX "INTERFACE specified as \"%s\", but no active interface by that name found\n", iface);
	}
	    
	first_call = 0;
	gethostname(hostname_buf, sizeof(hostname_buf));
	if (index(hostname_buf, '.') != NULL) {
	    /* don't even check for host if not fully qualified */
	    host = gethostbyname(hostname_buf);
	}
	if (host != NULL) {
	    for (p = host->h_addr_list; *p != 0; p++) {
		struct in_addr *in = *(struct in_addr **) p;
		if (!ipv4_is_loopback(ntohl(in->s_addr))) {
		    char str[INET_ADDRSTRLEN];
					      
		    inet_ntop(AF_INET, &(in->s_addr), str, sizeof(str));
		    trace_func(trace_data, "CM<IP_CONFIG> Prefer IP associated with hostname net -> %s", str);
		    free(if_addrs);
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
	    trace_func(trace_data, "CM<IP_CONFIG> get_self_ip_addr returning first avail -> %s : %s",
			       if_addr->ifa_name,
			       inet_ntop(family, tmp, buf, sizeof(buf)));
	    ret_ip = (ntohl(*(uint32_t*)tmp));
	    free(if_addrs);
	    return ret_ip;
	}
    }
#endif	
    gethostname(hostname_buf, sizeof(hostname_buf));
    if (strchr(hostname_buf, '.') != NULL) {
	/* don't even check for host if not fully qualified */
	host = gethostbyname(hostname_buf);
    }
    if (host != NULL) {
	for (p = host->h_addr_list; *p != 0; p++) {
	    struct in_addr *in = *(struct in_addr **) p;
	    if (!ipv4_is_loopback(ntohl(in->s_addr))) {
	        char str[INET_ADDRSTRLEN];
					      
		inet_ntop(AF_INET, &(in->s_addr), str, sizeof(str));
		trace_func(trace_data, "CM<IP_CONFIG> - Get self IP addr %lx, net %s",
			   ntohl(in->s_addr), str);
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
	    char str[INET_ADDRSTRLEN];
	    /*
	     * Basically we'll take the first interface satisfying 
	     * the following: 
	     *   up, running, not loopback, address family is INET.
	     */
	    ioctl(ss, SIOCGIFFLAGS, ifr);
	    sai = (struct sockaddr_in *) &(ifr->ifr_addr);
	    if (ifr->ifr_flags & IFF_LOOPBACK) {
		trace_func(trace_data, "CM<IP_CONFIG> - Get self IP addr %s, rejected, loopback",
			   inet_ntoa(sai->sin_addr));
		continue;
	    }
	    if (!(ifr->ifr_flags & IFF_UP)) {
		trace_func(trace_data, "CM<IP_CONFIG> - Get self IP addr %lx, rejected, not UP",
			   ntohl(sai->sin_addr.s_addr));
		continue;
	    }
	    if (!(ifr->ifr_flags & IFF_RUNNING)) {
		trace_func(trace_data, "CM<IP_CONFIG> - Get self IP addr %lx, rejected, not RUNNING",
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
					      
	    inet_ntop(AF_INET, &(sai->sin_addr.s_addr), str, sizeof(str));
	    trace_func(trace_data, "CM<IP_CONFIG> - Get self IP addr DHCP %lx, net %s",
		       ntohl(sai->sin_addr.s_addr), str);
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
	char *c = getenv(IPCONFIG_ENVVAR_PREFIX "LAST_RESORT_IP_ADDR");
	trace_func(trace_data, "CM<IP_CONFIG> - Get self IP addr at last resort");
	if (c != NULL) {
	    trace_func(trace_data, "CM<IP_CONFIG> - Translating last resort %s", c);
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
get_self_ip_addr(CMTransport_trace trace_func, void* trace_data)
{
    return get_self_ip_iface(trace_func, trace_data, NULL);
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
get_qual_hostname(char *buf, int len, attr_list attrs,
		  int *network_p, CMTransport_trace trace_func, void *trace_data)
{
    struct hostent *host = NULL;

    char *network_string = getenv(IPCONFIG_ENVVAR_PREFIX "NETWORK");
    char *hostname_string = getenv(IPCONFIG_ENVVAR_PREFIX "HOSTNAME");
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
		tmp_name = host->h_name;
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
    trace_func(trace_data, "CM<IP_CONFIG> - Tentative Qualified hostname %s", buf);
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
	        char str[INET_ADDRSTRLEN];
		good_addr++;
					      
		inet_ntop(AF_INET, &(in->s_addr), str, sizeof(str));
		trace_func(trace_data, "CM<IP_CONFIG> - Hostname gets good addr %lx, %s",
			   ntohl(in->s_addr), str);
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
	char *iface;
	if (get_string_attr(attrs, CM_IP_INTERFACE, &iface)){
	    IP.s_addr = htonl(get_self_ip_iface(trace_func, trace_data, iface));
	} else {
	    IP.s_addr = htonl(get_self_ip_addr(trace_func, trace_data));
	}
	trace_func(trace_data, "CM<IP_CONFIG> - No hostname yet, trying gethostbyaddr on IP %lx", IP);
	if (!is_private_IP(ntohl(IP.s_addr))) {
	    host = gethostbyaddr((char *) &IP, sizeof(IP), AF_INET);
	    if (host != NULL) {
		trace_func(trace_data, "     result was %s", host->h_name);
		strncpy(buf, host->h_name, len);
	    } else {
		trace_func(trace_data, "     FAILED, errno %d", h_errno);
	    }
	}
    }
    if (network_string == NULL) {
	static atom_t CM_NETWORK_POSTFIX = -1;
	if (CM_NETWORK_POSTFIX == -1) {
	    CM_NETWORK_POSTFIX = attr_atom_from_string(IPCONFIG_ENVVAR_PREFIX "NETWORK_POSTFIX");
	}
	if (!get_string_attr(attrs, CM_NETWORK_POSTFIX, &network_string)) {
	    trace_func(trace_data, "TCP/IP transport found no NETWORK POSTFIX attribute");
	} else {
	    trace_func(trace_data, "TCP/IP transport found NETWORK POSTFIX attribute %s", network_string);
	}
    }
    if (network_string != NULL) {
	size_t name_len = strlen(buf) + 2 + strlen(network_string);
	char *new_name_str = malloc(name_len);
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
	free(new_name_str);
    }

    if ((buf[0] == 0) ||
	((host = gethostbyname(buf)) == NULL) ||
	(memchr(buf, '.', strlen(buf)) == NULL)) {
	/* just use the bloody IP address */
	struct in_addr IP;
	IP.s_addr = htonl(get_self_ip_addr(trace_func, trace_data));
	if (IP.s_addr != 0) {
	    struct in_addr ip;
	    ip.s_addr = htonl(get_self_ip_addr(trace_func, trace_data));
	    if (!inet_ntop(AF_INET, &(ip.s_addr), buf, len)) {
	        trace_func(trace_data, "inet_ntop failed\n");
	    }
	} else {
	    static int warn_once = 0;
	    if (warn_once == 0) {
		warn_once++;
		trace_func(trace_data, "Attempts to establish your fully qualified hostname, or indeed any\nuseful network name, have failed horribly.  using localhost.\n");
	    }
	    strncpy(buf, "localhost", len);
	}
    }
    trace_func(trace_data, "CM<IP_CONFIG> - GetQualHostname returning %s", buf);
}


static char *IP_config_diagnostics = NULL;
static int IP_config_output_len = -1;

extern char*
IP_get_diagnostics(CManager cm, CMTransport_trace trace_out)
{
    char *output;
    IP_config_output_len = 0;   /* setup output */
    get_IP_config(NULL, 0, NULL, NULL, NULL, NULL, NULL,
		  trace_out, cm);
    IP_config_output_len = -1;   /* disable output */
    output = IP_config_diagnostics;
    IP_config_diagnostics = NULL;  /* not ours anymore */
    return output;
}

static void
dump_output(int length_estimate, char *format, ...)
{
    char buf[1024];
    char *tmp = &buf[0];
    va_list ap;
    int free_tmp = 0;

    if (IP_config_output_len == -1) return;

    IP_config_diagnostics = realloc(IP_config_diagnostics, IP_config_output_len + length_estimate + 1);
    tmp = IP_config_diagnostics + IP_config_output_len;

    if (length_estimate > 1024) {
	tmp = malloc(length_estimate + 1);
	free_tmp++;
    }
#ifdef STDC_HEADERS
    va_start(ap, format);
#else
    va_start(ap);
#endif
    vsprintf(tmp, format, ap);
    va_end(ap);
    IP_config_output_len += (int)strlen(tmp);
    if (free_tmp) free(tmp);
}

#ifndef HOST_NAME_MAX
#define HOST_NAME_MAX 255
#endif
#ifdef _MSC_VER
static int inet_aton(const char* cp, struct in_addr* addr)
{
    addr->s_addr = inet_addr(cp);
    return (addr->s_addr == INADDR_NONE) ? 0 : 1;
}
#endif


extern void
get_IP_config(char *hostname_buf, int len, int* IP_p, int *port_range_low_p, int *port_range_high_p, 
	      int *use_hostname_p, attr_list attrs, CMTransport_trace trace_func, void *trace_data)
{
    static int first_call = 1;
    static char determined_hostname[HOST_NAME_MAX+1];
    static int determined_IP = -1;
    static int port_range_low = 0, port_range_high = 0;
    static int use_hostname = 0;
    char hostname_to_use[HOST_NAME_MAX+1];
    int IP_to_use;
    char *iface = NULL;

    if (first_call) {
	char *preferred_hostname = getenv(IPCONFIG_ENVVAR_PREFIX "HOSTNAME");
	char *preferred_IP = getenv(IPCONFIG_ENVVAR_PREFIX "IP");
	char *port_range = getenv(IPCONFIG_ENVVAR_PREFIX "PORT_RANGE");
	CM_IP_INTERFACE = attr_atom_from_string("IP_INTERFACE");
	CM_IP_PORT = attr_atom_from_string("IP_PORT");
	(void) CM_IP_PORT;
	atom_init++;
	first_call = 0;
	determined_hostname[0] = 0;
	
	if (preferred_IP != NULL) {
	    struct in_addr addr;
	    if (preferred_hostname) printf("Warning, " IPCONFIG_ENVVAR_PREFIX "HOSTNAME and " IPCONFIG_ENVVAR_PREFIX "IP are both set, preferring " IPCONFIG_ENVVAR_PREFIX "IP\n");
	    if (inet_aton(preferred_IP, &addr) == 0) {
		fprintf(stderr, "Invalid address %s specified for " IPCONFIG_ENVVAR_PREFIX "IP\n", preferred_IP);
	    } else {
		trace_func(trace_data, "CM IP_CONFIG Using IP specified in " IPCONFIG_ENVVAR_PREFIX "IP, %s", preferred_IP);
		determined_IP =  (ntohl(addr.s_addr));
		dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "IP environment variable found, preferring IP %s\n", preferred_IP);
	    }
	} else if (preferred_hostname != NULL) {
	    struct hostent *host;
	    use_hostname = 1;
	    trace_func(trace_data, "CM<IP_CONFIG> CM_HOSTNAME set to \"%s\", running with that.", preferred_hostname);
	    dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "HOSTNAME environment variable found, trying \"%s\"\n", preferred_hostname);
	    host = gethostbyname(preferred_hostname);
	    strcpy(determined_hostname, preferred_hostname);
	    if (!host) {
		printf("Warning, " IPCONFIG_ENVVAR_PREFIX "HOSTNAME is \"%s\", but gethostbyname fails for that string.\n", preferred_hostname);
		dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "HOSTNAME \"%s\" fails to translate to IP address.\n", preferred_hostname);
	    } else {
		char **p;
		for (p = host->h_addr_list; *p != 0; p++) {
		    struct in_addr *in = *(struct in_addr **) p;
		    if (!ipv4_is_loopback(ntohl(in->s_addr))) {
		        char str[INET_ADDRSTRLEN];
					      
			inet_ntop(AF_INET, &(in->s_addr), str, sizeof(str));
			trace_func(trace_data, "CM IP_CONFIG Prefer IP associated with hostname net -> %s", str);
			dump_output(1023, "\t" "HOSTNAME \"%s\" translates to IP address %s.\n", preferred_hostname, str);
			determined_IP  = (ntohl(in->s_addr));
		    }
		}
		if (determined_IP == -1) {
		    dump_output(1023, "\t" "No non-loopback interfaces found for hostname \"%s\", rejected for IP use.\n", preferred_hostname);
		}
	    }
	} else {
	    get_qual_hostname(determined_hostname, sizeof(determined_hostname) - 1, NULL /* attrs */, NULL, trace_func, trace_data);
	    dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "IP_CONFIG best guess hostname is \"%s\"\n", determined_hostname);
	}
	if (determined_IP == -1) {
	    /* I.E. the specified hostname didn't determine what IP we should use */
	    char str[INET_ADDRSTRLEN];
	    struct in_addr addr;
	    determined_IP = get_self_ip_addr(trace_func, trace_data);
					      
	    addr.s_addr = ntohl(determined_IP);
	    inet_ntop(AF_INET, &(addr.s_addr), str, sizeof(str));
	    dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "IP_CONFIG best guess IP is \"%s\"\n", str);
	}
	if (port_range == NULL) {
	    // no getenv
	    port_range = EVPATH_DEFAULT_PORT_RANGE;
	}
	if (port_range != NULL) {
	    if (isalpha(port_range[0])) {
		char *t = strdup(port_range);
		char *lower = t;
		for ( ; *lower; ++lower) *lower = tolower(*lower);
		if (strcmp(t, "any") == 0) {
		    port_range_high = -1;
		    port_range_low = -1;
		} else {
		    printf(IPCONFIG_ENVVAR_PREFIX "PORT_RANGE spec not understood \"%s\"\n", port_range);
		}
		free(t);
	    } else {
		if (sscanf(port_range, "%d:%d", &port_range_high, &port_range_low) != 2) {
		    printf(IPCONFIG_ENVVAR_PREFIX "PORT_RANGE spec not understood \"%s\"\n", port_range);
		} else {
		    if (port_range_high < port_range_low) {
			int tmp = port_range_high;
			port_range_high = port_range_low;
			port_range_low = tmp;
		    }
		}
	    }
	}
	if (port_range_low == -1) {
	    dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "IP_CONFIG specified port range is \"ANY\" (unspecified)\n");
	} else {
	    dump_output(1023, "\t" IPCONFIG_ENVVAR_PREFIX "IP_CONFIG specified port range is %d:%d\n", port_range_low, port_range_high);
	}	    
    }


    if (get_string_attr(attrs, CM_IP_INTERFACE, &iface)) {
	/* don't use predetermined stuff ! */
	get_qual_hostname(hostname_to_use, sizeof(hostname_to_use) - 1 , attrs, NULL, trace_func, trace_data);
	IP_to_use = get_self_ip_iface(trace_func, trace_data, iface);
    } else {
	strcpy(hostname_to_use, determined_hostname);
	IP_to_use = determined_IP;
    }
    if (hostname_buf && (len > strlen(determined_hostname))) {
	strcpy(hostname_buf, hostname_to_use);
    }
    if (IP_p && (determined_IP != -1)) {
	*IP_p = IP_to_use;
    }
    
    if (port_range_low_p) {
	*port_range_low_p = port_range_low;
    }
    if (port_range_high_p) {
	*port_range_high_p = port_range_high;
    }
    if (use_hostname_p) {
	*use_hostname_p = use_hostname;
    }
    {
	char buf[256];
	int net_byte_order = htonl(IP_to_use);
	trace_func(trace_data, "CM<IP_CONFIG> returning hostname \"%s\", IP %s, use_hostname = %d, port range %d:%d",
		   hostname_to_use, inet_ntop(AF_INET, &net_byte_order, &buf[0], 256), use_hostname, port_range_low, port_range_high);
    }
}
