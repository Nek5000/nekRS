#include "config.h"
#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#include <sys/stat.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_SYSINFO
#include <sys/sysinfo.h>
#endif

#ifdef _MSC_VER
#define open _open
#define read _read
#define close _close
#include <sys/timeb.h>
#include <time.h>
#endif

#ifdef HAVE_MAC_SYSCTL
#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/socketvar.h>
#include <netinet/ip.h>
#include <netinet/ip_var.h>
#include <netinet/tcp.h>
#include <netinet/tcp_var.h>
#include <net/if.h>
#endif

#ifdef HAVE_UNAME
#include <sys/utsname.h>
#endif

#include "evpath.h"
#include "cm_internal.h"
#ifdef HAVE_COD_H
#include "cod.h"
#else
#define cod_assoc_externs(x,y)
#define cod_parse_for_context(x,y)
#endif
#include <fcntl.h>

/*#define CPU_FREQ_SCALING_MAX_FREQ "/sys/devices/system/cpu/cpu0/cpufreq/scaling_max_freq"*/
/*#define CPU_FREQ_SCALING_MIN_FREQ "/sys/devices/system/cpu/cpu0/cpufreq/scaling_min_freq"*/
/*#define CPU_FREQ_SCALING_CUR_FREQ "/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq"*/
/*#define CPU_FREQ_SCALING_AVAILABLE_FREQ "/sys/devices/system/cpu/cpu0/cpufreq/scaling_available_frequencies"*/

/*#define CPU_FREQ_SCALING_GOVERNOR "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"*/
/*#define CPU_FREQ_SCALING_AVAILABLE_GOVERNORS "/sys/devices/system/cpu/cpu0/cpufreq/scaling_available_governors"*/

/*#define CPU_FREQ_AFFECTED_CPU "/sys/devices/system/cpu/cpu0/cpufreq/affected_cpus"*/

/*#define CPU_FREQ_ONDEMAND_SAMPLING_RATE_MIN "/sys/devices/system/cpu/cpu0/cpufreq/ondemand/sampling_rate_min"*/
/*#define CPU_FREQ_ONDEMAND_SAMPLING_RATE_MAX "/sys/devices/system/cpu/cpu0/cpufreq/ondemand/sampling_rate_max"*/
/*#define CPU_FREQ_ONDEMAND_SAMPLING_RATE "/sys/devices/system/cpu/cpu0/cpufreq/ondemand/sampling_rate"*/
/*#define CPU_FREQ_ONDEMAND_UP_THRESHOLD "/sys/devices/system/cpu/cpu0/cpufreq/ondemand/up_threshold"*/
/*#define CPU_FREQ_ONDEMAND_IGNORE_NICE_LOAD "/sys/devices/system/cpu/cpu0/cpufreq/ondemand/ignore_nice_load"*/
/*#define CPU_FREQ_ONDEMAND_POWERSAVE_BIAD "/sys/devices/system/cpu/cpu0/cpufreq/ondemand/powersave_bias"*/

typedef struct sensor_slurp {
	const char *name;
	char buffer[8192];
} sensor_slurp, *sensor_slurp_ptr;

/*sensor_slurp proc_stat    = { "/proc/stat" };*/
/*sensor_slurp proc_loadavg = { "/proc/loadavg" };*/
/*sensor_slurp proc_meminfo = { "/proc/meminfo" };*/
/*sensor_slurp proc_net_dev = { "/proc/net/dev" };*/

char * skip_whitespace (const char *p)
{
    while (isspace((unsigned char)*p)) p++;
    return (char *)p;
}
 
char * skip_token (const char *p)
{
    while (isspace((unsigned char)*p)) p++;
    while (*p && !isspace((unsigned char)*p)) p++;
    return (char *)p;
}         

int slurpfile(const char *filename, char *buffer, int buflen)
{
#ifndef _MSC_VER
    int fd, read_len;
    fd = open(filename, O_RDONLY);
    if(fd < 0) {
	printf("open() error on file %s \n", filename); 
	exit(0);
    }
    read_len = read(fd, buffer, buflen);
    if(read_len <= 0) {
	printf("read() error on file %s \n", filename); 
	exit(0);
    }
    if (read_len == buflen) {
	--read_len;
	printf("slurpfile() read() buffer overflow on file %s", filename);
    }
    buffer[read_len] = '\0';
    close(fd);
    return read_len;
#else
    return 0;
#endif
}

char *update_file(struct sensor_slurp *sf) {
    (void) slurpfile(sf->name, &sf->buffer[0], 8192);
    return sf->buffer;
}

unsigned int num_cpustates_func ( void ) {
   char *p;
   unsigned int i=0;

   sensor_slurp proc_stat    = { "/proc/stat" };
   p = update_file(&proc_stat);

/**
** Skip initial "cpu" token
**/
   p = skip_token(p);
   p = skip_whitespace(p);
/*
** Loop over file until next "cpu" token is found.
** i=4 : Linux 2.4.x
** i=7 : Linux 2.6.x
*/
   while (strncmp(p,"cpu",3)) {
     p = skip_token(p);
     p = skip_whitespace(p);
     i++;
     }

   return i;
}

unsigned long total_jiffies_func ( void ) {

   char *p;
   unsigned long user_jiffies, nice_jiffies, system_jiffies, idle_jiffies,
                 wio_jiffies, irq_jiffies, sirq_jiffies;
	
   sensor_slurp proc_stat = { "/proc/stat" };
   p = update_file(&proc_stat);
   p = skip_token(p);
   p = skip_whitespace(p);
   user_jiffies = (unsigned int)strtod( p, &p );
   p = skip_whitespace(p);
   nice_jiffies = (unsigned int)strtod( p, &p ); 
   p = skip_whitespace(p);
   system_jiffies = (unsigned int)strtod( p , &p ); 
   p = skip_whitespace(p);
   idle_jiffies = (unsigned int)strtod( p , &p );

   int num_cpustates = num_cpustates_func();
   if(num_cpustates == 4) /*NUM_CPUSTATES_24X*/
       return user_jiffies + nice_jiffies + system_jiffies + idle_jiffies;

   p = skip_whitespace(p);
   wio_jiffies = (unsigned int)strtod( p , &p );
   p = skip_whitespace(p);
   irq_jiffies = (unsigned int)strtod( p , &p );
   p = skip_whitespace(p);
   sirq_jiffies = (unsigned int)strtod( p , &p );
  
   return user_jiffies + nice_jiffies + system_jiffies + idle_jiffies +
          wio_jiffies + irq_jiffies + sirq_jiffies; 
}   

void cpu_and_core_usage_func (double usage[])
{
#ifndef _MSC_VER
   char *p;
   int n, numcores, i;
   numcores = num_cpustates_func();
   n = numcores + 1;
   unsigned long user_jiffies[n], last_user_jiffies[n],  nice_jiffies[n], last_nice_jiffies[n], system_jiffies[n], last_system_jiffies[n], idle_jiffies[n], last_idle_jiffies[n],  wio_jiffies[n], last_wio_jiffies[n], irq_jiffies[n], last_irq_jiffies[n], sirq_jiffies[n], last_sirq_jiffies[n], delta_busy[n], delta_total[n];
      sensor_slurp proc_stat = { "/proc/stat" };
   p = update_file(&proc_stat);
   p = skip_token(p);
   p = skip_whitespace(p);
   for(i=0;i<n;++i)
   {
	user_jiffies[i] = (unsigned int)strtod( p, &p );
	p = skip_whitespace(p);
	nice_jiffies[i] = (unsigned int)strtod( p, &p );
	p = skip_whitespace(p);
	system_jiffies[i] = (unsigned int)strtod( p , &p );
	p = skip_whitespace(p);
	idle_jiffies[i] = (unsigned int)strtod( p , &p );
	p = skip_whitespace(p);
	wio_jiffies[i] = (unsigned int)strtod( p , &p );
	p = skip_whitespace(p);
	irq_jiffies[i] = (unsigned int)strtod( p , &p );
	p = skip_whitespace(p);
	sirq_jiffies[i] = (unsigned int)strtod( p , &p );
	p = skip_whitespace(p);
	delta_busy[i] = (user_jiffies[i] - last_user_jiffies[i]) + (nice_jiffies[i] - last_nice_jiffies[i]) + (system_jiffies[i] - last_system_jiffies[i]);
	delta_total[i] = delta_busy[i] + (idle_jiffies[i] - last_idle_jiffies[i]) + (wio_jiffies[i] - last_wio_jiffies[i]) + (irq_jiffies[i] - last_irq_jiffies[i]) + (sirq_jiffies[i] - last_sirq_jiffies[i]);
	usage[i] = (double) delta_busy[i] / (double) delta_total[i] ;
	last_user_jiffies[i] = user_jiffies[i];
	last_nice_jiffies[i] = nice_jiffies[i];
	last_system_jiffies[i] = system_jiffies[i];
	last_idle_jiffies[i] = idle_jiffies[i];
	last_wio_jiffies[i] = wio_jiffies[i];
	last_irq_jiffies[i] = irq_jiffies[i];
	last_sirq_jiffies[i] = sirq_jiffies[i];
	while (strncmp(p,"cpu",3)) {
	p = skip_token(p);
	p = skip_whitespace(p);
	}
	p = skip_token(p);
        p = skip_whitespace(p);

   }
#endif
}

double cpu_user_func ( void )
{
    char *p;
    double val;
    static double last_user_jiffies, last_total_jiffies;
    double user_jiffies, total_jiffies, diff;
   
      sensor_slurp proc_stat    = { "/proc/stat" };
    p = update_file(&proc_stat);

    p = skip_token(p);
    user_jiffies  = strtod( p , (char **)NULL );
    total_jiffies = (double) total_jiffies_func();

    diff = user_jiffies - last_user_jiffies; 
    
    if ( diff )
	val = (diff/(total_jiffies - last_total_jiffies))*100;
    else
	val = 0.0;

    last_user_jiffies  = user_jiffies;
    last_total_jiffies = total_jiffies;
    
    return val;
}

double cpu_nice_func ( void )
{
    char *p;
    double val;
    static double last_nice_jiffies, last_total_jiffies;
    double nice_jiffies, total_jiffies, diff;
    
      sensor_slurp proc_stat    = { "/proc/stat" };
    p = update_file(&proc_stat);
 
    p = skip_token(p);
    p = skip_token(p);
    nice_jiffies  = strtod( p , (char **)NULL );
    total_jiffies = (double) total_jiffies_func();

    diff = (nice_jiffies  - last_nice_jiffies);

    if ( diff )
	val = (diff/(total_jiffies - last_total_jiffies))*100;
    else
	val = 0.0;

    last_nice_jiffies  = nice_jiffies;
    last_total_jiffies = total_jiffies;
    
    return val;
}

double cpu_system_func ( void )
{
    char *p;
    static double val;
    static double last_system_jiffies,  system_jiffies,
	last_total_jiffies, total_jiffies, diff;
 
      sensor_slurp proc_stat    = { "/proc/stat" };
    p = update_file(&proc_stat);
    p = skip_token(p);
    p = skip_token(p);
    p = skip_token(p);
    system_jiffies = strtod( p , (char **)NULL );
	int num_cpustates = num_cpustates_func();
    if (num_cpustates > 4 ) /*NUM_CPUSTATES_24X)*/ {
	p = skip_token(p);
	p = skip_token(p);
	p = skip_token(p);
	system_jiffies += strtod( p , (char **)NULL ); /* "intr" counted in system */
	p = skip_token(p);
	system_jiffies += strtod( p , (char **)NULL ); /* "sintr" counted in system */
    }
    total_jiffies  = (double) total_jiffies_func();

    diff = system_jiffies  - last_system_jiffies;

    if ( diff )
	val = (diff/(total_jiffies - last_total_jiffies))*100;
    else
	val = 0.0;

    last_system_jiffies  = system_jiffies;
    last_total_jiffies = total_jiffies;   
	
    return val;
}
 
double cpu_idle_func ( void )
{
    char *p;
    static double val; 
    static double last_idle_jiffies,  idle_jiffies,
	last_total_jiffies, total_jiffies, diff;
									   
      sensor_slurp proc_stat    = { "/proc/stat" };
    p = update_file(&proc_stat);
		      
    p = skip_token(p);
    p = skip_token(p);
    p = skip_token(p);
    p = skip_token(p);
    idle_jiffies  = strtod( p , (char **)NULL );
    total_jiffies = (double) total_jiffies_func();
		 
    diff = idle_jiffies - last_idle_jiffies;

    if ( diff ) 
	val = (diff/(total_jiffies - last_total_jiffies))*100;
    else 
	val = 0.0; 
    last_idle_jiffies = idle_jiffies;
    last_total_jiffies = total_jiffies;
    return val; 
}

/**************NET FUNCTIONS*****************/
long received_bytes(char *iface)
{
	static long r_bytes=0;
#ifdef HAVE_MAC_SYSCTL
	char *buf = NULL, *lim, *next;
	static int mib[6];
  	size_t mlen,vlen;
	struct if_msghdr *ifm;
	unsigned int ifindex;
  	mib[0] = CTL_NET;
  	mib[1] = PF_ROUTE;
	mib[2] = 0;
	mib[3] = 0;
	mib[4] = NET_RT_IFLIST;
	mib[5] = 0;
	
	mlen = 6;
	ifindex = if_nametoindex(iface);
	sysctl(mib,mlen,NULL,&vlen,NULL,0);
	buf = malloc(vlen);
	sysctl(mib,mlen,buf,&vlen,NULL,0);
	lim = buf+vlen;
	for(next = buf; next < lim; )
	{	
		ifm = (struct if_msghdr *)next;
		next += ifm->ifm_msglen;
	//	printf("Index %d, ifm index %d\n",ifindex,ifm->ifm_index);
	//	printf("ibytes %ld\n",ifm->ifm_data.ifi_ibytes);
		if(ifindex == ifm->ifm_index)
		{
			//printf("Correct interface\n");
			r_bytes = (long) ifm->ifm_data.ifi_ibytes;
			break;
		}
	}
	//printf("Received bytes %ld\n",r_bytes);
	free(buf);
	
#else
	char *temp_if;
	char *p;
	if(iface == NULL)
	   temp_if = strdup("eth0");
	else
	   temp_if = strdup((const char *)iface);
	sensor_slurp proc_net = { "/proc/net/dev" };
	p = update_file(&proc_net);
	while (p && strncmp(p,temp_if,strlen(temp_if)))
	{
		p = skip_token(p);
		p = skip_whitespace(p);
	}

	p = skip_token(p);
	p = skip_whitespace(p);
	r_bytes = strtol( p ,NULL,10);
	free(temp_if);
#endif // end HAVE_MAC_SYSCTL
	return r_bytes;
}

long sent_bytes(char *iface)
{
	long s_bytes;
#ifdef HAVE_MAC_SYSCTL

	char *buf = NULL, *lim, *next;
	static int mib[6];
  	size_t mlen,vlen;
	struct if_msghdr *ifm;
	unsigned int ifindex;
  	mib[0] = CTL_NET;
  	mib[1] = PF_ROUTE;
	mib[2] = 0;
	mib[3] = 0;
	mib[4] = NET_RT_IFLIST;
	mib[5] = 0;
	
	mlen = 6;
	ifindex = if_nametoindex(iface);
	sysctl(mib,mlen,NULL,&vlen,NULL,0);
	buf = malloc(vlen);
	sysctl(mib,mlen,buf,&vlen,NULL,0);
	lim = buf+vlen;
	for(next = buf; next < lim; )
	{	
		ifm = (struct if_msghdr *)next;
		next += ifm->ifm_msglen;
		if(ifindex == ifm->ifm_index)
		{
		//	printf("Correct interface\n");
			s_bytes = (long) ifm->ifm_data.ifi_obytes;
			break;
		}
	}
//	printf("Sent bytes %ld\n",s_bytes);
	free(buf);
#else

        char *temp_if;
        char *p;
	int  i=0;
        if(iface == NULL)
           temp_if = strdup("eth0");
        else
           temp_if = strdup((const char *)iface);
        sensor_slurp proc_net = { "/proc/net/dev" };
        p = update_file(&proc_net);
	while (p && strncmp(p,temp_if,strlen(temp_if)))
        {
                p = skip_token(p);
                p = skip_whitespace(p);
        }
	while(i <= 8) {
	 	p = skip_token(p);
		p = skip_whitespace(p);
		++i;
	}
	s_bytes = strtol (p, NULL,10);
	free(temp_if);
#endif // end HAVE_MAC_SYSCTL
	return s_bytes;
	
	
}
/* 
 * This function returns the bandwidth of the network(in bits per second) as sampled over a certain window of time.
 * stage can be "start" where the network sampling begins, and at stage end the monitored network 
 * bandwidth is returned. On OS X the two calls have to be atleast 425 microseconds apart.
 * */

double net_bw(char *iface, char *stage)
{
	static long old_s_bytes, new_s_bytes, old_r_bytes, new_r_bytes, start_time, end_time;
	long s_bytes, r_bytes;
	struct timeval t;
	double net_bw;
#ifdef HAVE_MAC_SYSCTL
	char *buf = NULL, *lim, *next;
	static int mib[6];
  	size_t mlen,vlen;
	struct if_msghdr *ifm;
	unsigned int ifindex;
  	mib[0] = CTL_NET;
  	mib[1] = PF_ROUTE;
	mib[2] = 0;
	mib[3] = 0;
	mib[4] = NET_RT_IFLIST;
	mib[5] = 0;
	
	mlen = 6;
	ifindex = if_nametoindex(iface);
	sysctl(mib,mlen,NULL,&vlen,NULL,0);
	buf = malloc(vlen);
	sysctl(mib,mlen,buf,&vlen,NULL,0);
	lim = buf+vlen;
	for(next = buf; next < lim; )
	{	
		ifm = (struct if_msghdr *)next;
		next += ifm->ifm_msglen;
	//	printf("Index %d, ifm index %d\n",ifindex,ifm->ifm_index);
	//	printf("ibytes %ld\n",ifm->ifm_data.ifi_ibytes);
		if(ifindex == ifm->ifm_index)
		{
	//		printf("Correct interface\n");
			r_bytes = (long) ifm->ifm_data.ifi_ibytes;
			s_bytes = (long) ifm->ifm_data.ifi_obytes;
			break;
		}
	}
//	printf("Received bytes %ld Sent bytes %ld\n",r_bytes,s_bytes);
	free(buf);
	
#else

	int i=0;
	char *temp_if;
	char *p;

	if(iface == NULL)
		temp_if = strdup("eth0");
	else
		temp_if = strdup((const char*)iface);
	sensor_slurp proc_net = { "/proc/net/dev" };
        p = update_file(&proc_net);
#ifndef HAVE_WINDOWS_H
	gettimeofday(&t, NULL);
#else
	/* GSE...  No gettimeofday on windows.
	 * Must use _ftime, get millisec time, convert to usec.  Bleh.
	 */

	struct _timeb nowb;
	_ftime(&nowb);
	t.tv_sec = (long)nowb.time;
	t.tv_usec = nowb.millitm * 1000;
#endif
        while (p && strncmp(p,temp_if,strlen(temp_if)))
        {
                p = skip_token(p);
                p = skip_whitespace(p);
        }
	
        p = skip_token(p);
        p = skip_whitespace(p);
        r_bytes = strtol( p ,NULL,10);
	while(i <= 7) {
                p = skip_token(p);
                p = skip_whitespace(p);
                ++i;
        }
	s_bytes = strtol (p, NULL,10);
	free(temp_if);
#endif // end HAVE_MAC_SYSCTL
	if(strncmp(stage,"start",3)==0)
	{
		printf("Start %ld %ld \n",s_bytes,r_bytes);
		
		old_s_bytes = s_bytes;
		old_r_bytes = r_bytes;
		start_time = t.tv_sec * 1000000 + t.tv_usec;

		return 0.0;
	}
	printf("End %ld %ld \n",s_bytes,r_bytes);
	new_s_bytes = s_bytes;
	new_r_bytes = r_bytes;
	end_time = t.tv_sec * 1000000 + t.tv_usec;
	net_bw = (double) (((new_s_bytes - old_s_bytes) + (new_r_bytes - old_r_bytes)) / (double)(end_time - start_time))*1000000*8;
	printf("Bandwidth = %f\n",net_bw);
	return net_bw;			
	
}

/**************TIMING FUNCTIONS**************/

double dgettimeofday( void )
{
#ifdef HAVE_GETTIMEOFDAY
    double timestamp;
    struct timeval now;
    gettimeofday(&now, NULL);
    timestamp = now.tv_sec + now.tv_usec* 1.0e-6 ;
    return timestamp;
#else
    return -1;
#endif
}

/**************OS FUNCTIONS**************/
char*  os_type() {
#ifndef _MSC_VER
  static struct utsname output;
  static int first = 1;
  if (first) {
    uname(&output);
    first = 0;
  }
  return strdup(output.sysname);
#else
    return strdup("Windows");
#endif
}

char*  os_release() {
#ifndef _MSC_VER
  static struct utsname output;
  static int first = 1;
  if (first) {
    uname(&output);
    first = 0;
  }
  return strdup(output.release);
#else
    return strdup("Windows Release");
#endif
}

/* Should probably test if gethostname & uname exist on box before using them.... */
char* hostname() {
  static char val[256];
  static int first = 1;
  if (first) {
    gethostname(&val[0],256);
    first = 0;
  }
  return strdup(&val[0]);
}
    

/**************Stat FUNCTIONS**************/
double  stat_uptime() {
  double val=0;
#ifdef HAVE_SYSINFO
  struct sysinfo info;
  sysinfo(&info);
  val = (double) info.uptime;
#else
#ifdef HAVE_MAC_SYSCTL
  static int mib[4];
  struct timeval boottime;
  size_t mlen,vlen;
  mlen = 2;
  mib[0] = CTL_KERN;
  mib[1] = KERN_BOOTTIME;
  vlen = sizeof(struct timeval);
  sysctl(mib,mlen,&boottime,&vlen,NULL,0);  
  val =  dgettimeofday() - (double) boottime.tv_sec - boottime.tv_usec* 1.0e-6;

#endif // end HAVE_MAC_SYSCTL
#endif // end HAVE_SYSINFO

  return val;
 
}

/* double stat_cpu_user() { */
/* } */

/* double stat_cpu_system() { */
/* } */

double stat_loadavg_one ( void )
{
  double val;
#ifdef HAVE_GETLOADAVG
  double loadavg[3];
  getloadavg(loadavg,3);
  val = loadavg[0];
  
#else      
  sensor_slurp proc_loadavg = { "/proc/loadavg" };
  val = strtod( update_file(&proc_loadavg), (char **)NULL);
#endif //end HAVE_GETLOADAVG
    return val;
}

double stat_loadavg_five ( void )
{
  double val;
#ifdef HAVE_GETLOADAVG
  double loadavg[3];
  getloadavg(loadavg,3);
  val = loadavg[1];
  
#else      
    char *p;
 
      sensor_slurp proc_loadavg = { "/proc/loadavg" };
    p = update_file(&proc_loadavg);
    p = skip_token(p);
    val = strtod( p, (char **)NULL);
#endif //end HAVE_GETLOADAVG
    return val;
}

double stat_loadavg_fifteen ( void )
{
  double val;  
#ifdef HAVE_GETLOADAVG
  double loadavg[3];
  getloadavg(loadavg,3);
  val = loadavg[2];
  
#else    
  char *p;
  sensor_slurp proc_loadavg = { "/proc/loadavg" };
  p = update_file(&proc_loadavg);
  
  p = skip_token(p);
  p = skip_token(p);
  val = strtod( p, (char **)NULL);
#endif //end HAVE_GETLOADAVG
  return val;
}

/**************MEMORY FUNCTIONS**************/

unsigned long vm_mem_total() {
  unsigned long val = 0;
#ifdef HAVE_MAC_SYSCTL
  /* static int flag=0; */
  static int mib[4];
  size_t mlen,vlen;
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  mlen = 2;
  /* if (flag == 0) { */
  /*   sysctlnametomib("hw.memsize",mib,&mlen); */
  /*   flag=1; */
  /* } */
  vlen = sizeof(unsigned long);
  sysctl(mib,mlen,&val,&vlen,NULL,0);
  
#else
  char *p;
  
  sensor_slurp proc_meminfo = { "/proc/meminfo" };
  p = strstr( update_file(&proc_meminfo), "MemTotal:" );
  if(p) {
    p = skip_token(p);
    val = strtoul( p, NULL, 10 );
  } else {
    val = 0;
  }
#endif //end of HAVE_MAC_SYSCTL  
  return val;
}


unsigned long vm_mem_free() {
  unsigned long val = 0;
#ifdef HAVE_MAC_SYSCTL
  /* static int flag=0; */
  static int mib[2] = {CTL_VM, VM_SWAPUSAGE};
  struct xsw_usage swapval;
  size_t vlen;
  /* mlen = 4; */
  /* if (flag == 0) { */
  /*   sysctlnametomib("vm.swapusage",mib,&mlen); */
  /*   flag=1; */
  /* } */
  vlen = sizeof(struct xsw_usage);
  sysctl(mib,2,&swapval,&vlen,NULL,0);
  val = (unsigned long) swapval.xsu_avail;
  
#else
    char *p;

    sensor_slurp proc_meminfo = { "/proc/meminfo" };
    p = strstr( update_file(&proc_meminfo), "MemFree:" );
    if(p) {
	p = skip_token(p);
	val = strtoul( p, NULL, 10 );
    } else {
	val = 0;
    }
#endif //end of HAVE_MAC_SYSCTL
    
    return val;
}


unsigned long vm_swap_total ( void )
{
  unsigned long val = 0;
#ifdef HAVE_MAC_SYSCTL
  static int mib[2] = {CTL_VM, VM_SWAPUSAGE};
  struct xsw_usage swapval;
  size_t vlen;
  vlen = sizeof(struct xsw_usage);
  sysctl(mib,2,&swapval,&vlen,NULL,0);
  val = (unsigned long) swapval.xsu_total;
  
#else
  char *p;
   
  sensor_slurp proc_meminfo = { "/proc/meminfo" };
  p = strstr( update_file(&proc_meminfo), "SwapFree:" );
  if(p) {
    p = skip_token(p);
    val =  strtoul( p, NULL, 10 ); 
  } else {
    val = 0;
  }
#endif //end HAVE_MAC_SYSCTL
  return val;
}

unsigned long vm_swap_free ( void )
{
  unsigned long val = 0;
#ifdef HAVE_MAC_SYSCTL
  static int mib[2] = {CTL_VM, VM_SWAPUSAGE};
  struct xsw_usage swapval;
  size_t vlen;
  vlen = sizeof(struct xsw_usage);
  sysctl(mib,2,&swapval,&vlen,NULL,0);
  val = (unsigned long) swapval.xsu_avail;
  
#else
    char *p;
   
      sensor_slurp proc_meminfo = { "/proc/meminfo" };
    p = strstr( update_file(&proc_meminfo), "SwapFree:" );
    if(p) {
	p = skip_token(p);
	val = strtoul( p, NULL, 10 ); 
    } else {
	val = 0;
    }
#endif //end HAVE_MAC_SYSCTL
    return val;
}

/**************HARDWARE FUNCTIONS**************/
int hw_cpus () {
  long val = -1;
#ifdef HAVE_MAC_SYSCTL
  static int mib[2] = {CTL_HW, HW_NCPU};
  size_t vlen;
  vlen = sizeof(long);
  sysctl(mib,2,&val,&vlen,NULL,0);
#else
#if defined(HAVE_SYSCONF)
  val = (long) sysconf(_SC_NPROCESSORS_ONLN);
#else
  val = 1; // if we don't know any better...
#endif // end HAVE_SYSCONFIG
#endif //end of HAVE_MAC_SYSCTL
  return (int) val;
  
}


long hw_cpu_max_freq( void )
{
  long val = 0;
#ifdef HAVE_MAC_SYSCTL
  size_t vlen;
  vlen = sizeof(long);
  sysctlbyname("hw.cpufrequency_max",&val,&vlen,NULL,0);
  val = val/1e6;
#else
  char *p;
  
  struct stat struct_stat;
  char sys_devices_system_cpu[32];
  const char *CPU_FREQ_SCALING_MAX_FREQ = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_max_freq";
  if ( stat(CPU_FREQ_SCALING_MAX_FREQ, &struct_stat) == 0 ) {
    if(slurpfile(CPU_FREQ_SCALING_MAX_FREQ, sys_devices_system_cpu, 32)) {
      p = sys_devices_system_cpu;
      val = (strtol( p, (char **)NULL , 10 ) / 1000 );
    }
  }
#endif // end HAVE_MAC_SYSCTL
  return val;
}

long hw_cpu_min_freq( void )
{
  long val = -1;
#ifdef HAVE_MAC_SYSCTL
  size_t vlen;
  vlen = sizeof(long);
  sysctlbyname("hw.cpufrequency_min",&val,&vlen,NULL,0);
  val = val/1e6;

#else
  char *p;

  struct stat struct_stat;
  char sys_devices_system_cpu[32];
  const char *CPU_FREQ_SCALING_MIN_FREQ = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_min_freq";
  if ( stat(CPU_FREQ_SCALING_MIN_FREQ, &struct_stat) == 0 ) {
    if(slurpfile(CPU_FREQ_SCALING_MIN_FREQ, sys_devices_system_cpu, 32)) {
      p = sys_devices_system_cpu;
      val = (strtol( p, (char **)NULL , 10 ) / 1000 );
    }
  }
#endif //end HAVE_MAC_SYSCTL
  return val;
}

long hw_cpu_curr_freq( void )
{
  long val = 0;
#ifdef HAVE_MAC_SYSCTL
  /* static int flag=0; */
  static int mib[2] = {CTL_HW, HW_CPU_FREQ};
  size_t vlen;
  /* if (flag == 0) { */
  /*   sysctlnametomib("hw.cpufrequency",mib,&mlen); */
  /*   printf("\nmib = { %i, %i, %i, %i};\n",mib[0],mib[1],mib[2],mib[3]); */
  /*   flag=1; */
  /* } */
  vlen = sizeof(long);
  sysctl(mib,2,&val,&vlen,NULL,0);
  val = val/1e6;
  
#else
  char *p;
  
  struct stat struct_stat;
  char sys_devices_system_cpu[32];
  const char *CPU_FREQ_SCALING_CUR_FREQ = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq";
  if ( stat(CPU_FREQ_SCALING_CUR_FREQ, &struct_stat) == 0 ) {
    if(slurpfile(CPU_FREQ_SCALING_CUR_FREQ, sys_devices_system_cpu, 32)) {
      p = sys_devices_system_cpu;
      val = (unsigned int) (strtoul( p, (char **)NULL , 10 ) / 1000 );
    }
  }
#endif //end HAVE_MAC_SYSCTL
  return val;
}


void
add_metrics_routines(stone_type stone, cod_parse_context context)
{
    static char extern_string[] = "\
       double    dgettimeofday();         \n		\
/* number of cpus */  \n				\
       int           hw_cpus();           \n		\
/* minimum allowed frequency in MHz */ \n	\
       long           hw_cpu_min_freq();    \n \
/* maximum allowed frequency in MHz  */ \n	\
       long           hw_cpu_max_freq();   \n	\
/* current frequency  in Mhz */ \n	\
       long           hw_cpu_curr_freq();  \n	\
/* a string to identify the local OS type -- ie Linux  */ \n 	\
       char*       os_type();              \n	\
/* a string to identify the current release -- ie FC14 */  \n	\
       char*       os_release();          \n\
/* a string to identify the hostname -- ie maquis1 */  \n	\
       char*       hostname();          \n\
/* time in seconds that the computer has been up  */  \n	\
       double    stat_uptime();           \n\
/* load average over the last one minute  */ \n	\
       double    stat_loadavg_one();       \n	\
/* load average over the last five minutes  */ \n	\
       double    stat_loadavg_five();      \n	\
/* load average over the last fifteen minute  */ \n	\
       double    stat_loadavg_fifteen();   \n	\
/* total physical memory  */  \n	\
       unsigned long    vm_mem_total();   \n	\
/* free physical memory   */ \n	\
       unsigned long    vm_mem_free();    \n	\
/* total swap available */  \n	\
       unsigned long    vm_swap_total();   \n	\
/* free swap  */  \n	\
       unsigned long    vm_swap_free();     \n";

    static cod_extern_entry externs[] = {
	{"dgettimeofday", (void *) 0},	        // 0
	{"hw_cpus", (void *) 0},		// 1
	{"hw_cpu_min_freq", (void *) 0},	// 2
	{"hw_cpu_max_freq", (void *) 0},	// 3
	{"hw_cpu_curr_freq", (void *) 0},	// 4
	{"os_type", (void *) 0},		// 5
	{"os_release", (void *) 0},		// 6
	{"hostname", (void *) 0},               // 7
	{"stat_uptime", (void *) 0},		// 8
	{"stat_loadavg_one", (void *) 0},	// 9
	{"stat_loadavg_five", (void *) 0},	// 10
	{"stat_loadavg_fifteen", (void *) 0},	// 11
	{"vm_mem_total", (void *) 0},		// 12
	{"vm_mem_free", (void *) 0},		// 13
	{"vm_swap_total", (void *) 0},		// 14
	{"vm_swap_free", (void*) 0},		// 15
	{(void *) 0, (void *) 0}
    };

    (void) stone;
    /*
     * some compilers think it isn't a static initialization to put this
     * in the structure above, so do it explicitly.
     */
    externs[0].extern_value = (void *) (intptr_t) dgettimeofday;
    externs[1].extern_value = (void *) (intptr_t) hw_cpus;
    externs[2].extern_value = (void *) (intptr_t) hw_cpu_min_freq;
    externs[3].extern_value = (void *) (intptr_t) hw_cpu_max_freq;
    externs[4].extern_value = (void *) (intptr_t) hw_cpu_curr_freq;
    externs[5].extern_value = (void *) (intptr_t) os_type;
    externs[6].extern_value = (void *) (intptr_t) os_release;
    externs[7].extern_value = (void *) (intptr_t) hostname;
    externs[8].extern_value = (void *) (intptr_t) stat_uptime;
    externs[9].extern_value = (void *) (intptr_t) stat_loadavg_one;
    externs[10].extern_value = (void *) (intptr_t) stat_loadavg_five;
    externs[11].extern_value = (void *) (intptr_t) stat_loadavg_fifteen;
    externs[12].extern_value = (void *) (intptr_t) vm_mem_total;
    externs[13].extern_value = (void *) (intptr_t) vm_mem_free;
    externs[14].extern_value = (void *) (intptr_t) vm_swap_total;
    externs[15].extern_value = (void *) (intptr_t) vm_swap_free;

    cod_assoc_externs(context, externs);
    cod_parse_for_context(extern_string, context);
}


/*int main(int argc, char **argv) {
	//num_cpustates_func();
	printf (" CPU USER %  :  %lf \n", cpu_user_func());
	printf (" CPU NICE %  :  %lf \n", cpu_nice_func());
	printf (" CPU SYSTEM %  :  %lf \n", cpu_system_func());
	printf (" CPU IDLE %  :  %lf \n", cpu_idle_func());
	
	printf (" LOAD 1 %  :  %lf \n", load_one_func());
	printf (" LOAD 5 %  :  %lf \n", load_five_func());
	printf (" LOAD 15 %  :  %lf \n", load_fifteen_func());

	printf (" MEM FREE %  :  %lf \n", mem_free_func());
	printf (" MEM BUFF %  :  %lf \n", mem_buffers_func());
	printf (" MEM CACHE %  :  %lf \n", mem_cached_func());

	printf (" SWAP FREE %  :  %lf \n", swap_free_func());
	printf (" GET TIMESTAMP %  :  %lf \n", gettimeofday_func());

	printf (" CPU_MAX_FREQ %  :  %d \n", cpu_max_freq_func());
	printf (" CPU_MIN_FREQ %  :  %d \n", cpu_min_freq_func());
	printf (" CPU_CUR_FREQ %  :  %d \n", cpu_cur_freq_func());
	int *val;
	val = cpu_available_freq_func();
	printf (" CPU_AVAILABLE_FREQ %  :  %d , %d , %d \n", val[0], val[1], val[2]);

	printf (" CPU_FREQ_GOVERNOR %  :  %s \n", cpu_scaling_governor_func());
	char **avail;
	avail = cpu_scaling_available_governors_func();
	printf (" CPU_FREQ_AVAILABLE_GOVERNORS %  :  %s , %s , %s \n", avail[0], avail[1], avail[2]);
}*/

/* Dead functions.

char *cpu_scaling_governor_func( void )
{
    char *p = NULL;
	struct stat struct_stat;
	char sys_devices_system_cpu[32];
	const char *CPU_FREQ_SCALING_GOVERNOR = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor";
    if ( stat(CPU_FREQ_SCALING_GOVERNOR, &struct_stat) == 0 ) {
		if(slurpfile(CPU_FREQ_SCALING_GOVERNOR, sys_devices_system_cpu, 32)) {
			p = strdup(sys_devices_system_cpu);
		}
    }
    return p;
}

char **cpu_scaling_available_governors_func( void )
{
    static char *val[3];
    char *p;

	struct stat struct_stat;
	char sys_devices_system_cpu_available[128];
	const char *CPU_FREQ_SCALING_AVAILABLE_GOVERNORS  = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_available_governors";
    if ( stat(CPU_FREQ_SCALING_AVAILABLE_GOVERNORS, &struct_stat) == 0 ) {
		if(slurpfile(CPU_FREQ_SCALING_AVAILABLE_GOVERNORS, sys_devices_system_cpu_available, 128)) {
			p = strdup(sys_devices_system_cpu_available);
			val[0] = strtok(p, " ");
			val[1] = strtok(NULL, " ");
			val[2] = strtok(NULL, " ");
		}
    }
    return val;
}

int *cpu_available_freq_func( void )
{
    char *p;
    static int val[3];

	struct stat struct_stat;
	char sys_devices_system_cpu_available[128];
	const char *CPU_FREQ_SCALING_AVAILABLE_FREQ = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_available_frequencies";
    if ( stat(CPU_FREQ_SCALING_AVAILABLE_FREQ, &struct_stat) == 0 ) {
		if(slurpfile(CPU_FREQ_SCALING_AVAILABLE_FREQ, sys_devices_system_cpu_available, 128)) {
			p = sys_devices_system_cpu_available;
			val[0] = (strtol( p, (char **)NULL , 10 ) / 1000 );
			p = skip_token(p);
			val[1] = (strtol( p, (char **)NULL , 10 ) / 1000 );
			p = skip_token(p);
			val[2] = (strtol( p, (char **)NULL , 10 ) / 1000 );
		}
    }
    return val;
}

double mem_cached_func ( void )
{
    char *p;
    double val;

      sensor_slurp proc_meminfo = { "/proc/meminfo" };
    p = strstr( update_file(&proc_meminfo), "Cached:");
    if(p) {
	p = skip_token(p);
	val = atof( p );
    } else {
	val = 0;
    }
   return val;
}

double mem_buffers_func ( void )
{
    char *p;
    double val;

      sensor_slurp proc_meminfo = { "/proc/meminfo" };
    p = strstr( update_file(&proc_meminfo), "Buffers:" );
    if(p) {
	p = skip_token(p);
	val = atof( p ); 
    } else {
	val = 0;
    }

    return val;
}

*/

#ifdef NVML_FOUND
#include "nvml.h"

void gpu_mon_init_func()
{
   nvmlInit();
}

void gpu_mon_cleanup_func()
{
   nvmlShutdown();
}

unsigned int num_gpus_func()
{
   unsigned int ret;
   nvmlDeviceGetCount(&ret);
   return ret;
}

void gpu_utilization_func(unsigned int *gpu_util_rates)
{
   unsigned int num_dev, i;
   unsigned int *util;
   nvmlUtilization_t nvmlUtil;
   nvmlDevice_t dev;
   num_dev = num_gpus_func();
 //  gpu_util_rates=(unsigned int *)malloc(sizeof(unsigned int)*num_dev);
   for(i=0;i<num_dev;++i)
   {
	nvmlDeviceGetHandleByIndex(i,&dev);
	nvmlDeviceGetUtilizationRates(dev,&nvmlUtil);
	*(gpu_util_rates+i)=nvmlUtil.gpu;
   }
      
}

void gpu_free_mem_func(unsigned long long *gpu_free_mem)
{
   unsigned int num_dev,i;
   nvmlDevice_t dev;
   nvmlMemory_t mem;
   num_dev=num_gpus_func();
  // gpu_free_mem=(unsigned long long *)malloc(sizeof(unsigned long long)*num_dev);
   for(i=0;i<num_dev;++i)
   {
	nvmlDeviceGetHandleByIndex(i,&dev);
	nvmlDeviceGetMemoryInfo(dev,&mem);
	*(gpu_free_mem+i)=mem.free;
   }
}

void gpu_total_mem_func(unsigned long long *gpu_total_mem)
{
   unsigned int num_dev,i;
   nvmlDevice_t dev;
   nvmlMemory_t mem;
   num_dev=num_gpus_func();
 //  gpu_total_mem=(unsigned long long *)malloc(sizeof(unsigned long long)*num_dev);
   for(i=0;i<num_dev;++i)
   {
        nvmlDeviceGetHandleByIndex(i,&dev);
        nvmlDeviceGetMemoryInfo(dev,&mem);
        *(gpu_total_mem+i)=mem.total;
   }
}

void gpu_used_mem_func(unsigned long long *gpu_used_mem)
{
   unsigned int num_dev,i;
   nvmlDevice_t dev;
   nvmlMemory_t mem;
   num_dev=num_gpus_func();
 //  gpu_used_mem=(unsigned long long *)malloc(sizeof(unsigned long long)*num_dev);
   for(i=0;i<num_dev;++i)
   {
        nvmlDeviceGetHandleByIndex(i,&dev);
        nvmlDeviceGetMemoryInfo(dev,&mem);
        *(gpu_used_mem+i)=mem.used;
   }
}

void gpu_temp_func(unsigned int *temp)
{
   unsigned int num_dev,i,t;
   nvmlDevice_t dev;
   num_dev=num_gpus_func();
//   temp=(unsigned int *)malloc(sizeof(unsigned int)*num_dev);
   for(i=0;i<num_dev;++i)
   {
	nvmlDeviceGetHandleByIndex(i,&dev);
	nvmlDeviceGetTemperature(dev,NVML_TEMPERATURE_GPU,&t);
	*(temp+i)=t;
   }
}

void gpu_power_util_func(unsigned int *power)
{
   unsigned int num_dev,i,p;
   nvmlDevice_t dev;
   num_dev=num_gpus_func();
//   power=(unsigned int*)malloc(sizeof(unsigned int)*num_dev);
   for(i=0;i<num_dev;++i)
   {
	nvmlDeviceGetHandleByIndex(i,&dev);
	nvmlDeviceGetPowerUsage(dev,&p);
	*(power+i)=p;
   }

}

#endif
