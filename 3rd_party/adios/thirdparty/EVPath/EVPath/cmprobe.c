#include "config.h"

#include <stdio.h>
#include <atl.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>
#include "evpath.h"
#include <errno.h>
#ifdef _MSC_VER
#define sleep(x) Sleep(x*1000)
#endif

int
main(int argc, char **argv)
{
    CManager cm;
    CMConnection conn = NULL;
    static int atom_init = 0;

    cm = CManager_create();
    (void) CMfork_comm_thread(cm);

    atom_t CM_REBWM_RLEN=-1, CM_REBWM_REPT=-1, CM_BW_MEASURE_INTERVAL=-1, CM_BW_MEASURE_SIZE=-1, CM_BW_MEASURE_SIZEINC=-1, CM_TRANSPORT=-1;

    if (atom_init == 0) {
	CM_REBWM_RLEN = attr_atom_from_string("CM_REBWM_RLEN");
	CM_REBWM_REPT = attr_atom_from_string("CM_REBWM_REPT");
	CM_BW_MEASURE_INTERVAL = attr_atom_from_string("CM_BW_MEASURE_INTERVAL");
	CM_BW_MEASURE_SIZE = attr_atom_from_string("CM_BW_MEASURE_SIZE");
	CM_BW_MEASURE_SIZEINC = attr_atom_from_string("CM_BW_MEASURE_SIZEINC");
	CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
	atom_init++;
    }


    if (argc == 1) {
	attr_list contact_list, listen_list = NULL;
	char *transport = NULL;
	if ((transport = getenv("CMTransport")) != NULL) {

	    listen_list = create_attr_list();
	    add_string_attr(listen_list, CM_TRANSPORT, strdup(transport));
	}
	CMlisten_specific(cm, listen_list);
	contact_list = CMget_contact_list(cm);
	printf("Contact list \"%s\"\n", attr_list_to_string(contact_list));
	CMsleep(cm, 1200);
    } else {
	int size = 100;
	int i;
	int N, repeat_time, size_inc;

	attr_list contact_list = NULL;
	attr_list bw_list;

	for (i = 1; i < argc; i++) {
	    char *final;
	    long value;
	    errno = 0;
	    value = strtol(argv[i], &final, 10);
	    if ((errno == 0) && (final == (argv[i] + strlen(argv[i])))) {
		/* valid number as an argument, must be byte size */
		size = (int) value;
	    } else {
		contact_list = attr_list_from_string(argv[i]);
		if (contact_list == NULL) {
		    printf("Argument \"%s\" not recognized as size or contact list\n",
			   argv[i]);
		}
	    }
	}
	if (contact_list == NULL) {
	    exit(1);
	}
	conn = CMinitiate_conn(cm, contact_list);
	if (conn == NULL) {
	    printf("No connection\n");
	    exit(1);
	}


	N=3;
	repeat_time=2;
	size_inc=200;
    
	bw_list = create_attr_list();	    
	    
	/*Each measurement done by CMregressive_probe_bandwidth uses N streams of different size, each stream is sent out for repeat_time times. */
	/*For scheduled measurment, Each measurement is done every 2 seconds*/
	add_int_attr(bw_list, CM_REBWM_RLEN,  N);
	add_int_attr(bw_list, CM_REBWM_REPT,  repeat_time);
	add_int_attr(bw_list, CM_BW_MEASURE_INTERVAL, 2);
	add_int_attr(bw_list, CM_BW_MEASURE_SIZE, size); 
	add_int_attr(bw_list, CM_BW_MEASURE_SIZEINC, size_inc);
		
/* Example when invoking CMregressive_probe_bandwidth on demand is needed: */
	for(i=1; i<120; i++)
	{
	    double bandwidth;

	    bandwidth=CMprobe_bandwidth(conn, size, bw_list);
	    printf("Estimated bandwidth at size %d is %f Mbps\n", size, bandwidth);
	    sleep(1);
	    if(bandwidth>0)
		size+=size_inc;
	}
	    
    }
    CManager_close(cm);
    return 0;
}
