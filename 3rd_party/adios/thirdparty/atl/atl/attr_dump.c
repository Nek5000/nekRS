#include "config.h"
#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdint.h>
#include "atl.h"

int
main(int argc, char **argv)
{
    attr_list list;

    if (argc != 2) {
	printf("Usage:  %s <stringified attributes>\n", argv[0]);
	exit(1);
    }
	    
    list = attr_list_from_string(argv[1]);
    dump_attr_list(list);
    printf("\n");

    return 0;
}
