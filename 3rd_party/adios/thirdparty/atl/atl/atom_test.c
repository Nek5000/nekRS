#include <stdlib.h>
#include <stdio.h>
#include <unix_defs.h>
#include <stdint.h>
#include <atl.h>
#include <ctype.h>
#include <string.h>

extern
void
set_string_and_atom(atom_server as, char *str, atom_t atom);

void
print_atom(int atom)
{
    unsigned char atomc[4];
    int i;
    memcpy(&atomc[0], &atom, sizeof(atomc));
    for (i=0; i < sizeof(atomc); i++) {
	if (!isprint((int)atomc[i])) atomc[i] = ' ';
    }
    printf("%d, 0x%x, '%c%c%c%c'", atom, atom, atomc[0], atomc[1], 
	   atomc[2], atomc[3]);
}

int
main()
{
    atom_server as;
    as = init_atom_server(prefill_atom_cache);
    printf("Connected to server \"%s\"\n", get_server_id(as));
    printf("Usage:   D<number> to translate a decimal value to a string\n");
    printf("         X<number> to translate a hex value to a string\n");
    printf("         S<string> to translate a string to a value\n");
    printf("         A<decimal number> <string> to set a string/value association\n");
    printf("         a<hex number> <string> to set a string/value association\n");
    while (1) {
	char line[80];
	long atom = 0;
	int len;
	int radix = 10;

	if (fgets(line, sizeof(line), stdin) == NULL)
	    return 0;
	len = (int) strlen(line);
	if (len > 0)
	    line[len - 1] = 0;	/* truncate \n */
	if (line[0] == 0)
	    continue;
	switch(line[0]) {
	case 'X':
	    radix = 16;
	    /* falling through */
	case 'D': {
	    char *str, *end;
	    atom = strtol(&line[1], &end, radix);
	    if (end == &line[1]) {
		printf("\tMissing number\n");
		continue;
	    }
	    while (isspace((int)*end)) end++;
	    if (*end != 0) printf("\tskipped %s at end of number\n", end);
	    printf("inquiring for atom \n");
	    print_atom(atom); printf("\n");
	    str = string_from_atom(as, atom);
	    if (str) {
		printf("	string translation is \"%s\"\n", str);
		free(str);
	    } else {
		printf("	no atom %ld\n", atom);
	    }
	    break;
	}
	case 'S': {
	    char *str = &line[1];
	    int atom;
	    while (isspace((int)*str)) str++;
	    printf("inquiring for string \"%s\"\n", str);
	    atom = atom_from_string(as, str);
	    printf("\tatom translation is ");print_atom(atom);
	    printf("\n");
	    break;
	}
	case 'a':
	    radix = 16;
	    /* falling through */
	case 'A': {
	    long atom;
	    char *end;
	    atom = strtol(&line[1], &end, radix);
	    if (end == &line[1]) {
		printf("\tMissing number\n");
		continue;
	    }
	    while (isspace((int)*end)) end++;
	    if (*end == 0) {
		printf("\tMissing string\n");
		continue;
	    }
	    set_string_and_atom(as, end, atom);
	    break;
	case 'e':
	    exit(0);
	    break;
	}
	default:
	    printf("Command not understood \"%s\"\n", line);
	}
    }
}
