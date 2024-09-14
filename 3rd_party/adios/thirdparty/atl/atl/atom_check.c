#include <stdlib.h>
#include <stdio.h>
#include <unix_defs.h>
#include <stdint.h>
#include <atl.h>
#include <ctype.h>
#include <string.h>

int
main()
{
    atom_t atom = attr_atom_from_string("ECho_attr_test_atom");
    if (atom != 0x732a7ade) {
	printf("Atom test for ECho_Attr_test_atom is %x, expected 732a7ade\n", atom);
	return 1;
    }
    return 0;
}
