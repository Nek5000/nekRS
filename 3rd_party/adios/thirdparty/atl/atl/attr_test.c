#include "config.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdio.h>
#include <stdint.h>

#include "atl.h"
#ifdef _MSC_VER
    #define strdup _strdup
    #include <io.h>
#pragma warning(disable: 4996)
#endif

int
main()
{
    attr_list list = create_attr_list();
    attr_list list2 = create_attr_list();
    char *string;
    char *string2;
    atom_t red_atom = -1;
    atom_t ip_atom = -1;
    atom_t blue_atom = -1;
    atom_t green_atom = -1;
    atom_t yellow_atom = -1;
    atom_t magenta_atom = -1;
    atom_t cyan_atom = -1;
    char buffer[48];
    int i;

    attr_list al = create_attr_list();
    attr_list al2 = create_attr_list();

    set_attr (al, 2000, Attr_String, (attr_value)"test string");
    set_attr (al, 2001, Attr_Int4, (attr_value)2001);

    attr_add_list (al, al2);
    attr_list_to_string (al);

    dump_attr_list (al);
    free_attr_list(al);
    free_attr_list(al2);
    
    for (i=0; i < sizeof(buffer); i++) {
        buffer[i] = 0x30 + i;
    }
    red_atom = attr_atom_from_string("Red");
    ip_atom = attr_atom_from_string("IP_ADDR");
    blue_atom = attr_atom_from_string("Red Stripe");
    green_atom = attr_atom_from_string("Green");
    yellow_atom = attr_atom_from_string("Yellow");
    magenta_atom = attr_atom_from_string("Magenta");
    cyan_atom = attr_atom_from_string("Cyan");

    add_attr(list, ip_atom, Attr_Int4, (attr_value) 
	     (((unsigned int)130<<24) + ((unsigned int) 207<<16) +((unsigned int)5 << 8) + 68));
    add_attr(list, red_atom, Attr_Int4, (attr_value) 2);
    add_attr(list, blue_atom, Attr_String, (attr_value) strdup("The sky"));
    add_double_attr(list, magenta_atom, 3.14159);
    add_opaque_attr(list, cyan_atom, 48, &buffer[0]);
    dump_attr_list(list);
/*    add_attr(list2, RED_ATOM, Attr_Int4, (attr_value) 4);*/
    attr_add_list(list, list2);
    dump_attr_list(list);
    printf("stringified is %s\n", attr_list_to_string(list2));
    putenv("ATL_BASE64_STRINGS=1");
    string = attr_list_to_string(list);
    printf("stringified version is >%s<\n", string);
    free_attr_list(list);
    list2 = attr_list_from_string(string);
    free(string);
    dump_attr_list(list2);
    printf("\n");
    list = attr_copy_list(list2);
    free_attr_list(list2);
    dump_attr_list(list);
    printf("\n");
    free_attr_list(list);

    list = create_attr_list();
    list2 = create_attr_list();

    add_attr(list, red_atom, Attr_Int4, (attr_value) 2);
    add_attr(list, green_atom, Attr_Int4, (attr_value) 20);
    add_attr(list, blue_atom, Attr_String, (attr_value) strdup("The sky"));
    add_attr(list, yellow_atom, Attr_String, (attr_value) strdup("The sun"));

    add_attr(list2, green_atom, Attr_Int4, (attr_value) 20);
    add_attr(list2, red_atom, Attr_Int4, (attr_value) 2);
    add_attr(list2, yellow_atom, Attr_String, (attr_value) strdup("The sun"));
    add_attr(list2, blue_atom, Attr_String, (attr_value) strdup("The sky"));

    dump_attr_list(list);
    printf("\n");
    dump_attr_list(list2);
    printf("\n");

    string = attr_list_to_string(list);
    string2 = attr_list_to_string(list2);
    
    printf ("string [ %s ]\n", string);
    printf ("string2 [ %s ]\n", string2);

    printf ("strcmp returned %d\n", strcmp(string,string2));
    free_attr_list(list2);
    free_attr_list(list);
    free(string); free(string2);

    list = create_attr_list();
    list2 = create_attr_list();

    add_attr(list, 130, Attr_Int4, (attr_value) 130);
    add_attr(list, 120, Attr_Int4, (attr_value) 120);
    add_attr(list, 110, Attr_Int4, (attr_value) 110);
    add_attr(list, 230, Attr_String, (attr_value) strdup("The sky"));
    add_attr(list, 220, Attr_String, (attr_value) strdup("The sun"));
    add_attr(list, 210, Attr_String, (attr_value) strdup("The clouds"));

    add_attr(list2, 130, Attr_Int4, (attr_value) 130);
    add_attr(list2, 110, Attr_Int4, (attr_value) 110);
    add_attr(list2, 120, Attr_Int4, (attr_value) 120);

    add_attr(list2, 230, Attr_String, (attr_value) strdup("The sky"));
    add_attr(list2, 210, Attr_String, (attr_value) strdup("The clouds"));
    add_attr(list2, 220, Attr_String, (attr_value) strdup("The sun"));

    dump_attr_list(list);
    printf("\n");
    dump_attr_list(list2);
    printf("\n");

    string = attr_list_to_string(list);
    string2 = attr_list_to_string(list2);
    
    printf ("string [ %s ]\n", string);
    printf ("string2 [ %s ]\n", string2);

    printf ("strcmp returned %d\n", strcmp(string,string2));
    free_attr_list(list2);
    free_attr_list(list);

    /*    list = create_attr_list();
    list2 = create_attr_list();
    */
    list = attr_list_from_string(string);
    list2 = attr_list_from_string(string2);

    dump_attr_list(list);
    printf("\n");
    dump_attr_list(list2);
    printf("\n");

    free_attr_list(list2);
    free_attr_list(list);
    free(string); free(string2);
    return 0;
}
