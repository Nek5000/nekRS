extern double drand48();
extern long lrand48();
#include "../../ffs.h"
typedef struct _first_rec {
    int         i;
    long        j;
    double      d;
    char        c;
} first_rec, *first_rec_ptr;

static FMField field_list[] = {
    {"i", "integer", sizeof(int), FMOffset(first_rec_ptr, i)},
    {"j", "integer", sizeof(long), FMOffset(first_rec_ptr, j)},
    {"d", "float",   sizeof(double), FMOffset(first_rec_ptr, d)},
    {"c", "integer", sizeof(char), FMOffset(first_rec_ptr, c)},
    {NULL, NULL, 0, 0},
};

typedef struct _second_rec {
    int         *items;
    long        count;
    char *      name;
} second_rec, *second_rec_ptr;

static FMField second_field_list[] = {
    {"items", "integer[count]", sizeof(int), FMOffset(second_rec_ptr, items)},
    {"count", "integer", sizeof(long), FMOffset(second_rec_ptr, count)},
    {"name", "string",   sizeof(char*), FMOffset(second_rec_ptr, name)},
    {NULL, NULL, 0, 0},
};

int main(int argc, char **argv)
{
    FFSFile file = open_FFSfile("test_output2", "w");
    FMContext fmc = FMContext_of_file(file);
    FMFormat rec_ioformat, rec_ioformat2;
    first_rec rec1;
    second_rec rec2;
    int i;

    srand48(time());
    rec_ioformat = register_simple_format(fmc, "first format", field_list, sizeof(first_rec));
    rec_ioformat2 = register_simple_format(fmc, "second format", second_field_list, sizeof(second_rec));
    for (i=0; i < 10; i++) {
	if ((lrand48() % 2) == 1) {
	    rec1.i = lrand48()%100;
	    rec1.j = lrand48()%100;
	    rec1.d = drand48();
	    rec1.c = lrand48()%26+'A';
	    write_FFSfile(file, rec_ioformat, &rec1);
	    printf("Wrote i=%d, j=%ld, d=%g, j=%c\n", rec1.i, rec1.j, rec1.d, rec1.c);
	} else {
	    char str[32];
	    int array[10];
	    int j;
	    sprintf(str, "random string %d", (int)lrand48()%100);
	    rec2.name = &str[0];
	    rec2.count = lrand48()%10;
	    rec2.items = &array[0];
	    for(j=0; j<rec2.count; j++) rec2.items[j] = lrand48()%5+2;
	    write_FFSfile(file, rec_ioformat2, &rec2);
	    printf("Wrote name=%s, count=%ld, items = [", rec2.name, rec2.count);
	    for(j=0; j<rec2.count; j++) printf("%d, ", rec2.items[j]);
	    printf("]\n");
	}
    }	    
    close_FFSfile(file);
}
