#include "fm.h"
#ifdef _MSC_VER
#define strdup(s) _strdup(s)
#endif
extern void
write_buffer(char *filename, FMStructDescList desc, void *data, 
             int test_num);
extern char *read_buffer(FMContext c, char *read_file, int test_num);

