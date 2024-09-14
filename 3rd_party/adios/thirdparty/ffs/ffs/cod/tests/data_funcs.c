#include "config.h"
#include "fm.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _MSC_VER
#define read(x,y,z) _read(x,y,z)
#define open(x,y,z) _open(x,y,z)
#include <io.h>
#endif

#define MAGIC 0x4356ffaa	/* random magic */
#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))
void
write_buffer(char *filename, FMStructDescList desc, void *data, 
             int test_num)
{
/*     static FMContext c = NULL; */
/*     FMFormat ioformat = NULL; */
/*     char *enc_buf; */
/*     int buf_len; */
/*     char *server_id; */
/*     int id_len; */
/*     char *server_rep; */
/*     int rep_len; */

/*     int i = 0; */
/*     int magic = MAGIC; */
/*     int tmp = test_num; */
/*     static int f = 0; */
/*     if (f == 0) { */
/* 	f = open(filename, O_WRONLY|O_CREAT, 0777); */
/*         write(f, &magic, 4); */
/*     } */
/*     write(f, &tmp, 4); */
/*     if (c == NULL) { */
/* 	c = create_local_IOcontext(); */
/*     } */
/*     ioformat = register_data_format(c, desc); */
/*     server_id = get_server_ID_FMformat(ioformat, &id_len); */
/*     server_rep = get_server_rep_FMformat(ioformat, &rep_len); */
/*     write(f, &id_len, 4); */
/*     write(f, &server_id[0], id_len); */
/*     write(f, &rep_len, 4); */
/*     write(f, &server_rep[0], rep_len); */
/*     top_ioformat = ioformat; */

/*     enc_buf = encode_IOcontext_buffer(c, top_ioformat, data, &buf_len); */
/*     write(f, &buf_len, 4); */
/*     write(f, &enc_buf[0], buf_len); */
}

static void
byte_swap(char *data, int size)
{
    int i;
    assert((size % 2) == 0);
    for (i = 0; i < size / 2; i++) {
	char tmp = data[i];
	data[i] = data[size - i - 1];
	data[size - i - 1] = tmp;
    }
}

static
int read_data_block(int f, FMContext c, int swap, char **buf_p)
{
/*     int test_num = 0; */
/*     int format_count = 0; */
/*     int len = 0; */
/*     int i = 0; */
/*     if (read(f, &test_num, 4) != 4) return -1; */
/*     if (swap) byte_swap(&test_num, 4); */
/*     if (read(f, &format_count, 4) != 4) return -1; */
/*     if (swap) byte_swap(&format_count, 4); */
/*     for (i=0; i < format_count; i++) { */
/* 	int id_len; */
/* 	IOFormat ioformat; */
/* 	char *id; */
/* 	char *body; */
/* 	read(f, &id_len, 4); */
/* 	if (swap) byte_swap(&id_len, 4); */
/* 	id = malloc(id_len); */
/* 	read(f, &id[0], id_len); */
/* 	read(f, &len, 4); */
/* 	if (swap) byte_swap(&len, 4); */
/* 	body = malloc(len); */
/* 	read(f, &body[0], len); */
/* 	ioformat = load_external_format_IOcontext(c, id, id_len, body); */
/*     } */
/*     read(f, &len, 4); */
/*     if (swap) byte_swap(&len, 4); */
/*     *buf_p = malloc(len); */
/*     read(f, *buf_p, len); */
/*     return test_num; */
  return 0;
}
	
extern
char *read_buffer(FMContext c, char *read_file, int test_num)
{
    static int f = 0;
    static int swap = 0;
    static int data_test_num = -1;
    char *buf;
    if (f == 0) {
	int in_magic = 0;
	f = open(read_file, O_RDONLY, 0777);
	if(read(f, &in_magic, 4) != 4) exit(1);
	if (in_magic != MAGIC) {
	    byte_swap((char*)&in_magic, 4);
	    swap++;
	}
	if (in_magic != MAGIC) {
	    printf("Bad magic number in file, %x\n", in_magic);
	    exit(1);
	}
    }
    data_test_num = read_data_block(f, c, swap, &buf);
    while (data_test_num != test_num) {
	if (data_test_num == -1) {
	    printf("Data for test %d not in file\n", test_num);
	    exit(1);
	}
	data_test_num = read_data_block(f, c, swap, &buf);
	if (data_test_num == -1 ) {
	    printf("End of file error\n");
	}
    }
    return buf;
}
