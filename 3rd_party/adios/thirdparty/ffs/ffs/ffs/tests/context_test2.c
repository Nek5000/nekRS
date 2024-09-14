
#include "config.h"
#include <fcntl.h>
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>
#include <assert.h>
#include "ffs.h"

#include "test_funcs.h"
#include <sys/types.h>
#include <sys/stat.h>

static void test_receive(char *buffer, size_t buf_size, int finished,
			       int test_level);
static void test_all_receive(char *buffer, size_t buf_size, int finished);
static void write_buffer(char *buf, size_t size);
static void read_test_only();

static int write_output = 0;
static char *output_file = NULL;
static char *read_file = NULL;
static int fail = 0;
static char *test_only = NULL;

static FFSContext rcv_context = NULL;
static int verbose = 0;
static FMContext loaded_FMcontext = NULL;

int
main(int argc, char **argv)
{
    FMContext src_context;
    FFSBuffer encode_buffer;
    char *xfer_buffer;
    size_t buf_size;
    int i;
    struct node nodes[10];
    struct visit_table v;
    FMFormat node_ioformat, psa_ioformat, reader_register_ioformat;

    init_written_data();

    for (i=1; i<argc; i++) {
	if (strcmp(argv[i], "-w") == 0) {
	    output_file = argv[++i];
	    write_output++;
	} else if (strcmp(argv[i], "-r") == 0) {
	    read_file = argv[++i];
	} else if (strcmp(argv[i], "-v") == 0) {
	    verbose++;
	} else if (strcmp(argv[i], "-t") == 0) {
	    test_only = argv[++i];
	} else {
	    printf("Unknown argument %s\n", argv[i]);
	    printf("Usage:\n\t-w\t  write test output\n");
	    printf("\t-r file\t  rest/process test data in <file>\n");
	    printf("\t-t test\t  test only format <test>\n");
	    exit(1);
	}
    }

    if (read_file) {
	read_test_only();
	free_written_data();
	if (rcv_context != NULL) {
	    free_FFSContext(rcv_context);
	    rcv_context = NULL;
	}
	if (loaded_FMcontext) {
	    free_FMcontext(loaded_FMcontext);
	}
	if (fail) exit(1);
	exit(0);
    }
    src_context = create_FMcontext();
    encode_buffer = create_FFSBuffer();

    node_ioformat = register_data_format(src_context, node_format_list);
    psa_ioformat = register_data_format(src_context, pointer_to_static_format_list);
    reader_register_ioformat = register_data_format(src_context, reader_register_format_list);

    for (i = 0; i < sizeof(nodes)/sizeof(nodes[0]); i++) {
	nodes[i].node_num = i;
	nodes[i].link1 = nodes[i].link2 = NULL;
    }

    for (i=0; i <  sizeof(nodes)/sizeof(nodes[0]) - 1; i++) {
	nodes[i].link1 = &nodes[i+1];
    }
/*    nodes[0].link2 = &nodes[sizeof(nodes)/sizeof(nodes[0])-1];*/
    v.node_count = 0;
    nodes[0].node_num = calc_signature(&nodes[0], &v);
    xfer_buffer = FFSencode(encode_buffer, node_ioformat,
					  &nodes[0], &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);

    nodes[0].link2 = NULL;
    nodes[sizeof(nodes)/sizeof(nodes[0]) - 1].link1 = &nodes[2];
    v.node_count = 0;
    nodes[0].node_num = 0;
    nodes[0].node_num = calc_signature(&nodes[0], &v);
    xfer_buffer = FFSencode(encode_buffer, node_ioformat,
					  &nodes[0], &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);

    for (i=0; i <  sizeof(nodes)/sizeof(nodes[0]) - 1; i++) {
	nodes[i].link1 = nodes[i].link2 = NULL;
    }
    nodes[0].link1 = &nodes[1];
    nodes[0].link2 = &nodes[2];
    nodes[1].link1 = &nodes[3];
    nodes[1].link2 = &nodes[4];
    nodes[2].link1 = &nodes[5];
    nodes[2].link2 = &nodes[6];
    nodes[3].link1 = &nodes[7];
    nodes[3].link2 = &nodes[8];
    nodes[4].link1 = &nodes[9];

    v.node_count = 0;
    nodes[0].node_num = 0;
    nodes[0].node_num = calc_signature(&nodes[0], &v);
    xfer_buffer = FFSencode(encode_buffer, node_ioformat,
					  &nodes[0], &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    
    xfer_buffer = FFSencode(encode_buffer, psa_ioformat,
					  &psa, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);

    xfer_buffer = FFSencode(encode_buffer, reader_register_ioformat,
					  &reader_register, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    

    free_FMcontext(src_context);
    free_FFSBuffer(encode_buffer);
    src_context = NULL;
    test_all_receive(NULL, 0, 1);
    write_buffer(NULL, 0);

    free_written_data();
    if (rcv_context != NULL) {
	free_FFSContext(rcv_context);
	rcv_context = NULL;
    }
    if (loaded_FMcontext) {
	free_FMcontext(loaded_FMcontext);
    }
    if (fail) exit(1);
    return 0;
}

/* NT needs O_BINARY, but it doesn't exist elsewhere */
#ifndef O_BINARY
#define O_BINARY 0
#endif

static char *
get_buffer(size_t *size_p)
{
    static int file_fd = 0;
    static char *buffer = NULL;
    char *tmp_buffer;
    static int last_size = -1;
    unsigned short ssize;
    int to_read;
    int tmp_size;
    unsigned char csize;
    unsigned int beef = 0xdeadbeef;

    if (read_file == NULL) exit(1);

    if (loaded_FMcontext == NULL) {
	loaded_FMcontext = create_local_FMcontext();
    }
    if (file_fd == 0) {
	file_fd = open(read_file, O_RDONLY|O_BINARY, 0777);
	buffer = malloc(1);
    }
    if (last_size != -1) {
	if (memcmp(buffer+last_size, &beef, 4) != 0) {
	    printf("memory overwrite error\n");
	}
    }
    if (read(file_fd, &csize, 1) != 1) exit(1);	/* low byte of 2-byte size */
    ssize = csize;
    if (read(file_fd, &csize, 1) != 1) exit(1);	/* high byte of 2-byte size */
    ssize += ((csize << 8) & 0xff00);
    to_read = ssize;
    buffer = realloc(buffer, to_read+4);
    tmp_buffer = buffer;
    while((tmp_size = read(file_fd, tmp_buffer, to_read)) != to_read) {
	if (tmp_size == 0) {
	    free(buffer);
	    return NULL;
	} else if (tmp_size == -1) {
	    perror("Read failure");
	    free(buffer);
	    return NULL;
	}
	to_read -= tmp_size;
	tmp_buffer += tmp_size;
    } 
    last_size = ssize;
    memcpy(buffer+last_size, &beef, 4);
    if (ssize == 0) {
	free(buffer);
	close(file_fd);
	file_fd = 0;
	return NULL;
    } else {
	*size_p = ssize;
	return buffer;
    }
}

static void
read_test_only()
{
    char *input;
    size_t size;
    while ((input = get_buffer(&size)) != NULL) {
	test_all_receive(input, size, 0);
    }
    test_all_receive(NULL, 0, 1);
    free(input);
    if (rcv_context != NULL) {
	free_FFSContext(rcv_context);
	rcv_context = NULL;
    }
}

static FFSTypeHandle node_ioformat;
static FFSTypeHandle psa_ioformat;
static FFSTypeHandle reader_register_ioformat;

static void
set_targets(FFSContext context)
{
    node_ioformat = FFSset_fixed_target(context, node_format_list);
    psa_ioformat = FFSset_fixed_target(context, pointer_to_static_format_list);
    reader_register_ioformat = FFSset_fixed_target(context, reader_register_format_list);
}

size_t base_size_func(FFSContext context, char *src, size_t rec_len,
		   size_t native_struct_size)
{
    return native_struct_size;
}

size_t total_size_func(FFSContext context, char *src, size_t rec_len, 
		    size_t native_struct_size)
{
    return FFS_est_decode_length(context, src, rec_len);
}

static int 
decode_in_place(FFSContext context, char *src, size_t src_len, void *dest)
{
    if (FFSdecode_in_place_possible(FFSTypeHandle_from_encode(context, src))) {
	int ret, header_len;
	char *real_dest;
	ret = FFSdecode_in_place(context, src, (void**)&real_dest);
	header_len = (int)(intptr_t)(real_dest - src);
	memcpy(dest, real_dest, src_len - header_len);
	return ret;
    } else {
	return FFSdecode_to_buffer(context, src, dest);
    }
}

static int
decode_IOcontext_wrapper(FFSContext context, char *src, size_t src_len, void *dest)
{
    return FFSdecode(context, src, dest);
}

static int
decode_to_buffer_IOcontext_wrapper(FFSContext context, char *src, size_t src_len,
				   void *dest)
{
    return FFSdecode_to_buffer(context, src, dest);
}

typedef size_t (*size_func_t)(FFSContext context, char *src, size_t buf_size, 
				   size_t nativ_struct);

typedef int (*decode_func_t)(FFSContext context, char *src, size_t src_len, 
				   void *dest);

size_func_t size_funcs[] = {base_size_func, total_size_func, total_size_func};
decode_func_t decode_funcs[] = {decode_IOcontext_wrapper, 
				decode_to_buffer_IOcontext_wrapper,
				decode_in_place};


#define NUM_TESTS 3

static void
test_all_receive(char *buffer, size_t buf_size, int finished)
{
    int test_type = 0;
    char *tmp_buffer = malloc(buf_size);
    for ( test_type = 0; test_type < NUM_TESTS; test_type++) {
	memcpy(tmp_buffer, buffer, buf_size);
	test_receive(tmp_buffer, buf_size, finished, test_type);
    }
    free(tmp_buffer);
}
	
static void*
get_mem(size_t size)
{
    char *buffer;
    unsigned int beef = 0xdeadbeef;

    buffer = malloc(size+4);
    memcpy(buffer+size, &beef, 4);
    return buffer;
}

static void
check_mem(size_t size, char *buffer)
{
    unsigned int beef = 0xdeadbeef;
    if (memcmp(buffer+size, &beef, 4) != 0) {
	printf("memory overwrite error\n");
    }
}

    
static void
test_receive(char *buffer, size_t buf_size, int finished, int test_level)
{
    static FFSContext c = NULL;
    static int node_rec_count[NUM_TESTS] = {0,0,0};

    static int unknown_rec_count[NUM_TESTS] = {0,0,0};
    size_func_t size_func = size_funcs[test_level];
    decode_func_t decode_func = decode_funcs[test_level];
    if (c == NULL) {
	c = create_FFSContext();
	rcv_context = c;
	set_targets(rcv_context);
    }
    if (!finished) {
/*	char *comment;*/
	FFSTypeHandle buffer_format = FFS_target_from_encode(rcv_context, buffer);

	if (buffer_format == NULL) {
	    printf("No format!\n");
	    exit(1);
	}
	if (((test_only == NULL) || (strcmp(test_only, "node") == 0)) &&
		   (buffer_format == node_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, 
				 sizeof(add_rec));
	    node_ptr read_data = get_mem(size);
	    struct visit_table v;
	    int signature, calculated_sig;
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, decode msg format\n");
	    signature = read_data->node_num;
	    read_data->node_num = 0;
	    v.node_count = 0;
	    calculated_sig = calc_signature(read_data, &v);
	    if (signature != calculated_sig) {
		printf("Node sig not right, calculated %d, expected %d\n", 
		       calculated_sig, signature);
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	    node_rec_count[test_level]++;
	} else 	if (((test_only == NULL) || (strcmp(test_only, "psa") == 0)) &&
		   (buffer_format == psa_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, 
				 sizeof(add_rec));
	    pointer_to_static_array_ptr read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, pointer_to_static_array data failed\n");
	    if (!pointer_to_static_rec_eq(read_data, &psa)) {
		printf("Pointer to static array failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else 	if (((test_only == NULL) || (strcmp(test_only, "reader_register") == 0)) &&
		   (buffer_format == reader_register_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, 
				 sizeof(add_rec));
	    struct _reader_register_msg *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, pointer_to_static_array data failed\n");
	    if (!reader_register_rec_eq(read_data, &reader_register)) {
		printf("Pointer to reader_register failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (test_only == NULL) {
	    printf("discarding a record\n");
	    unknown_rec_count[test_level]++;
	}
    } else {
	/* finished */
	if (test_level == 0) {
	    free_FFSContext(rcv_context);
	    rcv_context = NULL;
	}
	if (test_only != NULL) return;
	if (node_rec_count[test_level] != 3) {
	    printf("Missed node, level %d\n", test_level);
	    fail++;
	}
	if (unknown_rec_count[test_level] != 0) {
	    printf("Got unknown\n");
	    fail++;
	}
    }
}

static void
write_buffer(char *buf, size_t size)
{
    static int file_fd = 0;
    unsigned short ssize;
    unsigned char csize;
    if (output_file == NULL) return;

    if (file_fd == 0) {
	file_fd = open(output_file, O_WRONLY|O_CREAT|O_TRUNC|O_BINARY, 0777);
    }
    printf("Writing buffer of size %zd\n", size);
    ssize = (unsigned short) size;
    csize = ssize & 0xff;
    if (write(file_fd, &csize, 1) != 1) exit(1);/* low byte of 2-byte size */
    csize = ((ssize >> 8) & 0xff);
    if (write(file_fd, &csize, 1) != 1) exit(1);/* high byte of 2-byte size */
    if (write(file_fd, buf, (int)size) != size) exit(1);
}
