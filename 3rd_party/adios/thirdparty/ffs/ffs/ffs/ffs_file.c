
#include "config.h"

#include "assert.h"
#include "atl.h"
#include "ffs.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#ifndef _MSC_VER
#include "unistd.h"
#include <arpa/inet.h>
#else
#include <winsock.h>
#include <io.h>
#endif
#include "errno.h"
#include "string.h"
#include "sys/types.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "fm_internal.h"
#include "ffs_internal.h"
#include "io_interface.h"

typedef struct _CDLLnode {
    void *data;
    struct _CDLLnode *next;
    struct _CDLLnode *prev;
}CDLLnode;

typedef enum {
    OpenForRead, OpenForWrite, Closed
} Status;

typedef enum {
    Simple, Indexed
} Organization;

typedef struct format_info {
    int written_to_file;
} format_info;

typedef struct write_index_info {
    off_t base_file_pos;
    int data_index_start;
    int data_index_end;
    int index_block_size;
    int next_item_offset;
    unsigned char *index_block;
} write_index_info;

typedef union {
    write_index_info write_info;
} *FFSIndexType;

typedef struct _FFSIndexElement {
    FFSRecordType type;
    long fpos;
    char *format_id;
    int fid_len;
    attr_list attrs;
} FFSIndexElement;

typedef struct _FFSIndexItem {
    long next_index_offset;
    long this_index_fpos;
    int start_data_count;
    int last_data_count;
    int end_item_offset;
    int elem_count;
    struct _FFSIndexElement *elements;
    struct _FFSIndexItem *next;
    struct _FFSIndexItem *prev;
} FFSIndexItemStruct;

struct _FFSFile {
    FFSContext c;
    FMContext fmc;

    Organization  file_org;
    FFSBuffer tmp_buffer;
    void *file_id;
    format_info *info;
    int info_size;
    int next_fid_len;
    DATA_LEN_TYPE next_data_len;
    FFSBuffer buf;
    int read_ahead;
    int errno_val;
    int raw_flag;
    FFSRecordType next_record_type;
    FFSTypeHandle next_data_handle;
    FFSTypeHandle next_actual_handle;

    CDLLnode *attr_list;
    long data_block_no;

    off_t fpos;
    int data_count;
    FFSIndexType cur_index;
    FFSIndexItem read_index;
    int visible_items_bitmap;

    FFSIndexItem index_head;
    FFSIndexItem index_tail;

    Status status;
    IOinterface_func write_func;
    IOinterface_func read_func;
    int max_iov;
    IOinterface_funcv writev_func;
    IOinterface_funcv readv_func;
    IOinterface_close close_func;

};

#define INDEX_BLOCK_SIZE 256

#define NUM_OFFSETS 10

typedef struct _offset_node {
    int data_offset[NUM_OFFSETS];
}offset_node;

typedef struct _attr_node {
    char *attr;
    CDLLnode *offsetLL;
    int len_offsetLL;
    int num_offsets;
}attr_node;


/*
 * Generic circular doubly linked list routine
 * returns 0 on success and 1 otherwise
 */
int
CDLLadd(CDLLnode **start, void *data, unsigned int data_size)
{
    CDLLnode *nodeP;
    CDLLnode *startP = *start;
    int ret = 1;

    nodeP = (CDLLnode *)malloc(sizeof(CDLLnode));
    if (!nodeP)
        goto exit;

    nodeP->data = malloc(data_size);
    if (!(nodeP->data))
        goto free_nodeP;

    memcpy(nodeP->data, data, data_size);
    nodeP->prev = NULL;
    nodeP->next = NULL;

    if (!startP) {
        nodeP->next = nodeP;
        nodeP->prev = nodeP;
        *start = nodeP;
    }
    else {
        nodeP->next = startP;
        startP->prev->next = nodeP;
        nodeP->prev = startP->prev;
        startP->prev = nodeP;
    }

    ret = 0;
    goto exit;

free_nodeP:
    free(nodeP);

exit:
    return ret;
}


/*
 * Generic circular doubly linked list routine
 * returns 0 on success and 1 otherwise
 */
int
CDLLremove(CDLLnode **start, CDLLnode *remove_nodeP)
{
    CDLLnode *startP = *start;
    CDLLnode *prev_rem_nodeP = remove_nodeP->prev;
    CDLLnode *next_rem_nodeP = remove_nodeP->next;

    if (!startP || !remove_nodeP)
        goto exit;

    if (remove_nodeP->data)
        free(remove_nodeP->data);

    if (startP == remove_nodeP) {
        if (startP->next == startP) {
            free(remove_nodeP);
            *start = NULL;
        }
        else {
            free(remove_nodeP);
            prev_rem_nodeP->next = next_rem_nodeP;
            next_rem_nodeP->prev = prev_rem_nodeP;
            *start = next_rem_nodeP;
        }
    }
    else {
        free(remove_nodeP);
        prev_rem_nodeP->next = next_rem_nodeP;
        next_rem_nodeP->prev = prev_rem_nodeP;
    }

exit:
    return 0;
}

static void read_all_index_and_formats(FFSFile file);
static void convert_last_index_block(FFSFile file);


static FFSRecordType next_record_type(FFSFile ffsfile);

static void
update_fpos(FFSFile f)
{
    f->fpos = ffs_file_lseek_func(f->file_id, 0, SEEK_CUR);
}

static void
free_FFSIndexItem(FFSIndexItemStruct *item)
{
    int i;
    for (i = 0; i < item->elem_count; i++) {
	if (item->elements[i].format_id) {
	    free(item->elements[i].format_id);
	}
	if (item->elements[i].attrs) {
	    free_attr_list(item->elements[i].attrs);
	}
    }
    if (item->elements) free(item->elements);
    free(item);
}

extern
void
free_FFSfile(FFSFile f)
{
    FFSIndexItemStruct *i;
    free(f->info);
    f->info = NULL;
    f->info_size = 0;
    if (f->buf) free_FFSBuffer(f->buf);
    if (f->tmp_buffer) free_FFSBuffer(f->tmp_buffer);
    if (f->cur_index) {
	if (f->cur_index->write_info.index_block) {
	    free(f->cur_index->write_info.index_block);
	}
	free(f->cur_index);
    }
    f->buf = NULL;
    i = f->index_head;
    while (i != NULL) {
	FFSIndexItemStruct *next = i->next;
	free_FFSIndexItem(i);
	i = next;
    }
    free_FFSContext(f->c);
    free_FMcontext(f->fmc);
    free(f);
}


static void
parse_flags(const char *flags, int *allow_input_p, int *allow_output_p, 
	    int *raw_p, int *index_p)
{
    int input = -1;
    int output = -1;
    int raw = -1;
    int index = -1;

    if (flags == NULL) return;

    while (flags[0] != 0) {
	switch(flags[0]) {
	case 'R':
	    raw = 1;
	    /* falling through */
	case 'r':
	    if (output == 1) {
		printf("Warning, read flag specified after write flag\n");
	    }
	    if (index != -1) {
		printf("Cannot specify index on input files\n");
		index = 0;
	    }
	    input = 1;
	    output = 0;
	    break;
	case 'i':
	    index = 1;
	    break;
	case 'w':
	    if (input == 1) {
		printf("Warning, write flag specified after read flag\n");
	    }
	    output = 1;
	    input = 0;
	    break;
	case 'a':
	    if ((input == 1) || (output == 1)) {
		printf("Warning, append flag specified after read flag\n");
	    }
	    output = 1;
	    input = 1;
	    break;
	default:
	    printf("Warning, unknown flag in FFS open call, \'%c\'\n", 
		   flags[0]);
	}
	flags++;
    }
    if (index == -1) index = 0;
    *index_p = index;
    if (raw == -1) raw = 0;
    *raw_p = raw;
    if (input == -1) input = 0;
    *allow_input_p = input;
    if (output == -1) output = 0;
    *allow_output_p = output;
}

extern FFSFile
open_FFSfd(void *fd, const char *flags)
{
    void *file = fd;
    FFSFile f;
    int allow_input = 0, allow_output = 0;
    int raw = 0, index = 0;

    f = malloc(sizeof(struct _FFSFile));
    memset(f, 0, sizeof(*f));
    f->file_id = file;
    f->visible_items_bitmap = FFSend|FFSerror|FFSdata|FFSformat|FFScomment;

    parse_flags(flags, &allow_input, &allow_output, &raw, &index);

    f->file_org = Simple;
    if (index) f->file_org = Indexed;
    if (raw) f->raw_flag = 1;

    f->attr_list = NULL;
    f->data_block_no = 0;

    set_interface_FFSFile(f, ffs_file_write_func, ffs_file_read_func,
                          ffs_file_writev_func, ffs_file_readv_func,
                          ffs_max_iov, ffs_close_func);

    f->buf = create_FFSBuffer();
    f->status = OpenForRead;
    f->fmc = create_local_FMcontext();
    f->c = create_FFSContext_FM(f->fmc);
    if (allow_input) {
	int magic_number;
	int count;

	count = f->read_func(f->file_id, &magic_number, 4, NULL, NULL);
	if ((count <= 0) && allow_input && allow_output) {
	    /* newly created file opened for append, just do write */
	    allow_input = 0;
	} else if ((count != 4) || (magic_number != htonl(MAGIC_NUMBER))) {
	    printf("read headers failed\n");
	    return NULL;
	}
	f->status = OpenForRead;
    }
    if (allow_input && allow_output) {   /* append mode! */
	read_all_index_and_formats(f);
	convert_last_index_block(f);
	f->status = OpenForWrite;
    } else if (allow_output) {
	int magic_number = htonl(MAGIC_NUMBER);
	assert(sizeof(int) == 4);	/* otherwise must do other stuff */

	if (f->write_func(file, &magic_number, 4, NULL, NULL) != 4) {
	    printf("write headers failed\n");
	    return NULL;
	}
	update_fpos(f);
	f->status = OpenForWrite;
    }
    return f;
}

extern FFSFile
open_FFSfile(const char *path, const char *flags)
{
    void *file;
    FFSFile f;
    int allow_input = 0, allow_output = 0;
    int raw, index;

    parse_flags(flags, &allow_input, &allow_output, &raw, &index);

    if (allow_input && allow_output) {
	file = (FFSFile) (ffs_file_open_func)(path, "a", NULL, NULL);
	if (file == (void*)0) {
	    /* if open for append failed, try creating it */
	    file = (FFSFile)(ffs_file_open_func)(path, "w", NULL, NULL);
	}
    } else if (allow_input) {
	file = (FFSFile)(ffs_file_open_func)(path, "r", NULL, NULL);
    } else {
	file = (FFSFile)(ffs_file_open_func)(path, "w", NULL, NULL);
    }

    if (file == NULL) {
	char msg[128];
	(void) snprintf(msg, sizeof(msg), "open_FFSfile failed for %s :", path);
	perror(msg);
	return NULL;
    }
    f = open_FFSfd(file, flags);
    return f;
}

FFSContext
FFSContext_of_file(FFSFile f)
{
    return f->c;
}

FMContext
FMContext_of_file(FFSFile f)
{
    return f->fmc;
}

extern void
set_interface_FFSFile(FFSFile f, IOinterface_func write_func, IOinterface_func read_func, IOinterface_funcv writev_func,
                      IOinterface_funcv readv_func, int max_iov, IOinterface_close close_func)
{
    f->write_func = write_func;
    f->read_func = read_func;
    f->max_iov = max_iov;
    f->writev_func = writev_func;
    f->readv_func = readv_func;
    f->close_func = close_func;
}

static void dump_index_block(FFSFile f);

extern void
close_FFSfile(FFSFile file)
{
    if (file->file_org == Indexed) dump_index_block(file);
    file->close_func(file->file_id);
}

void
init_format_info(FFSFile f, int index)
{
    if (f->info == NULL) {
	f->info = malloc(sizeof(f->info[0]) * (index+1));
	memset(f->info, 0, sizeof(f->info[0]) * (index+1));
	f->info_size = index + 1;
    } else if (f->info_size <= index) {
	f->info = realloc(f->info,
				 (index+1)*sizeof(f->info[0]));
	memset(&f->info[f->info_size], 0, 
	       sizeof(f->info[0]) * ((index+1) - f->info_size));
	f->info_size = index + 1;
    }
}

static void
init_write_index_block(FFSFile f)
{
    int data_index_start = 0;
    off_t end_of_index;
    if (f->read_index == NULL) { /* if not append */
	end_of_index = ffs_file_lseek_func(f->file_id, INDEX_BLOCK_SIZE, SEEK_CUR);
	if (f->cur_index) {
	    data_index_start = f->cur_index->write_info.data_index_end;
	} else {
	    f->cur_index = malloc(sizeof(*(f->cur_index)));
	    memset(f->cur_index, 0, sizeof(*(f->cur_index)));
	}
    
	f->cur_index->write_info.base_file_pos = end_of_index - INDEX_BLOCK_SIZE;
	f->cur_index->write_info.data_index_start = data_index_start;
	f->cur_index->write_info.data_index_end = f->cur_index->write_info.data_index_start;
	f->cur_index->write_info.index_block_size = INDEX_BLOCK_SIZE;
	if (!f->cur_index->write_info.index_block) {
	    f->cur_index->write_info.index_block = malloc(INDEX_BLOCK_SIZE);
	    memset(f->cur_index->write_info.index_block, 0, INDEX_BLOCK_SIZE);
	}
	f->cur_index->write_info.next_item_offset = 16;   /* number of bytes written below */
	f->fpos = end_of_index;
    } else {
	f->read_index = NULL;
    }
}

static void output_index_end(FFSFile f);

static void
dump_index_block(FFSFile f)
{
    off_t end = ffs_file_lseek_func(f->file_id, 0, SEEK_CUR);
    int ret;

    size_t size =  f->cur_index->write_info.index_block_size;
    unsigned char *index_base = f->cur_index->write_info.index_block;

    output_index_end(f);
    f->cur_index->write_info.data_index_end = f->data_count-1;
    ffs_file_lseek_func(f->file_id, f->cur_index->write_info.base_file_pos, SEEK_SET);
    /*
     * next_data indicator is a 2 4-byte chunks in network byte order.
     * In the first chunk, 
     *    the top byte is 0x4, next three are length of the index block.
     */
    *((int*)index_base) = (int) htonl((0x4<<24) | (long)size);
    *((int*)(index_base + 4)) = (int)htonl(end);  /* link to next index */
    *((int*)(index_base+8)) = (int)htonl(f->cur_index->write_info.data_index_start); /* data_index_start); */
    *((int*)(index_base+12)) = (int)htonl(f->cur_index->write_info.data_index_end); /* data_index_end); */
    ret = f->write_func(f->file_id, index_base, size, NULL, NULL);
    if (ret != size) {
	printf("Index write failed errno %d\n", errno);
    }
    ffs_file_lseek_func(f->file_id, end, SEEK_SET);
    init_write_index_block(f);
}

static void
prepare_index_item(FFSFile f, int item_len)
{
    int next_item_end_offset;
    if (!f->cur_index) {
	init_write_index_block(f);
    }
    next_item_end_offset = f->cur_index->write_info.next_item_offset + item_len;
    if (next_item_end_offset >= (f->cur_index->write_info.index_block_size - 8)) {
	dump_index_block(f);
	f->cur_index->write_info.data_index_start = f->data_count;
    }
}

typedef enum {
    Format_Item = 0, Data_Item = 1, End_Item = 0xff
} IndexItems;

/*
 *  Format item in Index is:
 *  int4   Format_Item + id_size << 8;  (Top 2 bytes reserved, network byte order)
 *  int8   fpos
 *  char *format_id;
 */

static void
output_index_end(FFSFile f)
{
    unsigned char *item_base;
    if (f->file_org != Indexed) return;
    item_base = f->cur_index->write_info.index_block;
    item_base += f->cur_index->write_info.next_item_offset;
    
    *(unsigned int *) item_base = htonl(End_Item);
}

static void
output_format_index(FFSFile f, char *server_id, int id_len)
{
    unsigned char *item_base;
    if (f->file_org != Indexed) return;
    prepare_index_item(f, 12 + id_len);
    item_base = f->cur_index->write_info.index_block;
    item_base += f->cur_index->write_info.next_item_offset;
    
    *(unsigned int *) item_base= htonl(Format_Item + (id_len << 8));
    f->cur_index->write_info.next_item_offset += 4;
    item_base += 4;

    /* upper 32 bits */
#if SIZEOF_OFF_T == 4
    *(unsigned int *) item_base= 0;
#else
    *(unsigned int *) item_base= htonl((f->fpos >> 32));
#endif
    f->cur_index->write_info.next_item_offset += 4;
    item_base += 4;

    /* lower 32 bits */
    *(unsigned int *) item_base= htonl(f->fpos & 0xffffffff);
    f->cur_index->write_info.next_item_offset += 4;
    item_base += 4;

    memcpy(item_base, server_id, id_len);
    f->cur_index->write_info.next_item_offset += id_len;
    item_base += (id_len + 3) & -4;
}

/*
 *  Format item in Index is:
 *  int4   Data_Item + id_size << 8;  (Top 2 bytes reserved, network byte order)
 *  int8   fpos
 *  char   format_id;
 */

static void
output_data_index(FFSFile f, char *server_id, int id_len, char *attr_list,
    int attr_len)
{
    unsigned char *item_base;
    int id_len_nl = htonl(id_len);
    int total_len = (!attr_list)?0:attr_len;
    total_len += id_len + 4;

    if (f->file_org != Indexed) return;
    prepare_index_item(f, 12 + total_len);
    item_base = f->cur_index->write_info.index_block;
    item_base += f->cur_index->write_info.next_item_offset;
    
    *(unsigned int *) item_base= htonl(Data_Item + (total_len << 8));
    f->cur_index->write_info.next_item_offset += 4;
    item_base += 4;

    /* upper 32 bits */
#if SIZEOF_OFF_T == 4
    *(unsigned int *) item_base= 0;
#else
    *(unsigned int *) item_base= htonl((f->fpos >> 32));
#endif
    f->cur_index->write_info.next_item_offset += 4;
    item_base += 4;

    /* lower 32 bits */
    *(unsigned int *) item_base= htonl(f->fpos & 0xffffffff);
    f->cur_index->write_info.next_item_offset += 4;
    item_base += 4;

    /*
     * Entry is organized as |id_len|server_id|attribute_list|
     * id_len gives the length of server_id
     * total_len - id_len gives length of attribute list
     * length of id_length is 4 bytes
     */

    memcpy(item_base, &id_len_nl, 4);
    memcpy(item_base + 4, server_id, id_len);
    memcpy(item_base + 4 + id_len, attr_list, total_len - id_len - 4);

    f->cur_index->write_info.next_item_offset += ((total_len+3) & -4);
    item_base += (total_len + 3) & -4; // alligning on boundary of 4
}

FFSIndexItem
parse_index_block(char *index_base)
{
    FFSIndexItem item;
    int item_count = 0, block_size;
    int cur_offset;
    int done = 0;
    item = malloc(sizeof(FFSIndexItemStruct));
    block_size = htonl(*((int*)(index_base+4))) & 0xffffff;
    (void) block_size; /* avoid warning */
    item->next_index_offset = htonl(*((int*)(index_base+4)));
    item->start_data_count = htonl(*((int*)(index_base+8)));
    item->last_data_count = htonl(*((int*)(index_base+12)));
    item->elements = malloc(sizeof(item->elements[0]));
    cur_offset = 16;
    while (!done) {
	char *elem = (char*)index_base + cur_offset;
	IndexItems item_type = 	(IndexItems) (htonl(*(unsigned int *)elem) & 0xff);
	item_count++;
	item->elements = realloc(item->elements,
				 item_count * sizeof(item->elements[0]));
	memset(&item->elements[item_count-1], 0, sizeof(item->elements[0]));
	switch (item_type) {
	case Format_Item: {
	    unsigned int *ielem = (unsigned int *)((char*) index_base + cur_offset);
	    int id_len = htonl(ielem[0]) >> 8;
#if SIZEOF_OFF_T == 4
	    off_t fpos = 0;
#else
	    off_t fpos = ((off_t)htonl(ielem[1])) << 32;
#endif
	    fpos += htonl(ielem[2]);
	    item->elements[item_count-1].type = FFSformat;
	    item->elements[item_count-1].fpos = fpos;
	    item->elements[item_count-1].format_id = malloc(id_len);
	    item->elements[item_count-1].fid_len = id_len;
	    memcpy(item->elements[item_count-1].format_id,
		   elem + 12, id_len);
	    cur_offset += (id_len + 12 + 3) & -4;
	    break;
	}
	case Data_Item: {
	    unsigned int *ielem = (unsigned int *)((char*) index_base + cur_offset);
	    int total_len = htonl(*(unsigned int *)elem) >> 8;
	    int id_len;
	    int attr_len;

#if SIZEOF_OFF_T == 4
	    off_t fpos = 0;
#else
	    off_t fpos = ((off_t)htonl(ielem[1])) << 32;
#endif
	    fpos += htonl(ielem[2]);
	    id_len = htonl(ielem[3]);
	    attr_len = total_len - id_len - 4;
	    item->elements[item_count-1].type = FFSdata;
	    item->elements[item_count-1].fpos = fpos;
	    item->elements[item_count-1].format_id = malloc(id_len);
	    item->elements[item_count-1].fid_len = id_len;
            if (attr_len <= 1) {
                item->elements[item_count-1].attrs = NULL;
	    } else {
                item->elements[item_count-1].attrs = decode_attr_from_xmit(elem + 16 + id_len);
	    }
	    memcpy(item->elements[item_count-1].format_id,
		   elem + 12, id_len);
	    cur_offset += (id_len + attr_len + 12 + 4) & -4;
	    break;
	}
	case End_Item:
	    item->end_item_offset = cur_offset;
	    done++;
	    break;
	default:
	    printf("Unknown format element\n");
	    /* hope for the best, recovery.... */
	    cur_offset++;
	}
	
    }
    item->elem_count = item_count - 1;
    return item;
}

static
int
write_format_to_file(FFSFile f, FMFormat format)
{
    struct {
	int indicator;
	int format_len;
    } format_header;
    struct iovec vec[4];
    char *server_id;
    int id_len;
    char *server_rep;
    int rep_len;

    server_id = get_server_ID_FMformat(format, &id_len);
    server_rep = get_server_rep_FMformat(format, &rep_len);

    output_format_index(f, server_id, id_len);

    /*
     * next_data indicator is a 2 4-byte chunks in network byte order.
     * In the first chunk, 
     *    the top byte is 0x2, middle 2 are unused and
     *    the bottom byte are the size of the format id;
     * The second chunk holds the length of the format rep;
     */

    format_header.indicator = htonl((id_len & 0xff) | 0x2 << 24);
    format_header.format_len = htonl(rep_len);

    vec[0].iov_len = 8;
    vec[0].iov_base = &format_header;
    vec[1].iov_len = id_len;
    vec[1].iov_base = server_id;
    vec[2].iov_len = rep_len;
    vec[2].iov_base = server_rep;
    vec[3].iov_len = 0;
    vec[3].iov_base = NULL;
    if (f->writev_func(f->file_id, vec, 3, NULL, NULL) != 3) {
	printf("Write failed errno %d\n", errno);
	return 0;
    }
    f->info[format->format_index].written_to_file = 1;
    update_fpos(f);
    return 1;
}

extern int
write_comment_FFSfile(FFSFile f, const char *comment)
{
    struct iovec vec[2];

    size_t byte_size = strlen(comment) + 1;
    int indicator;
    /*
     * next_comment indicator is a 4-byte chunk in network byte order.
     * The top byte is 0x1.  The bottom 3 bytes are the size of the data.
     */
    indicator = htonl((byte_size & 0xffffff) | 0x1 << 24);

    vec[0].iov_len = 4;
    vec[0].iov_base = &indicator;
    vec[1].iov_len = byte_size;
    vec[1].iov_base = (void*)comment;
    if (f->writev_func(f->file_id, vec, 2, NULL, NULL) != 2) {
	printf("Write failed errno %d\n", errno);
	return 0;
    }
    update_fpos(f);
    return 1;
}

extern int
write_encoded_FFSfile(FFSFile f, void *data, DATA_LEN_TYPE byte_size, FFSContext c,
		      attr_list attrs)
{
    FFSTypeHandle h = FFSTypeHandle_from_encode(c, (char*)data);
    FMFormat cf = FMFormat_of_original(h);
    int id_len = 0;
    char *id = get_server_ID_FMformat(cf, &id_len);
    int rep_len = 0;
    char *rep = get_server_rep_FMformat(cf, &rep_len);

    FMFormat f2 = load_external_format_FMcontext(f->fmc, id, id_len, rep);
    int index = f2->format_index;

    struct FFSEncodeVec vec[2];
    int indicator[2];
    int attr_len = 0;
    char *attr_block;
    AttrBuffer b = NULL;
    
    if (attrs) {
	b = create_AttrBuffer();
	attr_block = encode_attr_for_xmit(attrs, b, &attr_len);
    }

    attr_len = (f->file_org == Indexed)?attr_len:0;

    init_format_info(f, index);
    if (!f->info[index].written_to_file) {
	if (write_format_to_file(f, f2) != 1) return 0;
    }

    output_data_index(f, id, id_len, attr_block, attr_len);
    /*
     * next_data indicator is two 4-byte chunks in network byte order.
     * The top byte is 0x3.  The next byte is reserved for future use.
     * The following 6-bytes are the size of the data -- assume size fits
     * in 6 bytes.
     */
    indicator[0] = htonl((0x3 << 24) + (byte_size >> 32));
    indicator[1] = htonl(byte_size & 0xffffffff); 

    vec[0].iov_len = 8;
    vec[0].iov_base = indicator;
    vec[1].iov_len = byte_size;
    vec[1].iov_base = data;
    if (f->writev_func(f->file_id, (struct iovec *)vec, 2, 
		       NULL, NULL) != 2) {
	printf("Write failed, errno %d\n", errno);
	return 0;
    }
    if (b) free_AttrBuffer(b);
    f->data_count++;
    update_fpos(f);
    return 1;
}

				   
extern int
write_FFSfile(FFSFile f, FMFormat format, void *data)
{
    return write_FFSfile_attrs(f, format, data, NULL);
}

extern int
write_FFSfile_attrs(FFSFile f, FMFormat format, void *data, attr_list attrs)
{
    DATA_LEN_TYPE byte_size;
    int index = format->format_index;
    int vec_count;
    FFSEncodeVector vec;
    int indicator[2];
    int id_len = 0;
    char *id = get_server_ID_FMformat(format, &id_len);
    int attr_len = 0;
    char *attr_block = NULL;
    AttrBuffer b = NULL;

    if (attrs) {
	b = create_AttrBuffer();
	attr_block = encode_attr_for_xmit(attrs, b, &attr_len);
    }

    attr_len = (f->file_org == Indexed)?attr_len:0;

    init_format_info(f, index);
    if (!f->info[index].written_to_file) {
	if (write_format_to_file(f, format) != 1) return 0;
    }

    output_data_index(f, id, id_len, attr_block, attr_len);
    //build_reverse_index(f, attr_list, attr_len);
    // Above ... deepak u need to replace attr_list with attr_block
    vec = FFSencode_vector(f->buf, format, data);
    vec_count = 0;
    byte_size = 0;
    while (vec[vec_count].iov_base != NULL) {
	byte_size += vec[vec_count].iov_len;
	vec_count++;
    }

    /*
     * next_data indicator is two 4-byte chunks in network byte order.
     * The top byte is 0x3.  The next byte is reserved for future use.
     * The following 6-bytes are the size of the data -- assume size fits
     * in a signed int.
     */
    indicator[0] = htonl((0x3 << 24) + (byte_size >> 32));
    indicator[1] = htonl(byte_size & 0xffffffff); 

    /* 
     *  utilize the fact that we have some blank vec entries *before the 
     *  beginning* of the vec we are given.  see ffs.c 
     */
    vec--;
    vec_count++;
    vec[0].iov_len = 8;
    vec[0].iov_base = indicator;
    while (vec_count > f->max_iov) {
	/* 
	 * if iovcnt is more than the number of chunks we can write in a 
	 * single AtomicWriteV, loop to write out max allowed.
	 */
	if (f->writev_func(f->file_id, (struct iovec *)vec, f->max_iov, 
			   NULL, NULL) != f->max_iov) {
	    printf("Write failed, errno %d\n", errno);
	    return 0;
	}
	vec_count -= f->max_iov;
	vec += f->max_iov;
    }
    if (f->writev_func(f->file_id, (struct iovec *)vec, vec_count, 
		       NULL, NULL) != vec_count) {
	printf("Write failed, errno %d\n", errno);
	return 0;
    }
    if (b) free_AttrBuffer(b);

    f->data_count++;
    update_fpos(f);
    return 1;
}

extern FFSTypeHandle
FFSread_format(FFSFile ffsfile)
{
    char *id;
    char *rep;
    FMFormat format;
    FFSTypeHandle handle;
    if (ffsfile->read_ahead == FALSE) {
	(void) next_record_type(ffsfile);
    }
    while (ffsfile->next_record_type != FFSformat) {
	switch (ffsfile->next_record_type) {
	case FFScomment:
	    if (ffsfile->tmp_buffer == NULL) {
		ffsfile->tmp_buffer = create_FFSBuffer();
	    }
	    (void) FFSread_comment(ffsfile);
	    (void) next_record_type(ffsfile);
	    break;
	case FFSdata:
	    if (ffsfile->tmp_buffer == NULL) {
		ffsfile->tmp_buffer = create_FFSBuffer();
	    }
	    (void) FFSread(ffsfile, NULL);
	    (void) next_record_type(ffsfile);
	    break;
	default:
	    return NULL;
	}
    }
    
    id = malloc(ffsfile->next_fid_len);
    rep = malloc(ffsfile->next_data_len);
    if (ffsfile->read_func(ffsfile->file_id, id, 
			   ffsfile->next_fid_len, NULL, NULL)
	!= ffsfile->next_fid_len) {
	printf("Read failed, errno %d\n", errno);
	return NULL;
    }
    if (ffsfile->read_func(ffsfile->file_id, rep, 
			   ffsfile->next_data_len, NULL, NULL)
	!= ffsfile->next_data_len) {
	printf("Read failed, errno %d\n", errno);
	return NULL;
    }
    ffsfile->read_ahead = FALSE;
    format = load_external_format_FMcontext(ffsfile->c->fmc, id, 
					    ffsfile->next_fid_len, rep);
    free(id);
    handle = FFSTypeHandle_by_index(ffsfile->c, format->format_index);

    /* in case we're in append mode, mark this as ready in the file */
    init_format_info(ffsfile, format->format_index);
    ffsfile->info[format->format_index].written_to_file = 1;
    return handle;
}

static int
FFSconsume_next_item(FFSFile ffsfile)
{
    switch(next_record_type(ffsfile)) {
    case FFScomment:
	if (ffsfile->tmp_buffer == NULL) {
	    ffsfile->tmp_buffer = create_FFSBuffer();
	}
	(void) FFSread_comment(ffsfile);
	(void) next_record_type(ffsfile);
	return 1;
    case FFSdata:
	if (ffsfile->tmp_buffer == NULL) {
	    ffsfile->tmp_buffer = create_FFSBuffer();
	}
	(void) FFSread(ffsfile, NULL);
	(void) next_record_type(ffsfile);
	return 1;
    case FFSformat:
	(void) FFSread_format(ffsfile);
	(void) next_record_type(ffsfile);
	return 1;
    case FFSindex:
	(void) FFSread_index(ffsfile);
	(void) next_record_type(ffsfile);
	return 1;
    default:
	return 0;
    }
}

extern void
FFSdump_index(FFSIndexItem index_item);

extern FFSIndexItem
FFSread_index(FFSFile ffsfile)
{
    char *index_data;
    FFSIndexItem index_item;
    off_t index_fpos;
    size_t index_size;
    int currentPos = ffs_file_lseek_func(ffsfile->file_id, (size_t)0, SEEK_CUR);
    int end = ffs_file_lseek_func(ffsfile->file_id, (size_t)0, SEEK_END);
    ffs_file_lseek_func(ffsfile->file_id, currentPos, SEEK_SET);   // seek back to original spot

    if (ffsfile->read_ahead == FALSE) {
	(void) next_record_type(ffsfile);
    }
    while (ffsfile->next_record_type != FFSindex) {
	if (!FFSconsume_next_item(ffsfile)) return NULL;
    }
    index_data = malloc(ffsfile->next_data_len);
    index_size = ffsfile->next_data_len;
    update_fpos(ffsfile);
    index_fpos = ffsfile->fpos - 4;
    if (ffsfile->read_func(ffsfile->file_id, index_data+4, 
			   ffsfile->next_data_len-4, NULL, NULL)
	!= ffsfile->next_data_len-4) {
	printf("Read failed, errno %d\n", errno);
	return NULL;
    }
    ffsfile->read_ahead = FALSE;
    index_item = parse_index_block(index_data);
    ffsfile->read_index = index_item;
    if (index_item->next_index_offset == end) {
	ffsfile->cur_index = malloc(sizeof(*(ffsfile->cur_index)));
	memcpy(ffsfile->cur_index, index_item, sizeof(*(ffsfile->cur_index)));
	ffsfile->cur_index->write_info.base_file_pos = index_fpos;
	ffsfile->cur_index->write_info.data_index_start = 0;
	ffsfile->cur_index->write_info.index_block_size = INDEX_BLOCK_SIZE;
	ffsfile->cur_index->write_info.data_index_end = 0;
	ffsfile->cur_index->write_info.next_item_offset = index_item->end_item_offset;
	ffsfile->cur_index->write_info.index_block = malloc(INDEX_BLOCK_SIZE);
	memset(	ffsfile->cur_index->write_info.index_block, 0, INDEX_BLOCK_SIZE);
	memcpy(ffsfile->cur_index->write_info.index_block, index_data, index_size);
	ffsfile->data_count = index_item->last_data_count+1;
    }
    free(index_data);
    index_item->this_index_fpos = index_fpos;
    if (ffsfile->index_head == NULL) {
	ffsfile->index_head = index_item;
	ffsfile->index_tail = index_item;
	index_item->prev = NULL;
	index_item->next = NULL;
    } else {
	ffsfile->index_tail->next = index_item;
	index_item->prev = ffsfile->index_tail;
	index_item->next = NULL;
	ffsfile->index_tail = index_item;
    }
    return index_item;
}

extern void
FFSdump_index(FFSIndexItem index_item)
{
    int i;
    printf("Index item : \n");
    printf(" Next index offset : %ld\n", index_item->next_index_offset);
    printf("  Start data count : %d\n", index_item->start_data_count);
    printf("  End data count   : %d\n", index_item->last_data_count);
    for (i=0; i < index_item->elem_count; i++) {
	FFSIndexElement *elem = &index_item->elements[i];
	int j;
	switch(elem->type) {
	case FFSdata:
	    printf("   Data item at fpos %ld, format ", elem->fpos);
	    for(j=0; j<elem->fid_len; j++) {
		printf("%02x", (unsigned char)(elem->format_id[j]));
	    }
	    if (elem->attrs) {
		printf(", attribute list:");
		dump_attr_list(elem->attrs);
	    } else {
		printf("\n");
	    }
	    break;
	case FFSformat:
	    printf("   Format item at fpos %ld, format ", elem->fpos);
	    for(j=0; j<elem->fid_len; j++) {
		printf("%02x", (unsigned char)(elem->format_id[j]));
	    }
	    printf("\n");
	    break;
	default:
	    break;
	}
    }
}

extern
FFSTypeHandle
FFSnext_type_handle(FFSFile ffsfile)
{
    if (ffsfile->status != OpenForRead)
	return NULL;

    if (ffsfile->read_ahead == FALSE) {
	(void) next_record_type(ffsfile);
    }
    while (ffsfile->next_record_type != FFSdata) {
	if (!FFSconsume_next_item(ffsfile)) return NULL;
    }
    if (ffsfile->next_actual_handle->status == not_checked) {
	FFS_determine_conversion(ffsfile->c, ffsfile->next_actual_handle);
	ffsfile->next_data_handle = ffsfile->next_actual_handle->conversion_target;
    }
    return ffsfile->next_data_handle;
}

extern size_t
FFSfile_next_decode_length(FFSFile iofile)
{
    FFSContext context = iofile->c;
    FFSTypeHandle th = FFSnext_type_handle(iofile);
    DATA_LEN_TYPE len = iofile->next_data_len;
    th = iofile->next_actual_handle;
    return FFS_decode_length_format(context, th, len);
}

extern
char *
FFSread_comment(FFSFile ffsfile)
{
    if (ffsfile->status != OpenForRead)
	return NULL;

    if (ffsfile->read_ahead == FALSE) {
	(void) next_record_type(ffsfile);
    }
    while (ffsfile->next_record_type != FFScomment) {
	if (!FFSconsume_next_item(ffsfile)) return NULL;
    }
    if (ffsfile->tmp_buffer == NULL) ffsfile->tmp_buffer = create_FFSBuffer();
    make_tmp_buffer(ffsfile->tmp_buffer, ffsfile->next_data_len);
    if (ffsfile->read_func(ffsfile->file_id, ffsfile->tmp_buffer->tmp_buffer, 
			   ffsfile->next_data_len, NULL, NULL)
	!= ffsfile->next_data_len) {
	printf("Read failed, errno %d\n", errno);
	return NULL;
    }
    ffsfile->read_ahead = FALSE;
    return ffsfile->tmp_buffer->tmp_buffer;
}

static int
FFSset_fpos(FFSFile file,  off_t fpos)
{
    /* dangerous to set FPOS if not indexed, but we'll allow it */
    if (file->file_org == Indexed) {
	/* 
	 *  We count on index entries on read-side files to be accurate,
	 *  so we don't want to set our FPOS forward past the point where 
	 *  we've read index entries.  If we're trying, read more indexes 
	 *  until we have indexes that cover the desired fpos.
	 */
	int last_element_index;
	struct _FFSIndexElement *last_element;

	last_element_index = file->index_tail->elem_count -1;
	last_element = &file->index_tail->elements[last_element_index];
	while(fpos > last_element->fpos) {
	    /* don't skip forward past index blocks without reading them */
	    if (ffs_file_lseek_func(file->file_id, file->index_tail->next_index_offset, SEEK_SET) == -1)
		return 0;
	    file->read_ahead = FALSE;
	    (void) FFSread_index(file);
	    /* this should update file->index_tail, so continue while loop */
	    last_element_index = file->index_tail->elem_count -1;
	    last_element = &file->index_tail->elements[last_element_index];
	}
    }
    if (ffs_file_lseek_func(file->file_id, fpos, SEEK_SET) == -1) return 0;
    file->read_ahead = FALSE;
    return 1;
}

/*
 * Return data_item if successful else returns 0
 * This routine allows user to seek to kth data block
 */
extern int
FFSseek(FFSFile file, int data_item)
{
    struct _FFSIndexItem *index;
    off_t fpos;
    int index_item;
    FFSIndexItem prev_index_tail = NULL;

    if (data_item < 0)
        /* Or should it be set to 0  */
        return 0;

    if (!file->index_head)
        FFSread_index(file);
    if (!file->index_head)
        // If file is not indexed then this condition will be activated.
        return 0;

    /* seek to the N'th data item in the file */
    while (data_item > file->index_tail->last_data_count &&
           file->index_tail != prev_index_tail) {
	/* don't skip forward past index blocks without reading them */
	if (ffs_file_lseek_func(file->file_id, file->index_tail->next_index_offset, SEEK_SET) == -1)
	    return 0;
    	file->read_ahead = FALSE;
        prev_index_tail = file->index_tail;
	(void) FFSread_index(file);
    }

    if (file->index_tail->last_data_count < data_item)
        // How can I can more than what has been written! :D
        return 0;

    index = file->index_head;
    while (data_item > index->last_data_count)
        index = index->next;

    data_item -= index->start_data_count;
    data_item++;
    index_item = 0;
    while (data_item > 0) {
	if (index->elements[index_item].type == FFSdata) {
	    data_item--;
	}
	index_item++;
    }
    index_item--;
    fpos = index->elements[index_item].fpos;
    FFSset_fpos(file, fpos);
    file->data_block_no = data_item;

    return data_item;
}
    
static
int
get_AtomicInt(FFSFile file, FILE_INT *file_int_ptr)
{
#if SIZEOF_INT == 4
    int tmp_value;
    if (file->read_func(file->file_id, &tmp_value, 4, NULL, NULL) != 4)
	return 0;

#else
    Baaad shit;
#endif
    *file_int_ptr = tmp_value;
    return 1;
}

static void
read_all_index_and_formats(FFSFile file)
{
    off_t fpos = 1;
    int currentPos = ffs_file_lseek_func(file->file_id, (size_t)0, SEEK_CUR);
    int end = ffs_file_lseek_func(file->file_id, (size_t)0, SEEK_END);
    ffs_file_lseek_func(file->file_id, currentPos, SEEK_SET);   // seek back to original spot

    if (!file->index_head)
        FFSread_index(file);
    if (!file->index_head) {
        /* If file is not indexed then this condition will be activated.
	   The FFSread_index() will have caused the entire file to have 
	   been read and fpos set to the EOF.  Return and continue writing.
	*/
	return;
    }
    file->file_org = Indexed;
    while (fpos != end) {
	/* Read all formats in this index block */
	int i;
	for (i = 0; i < file->index_tail->elem_count; i++) {
	    if (file->index_tail->elements[i].type == FFSformat) {
		fpos = file->index_tail->elements[i].fpos;
		if (ffs_file_lseek_func(file->file_id, fpos, SEEK_SET) == -1)
		    return;
		(void) FFSread_format(file);
	    }
	}
	/* skip to next index block */
	fpos = file->index_tail->next_index_offset;
	if (fpos != end) {
	    if (ffs_file_lseek_func(file->file_id, fpos, SEEK_SET) == -1)
		return;
	    FFSread_index(file);
	}
    }
    ffs_file_lseek_func(file->file_id, 0, SEEK_END);
    file->fpos = ffs_file_lseek_func(file->file_id, 0, SEEK_CUR);
}

static void
convert_last_index_block(FFSFile ffsfile)
{
    FFSIndexItem read_index = ffsfile->index_tail;

    init_write_index_block(ffsfile);

    unsigned char *index_data;

    if (read_index == NULL) return;

    FFSset_fpos(ffsfile,  read_index->this_index_fpos);
    
    (void) next_record_type(ffsfile);
    assert(ffsfile->next_record_type == FFSindex);

    index_data = ffsfile->cur_index->write_info.index_block;
    if (ffsfile->read_func(ffsfile->file_id, index_data+4, 
			   ffsfile->next_data_len-4, NULL, NULL)
	!= ffsfile->next_data_len-4) {
	printf("Read failed, errno %d\n", errno);
	return;
    }
    ffsfile->cur_index->write_info.data_index_start =  htonl(*((int*)(index_data+8)));;
    ffsfile->data_count = read_index->last_data_count + 1;
    if (ffs_file_lseek_func(ffsfile->file_id, 0, SEEK_END) == -1)
	return;
}

FFSRecordType
FFSnext_record_type(FFSFile ffsfile)
{
    FFSRecordType next = next_record_type(ffsfile);
    while ((next & ffsfile->visible_items_bitmap) != next) {
	FFSconsume_next_item(ffsfile);
	next = next_record_type(ffsfile);
    }
    return next;
} 

extern void
FFSset_visible(FFSFile ffsfile, int bitmap)
{
    bitmap |= FFSerror | FFSend;  /* always visible */
    ffsfile->visible_items_bitmap = bitmap;
}

static
FFSRecordType
next_record_type(FFSFile ffsfile)
{
    FILE_INT indicator_chunk = 0;
 restart:
    if (ffsfile->status != OpenForRead) {
	return FFSerror;
    }
    if (ffsfile->tmp_buffer == NULL) {
	ffsfile->tmp_buffer = create_FFSBuffer();
    }
    if (ffsfile->read_ahead == FALSE) {
	if (!get_AtomicInt(ffsfile, &indicator_chunk)) {
	    ffsfile->next_record_type = 
		(ffsfile->errno_val) ? FFSerror : FFSend;
	    return ffsfile->next_record_type;
	}
	
	indicator_chunk = ntohl(indicator_chunk);
	switch (indicator_chunk >> 24) {
	case 0x1: /* comment */
		ffsfile->next_record_type = FFScomment;
		ffsfile->next_data_len = indicator_chunk & 0xffffff;
		break;
	case 0x2: /* format */
		ffsfile->next_record_type = FFSformat;
		ffsfile->next_fid_len = indicator_chunk & 0xffffff;
		if (!get_AtomicInt(ffsfile, &indicator_chunk)) {
		    ffsfile->next_record_type = (ffsfile->errno_val) ? FFSerror : FFSend;
		    return ffsfile->next_record_type;
		}
		ffsfile->next_data_len = ntohl(indicator_chunk);
		break;
	case 0x3: /* data */ {
		char *tmp_buf;
		int header_size;
		DATA_LEN_TYPE top_data_len_bytes;
		ffsfile->next_record_type = FFSdata;
		top_data_len_bytes = (indicator_chunk & 0xffff);
		if (!get_AtomicInt(ffsfile, &indicator_chunk)) {
		    ffsfile->next_record_type = (ffsfile->errno_val) ? FFSerror : FFSend;
		    return ffsfile->next_record_type;
		}
		ffsfile->next_data_len = ntohl(indicator_chunk) + 
		    (top_data_len_bytes << 32);
		make_tmp_buffer(ffsfile->tmp_buffer, ffsfile->next_data_len);
		tmp_buf = ffsfile->tmp_buffer->tmp_buffer;
		/* first get format ID, at least 8 bytes */
		if (ffsfile->read_func(ffsfile->file_id, tmp_buf, 8, NULL,
                                       NULL) != 8) {
		    ffsfile->next_record_type = (ffsfile->errno_val) ?
                        FFSerror : FFSend;
		    return ffsfile->next_record_type;
		}
		ffsfile->next_fid_len = FMformatID_len(tmp_buf);
		if (ffsfile->next_fid_len > 8) {
		    int more = ffsfile->next_fid_len - 8;
		    if (ffsfile->read_func(ffsfile->file_id, tmp_buf + 8,
                                           more, NULL, NULL) != more) {
			ffsfile->next_record_type = (ffsfile->errno_val) ?
                            FFSerror : FFSend;
			return ffsfile->next_record_type;
		    }
		}

		ffsfile->next_actual_handle = 
		    FFSTypeHandle_from_encode(ffsfile->c, tmp_buf);

                if (!ffsfile->next_actual_handle && ffsfile->index_head) {

                    struct _FFSIndexItem *index = NULL;
                    off_t fpos_bak = ffs_file_lseek_func(ffsfile->file_id, 0, SEEK_CUR);
		    int fid_len = ffsfile->next_fid_len;
		    char tmp_fid_storage[64];
		    size_t tmp_data_len;
		    int done = 0;
		    assert(sizeof(tmp_fid_storage) > fid_len);
		    /* store away the format ID we've read */
		    memcpy(tmp_fid_storage, tmp_buf, fid_len);
		    tmp_data_len = ffsfile->next_data_len;

                    index = ffsfile->index_head;
                    while (!done && index) {
			int i;
			for (i=0 ; i < index->elem_count; i++) {
			    struct _FFSIndexElement *elem;
			    elem = &index->elements[i];
			    if (elem->type == FFSformat &&
				elem->fid_len == ffsfile->next_fid_len &&
				!(memcmp(elem->format_id,
					 tmp_buf, ffsfile->next_fid_len))) {

				if (ffs_file_lseek_func(ffsfile->file_id, elem->fpos, SEEK_SET) != -1) {
				    ffsfile->read_ahead = FALSE;
				    FFSread_format(ffsfile);
				    ffs_file_lseek_func(ffsfile->file_id, fpos_bak, SEEK_SET);
				    ffsfile->read_ahead = TRUE;
				    ffsfile->next_record_type = FFSdata;
				    /* tmp_buf might have changed */
				    tmp_buf = ffsfile->tmp_buffer->tmp_buffer;
				    /* put back the format ID we read earlier */
				    memcpy(tmp_buf, tmp_fid_storage, fid_len);
				    ffsfile->next_data_len = tmp_data_len;

				    ffsfile->next_actual_handle = 
					FFSTypeHandle_from_encode(ffsfile->c,
								  tmp_buf);
				    done++;
				    break;
				}
			    }
                        }

                        index = index->next;
                    }
                }
		/* GSE
		 * If ffsfile->next_actual_handle is NULL here, we have 
		 * a problem.  The only way this should happen is if we 
		 * have used FFSset_fpos() to seek forward in a file, 
		 * skipping over the Format record associated with this
		 * data item.  But now we NEED IT.  We must :
		 *  1) Save our current FPOS
		 *  2) find the format in the index
		 *  3) lseek the its location (set read_ahead to false)
		 *  4) read the format
		 *  5) lseek back to our saved FPOS
		 *  6) call FFSTypeHandle_from_encode() again to set
		 *     next_actual_handle correctly.
		 */
		ffsfile->next_data_handle = 
		    FFS_target_from_encode(ffsfile->c, tmp_buf);

		if ((ffsfile->next_data_handle == NULL) &&
                    (!ffsfile->raw_flag)) {
		    /* no target for this format, discard */
		    size_t more = ffsfile->next_data_len - ffsfile->next_fid_len;
		    if (ffsfile->read_func(ffsfile->file_id, tmp_buf +
                                           ffsfile->next_fid_len, more, NULL,
                                           NULL) != more) {
			ffsfile->next_record_type = (ffsfile->errno_val) ?
                            FFSerror : FFSend;
			return ffsfile->next_record_type;
		    }
		    ffsfile->read_ahead = FALSE;
		    goto restart;
		    
		}
		header_size = FFSheader_size(ffsfile->next_actual_handle);
		if (header_size > ffsfile->next_fid_len) {
		    int more = header_size - ffsfile->next_fid_len;
		    if (ffsfile->read_func(ffsfile->file_id, tmp_buf +
                                           ffsfile->next_fid_len, more, NULL,
                                           NULL) != more) {
			ffsfile->next_record_type = (ffsfile->errno_val) ?
                            FFSerror : FFSend;
			return ffsfile->next_record_type;
		    }
		}
                ffsfile->data_block_no++;
		break;
	case 0x4: /* index */
		ffsfile->next_record_type = FFSindex;
		ffsfile->next_data_len = indicator_chunk & 0xffffff;
/*		if (!ffsfile->expose_index) {
		    ffs_file_lseek_func((int)(intptr_t)ffsfile->file_id, INDEX_BLOCK_SIZE-4, SEEK_CUR);
		    ffsfile->read_ahead = FALSE;
		    return next_record_type(ffsfile);
		    }*/
		break;
	default:
	    printf("CORRUPT FFSFILE\n");
	    exit(1);
	}
	}
	ffsfile->read_ahead = TRUE;
    }
    return ffsfile->next_record_type;
}

extern size_t
FFSnext_data_length(FFSFile file)
{
    if (file->status != OpenForRead)
	return 0;

    if (file->read_ahead == FALSE) {
	(void) next_record_type(file);
    }
    while (file->next_record_type != FFSdata) {
	if (!FFSconsume_next_item(file)) return 0;
    }
    return file->next_data_len;
}

extern int
FFSread(FFSFile file, void *dest)
{
    int header_size;
    size_t read_size;
    char *tmp_buf;

    if (file->status != OpenForRead)
	return 0;

    if (file->read_ahead == FALSE) {
	(void) next_record_type(file);
    }
    while (file->next_record_type != FFSdata) {
	if (!FFSconsume_next_item(file)) return 0;
    }

    header_size = FFSheader_size(file->next_actual_handle);
    read_size = file->next_data_len - header_size;
    tmp_buf = file->tmp_buffer->tmp_buffer;
    /* should have buffer optimization logic here.  
     * I.E. if decode_in_place_possible() handle differently.  later
     */
    /* read into temporary memory */
    if (file->read_func(file->file_id, tmp_buf + header_size, read_size, NULL, NULL) != read_size) {
	file->next_record_type = (file->errno_val) ? FFSerror : FFSend;
	return 0;
    }
    if (dest != NULL) {
	FFSdecode(file->c, file->tmp_buffer->tmp_buffer, dest);
    }
    file->read_ahead = FALSE;
    return 1;
}

/*
 * In future can be exposed to user so that user can get attr_list
 * of kth data block.
 * Currently this function does not read the unread index blocks to satify
 * the request when index item of data_item is not read.
 *
 * CAUTION:
 * attr_list belongs to index and should not be manipulated/freed by user.
 */
static attr_list
FFSread_data_attr(FFSFile file, int data_item)
{
    struct _FFSIndexItem *index;
    int index_item;
    attr_list atl = NULL;

    if (data_item < 0 || !file->index_head ||
        file->index_tail->last_data_count < data_item)
        goto exit;

    index = file->index_head;
    while (data_item > index->last_data_count)
        index = index->next;

    data_item -= index->start_data_count;
    data_item++;
    index_item = 0;
    while (data_item > 0) {
	if (index->elements[index_item].type == FFSdata) {
	    data_item--;
	}
	index_item++;
    }
    index_item--;

    atl = index->elements[index_item].attrs;

exit:
    return atl;
}

extern attr_list
FFSattrs_from_last_read(FFSFile file)
{
    return FFSread_data_attr(file, file->data_block_no-1);
}

extern int
FFSread_attr(FFSFile file, void *dest, attr_list *attr)
{
    int ret = FFSread(file, dest);
    attr_list atl = FFSread_data_attr(file, file->data_block_no-1);
    if (attr) *attr = atl;

    return ret;
}

extern int
FFSread_raw(FFSFile file, void *dest, int buffer_size, FFSTypeHandle *fp)
{
    FFSTypeHandle f;
    int header_size;
    size_t read_size;

    if (file->status != OpenForRead)
	return 0;

    if (file->read_ahead == FALSE) {
	(void) next_record_type(file);
    }
    while (file->next_record_type != FFSdata) {
	if (!FFSconsume_next_item(file)) return 0;
    }

    f = file->next_actual_handle;
    *fp = f;
    header_size = FFSheader_size(f);
    read_size = file->next_data_len - header_size;

    if (file->read_func(file->file_id, dest, read_size, NULL, NULL) != read_size) {
	file->next_record_type = (file->errno_val) ? FFSerror : FFSend;
	return 0;
    }

    file->read_ahead = FALSE;
    return 1;
}

extern int
FFSread_raw_header(FFSFile file, void *dest, int buffer_size, FFSTypeHandle *fp)
{
    FFSTypeHandle f;
    int header_size;
    size_t read_size;

    if (file->status != OpenForRead)
	return 0;

    if (file->read_ahead == FALSE) {
	(void) next_record_type(file);
    }
    while (file->next_record_type != FFSdata) {
	if (!FFSconsume_next_item(file)) return 0;
    }

    f = file->next_actual_handle;
    *fp = f;
    header_size = FFSheader_size(f);
    read_size = file->next_data_len - header_size;
    /* should have buffer optimization logic here.  
     * I.E. if decode_in_place_possible() handle differently.  later
     */
    /* read into temporary memory */
    memset(dest, 0, header_size);
    memcpy(dest, f->body->server_ID.value, f->body->server_ID.length);
    if (file->read_func(file->file_id, (char*)dest+header_size, read_size, NULL, NULL) != read_size) {
	file->next_record_type = (file->errno_val) ? FFSerror : FFSend;
	return 0;
    }

    file->read_ahead = FALSE;
    return 1;
}

extern int
FFSread_to_buffer(FFSFile file, FFSBuffer b,  void **dest)
{
    int header_size;
    size_t read_size;
    char *tmp_buf;

    if (file->status != OpenForRead)
	return 0;

    if (file->read_ahead == FALSE) {
	(void) next_record_type(file);
    }
    while (file->next_record_type != FFSdata) {
	if (!FFSconsume_next_item(file)) return 0;
    }

    header_size = FFSheader_size(file->next_actual_handle);
    read_size = file->next_data_len - header_size;
    tmp_buf = file->tmp_buffer->tmp_buffer;
    /* should have buffer optimization logic here.  
     * I.E. if decode_in_place_possible() handle differently.  later
     */
    /* read into temporary memory */
    if (file->read_func(file->file_id, tmp_buf + header_size, read_size, NULL, NULL) != read_size) {
	file->next_record_type = (file->errno_val) ? FFSerror : FFSend;
	return 0;
    }
    FFSdecode_to_buffer(file->c, file->tmp_buffer->tmp_buffer, b->tmp_buffer);
    file->read_ahead = FALSE;
    if (dest) *dest = b->tmp_buffer;
    return 1;
}

