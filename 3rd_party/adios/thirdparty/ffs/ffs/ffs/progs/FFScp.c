#include "config.h"
#include <fcntl.h>
#include <string.h>
#include <assert.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "ffs.h"

char *usage = "\
Usage: FFScp [<option>] <in filename> <out filename>\n\
       FFScp copies (rewrites) FFS files. \n";

extern char *
dump_raw_FMrecord_to_string(FMContext fmc, FMFormat format, void *data);

extern int
FFSread_raw_header(FFSFile file, void *dest, int buffer_size, FFSTypeHandle *fp);
extern int
write_encoded_FFSfile(FFSFile f, void *data, DATA_LEN_TYPE byte_size, FFSContext c,
		      attr_list attrs);

int
main(int argc, char **argv)
{
    FFSFile in_file = NULL, out_file = NULL;
    int buffer_size = 1024;
    char *buffer = NULL;
    int i;
    int bitmap, format_count;
    FFSTypeHandle *in_formats;
    FMFormat *out_formats;
    char *in_filename = NULL, *out_filename = NULL;
    int rewrite = 0;

    for (i = 1; i < argc; i++) {
	if (argv[i][0] == '-') {
	    if (strcmp(argv[i], "-rewrite") == 0) {
		rewrite++;
	    } else {
		printf("Unknown argument \"%s\"\n", argv[i]);
	    }
	} else {
	    if (in_filename == NULL) {
		in_filename = argv[i];
	    } else if (out_filename == NULL) {
		out_filename = argv[i];
	    } else {
		fprintf(stderr, "Extra argument specified \"%s\"\n%s\n",
			argv[i], usage);
		exit(1);
	    }
	}
    }
    if (!in_filename || !out_filename) {
	fprintf(stderr, "%s", usage);
	exit(1);
    }
    if (rewrite) {
	in_file = open_FFSfile(in_filename, "r");
    } else {
	in_file = open_FFSfile(in_filename, "R");
    }
    out_file = open_FFSfile(out_filename, "w");

    if (in_file == NULL) {
	printf("File Open Failure \"%s\"", in_filename);
	perror("Opening input file");
	exit(1);
    }

    if (out_file == NULL) {
	printf("File Open Failure \"%s\"", in_filename);
	perror("Opening output file");
	exit(1);
    }

    bitmap = FFSdata;
    if (rewrite) bitmap |= FFSformat;
    FFSset_visible(in_file, bitmap);

    buffer = malloc(1024);
    buffer_size = 1024;
    format_count = 0;
    in_formats = malloc(sizeof(FMFormat));
    out_formats = malloc(sizeof(FMFormat));
    while (1) {
        switch (FFSnext_record_type(in_file)) {
        case FFSformat:
	  {
	    /* this only happens when we are rewriting the file */
            FFSTypeHandle format = FFSread_format(in_file);
	    FMFormat original = FMFormat_of_original(format);
	    FMStructDescList local = get_localized_formats(original);
	    FFSContext in_context = FFSContext_of_file(in_file);
	    FFSContext out_context = FFSContext_of_file(out_file);
	    in_formats[format_count] = FFSset_fixed_target(in_context, local);
	    out_formats[format_count] = register_data_format(FMContext_from_FFS(out_context), local);
	    format_count++;
	    in_formats = realloc(in_formats, sizeof(FMFormat)*(format_count+1));
	    out_formats = realloc(out_formats, sizeof(FMFormat)*(format_count+1));
	    in_formats[format_count] = NULL;
            break;
        }
        case FFSdata:{
            if (buffer_size < FFSnext_data_length(in_file)) {
                buffer_size = FFSnext_data_length(in_file);
                buffer = realloc(buffer, buffer_size);
            }
	    if (rewrite) {
		FFSTypeHandle format = FFSnext_type_handle(in_file);
		FMFormat out_format = NULL;
		int i;
		for (i=0 ; i < format_count; i++) {
		    if (in_formats[i] == format) {
			out_format = out_formats[i];
		    }
		}
		assert(out_format);
		FFSread(in_file, buffer);
		write_FFSfile(out_file, out_format, buffer);
	    } else {
		FFSTypeHandle format;
	        FFSread_raw_header(in_file, buffer, buffer_size, &format);
		write_encoded_FFSfile(out_file, buffer, buffer_size, 
				      FFSContext_of_file(in_file), NULL);
	    }
	    break;
        }
        case FFSerror:
        case FFSend:
            close_FFSfile(in_file);
            close_FFSfile(out_file);
            exit(0);
	default:
	    break;
        }
    }
}
