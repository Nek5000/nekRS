#include "config.h"
#include <fcntl.h>
#include <string.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "ffs.h"
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#endif

extern int FMdumpVerbose;
extern void FFSprint_version();
extern void FFSdump_index(FFSIndexItem index_item);
extern void dump_FMFormat(FMFormat ioformat);
extern int FFSread_raw(FFSFile file, void *dest, int buffer_size, FFSTypeHandle *fp);

char *usage = "\
Usage: FFSdump [<option>] <filename>\n\
       FFSdump dumps records from FFS files.  Records may be\n\
       file headers, comments, record formats and data.\n\
       By default, record formats are not printed.  Printing for\n\
       any record type can be turned on or off with an option.\n\
          <option> is of the form {+,-}{<record type>}\n\
       where <record type> is one of : header, comments, index, formats, data\n";

extern char *
dump_raw_FMrecord_to_string(FMContext fmc, FMFormat format, void *data);

int
main(int argc, char **argv)
{
    FFSFile ffsfile = NULL;
    int buffer_size = 1024;
    char *buffer = NULL;
    int i;
    char *filename = NULL;
    int dump_formats = 0, dump_index = 0, dump_header = 1, dump_comments = 1, dump_data = 1;
    int bitmap;
#ifdef HAVE_WINDOWS_H
    /* set up a scrolling console screen */
    HANDLE hConsole=GetStdHandle(STD_OUTPUT_HANDLE);
    COORD xy;
    SMALL_RECT rect;

    xy.X=100;
    xy.Y=300;
    rect.Left=0;
    rect.Top=0;
    rect.Right=90;
    rect.Bottom=50;

    SetConsoleScreenBufferSize(hConsole,xy);
    SetConsoleWindowInfo(hConsole,TRUE,&rect);
#endif

    printf("Dumping program is ");
    FFSprint_version();
    for (i = 1; i < argc; i++) {
	if ((argv[i][0] == '-') || (argv[i][0] == '+')) {
	    int value = (argv[i][0] == '-') ? 0 : 1;
	    if (strcmp(&argv[i][1], "formats") == 0) {
		dump_formats = value;
	    } else if (strcmp(&argv[i][1], "index") == 0) {
		dump_index = value;
	    } else if (strcmp(&argv[i][1], "header") == 0) {
		dump_header = value;
	    } else if (strcmp(&argv[i][1], "comments") == 0) {
		dump_comments = value;
	    } else if (strcmp(&argv[i][1], "data") == 0) {
		dump_data = value;
	    } else {
		fprintf(stderr, "Unknown option \"%s\"\n%s\n", argv[i], usage);
		exit(1);
	    }
	} else {
	    if (filename == NULL) {
		filename = argv[i];
	    } else {
		fprintf(stderr, "Extra argument specified \"%s\"\n%s\n",
			argv[i], usage);
		exit(1);
	    }
	}
    }

    if (filename == NULL) {
	fprintf(stderr, "%s", usage);
	exit(1);
    }
    ffsfile = open_FFSfile(filename, "R");

    if (ffsfile == NULL) {
	printf("File Open Failure \"%s\"", filename);
	perror("Opening input file");
	exit(1);
    }
    FMdumpVerbose = 1;
    bitmap = 0;
    bitmap |= (dump_formats?FFSformat:0);
    bitmap |= (dump_index?FFSindex:0);
    bitmap |= (dump_data?FFSdata:0);
    bitmap |= (dump_comments?FFScomment:0);

    FFSset_visible(ffsfile, bitmap);

    if (dump_header)
/*	dump_FFSFile(ffsfile);*/

    buffer = malloc(1024);
    buffer_size = 1024;

    while (1) {
        switch (FFSnext_record_type(ffsfile)) {
        case FFSindex: {
            FFSIndexItem item = FFSread_index(ffsfile);
            FFSdump_index(item);
            break;
        }
        case FFScomment: {
            char *comment = FFSread_comment(ffsfile);
            printf("Got a comment\n");
            if (dump_comments) {
                if (comment)
                    printf("Comment -> \"%s\"\n", comment);
            }
            break;
        }
        case FFSformat:
        {
            FFSTypeHandle format = FFSread_format(ffsfile);
            dump_formats = 1;
            if (dump_formats) {
                dump_FMFormat(FMFormat_of_original(format));
            }
            break;
        }
        case FFSdata:{
            FFSTypeHandle format;
            if (buffer_size < FFSnext_data_length(ffsfile)) {
                buffer_size = FFSnext_data_length(ffsfile);
                buffer = realloc(buffer, buffer_size);
            }
            FFSread_raw(ffsfile, buffer, buffer_size, &format);
            if (dump_data) {
/*              dump_raw_FMrecord(ffsfile, FMFormat_of_original(format), buffer);*/
                char *str = dump_raw_FMrecord_to_string(NULL, FMFormat_of_original(format), buffer);
                printf("%s", str);
                free(str);
            }
            break;
        }
        case FFSerror:
        case FFSend:
            close_FFSfile(ffsfile);
            exit(0);
        }
    }
}
