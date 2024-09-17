#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "drisc.h"
#include "ecl.h"

/* image size */
#define XSZ (640)
#define YSZ (320)

/* input/output image buffers: raw RGB color data */
unsigned char inimg [3 * YSZ * XSZ];
unsigned char outimg[3 * YSZ * XSZ];

/* in: half-white/half-black, out: pure black */
void init()
{
    memset((void *)inimg,                     0, XSZ*YSZ*3);
    memset((void *)(inimg + XSZ*(YSZ/2)*3), 255, XSZ*(YSZ/2)*3);
    memset((void *)outimg,                    0, XSZ*YSZ*3);
}


/* macro for 2D index expansion & access to color component */
#define pix(ptr,j,i,y,x,c)     (*((ptr)+3*((j)*(x)+(i))+(c)))

void edge(unsigned char * pin, unsigned char * pout, int y, int x)
{
    int i;
    int j;
    
    int jp, jn;
    int ip, in;
    
    for(i=0; i<x; i++) {
	for(j=0; j<y; j++) {

	    ip = i - 1; in = i + 1;
	    jp = j - 1; jn = j + 1;

	    /* Clamp (i,j) at (0,0) and (x-1,y-1) */
	    if(i-1 <  0) ip = 0; if(i+1 >= x) in = x-1;
	    if(j-1 <  0) jp = 0; if(j+1 >= y) jn = y-1;

	    /* red */
	    pix(pout,j,i,y,x,0) =
		- pix(pin,jp,ip,y,x,0) -     pix(pin,jp,i,y,x,0) - pix(pin,jp,in,y,x,0)
		- pix(pin,j ,ip,y,x,0) + 8 * pix(pin,j, i,y,x,0) - pix(pin,j, in,y,x,0)
		- pix(pin,jn,ip,y,x,0) -     pix(pin,jn,i,y,x,0) - pix(pin,jn,in,y,x,0);

	    /* green */
	    pix(pout,j,i,y,x,1) =
		- pix(pin,jp,ip,y,x,1) -     pix(pin,jp,i,y,x,1) - pix(pin,jp,in,y,x,1)
		- pix(pin,j ,ip,y,x,1) + 8 * pix(pin,j, i,y,x,1) - pix(pin,j, in,y,x,1)
		- pix(pin,jn,ip,y,x,1) -     pix(pin,jn,i,y,x,1) - pix(pin,jn,in,y,x,1);


	    /* blue */
	    pix(pout,j,i,y,x,2) =
		- pix(pin,jp,ip,y,x,2) -     pix(pin,jp,i,y,x,2) - pix(pin,jp,in,y,x,2)
		- pix(pin,j ,ip,y,x,2) + 8 * pix(pin,j, i,y,x,2) - pix(pin,j, in,y,x,2)
		- pix(pin,jn,ip,y,x,2) -     pix(pin,jn,i,y,x,2) - pix(pin,jn,in,y,x,2);

	}
    }
    
}

/* same code as edge() above, only w/ macros expanded */
char edge_proto[] = "void edge(unsigned char * pin, unsigned char * pout, int y, unsigned char x)";
char edge_body [] = "{ \
x =  - (*pin);\
 }";\


/* -------------------------------------------------------- */

/* PNG header for writing to disk */
char header[] = "P6 640 320 255\n";

int main()
{
    ecl_parse_context context = new_ecl_parse_context();
    ecl_code gen_code;
    void (*func)(unsigned char *, unsigned char *, int, int);
    long result;


    ecl_subroutine_declaration(edge_proto, context);
    gen_code = ecl_code_gen(edge_body, context);
    func = (void (*)(unsigned char *, unsigned char *, int, int)) (long) gen_code->func;
 

    init();
    edge(inimg, outimg, YSZ, XSZ);

    {
	int fd;

	fd = creat("gcc.png", O_RDWR | 0644);
	if(fd < 0) { perror("creat"); return -1; }
	
	write(fd, header, strlen(header));
	write(fd, outimg, XSZ*YSZ*3);
	
	close(fd);
    }

    init();
    func(inimg, outimg, YSZ, XSZ);

    {
	int fd;

	fd = creat("ecl.png", O_RDWR | 0644);
	if(fd < 0) { perror("creat"); return -1; }
	
	write(fd, header, strlen(header));
	write(fd, outimg, XSZ*YSZ*3);
	
	close(fd);
    }

    ecl_code_free(gen_code);
    ecl_free_parse_context(context);

    return 0;
}
