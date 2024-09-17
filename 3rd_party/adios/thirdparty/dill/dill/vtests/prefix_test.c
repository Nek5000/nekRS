#include "stdio.h"
#include "dill.h"

static int verbose = 0;

int gg(int a, int b) {
     return a+b;
}
int ff(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j) {
    if (verbose) printf("In ff  a=%d, b=%d, c=%d, d=%d, e=%d, f=%d, g=%d, h=%d, i=%d, j=%d\n",
	   a, b, c, d, e, f, g, h, i, j);
     return a+b+c+d+e+f+g+h+i+j;
}

int a () {
    dill_stream s = dill_create_stream();
     dill_exec_handle handle;
     dill_reg a,b,p3,d,e,f,g,h,i,j,w,z,cnt;
/*     dill_reg func;*/
     int L1;
     int (*ip)();
     int result;

     dill_start_simple_proc(s, "a_gen", DILL_I);
     L1 = dill_alloc_label(s, "before call");
/*     func = dill_getreg(s, DILL_P);*/
     cnt = dill_getreg(s, DILL_P);
     a = dill_getreg(s, DILL_I);
     b = dill_getreg(s, DILL_I);
     p3 = dill_getreg(s, DILL_I);
     d = dill_getreg(s, DILL_I);
     e = dill_getreg(s, DILL_I);
     f = dill_getreg(s, DILL_I);
     g = dill_getreg(s, DILL_I);
     h = dill_getreg(s, DILL_I);
     i = dill_getreg(s, DILL_I);
     j = dill_getreg(s, DILL_I);
     z = dill_getreg(s, DILL_I);
     w = dill_getreg(s, DILL_I);

/*     dill_setp(s, func, ff);*/
     dill_seti(s, h, 8);
     dill_seti(s, i, 9);
     dill_seti(s, j, 0);
     dill_seti(s, z, 0);
     dill_seti(s, cnt, 0);
     dill_mark_label(s, L1);
     dill_push_init(s);
     if (!dill_do_reverse_vararg_push(s)) {
	 dill_push_argi(s, a);
	 dill_push_argi(s, b);
	 dill_push_argi(s, p3);
	 dill_push_argi(s, d);
	 dill_push_argi(s, e);
	 dill_push_argi(s, f);
	 dill_push_argi(s, g);
	 dill_push_argi(s, h);
	 dill_push_argi(s, i);
	 dill_push_argi(s, j);
     } else {
	 dill_push_argi(s, j);
	 dill_push_argi(s, i);
	 dill_push_argi(s, h);
	 dill_push_argi(s, g);
	 dill_push_argi(s, f);
	 dill_push_argi(s, e);
	 dill_push_argi(s, d);
	 dill_push_argi(s, p3);
	 dill_push_argi(s, b);
	 dill_push_argi(s, a);
     }	 
     w = dill_calli(s, (void*)ff, "ff");
     dill_addi(s, z,z,w);
     dill_addii(s, cnt,cnt,1);
     dill_bltii(s, cnt,2,L1);
     dill_addi(s, z, z, a);
     dill_addi(s, z, z, b);
     dill_addi(s, z, z, p3);
     dill_addi(s, z, z, d);
     dill_addi(s, z, z, e);
     dill_addi(s, z, z, f);
     dill_addi(s, z, z, g);
     dill_addi(s, z, z, h);
     dill_addi(s, z, z, i);
     dill_addi(s, z, z, j);
     dill_reti(s, z);			/* (9*10/2)*3 = 135 */
     dill_begin_prefix_code(s);
     dill_seti(s, a, 1);
     dill_seti(s, b, 2);
     dill_seti(s, p3, 3);
     dill_seti(s, d, 4);
     dill_seti(s, e, 5);
     dill_seti(s, f, 6);
     dill_seti(s, g, 7);

     handle = dill_finalize(s);
     ip = (int(*)())dill_get_fp(handle);

     if (verbose) dill_dump(s);

     result = (*ip)();
     dill_free_stream(s);
     dill_free_handle(handle);
     if (result != 135) printf("Result not 135, %d\n",result);
     return result == 135;
}

int c () 
{
     dill_reg a,b,c, f;
     int (*ip)();
     int ret;
     dill_exec_handle h;
     
     dill_stream s = dill_create_stream();

     dill_start_simple_proc(s, "c_gen", DILL_I);
     a = dill_getreg(s, DILL_I);
     c = dill_getreg(s, DILL_I);
     b = dill_getreg(s, DILL_I);
     f = dill_getreg(s, DILL_I);

     dill_reti(s, c);
     dill_begin_prefix_code(s);
     dill_setp(s, f, (void*)gg);
     dill_seti(s, a, 1);
     dill_seti(s, b, 2);  
     dill_push_init(s);
     if (!dill_do_reverse_vararg_push(s)) {
	 dill_push_argi(s, a);
	 dill_push_argi(s, b);
     } else {
	 dill_push_argi(s, b);
	 dill_push_argi(s, a);
     }
     a = dill_callri(s, f);
     dill_movi(s, c, a);
     h = dill_finalize(s);
     ip = (int(*)())dill_get_fp(h);

     if (verbose) dill_dump(s);

     ret = (*ip)();
     dill_free_stream(s);
     dill_free_handle(h);
     return ret == 3;
}


int 
main(int argc, char **argv)
{
    if (argc > 1) verbose++;
    if (a() == 0) return 1;
    if (c() == 0) return 1;
    return 0;
}
