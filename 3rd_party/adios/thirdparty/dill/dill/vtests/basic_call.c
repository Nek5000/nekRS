#include "stdio.h"
#include "dill.h"

static int verbose = 0;

int gg(int a, int b) {
    printf("In gg  a=%d, b=%d\n", a, b);
     return a+b;
}
int ff(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j) {
    printf("In ff  a=%d, b=%d, c=%d, d=%d, e=%d, f=%d, g=%d, h=%d, i=%d, j=%d\n",
	   a, b, c, d, e, f, g, h, i, j);
     return a+b+c+d+e+f+g+h+i+j;
}

void a () {
    dill_stream s = dill_create_stream();
     dill_exec_handle handle;
     dill_reg a,b,p3,d,e,f,g,h,i,j,w,z,cnt;
/*     dill_reg func;*/
     int L1;
     int (*ip)();

     dill_start_simple_proc(s, "a_gen", DILL_I);
     L1 = dill_alloc_label(s, NULL);
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
     dill_seti(s, a, 1);
     dill_seti(s, b, 2);
     dill_seti(s, p3, 3);
     dill_seti(s, d, 4);
     dill_seti(s, e, 5);
     dill_seti(s, f, 6);
     dill_seti(s, g, 7);
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
     handle = dill_finalize(s);
     ip = (int(*)())dill_get_fp(handle);

     if (verbose) dill_dump(s);

     printf("**135=%d\n", (*ip)());
     dill_free_handle(handle);
     dill_free_stream(s);
}

void b () {
     dill_reg f;
     void *(*pp)();
     dill_exec_handle h;
     
     dill_stream s = dill_create_stream();

     dill_start_proc(s, "b_gen", DILL_P, "%i%i");
     f = dill_getreg(s, DILL_P);
     dill_setp(s, f, (void *)gg);
     dill_retp(s, f);
     h = dill_finalize(s);
     pp = (void *(*)())dill_get_fp(h);

     if (verbose) dill_dump(s);

     printf("**3=%d\n", (*(int (*)(int, int))(IMM_TYPE)((*pp)()))(1,2));
     dill_free_handle(h);
     dill_free_stream(s);
}

void c () {
     dill_reg a,b,f;
     int (*ip)();
     dill_exec_handle h;
     
     dill_stream s = dill_create_stream();

     dill_start_simple_proc(s, "c_gen", DILL_I);
     a = dill_getreg(s, DILL_I);
     b = dill_getreg(s, DILL_I);
     f = dill_getreg(s, DILL_I);

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
     dill_reti(s, a);
     h = dill_finalize(s);
     ip = (int(*)())dill_get_fp(h);

     if (verbose) dill_dump(s);

     printf("**3=%d\n", (*ip)());
     dill_free_handle(h);
     dill_free_stream(s);
}

void d () {
     dill_reg a,b;
     int (*ip)();
     dill_exec_handle h;
     
     dill_stream s = dill_create_stream();

     dill_start_simple_proc(s, "d_gen", DILL_I);
     a = dill_getreg(s, DILL_I);
     b = dill_getreg(s, DILL_I);

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
     a = dill_calli(s, (void*)gg, "gg");
     dill_reti(s, a);
     h = dill_finalize(s);
     ip = (int(*)())dill_get_fp(h);

     if (verbose) dill_dump(s);

     printf("**3=%d\n", (*ip)());
     dill_free_handle(h);
     dill_free_stream(s);
}

void e () {
     dill_reg a,b,p3,d,e,f,g,h,i,j,k,l;
     int (*ip)(int (*)(int,int),int,int,int,int,int,int,int,int,int);
     dill_exec_handle handle;
     
     dill_stream s = dill_create_stream();

     dill_start_proc(s, "e_gen", DILL_I, "%p%i%i%i%i%i%i%i%i%i");
     f = dill_vparam(s, 0);
     a = dill_vparam(s, 1);
     b = dill_vparam(s ,2);
     p3 = dill_vparam(s, 3);
     d = dill_vparam(s, 4);
     e = dill_vparam(s, 5);
     g = dill_vparam(s, 6);
     h = dill_vparam(s, 7);
     i = dill_vparam(s, 8);
     j = dill_vparam(s, 9);
     k = dill_getreg(s, DILL_I);
     l = dill_getreg(s, DILL_I);

     dill_push_init(s);
     if (!dill_do_reverse_vararg_push(s)) {
	 dill_push_argi(s, a);
	 dill_push_argi(s, b);
     } else {
	 dill_push_argi(s, b);
	 dill_push_argi(s, a);
     }
     k = dill_callri(s, f);
     dill_push_init(s);
     if (!dill_do_reverse_vararg_push(s)) {
	 dill_push_argi(s, p3);
	 dill_push_argi(s, d);
     } else {
	 dill_push_argi(s, d);
	 dill_push_argi(s, p3);
     }
     d = dill_callri(s, f);
     dill_addi(s, k,k,d);
     dill_push_init(s);
     if (!dill_do_reverse_vararg_push(s)) {
	 dill_push_argi(s, e);
	 dill_push_argi(s, g);
     } else {
	 dill_push_argi(s, g);
	 dill_push_argi(s, e);
     }
     g = dill_callri(s, f);
     dill_addi(s, k,k,g);
     dill_push_init(s);
     if (!dill_do_reverse_vararg_push(s)) {
	 dill_push_argi(s, h);
	 dill_push_argi(s, i);
     } else {
	 dill_push_argi(s, i);
	 dill_push_argi(s, h);
     }
     i = dill_callri(s, f);
     dill_addi(s, k,k,i);
     dill_addi(s, k,k,j);
     dill_addii(s, l,k,3);
     dill_reti(s, l);
     handle = dill_finalize(s);
     ip = (int(*)(int (*)(int, int), int, int, int, int, int, int, int, int, int))dill_get_fp(handle);

     if (verbose) dill_dump(s);

     printf("**48=%d\n", (*ip)(gg,1,2,3,4,5,6,7,8,9));
     dill_free_handle(handle);
     dill_free_stream(s);
}

void f () {
     dill_reg a, b, d;
     double dp = 3.14159;
     int ip = 5;
     char *pp = "hello!";
     void *(*proc)();
     dill_exec_handle h;
     
     dill_stream s = dill_create_stream();

     dill_start_proc(s, "f_gen", DILL_I, "");
     a = dill_getreg(s, DILL_D);
     b = dill_getreg(s, DILL_I);
     d = dill_getreg(s, DILL_P);
     dill_setd(s, a, dp);
     dill_seti(s, b, ip);
     dill_setp(s, d, pp);
     dill_push_varidiac_init(s, 1);
     if (!dill_do_reverse_vararg_push(s)) {
	 dill_push_argpi(s, (void*)"values are %d, %g, %s\n");
	 dill_push_argi(s, b);
	 dill_push_argd(s, a);
	 dill_push_argp(s, d);
     } else {
	 dill_push_argp(s, d);
	 dill_push_argd(s, a);
	 dill_push_argi(s, b);
	 dill_push_argpi(s, (void*)"values are %d, %g, %s\n");
     }
     a = dill_calli(s, (void*)printf, "printf");
     dill_reti(s, a);
     h = dill_finalize(s);
     proc = (void *(*)())dill_get_fp(h);

     if (verbose) dill_dump(s);

     printf("expect: values are %d, %g, %s\n", ip, dp, pp);
     proc();
     dill_free_handle(h);
     dill_free_stream(s);
}

#ifdef NOTDEF
void g () {
     dill_reg a, b, d, e, f, g, h, i;
     double da = 3.1;
     double db = 4.1;
     double dd = 5.1;
     double de = 6.1;
     double df = 7.1;
     double dg = 8.1;
     double dh = 9.1;
     double di = 10.1;
     void *(*proc)();
     dill_exec_handle h;
     
     dill_stream s = dill_create_stream();

     dill_start_proc(s, "g_gen", "");
     a = dill_getreg(s, DILL_D);
     b = dill_getreg(s, DILL_D);
     d = dill_getreg(s, DILL_D);
     e = dill_getreg(s, DILL_D);
     f = dill_getreg(s, DILL_D);
     g = dill_getreg(s, DILL_D);
     h = dill_getreg(s, DILL_D);
     i = dill_getreg(s, DILL_D);
     dill_setd(s, a, da);
     dill_setd(s, b, db);
     dill_setd(s, d, dd);
     dill_setd(s, e, de);
     dill_setd(s, f, df);
     dill_setd(s, g, dg);
     dill_setd(s, h, dh);
     dill_setd(s, i, di);
     dill_push_init(s);
     if (!dill_do_reverse_vararg_push(s)) {
	 dill_push_argpi(s, (void*)"values are %g, %g, %g, %g, %g, %g, %g, %g, %g\n");
	 dill_push_argd(s, a);
	 dill_push_argd(s, b);
	 dill_push_argd(s, d);
	 dill_push_argd(s, e);
	 dill_push_argd(s, f);
	 dill_push_argd(s, g);
	 dill_push_argd(s, h);
	 dill_push_argd(s, i);
     }
     a = dill_calli(s, (void*)printf, "printf");
     dill_reti(s, a);
     h = dill_finalize(c);
     proc = (void *(*)())dill_get_fp(h);

     if (verbose) dill_dump(s);

     printf("expect: values are %g, %g, %g, %g, %g, %g, %g, %g\n", da, db, dd, de, df, dg, dh, di);
     proc();
     dill_free_handle(h);
     dill_free_stream(s);
}
#endif

int 
main(int argc, char **argv)
{
    if (argc > 1) verbose++;
    printf("########## A\n");
    a();
    printf("########## B\n");
    b();
    printf("########## C\n");
    c();
    printf("########## D\n");
    d();
    printf("########## E\n");
    e();
    printf("########## F\n");
    f();
    printf("########## end\n");
    return 0;
}
