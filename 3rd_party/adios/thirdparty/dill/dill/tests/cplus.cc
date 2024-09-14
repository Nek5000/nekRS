#include "config.h"
#include "dill.h"

#include <iostream>

static int verbose = 0;

class C
{
  int state;
public:
  C (int);
  void sub(int in);
  virtual void add(int add);
  virtual int get();
};

C::C(int init) {
  state = init;
}

void
C::sub(int in) {
  state -= in;
  if (verbose) {
    std::cout << "Called sub " << in << std::endl;
  }
}
void
C::add(int add) {
  state += add;
  if (verbose) {
    std::cout << "Called add " << add << std::endl;
  }
}
int 
C::get() {
  if (verbose) {
    std::cout << "Called get" << std::endl;
  }
  return state;
}

int c_based_calls(C*c)
{
  void (C::*p)(int);
  void (*x)(void *th, int in);
  int (C::*p2)();
  int (*x2)(void *th);

  p = &C::sub; // prepare to call C::sub
  void * th = get_this_ptr(&p, c);
  x = (void(*)(void*,int)) get_xfer_ptr(&p, c);
  x(th, 5);

  p = &C::add; // prepare to call C::add
  th = get_this_ptr(&p, c);
  x = (void(*)(void*,int)) get_xfer_ptr(&p, c);
  x(th, 10);

  p2 = &C::get;   // prepare to call C::get
  th = get_this_ptr(&p2, c);
  x2 = (int(*)(void*)) get_xfer_ptr(&p2, c);
  return x2(th);
};

int dcg_calls(C*c)
{
  dill_exec_handle h;
  void (C::*p)(int);
  void *x;
  int (C::*p2)();
  void *x2;
  dill_stream s = dill_create_raw_stream();
  dill_reg this_reg;
  dill_reg ret_reg;
  int (*proc)();
  int result;

  p = &C::sub; // prepare to call C::sub
  void * th = get_this_ptr(&p, c);
  x = get_xfer_ptr(&p, c);

  dill_start_simple_proc(s, "foo", DILL_V);
  dill_raw_getreg(s, &this_reg, DILL_P, DILL_TEMP);

  dill_setp(s, this_reg, th);
  dill_scallv(s, x, "C::sub", "%p%I", this_reg, 5);
  //  x(th, 5);

  p = &C::add; // prepare to call C::add
  th = get_this_ptr(&p, c);
  x = get_xfer_ptr(&p, c);
  dill_setp(s, this_reg, th);
  dill_scallv(s, x, "C::add", "%p%I", this_reg, 10);
  //  x(th, 10);

  p2 = &C::get;   // prepare to call C::get
  th = get_this_ptr(&p2, c);
  x2 = get_xfer_ptr(&p2, c);
  dill_setp(s, this_reg, th);
  ret_reg = dill_scalli(s, x2, "C::add", "%p", this_reg);
  dill_reti(s, ret_reg);
  h = dill_finalize(s);
  proc = (int(*)()) dill_get_fp(h);
  result = proc();
  dill_free_handle(h);
  dill_free_stream(s);
  return result;
};

int main(int argc, char **argv)
{
  C c(25);
  int no_dcg = 0;
  int result;
  for (int i=1; i < argc; i++) {
    std::string s = argv[i];
    if (s == "-no_dcg") {
      no_dcg++;
    } else if (s == "-v") {
      verbose++;
    } else if (s == "-verbose") {
      verbose++;
    } else {
      std::cout << "Unknown argument \"" << s << "\"" << std::endl;
    }
  }

  if (no_dcg) {
    result = c_based_calls(&c);
  } else {
    result = dcg_calls(&c);
  }    
  if (result != 30) std::cout << "FAILURE!" << std::endl;
};
