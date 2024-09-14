#include <assert.h>
#include "vcode.h"

int main() { 
	static unsigned insn[1000];
	struct v_reg al[10];
	v_reg_type r;
	struct v_cstate c;

	v_clambda("foo", V_NLEAF, insn);
	v_getreg(&r, V_I, V_TEMP);

	v_push_init(&c);

	v_push_argpi(c, "Hello: %d %d\n");
	v_seti(r, 10);
	v_push_argi(c, r);
	v_seti(r,20);
	v_push_argi(c, r);
	
	v_ccallv(c, printf);
	v_end().v();
	return 0;
}
