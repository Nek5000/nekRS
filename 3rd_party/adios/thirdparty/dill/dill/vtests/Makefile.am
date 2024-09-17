

AUTOMAKE_OPTIONS = no-dependencies foreign

if BUILD_LIBRARY
EXTRA_PROGRAMS = basic_call general t1 opt branch pkg_test multi_test prefix_test

TESTS = general$(EXEEXT) t1$(EXEEXT) opt branch pkg_test multi_test prefix_test

check-local: basic_call.stamp

valgrind:  basic_call.valgrind general.valgrind t1.valgrind opt.valgrind branch.valgrind pkg_test.valgrind multi_test.valgrind
endif BUILD_LIBRARY

VALGRIND= ../libtool --mode=execute valgrind --tool=memcheck --error-exitcode=1 --quiet --trace-children=yes

LDADD = ../libdill.la

general_LDADD = ../libdill.la ${MATH_LIB} 

NOOPT_COMPILE = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
	$(CPPFLAGS) $(AM_CFLAGS) $(NOOPT_CFLAGS)

general.o: general.c
	$(NOOPT_COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

general.c:  general.ops
	perl $(srcdir)/general.ops

basic_call.stamp:  basic_call@EXEEXT@
	./basic_call | head -200 1>tmp.out  2>&1
	@grep -v "lib[vi]code.a" tmp.out > tmp1.out
	@mv tmp1.out tmp.out
	diff -w tmp.out $(srcdir)/basic_call.out && touch basic_call.stamp && echo PASS: basic_call
	@rm tmp.out

basic_call.valgrind:  basic_call@EXEEXT@
	$(VALGRIND) ./basic_call 1>/dev/null  2>&1

general.valgrind:  general@EXEEXT@
	$(VALGRIND) ./general 1>/dev/null  2>&1

t1.valgrind:  t1@EXEEXT@
	$(VALGRIND) ./t1 1>/dev/null  2>&1

opt.valgrind:  opt@EXEEXT@
	$(VALGRIND) ./opt 1>/dev/null  2>&1

branch.valgrind:  branch@EXEEXT@
	$(VALGRIND) ./branch 1>/dev/null  2>&1

pkg_test.valgrind:  pkg_test@EXEEXT@
	$(VALGRIND) ./pkg_test 1>/dev/null  2>&1

multi_test.valgrind:  multi_test@EXEEXT@
	$(VALGRIND) ./multi_test 1>/dev/null  2>&1

# t2.stamp:  t2@EXEEXT@
# 	./t2 1>tmp.out 2>&1
# 	@grep -v "lib[vi]code.a" tmp.out > tmp1.out
# 	@mv tmp1.out tmp.out
# 	diff -w tmp.out t2.out && touch t2.stamp
# 	@rm tmp.out

clean-local:
	rm -f t4.stamp t2.stamp ${EXTRA_PROGRAMS} general.c

