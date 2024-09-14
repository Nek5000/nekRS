
AUTOMAKE_OPTIONS = no-dependencies foreign

if NATIVE_CG
EXTRA_PROGRAMS = regress ctest stest call-test t1 pkg_test cplus
check-local: test
endif NATIVE_CG

noinst_PROGRAMS = 

LDADD = ../libdill.la

TESTS = cplus

cplus_SOURCES= cplus.cc

VALGRIND= ../libtool --mode=execute valgrind --tool=memcheck --error-exitcode=1 --quiet --trace-children=yes

test: call-test.stamp stest.stamp ctest.stamp regress.stamp t1.stamp pkg_test.stamp

valgrind: call-test.valgrind stest.valgrind ctest.valgrind regress.valgrind t1.valgrind pkg_test.valgrind

regress_LDADD = ../libdill.la ${MATH_LIB} 

regress.c:  test-gen
	perl -w $(srcdir)/test-gen @TEST_PERL_FLAGS@ > regress.c

call-test.c:  call-gen
	perl -w $(srcdir)/call-gen @TEST_PERL_FLAGS@ > call-test.c

call-test.stamp:  call-test@EXEEXT@
	./call-test 1>tmp.out  2>&1
	diff -w tmp.out $(srcdir)/call-test.out && touch call-test.stamp
	@rm tmp.out

call-test.valgrind:  call-test@EXEEXT@
	$(VALGRIND) ./call-test 1>/dev/null  2>&1

stest.stamp:  stest@EXEEXT@
	./stest 1>tmp.out  2>&1
	diff -w tmp.out $(srcdir)/stest.out && touch stest.stamp
	@rm tmp.out

stest.valgrind:  stest@EXEEXT@
	$(VALGRIND) ./stest 1>/dev/null  2>&1

ctest.stamp:  ctest@EXEEXT@
	./ctest @TEST_PERL_FLAGS@ 1>tmp.out  2>&1
	diff -w tmp.out $(srcdir)/ctest.out && touch ctest.stamp
	@rm tmp.out

ctest.valgrind:  ctest@EXEEXT@
	$(VALGRIND) ./ctest @TEST_PERL_FLAGS@ 1>/dev/null  2>&1

regress.stamp:  regress@EXEEXT@
	./regress 1>tmp.out  2>&1
	diff -w tmp.out $(srcdir)/regress.out && touch regress.stamp
	@rm tmp.out

regress.valgrind:  regress@EXEEXT@
	$(VALGRIND) ./regress 1>/dev/null  2>&1

t1.stamp: t1@EXEEXT@
	./t1 && touch t1.stamp

t1.valgrind:  t1@EXEEXT@
	$(VALGRIND) ./t1 1>/dev/null  2>&1

pkg_test.stamp: pkg_test@EXEEXT@
	./pkg_test && touch pkg_test.stamp

pkg_test.valgrind:  pkg_test@EXEEXT@
	$(VALGRIND) ./pkg_test 1>/dev/null  2>&1

# t2.stamp:  t2@EXEEXT@
# 	./t2 1>tmp.out 2>&1
# 	@grep -v "lib[vi]code.a" tmp.out > tmp1.out
# 	@mv tmp1.out tmp.out
# 	diff -w tmp.out t2.out && touch t2.stamp
# 	@rm tmp.out

clean-local:
	rm -f *.stamp ${EXTRA_PROGRAMS} regress.c call-test.c

