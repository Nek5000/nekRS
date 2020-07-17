MPI ?= 1
MPIIO ?= 1
ADDUS ?= 1
USREXIT ?= 0
LIBNAME ?= gs
BLAS ?= 0
CFLAGS ?= $(LIBP_OPT_FLAGS)
FFLAGS ?= $(LIBP_OPT_FLAGS)

SRCROOT=.
TESTDIR=$(SRCROOT)/tests
FTESTDIR=$(TESTDIR)/fortran
SRCDIR=$(SRCROOT)/src
INCDIR=$(SRCROOT)/src
LIBDIR=$(SRCROOT)/lib

ifneq (,$(strip $(PREFIX)))
INSTALL_ROOT = $(PREFIX)
else
INSTALL_ROOT = $(LIBDIR)
endif

ifneq (0,$(MPI))
  G+=-DMPI
endif

ifneq (0,$(MPIIO))
  ifneq (0,$(MPI))
    G+=-DUSEMPIIO
  endif
endif

ifneq (0,$(ADDUS))
  G+=-DUNDERSCORE
endif

ifneq (0,$(USREXIT))
  G+=-DUSE_USR_EXIT
endif

ifeq (0,$(BLAS))
  G+=-DUSE_NAIVE_BLAS
endif

ifeq (1,$(BLAS))
  G+=-DUSE_CBLAS
endif

ifneq ($(PREFIX),)
  G+=-DPREFIX=$(PREFIX)
endif

ifneq ($(FPREFIX),)
  G+=-DFPREFIX=$(FPREFIX)
endif

G+=-DGLOBAL_LONG_LONG
#G+=-DPRINT_MALLOCS=1
#G+=-DGS_TIMING -DGS_BARRIER

CCCMD=$(CC) $(CFLAGS) -I$(INCDIR) $(G)
FCCMD=$(FC) $(FFLAGS)

LINKCMD=$(CC) $(CFLAGS) -I$(INCDIR) $(G) $^ -o $@ -L$(SRCDIR) \
        -l$(LIBNAME) -lm $(LDFLAGS)

TESTS=$(TESTDIR)/sort_test $(TESTDIR)/sort_test2 $(TESTDIR)/sarray_sort_test \
      $(TESTDIR)/comm_test $(TESTDIR)/crystal_test \
      $(TESTDIR)/sarray_transfer_test $(TESTDIR)/gs_test \
      $(TESTDIR)/gs_test_gop_blocking $(TESTDIR)/gs_test_gop_nonblocking \
      $(TESTDIR)/gs_unique_test $(TESTDIR)/gs_test_old \
      $(TESTDIR)/findpts_el_2_test \
      $(TESTDIR)/findpts_el_2_test2 $(TESTDIR)/findpts_el_3_test \
      $(TESTDIR)/findpts_el_3_test2 $(TESTDIR)/findpts_local_test \
      $(TESTDIR)/findpts_test $(TESTDIR)/poly_test \
      $(TESTDIR)/lob_bnd_test $(TESTDIR)/obbox_test

FTESTS=$(FTESTDIR)/f-igs

GS=$(SRCDIR)/gs.o $(SRCDIR)/sort.o $(SRCDIR)/sarray_transfer.o \
   $(SRCDIR)/sarray_sort.o $(SRCDIR)/gs_local.o $(SRCDIR)/fail.o \
   $(SRCDIR)/crystal.o $(SRCDIR)/comm.o $(SRCDIR)/tensor.o

FWRAPPER=$(SRCDIR)/fcrystal.o $(SRCDIR)/findpts.o
INTP=$(SRCDIR)/findpts_local.o $(SRCDIR)/obbox.o $(SRCDIR)/poly.o \
     $(SRCDIR)/lob_bnd.o $(SRCDIR)/findpts_el_3.o $(SRCDIR)/findpts_el_2.o

.PHONY: all lib install deps tests clean objects odepinfo

all : lib tests

lib: $(GS) $(FWRAPPER) $(INTP) $(SRCDIR)/rand_elt_test.o
	@$(AR) cr $(SRCDIR)/lib$(LIBNAME).a $?
	@ranlib $(SRCDIR)/lib$(LIBNAME).a

install: lib
	@mkdir -p $(INSTALL_ROOT) 2>/dev/null
	@cp -v $(SRCDIR)/lib$(LIBNAME).a $(INSTALL_ROOT) 2>/dev/null

tests: $(TESTS) $(FTESTS)

clean: ; @$(RM) $(SRCDIR)/*.o $(SRCDIR)/*.s $(SRCDIR)/*.a $(TESTS) $(TESTS)/*.o $(FTESTS) $(FTESTS)/*.o

cmds: ; @echo CC = $(CCCMD); echo LINK = $(LINKCMD);

deps: ; ./cdep.py *.c > makefile.cdep;

odepinfo: deps objects; @./odep_info.py *.o

$(TESTS): % : %.o | lib
	$(LINKCMD)

$(FTESTS): % : %.o | lib
	$(FCCMD) $^ -o $@ -L$(SRCDIR) -l$(LIBNAME)

-include makefile.cdep

%.o: %.c ; $(CCCMD) -c $< -o $@
%.o: %.f ; $(FCCMD) -c $< -o $@
%.s: %.c ; $(CCCMD) -S $< -o $@
objects: $(OBJECTS) ;

#poly_imp.h: gen_poly_imp.c
#	$(RM) poly_imp.h;
#	$(CC) -lgmp -lm gen_poly_imp.c -o gen_poly_imp;
#	./gen_poly_imp > poly_imp.h;
#	$(RM) gen_poly_imp
