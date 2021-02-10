MPI ?= 1
ADDUS ?= 1
USREXIT ?= 0
NBC ?= 0
LIBNAME ?= gs
BLAS ?= 0
DEBUG ?= 0
CFLAGS ?= -O2
FFLAGS ?= -O2
CPREFIX ?= gslib_
FPREFIX ?= fgslib_

SRCROOT=.
TESTDIR=$(SRCROOT)/tests
FTESTDIR=$(TESTDIR)/fortran
SRCDIR=$(SRCROOT)/src
INCDIR=$(SRCROOT)/src
LIBDIR=$(SRCROOT)/lib

ifneq (,$(strip $(DESTDIR)))
INSTALL_ROOT = $(DESTDIR)
else
INSTALL_ROOT = $(SRCROOT)/build
endif

$(shell >config.h)
ifneq (0,$(MPI))
  SN=MPI
  G:=$(G) -D$(SN)
  ifeq ($(origin CC),default)
    CC = mpicc
  endif
  ifeq ($(origin FC),default)
    FC = mpif77
  endif
  $(shell printf "#ifndef ${SN}\n#define ${SN}\n#endif\n" >>config.h)
endif

ifneq (0,$(ADDUS))
  SN=UNDERSCORE
  G:=$(G) -D$(SN)
  $(shell printf "#ifndef ${SN}\n#define ${SN}\n#endif\n" >>config.h)
endif

SN=GLOBAL_LONG_LONG
G:=$(G) -D$(SN)
$(shell printf "#ifndef ${SN}\n#define ${SN}\n#endif\n" >>config.h)

SN=PREFIX
G:=$(G) -D$(SN)=$(CPREFIX)
$(shell printf "#ifndef ${SN}\n#define ${SN} ${CPREFIX}\n#endif\n" >>config.h)

SN=FPREFIX
G:=$(G) -D$(SN)=$(FPREFIX)
$(shell printf "#ifndef ${SN}\n#define ${SN} ${FPREFIX}\n#endif\n" >>config.h)

ifneq (0,$(USREXIT))
  G+=-DUSE_USR_EXIT
endif

ifneq (0,$(NBC))
  G+=-DUSE_NBC
endif

ifeq (0,$(BLAS))
  SN=USE_NAIVE_BLAS
  G:=$(G) -D$(SN)
  $(shell printf "#ifndef ${SN}\n#define ${SN}\n#endif\n" >>config.h)
endif

ifeq (1,$(BLAS))
  G+=-DUSE_CBLAS
endif

ifneq (0,$(DEBUG))
  G+=-DGSLIB_DEBUG
  CFLAGS+=-g
endif

CCCMD=$(CC) $(CFLAGS) -I$(INCDIR) $(G)
FCCMD=$(FC) $(FFLAGS) -I$(INCDIR) $(G)

TESTS=$(TESTDIR)/sort_test $(TESTDIR)/sort_test2 $(TESTDIR)/sarray_sort_test \
      $(TESTDIR)/comm_test $(TESTDIR)/crystal_test \
      $(TESTDIR)/sarray_transfer_test $(TESTDIR)/gs_test \
      $(TESTDIR)/gs_test_gop_blocking $(TESTDIR)/gs_test_gop_nonblocking \
      $(TESTDIR)/gs_unique_test \
      $(TESTDIR)/findpts_el_2_test \
      $(TESTDIR)/findpts_el_2_test2 $(TESTDIR)/findpts_el_3_test \
      $(TESTDIR)/findpts_el_3_test2 $(TESTDIR)/findpts_local_test \
      $(TESTDIR)/findpts_test $(TESTDIR)/findpts_test_ms $(TESTDIR)/poly_test \
      $(TESTDIR)/lob_bnd_test $(TESTDIR)/obbox_test

FTESTS=$(FTESTDIR)/f-igs

GS=$(SRCDIR)/gs.o $(SRCDIR)/sort.o $(SRCDIR)/sarray_transfer.o \
   $(SRCDIR)/sarray_sort.o $(SRCDIR)/gs_local.o $(SRCDIR)/fail.o \
   $(SRCDIR)/crystal.o $(SRCDIR)/comm.o $(SRCDIR)/tensor.o

FWRAPPER=$(SRCDIR)/fcrystal.o $(SRCDIR)/findpts.o
INTP=$(SRCDIR)/findpts_local.o $(SRCDIR)/obbox.o $(SRCDIR)/poly.o \
     $(SRCDIR)/lob_bnd.o $(SRCDIR)/findpts_el_3.o $(SRCDIR)/findpts_el_2.o

.PHONY: all lib install tests clean objects

all : lib install

lib: $(GS) $(FWRAPPER) $(INTP)
	@$(AR) cr $(SRCDIR)/lib$(LIBNAME).a $?
	@ranlib $(SRCDIR)/lib$(LIBNAME).a

install: lib
	@mkdir -p $(INSTALL_ROOT)/lib 2>/dev/null
	@cp -v $(SRCDIR)/lib$(LIBNAME).a $(INSTALL_ROOT)/lib 2>/dev/null
	@mkdir -p $(INSTALL_ROOT)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(INSTALL_ROOT)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(INSTALL_ROOT)/include 2>/dev/null
	@mv config.h $(INSTALL_ROOT)/include 2>/dev/null

tests: $(TESTS)

clean: ; @$(RM) config.h $(SRCDIR)/*.o $(SRCDIR)/*.s $(SRCDIR)/*.a $(TESTDIR)/*.o $(FTESTDIR)/*.o $(TESTS)

$(TESTS): % : %.c | lib install
	$(CC) $(CFLAGS) -I$(INSTALL_ROOT)/include $< -o $@ -L$(INSTALL_ROOT)/lib -l$(LIBNAME) -lm $(LDFLAGS) 

$(FTESTS): % : %.o | lib install
	$(FCCMD) $^ -o $@ -L$(SRCDIR) -l$(LIBNAME)

%.o: %.c ; $(CCCMD) -c $< -o $@
%.o: %.f ; $(FCCMD) -c $< -o $@
%.s: %.c ; $(CCCMD) -S $< -o $@
objects: $(OBJECTS) ;
