MPI?=1
VALGRIND?=1
DEBUG?=0
UNDERSCORE?=1
CC?=mpicc
FC?=mpif77

SRCROOT=.

GSLIBDIR=$(GSLIBPATH)

SRCDIR  =$(SRCROOT)/src
INCDIR  =$(SRCDIR)
BUILDDIR=$(SRCROOT)/build
TESTDIR =$(SRCROOT)/tests

TARGET=parRSB
TESTS=$(TESTDIR)/gmsh/gmsh-test $(TESTDIR)/con/con-test
LIB=src/lib$(TARGET).a

INCFLAGS=-I$(INCDIR) -I$(GSLIBDIR)/include

ifneq (,$(strip $(DESTDIR)))
INSTALL_ROOT = $(DESTDIR)
else
INSTALL_ROOT = $(SRCROOT)/build
endif

TESTLDFLAGS:=-L$(SRCDIR) -l$(TARGET) -L $(GSLIBDIR)/lib -lgs -lm $(LDFLAGS)

# Main source
CSRCS:= $(SRCDIR)/genmap-vector.c $(SRCDIR)/genmap-algo.c \
	$(SRCDIR)/genmap-io.c $(SRCDIR)/genmap-comm.c $(SRCDIR)/genmap.c \
	$(SRCDIR)/parRSB.c $(SRCDIR)/genmap-quality.c \
	$(SRCDIR)/genmap-chelpers.c $(SRCDIR)/genmap-gmsh.c

COBJS:=$(CSRCS:.c=.o)

FSRCS:=
FOBJS:=$(FSRCS:.f=.o)

SRCOBJS :=$(COBJS) $(FOBJS)

PP= -DGENMAP_LONG_LONG
ifneq ($(MPI),0)
  PP += -DGENMAP_MPI
endif

ifneq ($(DEBUG),0)
  PP += -g -DGENMAP_DEBUG
endif

ifneq ($(UNDERSCORE),0)
  PP += -DGENMAP_UNDERSCORE
endif

.PHONY: default
default: check lib install

.PHONY: all
all: check lib tests install

.PHONY: install
install: lib
	@mkdir -p $(INSTALL_ROOT)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALL_ROOT)/lib 2>/dev/null
	@mkdir -p $(INSTALL_ROOT)/include 2>/dev/null
	@cp $(SRCDIR)/parRSB.h $(INSTALL_ROOT)/include 2>/dev/null


.PHONY: $(TARGET)
lib: $(SRCOBJS)
	@$(AR) cr $(LIB) $(SRCOBJS)
	@ranlib $(LIB)

.PHONY: check
check: 
ifeq ($(GSLIBPATH),)
	$(error Specify GSLIBPATH=<path to gslib>)
endif

$(COBJS): %.o: %.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

$(FOBJS): %.o: %.f
	$(FC) $(FFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

.PHONY: tests
tests: $(TESTS)

$(TESTS): lib
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $@.c -o $@ $(TESTLDFLAGS)

.PHONY: clean
clean:
	@rm -f $(SRCOBJS) $(LIB) $(TESTS) $(TESTS).o

.PHONY: astyle
astyle:
	astyle --style=google --indent=spaces=2 --max-code-length=72 \
	    --keep-one-line-statements --keep-one-line-blocks --lineend=linux \
            --suffix=none --preserve-date --formatted --pad-oper \
	    --unpad-paren tests/*.[ch] src/*.[ch] inc/*.[ch] io/*/*.[ch]

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true
