MPI?=1
DEBUG?=0
UNDERSCORE?=1
CC?=mpicc
FC?=mpif77
PAUL?=1

SRCROOT=.
GSLIBDIR=$(GSLIBPATH)

SRCDIR  =$(SRCROOT)/src
BUILDDIR=$(SRCROOT)/build
TESTDIR =$(SRCROOT)/example

TARGET=parRSB
TESTS=$(TESTDIR)/example
LIB=src/lib$(TARGET).a

INCFLAGS=-I$(SRCDIR) -I$(GSLIBDIR)/include

TESTLDFLAGS:=-L$(SRCDIR) -l$(TARGET) -L $(GSLIBDIR)/lib -lgs -lm $(LDFLAGS)

ifneq (,$(strip $(DESTDIR)))
INSTALL_ROOT = $(DESTDIR)
else
INSTALL_ROOT = $(SRCROOT)/build
endif

CSRCS:= $(SRCDIR)/genmap.c \
  	$(SRCDIR)/genmap-vector.c $(SRCDIR)/genmap-handle.c $(SRCDIR)/genmap-comm.c \
	$(SRCDIR)/genmap-eigen.c $(SRCDIR)/genmap-laplacian.c $(SRCDIR)/genmap-lanczos.c \
	$(SRCDIR)/genmap-rsb.c \
	$(SRCDIR)/genmap-quality.c $(SRCDIR)/genmap-chelpers.c \
	$(SRCDIR)/parRSB.c 
COBJS:=$(CSRCS:.c=.o)

SRCOBJS:=$(COBJS)

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

ifneq ($(PAUL),0)
  PP += -DGENMAP_PAUL
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
	$(error Specify GSLIBPATH=<path to gslib>/build)
endif

$(COBJS): %.o: %.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

.PHONY: tests
tests: $(TESTS)

$(TESTS): lib
	$(CC) $(CFLAGS) -DGLOBAL_LONG_LONG -DMPI -DPREFIX=gslib_ $(INCFLAGS) $@.c -o $@ $(TESTLDFLAGS)

.PHONY: clean
clean:
	@rm -f $(SRCOBJS) $(LIB) $(TESTS) $(TESTS).o

.PHONY: astyle
astyle:
	astyle --style=google --indent=spaces=2 --max-code-length=80 \
	    --keep-one-line-statements --keep-one-line-blocks --lineend=linux \
            --suffix=none --preserve-date --formatted --pad-oper \
	    --unpad-paren tests/*/*.[ch] src/*.[ch]

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true
