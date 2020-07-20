DEBUG ?= 0
PAUL ?= 1
CC ?= mpicc
CFLAGS ?= -O2

MKFILEPATH = $(abspath $(lastword $(MAKEFILE_LIST)))
SRCROOT_ ?= $(patsubst %/,%,$(dir $(MKFILEPATH)))
SRCROOT=$(realpath $(SRCROOT_))

GSLIBDIR=$(GSLIBPATH)

SRCDIR    =$(SRCROOT)/src
SORTDIR   =$(SRCROOT)/src/sort
BUILDDIR  =$(SRCROOT)/build
EXAMPLEDIR=$(SRCROOT)/example
TESTDIR   =$(SRCROOT)/tests

TARGET=parRSB
LIB=$(BUILDDIR)/lib/lib$(TARGET).a
EXAMPLE=$(EXAMPLEDIR)/example

INCFLAGS=-I$(SRCDIR) -I$(SORTDIR) -I$(GSLIBDIR)/include
LDFLAGS:=-L$(BUILDDIR)/lib -l$(TARGET) -L $(GSLIBDIR)/lib -lgs -lm

SRCS    =$(wildcard $(SRCDIR)/*.c)
SORTSRCS=$(wildcard $(SORTDIR)/*.c)
TESTSRCS=$(wildcard $(TESTDIR)/*.c)

SRCOBJS =$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(SRCS))
SRCOBJS+=$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(SORTSRCS))
TESTOBJS=$(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%,$(TESTSRCS))

PP=

ifneq ($(DEBUG),0)
  PP += -g -DGENMAP_DEBUG
endif

ifneq ($(PAUL),0)
  PP += -DGENMAP_PAUL
endif

INSTALLDIR=
ifneq (,$(strip $(DESTDIR)))
	INSTALLDIR=$(realpath $(DESTDIR))
endif

.PHONY: default
default: check lib install

.PHONY: all
all: check lib tests example install

.PHONY: install
install: lib
ifneq ($(INSTALLDIR),)
	@mkdir -p $(INSTALLDIR)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALLDIR)/lib 2>/dev/null
	@mkdir -p $(INSTALLDIR)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(SORTDIR)/*.h $(INSTALLDIR)/include 2>/dev/null
endif

.PHONY: lib
lib: $(SRCOBJS)
	@mkdir -p $(BUILDDIR)/lib
	@$(AR) cr $(LIB) $(SRCOBJS)
	@ranlib $(LIB)

.PHONY: check
check: 
ifeq ($(GSLIBPATH),)
  $(error Specify GSLIBPATH=<path to gslib>/build)
endif

$(BUILDDIR)/src/%.o: $(SRCROOT)/src/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

.PHONY: examples
examples: $(EXAMPLE)

$(EXAMPLE): install
	$(CC) $(CFLAGS) -I$(GSLIBDIR)/include -I$(SRCDIR) -I$(SORTDIR) $@.c -o $@ $(LDFLAGS)

.PHONY: tests
tests: install $(TESTOBJS)
	@cp $(TESTDIR)/run-tests.sh $(BUILDDIR)/tests/
	@cd $(BUILDDIR)/tests && ./run-tests.sh

$(BUILDDIR)/tests/%: $(SRCROOT)/tests/%.c
	$(CC) $(CFLAGS) -I$(GSLIBDIR)/include -I$(SRCDIR) -I$(SORTDIR) $< -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	@rm -rf $(BUILDDIR) $(EXAMPLE) $(EXAMPLE).o

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(shell mkdir -p $(BUILDDIR)/src/sort)
$(shell mkdir -p $(BUILDDIR)/tests)
