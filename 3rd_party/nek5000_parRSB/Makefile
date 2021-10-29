## General build parameters ##
DEBUG ?= 0
MPI ?= 1
CC ?= mpicc
CFLAGS ?= -g -O0
UNDERSCORE ?= 1
SYNC_BY_REDUCTION ?= 1
BLAS ?= 0
BLASDIR ?=
BLASFLAGS ?= -lblas -llapack

## Don't touch what follows ##
MKFILEPATH = $(abspath $(lastword $(MAKEFILE_LIST)))
SRCROOT_ ?= $(patsubst %/,%,$(dir $(MKFILEPATH)))
SRCROOT = $(realpath $(SRCROOT_))

BUILDDIR = $(SRCROOT)/build

SRCDIR = $(SRCROOT)/src
SRCS = $(wildcard $(SRCDIR)/*.c)
SRCOBJS = $(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(SRCS))

SORTDIR = $(SRCROOT)/src/sort
SORTSRCS = $(wildcard $(SORTDIR)/*.c)
SRCOBJS += $(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(SORTSRCS))

PRECONDDIR = $(SRCROOT)/src/precond
PRECONDSRCS = $(wildcard $(PRECONDDIR)/*.c)
SRCOBJS += $(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(PRECONDSRCS))

GENCONDIR = $(SRCROOT)/src/gencon
GENCONSRCS = $(wildcard $(GENCONDIR)/*.c)
SRCOBJS += $(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%.o,$(GENCONSRCS))

EXAMPLEDIR = $(SRCROOT)/examples
EXAMPLESRCS = $(wildcard $(EXAMPLEDIR)/*.c)
EXAMPLEOBJS = $(patsubst $(SRCROOT)/%.c,$(BUILDDIR)/%,$(EXAMPLESRCS))

TARGET = parRSB
LIB = $(BUILDDIR)/lib/lib$(TARGET).a

INCFLAGS = -I$(SRCDIR) -I$(SORTDIR) -I$(PRECONDDIR) -I$(GENCONDIR) -I$(GSLIBPATH)/include
LDFLAGS = -L$(BUILDDIR)/lib -l$(TARGET) -L$(GSLIBPATH)/lib -lgs -lm

PP=

ifneq ($(DEBUG),0)
  PP += -g -DGENMAP_DEBUG
endif

ifneq ($(UNDERSCORE),0)
  PP += -DGENMAP_UNDERSCORE
endif

ifneq ($(BLAS),0)
  PP += -DGENMAP_BLAS
  ifneq ($(BLASDIR),)
    LDFLAGS+= -L$(BLASDIR)
  endif
  LDFLAGS += $(BLASFLAGS)
endif

ifneq ($(MPI),0)
  PP += -DMPI
endif

ifneq ($(SYNC_BY_REDUCTION),0)
  PP += -DGENMAP_SYNC_BY_REDUCTION
endif

INSTALLDIR=
ifneq (,$(strip $(DESTDIR)))
  INSTALLDIR = $(realpath $(DESTDIR))
endif

.PHONY: all
all: check lib examples install

.PHONY: check
check:
ifeq ($(GSLIBPATH),)
  $(error Specify GSLIBPATH=<path to gslib>/build)
endif

.PHONY: install
install: lib
ifneq ($(INSTALLDIR),)
	@mkdir -p $(INSTALLDIR)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALLDIR)/lib 2>/dev/null
	@mkdir -p $(INSTALLDIR)/include 2>/dev/null
	@cp $(SRCDIR)/*.h $(SORTDIR)/*.h $(PRECONDDIR)/*.h \
		$(INSTALLDIR)/include 2>/dev/null
endif

.PHONY: lib
lib: $(SRCOBJS)
	@mkdir -p $(BUILDDIR)/lib
	@$(AR) cr $(LIB) $(SRCOBJS)
	@ranlib $(LIB)

$(BUILDDIR)/src/%.o: $(SRCROOT)/src/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) -c $< -o $@

.PHONY: examples
examples: install $(EXAMPLEOBJS)

$(BUILDDIR)/examples/%: $(SRCROOT)/examples/%.c
	$(CC) $(CFLAGS) $(PP) $(INCFLAGS) $< -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	@rm -rf $(BUILDDIR)

.PHONY: format
format:
	find . -iname *.h -o -iname *.c | xargs clang-format -i

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(shell mkdir -p $(BUILDDIR)/src/sort)
$(shell mkdir -p $(BUILDDIR)/src/precond)
$(shell mkdir -p $(BUILDDIR)/src/gencon)
$(shell mkdir -p $(BUILDDIR)/examples)
