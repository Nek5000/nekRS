ifndef OCCA_DIR
  $(error "Error, environment variable [OCCA_DIR] is not set")
endif

ifndef NEKRS_INSTALL_DIR
  $(error "Error, environment variable [NEKRS_INSTALL_DIR] is not set")
endif

ifndef NEKRS_LIBP_DEFINES
  $(error "Error, environment variable [NEKRS_LIBP_DEFINES] is not set")
endif

# define variables
HDRDIR      = $(NEKRS_BIN_THIRD_PARTY_DIR)/libparanumal
OGSDIR      = $(HDRDIR)/libs/gatherScatter
GSDIR       = $(HDRDIR)/3rdParty/gslib
ALMONDDIR   = $(HDRDIR)/libs/parAlmond
ELLIPTICDIR = $(HDRDIR)/solvers/elliptic
BLASDIR     = $(HDRDIR)/3rdParty/BlasLapack

# define LIBP_OPT_FLAGS
include ${NEKRS_SRC_THIRD_PARTY_DIR}/libparanumal/make.inp
include ${OCCA_DIR}/scripts/Makefile

INCLUDES = ins.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh2D.h \
$(HDRDIR)/mesh3D.h \
$(OGSDIR)/ogs.hpp \
$(ALMONDDIR)/parAlmond.hpp \
$(ELLIPTICDIR)/elliptic.h \
$(ELLIPTICDIR)/ellipticPrecon.h \
$(ELLIPTICDIR)/ellipticMultiGrid.h \

# types of files we are going to construct rules for
.SUFFIXES: .c

# rule for .c files
CFLAGS_:=$(CFLAGS)

.c.o: $(DEPS)
	$(CXX) $(CFLAGS_) -D DHOLMES=\"$(NEKRS_INSTALL_DIR)/libparanumal\" \
	-o $*.o -c $*.c $(paths)

## We need separate DHOLMES and DINS when CFLAGS gets quoted
CFLAGS:= $(CFLAGS_) -D DHOLMES="\\\"$(NEKRS_INSTALL_DIR)/libparanumal\\\"" \
	-D DOGS="\\\"$(NEKRS_INSTALL_DIR)/gatherScatter\\\"" \
	-D DELLIPTIC="\\\"$(NEKRS_INSTALL_DIR)/elliptic\\\"" \
	-D DPARALMOND="\\\"$(NEKRS_INSTALL_DIR)/parAlmond\\\"" \

# list of objects to be compiled

# library objects
LOBJS = \
$(HDRDIR)/src/meshConnectPeriodicFaceNodes3D.o \
$(HDRDIR)/src/meshConnectPeriodicFaceNodes2D.o \
$(HDRDIR)/src/meshSetupBoxHex3D.o \
$(HDRDIR)/src/meshSetupBoxQuad2D.o \
$(HDRDIR)/src/meshApplyElementMatrix.o \
$(HDRDIR)/src/meshConnect.o \
$(HDRDIR)/src/meshConnectFaceNodes2D.o \
$(HDRDIR)/src/meshConnectFaceNodes3D.o \
$(HDRDIR)/src/meshGeometricFactorsTet3D.o \
$(HDRDIR)/src/meshGeometricFactorsHex3D.o \
$(HDRDIR)/src/meshGeometricFactorsTri2D.o \
$(HDRDIR)/src/meshGeometricFactorsTri3D.o \
$(HDRDIR)/src/meshGeometricFactorsQuad2D.o \
$(HDRDIR)/src/meshGeometricFactorsQuad3D.o \
$(HDRDIR)/src/meshGeometricPartition2D.o \
$(HDRDIR)/src/meshGeometricPartition3D.o \
$(HDRDIR)/src/meshPartitionStatistics.o \
$(HDRDIR)/src/meshHaloExchange.o \
$(HDRDIR)/src/meshHaloExtract.o \
$(HDRDIR)/src/meshHaloSetup.o \
$(HDRDIR)/src/meshLoadReferenceNodesTri2D.o \
$(HDRDIR)/src/meshLoadReferenceNodesQuad2D.o \
$(HDRDIR)/src/meshLoadReferenceNodesTet3D.o \
$(HDRDIR)/src/meshOccaSetup2D.o \
$(HDRDIR)/src/meshOccaSetup3D.o \
$(HDRDIR)/src/meshOccaSetupQuad3D.o \
$(HDRDIR)/src/meshParallelConnectNodes.o \
$(HDRDIR)/src/meshParallelConnectOpt.o \
$(HDRDIR)/src/meshParallelConsecutiveGlobalNumbering.o\
$(HDRDIR)/src/meshParallelGatherScatterSetup.o \
$(HDRDIR)/src/meshParallelReaderTri2D.o \
$(HDRDIR)/src/meshParallelReaderHex3D.o \
$(HDRDIR)/src/meshParallelReaderTri3D.o \
$(HDRDIR)/src/meshParallelReaderQuad2D.o \
$(HDRDIR)/src/meshParallelReaderQuad3D.o \
$(HDRDIR)/src/meshParallelReaderTet3D.o \
$(HDRDIR)/src/meshPhysicalNodesTri2D.o \
$(HDRDIR)/src/meshPhysicalNodesTri3D.o \
$(HDRDIR)/src/meshPhysicalNodesQuad2D.o \
$(HDRDIR)/src/meshPhysicalNodesQuad3D.o \
$(HDRDIR)/src/meshPhysicalNodesTet3D.o \
$(HDRDIR)/src/meshPlotVTU2D.o \
$(HDRDIR)/src/meshPlotVTU3D.o \
$(HDRDIR)/src/meshPrint2D.o \
$(HDRDIR)/src/meshPrint3D.o \
$(HDRDIR)/src/meshSetup.o \
$(HDRDIR)/src/meshSetupTri2D.o \
$(HDRDIR)/src/meshSetupTri3D.o \
$(HDRDIR)/src/meshSetupQuad2D.o \
$(HDRDIR)/src/meshSetupQuad3D.o \
$(HDRDIR)/src/meshSetupHex3D.o \
$(HDRDIR)/src/meshSetupTet3D.o \
$(HDRDIR)/src/meshSurfaceGeometricFactorsTri2D.o \
$(HDRDIR)/src/meshSurfaceGeometricFactorsTri3D.o \
$(HDRDIR)/src/meshSurfaceGeometricFactorsQuad2D.o \
$(HDRDIR)/src/meshSurfaceGeometricFactorsQuad3D.o \
$(HDRDIR)/src/meshSurfaceGeometricFactorsTet3D.o \
$(HDRDIR)/src/meshSurfaceGeometricFactorsHex3D.o \
$(HDRDIR)/src/meshLoadReferenceNodesHex3D.o \
$(HDRDIR)/src/meshConnectBoundary.o \
$(HDRDIR)/src/meshVTU2D.o \
$(HDRDIR)/src/meshVTU3D.o \
$(HDRDIR)/src/matrixInverse.o \
$(HDRDIR)/src/matrixConditionNumber.o \
$(HDRDIR)/src/mysort.o \
$(HDRDIR)/src/parallelSort.o\
$(HDRDIR)/src/hash.o\
$(HDRDIR)/src/setupAide.o \
$(HDRDIR)/src/readArray.o\
$(HDRDIR)/src/occaHostMallocPinned.o \
$(HDRDIR)/src/timer.o

ifeq ($(OS),Windows_NT)
    detected_OS := Windows
else
    detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

SONAME=
WHOLE=
NOWHOLE=
SHARED=
EXT=

ifeq ($(detected_OS),Linux)
	SONAME=-soname
	WHOLE=--whole-archive
	NOWHOLE=--no-whole-archive
	SHARED=-shared
	EXT=so
endif
ifeq ($(detected_OS),Darwin)
	SONAME=-install_name
	WHOLE=-all_load
	NOWHOLE=-noall_load
	SHARED=-dynamiclib -undefined dynamic_lookup
	EXT=dylib
endif

libelliptic: libogs libP libparAlmond
	$(MAKE) -C $(ELLIPTICDIR) cc="$(CC)" FC="$(FC)" CC="$(CXX)" LD="$(CXX)" \
	  CFLAGS="$(CFLAGS)" GSDIR="$(GSDIR)" LDFLAGS="$(LDFLAGS)" sharedlib

libP: libblas libogs $(LOBJS)
	cd $(HDRDIR) && $(CXX) $(SHARED) -o libP.$(EXT) $(LOBJS) -Wl,$(WHOLE) \
        -L$(BLASDIR) -lBlasLapack -Wl,$(NOWHOLE) \
	-Wl,$(SONAME),$(NEKRS_INSTALL_DIR)/libparanumal/libP.$(EXT)

libogs:
	$(MAKE) -C $(OGSDIR) CC="$(CC)" CXX="$(CXX)" LD="$(CXX)" \
	  CFLAGS="$(CFLAGS)" GSDIR="$(GSDIR)" LDFLAGS="$(LDFLAGS)" sharedlib

libblas:
	$(MAKE) -C $(BLASDIR) FC="$(FC)" FCFLAGS="$(FFLAGS)" lib

libparAlmond:
	$(MAKE) -C $(ALMONDDIR) CC="$(CC)" CXX="$(CXX)" LD="$(CXX)" \
	CFLAGS="$(CFLAGS)" GSDIR="$(GSDIR)" LDFLAGS="$(LDFLAGS)" \
	ENABLE_HYPRE="$(ENABLE_HYPRE)" lib


all: libP libelliptic

clean:
	cd $(ELLIPTICDIR); make clean; cd -
	cd ${HDRDIR}/src; rm *.o; cd -

realclean:
	cd $(ELLIPTICDIR); make realclean; cd -
	cd ${HDRDIR}/src; rm *.o; cd -
