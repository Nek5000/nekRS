FC = gfortran
AR = ar

SDIR = ./
INCLUDES = -I.
OPTFLAGS = -O3
ARFLAGS = -cr

.SUFFIXES: .f

SRC = $(wildcard $(SDIR)/*f)
OBJS = $(subst $(SDIR)/,$(SDIR),$(SRC:.f=.o))

FCFLAGS = $(FCOPTIONS) $(OPTFLAGS)

.f.o:
	$(FC) $(FCFLAGS) -c $*.f 

lib: $(OBJS)
	$(AR) $(ARFLAGS) libBlasLapack.a *.o

clean:
	rm -f *.o *.a
