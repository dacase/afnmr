#  Simple hard-wired config.h file for gnu compilers

#set AFNMRHOME here, or use an environment variable:

BINDIR=$(AFNMRHOME)/bin
LIBDIR=$(AFNMRHOME)/lib
INCDIR=$(AFNMRHOME)/include
DATDIR=$(AFNMRHOME)/dat
LOGDIR=$(AFNMRHOME)/logs

OS=$(shell uname -s)
ifeq "$(OS)" "Darwin"
   # MacOSX rules:
   SHARED_SUFFIX=.dylib
   MAKE_SHARED=-dynamiclib
   LIBGFORTRAN=$(shell $(AFNMRHOME)/src/libgfortran.sh)
   LDFLAGS=
   LM=
else
   # Linux rules:
   SHARED_SUFFIX=.so
   MAKE_SHARED=-shared
   LDFLAGS=
   LM=-lm
endif

ifeq "$(MKLROOT)" ""
   FLIBS=-lsff -llapack -lblas $(LIBGFORTRAN) -lgfortran
   RISM=skip
   RISMSFF=""
   BLAS=install
   LAPACK=install
else
   FLIBS=-lsff -llapack -lblas -lrism -lmpi
   RISM=install
   BLAS=skip
   LAPACK=skip
   RISMSFF=-DRISMSFF
endif

CC=icc
CFLAGS=
COPTFLAGS=-O3 

CXX=icpc
CXXFLAGS=
CXXOPTFLAGS=-O3

FC=ifort
ifeq "$(MKLROOT)" ""
   FFLAGS=-I$(INCDIR)
else
   FFLAGS=-I$(INCDIR) -I$(MKLROOT)/include -I$(MKLROOT)/include/fftw
endif
FOPTFLAGS=-O3 -mtune=native

AR=ar rv
RANLIB=ranlib
FLEX=flex
YACC=bison -y
VB=@
