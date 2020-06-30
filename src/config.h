#  Simple hard-wired config.h file for gnu compilers
#    If MKLROOT is defined, use MKL; otherwise not.  Note that (for now)
#    MKL is required for the 3D-RISM option.

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
   FFLAGS=-I$(INCDIR)
else
   FLIBS=-lsff -lrism $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_sequential.a $(MKLROOT)/lib/libmkl_core.a -lpthread -ldl $(LIBGFORTRAN) -lgfortran
   RISM=install
   RISMSFF=-DRISMSFF
   BLAS=skip
   LAPACK=skip
   FFLAGS=-I$(INCDIR) -I$(MKLROOT)/include -I$(MKLROOT)/include/fftw
endif

CC=gcc
CFLAGS=
COPTFLAGS=-O3 -mtune=native

CXX=g++
CXXFLAGS=
CXXOPTFLAGS=-O3 -mtune=native

FC=gfortran
FOPTFLAGS=-O3 -mtune=native

AR=ar rv
RANLIB=ranlib
FLEX=flex
YACC=bison -y
VB=@
