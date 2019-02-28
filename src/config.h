#  Simple hard-wired config.h file for gnu compilers

#set AFNMRHOME here, or use an environment variable:

BINDIR=$(AFNMRHOME)/bin
LIBDIR=$(AFNMRHOME)/lib
INCDIR=$(AFNMRHOME)/include
DATDIR=$(AFNMRHOME)/dat
LOGDIR=$(AFNMRHOME)/logs

OS=$(shell uname -s)
ifeq "$(OS)" "Darwin"
   # MacOSX rules for making shared objects
   SHARED_SUFFIX=.dylib
   MAKE_SHARED=-dynamiclib
   LIBGFORTRAN=$(shell $(AFNMRHOME)/src/libgfortran.sh)
else
   # Linux rules for shared libraries:
   SHARED_SUFFIX=.so
   MAKE_SHARED=-shared
endif

FLIBS=-lsff -llapack -lblas $(LIBGFORTRAN) -lgfortran
LM=-lm

CC=gcc
CFLAGS=
COPTFLAGS=-O3 -mtune=native

CXX=g++
CXXFLAGS=
CXXOPTFLAGS=-O3

FC=gfortran
FFLAGS=-I$(INCDIR)
FOPTFLAGS=-O3 -mtune=native

AR=ar rv
RANLIB=ranlib
LEX=flex
YACC=bison -y
VB=@
