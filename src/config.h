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

FLIBS=-lsff -llapack -lblas $(LIBGFORTRAN) -lgfortran

CC=gcc
CFLAGS=-g
COPTFLAGS=-O3 -mtune=native

CXX=g++
CXXFLAGS=-g
CXXOPTFLAGS=-O3

FC=gfortran
FFLAGS=-I$(INCDIR) -g
FOPTFLAGS=-O3 -mtune=native

AR=ar rv
RANLIB=ranlib
FLEX=flex
YACC=bison -y
VB=@
