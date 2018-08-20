#  Simple hard-wired config.h file for gnu compilers

#set AFNMRHOME here, or use an environment variable:

BINDIR=$(AFNMRHOME)/bin
LIBDIR=$(AFNMRHOME)/lib
INCDIR=$(AFNMRHOME)/include
DATDIR=$(AFNMRHOME)/dat
LOGDIR=$(AFNMRHOME)/logs

FLIBS=-lsff -llapack -lblas -lgfortran

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
YACC=yacc
