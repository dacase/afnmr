Simple hard-wired config.h file for gnu compilers

AMBERHOME=/home/case/afnmr     # edit to match this installation!
BINDIR=$(AFNMRHOME)/bin
LIBDIR=$(AFNMRHOME)/lib
INCDIR=$(AFNMRHOME)/include
DATDIR=$(AFNMRHOME)/dat
LOGDIR=$(AFNMRHOME)/logs

FLIBS=-lsff -larpack -llapack -lblas

CC=gcc
CFLAGS=
COPTFLAGS=-O3 -mtune=native

CXX=g++
CXXFLAGS=
CXXOPTFLAGS=-O3

FC=gfortran
FFLAGS=-I$(INCDIR)
FOPTFLAGS=-O3 -mtune=native
