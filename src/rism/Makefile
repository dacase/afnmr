SHELL=/bin/sh
include ../config.h

# trick to use $(AMBER_PREFIX)/bin/nab
export AMBERHOME=$(AMBER_PREFIX)

# rism1d Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

# ------- rism1d information: ----------------------------------------------

OBJ=rism1d.o
SRC=rism1d.F90

.PHONY:  yes no \
	librism librism_mpi rism1d$(SFX) \
	rism3d.snglpnt$(SFX) rism3d.snglpnt.MPI$(SFX) \
	install install_mpi uninstall clean


rism1d$(SFX): $(OBJ) librism
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o rism1d$(SFX) $(OBJ) \
	        -L$(LIBDIR) -lrism $(FLIBS_FFTW3) \
	        $(FLIBSF)  $(LDFLAGS) $(AMBERLDFLAGS)  $(LIBDIR)/libxblas-amb.a

# ------- rism3d information: ----------------------------------------------

LIBOBJ = rism1d_c.o solvmdl_c.o rism1d_potential_c.o\
	rism1d_closure_c.o rism1d_kh_c.o rism1d_hnc_c.o rism1d_py_c.o \
	rism1d_mv0_c.o rism1d_psen_c.o rism1d_nub_c.o rism1d_polyt_c.o \
        rism1d_devb_c.o \
	quaternion.o array_util.o rism_util.o rism_report_c.o \
	rism3d_grid_c.o rism3d_closure_c.o \
	rism3d_kh_c.o rism3d_hnc_c.o rism3d_psen_c.o \
	rism3d_c.o rism3d_potential_c.o rism3d_restart.o \
	rism3d_xyzv.o rism3d_opendx.o rism3d_ccp4.o \
	rism3d_solvent_c.o rism3d_solute_c.o \
	pubfft.o rism3d_fft.o rism_parm.o rism_io.o \
	mdiis_orig_c.o mdiis_blas_c.o mdiis_blas2_c.o mdiis_c.o \
	fce_c.o erfcfun.o safemem.o blend.o rism_timer_c.o constants.o \
	getopts_c.o bspline.o  fftw3.o #pmepot/pmepot.o #rism3d_debug_c.o

LIBOBJ_MPI = quaternion.o array_util.o rism_util.o rism_report_c.o \
	rism3d_grid_c.o rism3d_closure_c.o \
	rism3d_kh_c.o rism3d_hnc_c.o rism3d_psen_c.o \
	rism3d_c.o rism3d_potential_c.o rism3d_restart.o \
	rism3d_xyzv.o rism3d_opendx.o rism3d_ccp4.o \
	rism3d_solvent_c.o rism3d_solute_c.o \
	pubfft.o rism3d_fft.o rism_parm.o rism_io.o \
	mdiis_orig_c.o mdiis_blas_c.o mdiis_blas2_c.o mdiis_c.o \
	fce_c.o erfcfun.o safemem.o blend.o rism_timer_c.o constants.o \
	bspline.o fftw3.o #pmepot/pmepot.o #rism3d_debug_c.o


librism: $(LIBOBJ)
	$(AR) $(LIBDIR)/$@.a $(LIBOBJ)
	$(RANLIB) $(LIBDIR)/$@.a

librism_mpi: $(LIBOBJ_MPI)
	$(AR) $(LIBDIR)/$@.a $(LIBOBJ_MPI)
	$(RANLIB) $(LIBDIR)/$@.a

# ------ volslice information ------------------------------------------
VSSRC=volslice.F90
VSOBJ=volslice.o
volslice$(SFX): $(VSOBJ) librism
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $(BINDIR)/volslice$(SFX) $(VSOBJ) \
	        -lrism $(FLIBSF) $(FLIBS_FFTW3) \
	        $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS) $(LIBDIR)/libxblas-amb.a -L$(LIBDIR)

# ------- orientational averaging information: ----------------------------------------------
ORAVE_OBJ=rism3d.orave.o getopts_c.o array_util.o
#create an include file that has all of the quadrature files in it
orave.h:

rism3d.orave$(SFX): $(ORAVE_OBJ) librism orave.h
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o rism3d.orave$(SFX) $(ORAVE_OBJ) \
	        -L$(LIBDIR) -lrism  $(FLIBS_FFTW3) \
	        $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

# ------- thermo information: ----------------------------------------------
THERMOOBJ=rism3d.thermo.o\
	getopts_c.o\
	array_util.o\
	rism_parm.o
THERMOSRC=rism3d.thermo.F90 getopts_c.F90 array_util.F90

rism3d.thermo$(SFX): $(THERMOOBJ) librism
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o rism3d.thermo$(SFX) \
	      -L$(LIBDIR) $(THERMOOBJ) -lrism  $(FLIBS_FFTW3) \
	      $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)
# ------ Single-Point information ------------------------------------------
SPSRC=rism3d.snglpnt.c
rism3d.snglpnt$(SFX): $(SPSRC)
	$(CC) $(CFLAGS) -I$(INCDIR) $(SPSRC) -L$(LIBDIR) $(FLIBS) $(LM) \
       -o $(BINDIR)/rism3d.snglpnt$(SFX)

rism3d.snglpnt.MPI$(SFX): $(SPSRC)
	$(CC) $(CFLAGS)  -I$(INCDIR) $(SPSRC) -L$(LIBDIR) $(FLIBS) $(LM) \
         -o $(BINDIR)/rism3d.snglpnt.MPI$(SFX)
# --------------------------------------------------------------------------

yes: install

no:

install: librism rism1d$(SFX) rism3d.thermo$(SFX) rism3d.orave$(SFX) \
	volslice$(SFX)
	/bin/mv rism1d$(SFX) rism3d.thermo$(SFX)  rism3d.orave$(SFX) \
	$(BINDIR)

install_rism3d.snglpnt: librism rism3d.snglpnt$(SFX)
install_rism3d.snglpnt.MPI: install_mpi rism3d.snglpnt.MPI$(SFX)

install_mpi: librism_mpi

uninstall:
	-rm -f $(BINDIR)/rism1d$(SFX) \
	$(BINDIR)/rism1d.orig$(SFX) \
	$(BINDIR)/rism3d.thermo$(SFX) \
	$(BINDIR)/rism3d.orave$(SFX) \
	$(BINDIR)/rism3d.snglpnt$(SFX) \
	$(BINDIR)/rism3d.snglpnt.MPI$(SFX) \
	$(BINDIR)/volslice$(SFX) \
	$(LIBDIR)/librism.a \
	$(LIBDIR)/librism_mpi.a

clean:
	-/bin/rm -f *.o  pmepot/*.o *.mod rism1d$(SFX) rism3d.thermo$(SFX) \
	rism3d.snglpnt$(SFX) rism3d.snglpnt.MPI$(SFX)

doc:
	doxygen rism.Doxyfile

depend::
	./makedepend > depend

include depend

# pmepot.o: pmepot/pmepot.cpp pmepot/pmepot.hpp
# 	$(CXX) -c -std=c++11 -O3 -I$(INCDIR) pmepot/pmepot.cpp
