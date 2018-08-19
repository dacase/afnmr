
#************************************************************************
#
#************************************************************************
#
include ../config.h

CONFIG_FILE=../config.h
CONFIG_COMMAND=configure
PARALLEL_TOKEN=DMPI

CUDA_TOKEN=DCUDA
CUDA_INC=-I$(CUDA_HOME)/include
CUSP_LIB=-I../cusplibrary-cuda9

# --- Begin compilation macros for CUDA CG ---
#  Caution! No space after variable value!
CG_CU_LIB=CUSP# [CUSPARSE, CUSP]
CG_METHOD=PCG# [CG, PCG]
PCG_OPT=Jacobi# [ILU0, IC02, Jacobi, Smooth]
MATRIX_FORMAT=DIA# [CSR, DIA, ELL, COO, HYB]

# Check combinations
.PHONY: .cuda_precheck .UPDATE
.cuda_precheck:
	@if [ "$(CG_CU_LIB)" = "CUSPARSE" ]; then \
		if [ "$(CG_METHOD)" = "PCG" ] && [ "$(PCG_OPT)" != "ILU0" ] && [ "$(PCG_OPT)" != "IC02" ]; then \
			echo "Error: Currently CUSPARSE only supports PCG ILU0 & IC02."; \
			exit 2; \
		fi; \
		if [ "$(MATRIX_FORMAT)" != "CSR" ]; then \
			echo "Error: Currently CUSPARSE only supports matrix format CSR."; \
			exit 2; \
		fi; \
	elif [ "$(CG_CU_LIB)" = "CUSP" ]; then \
		if [ "$(CG_METHOD)" = "PCG" ]; then \
			if [ "$(PCG_OPT)" != "Jacobi" ] && [ "$(PCG_OPT)" != "Smooth" ]; then \
				echo "Error: Currently CUSP only supports PCG Jacobi & Smooth."; \
				exit 2; \
			fi; \
			if [ "$(PCG_OPT)" = "Smooth" ] && [ "$(MATRIX_FORMAT)" = "DIA" ]; then \
				echo "Error: CUSP Smooth PCG does not support matrix format DIA."; \
				exit 2; \
			fi; \
		fi; \
	else \
		echo "Error: CG_CU_LIB is not properly set."; \
		exit 2; \
	fi; \
	if [ "$(CG_METHOD)" = "CG" ]; then \
		echo "You have selected 'CG' method, now ignoring PCG option $(PCG_OPT)."; \
	fi
# --- End compilation macros for CUDA CG ---

# ---- Begin file definitions ----
# PBSA Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

SHELL=/bin/sh
VPATH=

.SUFFIXES:
.SUFFIXES: .c .F90 .o

OBJ=\
	dsvdc.o  pbsa.o pb_read.o pb_write.o getcoor.o \
	runmd.o runmin.o force.o pbsa_lib.o pb_init.o sa_driver.o density.o np_force.o pb_force.o \
        pb_p3m.o pb_fdfrc.o pb_crgview.o pb_fddrv.o pb_iimdrv.o pb_anaiim.o pb_nlsolver.o \
	pb_bldsys.o pb_mpfrc.o pb_direct.o pb_list.o timer.o \
	egb.o ene.o locmem.o myopen.o rdparm.o pb_lsolver.o \
	decomp.o rgroup.o rfree.o debug.o \
	irre32.o  phidev.o  project.o  qld.o \
	irre31.o  miniop.o \
	transf.o curv.o GrToPr.o jumps.o \
	regular.o gmres.o \
	daxpy.o  ddot.o    dheqr.o  dpigmr.o  dscal.o   dslui2.o  dxlcal.o \
	dchkw.o  dgmres.o  dnrm2.o  drlcal.o  dsilus.o  dslui.o   isdgmr.o \
	dcopy.o  dhels.o   dorth.o  ds2y.o    dslugm.o  dsmv.o    qs2i1d.o \
	d1mach.o  fdump.o   i1mach.o  j4save.o  xercnt.o  xerhlt.o  xermsg.o  \
        xerprn.o  xersve.o  xgetua.o \
	bicg.o dbcg.o dslui4.o dsluti.o dsmtv.o dslubc.o isdbcg.o \
        pb_chkfrc.o svdcmp.o svbksb.o pythag.o\
        pb_augdrv.o interpX.o gen_dx_file.o \
        parms.o memory_module.o rdpqr.o

CUDAOBJ=\
	cuda_mg_wrapper.o cuda_cg_wrapper.o band2csr.o

# GM: when updating this variable, please don't omit to update makedepend
LIBPBSAOBJS= \
        interpX.LIBPBSA.o parms.LIBPBSA.o pb_fdfrc.LIBPBSA.o pb_force.LIBPBSA.o pb_p3m.LIBPBSA.o \
        pb_read.LIBPBSA.o sa_driver.LIBPBSA.o \
        curv.LIBPBSA.o decomp.LIBPBSA.o density.LIBPBSA.o dsvdc.LIBPBSA.o gen_dx_file.LIBPBSA.o \
        GrToPr.LIBPBSA.o interface.o memory_module.LIBPBSA.o miniop.LIBPBSA.o myopen.LIBPBSA.o \
        np_force.LIBPBSA.o pb_bldsys.LIBPBSA.o pb_crgview.LIBPBSA.o pb_direct.LIBPBSA.o pb_fddrv.LIBPBSA.o \
        pb_init.LIBPBSA.o pb_list.LIBPBSA.o pb_lsolver.LIBPBSA.o pb_mpfrc.LIBPBSA.o pb_nlsolver.LIBPBSA.o  \
        pbsaapi.LIBPBSA.o pbsa_lib.LIBPBSA.o pb_write.LIBPBSA.o phidev.LIBPBSA.o project.LIBPBSA.o \
        pythag.LIBPBSA.o rfree.LIBPBSA.o rgroup.LIBPBSA.o svbksb.LIBPBSA.o svdcmp.LIBPBSA.o \
        timer.LIBPBSA.o transf.LIBPBSA.o

# ---- End file definitions ----

# ---- Begin installation options ----
install: install.$(INSTALLTYPE)

install.serial: serial
	/bin/mv pbsa$(SFX) $(BINDIR)

install.cuda_serial: cuda
	/bin/mv pbsa.cuda$(SFX) $(BINDIR)

serial: configured_serial pbsa$(SFX)

cuda: configured_cuda pbsa.cuda$(SFX)

uninstall:
	-/bin/rm -f $(BINDIR)/pbsa$(SFX) pbsa.cuda$(SFX) $(LIBDIR)/libpbsa.a

libinstall: libinstall.$(INSTALLTYPE)

libinstall.openmp: libinstall.serial

libinstall.serial: libpbsa.a
	/bin/mv libpbsa.a $(LIBDIR)

libinstall.parallel: libpbsa.a
	/bin/mv libpbsa.a $(LIBDIR)
# ---- End installation options ----

# ---- Begin configuration check ----
configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
		echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE) !" ;\
		exit 2 ;\
		fi ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial !" ;\
		echo "  Rerun $(CONFIG_COMMAND) and do NOT specify -mpi." ;\
		exit 2 ;\
	fi ;\
	)

configured_cuda: configured_serial .cuda_precheck
	        @(grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null || \
	                { echo "Error: $(CONFIG_FILE) is not of type cuda !" ;\
	                        echo " Rerun $(CONFIG_COMMAND) and specify -cuda." ;\
	                        exit 2 ;\
	                };\
	        if [ -z $(CUDA_HOME) ]; then \
	                echo "Error: CUDA_HOME is not set. This must point to your NVIDIA tools installation"; \
	                exit 2; \
	                fi; \
	        )
# ---- End configuration check ----

# ---- Begin cuda object files ----
.UPDATE:
	$(info Always rebuild [ pb_lsolver.o pb_fddrv.o [ $(CUDAOBJ) ]] to make predefined macros effective)

$(CUDAOBJ) pb_fddrv.o pb_lsolver.o: .UPDATE

cuda_mg_wrapper.o: cuda_mg_wrapper.cu
	$(NVCC) -ccbin $(CXX) -o $@ -c $<
cuda_cg_wrapper.o: cuda_cg_wrapper.cu
	$(NVCC) -ccbin $(CXX) $(CUSP_LIB) -o $@ -c $< -D$(CG_CU_LIB) -D$(CG_METHOD) -D$(PCG_OPT) -D$(MATRIX_FORMAT)
band2csr.o: band2csr.cpp
	$(CXX) $(CXXOPTFLAGS) -o $@ -c $<
# ---- End buda object files ----

# ---- Begin build pbsa options ----
pbsa$(SFX): $(OBJ) syslib netlib c9x-complex configured_serial
	$(FC) $(PBSAFLAG) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o pbsa$(SFX) $(OBJ) \
	        ../lib/nxtsec.o ../lib/random.o \
	        -L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

pbsa.cuda$(SFX): $(OBJ) syslib netlib c9x-complex configured_serial $(CUDAOBJ)
	$(FC) $(PBSAFLAG) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o pbsa.cuda$(SFX) \
		$(OBJ) $(CUDAOBJ) $(PBSA_CU_LIBS) $(CUDA_INC) \
		../lib/nxtsec.o ../lib/random.o \
		-L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)
# ---- End build pbsa options ----

# ---- Begin build other options ----

simple:  $(BINDIR)/simplepbsa$(SFX)

$(BINDIR)/simplepbsa$(SFX): simplepbsa.o gopt.o sfflib syslib netlib c9x-complex
	$(CC) -o $(BINDIR)/simplepbsa$(SFX) simplepbsa.o gopt.o -L$(LIBDIR) $(FLIBS) 

libpbsa.a: $(LIBPBSAOBJS) syslib
	-rm -f libpbsa.a
	$(AR) libpbsa.a $(LIBPBSAOBJS) ../lib/random.o
	$(RANLIB) libpbsa.a

# ---- End build other options ----

syslib:
	cd ../lib && $(MAKE) nxtsec.o random.o

sfflib:
	cd ../sff && $(MAKE) install

netlib:
	cd ../lapack && $(MAKE) $(LAPACK)
	cd ../blas && $(MAKE) $(BLAS)
	cd ../arpack && $(MAKE) install

c9x-complex:
	@if test $(C9XCOMPLEX) != "skip"; then \
	   cd ../c9x-complex && $(MAKE) libmc.a; \
	 fi

# --- Predefine -DCUDA for pb_lsolver.F90/pb_fddrv.F90 ---
ifeq ($(INSTALLTYPE), cuda_serial)
    CU_define=-$(CUDA_TOKEN)
else
    CU_define=
endif

pb_init.o: pb_init.F90
	$(FC) $(PBSAFLAG) $(FPPFLAGS) -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

pb_init.LIBPBSA.o: pb_init.F90
	$(FC) $(PBSAFLAG) $(FPPFLAGS) -DLIBPBSA -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

np_force.LIBPBSA.o: np_force.F90
	$(FC) $(PBSAFLAG) $(FPPFLAGS) -DLIBPBSA -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

%.LIBPBSA.o: %.F90
	$(FC) $(PBSAFLAG) $(FPPFLAGS) -DLIBPBSA -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

interface.o: interface.c
	$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERFFLAGS) -o $@ $<

.F90.o:   $<
	$(FC) $(PBSAFLAG) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $< $(CU_define)

.c.o: $<
	$(CC) $(CFLAGS) $(COPTFLAGS) $(AMBERFLAGS) -c -o $@ $<

depend::
	./makedepend > depend

clean:
	/bin/rm -f pbsa$(SFX) libpbsa.a
	/bin/rm -f *.o
	/bin/rm -f *.mod *__genmod.F90
	/bin/rm -f *.d

distclean: clean



#----------HEADER DEPENDENCIES (currently from makedepend)
include depend
# DO NOT DELETE
