#
#************************************************************************
#                             AMBER                                    **
#                                                                      **
#                Copyright (c) 1986, 1991, 1995                        **
#             Regents of the University of California                  **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#

.SUFFIXES : .c .F .F90 .o
SHELL=/bin/sh

AOBJ=   ambpdb.o forFIRSTtype.o forFIRSTbond.o forFIRSThbond.o \
        forFIRSTteth.o forFIRSToutp.o forFIRSTsugar.o nxtsec.o \
        amopen.o mexit.o

ROBJ=   resp.o lapack.o amopen.o mexit.o

include ../config.h

PROGS= $(BINDIR)/elsize$(SFX) $(BINDIR)/resp$(SFX) $(BINDIR)/ambpdb$(SFX)

PERLPROGS= process_mdout.perl process_minout.perl mdout2pymbar.pl

PYTHONPROGS=cpinutil.py ceinutil.py mdout_analyzer.py finddgref.py fitpkaeo.py fixremdcouts.py genremdinputs.py softcore_setup.py charmmlipid2amber.py

SHELLPROGS= sgldinfo.sh sgldwt.sh

$(BINDIR)/ambpdb$(SFX): 	$(AOBJ)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $(BINDIR)/ambpdb$(SFX) \
		$(AOBJ) $(LDFLAGS) $(AMBERLDFLAGS)

$(BINDIR)/add_pdb$(SFX): add_pdb.o nxtsec.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $(BINDIR)/add_pdb$(SFX) \
		add_pdb.o nxtsec.o $(LDFLAGS) $(AMBERLDFLAGS)

$(BINDIR)/elsize$(SFX): elsize.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) -o $(BINDIR)/elsize$(SFX) elsize.o $(LM)

$(BINDIR)/resp$(SFX):	$(ROBJ)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	-o $(BINDIR)/resp$(SFX) $(ROBJ)

install:	$(PROGS)

parallel:
	$(FC) $(FPPFLAGS) -o numprocs numprocs.F90
	/bin/cp numprocs ../../test

nxtsec.o:  ../sander/nxtsec.F90
	$(FC) $(FPPFLAGS) -c $(FREEFORMAT_FLAG) $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

amopen.o:  ../sander/amopen.F90
	$(FC) $(FPPFLAGS) -c $(FREEFORMAT_FLAG) $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

mexit.o:  ../sander/mexit.F90
	$(FC) $(FPPFLAGS) -c $(FREEFORMAT_FLAG) $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

python:
	(if [ "$(SKIP_PYTHON)" = "no" ]; then \
		$(PYTHON) setup.py install -f $(PYTHON_INSTALL) --install-scripts=$(BINDIR) ;\
	fi;\
	)

clean:
	/bin/rm -f *.o 90 *.mod  $(PROGS) molsurf.c
	/bin/rm -rf build

uninstall:
	-(for file in $(PROGS) $(PERLPROGS) $(PYTHONPROGS) $(SHELLPROGS) ; do \
		/bin/rm -fr $(BINDIR)/$${file} ; \
	done)

#  default rules for Fortran and C compilation:

.F90.o: $<
	$(FC) $(FPPFLAGS) -c $(FREEFORMAT_FLAG) $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

# lapack.F still uses fixed-format F77
.F.o:   $<
	$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

.c.o:
	$(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

include depend
