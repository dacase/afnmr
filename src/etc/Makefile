.SUFFIXES : .F90 .o

AOBJ=   ambpdb.o forFIRSTtype.o forFIRSTbond.o forFIRSThbond.o \
        forFIRSTteth.o forFIRSToutp.o forFIRSTsugar.o nxtsec.o \
        amopen.o mexit.o

include ../config.h

PROGS= $(BINDIR)/ambpdb $(BINDIR)/add_pdb

$(BINDIR)/ambpdb: $(AOBJ)
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(AOBJ)

$(BINDIR)/add_pdb: add_pdb.o nxtsec.o
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ add_pdb.o nxtsec.o

install:	$(PROGS)

clean:
	/bin/rm -f *.o *.mod

#  default rules for Fortran and C compilation:

.F90.o: $<
	@echo "[ETC]  FC $<"
	$(VB)$(FC) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

