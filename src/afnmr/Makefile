include ../config.h

.SUFFIXES: .o .c 

install:  $(BINDIR)/lsq $(BINDIR)/ms $(BINDIR)/afnmr.x  \
          $(BINDIR)/gridprune $(BINDIR)/pair $(BINDIR)/txmin
	$(VB)/bin/cp afnmr getshifts-* make-shiftsrdb gettensor make-tensorrdb \
          multiafnmr.sh log_to_hlypqr get_hlypqr $(BINDIR)

$(BINDIR)/lsq: lsq.o conjgrad.o rand2.o memutil.o
	@echo "[AFNMR]  CC $@"
	$(VB)$(CC) $(LDFLAGS) -o $@ lsq.o conjgrad.o rand2.o memutil.o $(LM)

$(BINDIR)/ms: ms.f
	@echo "[AFNMR]  FC $@"
	$(VB)$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $<

$(BINDIR)/pair: pair.o pof.o strings.o
	@echo "[AFNMR]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ pair.o pof.o strings.o \
            $(LM)

$(BINDIR)/afnmr.x: afnmr.F90
	@echo "[AFNMR]  FC $@"
	$(VB)$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $<

$(BINDIR)/gridprune: gridprune.F90
	@echo "[AFNMR]  FC $@"
	$(VB)$(FC) $(FFLAGS) $(FOPTFLAGS) $(LDFLAGS) -o $@ $< -L$(LIBDIR) $(FLIBS)

$(BINDIR)/txmin: txmin.c
	@echo "[AFNMR]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< -L$(LIBDIR) $(FLIBS) $(LM)

.c.o:
	@echo "[AFNMR]  CC $<"
	$(VB)$(CC) -c $(CFLAGS) -o $@ $<

clean:
	/bin/rm -f *.o *.mod 
