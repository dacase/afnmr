include ../config.h

install:  $(BINDIR)/lsq $(BINDIR)/ms $(BINDIR)/afnmr.x  \
          $(BINDIR)/gridprune $(BINDIR)/pair $(BINDIR)/txmin
	/bin/cp afnmr getshifts-* multiafnmr.sh $(BINDIR)

$(BINDIR)/lsq: lsq.o conjgrad.o rand2.o memutil.o
	$(CC) -o $@ lsq.o conjgrad.o rand2.o memutil.o $(LM)

$(BINDIR)/ms: ms.f
	$(FC) $(FFLAGS) -o $@ $<

$(BINDIR)/pair: pair.o number.o pof.o strings.o checkio.o 
	$(CC) $(CFLAGS) -o $@ pair.o number.o pof.o strings.o \
            checkio.o $(LM)

$(BINDIR)/afnmr.x: afnmr.F90
	$(FC) $(FFLAGS) -o $@ $<

$(BINDIR)/gridprune: gridprune.F90
	$(FC) $(FFLAGS) $(FOPTFLAGS) -o $@ $< -L$(LIBDIR) $(FLIBSF)

$(BINDIR)/txmin: txmin.c
	$(CC) $(CFLAGS) -o $@ $< -L$(LIBDIR) $(FLIBS)

clean:
	/bin/rm -f *.o *.mod 
