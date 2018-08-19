include ../config.h

.c.o:
	$(CC) -c -Dflex $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(RISMSFF) $(NETCDFINC) -o $@ $<

.cpp.o:
	$(CXX) -c $(CXXOPTFLAGS) -o $@ $<
	
OBJS = conjgrad.o lmodC.o memutil.o nblist.o newton.o nmode.o \
	   prm.o rand2.o sasad.o sff.o xminC.o ptrajmask.o AmberNetcdf_dummy.o

install: libsff 

libsff: $(OBJS)
	-rm -f $(LIBDIR)/libsff.a
	$(AR) $(LIBDIR)/libsff.a $(OBJS)	
	$(RANLIB) $(LIBDIR)/libsff.a
	cp sff.h AmberNetcdf.h $(INCDIR)

libsff_mpi: $(OBJS)
	-rm -f $(LIBDIR)/libsff_mpi.a
	$(AR) $(LIBDIR)/libsff_mpi.a $(OBJS)	
	$(RANLIB) $(LIBDIR)/libsff_mpi.a
	cp sff.h AmberNetcdf.h $(INCDIR)

uninstall:
	-rm -f $(LIBDIR)/libsff.a
	-rm -f $(LIBDIR)/libsff_mpi.a
	-rm -f $(INCDIR)/sff.h $(INCDIR)/AmberNetcdf.h

clean:
	-/bin/rm -f $(OBJS) lex.mm_options.c

lex.mm_options.c: mm_options.l
	$(LEX) -t mm_options.l | sed -e 's/yy/mmo/g' -e 's/stdout/nabout/' > lex.mm_options.c

#  dependencies (hand-created for now):

binpos.o: memutil.h

conjgrad.o: memutil.h

lmodC.o: sff.h

memutil.o: sff.h

nblist.o: nblist.c
	$(CC) -c -Dflex $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(RISMSFF) \
	    $(NETCDFINC) $(SFF_INTEL_BUG361) -o $@ $<

newton.o: sff.h memutil.h

nmode.o: sff.h memutil.h

prm.o:   sff.h

ptrajmask.o:  ptrajmask.h

rand2.c:  sff.h

sasad.o: sff.h

sff.o: sff.h memutil.h timer.h\
       debug.h eff_box.c dssp.c eff.c \
       gbsa.c egbr6.c sff2.c mask.c rattle.c lex.mm_options.c \
       hcp.c hcp_gb.c timer.c AmberNetcdf.h

xminC.o:  sff.h

AmberNetcdf.o:  AmberNetcdf.h

hcp_getpdb$(SFX): hcp_getpdb.c
	$(CC) -o $(BINDIR)/hcp_getpdb$(SFX) hcp_getpdb.c

.PHONY: ../rism/amber_rism_interface.NAB.o
../rism/amber_rism_interface.NAB.o:
	cd ../rism && \
	$(MAKE) amber_rism_interface.NAB.o


