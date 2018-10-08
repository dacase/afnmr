include ../config.h

.c.o:
	$(CC) -c $(COPTFLAGS) $(CFLAGS) -o $@ $<

OBJS = conjgrad.o lmodC.o memutil.o nblist.o newton.o nmode.o wallclock.o \
	   prm.o rand2.o sasad.o sff.o xminC.o ptrajmask.o

install: $(LIBDIR)/libsff.a
	cp sff.h $(INCDIR)

$(LIBDIR)/libsff.a: $(OBJS)
	-rm -f $(LIBDIR)/libsff.a
	$(AR) $(LIBDIR)/libsff.a $(OBJS)	
	$(RANLIB) $(LIBDIR)/libsff.a

uninstall:
	-rm -f $(LIBDIR)/libsff.a
	-rm -f $(LIBDIR)/libsff_mpi.a
	-rm -f $(INCDIR)/sff.h

clean:
	-/bin/rm -f $(OBJS) lex.mm_options.c

lex.mm_options.c: mm_options.l
	$(LEX) -t mm_options.l | sed -e 's/stdout/nabout/' > lex.mm_options.c

depend: $(OBJS:.o=.c) lex.mm_options.c 
	$(CC) -MM $^ > $@

include depend

