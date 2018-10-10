include ../config.h

.SUFFIXES:
.SUFFIXES:  .c .o

.c.o:
	$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $<

CIFOBS=	\
	lex.cifp.o \
	cifp.tab.o \
	cifparse.o

$(LIBDIR)/libcifparse.a: $(CIFOBS)
	/bin/rm -f $(LIBDIR)/libcifparse.a
	$(AR) $(LIBDIR)/libcifparse.a $(CIFOBS)
	-$(RANLIB) $(LIBDIR)/libcifparse.a

cifp.tab.h: cifp.tab.c

cifp.tab.c:  cifparse.y lex.cifp.c
	$(YACC) -v -d cifparse.y

cifp.tab.o: cifp.tab.c
	$(CC) -c -D$(LEX) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) cifp.tab.c

lex.cifp.c: cifparse.l 
	$(LEX) cifparse.l

install:  $(LIBDIR)/libcifparse.a

clean:
	/bin/rm -f *.o  cifp.output
	/bin/rm -f lex.cifp.c cifp.tab.c cifp.tab.h

uninstall: clean
	/bin/rm -f $(LIBDIR)/libcifparse.a

cif2cif:  cif2cif.c $(LIBDIR)/libcifparse.a
	$(CC) -o cif2cif cif2cif.c $(LIBDIR)/libcifparse.a

nef2nef:  nef2nef.c $(LIBDIR)/libcifparse.a
	$(CC) -o nef2nef nef2nef.c $(LIBDIR)/libcifparse.a

lex.cifp.o:	\
	cifp.tab.h
