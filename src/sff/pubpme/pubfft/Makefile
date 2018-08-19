.f.o:  $<
	gfortran -c -O2 -o $@ $<

.c.o: $<
	gcc -c -O2 -o $@ $<

OBJECTS = cfftb.o cfftf1.o passb.o passb4.o passf2.o\
passf5.o cfftb1.o cffti.o passb2.o passb5.o passf3.o\
cfftf.o cffti1.o  passb3.o  passf.o passf4.o pub3dfft.o
LIB = libpubfft.a

$(LIB) : $(OBJECTS)
	ar rlvs $(LIB) $(OBJECTS)

clean:
	/bin/rm -f $(OBJECTS) $(LIB)
