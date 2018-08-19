.f.o: $<
	gfortran -c -O2 -o $@ $<

OBJECTS = reg_kspace.o

LIB = reg_recip.a

$(LIB) : $(OBJECTS)
	ar rlvs $(LIB) $(OBJECTS)
