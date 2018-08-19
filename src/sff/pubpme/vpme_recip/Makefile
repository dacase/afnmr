
.f.o: $<
	cpp -traditional -DPUBFFT $< > _$<
	gfortran -c -O2 -o $@ _$<

#.c.o: $<
#	gcc -c -O2 -o $@ $<

OBJECTS = charge_grid.o fftcalls.o bspline.o pmesh_kspace.o

LIB = pme_recip.a

$(LIB) : $(OBJECTS)
	ar rlvs $(LIB) $(OBJECTS)

clean:
	/bin/rm -f $(OBJECTS) $(LIB)
