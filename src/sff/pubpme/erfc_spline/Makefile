FFLAGS = -O -mips2
#FFLAGS = -O 

OBJECTS = erfc_locmem.o erfc_spline.o erfcfun.o 

LIB = ../../lib/erf_table.a

$(LIB) : $(OBJECTS)
	ar rlvs $(LIB) $(OBJECTS)
