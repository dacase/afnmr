.f.o:  $<
	gfortran -c -O0 -g -o $@ $<

../exe/filter:	filter_pdb.o
	gfortran -o ../exe/filter filter_pdb.o

../exe/exact:	exact.o utility.o erfcfun.o ../reg_recip/reg_recip.a
	gfortran -o ../exe/exact exact.o utility.o ../reg_recip/reg_recip.a \
erfcfun.o -lm

../exe/chk_exact: chk_exactness.o utility.o
	gfortran -o ../exe/chk_exact chk_exactness.o utility.o erfcfun.o

../exe/reg_recip:	dr_reg_recip.o utility.o ../reg_recip/reg_recip.a
	gfortran -o ../exe/reg_recip dr_reg_recip.o utility.o \
../reg_recip/reg_recip.a erfcfun.o -lm

../exe/pme_recip:	dr_pme_recip.o utility.o ../pme_recip/pme_recip.a ../pubfft/libpubfft.a
	gfortran -o ../exe/pme_recip dr_pme_recip.o utility.o erfcfun.o \
../pme_recip/pme_recip.a ../pubfft/libpubfft.a -lm

../exe/vpme_recip:	dr_pme_recip.o utility.o ../vpme_recip/pme_recip.a ../pubfft/libpubfft.a
	gfortran -o ../exe/vpme_recip dr_pme_recip.o utility.o erfcfun.o \
../vpme_recip/pme_recip.a ../pubfft/libpubfft.a -lm

all:	../exe/filter ../exe/exact ../exe/chk_exact ../exe/reg_recip \
		../exe/pme_recip ../exe/vpme_recip

clean:
	/bin/rm -f *.o
