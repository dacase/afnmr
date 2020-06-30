#   top-level Makefile for afnmr
include config.h

install: 
	(cd blas && make $(BLAS) )
	(cd lapack && make $(LAPACK) )
	(cd leap && make install )
	(cd etc && make install )
	(cd sff && make install )
	(cd MEAD && make install)
	(cd rism && make $(RISM) )
	(cd afnmr && make install )

clean::
	-(cd blas && make clean )
	-(cd lapack && make clean )
	-(cd leap && make clean )
	-(cd etc && make clean )
	-(cd sff && make clean )
	-(cd MEAD && make clean)
	-(cd rism && make clean)
	-(cd afnmr && make clean)

