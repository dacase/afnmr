#  Simple top-level Makefile to point users to those hidden below:

install: is_afnmrhome_defined
	./mkdirs
	cd src && make install

clean:
	-(cd src && make clean)

test::
	cd test && make test

is_afnmrhome_defined:
	@(if [ -z "$(AFNMRHOME)" ] ; then \
	    echo "Error: AFNMRHOME is not defined !" ;\
	    exit 2 ;\
	fi ;\
	)

