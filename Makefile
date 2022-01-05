#  Simple top-level Makefile to point users to those hidden below:

install: is_afnmrhome_defined
	./mkdirs
	cd src && make install
	@echo "Finished installation of afnmr-1.5"

clean:
	-(cd src && make clean)

uninstall: clean
	/bin/rm -rf bin lib logs include share

test::
	cd test && make test

is_afnmrhome_defined:
	@(if [ -z "$(AFNMRHOME)" ] ; then \
	    echo "Error: AFNMRHOME is not defined !" ;\
	    exit 2 ;\
    elif [ ! "$(AFNMRHOME)" = `pwd` ]; then \
        echo "Error: AFNMRHOME should be `pwd`, " ; \
        echo "       but is currently $(AFNMRHOME); this needs to be fixed!" ; \
	    exit 2 ;\
	fi ;\
	)

