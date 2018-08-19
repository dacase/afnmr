include ../config.h

install:
	cd src/leap && $(MAKE) install
	cp tleap $(BINDIR)/tleap
	chmod +x $(BINDIR)/tleap

uninstall:
	cd src/leap && $(MAKE) uninstallteLeap
	-rm -f $(BINDIR)/tleap
	-cd src/leap && $(MAKE) uninstall

clean:
	cd src/leap && $(MAKE) clean

