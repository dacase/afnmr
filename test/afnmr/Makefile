test: 1d3z_test1 1d3z_test2 1d3z_test4 1d3z_test5 m6ab_test1 clean

1d3z_test1::
	@echo "Running test quantum fragmentation on ubiquitin: orca"
	@echo
	../../bin/afnmr -orca -nomin -list 3 1d3z1
	@../dacdif -a 1.e-5 1d3z1003.orcainp.save 1d3z1003.orcainp
	@../dacdif 1d3z1003.pos.save 1d3z1003.pos
	@/bin/rm -f 1d3z1.afnmr.log 1d3z1.pqr 1d31003.pqr

1d3z_test2::
	@echo "Running test quantum fragmentation on ubiquitin: demon"
	@echo
	../../bin/afnmr -demon5 -nomin -list 3 1d3z1
	@../dacdif -a 1.e-5 1d3z1003.inp.save 1d3z1003.inp
	@../dacdif 1d3z1003.pqr.save 1d3z1003.pqr
	@/bin/rm -f 1d3z1.afnmr.log 1d3z1.pqr

1d3z_test3::
	@echo "Running test quantum fragmentation on ubiquitin: terachem"
	@echo
	../../bin/afnmr -tc -qopt -list 2 1d3z1
	@../dacdif -a 1.e-5 1d3z1002.opt.save 1d3z1002.opt
	@../dacdif -a 1.e-5  1d3z1002.xyz.save 1d3z1002.xyz
	@/bin/rm -f 1d3z1.afnmr.log 1d3z1002.pqr 1d3z1.min.pqr 1d3z1.pqr \
                1d3z1002.pos

1d3z_test4::
	@echo "Running test quantum fragmentation on ubiquitin: gaussian"
	@echo
	../../bin/afnmr -gaussian -nomin -list 3 1d3z1
	@../dacdif -a 1.e-5 1d3z1003.com.save 1d3z1003.com
	@/bin/rm -f 1d3z1.afnmr.log 1d3z1003.pqr 1d3z1.min.pqr 1d3z1.pqr

1d3z_test5::
	@echo "Running test quantum fragmentation on ubiquitin: qchem"
	@echo
	../../bin/afnmr -qchem -nomin -list 5 1d3z1
	@../dacdif -a 1.e-5 1d3z1005.in.save 1d3z1005.in
	@/bin/rm -f 1d3z1.afnmr.log 1d3z1005.pqr 1d3z1.min.pqr 1d3z1.pqr

1d3z_test6::
	@echo "Running test quantum fragmentation on ubiquitin with xtb:"
	@echo
	../../bin/afnmr -demon5 -xtb -nomin -list 6 1d3z1
	@../dacdif -r 1.e-6 1d3z1006.inp.save 1d3z1006.inp
	@../dacdif 1d3z1006.pqr.save 1d3z1006.pqr
#	@/bin/rm -f 1d3z1.afnmr.log 1d3z1.pqr 1d3z1006.xtb.log

1d3z_test7::
	@echo "Running test quantum fragmentation on ubiquitin: jaguar"
	@echo
	../../bin/afnmr -jaguar -nomin -list 7 1d3z1
	@../dacdif -a 1.e-5 1d3z1007.in.save 1d3z1007.in
	@/bin/rm -f 1d3z1.afnmr.log 1d3z1007.pqr 1d3z1.min.pqr 1d3z1.pqr

m6ab_test1::
	@echo "Running test quantum fragmentation on RNA with demon"
	@echo
	../../bin/afnmr -demon5 -nobuild -list 6 m6ab.1
	@../dacdif -a 1.e-5 m6ab.1006.inp.save m6ab.1006.inp
	@../dacdif m6ab.1006.pqr.save m6ab.1006.pqr
	@/bin/rm -f m6ab.1.afnmr.log 

clean:
	/bin/rm -f 1d3z1.0.pdb

