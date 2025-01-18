------------------------------------------------------------------------------
   Automatic fragmentation (AFNMR) program for calculating chemical shifts

            Tong Zhu, Xiao He, Jason Swails and David A. Case
------------------------------------------------------------------------------

* Full instructions are in the doc/afnmr.pdf file

* License information is in the LICENSE file

* For the impatient:  'make install && make test'

* Best papers to read/cite: 

  *  J. Swails, T. Zhu, X. He and David A. Case. 
       AFNMR: Automated fragmentation quantum mechanical 
       calculation of NMR chemical shifts for biomolecules. 
       J. Biomol. NMR 63, 125-139 (2015). 

  *  D.A. Case.  Using quantum chemistry to estimate chemical 
       shifts in biomolecules. Biophys. Chem. 267, 106476 (2020).


* Version history:
  * Version 1.0 is the initial split from the "shifts" code base
  * Version 1.1 makes it easier to use snapshots for solvated molecular 
       dynamics simulations as inputs.
  * Version 1.2 is the first to be hosted on github, and changes the 
       default basis to pcSseg-0.
  * Version 1.3 updates the way in which reference shieldings are estimated,
       using DFT calculations on reference compounds. 
  * Version 1.3.1 fixes the ways in which waters and ligands are handled.
  * Version 1.4 adds support for jaguar, plus other small tweaks
  * Version 1.5 provides for fixes for non-sequential residue numbers, adds
       tweaks for quantum geometry optimization
  * Version 1.6 fixes a bug in how the external charges were written for ORCA
  * Version 1.6.1 fixes a mistake in the documentation about the -mol2 flag;
       tweaks for handling extra points
  * Version 1.7 renames basis functions, adds xtb optimization, etc
  * Version 1.8 tweaks how xtb is done, adds -J option for spin-spin coupling
