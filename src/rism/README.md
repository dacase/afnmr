Contributions
========

Please make sure all contributions adhere to the coding style and
naming conventions below, as well as the following standards:

* Doxygen documentation should be created or updated for all
  substantial additions. Please follow the existing documentation
  style.

Coding Style Conventions
=========

Please adhere to the following conventions when contributing code to
RISM.

If you see places where conventions are not followed, feel free to
correct them.  If the changes only affect whitespace and newlines,
they can be included in any commit. Anything more should be in its own
separate commit.

* Variable names should use camel case, even in languages which are
  not case sensitive (such as Fortran).
* Variable names should be descriptive enough to be understood without
  looking up the documentation or definition.
* Variable names should use plain English and not mathematical
  symbols.

* Binary operators (+, -, =, ==, !=, .not., etc.) should be surrounded
  by space on both sides.
* In Fortran, use end if and end do to terminate blocks, not endif and enddo.
* In C, block brackets should be on the same line as the block:
    if (true)
    { // YES
    }
    
    if (true) { // NO
    }
* 

Naming Conventions
==========

The RISM codebase uses several naming conventions for concision. These
are kept to a minimum to avoid obfuscation.

* Ends of spatial variables end in either R or K depending whether
  they are in the real space or Fourier k-space respectively.
* Due to their high frequency usage, correlation and distribution
  function names are abbreviated as follows:
  * DCF - direct correlation function
  * ICF - indirect correlation function
  * RDF - radial distribution function
  * TCF - total correlation function


