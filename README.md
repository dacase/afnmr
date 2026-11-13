------------------------------------------------------------------------------
   Automatic fragmentation (AFNMR) program for calculating chemical shifts

            Tong Zhu, Xiao He, Jason Swails and David A. Case
------------------------------------------------------------------------------

* Full instructions are in the doc/afnmr.pdf file

* License information is in the LICENSE file

* For the impatient:  have AFNMRHOME point to this folder;
                      then type 'make install && make test'

* Best papers to read/cite: 

  *  J. Swails, T. Zhu, X. He and David A. Case. 
       AFNMR: Automated fragmentation quantum mechanical 
       calculation of NMR chemical shifts for biomolecules. 
       J. Biomol. NMR 63, 125-139 (2015). 

  *  D.A. Case.  Using quantum chemistry to estimate chemical 
       shifts in biomolecules. Biophys. Chem. 267, 106476 (2020).

* Some papers than have used this code:

  * H. Zhang, G. Hou, M. Lu, J. Ahn, I.-J. Byeon, C.J. Langmead, J.R. Perilla, I. Hung, P.L. Gor'kov, Z. Gan, W. Brey, D.A. Case, K. Schulten, A.M. Gronenborn, and T. Polenova. HIV-1 Capsid Function is Regulated by Dynamics: Quantitative Atomic-Resolution Insights by Integrating Magic-Angle-Spinning NMR, QM/MM, and MD. J. Am. Chem. Soc. 138, 14066-14075 (2016).

   * M. Fritz, C.M. Quinn, M. Wang, G. Hou, X. Lu, L.M.I. Koharudin, J. Struppe, D.A. Case, T. Polenova and A.M. Gronenborn. Accurate determination of backbone chemical shift tensors in microcrystalline proteins by integrated MAS NMR and QM/MM. Phys. Chem. Chem. Phys. 20, 9543-9553 (2018).

   * H. Shi, M.C. Clay, A. Rangadurai, B. Sathyamoorthy, D.A. Case, and H.M. Al-Hashimi. Atomic Structures of Excited State A-T Hoogsteen Base Pairs in Duplex DNA by Combining NMR Relaxation Dispersion, Mutagenesis, and Chemical Shift Calculations. J. Biomol. NMR 70, 229-244 (2018).

   * H. Zhou, B. Sathyamoorthy, A. Stelling, Y. Xu, Y. Xue, Y.Z. Pigli, D.A. Case, P.A. Rice, and H.M. Al-Hashimi. Characterizing Watson–Crick versus Hoogsteen Base Pairing in a DNA–Protein Complex Using Nuclear Magnetic Resonance and Site-Specifically 13C-and 15N-Labeled DNA. Biochemistry 58, 1963-1974 (2019).

   * H. Shi, A. Rangadurai, H.A. Assi, R. Roy, D.A. Case, D. Herschlag, J.D. Yesselman, and H.M. Al-Hashimi. Rapid and accurate determination of atomistic RNA dynamic ensemble models using NMR and structure prediction. Nature Commun. 11, 5531 (2020).

   * B. Liu, H. Shi, A. Rangadurai, C.-C. Chu, F. Nussbaumer, K.A. Erharter, D.A. Case, C. Kreutz, and H.M. Al-Hashimi. A quantitative model predicts how m6A reshapes the kinetic landscape of nucleic acid hybridization and conformational transitions. Nature Commun. 12, 5201 (2021).

   * A. Geng, L. Ganser, R. Roy, S. Pratihar, D.A. Case, and H.M. Al-Hashimi. An RNA excited conformational state at atomic resolution. Nature Commun. 14, 8432 (2023).

* Version history: (see tags in git)
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
  * Version 1.7 is a significant update: incorporating pcS-seg2 basis, new
       mixed basis, tweaked fragmentation.
  * Version 1.8 add support for spin-spin couplings, updates to xtb optimization
  * Version 1.9 increases the cutoff for protein-water interactions; allows
       xtb optimization to go to completion or 500 cycles
  * Version 1.10 uses 3-step rattle minimization; avoids need for cpptraj
