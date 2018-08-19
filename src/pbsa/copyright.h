!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! PBSA
!
! An analysis program for solvent-mediated electrostatic and non-electrostatic
! interactions in biomolecules.
!
! Please acknowledge your use of PBSA by citing:
!
! Cai, Hsieh, Wang, and Luo, J. Chem. Comp. Theo. 6:203-211, 2009.
! Wang and Luo, J. Comp. Chem. 31:1689-1698, 2010.
! Wang, Cai, Xiang, and Luo, J. Chem. Comp. Theo. 8:2741-2751, 2012.
!
! Major Developers:
!
! Ruxi Qi, CUDA-enabled linear numerical solvers
! Li Xiao, full membrane protein support
! Wes Smith, periodic linear numerical solvers
! Jun Wang, linear Poisson-Boltzmann numerical solvers
! Qin Cai, nonlinear Poisson-Boltzmann numerical solvers
! Meng-Juei Hsieh, program interface and parallel implementations
! Xiang Ye, electrostatic energy and force numerical algorithms
! Chuck Tan, non-electrostatic energy and force numerical algorithms
! Ray Luo, coordinator of overall development
!
! Additional contributing authors are listed in the code documentation.
!
! Departments of Molecular Biology and Biochemistry, Biomedical Engineering,
! and Chemical Engineering and Materials Science, University of California,
! Irvine, California
!
! Copyright (c) 2004-2017. The Regents of the University of California.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Portion of the routine "pb_exmol" is modified from the UHBD program
! (Comp. Phys. Comm. 91:57-95, 1995), copyrighted by University of Houston,
! 1989-2009. See program documentation for more information.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Overview:
!
! PBSA models the electrostatic solvation interaction by the Poisson-Boltzmann
! equation. The implementation uses multiple finite-difference methods to solve
! the partial differential equation. Both linear and full nonlinear numerical
! solvers are implemented. Please refer to the following publications:
!
! Luo, David, and Gilson, J. Comp. Chem. 23:1244-1253, 2002.
! Cai, Wang, Zhao, and Luo, J. Chem. Phys. 130:14501, 2009.
! Wang and Luo, J. Comp. Chem. 31:1689-1698, 2010.
! Smith and Luo, J. Chem. Info. Model. 55:2187-2199, 2015.
!
! for implementation details of the linear solvers.
!
! Porting of linear solvers to the GPU platforms is documented in
!
! R. Qi and R. Luo, Submitted, 2017.
!
! The full nonlinear solvers are documented in:
!
! Cai, Hsieh, Wang, and Luo, J. Chem. Comp. Theo. 6:203-211, 2009.
!
! The electrostatic energy and forces are computed based on the finite-
! difference grid potentials as discussed in the following publications:
!
! Lu and Luo, J. Chem. Phys. 119:11035-11047, 2003.
! Cai, Ye, Wang, and Luo, Chem. Phys. Lett. 514:368-373, 2011.
! Cai, Ye, and Luo, Phys. Chem. Chem. Phys. 14:15917-15925, 2012.
! Xiao, Cai, Ye, and Luo, J. Chem. Phys. 139:094106, 2013.
! Xiao, Wang, Ye, and Luo, J. Phys. Chem. B. 120:8707-8721, 2016.
!
! The dielectric models and molecular surfaces and their discretization in the
! electrostatic solvation model are documented in:
!
! Ye, Wang, and Luo, J. Chem. Theo. Comp. 6:1157-1169, 2010.
! Wang, Cai, Ye, and Luo, J. Chem. Comp. Theo. 8:2741-2751, 2012.
!
! The parameters used in the electrostatic solvation model are optimized and
! validated in the following publications:
!
! Tan, Yang and Luo, J. Phys. Chem. 110:18680-18687, 2006.
! Wang, Tan, Chanco, and Luo, Phys. Chem. Chem. Phys. 12:1194-1202, 2007.
!
! PBSA models the non-electrostatic solvation interaction by two separate
! terms in this release: dispersion (or van der Waals or attractive) and cavity
! (or hydrophobic or repulsive). The dispersion term is computed by a numerical
! integration over the solvent accessible surface area. The cavity term is
! modeled by a term proportional to the molecular surface or volume with a
! single proportional constant. See the following publication for details:
!
! Tan, Tan and Luo, J. Phys. Chem. 111:12263-12274, 2007.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file is part of PBSA.
!
! PBSA is free software; you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free Software
! Foundation; either version 3 of the License, or (at your option) any later
! version.
!
! PBSA is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
! details.
!
! You should have received a copy of the GNU Lesser General Public License along
! with PBSA; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
! Suite 330, Boston, MA 02111-1307, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
