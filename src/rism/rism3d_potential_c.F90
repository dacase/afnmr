!<compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2011-2012 by
! Andriy Kovalenko, Tyler Luchko and David A. Case.
!
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU General Public License in the
! ../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
! Users of the 3D-RISM capability found here are requested to acknowledge
! use of the software in reports and publications.  Such acknowledgement
! should include the following citations:
!
! 1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999);
! ibid. 112:10391-10417 (2000).
!
! 2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by
! F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.
!
! 3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
! and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010).

#include "../include/dprec.fh"

!> A wrapper module to determine k-space long-range asymptotics cut
!! offs using Newton-Raphson root finding.
!!
!! DCF and TCF long-range asymptotics both have similar leading
!! factors that depend only on the wave vector squared, k2, and not on
!! the position of the atoms. The contribution of the atom positions
!! is no more that the total charge of the solute. A suitable cutoff
!! for k2 is when this factor, multiplied by the total charge, is very
!! small. We can determine this by solving for k2 that gives a
!! suitably small number.
!!
!! Since the TCF also has a factor of 1/(dielectric constant) and adds
!! the square of inverse Debye length to k2, it always has a longer
!! cut off and we can just use the cut off for the DCF.
!!
!! See Eqs. 9b and 10a of S. Genheden,
!! T. Luchko, S. Gusarov, A. Kovalenko, and U. Ryde, The Journal of
!! Physical Chemistry B 114, 8505 (2010).
!!
!! We can't pass internal functions as arguments until gfortran 4.6,
!! so wrapping this in a module is required.
module asympk_cut
  implicit none
  _REAL_ :: DCFTolerance, DCFCoefficient, SmearSq
  private :: DCFTolerance, DCFCoefficient, SmearSq
  contains

    !> Calculates the minimum cut off distance to achieve a maximum DCF error of tolerance.
    !! @param[in] totalCharge The total charge of the solute
    !! @param[in] solvent_charges Array of charges for each solvent site.
    !! @param[in] smear Gaussian charge smearing
    !! @param[in] boxVolume Volume of supercell
    !! @param[in] tolerance acceptable error due to the cutoff
    function asympck_cut_calc(soluteCharges, solventCharges, smear, boxVolume, tolerance) &
         result (cutoff)
      use rism_util, only: root_newton
      use constants, only : FOURPI
      implicit none
      _REAL_ :: cutoff
      _REAL_, intent(in) :: soluteCharges(:), solventCharges(:), tolerance, smear, boxVolume
      DCFTolerance = tolerance
      DCFCoefficient = maxval(abs(soluteCharges))*FOURPI*sqrt(2d0)/boxVolume*maxval(abs(solventCharges))
      smearSq = smear**2
      ! We need to that we approach the
      ! root from the left as approaching the root from the right can
      ! cause the Newton to overshoot and end-up with NaNs.  In the
      ! case of a NaN cutoff, the effect is to have no cutoff.
      cutoff = root_newton(asympck_cut, asympck_cut_deriv, 1d0, 1d-16)
    end function asympck_cut_calc
    
    !> Expression to find the k-space cut off for the DCF asymptotics.
    !! @param[in] k2 The k-space vector squared. This is a real number.
    !! @return A real number.  The value of this function.
    function asympck_cut(k2)
      implicit none
      _REAL_ :: asympck_cut
      _REAL_, intent(in) :: k2
      asympck_cut = DCFCoefficient/k2*exp(-0.25d0*k2*smearSq) - DCFTolerance
    end function asympck_cut

    !> Derivative of asympck_cut()
    !! @param[in] k2 The k-space vector squared. This is a real number.
    !! @return A real number.  The value of this function.
    function asympck_cut_deriv(k2)
      implicit none
      _REAL_ :: asympck_cut_deriv
      _REAL_, intent(in) :: k2
      asympck_cut_deriv = -DCFCoefficient/k2**2*exp(-0.25d0*k2*smearSq)*(1d0+0.25d0*smearSq*k2)
    end function asympck_cut_deriv

end module asympk_cut
  
!> A wrapper module to determine r-space Lennard-Jones cutoff using
!! Newton-Raphson root finding.
!!
!! We can't pass internal functions as arguments until gfortran 4.6,
!! so wrapping this in a module is required.
module lj_cut
  implicit none
  _REAL_ :: LJTolerance, LJepsilon, LJrmin, factor
  integer :: LJpower
  private :: LJTolerance, LJepsilon, LJrmin, LJpower, factor
  contains

    !> Calculates the minimum cut off distance to achieve a maximum LJ error of
    !> tolerance.
    !! @param[in] epsilon LJ epsilon parameter for solute-solvent pair
    !! @param[in] rmin LJ rmin parameter
    !! @param[in] tolerance acceptable error due to the cutoff
    function lennard_jones_cut_calc(epsilon, rmin, power, tolerance) &
         result (cutoff)
      use rism_util, only: root_newton
      use constants, only : FOURPI
      implicit none
      _REAL_ :: cutoff
      _REAL_, intent(in) :: epsilon, rmin, tolerance
      integer, intent(in) :: power
      LJTolerance = tolerance
      LJepsilon = epsilon
      LJrmin = rmin
      LJpower = power
      factor =1d0
      if (power==6) factor=2d0
      cutoff = root_newton(lennard_jones_cut, lennard_jones_cut_deriv, 1d0, 1d-16)
    end function lennard_jones_cut_calc
    
    !> Expression to find the r-space cut off for the LJ potential energy.
    !! @param[in] r The r-space distance. This is a real number.
    !! @return A real number.  The value of this function.
    function lennard_jones_cut(r)
      implicit none
      _REAL_ :: lennard_jones_cut
      _REAL_, intent(in) :: r
      _REAL_ :: rmin_r_power
      rmin_r_power = (LJrmin/r)**LJpower
      lennard_jones_cut = LJepsilon*factor*rmin_r_power - LJTolerance
    end function lennard_jones_cut

    !> Derivative of lennard_jones_cut()
    !! @param[in] r The r-space distance. This is a real number.
    !! @return A real number.  The value of this function.
    function lennard_jones_cut_deriv(r)
      implicit none
      _REAL_ :: lennard_jones_cut_deriv
      _REAL_, intent(in) :: r
      _REAL_ :: rmin_r_power
      rmin_r_power = (LJrmin/r)**(LJpower)
      lennard_jones_cut_deriv = -factor*(LJpower)*LJepsilon*rmin_r_power/r
    end function lennard_jones_cut_deriv

  end module lj_cut
  
!> Electrostatic potential class for 3D-RISM. Used to calculate/store
!! quantities that are potential dependent and do not change while
!! converging a solution calculation.
!!
!! Pointers to solute and solvent objects are maintained.  So, if
!! values change in these objects, they are automatically used in the
!! potential calculation.
!!
!! This class is generally MPI agnostic.  That is, MPI is only an
!! issue for setting the grid size where both information about the
!! size of the local slab and the total grid must be supplied.  There
!! is no MPI communication within the class.

module rism3d_potential_c
  use safemem
  use rism3d_solute_c
  use rism3d_solvent_c
  use rism3d_grid_c
  use rism_timer_c
  use rism3d_fft_c ! Only used for periodic code.
#include "def_time.h"
  type rism3d_potential
     !> Only called by periodic code.
     type(rism3d_fft), pointer :: fft => NULL()

     !> Cutoff for RISM potential calculations.
     _REAL_ :: cutoff
     !> Cutoff**2 for RISM LJ potential calculations.
     _REAL_ :: cutoff2
     _REAL_, pointer :: ljCutoffs2(:,:) => NULL()
     !> Cutoff for long range asymptotics of the total correlation function.
     _REAL_ :: cut_hlr
     !> Cutoff for long range asymptotics of the direct correlation function.
     _REAL_ :: cut_clr
     !> Cutoff for k-space long range asymptotics of the direct and total correlation functions.
     _REAL_ :: cut2_chlk

     !> Pointer to grid object.
     type(rism3d_grid), pointer :: grid => NULL()

     !> Pointer to solute object.
     type(rism3d_solute), pointer :: solute => NULL()
     !> Pointer to solvent object.
     type(rism3d_solvent), pointer :: solvent => NULL()

     !> Total time for potential and asymptotics.
     type(rism_timer) :: timer
     !> Timer for LJ potential.
     type(rism_timer) :: ljTimer
     !> Timer for Coulomb potential.
     type(rism_timer) :: coulombTimer
     !> Timer for long range portion Coulomb potential.
     type(rism_timer) :: coulombLongRangeTimer
     !> Timer for short range portion of Coulomb potential.
     type(rism_timer) :: coulombShortRangeTimer
     !> Timer for long range asymptotics.
     type(rism_timer) :: asympTimer
     !> Timer for long range asymptotics: TCF r-space.
     type(rism_timer) :: asympTcfRTimer
     !> Timer for long range asymptotics Dcf & TCF k-space.
     type(rism_timer) :: asympDcfTcfKTimer
     !> Timer for long range asymptotics DCF r-space.
     type(rism_timer) :: asympDcfRTimer

     !> Potential energy of the solvent about the solute.  This is
     !! recalculated for each solution but we want to reserve the
     !! memory and may want to change it if the box dimensions
     !! change. The fourth dimension is per solvent atom. [kT]
     _REAL_, pointer :: uuv(:,:,:,:) => NULL()

     !> The long range portion of the Ewald potential in k-space.
     !! Every two items are the real and imaginary components of the
     !! complex exponential. Thus even though only half the box points
     !! have their potential calculated, the array must have
     !! dimensions of the total number of k-space grid points. This
     !! also allows it to store the real space output which is for
     !! every grid point.
     _REAL_, pointer :: uuv1d(:,:) => NULL()

     !> Solute-solvent sigma interaction matrix.  Calculated once and
     !! used in later calls. [A]
     _REAL_, pointer :: ljSigmaUV(:,:) => NULL()
     !> Solute-solvent epsilon interaction matrix.  Calculated once and
     !! used in later calls. [kT]
     _REAL_, pointer :: ljEpsilonUV(:,:) => NULL()

     !> Solute-solvent LJ A coefficient
     _REAL_, pointer :: ljAUV(:,:) => NULL()
     !> Solute-solvent LJ B coefficient
     _REAL_, pointer :: ljBUV(:,:) => NULL()

     !> Supercell real space asymptotics for direct correlation function.
     _REAL_, pointer :: dcfLongRangeAsympR(:) => NULL()
     !> Supercell k-space asymptotics for direct correlation function.
     _REAL_, pointer :: dcfLongRangeAsympK(:) => NULL()
     !> Supercell real space asymptotics for total correlation function.
     _REAL_, pointer :: tcfLongRangeAsympR(:) => NULL()
     !> Supercell k-space asymptotics for total correlation function.
     _REAL_, pointer :: tcfLongRangeAsympK(:) => NULL()

     !> Long-range part of Huv(k) at k = 0 (2, solv%natom).
     _REAL_, pointer :: huvk0(:,:) => NULL()
     !> Long-range part of temperature derivative Huv(k) at k = 0 (2, solv%natom).
     _REAL_, pointer :: huvk0_dT(:,:) => NULL()

     !> If true, a periodic 3D-RISM calculation is performed. This
     !! mostly differs from simple 3D-RISM by using Ewald sum potential
     !! in place of Coulombic potential and invoking the minimum image
     !! convention while calculating potentials.
     logical :: periodic = .false.

     !> Specifies the periodic potential to use.
     !! Current valid values include:
     !!   'ewald' = Ewald sum potential
     !!   'pme'   = Particle Mesh Ewald potential
     character(len=256) :: periodicPotential = ''

     !> User specified uniform potential bias.  Currently only used
     !! for periodic solute. Applied equally to all solvent
     !! sites.
     _REAL_ :: biasPotential
     _REAL_, pointer :: phineut(:) => NULL()
     ! potential energy and long-range asymptotics options
     
     !> do treecode DCF
     logical :: treeDCF
     !> do treecode TCF
     logical :: treeTCF
     !> do treecode sum Coulomb
     logical :: treeCoulomb
     !> treecode multipole acceptance criterion parameter for DCF
     _REAL_ :: treeDCFMAC
     !> treecode multipole acceptance criterion parameter for TCF
     _REAL_ :: treeTCFMAC
     !> treecode multipole acceptance criterion parameter for Coulomb
     _REAL_ :: treeCoulombMAC
     !> treecode order of Taylor approximation for DCF
     integer :: treeDCFOrder
     !> treecode order of Taylor approximation for TCF
     integer :: treeTCFOrder
     !> treecode order of Taylor approximation for Coulomb
     integer :: treeCoulombOrder
     !> treecode maximum leaf size for DCF
     integer :: treeDCFN0
     !> treecode maximum leaf size for TCF
     integer :: treeTCFN0
     !> treecode maximum leaf size for Coulomb
     integer :: treeCoulombN0
     !> Charge smearing parameter for long-range
     !! asymtotics and Ewald, typically eta in the literature
     _REAL_ :: chargeSmear

     !> indicates if it is approriate to apply the LJ truncation
     !! correction for thermodynamic values
     logical :: applyLJCorrection = .false.
  end type rism3d_potential
  
  private mixSoluteSolventLJParameters
  private uvCoulombicPotential, uvEwaldSumPotential
  private uvLennardJonesPotentialWithCutoff, uvLennardJonesPotentialWithMinimumImage
  private getnojellywt, uvLJrEwaldPotentialWithMinimumImage, uvParticleMeshRecipEwaldPotential
contains


  !> Constructor.
  !! @param[in,out] this potential object
  !! @param[in] grid grid object.  A pointer to this will be retained.
  !! @param[in] solv solvent object.  A pointer to this will be retained.
  !! @param[in] solu solute object.  A pointer to this will be retained.
  !! @param[in] cut Cutoff.
  !! @param[in] fft Fast Fourier Transform object.
  !! @param[in] periodic True when calculating potentials for periodic solute.
  !! @param[in] biasPotential Uniform bias to the electrostatic potential (Coulomb or Ewald sum).
  !! @param[in] treeDCF :: perform treecode DCF
  !! @param[in] treeTCF :: perform treecode TCF
  !! @param[in] treeCoulomb :: perform treecode Coulomb potential
  !! @param[in] treeDCFMAC :: treecode multipole acceptance parameter for DCF
  !! @param[in] treeTCFMAC :: treecode multipole acceptance parameter for TCF
  !! @param[in] treeCoulombMAC :: treecode multipole acceptance parameter for Coulomb
  !! @param[in] treeDCFOrder :: treecode order parameter for DCF
  !! @param[in] treeTCFOrder :: treecode order parameter for TCF  
  !! @param[in] treeCoulombOrder :: treecode order parameter for Coulomb
  !! @param[in] treeDCFN0 :: treecode maximum leaf size for DCF
  !! @param[in] treeTCFN0 :: treecode maximum leaf size for TCF  
  !! @param[in] treeCoulombN0 :: treecode maximum leaf size for Coulomb
  !! @param[in] chargeSmear :: Charge smearing parameter for long-range
  !!       asymtotics and Ewald, typically eta in the literature
  subroutine rism3d_potential_new(this, grid, solv, solu, cut, fft, periodicPotential, &
       biasPotential,&
       treeDCF, treeTCF, treeCoulomb, &
       treeDCFMAC, treeTCFMAC, treeCoulombMAC, &
       treeDCFOrder, treeTCFOrder, treeCoulombOrder, &
       treeDCFN0, treeTCFN0, treeCoulombN0, &
       chargeSmear)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    type(rism3d_grid), target, intent(in) :: grid
    type(rism3d_solute), target, intent(in) :: solu
    type(rism3d_solvent), target, intent(in) :: solv
    _REAL_, intent(in):: cut
    type(rism3d_fft), target, intent(in) :: fft
    character(len=*), intent(in) :: periodicPotential
    _REAL_, intent(in) :: biasPotential
    logical, intent(in) :: treeDCF
    logical, intent(in) :: treeTCF
    logical, intent(in) :: treeCoulomb
    _REAL_, intent(in) :: treeDCFMAC
    _REAL_, intent(in) :: treeTCFMAC
    _REAL_, intent(in) :: treeCoulombMAC
    integer, intent(in) :: treeDCFOrder
    integer, intent(in) :: treeTCFOrder
    integer, intent(in) :: treeCoulombOrder
    integer, intent(in) :: treeDCFN0
    integer, intent(in) :: treeTCFN0
    integer, intent(in) :: treeCoulombN0
    _REAL_, intent(in) :: chargeSmear
    this%grid => grid
    this%solvent => solv
    this%solute => solu
    this%fft => fft
    this%periodicPotential = periodicPotential
    if (this%periodicPotential /= '') then
       this%periodic = .true.
    end if
    this%biasPotential = biasPotential
    this%treeDCF = treeDCF
    this%treeTCF = treeTCF
    this%treeCoulomb = treeCoulomb
    this%treeDCFMAC = treeDCFMAC
    this%treeTCFMAC = treeTCFMAC
    this%treeCoulombMAC = treeCoulombMAC
    this%treeDCFOrder = treeDCFOrder
    this%treeTCFOrder = treeTCFOrder
    this%treeCoulombOrder = treeCoulombOrder
    this%treeDCFN0 = treeDCFN0
    this%treeTCFN0 = treeTCFN0
    this%treeCoulombN0 = treeCoulombN0
    this%chargeSmear = chargeSmear
    this%ljCutoffs2 => safemem_realloc(this%ljCutoffs2, this%solute%numAtoms, this%solvent%numAtomTypes, .false.)
    this%ljSigmaUV => safemem_realloc(this%ljSigmaUV, this%solute%numAtoms, this%solvent%numAtomTypes, .false.)
    this%ljEpsilonUV => safemem_realloc(this%ljEpsilonUV, this%solute%numAtoms, this%solvent%numAtomTypes, .false.)
    this%ljAUV => safemem_realloc(this%ljAUV, this%solute%numAtoms, this%solvent%numAtomTypes, .false.)
    this%ljBUV => safemem_realloc(this%ljBUV, this%solute%numAtoms, this%solvent%numAtomTypes, .false.)
    call rism3d_potential_setCut_ljdistance(this, cut)
    call mixSoluteSolventLJParameters(this)
    call rism_timer_new(this%timer, "Potential")
    call rism_timer_new(this%ljTimer, "Lennard-Jones")
    call rism_timer_setParent(this%ljTimer, this%timer)
    call rism_timer_new(this%coulombTimer, "Coulomb")
    call rism_timer_setParent(this%coulombTimer, this%timer)
    if (this%periodic) then
       call rism_timer_new(this%coulombLongRangeTimer, "Long-Range Coulomb")
       call rism_timer_setParent(this%coulombLongRangeTimer, this%coulombTimer)
       call rism_timer_new(this%coulombShortRangeTimer, "Short-Range Coulomb")
       call rism_timer_setParent(this%coulombShortRangeTimer, this%coulombTimer)
    end if
    call rism_timer_new(this%asympTimer, "Asymptotics")
    call rism_timer_setParent(this%asympTimer, this%timer)
    call rism_timer_new(this%asympTcfRTimer, "Asymptotics TCF-R")
    call rism_timer_setParent(this%asympTcfRTimer, this%asympTimer)
    call rism_timer_new(this%asympDcfTcfKTimer, "Asymptotics Dcf/TCF-K")
    call rism_timer_setParent(this%asympDcfTcfKTimer, this%asympTimer)
    call rism_timer_new(this%asympDcfRTimer, "Asymptotics DCF-R")
    call rism_timer_setParent(this%asympDcfRTimer, this%asympTimer)
    
  end subroutine rism3d_potential_new


  !> Set parent for this timer
  !! @param[in,out] this rism3d_potential object.
  !! @param[in,out] parent Parent timer object.
  subroutine rism3d_potential_setTimerParent(this, parent)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    type(rism_timer), intent(inout) :: parent
    call rism_timer_start(this%timer)
    call rism_timer_setParent(this%timer, parent)
    call rism_timer_stop(this%timer)
  end subroutine rism3d_potential_setTimerParent

  !> Sets cut offs for the Lennard-Jones potential calculations.
  !! @param[in,out] this potential object.
  !! @param[in] ljTolerance The value below which the LJ potential can be neglected.
  subroutine rism3d_potential_setcut_ljtolerance(this, ljTolerance)
    use lj_cut
    implicit none
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(in) :: ljTolerance

    integer :: iu, iv

    !find the LJ cutoffs for this atom and each of the solvent sites
    if (ljTolerance == 0d0) then
       this%ljCutoffs2=HUGE(1d0)
    else
       do iv = 1, this%solvent%numAtomTypes
          do iu = 1, this%solute%numAtoms
             this%ljCutoffs2(iu,iv) = lennard_jones_cut_calc(this%ljEpsilonUV(iu, iv),&
                  this%ljSigmaUV(iu, iv),6,&
                  ljTolerance)
             this%ljCutoffs2(iu,iv) = this%ljCutoffs2(iu,iv)**2
          end do
       end do
    end if
    this%applyLJCorrection = .true.
  end subroutine rism3d_potential_setcut_ljtolerance
  
  !> Sets cut offs for the k-space long range asymptotics (LRA)
  !! calculations.
  !! @param[in,out] this potential object.
  !! @param[in] asympKTolerance The value below which the long range
  !!     asymptotics can be neglected.
  !! @param[in] boxVolume volume of the solvent box
  subroutine rism3d_potential_setcut_asympktolerance(this, asympKSpaceTolerance, boxVolume)
    use asympk_cut
#ifdef __PGI
    use ieee_arithmetic, only : ieee_is_nan
#endif
    implicit none
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(in) :: asympKSpaceTolerance
    _REAL_, intent(in) :: boxVolume
    ! asymptotics cutoffs
    if (asympKSpaceTolerance ==0d0) then
       this%cut2_chlk = HUGE(1d0)
    else
       this%cut2_chlk = asympck_cut_calc(this%solute%charge, this%solvent%charge,&
            this%chargeSmear, boxVolume, asympKSpaceTolerance)
      ! isnan is not standard, but GNU doesn't support Fortran 2003's
      ! ieee_arithmetic until version 5. Intel 15+ supports both; use
      ! the PGI define to enable the most support for old compilers.
#ifdef __PGI
      if (ieee_is_nan(this%cut2_chlk)) then
#else
      if (isnan(this%cut2_chlk)) then
#endif
         call rism_report_error('Could not converge k-space asymptotics cutoff.  '&
              //'Try using a smaller error tolerance or no cutoff.')

      end if
    end if
  end subroutine rism3d_potential_setcut_asympktolerance

  !> Directly a distance cut off for potential and force.
  !! @param[in,out] this potential object.
  !! @param[in] cut Distance cutoff for potential and force calculations.
  subroutine rism3d_potential_setcut_ljdistance(this, cut)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(in) :: cut

    ! Assign potential cutoff.
    ! Ensure that the square won't overflow.
    this%cutoff = min(sqrt(huge(1d0)), cut)
    this%cutoff2 = this%cutoff**2

    ! non-periodic code now uses this variable
    this%ljCutoffs2 = this%cutoff**2

  end subroutine rism3d_potential_setcut_ljdistance
  
  !> Calculates the potential on the grid.
  !! @param[in] ljTolerance :: Lennard-Jones potential tolerance.
  !!       Only grid points that have an approximate value greater
  !!       than this will be computed.
  subroutine rism3d_potential_calc(this, ljTolerance)
    use rism_util, only : checksum
    use rism3d_opendx, only : rism3d_opendx_write
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_potential), intent(inout) :: this !< potential object.

    _REAL_, intent(in) :: ljTolerance
    integer :: id

    integer :: ix, iy, iz, ierr
    character(len=30) :: filename
    integer :: iu 
    call rism_timer_start(this%timer)
    ! Ensure a grid size has been set.
    if (.not. associated(this%grid%waveVectors)) then
       call rism_report_error("rism3d_potential_calc: grid size not set")
       stop
    end if
    ! Check if the grid size has changed.
    if (ubound(this%uuv, 1) /= this%grid%localDimsR(1) .or. &
         ubound(this%uuv, 2) /= this%grid%localDimsR(2) .or. &
         ubound(this%uuv, 3) /= this%grid%localDimsR(3) .or. &
         .not. associated(this%uuv)) then
       this%uuv => safemem_realloc(this%uuv, &
            this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3),&
            this%solvent%numAtomTypes, .false.)
       this%uuv1d => safemem_realloc(this%uuv1d, this%grid%totalLocalPointsK, &
            this%solvent%numAtomTypes, o_preserve = .false., o_aligned = .true.)
    end if
     
    this%uuv = 0

    call timer_start(TIME_UCOULU)
    call rism_timer_start(this%coulombTimer)

    if (this%solute%charged) then
       if (this%periodic) then
          if (this%periodicPotential(1:3) == 'pme') then
             !write (6,*) "PME Electrostatics."
             call uvParticleMeshRecipEwaldPotential(this, this%uuv)
          else if (this%periodicPotential(1:5) == 'ewald') then
             !write (6,*) "Ewald Sum Electrostatics."
             call uvEwaldSumPotential(this, this%uuv) 
          end if
       else
          if (this%treeCoulomb) then
             ! write(0,*) 'tree coulomb'
             call uvCoulombicPotential_treecode(this, this%uuv)
          else
             ! write(0,*) 'notree coulomb'
             call uvCoulombicPotential(this, this%uuv)
          end if
       end if
    end if


    call rism_timer_stop(this%coulombTimer)
    call timer_stop(TIME_UCOULU)

    call timer_start(TIME_ULJUV)
    call rism_timer_start(this%ljTimer)


    if (this%periodic) then
         if (this%periodicPotential(1:5) == 'ewald' ) then
             call uvLennardJonesPotentialWithMinimumImage(this, this%uuv)
         else if (this%periodicPotential(1:3) == 'pme') then
             call uvLJrEwaldPotentialWithMinimumImage(this, this%uuv)
         end if
    else
       call uvLennardJonesPotentialWithCutoff(this, this%uuv, ljTolerance)
    end if
    call rism_timer_stop(this%ljTimer)
    call timer_stop(TIME_ULJUV)

    ! TODO: this routine should only need to be called once at the beginning:
    call getnojellywt(this)
    call rism_timer_stop(this%timer)

  end subroutine rism3d_potential_calc


  !> Set r-space supercell asymptotic of the total correlation function.
  !! @param[in,out] this Potential object.
  !! @param[in] tcfLongRangeAsympR Pointer to pre-calculated r-space
  !!                long range asymptotics.  This class will try to
  !!                deallocate this memory when destroyed.
  subroutine  rism3d_potential_setAsympTCF(this, tcfLongRangeAsympR)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    _REAL_, target, intent(in) :: tcfLongRangeAsympR(:)
    this%tcfLongRangeAsympR => tcfLongRangeAsympR
  end subroutine rism3d_potential_setAsympTCF


  !> Set r-space supercell asymptotic of the direct correlation function.
  !! @parm[in,out] this rism3d potential object.
  !! @param[in] dcfLongRangeAsympR Pointer to pre-calculated r-space
  !!                long range asymptotics.  This class will try to
  !!                deallocate this memory when destroyed.
  subroutine  rism3d_potential_setAsympDCF(this, dcfLongRangeAsympR)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    _REAL_, target, intent(in) :: dcfLongRangeAsympR(:)
    this%dcfLongRangeAsympR => dcfLongRangeAsympR
  end subroutine rism3d_potential_setAsympDCF


  !> Synthesizes supercell asymptotic function of C and H.
  !! Asymptotics are part of the supercell formalism and, as this
  !! may not be used, asymptotics grids are only allocated if/when
  !! the asymptotics calcluation is requested. In particular, if
  !! the solvent is not ionic, asymh[rk] is not allocated.  This
  !! should be invisible to the user.
  subroutine  rism3d_potential_dcf_tcf_long_range_asymptotics(this)
    use constants, only : PI, FOURPI
    use rism_util, only : checksum
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_potential), intent(inout) :: this !< potential object.

    _REAL_ :: uc1gc, uc1gh
    integer :: igx, igy, igz, ig, iu, iv, ig0
    _REAL_ :: phase, sumcos, sumsin, rx, ry, rz, delx, dely, delz, delyz2, ra, &
         uc1g, sr, uc, ucs, ucsc
    _REAL_ :: xappa2, smear2_4, asympk_const
    ! z-axis offset (important for spatially distributed MPI).
    _REAL_ :: offset

    _REAL_ :: solutePosition(3), totalSoluteCharge

    if (.not. this%solute%charged .and. .not. this%periodic) return
    
    ! logical :: printedImages = .false.

    call rism_timer_start(this%asympTimer)

    ! Ensure a grid size has been set.
    if (.not. associated(this%grid%waveVectors)) then
       call rism_report_error("rism3d_potential_dcf_tcf_long_range_asymptotics: grid size not set")
    end if

    ! Allocate memory if the grid size has changed.
    if (this%grid%totalLocalPointsK /= ubound(this%dcfLongRangeAsympK, 1) .or. &
         product(this%grid%localDimsR) /= ubound(this%dcfLongRangeAsympR, 1) .or. &
         .not. associated(this%dcfLongRangeAsympR) .or. &
         .not. associated(this%tcfLongRangeAsympR)) then
       this%dcfLongRangeAsympR => safemem_realloc(this%dcfLongRangeAsympR, product(this%grid%localDimsR), .false.)
       this%dcfLongRangeAsympK => safemem_realloc(this%dcfLongRangeAsympK, this%grid%totalLocalPointsK, .false.)
       if (this%solvent%ionic) then
          this%tcfLongRangeAsympR => safemem_realloc(this%tcfLongRangeAsympR, product(this%grid%localDimsR), .false.)
          this%tcfLongRangeAsympK => safemem_realloc(this%tcfLongRangeAsympK, this%grid%totalLocalPointsK, .false.)
       end if
    end if
    ! Allocate long range part of Huv(k) if not done already.
    if (.not. associated(this%huvk0)) then
       this%huvk0 => safemem_realloc(this%huvk0, 2, this%solvent%numAtomTypes)
       this%huvk0_dT => safemem_realloc(this%huvk0_dT, 2, this%solvent%numAtomTypes)
    end if

    ! Set some values for the long range asymptotics
    smear2_4 = this%chargeSmear**2 / 4d0
    xappa2 = this%solvent%xappa**2
    asympk_const = FOURPI/this%grid%boxVolume
    
    ! Initialize correlation function arrays.
    this%dcfLongRangeAsympK = 0d0
    if (this%solvent%ionic) then
       this%tcfLongRangeAsympK = 0d0
    end if

    if (.not. this%periodic) then
       !JAJ This is the real-space long range portion of the
       ! solute-solvent Ewald potential!

       ! Asymptotic values for the DCF in R-space.
       ! This is eq. 29 of Kovalenko/Hirata 2000.
       if (.not.this%treeDCF) then
          ! write(0,*) 'DCF R direct'

          call rism_timer_start(this%asympDcfRTimer)
          offset = this%grid%spacing(3) * (this%grid%offsetR(3))
          do igz = 0, this%grid%localDimsR(3) - 1
             rz = igz * this%grid%spacing(3) + offset
             do igy = 0, this%grid%localDimsR(2) - 1
                ry = igy * this%grid%spacing(2)
                do igx = 0, this%grid%localDimsR(1) - 1
                   ig = 1 + igx + igy * this%grid%localDimsR(1) + &
                        igz * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                   rx = igx * this%grid%spacing(1)
                   this%dcfLongRangeAsympR(ig) = &
                        dcf_long_range_asymptotics_R(this, (/rx, ry, rz/))
                   !if (this%solvent%ionic) then
                   !   this%tcfLongRangeAsympR(ig) = tcf_long_range_asymptotics_R(this, (/rx, ry, rz/))
                   !end if
                end do
             end do
          end do
          call rism_timer_stop(this%asympDcfRTimer)
       end if
       if (.not.this%treeTCF) then
          ! write(0,*) 'TCF R direct'
          ! Asymptotic values for  TCF in R-space.
          ! This eq. 32 of Kovalenko/Hirata 2000.
          if (this%solvent%ionic) then
             call rism_timer_start(this%asympTcfRTimer)
             offset = this%grid%spacing(3) * (this%grid%offsetR(3))
             do igz = 0, this%grid%localDimsR(3) - 1
                rz = igz * this%grid%spacing(3) + offset
                do igy = 0, this%grid%localDimsR(2) - 1
                   ry = igy * this%grid%spacing(2)
                   do igx = 0, this%grid%localDimsR(1) - 1
                      ig = 1 + igx + igy * this%grid%localDimsR(1) + &
                           igz * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                      rx = igx * this%grid%spacing(1)
                      ! this%dcfLongRangeAsympR(ig) = dcf_long_range_asymptotics_R(this, (/rx, ry, rz/))

                      this%tcfLongRangeAsympR(ig) = tcf_long_range_asymptotics_R(this, (/rx, ry, rz/))
                   end do
                end do
             end do
             call rism_timer_stop(this%asympTcfRTimer)
          end if
       end if

       if ((this%treeDCF) .or. (this%treeTCF)) then
          ! write(0,*) 'DCF/TCK R tree'
          call dcf_tcf_asymptotics_treecode(this)
       end if

    end if

    if ((.not. this%periodic) .or. (this%periodicPotential(1:5) == 'ewald')) then
       ! write(0,*) 'DCF/TCF k direct'

       ! Getting long-range part of the DCF and TCF in k-space.
       ! k = 0 is on master.  Set it to zero as it is handled separately.

       !JAJ This is the long range portion of the solute electrostatic 
       ! periodic field
       ! See Understanding Molecular Simulation 2E, eq. 12.1.15 on p. 297.
       ! Note: the result is divided by volume since FFTW does not do this
       ! during transforms, though technically this is not part of the 
       ! frequency space potential.
       call rism_timer_start(this%asympDcfTcfKTimer)
       if (this%grid%offsetK(3) == 0) then
          this%dcfLongRangeAsympK(1) = 0d0
          this%dcfLongRangeAsympK(2) = 0d0
          if (this%solvent%ionic) then
             this%tcfLongRangeAsympK(1) = 0d0
             this%tcfLongRangeAsympK(2) = 0d0
          end if
          ig0 = 2
       else
          ig0 = 1
       end if

       do ig = ig0, this%grid%totalLocalPointsK / 2
          ! If the wave vector squared is longer than the cut off,
          ! skip to the next value. The values are not ordered by
          ! size.
          if (this%grid%waveVectors2(ig) > this%cut2_chlk) then
             ! write(0,*) this%cut2_chlk, this%grid%waveVectors2(ig)
             this%dcfLongRangeAsympK(2 * ig - 1) = 0d0
             this%dcfLongRangeAsympK(2 * ig) = 0d0
             if (this%solvent%ionic) then
                this%tcfLongRangeAsympK(2 * ig - 1) = 0d0
                this%tcfLongRangeAsympK(2 * ig) = 0d0
             end if
             cycle
          end if
          sumcos = 0d0
          sumsin = 0d0
          if (this%periodic) then
             do iu = 1, this%solute%numAtoms
                solutePosition = this%solute%position(:, iu)
                solutePosition = minimumImage(this, solutePosition)
                phase = dot_product(this%grid%waveVectors(:, ig), solutePosition)
                sumcos = sumcos + this%solute%charge(iu) * cos(phase)
                sumsin = sumsin + this%solute%charge(iu) * sin(phase)
             end do
          else
             do iu = 1, this%solute%numAtoms
                solutePosition = this%solute%position(:, iu)
                phase = dot_product(this%grid%waveVectors(:, ig), solutePosition)
                sumcos = sumcos + this%solute%charge(iu) * cos(phase)
                sumsin = sumsin + this%solute%charge(iu) * sin(phase)
             end do
          end if

          ! DCF in k-space.
          uc1g = asympk_const * exp(-smear2_4 * this%grid%waveVectors2(ig))
          uc1gc = uc1g / this%grid%waveVectors2(ig)
          this%dcfLongRangeAsympK(2 * ig - 1) = uc1gc * sumcos
          this%dcfLongRangeAsympK(2 * ig) = uc1gc * sumsin

          ! TCF in k-space.
          if (this%solvent%ionic) then
             uc1gh = uc1g / ((this%grid%waveVectors2(ig) + xappa2) * this%solvent%dielconst)
             this%tcfLongRangeAsympK(2 * ig - 1) = uc1gh * sumcos
             this%tcfLongRangeAsympK(2 * ig) = uc1gh * sumsin
          end if
       end do
       call rism_timer_stop(this%asympDcfTcfKTimer)
    end if

    if (.not. this%periodic) then
       ! write(0,*) 'DCF/TCF k direct k=0'
       ! Getting the difference between long-range asymptotic functions
       ! of the TCF at k = 0.

       call rism_timer_start(this%asympDcfTcfKTimer)
       if (this%grid%offsetK(3) == 0) then
          sumcos = 0d0
          sumsin = 0d0
          do iu = 1, this%solute%numAtoms
             solutePosition = this%solute%position(:, iu)
             ! if (this%periodic) then
             !    ! Minimum image convention.
             !    !FIXME: This should be relative to the current grid point position!
             !    !FIXME: This should never be reached anyways, so it
             !    ! probably should be removed.
             !    solutePosition = minimumImage(this, solutePosition)
             ! end if
             phase = dot_product(this%grid%waveVectors(:, 1), solutePosition)
             sumcos = sumcos + this%solute%charge(iu) * cos(phase)
             sumsin = sumsin + this%solute%charge(iu) * sin(phase)
          end do
          do iv = 1, this%solvent%numAtomTypes
             this%huvk0(1, iv) = this%solvent%delhv0(iv) * sumcos / this%grid%boxVolume
             this%huvk0(2, iv) = this%solvent%delhv0(iv) * sumsin / this%grid%boxVolume
             ! Make sure we have read the information to do this.
             if (this%solvent%delhv0_dT(iv) /= huge(1d0)) then
                this%huvk0_dT(1, iv) = this%solvent%delhv0_dT(iv) * sumcos / this%grid%boxVolume
                this%huvk0_dT(2, iv) = this%solvent%delhv0_dT(iv) * sumsin / this%grid%boxVolume
             end if
          end do
       end if
       call rism_timer_stop(this%asympDcfTcfKTimer)
    end if
    call rism_timer_stop(this%asympTimer)
  end subroutine rism3d_potential_dcf_tcf_long_range_asymptotics



  !> Calculates direct and total correlation functions on grid
  !! using a treecode expansion
  !! @param[in,out] this rism3d potential object.
  subroutine dcf_tcf_asymptotics_treecode(this)
    use treecode_procedures
    implicit none
    type(rism3d_potential), intent(inout) :: this
    
    !Parameters found in/ calculated from rism3d_potential struct
    integer :: numSources, numTargets
    _REAL_ :: paramKappa, paramEta, paramDiel

    !Local variables
    integer :: iu
    _REAL_ :: offset
    _REAL_ :: targetGridLimits(6)


    !Setting values for treecode from members of rism3d_potential struct
    numSources = this%solute%numAtoms
    numTargets = product(this%grid%localDimsR)

    paramKappa = this%solvent%xappa
    paramEta = this%chargeSmear
    paramDiel = this%solvent%dielconst

    !Constructing grid arrays
    offset = this%grid%spacing(3) * this%grid%offsetR(3)

    targetGridLimits(1) = 0
    targetGridLimits(2) = (this%grid%localDimsR(1) - 1) * this%grid%spacing(1)
    targetGridLimits(3) = 0
    targetGridLimits(4) = (this%grid%localDimsR(2) - 1) * this%grid%spacing(2)
    targetGridLimits(5) = offset
    targetGridLimits(6) = (this%grid%localDimsR(3) - 1) * this%grid%spacing(3) + offset
        
    !Calling the driver routine for the treecode.
    if (this%solvent%ionic .and. this%treeTCF) then
        call rism_timer_start(this%asympTcfRTimer)
        call TREECODE(this%solute%position, &           !Array of source positions
             this%solute%charge, &                      !Array of source charges
             targetGridLimits, this%grid%localDimsR, &  !Target grid information
             numSources, numTargets, &                  !# of sources and targets
             this%tcfLongRangeAsympR, &                 !Asymptotic TCF grid
             this%treeTCFOrder, this%treeTCFMAC, & !Treecode tuning parameters
             this%treeTCFN0, 0, &                       !0 is TCF
             paramKappa, paramEta, paramDiel)           !Physical parameters
         call rism_timer_stop(this%asympTcfRTimer)
    end if

    if (this%treeDCF) then
       call rism_timer_start(this%asympDcfRTimer)
       call TREECODE(this%solute%position, &            !Array of source positions
            this%solute%charge, &                       !Array of source charges
            targetGridLimits, this%grid%localDimsR, &   !Target grid information
            numSources, numTargets, &                   !# of sources and targets
            this%dcfLongRangeAsympR, &                  !Asymptotic DCF grid
            this%treeDCFOrder, this%treeDCFMAC, &  !Treecode tuning parameters
            this%treeDCFN0, 1, &                        !1 is DCF
            paramKappa, paramEta, paramDiel)            !Pysical parameters
       call rism_timer_stop(this%asympDcfRTimer)
    end if

  end subroutine dcf_tcf_asymptotics_treecode

  !> Calculates Coulombic potential on grid
  !! using a treecode expansion
  !! @param[in] this rism3d potential object.
  !! @param[in,out] ucu grid on which potential is added
  subroutine uvCoulombicPotential_treecode(this, ucu)
    use treecode_procedures
    implicit none
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(inout) :: ucu(:,:,:,:) 
    
    !Values found in/ calculated from rism3d_potential struct
    integer :: numSources, numTargets

    !Local variables
    integer :: iv, igy, igz, igx, nn
    _REAL_ :: offset

    !Local arrays for grid locations and source locations
    _REAL_, pointer :: ucuFlatGrid(:) => NULL()
    _REAL_ :: targetGridLimits(6)


    !Setting parameters for treecode from members of rism3d_potential struct
    numSources = this%solute%numAtoms
    numTargets = product(this%grid%localDimsR)

    !Allocating space for source location arrays
    ucuFlatGrid => safemem_realloc(ucuFlatGrid, numTargets, .false.)

    !Constructing grid arrays
    offset = this%grid%spacing(3) * this%grid%offsetR(3)

    targetGridLimits(1) = 0
    targetGridLimits(2) = (this%grid%localDimsR(1) - 1) * this%grid%spacing(1)
    targetGridLimits(3) = 0
    targetGridLimits(4) = (this%grid%localDimsR(2) - 1) * this%grid%spacing(2)
    targetGridLimits(5) = offset
    targetGridLimits(6) = (this%grid%localDimsR(3) - 1) * this%grid%spacing(3) + offset

    !Calling the driver routine for the treecode.
    call TREECODE(this%solute%position, &           !Array of source charges
            this%solute%charge, &                   !Array of source charges
            targetGridLimits, &                     !Target grid information
            this%grid%localDimsR, &
            numSources,numTargets, &                !# of sources and targets
            ucuFlatGrid, &                          !Coulombic potential grid
            this%treeCoulombOrder, &                !Treecode tuning parameters
            this%treeCoulombMAC, &
            this%treeCoulombN0, 2, 0, 0, 0)         !2 is Coulombic potential

    nn = 0;
    do igz = 1, this%grid%localDimsR(3)
       do igy = 1, this%grid%localDimsR(2)
          do igx = 1, this%grid%localDimsR(1)
             nn = nn + 1;
             ucu(igx, igy, igz, 1) = -ucuFlatGrid(nn);
          end do
       end do
    end do

    do iv = this%solvent%numAtomTypes, 1, -1
       ucu(:,:,:,iv) =  ucu(:,:,:,1) * this%solvent%charge(iv)
    end do

    if (safemem_dealloc(ucuFlatGrid) /= 0) &
          call rism_report_error("ucuFlatGrid deallocation failed")

  end subroutine uvCoulombicPotential_treecode

  !> Calculates integrals of the long range asymptotics h**2 and h * c
  !! for all space (i.e. not just the grid).  These terms are used by
  !! the HNC term in the chemical potential so we calculate them
  !! together.
  !! @param[in] this the closure object
  !! @param[out] h2 integral of tcfLongRangeAsympR**2
  !! @param[out] hc integral of tcfLongRangeAsympR * dcfLongRangeAsympR
  subroutine rism3d_potential_int_h2_hc(this, h2, hc)
    use constants, only : PI
    use rism_util, only : gaussquad_legendre, checksum
    implicit none
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(out) :: h2(this%solvent%numAtomTypes), hc(this%solvent%numAtomTypes)

    _REAL_ :: sumhc, sumh2, sumb, xarg, x2arg, denom
    _REAL_ :: dx, dy , dz, r2, k, k2
    integer :: i, j, ik, iv , n0, nf
    integer, parameter :: Nintmx = 200
    _REAL_ :: argum(Nintmx), weight(Nintmx)

    if (.not. this%solute%charged .or. .not. this%solvent%ionic) then
       h2 = 0d0
       hc = 0d0
       return
    end if

    ! Get the numerical long-range contribution.  This uses
    ! Gauss-Legendre quadrature where the integration bounds have been
    ! converted from 0 to infinity to 0 to 1.
    ! See Genheden et al. J. Phys. Chem. B 114, 8505 (2010) Eq. 14

    ! Initialize Gauss-Legendre integration.
    call gaussquad_legendre (0d0, 1d0, argum, weight, Nintmx)
    sumhc = 0.d0
    sumh2 = 0.d0

    ! Parallelize over k-loop.  This is good for a maximum of Nintmx processes.
    ! Select range of sum for this process.
    ! Get the start point for this process.
    n0 = Nintmx / this%grid%mpisize * this%grid%mpirank +1

    ! The final point is the start point for one process
    ! higher.  Ensure we do not exceed the desired range of the
    ! sum
    nf = min(Nintmx / this%grid%mpisize * (this%grid%mpirank + 1), Nintmx)

    do ik = n0, nf

       k = argum(ik) / (1.d0 - argum(ik))

       ! Bessel part.
       sumb = 0.d0

       do i = 2, this%solute%numAtoms
          do j = 1, i - 1

             ! Site separation.
             dx = this%solute%position(1, i) - this%solute%position(1, j)
             dy = this%solute%position(2, i) - this%solute%position(2, j)
             dz = this%solute%position(3, i) - this%solute%position(3, j)

             r2 = dx * dx + dy * dy + dz * dz

             xarg = this%solvent%xappa * k * sqrt(r2)
             if (xarg == 0.d0) &
                  sumb = sumb + 1.d0 * this%solute%charge(i) * this%solute%charge(j)
             if (xarg /= 0.d0) &
                  sumb = sumb + sin(xarg) / xarg *this%solute%charge(i) * this%solute%charge(j)

          end do
       end do

       sumb = sumb * 2.d0

       do i = 1, this%solute%numAtoms
          sumb = sumb + 1.d0 * this%solute%charge(i) * this%solute%charge(i)
       end do

       ! End of Bessel part.

       x2arg = (k * this%solvent%xappa * this%chargeSmear)**2
       k2 = k**2

       denom = argum(ik)**2 + (1.0d0 - argum(ik))**2

       sumhc = sumhc - &
            exp(-x2arg / 2.d0) / denom * sumb * weight(ik)
       sumh2 = sumh2 + &
            exp(-x2arg / 2.d0) * argum(ik)**2 / (denom**2) * sumb * weight(ik)

    end do
    ! End of k-loop.

    ! Site xappa.
    do iv = 1, this%solvent%numAtomTypes
       h2(iv) = 8.d0 * PI / (this%solvent%dielconst) &
            *this%solvent%charge_sp(iv)**2
       hc(iv) = 8.d0 * PI / (this%solvent%dielconst) &
            *this%solvent%charge(iv) * this%solvent%charge_sp(iv)
    end do

    ! Unit for xappa.
    do iv = 1, this%solvent%numAtomTypes
       h2(iv) = h2(iv)
       hc(iv) = hc(iv)
    end do

    ! Divided by total xappa.
    if (this%solvent%xappa /= 0.d0) then
       do iv = 1, this%solvent%numAtomTypes
          h2(iv) = h2(iv) / this%solvent%xappa
          hc(iv) = hc(iv) / this%solvent%xappa
       end do
    else
       do iv = 1, this%solvent%numAtomTypes
          h2(iv) = 0.d0
          hc(iv) = 0.d0
       end do
    end if

    ! Multiplying xmulr by integral.
    do iv = 1, this%solvent%numAtomTypes
       h2(iv) =   h2(iv) * sumh2
       hc(iv) = - hc(iv) * sumhc
    end do

    ! J devided by KT.
    do iv = 1, this%solvent%numAtomTypes
       h2(iv) = h2(iv) / (PI * this%solvent%dielconst)
       hc(iv) = hc(iv) / PI
    end do
  end subroutine rism3d_potential_int_h2_hc


  !> Calculates integral of the long range asymptotics of h**n for all
  !! space (i.e. not just the grid).
  !! @param[in] this the closure object
  !! @param[in] power power to multiply t* to
  !! @return integral of tcfLongRangeAsympR**n
  function rism3d_potential_int_hn(this, power) result(integral)
    use constants, only : PI
    use rism_util, only : gaussquad_laguerre
    implicit none
    type(rism3d_potential), intent(in) :: this
    integer, intent(in) :: power
    _REAL_ :: integral(this%solvent%numAtomTypes)
    _REAL_ :: temp(this%solvent%numAtomTypes)
    _REAL_ :: total, shell, r(3)
    !      _REAL_, parameter :: err = 1d-8
    _REAL_, parameter :: err = 1d-5
    _REAL_ :: sm, tstar
    integer :: iu, iv, istep, ishell, ix, iy, iz, ig, n0, nf, npt
    ! nstep :: steps for Gauss-Laguerre quadrature
    integer, parameter :: nstep = 100, maxshell = 1000
    _REAL_ :: argum(nstep), weight(nstep)

    integral = 0d0
    if (.not. this%solvent%ionic .or. .not. this%solute%charged) return

    if (power == 2) then
       call rism3d_potential_int_h2_hc(this, integral, temp)
    else if (power == 1) then
       integral = rism3d_potential_int_h(this)
    else
       ! Since this integral converges faster than exp(r), we can use Gauss-Laguerre
       ! quadrature.
       !! $      call gaussquad_laguerre(argum, weight, nstep)

       ! Instead use a brute force grid based method for now.
       ! First sum the precalculated values.
       total = sum(this%tcfLongRangeAsympR**power)

       ! Add sum over shells around the original grid until the
       ! relative error falls below our selected threshold The shells
       ! are constructed as faces of the grid expanded by one point
       ! at each boundary. The faces are two points larger along one
       ! axis, such that the faces 'interlock' to cover the edges.
       ! Corners still need to be handled separately.

       ! NOTE: parallel version implemented here will not give
       ! identical results to the serial version since different
       ! processes may converge in different number of iterations.
       ! This is not a practical issue for the typical precision
       ! needed for 3D-RISM.
       do ishell = 1, maxshell

          shell = 0
          ig = 0
          ! xy plane.
          ! Select range of sum for this process.
          ! Total number of points in the outside loop to evaluate.
          npt = (this%grid%globalDimsR(2) + 2 * ishell)
          ! Get the start point for this process.
          n0 = npt / this%grid%mpisize * this%grid%mpirank - ishell
          ! The final point is the start point for one process
          ! higher.  Ensure we do not exceed the desired range of the
          ! sum.
          nf = min(npt / this%grid%mpisize * (this%grid%mpirank + 1) - ishell - 1&
               , this%grid%globalDimsR(2) - 1 + ishell)
          do iy = n0, nf
             r(2) = iy * this%grid%spacing(2)
             do ix = 1 - ishell, this%grid%globalDimsR(1)-2 + ishell
                r(1) = ix * this%grid%spacing(1)
                ! +z
                r(3) = (this%grid%globalDimsR(3) - 1 + ishell) * this%grid%spacing(3)
                shell = shell + tcf_long_range_asymptotics_R(this, r)**power
                ! -z
                r(3) = -ishell * this%grid%spacing(3)
                shell = shell + tcf_long_range_asymptotics_R(this, r)**power
                ig = ig + 2
             end do
          end do
          ! yz plane.
          npt = (this%grid%globalDimsR(3) + 2 * ishell)
          n0 = npt / this%grid%mpisize * this%grid%mpirank - ishell
          nf = min(npt / this%grid%mpisize * (this%grid%mpirank + 1) - ishell - 1&
               , this%grid%globalDimsR(3) - 1 + ishell)
          do iz = n0, nf
             r(3) = iz * this%grid%spacing(3)
             do iy = 1 - ishell, this%grid%globalDimsR(2) - 2 + ishell
                r(2) = iy * this%grid%spacing(2)
                ! +x
                r(1) = (this%grid%globalDimsR(1) - 1 + ishell) * this%grid%spacing(1)
                shell = shell + tcf_long_range_asymptotics_R(this, r)**power
                ! -x
                r(1) = -ishell * this%grid%spacing(1)
                shell = shell + tcf_long_range_asymptotics_R(this, r)**power
                ig = ig + 2
             end do
          end do
          ! zx plane.
          npt = (this%grid%globalDimsR(3) + 2 * ishell - 2)
          n0 = npt / this%grid%mpisize * this%grid%mpirank +1- ishell
          nf = min(npt / this%grid%mpisize * (this%grid%mpirank + 1) - ishell, &
               this%grid%globalDimsR(3) - 2 + ishell)
          do iz = n0, nf
             r(3) = iz * this%grid%spacing(3)
             do ix = -ishell, this%grid%globalDimsR(1) - 1 + ishell
                r(1) = ix * this%grid%spacing(1)
                !+y
                r(2) = (this%grid%globalDimsR(2) - 1 + ishell) * this%grid%spacing(2)
                shell = shell + tcf_long_range_asymptotics_R(this, r)**power
                ! - y
                r(2) = -ishell * this%grid%spacing(2)
                shell = shell + tcf_long_range_asymptotics_R(this, r)**power
                ig = ig + 2
             end do
          end do
          ! Corners.
          if (this%grid%mpirank == 0) then
             do ix = -ishell, this%grid%globalDimsR(1) - 1 + ishell, &
                  this%grid%globalDimsR(1) - 1 + 2 * ishell
                do iy = -ishell, this%grid%globalDimsR(2) - 1 + ishell, &
                     this%grid%globalDimsR(2) - 1 + 2 * ishell
                   do iz = -ishell, this%grid%globalDimsR(3) - 1 + ishell, &
                        this%grid%globalDimsR(3) - 1 + 2 * ishell
                      r(1) = ix * this%grid%spacing(1)
                      r(2) = iy * this%grid%spacing(2)
                      r(3) = iz * this%grid%spacing(3)
                      shell = shell + tcf_long_range_asymptotics_R(this, r)**power
                      ig = ig + 1
                   end do
                end do
             end do
          end if
          total = total + shell
          if (abs(shell / total) < err) then
             exit
          end if
       end do
       if (abs(shell / total) > err) &
            call rism_report_warn("Long range asymptotics integration failed to converge")
       integral = total * (this%solvent%charge_sp)**power * this%grid%voxelVolume
    end if
  end function rism3d_potential_int_hn


  !> Calculates integral of the long range asymptotics of h for all
  !! space (i.e. not just the grid).
  !! @param[in] this The potential object.
  !! @return Integral of tcfLongRangeAsympR for each solvent site.
  function rism3d_potential_int_h(this) result(h)
    use constants, only : PI
    implicit none
    type(rism3d_potential), intent(in) :: this
    _REAL_ :: h(this%solvent%numAtomTypes)
    _REAL_ :: qut

    if (this%grid%mpirank /= 0) then
       h = 0
       return
    end if

    ! We have an analytic expression for this integral so we only need
    ! to sum this over the solute and solvent sites.
    qut = sum(this%solute%charge)
    h = 0d0
    where (this%solvent%charge_sp /= 0d0)
       h = -4d0 * pi / (this%solvent%dielconst * this%solvent%xappa**2) &
            * this%solvent%charge_sp * qut
    end where
  end function rism3d_potential_int_h


  !> Frees all memory and resets values.
  !! @param[in,out] this potential object.
  subroutine rism3d_potential_destroy(this)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    call rism_timer_destroy(this%timer)
    call rism_timer_destroy(this%ljTimer)
    call rism_timer_destroy(this%coulombTimer)
    if (this%periodic) then
       call rism_timer_destroy(this%coulombLongRangeTimer)
       call rism_timer_destroy(this%coulombShortRangeTimer)
    end if
    call rism_timer_destroy(this%asympTimer)
    call rism_timer_destroy(this%timer)
    nullify(this%grid)
    nullify(this%solvent)
    nullify(this%solute)
    this%cutoff2 = 0
    if (safemem_dealloc(this%uuv) /= 0) &
         call rism_report_error("Uuv deallocation failed")
    if (safemem_dealloc(this%uuv1d,o_aligned=.true.) /= 0) &
         call rism_report_error("Uuv1d deallocation failed")
    if (safemem_dealloc(this%ljSigmaUV) /= 0) &
         call rism_report_error("LjSigmaUV deallocation failed")
    if (safemem_dealloc(this%ljEpsilonUV) /= 0) &
         call rism_report_error("EPSuv deallocation failed")
    if (safemem_dealloc(this%ljAUV) /= 0) &
         call rism_report_error("ljAUV deallocation failed")
    if (safemem_dealloc(this%ljBUV) /= 0) &
         call rism_report_error("ljBUV deallocation failed")
    if (safemem_dealloc(this%ljCutoffs2) /= 0) &
         call rism_report_error("ljCutoffs2 deallocation failed")
    if (safemem_dealloc(this%dcfLongRangeAsympR) /= 0) &
         call rism_report_error("asympcr deallocation failed")
    if (safemem_dealloc(this%dcfLongRangeAsympK) /= 0) &
         call rism_report_error("dcfLongRangeAsympK deallocation failed")
    if (safemem_dealloc(this%tcfLongRangeAsympR) /= 0) &
         call rism_report_error("tcfLongRangeAsympR deallocation failed")
    if (safemem_dealloc(this%tcfLongRangeAsympK) /= 0) &
         call rism_report_error("tcfLongRangeAsympK deallocation failed")
    if (safemem_dealloc(this%huvk0) /= 0) &
         call rism_report_error("huvk0 deallocation failed")
    if (safemem_dealloc(this%huvk0_dT) /= 0) &
         call rism_report_error("huvk0_dT deallocation failed")
  end subroutine rism3d_potential_destroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !> Mix per-site Lennard-Jones solute and solvent parameters to obtain LJ
  !! solute-solvent interaction parameters.
  !! @param[in,out] this potential object.
  subroutine mixSoluteSolventLJParameters(this)
    use safemem
    implicit none
    type(rism3d_potential), intent(inout) :: this
    integer :: iv, iu
    !GMG> where to write, read lj matrix
    integer, parameter :: out_unit=20
    integer, parameter :: sout_unit=22
    integer, parameter :: in_unit=21
    logical :: exist
    integer :: status
    _REAL_ :: sm,em
    !-------------------------------

    ! Read LJ-solute (uu) corrections in case ljsolute.mods.txt file is
    ! present.
    !--------------------------------------------------------------------
    inquire(file="ljsolute.mods.txt", exist=exist)
    if (exist) then
        open(in_unit, file="ljsolute.mods.txt", status="old", action="read")
        do
          read (in_unit,'(i6,2e16.8)',IOSTAT=status) iu,sm,em
          if (status /= 0) exit
          ! write (6,'(a5,i6,2e16.8)') "ljsolute-mods>",iu,sm,em
          ! call
          ! rism_report_message("(a20,i6,2e16.8)","ljsolute-mods>",iu,sm,em)
          this%solute%ljSigma(iu)   = sm
          this%solute%ljEpsilon(iu) = em
        enddo
        call rism_report_message("|-> Read LJ-solute modifications from ljsolute.mods.txt.")

       ! Write the LJ-solute (uu) parameters.
       !--------------------------------------------------------------------
       open (unit=sout_unit,file="ljsolute.orig.txt",action="write",&
             status="replace")
       do iu = 1, this%solute%numAtoms
          write (sout_unit,'(i6,2e16.8)') iu,this%solute%ljSigma(iu),&
                                             this%solute%ljEpsilon(iu)
       end do
       close(sout_unit)
    endif

    ! Compute the LJ-matrix (uv).
    ! ----------------------------------------------------------
    do iv = 1, this%solvent%numAtomTypes
       do iu = 1, this%solute%numAtoms
          this%ljSigmaUV(iu, iv) = this%solute%ljSigma(iu) + &
                                   this%solvent%ljSigma(iv)
          this%ljEpsilonUV(iu, iv) = sqrt(this%solute%ljEpsilon(iu) * &
                                          this%solvent%ljEpsilon(iv))
          this%ljAUV(iu, iv) = this%ljEpsilonUV(iu, iv)*this%ljSigmaUV(iu,iv)**12
          this%ljBUV(iu, iv) = 2d0*this%ljEpsilonUV(iu, iv)*this%ljSigmaUV(iu,iv)**6
       end do
    end do

    ! Read LJ-matrix (uv) corrections in case ljmatrix.mods.txt file is present.
    !--------------------------------------------------------------------
    !check file exist to trigger the LJ-matrix update
    inquire(file="ljmatrix.mods.txt", exist=exist)
    !if true, than read until eof has been reached
    if (exist) then
       open(in_unit, file="ljmatrix.mods.txt", status="old", action="read")
       do
         read (in_unit,'(2i6,2e16.8)',IOSTAT=status) iu,iv,sm,em
         if (status /= 0) exit
         !write (6,'(a5,2i6,2e16.8)') "mods>", iu,iv,sm,em
         !call rism_report_message("(a20,2i6,2e16.8)", "lj-matrix-mods>",
         !iu,iv,sm,em)
         this%ljSigmaUV(iu,iv) = sm
         this%ljEpsilonUV(iu,iv) = em
       enddo
       call rism_report_message("|-> Read LJ-matrix modifications from ljmatrix.mods.txt.")


       ! Write the LJ-matrix (uv) to ljmatrix.orig.txt.
       ! ----------------------------------------------------------
       open (unit=out_unit,file="ljmatrix.orig.txt",action="write",&
             status="replace")
       do iv = 1, this%solvent%numAtomTypes
          do iu = 1, this%solute%numAtoms
             write (out_unit,'(2i6,2e16.8)') iu,iv,this%ljSigmaUV(iu,iv), &
                                                   this%ljEpsilonUV(iu,iv)
          end do
       end do
       close (out_unit)
    endif
    flush(6)
  end subroutine mixSoluteSolventLJParameters


  !> Tabulate the solute-solvent (UV) 12-6 Lennard-Jones potential on
  !! box grid points. The potential is subject to a radial cutoff
  !! distance.
  !!TODO: Merge with minimum image variation.
  !! @param[in] this potential object.
  !! @param[in,out] ulj Grid to add potential to.
  subroutine uvLennardJonesPotentialWithCutoff(this, ulj, ljTolerance)
    use rism_util, only : checksum
    use constants, only : PI, KB
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(inout) :: ulj(:,:,:,:)
    _REAL_, intent(in) :: ljTolerance
    ! Amount to offset the grid in the z-axis (when using MPI).
    _REAL_ :: offset
    ! Number of gridpoints within cutoff for each dimension.
    integer :: cutoff(3)
    ! Closest grid point to solute atom center.
    integer :: closestPoint(3)
    ! First and last gridpoints which are within the cutoff distance
    ! from the solute atom while remaining in the box.
    integer :: first(3), last(3)
    ! Grid, solute, solvent, and dimension indices.
    integer ::  igx, igy, igz, ig, iu, iv, id
    ! Grid point position.
    _REAL_ ::  rx, ry, rz
    ! Distance of grid point from solute.
    _REAL_ ::  dx2, dy2, dz2, dy2dz2, r2
    ! LJ site term and its sixth power inverse.
    _REAL_ ::  rs2, rs2i, rs6i
    ! Minimum LJ term r / sigma and its square.
    _REAL_ ::  rcor, rcor2, rcor2i
    parameter (rcor = 0.2d0, rcor2 = rcor**2, rcor2i = 1d0/rcor2)
    logical :: warnNoCorrection
    ! No cutoff optimized code is about 30-40% faster when the cutoff is larger than the box.
    if (all(this%ljCutoffs2 .eq. HUGE(1d0))) then
       this%applyLJCorrection = .false.
       call uvLennardJonesPotentialWithoutCutoff(this, ulj)
       return
    end if

    warnNoCorrection = .False.
    
    ! ! Calculate the number of grid points necessary to cover this
    ! ! cutoff radius.  In case of a large cutoff, we need to protect
    ! ! against integer overflow.
    ! cutoff = nint(min(dble(huge(1) - 1), sqrt(this%cutoff2) / this%grid%spacing)) + 1
    ! Offset used for dividing the grid into slabs along the z-axis
    ! for MPI.
    offset = this%grid%spacing(3) * this%grid%offsetR(3)


    ! possible alternative for no cutoff LJ potential (untested)
    
    ! rx_v => safemem_reallocate(r, this%grid%localDimsR(1))
    ! ry_v => safemem_reallocate(r, this%grid%localDimsR(2))
    ! rz_v => safemem_reallocate(r, this%grid%localDimsR(3))

    ! do igx = 1, len(rx_v)
    !    rx_v(igx) = (igx-1) * this%grid%grdspc(1)
    ! end do
    ! do igy = 1, len(ry_v)
    !    ry_v(igy) = (igy-1) * this%grid%grdspc(2)
    ! end do
    ! do igz = 1, len(rz_v)
    !    rz_v(igx) = (igz-1) * this%grid%grdspc(3)    ! end do
    ! do iu = 1, this%solute%numAtoms
    !    rx_v = rx_v - this%solute%position(1, iu)
    !    ry_v = ry_v - this%solute%position(2, iu)
    !    rz_v = rz_v - this%solute%position(3, iu)
    !    do igz = 1, len(rz_v)
    !       do igy = 1, len(ry_v)
    !          do igx = 1, len(rx_v)
    !             r2(igx, igy, igx) = rx_v*rx_v + ry_v*ry_v + rz_v*rz_v
    !          end do
    !       end do
    !    end do
    !    r2 = 1d0/r2*r2*r2
    !    do iv = 1, this%solvent%numAtomTypes
    !        ulj(:,:,:,iu) = ulj(:,:,:,iu) + r2*&
    !(this%ljAUV(iu,iv)*r2- this%ljBUV(iu,iv))
    !    end do
    !    ulj(:,:,:,iu) = 
    !    rx_v = rx_v + this%solute%position(1, iu)
    !    ry_v = ry_v + this%solute%position(2, iu)
    !    rz_v = rz_v + this%solute%position(3, iu)
    ! end do

    ! Main optimizations are
    ! 1. Expanding exponentials into multiplications and avoiding all
    !     but one division
    ! 2. Moving the calculation of 1/r^6 outside of the inner loop
    ! 3. Precomputed LJ A and B coefficients.
    ! 4. loop over solute atoms is outside the grid point loop (not
    !     sure why this is faster)
    ! 5. r6 and r12 cutoffs (need to test impact)
    do iv = 1, this%solvent%numAtomTypes
       do iu = 1, this%solute%numAtoms
          ! Calculate the number of grid points necessary to cover this
          ! cutoff radius.  In case of a large cutoff, we need to protect
          ! against integer overflow.
          ! cutoff = nint(min(dble(huge(1) - 1), sqrt(maxval(this%ljCutoffs2)) / this%grid%spacing)) + 1
          ! write(0,*) 'box cut', sqrt(cutoff), sqrt(maxval(this%ljCutoffs2))
          ! Find closest grid point near solute atom center.
          ! do id = 1, 3
          !    closestPoint(id) = nint(this%solute%position(id, iu) / this%grid%spacing(id))
          ! end do
          ! Apply z-axis offset for MPI (if applicable).
          ! closestPoint(3) = closestPoint(3) - this%grid%offsetR(3)

          ! Find the first and last grid points which are within the
          ! cubic cutoff distance while remaining within the box.
          do id = 1, 3
             ! rounddown
             first(id) = floor( ( this%solute%position(id, iu) - sqrt(this%ljCutoffs2(iu,iv)) ) &
                  / this%grid%spacing(id)) + 1
             ! rounddown
             last(id) = floor( ( this%solute%position(id, iu) + sqrt(this%ljCutoffs2(iu,iv)) ) &
                  / this%grid%spacing(id)) 
             if(first(id) .lt. 1 .or. last(id) .gt. this%grid%globalDimsR(id)) then
                warnNoCorrection = .true.
             end if
             ! ensure we only access gridpoint inside the box
             first(id) = max(first(id), 1+this%grid%offsetR(id)) - this%grid%offsetR(id)
             last(id) =  min(last(id),  this%grid%offsetR(id)  + this%grid%localDimsR(id))&
                  -this%grid%offsetR(id)
          end do
          ! Calculate the solute-solvent Lennard-Jones potential at each
          ! point within the cutoff radius.
          do igz = first(3), last(3)
             rz = (igz - 1) * this%grid%spacing(3) + offset
             dz2 = (rz - this%solute%position(3, iu))*(rz - this%solute%position(3, iu))
             do igy = first(2), last(2)
                ry = (igy - 1) * this%grid%spacing(2)
                dy2 = (ry - this%solute%position(2, iu))*(ry - this%solute%position(2, iu))
                dy2dz2 = dz2 + dy2
                do igx = first(1), last(1)
                   rx = (igx - 1) * this%grid%spacing(1)

                   ! Find distance from grid point to solute site in real
                   ! space.
                   dx2 = (rx - this%solute%position(1, iu))*(rx - this%solute%position(1, iu))
                   r2 = dx2 + dy2dz2
                   ! Ensure grid point is within radial cutoff distance
                   ! from solute atom. This is necessary in order to
                   ! enforce a radial cutoff since 'first' and 'last'
                   ! grid points defined above treated cutoff as cubic.
                   if (r2 < this%ljCutoffs2(iu,iv)) then
                      rs2i = 1d0/r2
                      rs6i = (rs2i*rs2i*rs2i)
                      ulj(igx, igy, igz, iv) = ulj(igx, igy, igz, iv) &
                           + rs6i * &
                           (this%ljAUV(iu,iv)*rs6i - this%ljBUV(iu,iv))
                   end if
                end do
             end do
          end do
       end do
    end do
    if (warnNoCorrection) then
       this%applyLJCorrection = .false.
       call rism_report_warn(&
            'LJ tolerance extends beyond the solvent box and the LJ '//NEW_LINE('a')//&
            'correction will not be used. For more accurate calculations,'//NEW_LINE('a')//&
            'increase the tolerance, box dimensions, or use buffer=0')
    end if
    ! Enforce minimum LJ term value, rcor.
    where (abs(ulj) > HUGE(1d0) .or. ulj /= ulj)
       ulj = sqrt(HUGE(1d0))
    end where
  end subroutine uvLennardJonesPotentialWithCutoff


  !> Tabulate the solute-solvent (UV) 12-6 Lennard-Jones potential on
  !! box grid points. No cutoff is applied, making it reasonably fast.
  !! @param[in] this potential object.
  !! @param[in,out] ulj Grid to add potential to.
  subroutine uvLennardJonesPotentialWithoutCutoff(this, ulj)
    implicit none
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(inout) :: ulj(:,:,:,:)
    ! Amount to offset the grid in the z-axis (when using MPI).
    _REAL_ :: offset
    ! Grid, solute, solvent, and dimension indices.
    integer ::  igx, igy, igz,  iu, iv
    ! Grid point position.
    _REAL_ ::  rx, ry, rz
    ! Distance of grid point from solute.
    _REAL_ ::  dx2, dy2, dz2, dy2dz2, r2
    ! LJ site term and its sixth power inverse.
    _REAL_ ::  rs2, rs2i, rs6i

    ! Offset used for dividing the grid into slabs along the z-axis
    ! for MPI.
    offset = this%grid%spacing(3) * this%grid%offsetR(3)
    do iv = 1, this%solvent%numAtomTypes
       do iu = 1, this%solute%numAtoms
       ! Calculate the solute-solvent Lennard-Jones potential at each
       ! point within the cutoff radius.
       do igz = 1, this%grid%localDimsR(3)
          rz = (igz - 1) * this%grid%spacing(3) + offset
          dz2 = (rz - this%solute%position(3, iu))
          dz2 = dz2*dz2
          do igy = 1, this%grid%localDimsR(2)
             ry = (igy - 1) * this%grid%spacing(2)
             dy2 = (ry - this%solute%position(2, iu))
             dy2 = dy2*dy2
             dy2dz2 = dz2 + dy2
             do igx =  1, this%grid%localDimsR(1)
                rx = (igx - 1) * this%grid%spacing(1)

                ! Find distance from grid point to solute site in real
                ! space.
                dx2 = (rx - this%solute%position(1, iu))
                dx2 = dx2*dx2
                r2 = dx2 + dy2 + dz2
                ! reciprocal is the most expensive operation here
                rs2i = 1d0/r2
                rs6i = (rs2i*rs2i*rs2i)
                ulj(igx, igy, igz, iv) = ulj(igx, igy, igz, iv) &
                     + rs6i * &
                     (this%ljAUV(iu,iv)*rs6i - this%ljBUV(iu,iv))
                end do
             end do
          end do
       end do
    end do
    ! Enforce minimum LJ term value, rcor.
    where (abs(ulj) > HUGE(1d0) .or. ulj /= ulj)
       ulj = sqrt(HUGE(1d0))
    end where
  end subroutine uvLennardJonesPotentialWithoutCutoff

  !> Tabulate the solute-solvent 12-6 Lennard-Jones potential in the
  !! box subject to the minimum image convention.
  !! @param[in] this potential object
  !! @param[in,out] ulj grid to add potential to
  subroutine uvLennardJonesPotentialWithMinimumImage(this, ulj)
    use constants, only : PI
    use rism_util, only : checksum
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(inout) :: ulj(:,:,:,:)
    ! Amount to offset the grid in the z-axis (when using MPI).
    _REAL_ :: offset
    ! Grid, solute, slovent, and dimension indices.
    integer :: igx, igy, igz, ig, iu, iv, id, cnt
    ! Grid point position.
    _REAL_ :: rx(3), ry(3), rz(3), gridPoint(3)
    ! Distance of grid point from solute.
    _REAL_ :: soluteDistanceSquared
    ! Distance projected along unit cell z-axis.
    _REAL_ :: cellZ
    ! Base term in LJ equation (ratio of sigma and distance).
    _REAL_ :: ljBaseTerm
    ! Minimum LJ term r / sigma and its square.
    _REAL_ :: minLJBaseTerm
    parameter (minLJBaseTerm = 0.2d0**2)
    _REAL_, parameter :: minDistance = 0.002d0
    _REAL_, parameter :: minDistanceSquared = minDistance**2
    _REAL_ :: solutePosition(3)

    ! Offset used for dividing the box along the z-axis
    ! for MPI.
    !TODO:JAJ Need to adjust for MPI z-slab decomposition.
    offset = this%grid%spacing(3) * this%grid%offsetR(3)
    
    do iu = 1, this%solute%numAtoms
       ! Calculate the solute-solvent Lennard-Jones potential at each
       ! grid point in the box.
       do igz = 1, this%grid%localDimsR(3)
          rz = (igz - 1 + this%grid%offsetR(3)) * this%grid%voxelVectorsR(3, :)
          do igy = 1, this%grid%localDimsR(2)
             ry = (igy - 1) * this%grid%voxelVectorsR(2, :)
             do igx = 1, this%grid%localDimsR(1)
                rx = (igx - 1) * this%grid%voxelVectorsR(1, :)
                gridPoint = rx + ry + rz

                ! Vector from grid point to solute atom, in real space.
                solutePosition = gridPoint - this%solute%position(:, iu)

                ! Apply minimum image convention.
                solutePosition = minimumImage(this, solutePosition)

                soluteDistanceSquared = dot_product(solutePosition, solutePosition)
                if (soluteDistanceSquared < this%cutoff2) then
                   soluteDistanceSquared = max(soluteDistanceSquared, &
                                               minDistanceSquared )
                   do iv = 1, this%solvent%numAtomTypes
                      ljBaseTerm = soluteDistanceSquared / this%ljSigmaUV(iu, iv)**2
                      ljBaseTerm = 1d0 / ljBaseTerm**3
                      ulj(igx, igy, igz, iv) = ulj(igx, igy, igz, iv) &
                           + this%ljEpsilonUV(iu, iv) * ljBaseTerm * (ljBaseTerm - 2.d0)
                   end do
                end if
             end do
          end do
       end do
    end do
  end subroutine uvLennardJonesPotentialWithMinimumImage


  !> Tabulate the solute-solvent 12-6 Lennard-Jones potential, and
  !! the short-range Ewald electrostatic potential,  in the
  !! box subject to the minimum image convention.
  !! @param[in] this potential object
  !! @param[in,out] ulj grid to add potential to
  subroutine uvLJrEwaldPotentialWithMinimumImage(this, ulj)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(inout) :: ulj(:,:,:,:)
    ! Grid, solute, slovent, and dimension indices.
    integer :: igx, igy, igz, iu, iv
    ! Grid point position.
    _REAL_ :: rx(3), ry(3), rz(3), gridp(3)
    ! Distance of grid point from solute.
    _REAL_ :: sd2
    ! Base term in LJ equation (ratio of sigma and distance).
    _REAL_ :: ljBaseTerm(this%solvent%numAtomTypes)
    _REAL_ :: solutePosition(3)
    _REAL_ :: sd, sr
    _REAL_ :: sigma(this%solute%numAtoms,this%solvent%numAtomTypes), beta

    beta = 1.d0/this%chargeSmear
    do iu = 1, this%solute%numAtoms
       do iv = 1, this%solvent%numAtomTypes
          sigma(iu,iv) = 1.d0/this%ljSigmaUV(iu, iv)**2
       end do
    end do

#ifndef MPI
!$omp parallel do private (rx,ry,rz,solutePosition,sd2,sd,sr,ljBaseTerm) 
#endif
    do igz = 1, this%grid%localDimsR(3)
       rz = (igz - 1 + this%grid%offsetR(3)) * this%grid%voxelVectorsR(3, :)
       do igy = 1, this%grid%localDimsR(2)
          ry = (igy - 1) * this%grid%voxelVectorsR(2, :)
          do igx = 1, this%grid%localDimsR(1)
             rx = (igx - 1) * this%grid%voxelVectorsR(1, :)

             do iu = 1, this%solute%numAtoms

                solutePosition = rx + ry + rz - this%solute%position(:, iu)
                solutePosition = minimumImage(this, solutePosition)

                sd2 = max(4d-6, dot_product(solutePosition, solutePosition))
                if (sd2 < this%cutoff2) then

                   sd = sqrt(sd2)
                   sr = this%solute%charge(iu) * erfc(sd * beta) / sd
                   ljBaseTerm(:) = 1.d0 / (sd2 * sigma(iu,:))**3

                   ulj(igx,igy,igz,:) = ulj(igx,igy,igz,:) &
                           + this%ljEpsilonUV(iu, :) &
                             * ljBaseTerm(:) * (ljBaseTerm(:) - 2.d0) &
                           + sr * this%solvent%charge(:)
                end if
             end do
          end do
       end do
    end do
#ifndef MPI
!$omp end parallel do
#endif
  end subroutine uvLJrEwaldPotentialWithMinimumImage

  !> Synthesizing the Coulomb solute potential in a sparser box.
  !! The exact electrostatic potential is calculated for on a sparse
  !! grid (1/8 the number of points or 1/2 density) and for all grid
  !! points within cutoff of an atom.  To ensure that all interpolated
  !! points are enclosed by eight neighbours, the sparse grid is
  !! extended by one index in each dimension.
  !! @param[in] this potential object
  !! @param[in,out] ucu grid to add potential to
  subroutine uvCoulombicPotential(this, ucu)
    use rism_util, only : checksum
    use safemem
    implicit none
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(inout) :: ucu(:,:,:,:)

    ! For MPI spatial decomposition.
    ! offset :: distance in z-axis due to spatial decomp
    _REAL_ :: offset
    integer :: cutoff(3), closestPoint(3), id
    integer :: first(3), last(3)

    ! Temporary / intermediate variables.
    integer :: igx, igy, igz, ig, iu, iv
    _REAL_ :: r2, rx, ry, rz, dx2, dy2, dy2dz2, dz2
#ifdef RISM_DEBUG
    write(6,*) 'Tabulating fast solute Coulomb potential ...'
    call flush(6)
#endif

    ! Offset used for dividing the grid into slabs along the z-axis
    ! for MPI.
    offset = this%grid%spacing(3) * this%grid%offsetR(3)
    do igz = 1, this%grid%localDimsR(3)
       rz = (igz - 1) * this%grid%spacing(3) + offset
       do igy = 1, this%grid%localDimsR(2)
          ry = (igy - 1) * this%grid%spacing(2)
          do igx =  1, this%grid%localDimsR(1)
             rx = (igx - 1) * this%grid%spacing(1)
             do iu = 1, this%solute%numAtoms
                dz2 = (rz - this%solute%position(3, iu))
                dz2 = dz2*dz2
                dy2 = (ry - this%solute%position(2, iu))
                dy2 = dy2*dy2
                dy2dz2 = dz2 + dy2

                ! Find distance from grid point to solute site in real
                ! space.
                dx2 = (rx - this%solute%position(1, iu))
                dx2 = dx2*dx2
                r2 = dx2 + dy2 + dz2
                ! reciprocal is the most expensive operation here
                ucu(igx, igy, igz, 1) = ucu(igx, igy, igz, 1) + this%solute%charge(iu)/sqrt(r2)
             end do
          end do
       end do
    end do
    do iv = this%solvent%numAtomTypes, 1, -1
       ucu(:,:,:,iv) =  ucu(:,:,:,1) * this%solvent%charge(iv)
    end do
    ! Remove singularities (inf)
    where (abs(ucu) > HUGE(1d0) .or. ucu /= ucu)
       ucu = sqrt(HUGE(1d0))
    end where
  end subroutine uvCoulombicPotential


  !> Calculate the Ewald sum electric potential between the solute and
  !! solvent.  The minimum image convention is used to efficiently
  !! approximate infinitely periodic solute.
  !! References:
  !!   Understanding Molecular Simulations 2E, chapter on long range
  !!   interactions.
  !!
  !!   Lipkowitz, Kenny B., and Thomas R. Cundari, eds.
  !!   “Appendix F: Mathematical Aspects of Ewald Summation.”
  !!   Reviews in Computational Chemistry, 447–78. 2007.
  !!   doi/10.1002/9780470164112.app6
  !! @param[in] this Potential object.
  !! @param[in,out] ucu Grid to add potential to.
  subroutine uvEwaldSumPotential(this, ucu)
    use rism3d_fft_c
    use rism_util, only : checksum
    use safemem
    use constants, only : PI, FOURPI
    use rism3d_opendx, only : rism3d_opendx_write
    implicit none
#ifdef FFW_THREADS
    integer :: nthreads, totthreads
    integer, external :: OMP_get_max_threads, OMP_get_num_threads
    logical, external :: OMP_get_dynamic, OMP_get_nested
#endif
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(inout) :: ucu(:,:,:,:)

    ! Flag indicating that this index of the array has not been set.
    _REAL_ :: untouched = HUGE(0d0)

    ! Temporary / intermediate variables.
    integer :: igx, igy, igz
    integer :: iv
    integer :: ig

    _REAL_ :: total
    _REAL_ :: chargeCorrection


    call timer_start(TIME_UCOULULR)
    call rism_timer_start(this%coulombLongRangeTimer)

    ! Clearing potentials summators.
    !TODO: Make this debug-only.
    ucu(:,:,:,1) = untouched
    this%uuv1d(:,1) = untouched

    ! Long range (k-space) term in Ewald sum.

    ! In order to cancel the surface charge introduced by ionic or
    ! polar periodic unit cells, the crystal is assumed to be immersed
    ! in a medium of infinite dielectric constant (perfect conductor)
    ! which produces an infinitely large reaction field. This
    ! limiting case of infinite dielectric results in the k = 0
    ! contibution to the potential energy always being zero. See
    ! Understanding Molecular Simulation 2E, section 12.1.3.
    ! Note: Unclear where this assumption originated or how it is
    ! physically justified.
    do iv = 1, this%solvent%numAtomTypes
       this%uuv1d(:,iv) = this%dcfLongRangeAsympK(:)
    end do

    ! Convert from reciprocal space to Cartesian space.
#ifdef MPI
    this%uuv1d(2:this%grid%totalLocalPointsK:2, :) = &
         -this%uuv1d(2:this%grid%totalLocalPointsK:2, :)
    call rism3d_fft_bwd(this%fft, this%uuv1d)
#else
    call rism3d_fft_bwd(this%fft, this%uuv1d)
#endif

    ! Copy the long-range Gaussian charge portion of the Ewald sum to
    ! to electric potential.
    do igz = 0, this%grid%localDimsR(3) - 1
       do igy = 0, this%grid%localDimsR(2) - 1
          do igx = 0, this%grid%localDimsR(1) - 1
             ig = 1 + igx + (igy + igz * this%grid%localDimsR(2)) &
                  * this%grid%localDimsR(1)
             ucu(igx + 1, igy + 1, igz + 1, 1) = this%uuv1d(ig, 1)
          end do
       end do
    end do

    call rism_timer_stop(this%coulombLongRangeTimer)
    call timer_stop(TIME_UCOULULR)

    call timer_start(TIME_UCOULUSR)
    call rism_timer_start(this%coulombShortRangeTimer)

    ! Calculate short-range term of Ewald sum and combine with
    ! previously calculated long-range term.
    call uvEwaldSumShortRangePotential(this, ucu)

    call rism_timer_stop(this%coulombShortRangeTimer)
    call timer_stop(TIME_UCOULUSR)

    ! No correction for self-interaction of a point charge with its
    ! own Ewald gaussian charge is needed. Correction is only needed
    ! for intramolecular atomic interactions.  Solute-solvent is
    ! purely intermolecular.  See p. 298 of Understanding Molecular
    ! Simulation 2E by Frenkel & Smit (ISBN-13: 978-0-12-267351-1)
    ! for details.

    ! Additional term for charged systems
    chargeCorrection = - pi  * (this%solute%totalCharge / this%grid%boxVolume) &
                        *  this%chargeSmear * this%chargeSmear

    write (6,*) "Ewald charge correction: ", chargeCorrection
    ucu(:,:,:,1) = ucu(:,:,:,1) + chargeCorrection


    ! Calculate electrostatic potential energy on each grid point.
    do iv = this%solvent%numAtomTypes, 1, -1
          ucu(:,:,:,iv) = ucu(:,:,:,1) * this%solvent%charge(iv)
    end do
  end subroutine uvEwaldSumPotential

  !> Short-range portion of the Ewald sum electric potential.
  subroutine uvEwaldSumShortRangePotential(this, ucu)
    use constants, only: pi
    use rism3d_opendx, only : rism3d_opendx_write
    implicit none
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(inout) :: ucu(:,:,:,:)

    ! Minimum distance to prevent division by zero.
    _REAL_, parameter :: minDistance = 0.001
    _REAL_, parameter :: minDistance2 = minDistance**2
    integer :: igx, igy, igz
    integer :: ig
    integer :: iu, cnt

    _REAL_ :: rx(3), ry(3), rz(3)
    _REAL_ :: gridPoint(3), solutePosition(3)
    _REAL_ :: soluteDistance, soluteDistance2
    _REAL_ :: sr


    ! Calculate short-range term of Ewald sum and combine with
    ! previously calculated long-range term.

    do iu = 1, this%solute%numAtoms
       do igz = 1, this%grid%localDimsR(3)
          rz = (igz - 1 + this%grid%offsetR(3)) * this%grid%voxelVectorsR(3, :)
          do igy = 1, this%grid%localDimsR(2)
             ry = (igy - 1) * this%grid%voxelVectorsR(2, :)
             do igx = 1, this%grid%localDimsR(1)
                rx = (igx - 1) * this%grid%voxelVectorsR(1, :)
                gridPoint = rx + ry + rz

                ! Solute atom position relative to a grid point.
                solutePosition = gridPoint - this%solute%position(:, iu)

                ! Apply minimum image convention.
                solutePosition = minimumImage(this, solutePosition)

                ! Distance from solute atom to gridpoint.
                !soluteDistance = sqrt(dot_product(solutePosition, solutePosition))
                soluteDistance2 = dot_product(solutePosition, solutePosition)

                !if (soluteDistance2 < minDistance2) soluteDistance = minDistance
                if (soluteDistance2 < minDistance2) soluteDistance2 = minDistance2

                ! Short-range term of Ewald sum.
                !if (soluteDistance < this%cutoff) then
                if (soluteDistance2 < this%cutoff2) then
                   soluteDistance = sqrt(soluteDistance2)
                   sr = this%solute%charge(iu) &
                           * erfc(soluteDistance / this%chargeSmear) / soluteDistance
                   ucu(igx , igy , igz , 1) = ucu(igx , igy , igz , 1) + sr
                end if

             end do
          end do
       end do
    end do
  end subroutine uvEwaldSumShortRangePotential


  !> Long-range portion of the Particle Mesh Ewald (PME) electric potential.
  subroutine uvParticleMeshRecipEwaldPotential(this, ucu)
    use, intrinsic :: iso_c_binding
    use bspline
    use constants, only : pi
    use FFTW3
    use rism3d_opendx, only : rism3d_opendx_write
    use rism_util, only: r2c_pointer
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(inout) :: ucu(:,:,:,:)

    ! Ewald charge smear parameter.
    ! _REAL_, parameter :: smear = this%chargeSmear
    ! Order of b-spline interpolation.
    integer, parameter :: splineOrder = 6

    !_REAL_, parameter :: zeta = this%chargeSmear * this%chargeSmear
    !_REAL_, parameter :: t = -0.25 / zeta

    _REAL_ :: smear, zeta, t

    logical, parameter :: DEBUG = .false.

    integer :: i, j, k
    integer :: ix, iy, iz, ixyz
    integer :: igx, igy, igz
    integer :: igk, igk_fort, igk_cpp, igk_fort_global, igk_cpp_global
    integer :: id
    integer :: iu, iv

    integer :: gridPoints(splineOrder, 3)

    _REAL_, pointer :: kxi(:,:) => NULL()
    _REAL_, pointer :: kyi(:,:) => NULL()
    _REAL_, pointer :: kzi(:,:) => NULL()
    _REAL_ :: k2
    _REAL_ :: waveVector(3)
    integer :: lgx, lgy, lgz

    _REAL_ :: byz, bxyz

    _REAL_ :: weights(splineOrder, 3)
    _REAL_, pointer :: bsplineFourierCoeffX(:) => NULL()
    _REAL_, pointer :: bsplineFourierCoeffY(:) => NULL()
    _REAL_, pointer :: bsplineFourierCoeffZ(:) => NULL()
    _REAL_, pointer :: gaussianFourierCoeff(:) => NULL()
    _REAL_, pointer :: kernel(:) => NULL()

    type(C_PTR) :: uuv1d_r_cptr, uuv1d_c_cptr, uuv1d_final_cptr
    real(C_DOUBLE), pointer :: uuv1d_r(:,:,:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: uuv1d_c(:,:,:) => NULL()
    real(C_DOUBLE), pointer :: uuv1d_final(:,:,:) => NULL()
    integer(C_INTPTR_T) :: localPtsK
    real(C_DOUBLE), pointer :: outr_1d(:) => NULL()
    integer(C_INTPTR_T) :: L, M, N
    integer(C_INTPTR_T) :: local_N, local_k_offset

    type(C_PTR) :: planfwd = C_NULL_PTR, planbwd = C_NULL_PTR

    _REAL_ :: reciprocalPos(3)
    _REAL_ :: chargeCorrection
    ! integer :: numGridPointsK

    integer :: gridDimX_k

    character(len=120) :: filename, suffix
    integer :: ierr

    smear = this%chargeSmear
    zeta = this%chargeSmear * this%chargeSmear
    t = -0.25 * zeta

    L = this%grid%localDimsR(1)
    M = this%grid%localDimsR(2)
    N = this%grid%globalDimsR(3)

!    call timer_start(TIME_UCOULULR)
!    call rism_timer_start(this%coulombLongRangeTimer)

    gridDimX_k = this%grid%localDimsR(1) / 2 + 1
    kxi => safemem_realloc(kxi,  gridDimX_k, 3, .false.)
    kyi => safemem_realloc(kyi,  this%grid%localDimsR(2), 3, .false.)
    kzi => safemem_realloc(kzi,  this%grid%localDimsR(3), 3, .false.)

    bsplineFourierCoeffX => safemem_realloc(bsplineFourierCoeffX, this%grid%globalDimsR(1), .false.)
    bsplineFourierCoeffY => safemem_realloc(bsplineFourierCoeffY, this%grid%globalDimsR(2), .false.)
    bsplineFourierCoeffZ => safemem_realloc(bsplineFourierCoeffZ, this%grid%globalDimsR(3), .false.)

    gaussianFourierCoeff => safemem_realloc(gaussianFourierCoeff, &
         gridDimX_k * this%grid%localDimsR(2) * this%grid%localDimsR(3), .false.)

#ifdef MPI
    ! localPtsK = fftw_mpi_local_size_3d_transposed( &
    localPtsK = fftw_mpi_local_size_3d(N, M, L / 2 + 1, &
         this%grid%mpicomm, local_N, local_k_offset)
         ! localPtsK_y, localOffsetK_y)
#else
    localPtsK = gridDimX_k * this%grid%localDimsR(2) * this%grid%localDimsR(3)
    local_N = this%grid%localDimsR(3)
    local_k_offset = 0
#endif

    uuv1d_r_cptr = fftw_alloc_real(2 * localPtsK)
    uuv1d_c_cptr = fftw_alloc_complex(localPtsK)
#ifdef MPI
    call c_f_pointer(uuv1d_r_cptr, uuv1d_r, [2*(L/2+1), M, local_N])
    call c_f_pointer(uuv1d_r_cptr, outr_1d, [2*(L/2+1) * M * local_N])
    call c_f_pointer(uuv1d_c_cptr, uuv1d_c, [L/2+1, M, local_N])
    if (this%grid%mpirank == 0) then
       uuv1d_final_cptr = fftw_alloc_real(L * M * N)
       call c_f_pointer(uuv1d_final_cptr, uuv1d_final, [L, M, N])
    end if
#else
    call c_f_pointer(uuv1d_r_cptr, uuv1d_r, [L, M, local_N])
    call c_f_pointer(uuv1d_r_cptr, outr_1d, [L * M * local_N])
    call c_f_pointer(uuv1d_c_cptr, uuv1d_c, [L/2+1, M, local_N])
#endif

    kernel => safemem_realloc(kernel, &
         gridDimX_k * this%grid%localDimsR(2) * this%grid%localDimsR(3), .false.)

#if 0
    if (any(mod(this%grid%localDimsR, 2) > 0)) then
       !TODO: call RISM error functions
#ifdef MPI
       if( this%grid%mpirank == 0 ) &
#endif
       write(6,'(a,a)')  "| PME implementation prefers an even", &
            " number of grid points on each axis."
    end if
#endif

   ! Angular wave numbers.
   !FIXME: This needs to be removed and X should be the halved axis.
   !FIXME: Scrutinize grid dims in loops (K v. R, local vs. global).
   do igk = 0, gridDimX_k - 1
      kxi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(1,:) * igk
   end do
   do igk = 0, this%grid%localDimsR(2) - 1
      kyi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(2,:) &
           * merge(igk, -(this%grid%localDimsR(2) - igk), &
           igk < this%grid%localDimsR(2) / 2 + 1)
   end do
   do igk = 0, this%grid%localDimsR(3) - 1
      kzi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(3,:) &
           * merge(igk + this%grid%offsetR(3), &
           -(this%grid%globalDimsR(3) - (igk + this%grid%offsetR(3))), &
           (igk + this%grid%offsetR(3)) < this%grid%globalDimsR(3) / 2 + 1)
   end do

    ! Compute the discrete Fourier transform coefficients of the b-spline.
    call cardinal_bspline(merge(0d0, 0.5d0, mod(splineOrder, 2) == 0), &
         splineOrder, weights(:,1))
    !TODO: For triclinic case, b-spline may have axial interdependence
    ! and thus so will its Fourier coefficients.

    !FIXME: Make these functions MPI aware by passing in global +
    ! local dims and offset. For now, just compute all coeffs on every
    ! node and use the grid offset to index.
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%globalDimsR(1), splineOrder, &
         weights(:,1), bsplineFourierCoeffX, .true.)
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%globalDimsR(2), splineOrder, &
         weights(:,1), bsplineFourierCoeffY, .false.)
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%globalDimsR(3), splineOrder, &
         weights(:,1), bsplineFourierCoeffZ, .false.)

    !TODO: Much of this code can be linearized like the C++ counterpart.

    ! Reciprocal space kernel with b-spline discrete Fourier transform
    ! correction.
    !TODO: If b-spline Fourier coefficients are stored in a linearized
    ! array, then this could be simplified to a single loop.
    ixyz = 1
    do iz = 0, this%grid%localDimsR(3) - 1
       do iy = 0, this%grid%localDimsR(2) - 1
          byz = bsplineFourierCoeffY(iy + 1) * bsplineFourierCoeffZ(iz + 1 + this%grid%offsetR(3))
          do ix = 0, gridDimX_k - 1

!#if defined(MPI)
!             igk = 1 + ix + &
!                  iy * (this%grid%localDimsK(1) / 2) + &
!                  iz * this%grid%localDimsK(2) * (this%grid%localDimsK(1) / 2)
!#else
!
!             igk = 1 + ix + iy * (this%grid%localDimsR(1) / 2 ) +&
!                    iz * this%grid%localDimsR(2) * (this%grid%localDimsR(1) / 2 )
!             if (ix .eq. gridDimX_k - 1) then
!                 igk = 1 + iy + iz * this%grid%globalDimsR(2) &
!               + this%grid%globalDimsR(1) / 2 * this%grid%globalDimsR(2) * this%grid%globalDimsR(3)
!             end if
!#endif
             bxyz = bsplineFourierCoeffX(ix + 1) * byz
             waveVector = kxi(ix+1,:) + kyi(iy+1,:) + kzi(iz+1,:)

             k2 = dot_product(waveVector, waveVector)
             !k2 = this%grid%waveVectors2(igk)
!            if ( abs(k2-this%grid%waveVectors2(igk)) .gt. 1e-8   ) then
!               write (6,*) "-|: ", iz,iy,ix,ixyz,igk,k2-this%grid%waveVectors2(igk), this%grid%mpirank
!            end if
             kernel(ixyz) = (4 * pi / k2) * exp(k2 * t) / bxyz
             ixyz = ixyz + 1
          end do
       end do
    end do

    if (this%grid%offsetR(3) == 0) then
       ! Remove the k = 0 term (tinfoil boundary conditions).
       kernel(1) = 0d0
    end if

    ! Allocate the FFT.
    !TODO: Consider switching to in-place, transposed in/out FFT.
    ! Speed vs. memory usage tradeoff.
#ifdef MPI
    planfwd = fftw_mpi_plan_dft_r2c_3d(N, M, L, &
         uuv1d_r, uuv1d_c, this%grid%mpicomm, &
         !TODO: Enable FFTW_MPI_TRANSPOSED_OUT for speedup.
         ! ior(ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT), FFT_ALIGNED))
         ! ior(FFTW_ESTIMATE, FFTW_MPI_TRANSPOSED_OUT))
         FFTW_ESTIMATE)

    planbwd = fftw_mpi_plan_dft_c2r_3d(N, M, L, &
         uuv1d_c, uuv1d_r, this%grid%mpicomm, &
         ! ior(ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN), FFT_ALIGNED))
         FFTW_ESTIMATE)
#else
    planfwd = fftw_plan_dft_r2c_3d( &
         this%grid%globalDimsR(3), &
         this%grid%globalDimsR(2), &
         this%grid%globalDimsR(1), &
         uuv1d_r, uuv1d_c, FFTW_ESTIMATE)
    planbwd = fftw_plan_dft_c2r_3d( &
         this%grid%globalDimsR(3), &
         this%grid%globalDimsR(2), &
         this%grid%globalDimsR(1), &
         uuv1d_c, uuv1d_r, FFTW_ESTIMATE)
#endif

    ! Initialize grids.
    ucu = 0
    uuv1d_r = 0
    uuv1d_c = 0
    this%uuv1d = 0

    ! Spread the Gaussian charges to the grid using b-splines.
    do iu = 1, this%solute%numAtoms
#ifdef MPI
       !TODO: Perform a check earlier for whether the atom
       ! can contribute to this portion of the grid. If not, skip it.
#endif
#if 1
       ! Convert Cartesian position to reciprocal space by projecting
       ! to reciprocal unit cell vectors.
       do id = 1, 3
          reciprocalPos(id) = dot_product(this%solute%position(:,iu), &
               this%grid%unitCellVectorsK(id, :))
       end do

       ! Set the box length to unity since reciprocal space already
       ! divides positions by the box length.
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(1), 1d0, &
            this%grid%globalDimsR(1), &
            splineOrder, gridPoints(:,1), weights(:,1))
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(2), 1d0, &
            this%grid%globalDimsR(2), &
            splineOrder, gridPoints(:,2), weights(:,2))
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(3), 1d0, &
            this%grid%globalDimsR(3), &
            splineOrder, gridPoints(:,3), weights(:,3))
#else
       ! Assume orthorhombic unit cell.
       call cardinal_bspline_periodic_box( &
            this%solute%position(1,iu), this%grid%boxLength(1), &
            this%grid%globalDimsR(1), &
            splineOrder, weights(:,1), gridPoints(:,1))
       call cardinal_bspline_periodic_box( &
            this%solute%position(2,iu), this%grid%boxLength(2), &
            this%grid%globalDimsR(2), &
            splineOrder, weights(:,2), gridPoints(:,2))
       call cardinal_bspline_periodic_box( &
            this%solute%position(3,iu), this%grid%boxLength(3), &
            this%grid%globalDimsR(3), &
            splineOrder, weights(:,3), gridPoints(:,3))
#endif

       do k = 1, splineOrder
#ifdef MPI
          ! Check if grid point is within local grid.
          if (gridPoints(k,3) + 1 < this%grid%offsetR(3) + 1 &
               .or. gridPoints(k,3) + 1 > this%grid%offsetR(3) + this%grid%localDimsR(3)) then
             cycle
          end if
#endif
          do j = 1, splineOrder
             byz = weights(j,2) * weights(k,3)
             do i = 1, splineOrder
                uuv1d_r(1 + gridPoints(i,1), 1 + gridPoints(j,2), &
                     1 + gridPoints(k,3) - this%grid%offsetR(3)) = &
                     uuv1d_r(1 + gridPoints(i,1), 1 + gridPoints(j,2), &
                     1 + gridPoints(k,3) - this%grid%offsetR(3)) &
                     + weights(i,1) * byz * this%solute%charge(iu)
             end do
          end do
       end do
    end do


    ! Convert the charge density into reciprocal space.

#ifdef MPI
    call fftw_mpi_execute_dft_r2c(planfwd, uuv1d_r, uuv1d_c)
#else
    call fftw_execute_dft_r2c(planfwd, uuv1d_r, uuv1d_c)
#endif

    ! Convolute with the kernel.
    do igz = 0, local_N - 1
       do igy = 0, M - 1
          do igx = 0, (L / 2 + 1) - 1
             igk_fort = 1 + igx + (igy + igz * M) * (L / 2 + 1)
             uuv1d_c(igx + 1, igy + 1, igz + 1) &
                  = uuv1d_c(igx + 1, igy + 1, igz + 1) * kernel(igk_fort)
          end do
       end do
    end do

    ! Evaluate the recip-space potential at the grid points.
    ! outr = 0
    uuv1d_r = 0
#ifdef MPI
    call fftw_mpi_execute_dft_c2r(planbwd, uuv1d_c, uuv1d_r)
#else
    call fftw_execute_dft_c2r(planbwd, uuv1d_c, uuv1d_r)
#endif

    chargeCorrection = - pi * (this%solute%totalCharge/this%grid%boxVolume) &
            * zeta
    ! write (6,*) "PME charge correction: ", chargeCorrection, " with biasPotential: ", this%biasPotential
    do igz = 0, local_N - 1
       do igy = 0, M - 1
          do igx = 0, L - 1
             ucu(igx + 1, igy + 1, igz + 1, :) = uuv1d_r(igx + 1, igy + 1, igz + 1) &
                  / this%grid%boxVolume + chargeCorrection
          end do
       end do
    end do

!    call rism_timer_stop(this%coulombLongRangeTimer)
!    call timer_stop(TIME_UCOULULR)

    ! Deallocate the FFT plans.
    call fftw_destroy_plan(planfwd)
    call fftw_destroy_plan(planbwd)

    ! Calculate electrostatic potential energy on each grid point.
    ! No 1/2 term is required since solute affects solvent, but
    ! solvent does not affect solute, hence no double counting.
    do iv = this%solvent%numAtomTypes, 1, -1
          ucu(:,:,:,iv) = ucu(:,:,:,1) * this%solvent%charge(iv)
    end do

    !TODO: Still causes crashes in rare circumstances. Need to debug
    ! more thoroughly.
    if (safemem_dealloc(kxi) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(kyi) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(kzi) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    ! if (safemem_dealloc(k2i) /= 0) then
    !    call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    ! end if
    if (safemem_dealloc(bsplineFourierCoeffX) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(bsplineFourierCoeffY) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(bsplineFourierCoeffZ) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(gaussianFourierCoeff) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    call fftw_free(uuv1d_r_cptr)
    call fftw_free(uuv1d_c_cptr)
#ifdef MPI
    if (this%grid%mpirank == 0) then
       call fftw_free(uuv1d_final_cptr)
    end if
#endif
    if (safemem_dealloc(kernel) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
  end subroutine uvParticleMeshRecipEwaldPotential

  !> Returns the value of the k-space expression for the Lennard-Jones potential
  !! beyond the cutoff.
  !! @param[in] k wave number
  !! @param[in] iu solute atom number
  !! @param[in] iv solvent species number
  !! @return The value
  function rism3D_potential_ulj_gt_cut_k(this,k,iu,iv) result(ulj)
    use constants, only : pi
    implicit none
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(in) :: k
    integer, intent(in) :: iu, iv
    _REAL_ :: ulj
    _REAL_ :: cut, cut_inv
    _REAL_, external :: dsi
    cut = sqrt(this%ljCutoffs2(iu,iv))
    cut_inv = 1d0/cut
    ! could be improved with Horner's method

    ! write(6,*) 'ljCutoff',iu,iv,sqrt(this%ljCutoffs2(iu,iv))
    ! write(6,*) 'ljA ljB',iu,iv, this%ljAUV(iu,iv),this%ljBUV(iu,iv)
    ! write(6,*) 'dsi',cut, k,dsi(cut*k)
    ! r6
    ulj = -this%ljBUV(iu,iv) * pi/ (12d0 * k) &
         *( k**4 * (pi - 2*dsi(cut*k)) &
         
         -2*(cut_inv**2*k**2 - 6d0*cut_inv**4)*sin(cut*k) &
         
         -2*k*(cut_inv*k**2 - 2d0*cut_inv**3)*cos(cut*k))

    ulj = ulj + this%ljAUV(iu,iv) * pi/ (1814400d0 * k) &
         *(-k**10*(pi - 2*dsi(cut*k)) &
         
         +2*((((    cut_inv**2*k**2 &
              - 6d0*cut_inv**4)*k**2 &
            + 120d0*cut_inv**6)*k**2 &
           - 5040d0*cut_inv**8)*k**2 &
         + 362880d0*cut_inv**10)&
         *sin(cut*k) &
         
         +2*k*(((( cut_inv*k**2&
             - 2d0*cut_inv**3)*k**2 &
            + 24d0*cut_inv**5)*k**2 &
           - 720d0*cut_inv**7)*k**2 &
         + 40320d0*cut_inv**9)&
         *cos(cut*k))

  end function rism3D_potential_ulj_gt_cut_k
  
  !> Returns the values of the long range asymptotics of the total
  !! correlation function h(r) divided by the species charge at the
  !! given point.  A cut off is used to accelerate the calculation,
  !! which should be precomputed by _setcut().
  !! See equation 32 in Kovalenko & Hirata 2000.
  !! @param[in] this Potential object.
  !! @param[in] r Position.
  !! @return Long range asymptotics.
  function tcf_long_range_asymptotics_R(this, r) result(tcfLRasymp)
    use constants, only : pi
    implicit none
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(in) :: r(3)
    _REAL_ :: tcfLRasymp
    _REAL_ :: distance
    ! 1/2 * kappa * smear.
    _REAL_ :: half_xappa_smear
    integer :: iu

    half_xappa_smear = 0.5d0 * this%solvent%xappa * this%chargeSmear
    tcfLRasymp = 0
    if (.not. this%solvent%ionic) return
    do iu = 1, this%solute%numAtoms
       distance = sqrt(sum((r(:) - this%solute%position(:, iu))**2))
       if (distance > 0) then
          tcfLRasymp = tcfLRasymp &
               - this%solute%charge(iu) / distance &
               * (exp(-this%solvent%xappa * distance) &
               * erfc(half_xappa_smear - distance / this%chargeSmear) &
               - exp(this%solvent%xappa * distance) &
               * erfc(half_xappa_smear + distance / this%chargeSmear)) / 2d0
       else ! distance == 0d0
          tcfLRasymp = tcfLRasymp &
               - this%solute%charge(iu) &
               * (2d0 / (sqrt(PI) * this%chargeSmear) &
               - exp(half_xappa_smear**2) &
               * this%solvent%xappa * erfc(half_xappa_smear)) &
               / exp(half_xappa_smear**2)
          ! Outside the loop we multiply by the inverse of the last line.
          ! Since the statement is rarely executed, this should save
          ! a bit of time.
       end if
    end do
    tcfLRasymp = tcfLRasymp * exp(half_xappa_smear**2) / this%solvent%dielconst
  end function tcf_long_range_asymptotics_R


  !> Returns the long range asymptotic values of the direct
  !! correlation function c(r) divided by the species charge at the
  !! given point.  A cut off is used to accelerate the calculation,
  !! which should be precomputed by _setcut().
  !! See equation 29 in Kovalenko & Hirata 2000.
  !! @param[in] this Potential object.
  !! @param[in] r Position.
  !! @return The long range portion of the Ewald sum of the direct
  !! correlation function divided by the species charge at the given
  !! point.
  function dcf_long_range_asymptotics_R(this, r) result(dcfLRasymp)
    use constants, only : pi
    implicit none
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(in) :: r(3)
    _REAL_ :: dcfLRasymp
    _REAL_ :: distance
    integer :: iu

    !JAJ This returns the long range portion of the electrostatic
    ! periodic field by calculating the difference between the
    ! Coulombic potential and the short range Ewald term.
    dcfLRasymp = 0
    if (.not. this%solute%charged) return
    do iu = 1, this%solute%numAtoms
       distance = sqrt(sum((r(:) - this%solute%position(:, iu))**2))
       if (distance > 0) then
          dcfLRasymp = dcfLRasymp &
               - this%solute%charge(iu) / distance * erf(distance / this%chargeSmear)
       else ! distance == 0d0
          dcfLRasymp = dcfLRasymp &
               - this%solute%charge(iu) / (sqrt(PI) * this%chargeSmear) * 2d0
       end if
    end do
  end function dcf_long_range_asymptotics_R

  !> Applying minimum-image convention to find distance from grid
  !! point to nearest solute atom image, which may be in an adjacent
  !! cell. Closest image is based on solute unit cell dimensions.
  !! The closest image is whichever one is less than half a
  !! unit cell away, i.e. minimum image convention.
  !! This is only useful for periodic solute and requires a defined
  !! unit cell.
  !! Note that the minimum image convention implicitly implies a
  !! close-range interaction cutoff at half the unit cell length,
  !! which may or may not be physically realistic depending on the
  !! system.
  !! DAC note: might be worth having a special routine for
  !!   orthogonal boxes
  function minimumImage(this, position)
   implicit none
    _REAL_ :: minimumImage(3)
    type(rism3d_potential), intent(in) :: this !< potential object.
    _REAL_, intent(in) :: position(3) !< Position vector.
    !! The vector origin is usually the grid point a calculation is
    !! performed at.

    integer :: id

    ! Applying minimal-image convention to find distance from grid
    ! point to nearest solute atom image. Closest image is based on
    ! solute unit cell dimensions.

    _REAL_ :: f(3)
    do id = 1, 3
       ! 1. Transform to fractional coordinates (i.e., reciprocal space).
       f(id) = dot_product(position, this%grid%unitCellVectorsK(id, :))
       ! 2. Round to the nearest whole unit and subtract from
       ! coordinate to obtain the minimum image.
       f(id) = f(id) - anint(f(id))
    end do
    ! f = f - anint(f)
    do id = 1, 3
       ! 3. Transform back to Cartesian coordinates using the
       ! transpose of the matrix of Cartesian unit cell vectors.
       minimumImage(id) = dot_product(f, this%grid%unitCellVectorsR(:, id))
    end do
  end function minimumImage

  subroutine getnojellywt(this)
    implicit none
    type(rism3d_potential), intent(inout) :: this !< potential object.
    _REAL_ :: q1
    integer :: iv

    this%phineut => safemem_realloc(this%phineut,this%solvent%numAtomTypes, .false.)
    this%phineut(:) =  0.0

    if (this%solvent%ionic .and. this%periodic .and.  &
        (this%periodicPotential == 'pme' .or. this%periodicPotential == 'ewald') ) then
        q1 = 0
        do iv = 1,this%solvent%numAtomTypes
          if (this%solvent%atomName(iv) /= "O" .and. &
              this%solvent%atomName(iv) /= "H1" .and. &
              this%solvent%charge(iv) /= 0 ) then
                q1 = q1 + this%solvent%density(iv)*this%solvent%charge(iv)*this%solvent%charge(iv)
          end if
        end do

        q1 =  (this%solute%totalCharge)/(this%grid%boxVolume * q1 )

        do iv = 1,this%solvent%numAtomTypes
           if (this%solvent%atomName(iv) /= "O" .and. &
               this%solvent%atomName(iv) /= "H1") then
               this%phineut(iv)=this%solvent%charge(iv) * q1
           end if

        end do
    end if
    
  end subroutine getnojellywt

end module rism3d_potential_c

