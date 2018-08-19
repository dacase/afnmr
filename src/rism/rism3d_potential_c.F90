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
     !> Cutoff**2 for RISM potential calculations.
     _REAL_ :: cutoff2
     !> Cutoff for long range asymptotics of the total correlation function.
     _REAL_ :: cut_hlr
     !> Cutoff for long range asymptotics of the direct correlation function.
     _REAL_ :: cut_clr

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
  end type rism3d_potential

  private mixSoluteSolventLJParameters
  private uvCoulombicPotential, uvEwaldSumPotential !, uvParticleMeshEwaldPotential
  private uvParticleMeshEwaldPotential
  private uvLennardJonesPotentialWithCutoff, uvLennardJonesPotentialWithMinimumImage

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
  subroutine rism3d_potential_new(this, grid, solv, solu, cut, fft, periodicPotential, biasPotential)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    type(rism3d_grid), target, intent(in) :: grid
    type(rism3d_solute), target, intent(in) :: solu
    type(rism3d_solvent), target, intent(in) :: solv
    _REAL_, intent(in):: cut
    type(rism3d_fft), target, intent(in) :: fft
    character(len=*), intent(in) :: periodicPotential
    _REAL_, intent(in) :: biasPotential
    this%grid => grid
    this%solvent => solv
    this%solute => solu
    this%fft => fft
    this%periodicPotential = periodicPotential
    if (this%periodicPotential /= '') then
       this%periodic = .true.
    end if
    this%biasPotential = biasPotential
    call rism3d_potential_setCut(this, cut)
    this%ljSigmaUV => safemem_realloc(this%ljSigmaUV, this%solute%numAtoms, this%solvent%numAtomTypes, .false.)
    this%ljEpsilonUV => safemem_realloc(this%ljEpsilonUV, this%solute%numAtoms, this%solvent%numAtomTypes, .false.)
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


  !> Sets the cut off distance for potential, force and long range
  !! asymptotics calculations.  The long range asymptotics cut
  !! off is calculated from the solvent parameters.
  !! @param[in,out] this potential object.
  !! @param[in] cut Distance cutoff for potential and force calculations.
  subroutine rism3d_potential_setcut(this, cut)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(in) :: cut

    ! exp_quarter_xappa2_smear2 : exp(-1 / 4 * (kappa^2) * (smear^2))
    ! half_xappa_smear      : 1 / 2* kappa   * smear
    ! tcfLongRangeAsympR_coeff_at_cut   : long range h without charge or 1/r dependence (tcfLongRangeAsympR/(Qq) * r)
    ! dcfLongRangeAsympR_coeff_at_cut   : long range c without charge or 1/r dependence (dcfLongRangeAsympR/(Qq) * r)
    ! tcut                  : trial cutoff
    _REAL_ :: exp_quarter_xappa2_smear2, tcfLongRangeAsympR_coeff_at_cut, &
         dcfLongRangeAsympR_coeff_at_cut, half_xappa_smear, tcut


    ! Assign potential cutoff.
    ! Ensure that the square won't overflow.
    this%cutoff = min(sqrt(huge(1d0)), cut)
    this%cutoff2 = this%cutoff**2

    half_xappa_smear = 0.5d0 * this%solvent%xappa * this%solvent%smear
    exp_quarter_xappa2_smear2 = exp(half_xappa_smear**2)


    ! This is really a very simple minded approach to finding the cut off distance.
    ! Basically, it is a linear search for a value where the non-1/r part is less
    ! than an error value.  This really should be done with Newton-Raphson instead.
    tcut = 4d0
    this%cut_hlr = -1d0
    this%cut_clr = -1d0
    do while((this%cut_hlr <0 .or. this%cut_clr < 0) .and. tcut < 1000d0)
       tcfLongRangeAsympR_coeff_at_cut = exp_quarter_xappa2_smear2 &
            * (exp(-this%solvent%xappa * tcut) * erfc(half_xappa_smear - tcut / this%solvent%smear) &
            - exp(this%solvent%xappa * tcut) * erfc(half_xappa_smear + tcut / this%solvent%smear)) / 2d0
       dcfLongRangeAsympR_coeff_at_cut = (1d0 - erfc(tcut / this%solvent%smear))
       if (exp_quarter_xappa2_smear2 - tcfLongRangeAsympR_coeff_at_cut < 1d-7 .and. this%cut_hlr < 0) &
            this%cut_hlr = tcut
       if (1d0 - dcfLongRangeAsympR_coeff_at_cut < 1d-7 .and. this%cut_clr < 0) &
            this%cut_clr = tcut
       tcut = tcut + .1d0
    end do

    !! $      if (this%cut_hlr < 0) this%cut_hlr = HUGE(1d0)
    !! $      if (this%cut_clr < 0) this%cut_clr = HUGE(1d0)
    this%cut_hlr = HUGE(1d0)
    this%cut_clr = HUGE(1d0)
  end subroutine rism3d_potential_setcut


  !> Calculates the potential on the grid.
  subroutine rism3d_potential_calc(this)
    use rism_util, only : checksum
    use rism3d_opendx, only : rism3d_opendx_write
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif
    type(rism3d_potential), intent(inout) :: this !< potential object.

    integer :: id

    integer :: ix, iy, iz, ierr
    character(len=30) :: filename

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
          if (this%periodicPotential == 'pme') then
             call uvParticleMeshEwaldPotential(this, this%uuv)
          else
             call uvEwaldSumPotential(this, this%uuv)
          end if
       else
          call uvCoulombicPotential(this, this%uuv)
       end if
    end if
    call rism_timer_stop(this%coulombTimer)
    call timer_stop(TIME_UCOULU)

    call timer_start(TIME_ULJUV)
    call rism_timer_start(this%ljTimer)
    if (this%periodic) then
       call uvLennardJonesPotentialWithMinimumImage(this, this%uuv)
    else
       call uvLennardJonesPotentialWithCutoff(this, this%uuv)
    end if
    call rism_timer_stop(this%ljTimer)
    call timer_stop(TIME_ULJUV)

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
  !! should be transparent to the user.
  subroutine  rism3d_potential_dcf_tcf_long_range_asymptotics(this)
    use constants, only : PI, FOURPI
    use rism_util, only : checksum
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d_potential), intent(inout) :: this !< potential object.

    _REAL_ :: uc1gc, uc1gh
    integer :: igx, igy, igz, ig, iu, iv, ig0
    _REAL_ :: phase, sumcos, sumsin, rx, ry, rz, delx, dely, delz, delyz2, ra, &
         uc1g, sr, uc, ucs, ucsc
    _REAL_ :: xappa2, smear2_4
    ! z-axis offset (important for spatially distributed MPI).
    _REAL_ :: offset

    _REAL_ :: solutePosition(3), totalSoluteCharge

    !      if (.not. this%solute%charged) return
    
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

    ! Set Gaussian charge smearing coefficients.
    smear2_4 = this%solvent%smear**2 / 4d0
    xappa2 = this%solvent%xappa**2

    ! Initialize correlation function arrays.
    this%dcfLongRangeAsympK = 0d0
    if (this%solvent%ionic) then
       this%tcfLongRangeAsympK = 0d0
    end if

    if (.not. this%periodic) then
       !JAJ This is the real-space long range portion of the
       ! solute-solvent Ewald potential!

       ! Asymptotic values for the DCF and TCF in R-space.
       ! These are eq. 29 and 32 of Kovalenko/Hirata 2000.
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
                   if (this%solvent%ionic) then
                   this%tcfLongRangeAsympR(ig) = tcf_long_range_asymptotics_R(this, (/rx, ry, rz/))
                end if
             end do
          end do
       end do
    end if


    if ((.not. this%periodic) .or. (this%periodicPotential == 'ewald')) then

       ! Getting long-range part of the DCF and TCF in k-space.
       ! k = 0 is on master.  Set it to zero as it is handled separately.

       !JAJ This is the long range portion of the solute electrostatic 
       ! periodic field!
       ! See Understanding Molecular Simulation 2E, eq. 12.1.15 on p. 297.
       ! Note: the result is divided by volume since FFTW does not do this
       ! during transforms, though technically this is not part of the 
       ! frequency space potential.

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
          sumcos = 0d0
          sumsin = 0d0
          do iu = 1, this%solute%numAtoms
             solutePosition = this%solute%position(:, iu)
             if (this%periodic) then
                ! Minimum image convention.
                solutePosition = minimumImage(this, solutePosition)
             end if
             phase = dot_product(this%grid%waveVectors(:, ig), solutePosition)
             sumcos = sumcos + this%solute%charge(iu) * cos(phase)
             sumsin = sumsin + this%solute%charge(iu) * sin(phase)
          end do

          ! DCF in k-space.
          uc1g = FOURPI * exp(-smear2_4 * this%grid%waveVectors2(ig))
          uc1gc = uc1g / this%grid%waveVectors2(ig)
          this%dcfLongRangeAsympK(2 * ig - 1) = uc1gc * sumcos / this%grid%boxVolume
          this%dcfLongRangeAsympK(2 * ig) = uc1gc * sumsin / this%grid%boxVolume

          ! TCF in k-space.
          if (this%solvent%ionic) then
             uc1gh = uc1g / (this%grid%waveVectors2(ig) + xappa2)
             this%tcfLongRangeAsympK(2 * ig - 1) = uc1gh * sumcos / this%grid%boxVolume
             this%tcfLongRangeAsympK(2 * ig) = uc1gh * sumsin / this%grid%boxVolume
          end if
       end do

    end if

    if (.not. this%periodic) then
       ! Getting the difference between long-range asymptotic functions
       ! of the TCF at k = 0.

       !TODO: Understanding Molecular Simulation suggests setting the k = 0
       ! case to 0.  In that case, where does this analytic treatment arise?
       if (this%grid%offsetK(3) == 0) then
          sumcos = 0d0
          sumsin = 0d0
          do iu = 1, this%solute%numAtoms
             solutePosition = this%solute%position(:, iu)
             if (this%periodic) then
                ! Minimum image convention.
                !FIXME: This should be relative to the current grid point position!
                !FIXME: This should never be reached anyways, so it
                ! probably should be removed.
                solutePosition = minimumImage(this, solutePosition)
             end if
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
    end if
    call rism_timer_stop(this%asympTimer)
  end subroutine rism3d_potential_dcf_tcf_long_range_asymptotics


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

       x2arg = (k * this%solvent%xappa * this%solvent%smear)**2
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
  !! @param[in] this The closure object.
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
  subroutine uvLennardJonesPotentialWithCutoff(this, ulj)
    use rism_util, only : checksum
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(inout) :: ulj(:,:,:,:)
    ! Amount to offset the grid in the z-axis (when using MPI).
    _REAL_ :: offset
    ! Number of gridpoints within cutoff for each dimension.
    integer :: cutoff(3)
    ! Closest grid point to solute atom center.
    integer :: closestPoint(3)
    ! First and last gridpoints which are within the cutoff distance
    ! from the solute atom while remaining in the box.
    integer :: first(3), last(3)
    ! Grid, solute, slovent, and dimension indices.
    integer ::  igx, igy, igz, ig, iu, iv, id
    ! Grid point position.
    _REAL_ ::  rx, ry, rz
    ! Distance of grid point from solute.
    _REAL_ ::  dx2, dy2, dz2, dy2dz2, r2
    ! LJ site term and its sixth power inverse.
    _REAL_ ::  rs2, rs6i
    ! Minimum LJ term r / sigma and its square.
    _REAL_ ::  rcor, rcor2
    parameter (rcor = 0.2d0, rcor2 = rcor**2)

    ! Calculate the number of grid points necessary to cover this
    ! cutoff radius.  In case of a large cutoff, we need to protect
    ! against integer overflow.
    cutoff = nint(min(dble(huge(1) - 1), sqrt(this%cutoff2) / this%grid%spacing)) + 1
    ! Offset used for dividing the grid into slabs along the z-axis
    ! for MPI.
    offset = this%grid%spacing(3) * this%grid%offsetR(3)

    do iu = 1, this%solute%numAtoms
       ! Find closest grid point near solute atom center.
       do id = 1, 3
          closestPoint(id) = nint(this%solute%position(id, iu) / this%grid%spacing(id))
       end do
       ! Apply z-axis offset for MPI (if applicable).
       closestPoint(3) = closestPoint(3) - this%grid%offsetR(3)

       ! Find the first and last grid points which are within the
       ! cubic cutoff distance while remaining within the box.
       do id = 1, 3
          !TODO:JAJ Should we throw exception in event of cutoff going
          ! beyond the box edges?
          if (closestPoint(id) > 0) then
             first(id) = min(max(1, closestPoint(id) - cutoff(id)), &
               this%grid%localDimsR(id))
             last(id) = max(min(this%grid%localDimsR(id), &
               closestPoint(id) + min(huge(1) - closestPoint(id),cutoff(id))),1)
          else
             first(id) = min(max(1, closestPoint(id) - &
                  min(huge(1) + closestPoint(id), cutoff(id))), &
                  this%grid%localDimsR(id))
             last(id) = max(min(this%grid%localDimsR(id), &
                  closestPoint(id) + cutoff(id)), 1)
          end if
       end do

       ! Calculate the solute-solvent Lennard-Jones potential at each
       ! point within the cutoff radius.
       do igz = first(3), last(3)
          rz = (igz - 1) * this%grid%spacing(3) + offset
          dz2 = (rz - this%solute%position(3, iu))**2
          do igy = first(2), last(2)
             ry = (igy - 1) * this%grid%spacing(2)
             dy2 = (ry - this%solute%position(2, iu))**2
             dy2dz2 = dz2 + dy2
             do igx = first(1), last(1)
                rx = (igx - 1) * this%grid%spacing(1)

                ! Find distance from grid point to solute site in real
                ! space.
                dx2 = (rx - this%solute%position(1, iu))**2
                r2 = dx2 + dy2dz2

                ! Ensure grid point is within radial cutoff distance
                ! from solute atom. This is necessary in order to
                ! enforce a radial cutoff since 'first' and 'last'
                ! grid points defined above treated cutoff as cubic.
                if (r2 < this%cutoff2) then
                   do iv = 1, this%solvent%numAtomTypes
                      rs2 = r2 / this%ljSigmaUV(iu, iv)**2
                      ! Enforce minimum LJ term value, rcor.
                      if (rs2 < rcor2)  rs2 = rcor2
                      rs6i = 1d0 / rs2**3
                      ulj(igx, igy, igz, iv) = ulj(igx, igy, igz, iv) &
                           + this%ljEpsilonUV(iu, iv) * rs6i * (rs6i - 2.d0)
                   end do
                end if
             end do
          end do
       end do
    end do
  end subroutine uvLennardJonesPotentialWithCutoff


  !> Tabulate the solute-solvent 12-6 Lennard-Jones potential in the
  !! box subject to the minimum image convention. No cutoff is applied.
  !! @param[in] this potential object
  !! @param[in,out] ulj grid to add potential to
  subroutine uvLennardJonesPotentialWithMinimumImage(this, ulj)
    use constants, only : PI
    use rism_util, only : checksum
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(inout) :: ulj(:,:,:,:)
    ! Amount to offset the grid in the z-axis (when using MPI).
    _REAL_ :: offset
    ! Grid, solute, slovent, and dimension indices.
    integer :: igx, igy, igz, ig, iu, iv, id
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
          ! rz = (igz - 1) * this%grid%voxelVectorsR(3, :) + offset
          do igy = 1, this%grid%localDimsR(2)
             ry = (igy - 1) * this%grid%voxelVectorsR(2, :)
             do igx = 1, this%grid%localDimsR(1)
                rx = (igx - 1) * this%grid%voxelVectorsR(1, :)

                gridPoint = rx + ry + rz

                ! Vector from grid point to solute atom, in real space.
                solutePosition = gridPoint - this%solute%position(:, iu)

                ! Adjust distance to be to the closest solute atom
                ! image, which may be in an adjacent cell. (The
                ! closest image is whichever one is less than half a
                ! unit cell away, i.e. minimum image convention.)
                solutePosition = minimumImage(this, solutePosition)

                soluteDistanceSquared = dot_product(solutePosition, solutePosition)

                if (soluteDistanceSquared < this%cutoff2) then
                   do iv = 1, this%solvent%numAtomTypes
                      ljBaseTerm = soluteDistanceSquared / this%ljSigmaUV(iu, iv)**2
                      ! Enforce minimum LJ term to prevent division by zero below.
                      if (ljBaseTerm < minLJBaseTerm) ljBaseTerm = minLJBaseTerm
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

    ! Constants
    _REAL_ :: rcor, rcor2
    parameter (rcor = 0.002, rcor2 = rcor**2)

    ! Potential on sparse grid.
    _REAL_, pointer :: smucu(:,:,:) => NULL()
    ! Number of grid points on sparse grid.
    integer :: smng3(3)
    ! Ratio of coarse to fine grid dimensions.
    _REAL_ :: smratio(3)
    ! _REAL_ :: smratio
    ! For MPI spatial decomposition.
    ! offset :: distance in z-axis due to spatial decomp
    _REAL_ :: offset, smoffset
    ! smstart :: initial grid point where the sparse and fine grids have the same value
    integer :: smstart
    ! Linear spacing of the grid.
    _REAL_ :: smgrdspc(3)
    ! Number of gridpoints in each direction to use.
    ! Closest point to solute atom.
    integer :: cutoff(3), closestPoint(3), id, iu2
    ! Flag indicating that this index of the array has not been set
    _REAL_ :: untouched = HUGE(0d0)
    ! ???
    integer :: gxstep, gxstart, first(3), last(3)

    ! Temporary / intermediate variables.
    integer :: igx, igy, igz, ig, iu, iv
    _REAL_ :: r2, ra, rx, ry, rz, dx, dy, dz, dz2, dyz2
    integer :: iigx, iigy, iigz
    _REAL_ :: tmp(1)
    ! Center of (local?) r-space grid.
    integer :: rCenter(3)

    ! Screen output.
#ifdef RISM_DEBUG
    write(6,*) 'Tabulating fast solute Coulomb potential ...'
    call flush(6)
#endif /*RISM_DEBUG*/
    rCenter = ((this%grid%localDimsR + 1) / 2) + 1
    smucu => safemem_realloc(smucu, rCenter(1), rCenter(2), rCenter(3), .false.)

    ! Clearing summators.
    ucu(:,:,:,1) = untouched
    ! Calculate exact potential over the entire sparse grid.
    smucu = 0
    ! ng3 + 1 takes care of values that are odd
    smng3 = (this%grid%localDimsR + 1) / 2 + 1
    smgrdspc = this%grid%spacing * 2d0
    smratio = this%grid%spacing / smgrdspc
    ! smratio = (/ 2d0, 2d0, 2d0 /)
    ! smgrdspc = this%grid%spacing * smratio
    ! Calculate the number of grid point necessary to cover this
    ! cutoff range.  In case of a large cutoff, we need to protect
    ! against integer overflow.
    cutoff = nint(min(dble(huge(1) - 1), sqrt(this%cutoff2) / this%grid%spacing)) + 1
    offset = this%grid%spacing(3) * (this%grid%offsetR(3))
    smstart = mod(this%grid%offsetR(3), 2) + 1
    ! If we are start on an odd index (nz_start starts from 0) shift
    ! the small grid to a lower offset.
    smoffset = offset - (smstart - 1) * this%grid%spacing(3)

    ! Enumerating box grid points.
    do iu = 1, this%solute%numAtoms
       do igz = 1, smng3(3)
          rz = (igz - 1) * smgrdspc(3) + smoffset
          dz = rz - this%solute%position(3, iu)
          dz2 = dz * dz
          do igy = 1, smng3(2)
             ry = (igy - 1) * smgrdspc(2)
             dy = ry - this%solute%position(2, iu)
             dyz2 = dy * dy + dz2
             do igx = 1, smng3(1)
                !  ig = 1 + igx + igy * smng3(1) + igz * smng3(2) * smng3(1)
                rx = (igx - 1) * smgrdspc(1)
                ! Summing over solute sites.
                ! Getting site separation.
                dx = rx - this%solute%position(1, iu)

                ra = max(sqrt(dx * dx + dyz2), rcor)
                smucu(igx, igy, igz) = smucu(igx, igy, igz) + this%solute%charge(iu) / ra
             end do
          end do
       end do
    end do

    ! Transfer data from the sparse grid to the dense grid.
    do igz = smstart, this%grid%localDimsR(3), 2
       do igy = 1, this%grid%localDimsR(2), 2
          do igx = 1, this%grid%localDimsR(1) ,2
             ucu(igx, igy, igz, 1) &
                  = smucu((igx + 1) / 2, (igy + 1) / 2, smstart -1 + (igz + 1) / 2)
          end do
       end do
    end do

    ! Calculate exact potential for each grid point within the cutoff.
    do iu = 1, this%solute%numAtoms
       ! Determine grid points within both the cutoff and the closest solute site.
       do id = 1, 3
          ! Central point (grid point nearest atom).
          closestPoint(id) = nint(this%solute%position(id, iu) / this%grid%spacing(id))
       end do
       closestPoint(3) = closestPoint(3) - this%grid%offsetR(3)
       do id = 1, 3
          first(id) = min(max(1, closestPoint(id) - cutoff(id)), this%grid%localDimsR(id)) ! Smallest grid index within cutoff.
          ! Note: we have to protect against closestPoint(id) + cutoff(id) overflowing.
          ! Largest grid index within cutoff.
          last(id) = max(min(this%grid%localDimsR(id), &
               closestPoint(id) + min(huge(1) - closestPoint(id), cutoff(id))), 1)
       end do
       do igz = first(3), last(3)
          rz = (igz - 1) * this%grid%spacing(3) + offset
          dz = rz - this%solute%position(3, iu)
          dz2 = dz * dz
          do igy = first(2), last(2)
             ry = (igy - 1) * this%grid%spacing(2)
             dy = ry - this%solute%position(2, iu)
             dyz2 = dy * dy + dz2
             do igx = first(1), last(1)
                if (ucu(igx, igy, igz, 1) == untouched) then
                   rx = (igx - 1) * this%grid%spacing(1)

                   ! Site-site interaction subject to radial cutoff.
                   dx = rx - this%solute%position(1, iu)
                   r2 = dx * dx + dyz2
                   if (r2 < this%cutoff2) then
                      ucu(igx, igy, igz, 1) = 0d0
                      do iu2 = 1, this%solute%numAtoms
                         dx = rx - this%solute%position(1, iu2)
                         dy = ry - this%solute%position(2, iu2)
                         dz = rz - this%solute%position(3, iu2)

                         ra = max(sqrt(dx * dx + dy * dy + dz * dz), rcor)
                         ucu(igx, igy, igz, 1) = ucu(igx, igy, igz, 1) + &
                            this%solute%charge(iu2) / ra
                      end do
                   end if
                end if
             end do
          end do
       end do
    end do

    ! Interpolate missing values for remaining grid points.
    do igz = 1, this%grid%localDimsR(3)
       rz = (igz - 1) * this%grid%spacing(3) + offset
       first(3) = (igz + 1) / 2
       last(3) = first(3) + 1
       do igy = 1, this%grid%localDimsR(2)
          ry = (igy - 1) * this%grid%spacing(2)
          first(2) = (igy + 1) / 2
          last(2) = first(2) + 1
          if (mod(igy, 2) /= 0 .and. mod(igz + smstart - 1, 2) /= 0) then
             gxstep = 2
             gxstart = 2
          else
             gxstep = 1
             gxstart = 1
          endif
          do igx = gxstart, this%grid%localDimsR(1), gxstep
             rx = (igx - 1) * this%grid%spacing(1)
             first(1) = (igx + 1) / 2
             last(1) = first(1) + 1
             if (ucu(igx, igy, igz, 1) /= untouched) then
                cycle
             endif
             call blend_103(&
                  (dble(igx - 1) * smratio(1) - first(1) + 1) / (last(1) - first(1)),&
                  (dble(igy - 1) * smratio(2) - first(2) + 1) / (last(2) - first(2)),&
                  (dble(igz - 1 + smstart - 1) * smratio(3) - first(3) + 1) / (last(3) - first(3)),&
                  smucu(first(1), first(2), first(3)),&
                  smucu(first(1),  first(2), last(3)),&
                  smucu(first(1), last(2), first(3)),&
                  smucu(first(1), last(2), last(3)),&
                  smucu(last(1), first(2), first(3)),&
                  smucu(last(1), first(2), last(3)),&
                  smucu(last(1), last(2), first(3)),&
                  smucu(last(1), last(2), last(3)),&
                  ucu(igx, igy, igz, 1))
          end do
       end do
    end do

    do iv = this%solvent%numAtomTypes, 1, -1
       ucu(:,:,:,iv) =  ucu(:,:,:,1) * this%solvent%charge(iv)
    end do
    ucu = ucu

    if (safemem_dealloc(smucu) /= 0) &
         call rism_report_error("Failed to deallocate memory in uvCoulombicPotential")
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
#if defined(MPI)
    this%uuv1d(2:this%grid%totalLocalPointsK:2, :) = &
         -this%uuv1d(2:this%grid%totalLocalPointsK:2, :)
    call rism3d_fft_bwd(this%fft, this%uuv1d)
#else
    call rism3d_fft_bwd(this%fft, this%uuv1d)
#endif /*defined(MPI)*/

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

    ! Calculate electrostatic potential energy on each grid point.
    ! No 1/2 term is required since solute affects solvent, but
    ! solvent does not affect solute, hence no double counting.
    do iv = this%solvent%numAtomTypes, 1, -1
       if (index(this%solvent%atomName(iv), "-") == 0 .and. &
            index(this%solvent%atomName(iv), "+") == 0) then
          ucu(:,:,:,iv) = ucu(:,:,:,1) * this%solvent%charge(iv)
       else
          ucu(:,:,:,iv) = (ucu(:,:,:,1) + this%biasPotential) * this%solvent%charge(iv)
       end if
    end do
    ! TODO: Why is this done?
    ucu = ucu
  end subroutine uvEwaldSumPotential


  !> Short-range portion of the Ewald sum electric potential.
  !! TODO: Document me.
  subroutine uvEwaldSumShortRangePotential(this, ucu)
    use constants, only: pi
    use rism3d_opendx, only : rism3d_opendx_write
    implicit none
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(inout) :: ucu(:,:,:,:)

    ! Minimum distance to prevent division by zero.
    _REAL_, parameter :: minDistance = 0.002

    integer :: igx, igy, igz
    integer :: ig
    integer :: iu

    _REAL_ :: rx(3), ry(3), rz(3)
    _REAL_ :: gridPoint(3), solutePosition(3)
    _REAL_ :: soluteDistance
    _REAL_ :: sr
    ! Calculate short-range term of Ewald sum and combine with
    ! previously calculated long-range term.
    do igz = 0, this%grid%localDimsR(3) - 1
       rz = (igz + this%grid%offsetR(3)) * this%grid%voxelVectorsR(3, :)
       do igy = 0, this%grid%localDimsR(2) - 1
          ry = igy * this%grid%voxelVectorsR(2, :)
          do igx = 0, this%grid%localDimsR(1) - 1
             ig = 1 + igx + igy * this%grid%localDimsR(1) + &
                  igz * this%grid%localDimsR(2) * this%grid%localDimsR(1)
             rx = igx * this%grid%voxelVectorsR(1, :)

             gridPoint = rx + ry + rz

             do iu = 1, this%solute%numAtoms
                ! Solute atom position relative to a grid point.
                solutePosition = gridPoint - this%solute%position(:, iu)

                ! Apply minimum image convention.
                solutePosition = minimumImage(this, solutePosition)

                ! Distance from solute atom to gridpoint.
                soluteDistance = sqrt(dot_product(solutePosition, solutePosition))

                if (soluteDistance < minDistance) soluteDistance = minDistance

                ! Short-range term of Ewald sum.
                if (soluteDistance < this%cutoff) then
                   !TODO: Should this compare to some small value?
                   ! if (soluteDistance > 0) then
                      ! Second term in Eq. 16 of Kovalenko/Hirata 2000
                      ! (missing charge in paper is a typo). See
                      ! Understanding Molecular Simulations 2E chapter 12
                      ! for derivation (note the differing smear
                      ! coefficient definition).
                      sr = this%solute%charge(iu) &
                           * erfc(soluteDistance / this%solvent%smear) / soluteDistance
                   ! else
                   !    ! Limiting case approaching zero.
                   !    !TODO: Verify this limit using equation above.
                   !    sr = -2 * this%solute%charge(iu) / (sqrt(PI) * this%solvent%smear)
                   ! end if
                   ucu(igx + 1, igy + 1, igz + 1, 1) = ucu(igx + 1, igy + 1, igz + 1, 1) + sr
                end if
             end do
          end do
       end do
    end do
  end subroutine uvEwaldSumShortRangePotential

  !> Long-range portion of the Particle Mesh Ewald (PME) electric potential.
  !! TODO: Document me!
  subroutine uvParticleMeshEwaldPotential(this, ucu)
    use, intrinsic :: iso_c_binding
    use bspline
    use constants, only : pi
    use FFTW3
    use rism3d_opendx, only : rism3d_opendx_write
    use rism_util, only: r2c_pointer
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(inout) :: ucu(:,:,:,:)

    !TODO: Associate smear with that used for dilute 3D-RISM.
    ! Ewald charge smear parameter.
    _REAL_, parameter :: smear = 1
    ! Order of b-spline interpolation.
    integer, parameter :: splineOrder = 4

    _REAL_, parameter :: zeta = smear * smear
    _REAL_, parameter :: t = -0.25 / zeta

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

    L = this%grid%localDimsR(1)
    M = this%grid%localDimsR(2)
    N = this%grid%globalDimsR(3)

    call timer_start(TIME_UCOULULR)
    call rism_timer_start(this%coulombLongRangeTimer)

    gridDimX_k = this%grid%localDimsR(1) / 2 + 1
    kxi => safemem_realloc(kxi,  gridDimX_k, 3, .false.)
    kyi => safemem_realloc(kyi,  this%grid%localDimsR(2), 3, .false.)
    kzi => safemem_realloc(kzi,  this%grid%localDimsR(3), 3, .false.)

    bsplineFourierCoeffX => safemem_realloc(bsplineFourierCoeffX, this%grid%globalDimsR(1), .false.)
    bsplineFourierCoeffY => safemem_realloc(bsplineFourierCoeffY, this%grid%globalDimsR(2), .false.)
    bsplineFourierCoeffZ => safemem_realloc(bsplineFourierCoeffZ, this%grid%globalDimsR(3), .false.)

    gaussianFourierCoeff => safemem_realloc(gaussianFourierCoeff, &
         gridDimX_k * this%grid%localDimsR(2) * this%grid%localDimsR(3), .false.)

#if defined(MPI)
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
#if defined(MPI)
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

    if (any(mod(this%grid%localDimsR, 2) > 0)) then
       !TODO: call RISM error functions
       write(6,'(a,a)')  "| PME implementation prefers an even", &
            " number of grid points on each axis."
    end if

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

    ! Compute the discrete Fourier transform coefficients of a
    ! Gaussian monopole.
    ixyz = 1
    do iz = 0, this%grid%localDimsR(3) - 1
       do iy = 0, this%grid%localDimsR(2) - 1
          do ix = 0, gridDimX_k - 1
             waveVector = kxi(ix+1,:) + kyi(iy+1,:) + kzi(iz+1,:)
             k2 = dot_product(waveVector, waveVector)
             gaussianFourierCoeff(ixyz) = exp(k2 * t)
             ixyz = ixyz + 1
          end do
       end do
    end do

    ! Reciprocal space kernel with b-spline discrete Fourier transform
    ! correction.
    !TODO: If b-spline Fourier coefficients are stored in a linearized
    ! array, then this could be simplified to a single loop.
    ixyz = 1
    do iz = 0, this%grid%localDimsR(3) - 1
       do iy = 0, this%grid%localDimsR(2) - 1
          byz = bsplineFourierCoeffY(iy + 1) * bsplineFourierCoeffZ(iz + 1 + this%grid%offsetR(3))
          do ix = 0, gridDimX_k - 1
             igk = 1 + ix + (iy + iz * this%grid%localDimsR(2)) &
                  * gridDimX_k
             bxyz = bsplineFourierCoeffX(ix + 1) * byz
             waveVector = kxi(ix+1,:) + kyi(iy+1,:) + kzi(iz+1,:)
             k2 = dot_product(waveVector, waveVector)
             ! Note that bxyz is not squared here like it would normally be.
             kernel(ixyz) = (4 * pi / k2) * gaussianFourierCoeff(ixyz) / bxyz
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
#if defined(MPI)
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
#endif /*defined(MPI)*/

    ! Initialize grids.
    ucu = 0
    uuv1d_r = 0
    uuv1d_c = 0
    this%uuv1d = 0

    ! Spread the Gaussian charges to the grid using b-splines.
    do iu = 1, this%solute%numAtoms
#if defined(MPI)
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
#if defined(MPI)
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

#if defined(MPI)
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

    !TODO: The normalization, complex conjugate, and such done by the
    ! RISM FFT module causes problems. Need to use FFTW directly for
    ! now.

    ! Evaluate the recip-space potential at the grid points.
    ! outr = 0
    uuv1d_r = 0
#if defined(MPI)
    call fftw_mpi_execute_dft_c2r(planbwd, uuv1d_c, uuv1d_r)
#else
    call fftw_execute_dft_c2r(planbwd, uuv1d_c, uuv1d_r)
#endif


    chargeCorrection = - pi / zeta * this%solute%totalCharge / this%grid%boxVolume
    do igz = 0, local_N - 1
       do igy = 0, M - 1
          do igx = 0, L - 1
             igk_fort = 1 + igx + (igy + igz * M) * L
             igk_cpp = 1 + igz + (igy + igx * M) * local_N
             igk_fort_global = 1 + igx &
                  + (igy + (igz + this%grid%offsetR(3)) * M) * L
             igk_cpp_global = 1 + (igz + this%grid%offsetR(3)) &
                  + (igy + igx * M) * N
             ! ucu(igx + 1, igy + 1, igz + 1, :) = igk_fort_global
             ! ucu(igx + 1, igy + 1, igz + 1, :) = 5
             ucu(igx + 1, igy + 1, igz + 1, :) = uuv1d_r(igx + 1, igy + 1, igz + 1) &
                  / this%grid%boxVolume + chargeCorrection
          end do
       end do
    end do

    call rism_timer_stop(this%coulombLongRangeTimer)
    call timer_stop(TIME_UCOULULR)

    ! Deallocate the FFT plans.
    call fftw_destroy_plan(planfwd)
    call fftw_destroy_plan(planbwd)

    call timer_start(TIME_UCOULUSR)
    call rism_timer_start(this%coulombShortRangeTimer)

    ! For those grid points near the atoms, replace the Gaussian potential
    ! with the point-charge potential.
    !TODO: This can be shared with the EwaldSum routine.
    call uvEwaldSumShortRangePotential(this, ucu)

    call rism_timer_stop(this%coulombShortRangeTimer)
    call timer_stop(TIME_UCOULUSR)

    ! Calculate electrostatic potential energy on each grid point.
    ! No 1/2 term is required since solute affects solvent, but
    ! solvent does not affect solute, hence no double counting.
    do iv = this%solvent%numAtomTypes, 1, -1
       if (index(this%solvent%atomName(iv), "-") == 0 .and. &
            index(this%solvent%atomName(iv), "+") == 0) then
          ucu(:,:,:,iv) = ucu(:,:,:,1) * this%solvent%charge(iv)
       else
          ucu(:,:,:,iv) = (ucu(:,:,:,1) + this%biasPotential) * this%solvent%charge(iv)
       end if
    end do
    ! TODO: Why is this done?
    ucu = ucu

    !TODO: Allocate and deallocate memory outside this function to
    ! avoid waste during force calculations.

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
#if defined(MPI)
    if (this%grid%mpirank == 0) then
       call fftw_free(uuv1d_final_cptr)
    end if
#endif
    if (safemem_dealloc(kernel) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
  end subroutine uvParticleMeshEwaldPotential


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

    half_xappa_smear = 0.5d0 * this%solvent%xappa * this%solvent%smear
    tcfLRasymp = 0
    if (.not. this%solvent%ionic) return
    do iu = 1, this%solute%numAtoms
       distance = sqrt(sum((r(:) - this%solute%position(:, iu))**2))
       if (distance > this%cut_hlr) then
          tcfLRasymp = tcfLRasymp &
               - this%solute%charge(iu) / distance
       else if (distance > 0) then
          tcfLRasymp = tcfLRasymp &
               - this%solute%charge(iu) / distance &
               * (exp(-this%solvent%xappa * distance) &
               * erfc(half_xappa_smear - distance / this%solvent%smear) &
               - exp(this%solvent%xappa * distance) &
               * erfc(half_xappa_smear + distance / this%solvent%smear)) / 2d0
       else ! distance == 0d0
          tcfLRasymp = tcfLRasymp &
               - this%solute%charge(iu) &
               * (2d0 / (sqrt(PI) * this%solvent%smear) &
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
       if (distance > this%cut_clr) then
          dcfLRasymp = dcfLRasymp &
               - this%solute%charge(iu) / distance
       else if (distance > 0) then
          dcfLRasymp = dcfLRasymp &
               - this%solute%charge(iu) / distance * erf(distance / this%solvent%smear)
       else ! distance == 0d0
          dcfLRasymp = dcfLRasymp &
               - this%solute%charge(iu) / (sqrt(PI) * this%solvent%smear) * 2d0
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
  function minimumImage(this, position)
   implicit none
    _REAL_ :: minimumImage(3)
    type(rism3d_potential), intent(in) :: this !< potential object.
    _REAL_, intent(in) :: position(3) !< Position vector.
    !! The vector origin is usually the grid point a calculation is
    !! performed at.

    integer :: id

    !TODO: Calculating minimum image of each solute site at each grid
    ! point twice (once for Ewald sum and once LJ) seems inefficient,
    ! though storing and retrieving the data may be about as fast as
    ! calculating it again from scratch.

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

end module rism3d_potential_c
