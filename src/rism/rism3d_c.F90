!<compile=optimized>
! The 3D-RISM-KH software found here is copyright (c) 2010-2012 by
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

!> 3D-RISM solver.
!! This defines the 3D-RISM type and associated subroutines.  All type
!! elements are public.  In general, read but do not write these
!! variables.  This provides an object-orientiented interface without
!! needing a function to access every variable.
!! Features of this solver include:
!! o Multiple closures w/ temperature derivatives: KH, HNC, PSE-n
!! o Temperature derivative expressed as T*d/dT
!! o MDIIS accelerated solutions
!! o Optional cutoffs
!! o Supercell method for long range asymptotics
!! o Analytic forces
!! o Variable grid size and dynamic memory allocation
!! o MPI support
!! o units:  energy       [kT]
!!           distances    [A]         (Angstroms)
!!           site charges [sqrt(kT A)]
!!           temperature  [K]
!!           density      [#/A^3]
!!           mass         [au]
!! o To convert [e] to [sqrt(kT A)] * sqrt(COULOMB_CONST_E/ KB / temperature)
module rism3d_c
  use rism3d_solute_c
  use rism3d_solvent_c
  use rism3d_potential_c
  use rism3d_grid_c
  use rism3d_closure_c
  use rism_report_c
  use rism_timer_c
  use mdiis_c
  use rism3d_fft_c

  use rism3d_opendx
#ifdef RISM3D_DEBUG
  !    use rism3d_debug_c
#endif
  implicit none
#include "def_time.h"

  type rism3d
     !! Solute/solvent information.

     !> Solute object.
     type(rism3d_solute) :: solute
     !> Solvent object.
     type(rism3d_solvent) :: solvent
     !> Potential object.
     type(rism3d_potential) :: potential
     !> Grid object.
     type(rism3d_grid) :: grid
     !> Closure object.
     type(rism3d_closure) :: closure

     !> List of closure names to use in order.  Only the last closure
     !! is used for thermodynamic output.  This can be used to
     !! progressively increase the order of the closure to aid
     !! convergence.
     character(len = 8), pointer :: closureList(:) => NULL()

     ! TIMERS.  Subtimers only account for computation.  We ignore setup etc.
     ! timer :: timer for this class.  Activated for all public routines
     ! resizeTimer :: time to resize solvent box
     ! reorientTimer :: time to reorient solute
     ! cuvpropTimer :: time to propagate Cuv solution
     ! fftTimer :: specifically times FFT calculation
     ! solveTimer :: specifically times rism1d_solve calculation
     ! solve3DRISMTimer :: specifically times solve3DRISM calculation
     ! single3DRISMsolutionTimer :: specifically times single3DRISMsolution calculation
     ! thermoTimer :: specifically times thermodynamics calculations
     ! forceTimer :: specifically times force calculation
     ! excessChemicalPotentialTimer :: specificall times excess chemical potential calculation
     ! solve3DRISM_dTTimer :: specifically times solve3DRISM_dT calculation
     ! single3DRISM_dTsolutionTimer :: specifically times single3DRISM_dTsolution calculation
     ! fft_dTTimer :: specifically times FFT calculation for temperature derivatives
     type(rism_timer) :: timer, resizeTimer, reorientTimer, &
          cuvpropTimer, fftTimer, solventTimer, &
          solve3DRISMTimer, single3DRISMsolutionTimer, thermoTimer, &
          forceTimer, excessChemicalPotentialTimer, &
          solve3DRISM_dTTimer, single3DRISM_dTsolutionTimer, fft_dTTimer

     !! private !(should be)

     ! FFTW options

     ! FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
     integer :: fftw_planner = FFT_ESTIMATE
     ! .true.  - use aligned memory and to enable SIMD;
     ! .false. - don't use aligned memory
     logical :: fft_aligned = .true.
     ! Transpose site number and spatial data locally before and after FFT.
     logical :: fftw_localtrans = .true.

     !> Output verbosity.  Useful for debugging.
     !! 0 - no ouput
     !! 1 - memory allocation and steps for convergence
     !! 2 - 1 + convergence progress
     integer :: verbose = 0

     ! This is a bit ugly and there may be a better solution.  We
     ! need to keep track of the number of solutions for both charged
     ! and un-charged solutes.  When we change between the solutes we
     ! set the nsolutions pointer to the appropriate variable.
     ! However, the 'target' attribute is not allowed in type
     ! definitions so these variables have to be pointers and we have
     ! to allocate memory for them.

     !> Number of times full solutions have been calculated.
     integer, pointer :: nsolution => NULL()
     !> Number of times full solutions with a charged solute have been
     !! calculated.
     integer, pointer :: nsolutionChg => NULL()
     !> Number of times full solutions with an uncharged solute have
     !! been calculated.
     integer, pointer :: nsolutionNoChg => NULL()

     !> Center the solute in the solvation box.
     !! 0 - off
     !! 1 - center of mass
     !! 2 - center of geometry
     !! 3 - center of mass shifted to the nearest grid point
     !! 4 - center of geometry shifted to the nearest grid point
     !! For negative numbers the centering translation is only
     !! calculated the for the first solution and used for subsequent
     !! calculations.  This allows the solute to drift in the box.
     integer :: centering = 1

     !> Number of past direct correlation function time step saves.
     integer :: ncuvsteps ! numDCFsteps

     !> Buffer distance to the edge of the box for the solvent. [A]
     _REAL_ :: buffer = 12d0
     !> Fixed box size for 3D-RISM.
     _REAL_ :: fixedBoxDimensionsR(3)
     !> Number of Cartesian grid points in each dimension for a fixed box size.
     integer :: fixedNumGridPoints(3)

     !> Variable box size.
     logical :: varbox = .true.

     !> Number of vectors used for MDIIS (consequently, the number of
     !! copies of CUV we need to keep for MDIIS).
     integer :: NVec
     !> MDIIS implementation to use.
     integer :: mdiis_method
     type(mdiis) :: mdiis_o

     !> 'Step size' for MDIIS.
     _REAL_ :: deloz = 0.7d0
     !> Restart threshold factor. Ratio of the current residual to the
     !! minimum residual in the basis that causes a restart.
     _REAL_ :: mdiis_restart

     !! MPI Support !!
     integer :: mpirank = 0, mpicomm = 0, mpisize = 1

     !! LARGE ARRAYS !!
     !
     ! all arrays are declared as pointers to ensure we can reallocate them as necessary
     !

     ! xvva       :: solvent chi interpolated for our grid size
     ! guv        :: solvent distribution function
     ! huv        :: guv - 1
     ! cuv        :: solvent direct correlation function and points to the
     !              current active solution in cuvWRK
     ! cuvres     :: residual value for cuv calculation and points to the
     !              current active solution in cuvresWRK.
     ! cuvWRK     :: Working Cuv memory.  Holds Cuv from previous iterations.
     ! cuvresWRK  :: Working Cuvres memory.  Holds Cuvres from previous iterations.
     ! oldcuv     :: previous solutions of cuv. Points to oldcuvChg
     !              or oldcuvNoChg depending on the charge state of the
     !              calculation.
     ! oldcuvChg  :: previous solutions for the standard charged system
     ! oldcuvNoChg :: previous solutions for the chargeless system.  This is only allocated
     !                if _unsetCharges() is called
     ! electronMap :: smeared solvent electron density map.
     _REAL_, pointer :: xvva(:) => NULL(), &
          oldcuv(:, :, :, :, :) => NULL(), &
          oldcuvChg(:, :, :, :, :) => NULL(), &
          oldcuvNoChg(:, :, :, :, :) => NULL(), &
          cuv(:, :, :, :) => NULL(), cuvWRK(:, :, :, :, :) => NULL(), &
          cuvres(:, :) => NULL(), cuvresWRK(:, :, :) => NULL()


     _REAL_, pointer :: guv(:, :) => NULL(), huv(:, :) => NULL()

     _REAL_, pointer :: electronMap(:, :, :) => NULL()

     ! Temperature derivative memory
     ! cuvk        :: k-space Cuv solution from 3D-RISM
     !               solution. NOTE: we should consider using Huv or
     !               Guv memory instead.  However, it has to be
     !               checked first that it is not used for any thermodynamics calculations
     ! xvva_dT     :: solvent dT chi interpolated for our grid size
     ! guv_dT      :: solvent dT distribution function
     ! huv_dT      :: guvdT
     ! cuv_dT      :: solvent dT direct correlation function and points to the
     !              current active solution in cuvWRK
     ! cuvres_dT   :: residual value for dT cuv calculation and points to the
     !              current active solution in cuvresWRK.
     ! cuvWRK_dT   :: Working dT Cuv memory.  Holds Cuv from previous iterations.
     ! cuvresWRK_dT:: Working dT Cuvres memory.  Holds Cuvres from previous iterations.
     _REAL_, pointer ::  xvva_dT(:) => NULL(), &
          cuv_dT(:, :, :, :) => NULL(), cuvres_dT(:, :) => NULL(), &
          cuvWRK_dT(:, :, :, :, :) => NULL(), cuvresWRK_dT(:, :, :) => NULL()

     _REAL_, pointer :: cuvk(:, :) => NULL(), guv_dT(:, :) => NULL(), huv_dT(:, :) => NULL()

     ! fft :: fft object for standard 3D-RISM solution
     ! fft_dT :: fft object for temperature derivative 3D-RISM solution
     ! fft_cuv :: fft object for Cuv, part of the temperature derivative code
     type(rism3d_fft) :: fft, fft_dT, fft_cuv

     !> If true, a periodic 3D-RISM calculation is performed. This
     !! primarily differs from infinite dilution 3D-RISM by using
     !! Ewald sum potential in place of Coulombic potential and
     !! invoking the minimum image convention while calculating both
     !! the Ewald sum and Lennard-Jones potentials.
     logical :: periodic = .false.

     !> Lengths and interior angles of the unit cell. For aperiodic
     !! systems, the interior angles are always 90 degrees.
     _REAL_ :: unitCellDimensions(6)

     !> The abbreviated label of the periodic potential function used
     !! for periodic calculations. See rism3d_potential for valid values.
     character(len=255) :: periodicPotential = ""

  end type rism3d

  public :: rism3d_new, rism3d_destroy, rism3d_calculateSolution, rism3d_force, &
       rism3d_excessChemicalPotential_tot, rism3d_excessChemicalPotential, &
       rism3d_excessChemicalPotentialGF_tot, rism3d_excessChemicalPotentialGF, &
       rism3d_setbox_variable, rism3d_setbox_fixed, &
       rism3d_setclosure, rism3d_setverbosity, rism3d_setcut, rism3d_setmdiis

  private :: resizeBox, reallocateBox, reallocateBox_dT, &
       interpolateSolventSusceptibility, centerSolute, &
       solve3DRISM, single3DRISMsolution, &
       solve3DRISM_dT, single3DRISMsolution_dT, &
       lennardJonesForce, coulombicForce, guessDCF, updateDCFguessHistory

contains


  !> Constructor - precalculates the solute solvent terms that are not
  !! configuration dependent and sets box parameters.
  !!
  !! The solvation box may be fixed size or variable.  For fixed size,
  !! define o_boxlen and o_ng3.  For variable box size, define buffer
  !! and grdspc.  Do not mix these parameters as this will cause the
  !! program to halt.
  !!
  !! For periodic simulations, buffer can be set but will be ignored;
  !! (this should get fixed).  The unit cell parameters give the size of
  !! the box, and the grdspc(1:3) array gives an approximate grid
  !! spacing.  The actual spacing will be set that an exact number of 
  !! grids spans the box, and so that the number of grid points is even.  
  !! (In addtion, for MPI runs, the number of grids along y and z will 
  !! be adjusted to be a multiple of the number of MPI threads.)
  !!
  !! If this is an MPI run, supply the MPI communicator.  Only the rank
  !! 0 parameters will be used. However, due to the limitations of
  !! pre-Fortran2003, the closure must be the same length on all
  !! processes. The values and number of elements for the closure list
  !! on non-rank 0 processes still do not matter.
  !!
  !! IN:
  !!   this :: new rism3d object
  !!   solu :: 3D-RISM solute object
  !!   solv :: 3D-RISM solvent object
  !!   centering :: center the solute in the solvation box.
  !!   ncuvsteps :: number of past cuv time steps saves
  !!   closure :: list of closures. Closures may be KH, HNC or PSEn
  !!              where n is an integer. Ensure the length attribute is
  !!              the same on all processes.
  !!   cut     :: distance cutoff for potential and force calculations
  !!   mdiis_nvec :: number of MDIIS vectors (previous iterations) to keep
  !!   mdiis_del :: scaling factor (step size) applied to estimated gradient (residual)
  !!   mdiis_method :: which implementation of the algorithm
  !!   o_buffer :: (optional) shortest distance between solute and solvent box boundary
  !!   o_grdspc :: (optional) linear grid spacing for the solvent box in each dimension
  !!   o_boxlen :: (optional) solvent box size in each dimension [A]
  !!   o_ng3    :: (optional) number of grid points in each dimension
  !!   o_mpicomm :: (optional) MPI communicator
  !!   o_periodic :: (optional) periodic electric potential to use, if any
  !!   o_unitCellDimensions :: (optional) geometry of the system unit cell

  subroutine rism3d_new(this, solute, solvent, centering, ncuvsteps, &
       closure, cut, mdiis_nvec, mdiis_del, mdiis_method, mdiis_restart, &
       o_buffer, o_grdspc, o_boxlen, o_ng3, o_mpicomm, &
       o_periodic, o_unitCellDimensions, o_biasPotential)
    use rism3d_solute_c
    use rism3d_solvent_c
    use safemem
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    type(rism3d), intent(inout) :: this
    type(rism3d_solute), intent(in), target :: solute
    type(rism3d_solvent), intent(in), target :: solvent
    integer, intent(in) :: centering, ncuvsteps
    character(len = *), intent(in) :: closure(:)
    _REAL_, intent(in) :: cut
    integer, intent(in) :: mdiis_nvec, mdiis_method
    _REAL_, intent(in) :: mdiis_del, mdiis_restart
    _REAL_, optional, intent(in) :: o_buffer, o_grdspc(3)
    _REAL_, optional, intent(in) :: o_boxlen(3)
    integer, optional, intent(in) :: o_ng3(3)
    integer, optional, intent(in) :: o_mpicomm
    character(len = *), optional, intent(in) :: o_periodic
    _REAL_, optional, intent(in) :: o_unitCellDimensions(6)
    _REAL_, optional, intent(in) :: o_biasPotential
    ! temporary copies
    character(len = len(closure)), pointer :: t_closure(:)
    _REAL_ :: t_cut
    integer :: t_mdiis_nvec, t_mdiis_method
    _REAL_ :: t_mdiis_del, t_mdiis_restart
    _REAL_ :: t_buffer, t_grdspc(3)
    _REAL_ :: t_boxlen(3)
    integer :: t_ng3(3)
    _REAL_ :: t_unitCellDimensions(6)
    integer :: t_mpicomm
    integer :: nclosure
    integer :: err

    ! MPI set up starts by obtaining rank and size.  Temporary copies
    ! of input parameters that do not directly set object variables
    ! are made.  These are they broadcast to the rank > 0 processes.
    ! Then all processes complete the intitialization proceedure
    ! using the temporary copies.  This leaves the input parameters
    ! untouched.

    call rism_timer_new(this%timer, "3D-RISM")
    call rism_timer_start(this%timer)
    call rism_timer_new(this%thermoTimer, "Thermodynamics")
    call rism_timer_setParent(this%thermoTimer, this%timer)
    call rism_timer_new(this%forceTimer, "Force")
    call rism_timer_setParent(this%forceTimer, this%thermoTimer)
    call rism_timer_new(this%excessChemicalPotentialTimer, "Excess Chemical Potential")
    call rism_timer_setParent(this%excessChemicalPotentialTimer, this%thermoTimer)
    call rism_timer_new(this%solventTimer, "Solve 3D-RISM")
    call rism_timer_setParent(this%solventTimer, this%timer)
    call rism_timer_new(this%resizeTimer, "Solvation box resize")
    call rism_timer_setParent(this%resizeTimer, this%solventTimer)
    call rism_timer_new(this%reorientTimer, "Solute reorientation")
    call rism_timer_setParent(this%reorientTimer, this%solventTimer)
    call rism_timer_new(this%cuvpropTimer, "Cuv propagation")
    call rism_timer_setParent(this%cuvpropTimer, this%solventTimer)
    call rism_timer_new(this%solve3DRISMTimer, "RXRISM")
    call rism_timer_setParent(this%solve3DRISMTimer, this%solventTimer)
    call rism_timer_new(this%single3DRISMsolutionTimer, "R1RISM")
    call rism_timer_setParent(this%single3DRISMsolutionTimer, this%solve3DRISMTimer)
    call rism_timer_new(this%fftTimer, "FFT")
    call rism_timer_setParent(this%fftTimer, this%single3DRISMsolutionTimer)
    call rism_timer_new(this%solve3DRISM_dTTimer, "RXRISM_dT")
    call rism_timer_setParent(this%solve3DRISM_dTTimer, this%solventTimer)
    call rism_timer_new(this%single3DRISM_dTsolutionTimer, "R1RISM")
    call rism_timer_setParent(this%single3DRISM_dTsolutionTimer, this%solve3DRISM_dTTimer)
    call rism_timer_new(this%fft_dTTimer, "FFT")
    call rism_timer_setParent(this%fft_dTTimer, this%single3DRISM_dTsolutionTimer)

    nullify(t_closure)

    if (present(o_periodic)) then
       if (o_periodic /= '') then
          this%periodic = .true.
          this%periodicPotential = o_periodic
       end if
    end if

    ! GET RANK AND SIZE
    this%mpicomm = 0
    this%mpisize = 1
    this%mpirank = 0
#ifdef MPI
    if (present(o_mpicomm)) then
       this%mpicomm = o_mpicomm
       if (this%mpicomm == MPI_COMM_NULL) &
            call rism_report_error("RISM3D: received NULL MPI communicator")
       call mpi_comm_rank(this%mpicomm, this%mpirank, err)
       if (err /= 0) call rism_report_error &
            ("(a,i8)", "RISM3D: could not get MPI rank for communicator ", this%mpicomm)
       call mpi_comm_size(this%mpicomm, this%mpisize, err)
       if (err /= 0) call rism_report_error &
            ("(a,i8)", "RISM3D: could not get MPI size for communicator ", this%mpicomm)
       call rism_report_mpi(this%mpicomm)
    end if
#endif /*MPI*/
    ! MAKE TEMPORARY COPIES
    if (this%mpirank == 0) then
       call rism3d_solute_clone(solute, this%solute)
       call rism3d_solvent_clone(solvent, this%solvent)
       this%centering = centering
       this%ncuvsteps = ncuvsteps
       nclosure = size(closure)
       t_closure => safemem_realloc(t_closure, len(t_closure), nclosure)
       t_closure = closure
       t_cut = cut
       t_mdiis_nvec = mdiis_nvec
       t_mdiis_method = mdiis_method
       t_mdiis_del = mdiis_del
       t_mdiis_restart = mdiis_restart
       ! check box parameters
       if (present(o_buffer) .and. present(o_grdspc)) then
              !! (This branch should also be taken for periodic RISM)
          t_buffer = o_buffer
          t_grdspc = o_grdspc
          if (present(o_boxlen) .or. present(o_ng3)) &
               call rism_report_error("RISM3D: do not set BOXLEN or NG3 for variable box size")
       else if (present(o_boxlen) .and. present(o_ng3)) then
          t_boxlen = o_boxlen
          t_ng3 = o_ng3
          if (present(o_buffer) .or. present(o_grdspc)) &
               call rism_report_error("RISM3D: do not set BUFFER or GRDSPC for fixed box size")
       else
          call rism_report_error("RISM3D: not enough parameters for fixed or variable box size")
       end if
       if (present(o_unitCellDimensions)) then
          t_unitCellDimensions = o_unitCellDimensions
       end if
    end if
#ifdef MPI
    ! BROADCAST PARAMETERS
    ! set solu on all processes
    call rism3d_solute_mpi_clone(this%solute, this%mpirank, this%mpicomm)
    ! set solv on all processes
    call rism3d_solvent_mpi_clone(this%solvent, this%mpirank, this%mpicomm)
    ! set centering on all processes
    call mpi_bcast(this%centering, 1, mpi_integer, 0, this%mpicomm, err)
    if (err /= 0) call rism_report_error("RISM3D: broadcast CENTERING in constructor failed")
    ! set ncuvstpes on all processes
    call mpi_bcast(this%ncuvsteps, 1, mpi_integer, 0, this%mpicomm, err)
    if (err /= 0) call rism_report_error("RISM3D: broadcast NCUVSTEPS in constructor failed")
    call mpi_bcast(nclosure, 1, mpi_integer, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error &
         ("RISM3D interface: could not broadcast PROGRESS")
    if (this%mpirank/=0) &
         t_closure => safemem_realloc(t_closure, len(t_closure), nclosure)
    call mpi_bcast(t_closure, len(t_closure) * nclosure, mpi_character, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast CLOSURE in constructor failed")
    call mpi_bcast(t_cut, 1, mpi_double_precision, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast CUT in constructor failed")
    call mpi_bcast(t_mdiis_nvec, 1, mpi_integer, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast MDIIS_NVEC in constructor failed")
    call mpi_bcast(t_mdiis_del, 1, mpi_double_precision, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast MDIIS_DEL in constructor failed")
    call mpi_bcast(t_mdiis_restart, 1, mpi_double_precision, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast MDIIS_RESTART in constructor failed")
    call mpi_bcast(t_mdiis_method, 1, mpi_integer, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast MDIIS_METHOD in constructor failed")
    if (present(o_buffer) .and. present(o_grdspc)) then
       call mpi_bcast(t_buffer, 1, mpi_double_precision, 0, this%mpicomm, err)
       if (err /=0) call rism_report_error("RISM3D: broadcast BUFFER in constructor failed")
       call mpi_bcast(t_grdspc, 3, mpi_double_precision, 0, this%mpicomm, err)
       if (err /=0) call rism_report_error("RISM3D: broadcast GRDSPC in constructor failed")
    else
       call mpi_bcast(t_boxlen, 3, mpi_double_precision, 0, this%mpicomm, err)
       if (err /=0) call rism_report_error("RISM3D: broadcast BOXLEN in constructor failed")
       call mpi_bcast(t_ng3, 3, mpi_integer, 0, this%mpicomm, err)
       if (err /=0) call rism_report_error("RISM3D: broadcast NG3 in constructor failed")
    end if
    if (present(o_unitCellDimensions)) then
       call mpi_bcast(t_unitCellDimensions, 6, mpi_double_precision, 0, this%mpicomm, err)
       if (err /=0) call rism_report_error("RISM3D: broadcast UNITCELLDIMENSIONS in constructor failed")
    end if
#endif /*MPI*/
    ! INITIALIZE
    call rism3d_grid_new(this%grid, this%mpicomm)
    call rism_timer_stop(this%timer)
    call rism3d_setmdiis(this, t_mdiis_nvec, t_mdiis_del, t_mdiis_method, t_mdiis_restart)
    call rism_timer_start(this%timer)
    call rism3d_potential_new(this%potential, this%grid, this%solvent, this%solute, 0d0, &
         this%fft, this%periodicPotential, o_biasPotential)
    call rism3d_potential_setTimerParent(this%potential, this%solventTimer)

#ifdef MPI
    call mdiis_new_mpi(this%mdiis_o, this%mdiis_method, &
         this%deloz, 0d0, &
         this%MDIIS_restart, &
         this%mpirank, this%mpisize, this%mpicomm)
#else
    call mdiis_new(this%mdiis_o, this%mdiis_method, &
         this%deloz, 0d0, &
         this%MDIIS_restart)
#endif /*MPI*/

    call mdiis_setTimerParent(this%mdiis_o, this%single3DRISMsolutiontimer)
    call rism_timer_stop(this%timer)
    call rism3d_setcut(this, t_cut)
    call rism3d_setclosurelist(this, t_closure)
    if (present(o_buffer) .and. present(o_grdspc)) then
       call rism3d_setbox_variable(this, t_buffer, t_grdspc)
    else
       call rism3d_setbox_fixed(this, t_boxlen, t_ng3)
    end if

    if (present(o_unitCellDimensions)) then
       !TODO: This can probably be made a local variable.
       this%unitCellDimensions = t_unitCellDimensions
       call rism3d_grid_setUnitCellDimensions(this%grid, this%unitCellDimensions, this%periodic)
    end if

    allocate(this%nsolutionChg, this%nsolutionNoChg)
    this%nsolutionChg = 0
    this%nsolutionNoChg = 0
    this%nsolution => this%nsolutionChg

    call rism3d_fft_global_init()
    this%fftw_planner = FFT_MEASURE
#if defined(MPI)
    this%fft_aligned = .false.
#else
    this%fft_aligned = .true.
#endif
    this%fftw_localtrans = .true.

    ! Clean up locally allocated temporary memory.
    if (safemem_dealloc(t_closure) /= 0) &
         call rism_report_error("RISM3D:NEW: failed to deallocate t_closure")
#ifdef RISM3D_DEBUG
    call rism3d_debug_new(this%grid, this%solvent, this%mpirank, this%mpisize, this%mpicomm)
#endif

  end subroutine rism3d_new


  !> Check if we can perform a temperature derivative calculation
  !! (i.e. all the necessary information is available and the closure
  !! supports it).
  !! IN:
  !!   this : rism3d object
  !! OUT:
  !!    .true. if we can, .false. if we can't
  function rism3d_canCalc_DT(this) result(can_dT)
    implicit none
    type(rism3d), intent(in) :: this
    logical :: can_dT

    can_dT = rism3d_solvent_canCalc_dT(this%solvent) .and. &
         rism3d_closure_canCalc_dT(this%closure)
  end function rism3d_canCalc_DT


  !> Set parent for this timer
  !! IN:
  !!   this : rism3d object
  !!   parent : parent timer object
  subroutine rism3d_setTimerParent(this, parent)
    implicit none
    type(rism3d), intent(inout) :: this
    type(rism_timer), intent(inout) :: parent
    call rism_timer_start(this%timer)
    call rism_timer_setParent(this%timer, parent)
    call rism_timer_stop(this%timer)
  end subroutine rism3d_setTimerParent


  !> Sets the parameters for a variable solvation box.
  !! IN:
  !!  this :: rism3d object
  !!  buffer :: shortest distance between solute and solvent box boundary
  !!  grdspc :: linear grid spacing for the solvent box in each dimension
  subroutine rism3d_setbox_variable(this, buffer, grdspc)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: buffer, grdspc(3)
    call rism_timer_start(this%timer)
    this%varbox = .true.
    this%buffer = buffer
    call rism3d_grid_setSpacing(this%grid, grdspc)
    call rism_timer_stop(this%timer)
  end subroutine rism3d_setbox_variable


  !>Sets the parameters for a fixed solvation box.
  !! IN:
  !!  this :: rism3d object
  !!  boxlen :: solvent box size in each dimension in Angstroms
  !!  ng3    :: number of grid points in each dimension
  subroutine rism3d_setbox_fixed(this, boxlen, ng3)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: boxlen(3)
    integer, intent(in) :: ng3(3)
    call rism_timer_start(this%timer)
    this%varbox = .false.
    this%fixedBoxDimensionsR = boxlen
    this%fixedNumGridPoints = ng3
    call rism_timer_stop(this%timer)
  end subroutine rism3d_setbox_fixed


  !> Sets the closure list and sets the current closure to the first one
  !! in the list.  When there is no previous solution to work from, the
  !! solver will use each closure in the list in turn. By choosing the
  !! list to increase in order, it makes it possible to converge
  !! otherwise difficult closures. Only the last closure is used for
  !! thermodynamic output.
  !! IN:
  !!   this :: rism3d object
  !!   closure :: array of closure types (see closure enumeration).
  subroutine rism3d_setclosurelist(this, closure)
    implicit none
    type(rism3d), intent(inout) :: this
    character(len = *), intent(in) :: closure(:)
    call rism_timer_start(this%timer)
    this%closureList => safemem_realloc(this%closureList, len(this%closureList), &
         ubound(closure, 1))
    this%closureList = closure
    call rism_timer_stop(this%timer)
    call rism3d_setclosure(this, this%closureList(1))
  end subroutine rism3d_setclosurelist


  !> Sets the closure type.
  !! IN:
  !!   this :: rism3d object
  !!   closure :: closure type (see closure enumeration).
  subroutine rism3d_setclosure(this, closure)
    implicit none
    type(rism3d), intent(inout) :: this
    character(len = *), intent(in) :: closure
    call rism_timer_start(this%timer)
    call rism3d_closure_destroy(this%closure)
    call rism3d_closure_new(this%closure, closure, this%potential)
    call rism_timer_stop(this%timer)
  end subroutine rism3d_setclosure


  !> Sets verbosity of output.
  !! IN:
  !!   this :: rism3d object
  !!   verbosity :: 0 - no output
  !!                1 - memory allocation and steps for convergence
  !!                2 - 1 + convergence progress
  subroutine rism3d_setverbosity(this, verbosity)
    implicit none
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: verbosity
    call rism_timer_start(this%timer)
    this%verbose = verbosity
    call rism_timer_stop(this%timer)
  end subroutine rism3d_setverbosity


  !> Sets the cut off distance for potential and force calculations.
  !! IN:
  !!   this :: rism3d object
  !!   cut     :: distance cutoff for potential and force calculations
  subroutine rism3d_setcut(this, cut)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: cut
    call rism_timer_start(this%timer)
    call rism3d_potential_setCut(this%potential, cut)
    call rism_timer_stop(this%timer)
  end subroutine rism3d_setcut


  !> Sets MDIIS parameters
  !! IN:
  !!   this :: rism3d object!
  !!   nvec :: number of MDIIS vectors (previous iterations) to keep
  !!   del :: scaling factor (step size) applied to estimated gradient (residual)
  !!   method :: which implementation of the algorithm
  !!   restart :: restart threshold factor. Ratio of the current residual to the
  !!              minimum residual in the basis that causes a restart
  subroutine rism3d_setmdiis(this, nvec, del, method, restart)
    implicit none
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: nvec, method
    _REAL_, intent(in) :: del, restart
    call rism_timer_start(this%timer)
    this%NVec = nvec
    this%deloz = del
    this%mdiis_method = method
    this%mdiis_restart = restart
    call rism_timer_stop(this%timer)
  end subroutine rism3d_setmdiis


  !> Sets solute coordinates.
  !! IN:
  !!   this :: rism3d object
  !!   ratu :: coordinates
  subroutine rism3d_setCoord(this, solutePositions)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: solutePositions(:, :)
    call rism_timer_start(this%timer)
    call rism3d_solute_setCoord(this%solute, solutePositions)
    call rism_timer_stop(this%timer)
  end subroutine rism3d_setCoord


  !> Sets all solute partial charges to zero, resets MDIIS and wipes out
  !! working memory.
  !! IN:
  !!   this :: rism3d object
  subroutine rism3d_unsetCharges(this)
    implicit none
    type(rism3d), intent(inout) :: this
    integer :: i
    ! reset MDIIS.  This makes the working vector index 1
    call mdiis_reset(this%mdiis_o)
    this%cuv => this%cuvWRK(:, :, :, :, mdiis_getWorkVector(this%mdiis_o))
    this%cuvres => this%cuvresWRK(:, :, mdiis_getWorkVector(this%mdiis_o))
    ! turn off charges
    call rism3d_solute_unsetCharges(this%solute)
    ! Use the number of no charge solutions
    this%nsolution => this%nsolutionNoChg
    ! Use no charge previous soluitions
    this%oldcuv => this%oldcuvNoChg
    ! If we have run with no charges before, copy the previous solution
    if (associated(this%oldcuvNoChg)) &
         call dcopy(product(ubound(this%oldcuv)), this%oldcuv, 1,this%cuv, 1)
  end subroutine rism3d_unsetCharges


  !> Sets all solute partial charges to to their original
  !! values. (Undoes rism3d_unsetCharge().)
  !! IN:
  !!   this :: rism3d object
  subroutine rism3d_resetCharges(this)
    implicit none
    type(rism3d), intent(inout) :: this
    integer :: i
    ! reset MDIIS.  This makes the working vector index 1
    call mdiis_reset(this%mdiis_o)
    this%cuv => this%cuvWRK(:, :, :, :, mdiis_getWorkVector(this%mdiis_o))
    this%cuvres => this%cuvresWRK(:, :, mdiis_getWorkVector(this%mdiis_o))
    ! get back the charges
    call rism3d_solute_resetCharges(this%solute)
    ! restore the number of previous solutions
    this%nsolution => this%nsolutionChg
    ! point to previous charged solutions
    this%oldcuv => this%oldcuvChg
    ! If we have run with charges before, copy the previous solution
    if (associated(this%oldcuvNoChg)) &
         call dcopy(product(ubound(this%oldcuv)), this%oldcuv, 1,this%cuv, 1)
  end subroutine rism3d_resetCharges


  !> Calculates the full 3D-RISM solvent distribution.  This is required to
  !! calculate thermodynamic quantities.
  !! @param[in,out] this rism3d object.
  !! @param[in,out] ksave Save intermediate results every ksave
  !!            interations (0 means no saves).
  !! @param[in] kshow Print parameter for relaxation steps every kshow
  !!            iteration (0 means no print).
  !! @param[in] maxSteps Maximum number of rism relaxation steps.
  !! @param[in] tolerance Convergence tolerances. There should be one
  !!          tolerance per closure in the closure list.
  subroutine rism3d_calculateSolution(this, ksave, kshow, maxSteps, tolerance)
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: ksave, kshow, maxSteps
    _REAL_, intent(in) :: tolerance(:)
    ! iclosure :: counter for closures
    integer :: iclosure

    _REAL_ :: offset(3)
    integer :: id, iu

    call rism_timer_start(this%solventTimer)

    ! 0) Quick check that the tolerance list is of the correct length.
    if (ubound(tolerance, 1) /= ubound(this%closureList, 1)) &
         call rism_report_error("(a,i3,a,i3)", &
         "RISM3D_SOLVE: number of tolerances, ", &
         ubound(tolerance, 1), ", is not equal to numer of closures, ", &
         ubound(this%closureList, 1))


    ! 1) Reorient solute along the principal axis and resize the grids
    ! if necessary.

    ! 1b) Get the minimum box size for this frame.
    if (this%periodic .or. this%varbox .or. this%nsolution == 0) then
       call rism_timer_start(this%resizeTimer)
       call timer_start(TIME_RESIZE)
       call resizeBox(this)
       call timer_stop(TIME_RESIZE)
       call rism_timer_stop(this%resizeTimer)
    end if
    ! By default, center solute in box for the infinite dilution case,
    ! but not for the periodic case since the origin helps define
    ! rotation axes.
    call timer_start(TIME_REORIENT)
    call rism_timer_start(this%reorientTimer)
    call centerSolute(this)
    call rism_timer_stop(this%reorientTimer)
    call timer_stop(TIME_REORIENT)

    ! 2) Calculate long range asymptotics of the direct and total
    ! correlation functions about the solute.
    call timer_start(TIME_ASYMP)
    call rism3d_potential_dcf_tcf_long_range_asymptotics(this%potential)
    call timer_stop(TIME_ASYMP)

    ! 3) Calculate electrostatic and Lennard-Jones potential about the
    ! solute.
    call rism3d_potential_calc(this%potential)


    ! 4) Propagate previously saved solute-solvent DCF solutions to
    ! create an initial guess for this solution.
#ifdef RISM_DEBUG
    write(0,*) "(4)", this%grid%totalLocalPointsK, this%grid%waveNumberArraySize
    call flush(0)
#endif /*RISM_DEBUG*/
    call timer_start(TIME_CUVPROP)
    call rism_timer_start(this%cuvpropTimer)
    call guessDCF(this)
    call rism_timer_stop(this%cuvpropTimer)
    call timer_stop(TIME_CUVPROP)


    ! 5) Calculate 3D-RISM solution using MDIIS.
#ifdef RISM_DEBUG
    write(0,*) "(5) RXRISM"
    call flush(0)
#endif /*RISM_DEBUG*/
    ! If the user has to provide a list of closures, use it only if
    ! this is the first solution (nsolution == 0) or solution
    ! propagation is turned off (ncuvsteps == 0). Otherwise, the
    ! current closure will be the last one in the list.
    if (this%nsolution == 0 .or. this%ncuvsteps == 0) then
       do iclosure = 1, size(this%closureList)
          if (this%verbose >= 1) &
               call rism_report_message("|Switching to "// &
               trim(this%closureList(iclosure))//" closure")
          call rism_timer_stop(this%solventTimer)
          call rism3d_setClosure(this, this%closureList(iclosure))
          call rism_timer_start(this%solventTimer)
          call timer_start(TIME_RXRISM)
          call solve3DRISM(this, ksave, kshow, maxSteps, tolerance(iclosure))
          call timer_stop(TIME_RXRISM)
          ! Increment nsolution and ncuvsteps to ensure the previous
          ! closure solution is used.
          if (iclosure == 1) then
             this%nsolution = this%nsolution + 1
             this%ncuvsteps = this%ncuvsteps + 1
          end if
       end do
       this%nsolution = this%nsolution - 1
       this%ncuvsteps = this%ncuvsteps - 1
    else
       call timer_start(TIME_RXRISM)
       call solve3DRISM(this, ksave, kshow, maxSteps, tolerance(size(tolerance)))
       call timer_stop(TIME_RXRISM)
    end if

    ! 11) Update stored variables.
    call timer_start(TIME_CUVPROP)
    call rism_timer_start(this%cuvpropTimer)
    this%nsolution = this%nsolution + 1
    call updateDCFguessHistory(this)
    call rism_timer_stop(this%cuvpropTimer)
    call timer_stop(TIME_CUVPROP)

    call rism_timer_stop(this%solventTimer)
  end subroutine rism3d_calculateSolution


  !> Calculates the full 3D-RISM solvent distribution.  This is required to
  !! calculate solvation energies and entropies.
  !! IN:
  !!   this :: rism3d object
  !!   ksave  :: save itermediate results every ksave interations (0 means no saves)
  !!   kshow  :: print parameter for relaxation steps every kshow iteration (0 means no saves)
  !!   maxSteps :: maximum number of rism relaxation steps
  !!   tolerance    :: convergence tolerance
  subroutine rism3d_calculateSolution_dT(this, kshow, maxSteps, tolerance)
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: kshow, maxSteps
    _REAL_, intent(in) :: tolerance
    ! 1) Check for the existance of xvv_dT data.
    if (.not.rism3d_solvent_canCalc_dT(this%solvent)) then
       call rism_report_warn("cannot perform temperature derivative. No 1D-RISM temperature derivative data found")
       return
    end if


    ! 2) Check for the existance of a 3D-RISM solution.
    ! If memory has been allocated, the solution should already exist.
    if (.not.(associated(this%cuv) .and. associated(this%cuvres))) then
       call rism_report_error("a 3D-RISM solution must be present " &
            //"before a temperature derivative can be calculated")
    end if

    ! 3) Allocate memory.
    call reallocateBox_dT(this)

    ! 4) Calculate 3D-RISM DT solution using MDIIS.
    call timer_start(TIME_RXRISM)
    call  solve3DRISM_dT(this, kshow, maxSteps, tolerance)
    call timer_stop(TIME_RXRISM)

  end subroutine rism3d_calculateSolution_dT

  !> Calculates the forces on the solute contributed by the solvent according
  !! to 3D-RISM.  Just a wrapper for rism3d_closure_force().
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   ff   :: 3D-RISM forces
  subroutine rism3d_force(this, ff)
    implicit none
    type(rism3d):: this
    _REAL_, intent(out) :: ff(3, this%solute%numAtoms)
    call rism_timer_start(this%forceTimer)
    call timer_start(TIME_FF)
    call rism3d_closure_force(this%closure, ff, this%guv, this%periodicPotential)
    call timer_stop(TIME_FF)
    call rism_timer_stop(this%forceTimer)
  end subroutine rism3d_force


  !> Calculate the excess chemical potential for each solvent species
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
  !! OUT:
  !!    excess chemical potential of solvation for each solvent species
  function rism3d_excessChemicalPotential(this, o_lr) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: excessChemicalPotential(this%solvent%numAtomTypes)
    call rism_timer_start(this%excessChemicalPotentialTimer)

    lr = .true.
    if (present(o_lr)) lr = o_lr
    if (lr) then
       excessChemicalPotential = rism3d_closure_aexcessChemicalPotential(this%closure, this%huv, this%cuv(:, :, :, :))
    else
       excessChemicalPotential = rism3d_closure_excessChemicalPotential(this%closure, this%huv, this%cuv(:, :, :, :))
    end if
    call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_excessChemicalPotential


  !> Calculate the total excess chemical potential of solvation
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
  !! OUT:
  !!    total excess chemical potential of solvation
  function rism3d_excessChemicalPotential_tot(this, o_lr) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: excessChemicalPotential
    call rism_timer_start(this%excessChemicalPotentialTimer)

    lr = .true.
    if (present(o_lr)) lr = o_lr
    call rism_timer_stop(this%excessChemicalPotentialTimer)
    excessChemicalPotential = sum(rism3d_excessChemicalPotential(this, o_lr))
  end function rism3d_excessChemicalPotential_tot


  !> Returns a 3D map of the excess chemical potential in kT.
  !! Integrating this map gives the total excess chemical potential as
  !! returned from rism3d_excessChemicalPotential_tot().  Memory is allocated into a
  !! pointer and must be freed by the calling function. For MPI, only
  !! the grid points local to this process are allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the excess chemical potential contributions
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_excessChemicalPotential_tot_map(this) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    logical :: lr
    _REAL_, pointer :: excessChemicalPotential(:, :, :)
    call rism_timer_start(this%excessChemicalPotentialTimer)
    excessChemicalPotential => rism3d_closure_excessChemicalPotential_tot_map(this%closure, this%huv, this%cuv)
    call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_excessChemicalPotential_tot_map


  !> Returns a 3D map of the excess chemical potential in kT.
  !! Integrating this map gives the site excess chemical potential
  !! contribution as returned from rism3d_excessChemicalPotential_tot().  Memory is
  !! allocated into a pointer and must be freed by the calling
  !! function. For MPI, only the grid points local to this process are
  !! allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the excess chemical potential contributions (nx, ny, nz, nsite)
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_excessChemicalPotential_site_map(this) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    logical :: lr
    _REAL_, pointer :: excessChemicalPotential(:, :, :, :)
    call rism_timer_start(this%excessChemicalPotentialTimer)
    excessChemicalPotential => rism3d_closure_excessChemicalPotential_site_map(this%closure, this%huv, this%cuv)
    call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_excessChemicalPotential_site_map


  !> Returns a 3D map of the excess chemical potential in kT from the GF
  !! approximation.  Integrating this map gives the total excess
  !! chemical potential as returned from rism3d_excessChemicalPotential_tot().  Memory is
  !! allocated into a pointer and must be freed by the calling
  !! function. For MPI, only the grid points local to this process are
  !! allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the excess chemical potential contributions
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_excessChemicalPotentialGF_tot_map(this) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    logical :: lr
    _REAL_, pointer :: excessChemicalPotential(:, :, :)
    call rism_timer_start(this%excessChemicalPotentialTimer)
    excessChemicalPotential => rism3d_closure_excessChemicalPotentialGF_tot_map(&
         this%closure, this%huv, this%cuv)
    call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_excessChemicalPotentialGF_tot_map


  !> Returns a 3D map of the excess chemical potential in kT using the
  !! GF approximation.  Integrating this map gives the site excess
  !! chemical potential contribution as returned from
  !! rism3d_excessChemicalPotential_tot().  Memory is allocated into a pointer and must
  !! be freed by the calling function. For MPI, only the grid points
  !! local to this process are allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the excess chemical potential contributions (nx, ny, nz, nsite)
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_excessChemicalPotentialGF_site_map(this) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    logical :: lr
    _REAL_, pointer :: excessChemicalPotential(:, :, :, :)
    call rism_timer_start(this%excessChemicalPotentialTimer)
    excessChemicalPotential => rism3d_closure_excessChemicalPotentialGF_site_map(this%closure, this%huv, this%cuv)
    call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_excessChemicalPotentialGF_site_map

  !> Returns a 3D map of the excess chemical potential in kT from the
  !! PCPLUS.  Integrating this map gives the total excess chemical
  !! potential as returned from rism3d_excessChemicalPotentialPCPLUS_tot().
  !! Memory is allocated into a pointer and must be freed by the
  !! calling function. For MPI, only the grid points local to this
  !! process are allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the excess chemical potential contributions
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_excessChemicalPotentialPCPLUS_tot_map(this)&
       result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, pointer :: excessChemicalPotential(:, :, :)
    call rism_timer_start(this%excessChemicalPotentialTimer)
    excessChemicalPotential => rism3d_closure_excessChemicalPotentialPCPLUS_tot_map(&
         this%closure, this%huv, this%cuv)
    call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_excessChemicalPotentialPCPLUS_tot_map


  !> Returns a 3D map of the excess chemical potential in kT from the
  !! Universal Correction correction.  Integrating this map gives the
  !! total excess chemical potential as returned from
  !! rism3d_excessChemicalPotential_tot().
  !!
  !! Memory is allocated into a pointer and must be freed by the
  !! calling function. For MPI, only the grid points local to this
  !! process are allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   coeff: coefficients for correction.  For the original correction
  !!          (a, b) = coeff(1:2).
  !!          Extra coefficients are a1 and b1 from Johnson et al. 2016.
  !! OUT:
  !!   a 3D-grid of the excess chemical potential contributions
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_excessChemicalPotentialUC_tot_map(this, coeff) &
       result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    logical :: lr
    _REAL_, intent(in) :: coeff(:)
    _REAL_, pointer :: excessChemicalPotential(:, :, :)
    call rism_timer_start(this%excessChemicalPotentialTimer)
    
    excessChemicalPotential => rism3d_closure_excessChemicalPotentialUC_tot_map(&
         this%closure, this%huv, this%cuv, UC_temperature_coeff(this,coeff))
    call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_excessChemicalPotentialUC_tot_map


  !> Calculate the excess chemical potential of solvation for each solvent species
  !! with the Gaussian fluctuation correction
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
  !! OUT:
  !!    Gaussian fluctuation excess chemical potential of solvation for each
  !!    solvent species
  function rism3d_excessChemicalPotentialGF(this, o_lr) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: excessChemicalPotential(this%solvent%numAtomTypes)
    call rism_timer_start(this%excessChemicalPotentialTimer)

    lr = .true.
    if (present(o_lr)) lr = o_lr
    if (lr) then
       excessChemicalPotential = rism3d_closure_aexcessChemicalPotentialGF(this%closure, this%huv, this%cuv(:, :, :, :))
    else
       excessChemicalPotential = rism3d_closure_excessChemicalPotentialGF(this%closure, this%huv, this%cuv(:, :, :, :))
    end if
    call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_excessChemicalPotentialGF


  !> Calculate the total excess chemical potential of solvation with the Gaussian
  !! fluctuation correction
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
  !! OUT:
  !!    total Gaussian fluctuation excess chemical potential of solvation
  function rism3d_excessChemicalPotentialGF_tot(this, o_lr) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: excessChemicalPotential
    call rism_timer_start(this%excessChemicalPotentialTimer)

    lr = .true.
    if (present(o_lr)) lr = o_lr

    call rism_timer_stop(this%excessChemicalPotentialTimer)
    excessChemicalPotential = sum(rism3d_excessChemicalPotentialGF(this, o_lr))
  end function rism3d_excessChemicalPotentialGF_tot


  !>Calculate the total excess chemical potential of solvation with the
  !!PC+/3D-RISM Correction
  !!IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default=.true.) Apply asymptotic long range correction
  !!OUT:
  !!    total Universal Correction excess chemical potential of solvation
  function rism3d_excessChemicalPotentialPCPLUS (this,o_lr) result (excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: excessChemicalPotential

    lr = .true.
    if(present(o_lr)) lr = o_lr

    if(lr)then
       excessChemicalPotential = rism3d_closure_aexcessChemicalPotentialPCPLUS(this%closure,this%huv,&
            this%cuv(:,:,:,:))
    else
       excessChemicalPotential = rism3d_closure_excessChemicalPotentialPCPLUS(this%closure,this%huv,&
            this%cuv(:,:,:,:))
    end if
  end function rism3d_excessChemicalPotentialPCPLUS

  !>Calculate the total solvation energy with the PC+/3D-RISM
  !!Correction.
  !!
  !!IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default=.true.) Apply asymptotic long range correction
  !!OUT:
  !!    total PC+/3D-RISM Correction excess chemical potential of
  !!    solvation
  function rism3d_solvationEnergyPCPLUS (this,o_lr) result (solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: solvationEnergy

    lr = .true.
    if(present(o_lr)) lr = o_lr

    if(lr)then
       solvationEnergy = rism3d_closure_asolvationEnergyPCPLUS(this%closure,&
            this%huv_dT,this%cuv_dT,this%huv,this%cuv)
    else
       solvationEnergy = rism3d_closure_solvationEnergyPCPLUS(this%closure,&
            this%huv_dT,this%cuv_dT,this%huv,this%cuv)
    end if
  end function rism3d_solvationEnergyPCPLUS

  !> Calculate the total excess chemical potential of solvation with the
  !! Palmer et al. Universal Correction with optional temperature dependence.
  !!
  !! Uses the total molecular density of the solvent.  For pure water,
  !! this gives the original Palmer correction.  There is no definition
  !! for mixed solvents but this seems reasonable.
  !!
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   coeff: coefficients for correction.  For the original correction
  !!      (a, b) = coeff(1:2). Extra coefficients are a1, b1 from
  !!      Johson et al. 2016
  !!   o_lr :: (optional) (default = .true.) Apply asymptotic long
  !!       range correction
  !! OUT:
  !!    total Universal Correction excess chemical potential of solvation
  function rism3d_excessChemicalPotentialUC(this, coeff, o_lr) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: coeff(:)
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: excessChemicalPotential

    lr = .true.
    if (present(o_lr)) lr = o_lr

    if (lr) then
       excessChemicalPotential = rism3d_closure_aexcessChemicalPotentialUC(this%closure, this%huv, &
            this%cuv(:, :, :, :), UC_temperature_coeff(this,coeff))
    else
       excessChemicalPotential = rism3d_closure_excessChemicalPotentialUC(this%closure, this%huv, &
            this%cuv(:, :, :, :), UC_temperature_coeff(this,coeff))
    end if
  end function rism3d_excessChemicalPotentialUC


  !> Calculate the total solvation energy with the Palmer et
  !! al. Universal Correction with optional temperature dependence.
  !!
  !! Uses the total molecular density of the solvent.  For pure water,
  !! this gives the original Palmer correction.  There is no
  !! definition for mixed solvents but this seems reasonable.
  !!
  !! @param[in] this rism3d object with computed solution.
  !! @param[in] coeff coefficients for correction.  For the original
  !!     correction (a, b) = coeff(1:2). Extra coefficients are a1, b1
  !!     from Johnson et al. 2016
  !! @param[in] o_lr (optional) (default = .true.)
  !!                 Apply asymptotic long range correction.
  !! @return Total Universal Correction excess chemical potential of solvation.
  function rism3d_solvationEnergyUC(this, coeff, o_lr) result(solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: coeff(:)
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: solvationEnergy

    lr = .true.
    if (present(o_lr)) lr = o_lr

    if (lr) then
       solvationEnergy = rism3d_closure_asolvationEnergyUC(this%closure, &
            this%huv_dT, this%cuv_dT, this%huv, this%cuv, UC_temperature_coeff(this,coeff))
    else
       solvationEnergy = rism3d_closure_solvationEnergyUC(this%closure, &
            this%huv_dT, this%cuv_dT, this%huv, this%cuv, UC_temperature_coeff(this,coeff))
    end if
  end function rism3d_solvationEnergyUC


  !> Calculate the solvation interaction energy: de = density sum g*u for
  !! each solvent site.  I.e., the direct intection potential energy of
  !! solute and solvent and not the total solvation energy (see solvationEnergy).
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!    the contribution of each solvent site to the total solvation interaction
  !!    energy [kT]
  function rism3d_solventPotEne(this) result(ene)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: ene(this%solvent%numAtomTypes)
    call rism_timer_start(this%thermoTimer)
    ene = rism3d_closure_solvPotEne(this%closure, this%guv)
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_solventPotEne


  !> Calculate the total solvation interaction energy: de = density sum g*u for
  !! each solvent site.  I.e., the direct intection potential energy of
  !! solute and solvent and not the total solvation energy (see solvationEnergy).
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!    the total solvent-solute potential energy [kT]
  function rism3d_solventPotEne_tot(this) result(ene)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: ene
    call rism_timer_start(this%thermoTimer)
    ene = sum(rism3d_solventPotEne(this))
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_solventPotEne_tot


  !> Returns a 3D map of the solvent-solute potential energy in kT.
  !! Integrating this map gives the site solvent-solute potential energy
  !! contribution as returned from rism3d_closure_solvPotEne(, .false.).
  !! Memory is allocated into a pointer and must be freed by the calling
  !! function. For MPI, only the grid points local to this process are
  !! allocated and calculated.
  !! IN:
  !!    this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the solvent-solute energy contributions (nx, ny, nz, nsite)
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_solventPotEne_tot_map(this) result(solvPotEne)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, pointer :: solvPotEne(:, :, :)
    !      call rism_timer_start(this%excessChemicalPotentialTimer)
    solvPotEne => rism3d_closure_solvPotEne_tot_map(this%closure, &
         this%guv)
    !      call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_solventPotEne_tot_map


  !> Returns a 3D map of the solvent-solute potential energy in kT.
  !! Integrating this map gives the site solvent-solute potential energy
  !! contribution as returned from rism3d_closure_solvPotEne(, .false.).
  !! Memory is allocated into a pointer and must be freed by the calling
  !! function. For MPI, only the grid points local to this process are
  !! allocated and calculated.
  !! IN:
  !!    this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the solvent-solute energy contributions (nx, ny, nz, nsite)
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_solventPotEne_site_map(this) result(solvPotEne)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, pointer :: solvPotEne(:, :, :, :)
    !      call rism_timer_start(this%excessChemicalPotentialTimer)
    solvPotEne => rism3d_closure_solvPotEne_site_map(this%closure, &
         this%guv)
    !      call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_solventPotEne_site_map


  !> Calculate the total solvation energy dE = dmu + dTS which includes both
  !! solute-solvent interaction energy and the change in solvent
  !! self-energy due to rearranging around the solvent.
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
  !! OUT:
  !!    total solvation energy
  function rism3d_solvationEnergy_tot(this, o_lr) result(excessChemicalPotential_dT)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: excessChemicalPotential_dT

    lr = .true.
    if (present(o_lr)) lr = o_lr
    excessChemicalPotential_dT = sum(rism3d_solvationEnergy(this, lr))
  end function rism3d_solvationEnergy_tot


  !> Returns a 3D map of the solvation energy in kT.
  !! Integrating this map gives the total solvation energy as
  !! returned from rism3d_excessChemicalPotential_tot(, .false.).  Memory is allocated into a
  !! pointer and must be freed by the calling function. For MPI, only
  !! the grid points local to this process are allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the solvation energy contributions
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_solvationEnergy_tot_map(this) result(solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, pointer :: solvationEnergy(:, :, :)
    !      call rism_timer_start(this%excessChemicalPotentialTimer)
    solvationEnergy => rism3d_closure_solvationEnergy_tot_map(this%closure, &
         this%huv_dT, this%cuv_dT, this%huv, this%cuv)
    !      call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_solvationEnergy_tot_map


  !> Returns a 3D map of the solvation energy in kT.
  !! Integrating this map gives the site solvation energy contribution as
  !! returned from rism3d_excessChemicalPotential_tot(, .false.).  Memory is allocated into a
  !! pointer and must be freed by the calling function. For MPI, only
  !! the grid points local to this process are allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the solvation energy contributions (nx, ny, nz, nsite)
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_solvationEnergy_site_map(this) result(solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, pointer :: solvationEnergy(:, :, :, :)
    !      call rism_timer_start(this%excessChemicalPotentialTimer)
    solvationEnergy => rism3d_closure_solvationEnergy_site_map(this%closure, &
         this%huv_dT, this%cuv_dT, this%huv, this%cuv)
    !      call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_solvationEnergy_site_map


  !> Returns a 3D map of the solvation energy in kT from the GF
  !! approximation.  Integrating this map gives the total solvation
  !! energy as returned from rism3d_excessChemicalPotential_tot(, .false.).  Memory is
  !! allocated into a pointer and must be freed by the calling
  !! function. For MPI, only the grid points local to this process are
  !! allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the solvation energy contributions
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_solvationEnergyGF_tot_map(this) result(solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, pointer :: solvationEnergy(:, :, :)
    !      call rism_timer_start(this%excessChemicalPotentialTimer)
    solvationEnergy => rism3d_closure_solvationEnergyGF_tot_map(this%closure, &
         this%huv_dT, this%cuv_dT, this%huv, this%cuv)
    !      call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_solvationEnergyGF_tot_map


  !> Returns a 3D map of the solvation energy in kT from the GF approximation.
  !! Integrating this map gives the site solvation energy contribution as
  !! returned from rism3d_excessChemicalPotential_tot(, .false.).  Memory is allocated into a
  !! pointer and must be freed by the calling function. For MPI, only
  !! the grid points local to this process are allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the solvation energy contributions (nx, ny, nz, nsite)
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_solvationEnergyGF_site_map(this) result(solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, pointer :: solvationEnergy(:, :, :, :)
    !      call rism_timer_start(this%excessChemicalPotentialTimer)
    solvationEnergy => rism3d_closure_solvationEnergyGF_site_map(this%closure, &
         this%huv_dT, this%cuv_dT, this%huv, this%cuv)
    !      call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_solvationEnergyGF_site_map

  !! Returns a 3D map of the solvation energy in kT from the PCPLUS.
  !! Integrating this map gives the total solvation energy as returned
  !! from rism3d_excessChemicalPotential_tot(, .false.).  Memory is
  !! allocated into a pointer and must be freed by the calling
  !! function. For MPI, only the grid points local to this process are
  !! allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the solvation energy contributions
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_solvationEnergyPCPLUS_tot_map(this) result(solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, pointer :: solvationEnergy(:, :, :)
    !      call rism_timer_start(this%excessChemicalPotentialTimer)
    solvationEnergy => rism3d_closure_solvationEnergyPCPLUS_tot_map(this%closure, &
         this%huv_dT, this%cuv_dT, this%huv, this%cuv)
    !      call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_solvationEnergyPCPLUS_tot_map


  !! Returns a 3D map of the solvation energy in kT from the Universal
  !! Correction with optional temperature dependence.
  !!
  !! Integrating this map gives the total solvation energy as returned
  !! from rism3d_excessChemicalPotential_tot(, .false.).  Memory is
  !! allocated into a pointer and must be freed by the calling
  !! function. For MPI, only the grid points local to this process are
  !! allocated and calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   coeff: coefficients for correction.  For the original
  !!       correction (a, b) = coeff(1:2). Extra coefficients are a1,
  !!       b1 from Johnson et al. 2016
  !! OUT:
  !!   a 3D-grid of the solvation energy contributions
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_solvationEnergyUC_tot_map(this, coeff) result(solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: coeff(:)
    _REAL_, pointer :: solvationEnergy(:, :, :)
    !      call rism_timer_start(this%excessChemicalPotentialTimer)
    solvationEnergy => rism3d_closure_solvationEnergyUC_tot_map(this%closure, &
         this%huv_dT, this%cuv_dT, this%huv, this%cuv, UC_temperature_coeff(this,coeff))
    !      call rism_timer_stop(this%excessChemicalPotentialTimer)
  end function rism3d_solvationEnergyUC_tot_map


  !! Calculate the solvation energy dE = dmu + dTS which includes both
  !! solute-solvent interaction energy and the change in solvent
  !! self-energy due to rearranging around the solvent.
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
  !! OUT:
  !!    solvation energy
  function rism3d_solvationEnergy(this, o_lr) result(solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: solvationEnergy(this%solvent%numAtomTypes)

    lr = .true.
    if (present(o_lr)) lr = o_lr

    if (lr) then
       solvationEnergy = rism3d_closure_asolvationEnergy(this%closure, this%huv_dT, this%cuv_dT(:, :, :, :), &
            this%huv, this%cuv(:, :, :, :))
    else
       solvationEnergy = rism3d_closure_solvationEnergy(this%closure, this%huv_dT, this%cuv_dT(:, :, :, :), &
            this%huv, this%cuv(:, :, :, :))
    end if
  end function rism3d_solvationEnergy


  !! Calculate the solvation energy dE = dmu + dTS the Gaussian
  !! fluctuation correction, which includes both solute-solvent
  !! interaction energy and the change in solvent self-energy due to
  !! rearranging around the solvent.
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
  !! OUT:
  !!    solvation energy
  function rism3d_solvationEnergyGF(this, o_lr) result(solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: solvationEnergy(this%solvent%numAtomTypes)

    lr = .true.
    if (present(o_lr)) lr = o_lr

    if (lr) then
       solvationEnergy = rism3d_closure_asolvationEnergyGF(this%closure, this%huv_dT, this%cuv_dT(:, :, :, :), &
            this%huv, this%cuv(:, :, :, :))
    else
       solvationEnergy = rism3d_closure_solvationEnergyGF(this%closure, this%huv_dT, this%cuv_dT(:, :, :, :), &
            this%huv, this%cuv(:, :, :, :))
    end if
  end function rism3d_solvationEnergyGF


  !! Calculate the total solvation energy dE = dmu + dTSthe Gaussian
  !! fluctuation correction, which includes both solute-solvent
  !! interaction energy and the change in solvent self-energy due to
  !! rearranging around the solvent.
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
  !! OUT:
  !!    total solvation energy
  function rism3d_solvationEnergyGF_tot(this, o_lr) result(solvationEnergy)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: solvationEnergy

    lr = .true.
    if (present(o_lr)) lr = o_lr
    solvationEnergy = sum(rism3d_solvationEnergyGF(this, lr))
  end function rism3d_solvationEnergyGF_tot


  !! Calculating the partial molar volume of solute.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   partial molar volume
  function rism3d_partialMolarVolume(this) result(partialMolarVolume)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: partialMolarVolume
    call rism_timer_start(this%thermoTimer)
    partialMolarVolume = rism3d_closure_partialMolarVolume(this%closure, this%cuv(:, :, :, :))
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_partialMolarVolume


  !! Calculating the partial molar volume temperature derivative
  !! (T * d/dT) of solute.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   partial molar volume temperature derivative
  function rism3d_partialMolarVolume_dT(this) result(partialMolarVolume_dT)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: partialMolarVolume_dT
    call rism_timer_start(this%thermoTimer)
    partialMolarVolume_dT = rism3d_closure_partialMolarVolume_dT( &
         this%closure, this%cuv_dT, this%cuv)
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_partialMolarVolume_dT

  !> Returns a 3D map of the partial molar volume temperature derivative.
  !! Integrating this map gives the total PMV as returned from
  !! rism3d_partialMolarVolume_tot().  Memory is allocated into a
  !! pointer and must be freed by the calling function. For MPI, only
  !! the grid points local to this process are allocated and
  !! calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the partial molar volume contributions
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_partialMolarVolume_tot_map(this) result(partialMolarVolume)
    implicit none
    type(rism3d), intent(inout) :: this
    logical :: lr
    _REAL_, pointer :: partialMolarVolume(:, :, :)
    call rism_timer_start(this%thermoTimer)
    partialMolarVolume => rism3d_closure_partialMolarVolume_tot_map(this%closure, this%cuv)
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_partialMolarVolume_tot_map

  !> Returns a 3D map of the partial molar volume .
  !! Integrating this map gives the total PMV as returned from
  !! rism3d_partialMolarVolume_tot().  Memory is allocated into a
  !! pointer and must be freed by the calling function. For MPI, only
  !! the grid points local to this process are allocated and
  !! calculated.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   a 3D-grid of the partial molar volume contributions
  !! SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_partialMolarVolume_dT_tot_map(this) result(partialMolarVolume_dT)
    implicit none
    type(rism3d), intent(inout) :: this
    logical :: lr
    _REAL_, pointer :: partialMolarVolume_dT(:, :, :)
    call rism_timer_start(this%thermoTimer)
    partialMolarVolume_dT => rism3d_closure_partialMolarVolume_dT_tot_map(&
         this%closure, this%cuv_dT, this%cuv)
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_partialMolarVolume_dT_tot_map

  !> Calculating excess number of each solvent type associated with
  !! the solute.
  !! @param[in,out] this rism3d object with computed solution.
  !! @param[in] o_lr (optional) (default = .true.)
  !!                 Apply asymptotic long range correction.
  !! @return Excess number of each solvent type associated with the solute.
  function rism3d_excessParticles(this, o_lr) result(num)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: num(this%solvent%numAtomTypes)
    call rism_timer_start(this%thermoTimer)

    lr = .true.
    if (present(o_lr)) lr = o_lr

    if (lr) then
       num = rism3d_closure_aexcessParticles(this%closure, this%guv)
    else
       num = rism3d_closure_excessParticles(this%closure, this%guv)
    end if
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_excessParticles


  !> Calculating temperature derivative (T * d/dT) of the excess number of each
  !! solvent type associated with the solute.
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
  !! OUT:
  !!    excess number temperature derivative of each solvent type
  !!    associated with the solute
  function rism3d_excessParticles_dT(this, o_lr) result(num_dT)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: num_dT(this%solvent%numAtomTypes)
    call rism_timer_start(this%thermoTimer)

    lr = .true.
    if (present(o_lr)) lr = o_lr

    if (lr) then
       num_dT = rism3d_closure_aexcessParticles_dT(this%closure, this%guv_dT)
    else
       num_dT = rism3d_closure_excessParticles_dT(this%closure, this%guv_dT)
    end if
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_excessParticles_dT


  !! Calculate the Kirkwood-Buff integral for the solute. This is the
  !! all space integral of huv.
  !!
  !! J. G. Kirkwood; F. P. Buff. J. Chem. Phys. 1951, 19, 774-777
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    o_lr :: (optional) (default = .true.) Apply asymptotic long range
  !!            correction
  !! OUT:
  !!    Kirkwood-Buff integeral for each solvent site
  function rism3d_kirkwoodBuff(this, o_lr) result(kb)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: kb(this%solvent%numAtomTypes)
    call rism_timer_start(this%thermoTimer)

    lr = .true.
    if (present(o_lr)) lr = o_lr

    if (lr) then
       kb = rism3d_closure_akirkwoodbuff(this%closure, this%guv)
    else
       kb = rism3d_closure_kirkwoodBuff(this%closure, this%guv)
    end if
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_kirkwoodBuff


  !! Calculate the Kirkwood-Buff temperature derivative (T*d/dT) integral for the
  !! solute. This is the all space integral of huv_dT.
  !!
  !! J. G. Kirkwood; F. P. Buff. J. Chem. Phys. 1951, 19, 774-777
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    o_lr :: (optional) (default = .true.) Apply asymptotic long range
  !!            correction
  !! OUT:
  !!    Kirkwood-Buff integeral temperature derivative for each solvent site
  function rism3d_kirkwoodBuff_dT(this, o_lr) result(kb_dT)
    implicit none
    type(rism3d), intent(inout) :: this
    logical, optional, intent(in) :: o_lr
    logical :: lr
    _REAL_ :: kb_dT(this%solvent%numAtomTypes)
    call rism_timer_start(this%thermoTimer)

    lr = .true.
    if (present(o_lr)) lr = o_lr

    if (lr) then
       kb_dT = rism3d_closure_akirkwoodbuff_dT(this%closure, this%guv_dT)
    else
       kb_dT = rism3d_closure_kirkwoodBuff_dT(this%closure, this%guv_dT)
    end if
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_kirkwoodBuff_dT


  !! Calculates the direct correlation function integral for the solute. This is the
  !! all space integral of cuv.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!    DCF integeral for each solvent site
  function rism3d_DCFintegral(this) result(dcfi)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: dcfi(this%solvent%numAtomTypes)
    call rism_timer_start(this%thermoTimer)
    dcfi = rism3d_closure_DCFintegral(this%closure, this%cuv)
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_DCFintegral


  !! Calculates the direct correlation function integral temperature
  !! derivative (T*d/dT) for the solute. This is the all space integral of cuv.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!    DCF integeral temperature derivative for each solvent site
  function rism3d_DCFintegral_dT(this) result(dcfi_dT)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: dcfi_dT(this%solvent%numAtomTypes)
    call rism_timer_start(this%thermoTimer)
    dcfi_dT = rism3d_closure_DCFintegral_dT(this%closure, this%cuv_dT)
    call rism_timer_stop(this%thermoTimer)
  end function rism3d_DCFintegral_dT


  !! Perform some internal tests for correctness. Some tests require a
  !! converged standard [and temperature derivative] solution. A
  !! warning message will report when tests cannot be run.
  subroutine rism3d_selftest(this)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, pointer::workTot(:, :, :) => NULL()
    _REAL_, pointer::workSite(:, :, :, :) => NULL()
    _REAL_ :: temp(this%solvent%numAtomTypes)
    _REAL_ :: errorTolerance=1d-8
    _REAL_ :: uccoeff(4)
    integer :: iv
    integer :: totalTests, passedTests

    if (this%nsolution<1) then
       call rism_report_warn("SELFTEST cannot run without a RISM solution")
       return
    end if

    ! initialize counters
    totalTests=0
    passedTests=0

    ! UC coefficients are arbitrary.
    uccoeff=1d0
    ! test grid output against internal integral

    ! not tested: guv, uuv, asymp, entropy, entropyGF, entropyPCPLUS, entropyUC
    ! should be tested: cuv, huv, guv, chgdist, potUV

    !! *** Density grids *** 

    !! *** Standard Thermodynamics ***

    !! PMV
    ! total
    if (safemem_dealloc(workTot)/= 0) &
         call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
    workTot => rism3d_partialMolarVolume_tot_map(this)
    temp(1) = rism3d_partialMolarVolume(this)
    totalTests = totalTests + 1
    if (rism3d_integrateCompare(this,workTot,temp(1),errorTolerance,&
         "total partialMolarVolume grid ")) then
       passedTests = passedTests + 1
    end if

    !! *** exchem ***
    ! total
    if (safemem_dealloc(workTot)/= 0) &
         call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
    workTot => rism3d_excessChemicalPotential_tot_map(this)
    temp(1) = rism3d_excessChemicalPotential_tot(this,.false.)
    totalTests = totalTests + 1
    if (rism3d_integrateCompare(this,workTot,temp(1),errorTolerance,&
         "total excessChemicalPotential grid ")) then
       passedTests = passedTests + 1
    end if
    ! sites
    if (safemem_dealloc(workSite)/= 0) &
         call rism_report_error("RISM3D_SELFTEST: failed to deallocate workSite")
    workSite => rism3d_excessChemicalPotential_site_map(this)
    temp = rism3d_excessChemicalPotential(this,.false.)
    do iv=1, this%solvent%numAtomTypes
       totalTests = totalTests + 1
       if (rism3d_integrateCompare(this,workSite(:,:,:,iv),temp(iv),errorTolerance,&
            "excessChemicalPotential grid for site "//this%solvent%atomName(iv))) then
          passedTests = passedTests + 1
       end if
    end do

    !! *** exchemGF ***
    ! total
    if (safemem_dealloc(workTot)/= 0) &
         call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
    workTot => rism3d_excessChemicalPotentialGF_tot_map(this)
    temp(1) = rism3d_excessChemicalPotentialGF_tot(this,.false.)
    totalTests = totalTests + 1
    if (rism3d_integrateCompare(this,workTot,temp(1),errorTolerance,&
         "total excessChemicalPotentialGF grid ")) then
       passedTests = passedTests + 1
    end if
    ! sites
    if (safemem_dealloc(workSite)/= 0) &
         call rism_report_error("RISM3D_SELFTEST: failed to deallocate workSite")
    workSite => rism3d_excessChemicalPotentialGF_site_map(this)
    temp = rism3d_excessChemicalPotentialGF(this,.false.)
    do iv=1, this%solvent%numAtomTypes
       totalTests = totalTests + 1
       if (rism3d_integrateCompare(this,workSite(:,:,:,iv),temp(iv),errorTolerance,&
            "excessChemicalPotentialGF grid for site "//this%solvent%atomName(iv))) then
          passedTests = passedTests + 1
       end if
    end do

    !! *** exchemPCPLUS ***
    ! total
    if (safemem_dealloc(workTot)/= 0) &
         call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
    workTot => rism3d_excessChemicalPotentialPCPLUS_tot_map(this)
    temp(1) = rism3d_excessChemicalPotentialPCPLUS(this,.false.)
    totalTests = totalTests + 1
    if (rism3d_integrateCompare(this,workTot,temp(1),errorTolerance,&
         "total excessChemicalPotentialPCPLUS grid ")) then
       passedTests = passedTests + 1
    end if

    !! *** exchemUC ***
    ! total
    if (safemem_dealloc(workTot)/= 0) &
         call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
    workTot => rism3d_excessChemicalPotentialUC_tot_map(this,uccoeff)
    temp(1) = rism3d_excessChemicalPotentialUC(this,uccoeff,.false.)
    totalTests = totalTests + 1
    if (rism3d_integrateCompare(this,workTot,temp(1),errorTolerance,&
         "total excessChemicalPotentialUC grid ")) then
       passedTests = passedTests + 1
    end if

    !! *** entropic decomp ***

    if (associated(this%cuv_dT)) then
       !! PMV_dT
       ! total
       if (safemem_dealloc(workTot)/= 0) &
            call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
       workTot => rism3d_partialMolarVolume_dT_tot_map(this)
       temp(1) = rism3d_partialMolarVolume_dT(this)
       totalTests = totalTests + 1
       if (rism3d_integrateCompare(this,workTot,temp(1),errorTolerance,&
            "total partialMolarVolume_dT grid ")) then
          passedTests = passedTests + 1
       end if

       !! *** solvene ***
       ! total
       if (safemem_dealloc(workTot)/= 0) &
            call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
       workTot => rism3d_solvationEnergy_tot_map(this)
       temp(1) = rism3d_solvationEnergy_tot(this,.false.)
       totalTests = totalTests + 1
       if (rism3d_integrateCompare(this,workTot,temp(1),errorTolerance,&
            "total solvationEnergy grid ")) then
          passedTests = passedTests + 1
       end if
       ! sites
       if (safemem_dealloc(workSite)/= 0) &
            call rism_report_error("RISM3D_SELFTEST: failed to deallocate workSite")
       workSite => rism3d_solvationEnergy_site_map(this)
       temp = rism3d_solvationEnergy(this,.false.)
       do iv=1, this%solvent%numAtomTypes
          totalTests = totalTests + 1
          if (rism3d_integrateCompare(this,workSite(:,:,:,iv),temp(iv),errorTolerance,&
               "solvationEnergy grid for site "//this%solvent%atomName(iv))) then
             passedTests = passedTests + 1
          end if
       end do

       !! *** solveneGF ***
       ! total
       if (safemem_dealloc(workTot)/= 0) &
            call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
       workTot => rism3d_solvationEnergyGF_tot_map(this)
       temp(1) = rism3d_solvationEnergyGF_tot(this,.false.)
       totalTests = totalTests + 1
       if (rism3d_integrateCompare(this,workTot,temp(1),errorTolerance,&
            "total solvationEnergyGF grid ")) then
          passedTests = passedTests + 1
       end if
       ! sites
       if (safemem_dealloc(workSite)/= 0) &
            call rism_report_error("RISM3D_SELFTEST: failed to deallocate workSite")
       workSite => rism3d_solvationEnergyGF_site_map(this)
       temp = rism3d_solvationEnergyGF(this,.false.)
       do iv=1, this%solvent%numAtomTypes
          totalTests = totalTests + 1
          if (rism3d_integrateCompare(this,workSite(:,:,:,iv),temp(iv),errorTolerance,&
               "solvationEnergyGF grid for site "//this%solvent%atomName(iv))) then
             passedTests = passedTests + 1
          end if
       end do


       !! *** solvenePCPLUS ***
       ! total
       if (safemem_dealloc(workTot)/= 0) &
            call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
       workTot => rism3d_solvationEnergyPCPLUS_tot_map(this)
       temp(1) = rism3d_solvationEnergyPCPLUS(this,.false.)
       totalTests = totalTests + 1
       if (rism3d_integrateCompare(this,workTot,temp(1),errorTolerance,&
            "total solvationEnergyPCPLUS grid ")) then
          passedTests = passedTests + 1
       end if

       !! *** solveneUC ***
       ! total
       if (safemem_dealloc(workTot)/= 0) &
            call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
       workTot => rism3d_solvationEnergyUC_tot_map(this,uccoeff)
       temp(1) = rism3d_solvationEnergyUC(this,uccoeff,.false.)
       totalTests = totalTests + 1
       if (rism3d_integrateCompare(this,workTot,temp(1),errorTolerance,&
            "total solvationEnergyUC grid ")) then
          passedTests = passedTests + 1
       end if
    else
       call rism_report_warn("SELFTEST not performing temperature dependent checks.")
       return
    end if

    call rism_report_message("(a,i3,a,i3,a)", "3D-RISM selftest: ",passedTests,&
         " of ", totalTests, " passed.") 
    if (safemem_dealloc(workTot)/= 0) &
         call rism_report_error("RISM3D_SELFTEST: failed to deallocate workTot")
    if (safemem_dealloc(workSite)/= 0) &
         call rism_report_error("RISM3D_SELFTEST: failed to deallocate workSite")
  end subroutine rism3d_selftest

  !! Write to log file whether or not the test passes within the
  !! prescribed error. This is a common utility function for
  !! selftest().
  !!
  !! @param[in] this           rism3d object
  !! @param[in] grid           array to integrate.  Assumes gridspacing of
  !!                           this%grid%voxelVolume
  !! @param[in] reference      pre-integrated reference value
  !! @param[in] errorTolerance allowable error
  !! @param[in] description    description of the tested quantity
  !!
  !! @returns .true. if passed.
  function rism3d_integrateCompare(this,grid, reference, errorTolerance, description)&
       result (passed)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    type(rism3d) :: this
    _REAL_, intent(in) :: grid(:,:,:)
    _REAL_, intent(in) :: reference
    _REAL_, intent(in) :: errorTolerance
    character(len=*), intent(in) :: description
    logical :: passed
    _REAL_ :: gridsum, difference
    _REAL_ :: reduce_gridsum, reduce_reference
    integer :: err
    gridsum = sum(grid)
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
    reduce_gridsum = gridsum
    reduce_reference = reference
    if (this%mpirank == 0) then
       call mpi_reduce(MPI_IN_PLACE, reduce_gridsum, &
            1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, this%mpicomm, err)
       if (err/= 0) call rism_report_warn("RISM3D_SELFTEST: MPI_REDUCE failed.")
       call mpi_reduce(MPI_IN_PLACE, reduce_reference, &
            1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, this%mpicomm, err)
       if (err/= 0) call rism_report_warn("RISM3D_SELFTEST: MPI_REDUCE failed.")
    else
       call MPI_REDUCE(reduce_gridsum, reduce_gridsum, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0,this%mpicomm, err)
       if (err/= 0) call rism_report_warn("RISM3D_SELFTEST: MPI_REDUCE failed.")
       call MPI_REDUCE(reduce_reference, reduce_reference, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0,this%mpicomm, err)
       if (err/= 0) call rism_report_warn("RISM3D_SELFTEST: MPI_REDUCE failed.")
    end if
#  else
    call MPI_REDUCE(gridsum, reduce_gridsum, 1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, 0,this%mpicomm, err)
    if (err/= 0) call rism_report_warn("RISM3D_SELFTEST: MPI_REDUCE failed.")
    call MPI_REDUCE(reference, reduce_reference, 1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, 0,this%mpicomm, err)
    if (err/= 0) call rism_report_warn("RISM3D_SELFTEST: MPI_REDUCE failed.")
#  endif /*USE_MPI_IN_PLACE*/    
#else
    reduce_gridsum = gridsum
    reduce_reference = reference
#endif /*MPI*/    
    if (this%grid%mpirank == 0) then
       difference = abs(reduce_gridsum*this%grid%voxelVolume - reduce_reference)
       if ( difference < errorTolerance) then
          call rism_report_message("(a,e8.2,a,e8.2)","|SELFTEST of "//&
               description//" PASSED with a difference of ", difference,&
               " < ", errorTolerance)
          passed = .true.
       else
          call rism_report_warn("(a,e8.2,a,e8.2)","|SELFTEST of "//&
               description//" FAILED with a difference of ", difference,&
               " < ", errorTolerance)
          passed = .false.
       end if
       call rism_report_flush()
    end if
  end function rism3d_integrateCompare

!!!!!
  !! DEALLOCATE
!!!!!

!!!!!
  !! deconstructor - frees all memory
!!!!!
  subroutine rism3d_destroy(this)
    use safemem
    implicit none
    type(rism3d) :: this
    call rism_timer_destroy(this%fftTimer)
    call rism_timer_destroy(this%single3DRISMsolutionTimer)
    call rism_timer_destroy(this%solve3DRISMTimer)
    call rism_timer_destroy(this%fft_dTTimer)
    call rism_timer_destroy(this%single3DRISM_dTsolutionTimer)
    call rism_timer_destroy(this%solve3DRISM_dTTimer)
    call rism_timer_destroy(this%cuvpropTimer)
    call rism_timer_destroy(this%reorientTimer)
    call rism_timer_destroy(this%resizeTimer)
    call rism_timer_destroy(this%solventTimer)
    call rism_timer_destroy(this%excessChemicalPotentialTimer)
    call rism_timer_destroy(this%forceTimer)
    call rism_timer_destroy(this%thermoTimer)
    call rism_timer_destroy(this%timer)

    call rism3d_solvent_destroy(this%solvent)
    call rism3d_solute_destroy(this%solute)
    call rism3d_potential_destroy(this%potential)
    call rism3d_grid_destroy(this%grid)
    call rism3d_closure_destroy(this%closure)
    call mdiis_destroy(this%mdiis_o)

    if (safemem_dealloc(this%xvva) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate XVVA")
    if (safemem_dealloc(this%cuvWRK) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate CUVWRK")
    if (safemem_dealloc(this%oldcuvChg) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate OLDCUVCHG")
    if (safemem_dealloc(this%oldcuvNoChg) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate OLDCUVNOCHG")
    if (safemem_dealloc(this%cuvresWRK) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate CUVRESWRK")
    if (safemem_dealloc(this%xvva_dT) /= 0) call rism_report_error("RISM3D: failed to deallocate XVVDTA")
    if (safemem_dealloc(this%cuvWRK_dT) /= 0) call rism_report_error("RISM3D: failed to deallocate CUVDTWRK")
    if (safemem_dealloc(this%cuvresWRK_dT) /= 0) call rism_report_error("RISM3D: failed to deallocate CUVRESDTWRK")

    if (safemem_dealloc(this%closureList) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate CLOSURELIST")
    if (associated(this%nsolutionChg)) &
         deallocate(this%nsolutionChg)
    if (associated(this%nsolutionNoChg)) &
         deallocate(this%nsolutionNoChg)
    nullify(this%cuv)
    nullify(this%cuvres)
    nullify(this%oldcuv)
    nullify(this%nsolution)

    if (safemem_dealloc(this%guv, o_aligned = .true.) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate GUV")
    if (safemem_dealloc(this%huv, o_aligned = .true.) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate HUV")
    if (safemem_dealloc(this%cuvk, o_aligned = .true.) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate CUVK")
    if (safemem_dealloc(this%guv_dT, o_aligned = .true.) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate GUVDT")
    if (safemem_dealloc(this%huv_dT, o_aligned = .true.) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate HUVDT")
    if (safemem_dealloc(this%electronMap) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate electronMap")
    call rism3d_fft_destroy(this%fft)
    call rism3d_fft_destroy(this%fft_dT)
    call rism3d_fft_destroy(this%fft_cuv)
    call rism3d_fft_global_finalize()
  end subroutine rism3d_destroy


  !! PRIVATE



  !! Using the current orientation of the solute, define the minimum box size
  !! and resize all associated grids.
  !! IN:
  !!   this :: rism3d object
  subroutine resizeBox(this)
    use constants, only : PI
    use rism_util, only : isprime, lcm, isFactorable, largestPrimeFactor
    implicit none
#if defined(MPI)
    include 'mpif.h'
    integer :: ierr
#endif /*defined(MPI)*/
    type(rism3d), intent(inout) :: this
    integer :: ngr(3)
    integer :: id
    _REAL_ :: boxlen(3)
    integer :: primes(4) = (/2, 3, 5, 7/)
    logical :: cuv_dimension_changed
    integer :: ngr_size

    ! To properly distribute the run, y and z dimension must have a
    ! number of grid points that is a multiple of mpisize.  Thus,
    ! mpisize must be factorizable by small prime numbers.
    if (.not. isFactorable(this%mpisize, primes)) then
       call rism_report_error("(a,10i4)", &
            "Sorry, 3D-RISM requires that the number " &
            //"of processes be a product of ", primes)
    end if

    ! If we have a fixed box size, we simply retain all of the
    ! previously calculated box size values.

    if (this%periodic) then
       !TODO: For periodic case, specifying buffer or box dimensions
       ! is useless, so mention this in user documentation and print a
       ! warning if user attempts to combine periodicity with either.
       ! Also, automatically enable periodicity whenever a cell / box
       ! size is present.

       !TODO: Support both grid spacing or grid dimensions.

       boxlen(:) = this%grid%unitCellLengths(:)

       ! Ensure gridpoints fit in unit cell perfectly by adjusting
       ! grid spacing. Current approach treats user-specified grid
       ! spacing as a maximum spacing.

       ngr(:) = ceiling(boxlen(:) / this%grid%spacing(:))
       this%grid%spacing(:) = boxlen(:) / ngr(:)

       ! Determine if the number of grid points in each dimension
       ! product only has prime factors 2, 3, 5 or 7. If not
       ! increment the number of points (in that dimension) until
       ! this is true.

       ! Make sure that each dimension is divisible by 2 and that the
       ! y- and z-dimensions are divisible by this%mpisize if
       ! this%mpisize > 2.  This former is done for the sake of FFT
       ! libraries.
       ! dac note: this code does not ensure that the y and z dimesions
       !   will still be even after dividing things up among the MPI
       !   threads....
       ! dac note: why do both y and z need to be multiples of mpisize?
       !   is in not just z?

       ngr(:) = ngr(:) + mod(ngr(:), 2)
       if (this%mpisize > 2) then
          do id = 2, 3
             if (mod(ngr(id), this%mpisize) /= 0) then
                ngr(id) = ngr(id) + lcm(this%mpisize, 2) &
                     - mod(ngr(id), lcm(this%mpisize, 2))
             end if
          end do
       end if

       do id = 1, 3
          do while (largestPrimeFactor(ngr(id)) .gt. 9)
             if (this%mpisize > 1 .and. id > 1) then
                ngr(id) = ngr(id) + lcm(this%mpisize, 2)
             else
                ngr(id) = ngr(id) + 2
             end if
          end do
       end do
       this%grid%spacing = boxlen / ngr

#if 0
       if( this%mpirank == 0 ) then
          write(6,*) 'in resizeBox:'
          !write(6,*) 'angles:', this%grid%unitCellAngles, 'box calc:', &
          !     this%grid%unitCellLengths(3) * cos(this%grid%unitCellAngles(2))
          write(6,*) 'box length:', boxlen
          write(6,*) 'grid points:', ngr
          write(6,*) 'grid spacing:', this%grid%spacing
       endif
#endif

       cuv_dimension_changed = .false.
       do id = 1, 2, 3
          ngr_size = merge(ngr(id), ngr(id) / this%mpisize, id /= 3)
          if (ubound(this%cuv, id) /= ngr_size &
               .or. ubound(this%oldcuv, id) /= ngr_size) then
             cuv_dimension_changed = .true.
             exit
          end if
       end do
       if (.not. associated(this%cuv) &
            .or. .not. associated(this%oldcuv) &
            .or. cuv_dimension_changed) then
          call reallocateBox(this, ngr, this%grid%spacing)
       end if

       if (this%potential%cutoff >= this%grid%inscribedSphereRadius) then
          call rism_report_error('(a,f8.2,a)', 'solvcut must be < ', &
               this%grid%inscribedSphereRadius, &
               ', the largest inscribed sphere radius of the unit cell.')
       end if

    else if (this%varbox) then

       ! Get minimum box size defined by the buffer.
       do id = 1, 3
          boxlen(id) = maxval(this%solute%position(id, :)) &
               - minval(this%solute%position(id, :)) + 2 * this%buffer
       end do

       ! Round this box size using the prescribed linear density and
       ! get the number of grid points required.
       ngr = ceiling(boxlen / this%grid%spacing)
       boxlen = ngr * this%grid%spacing

       ! Determine if the number of grid points in each dimension
       ! product only has prime factors 2, 3, 5 or 7. If not
       ! increment the number of points (in that dimension) until
       ! this is true.

       ! Make sure that each dimension is divisible by 2 and that
       ! the y- and z-dimensions are divisible by this%mpisize if
       ! this%mpisize > 1.

       do id = 1, 3
          ngr(id) = ngr(id) + mod(ngr(id), 2)
       end do
       if (this%mpisize > 2) then
          do id = 2, 3
             if (mod(ngr(id), this%mpisize) /= 0) then
                ngr(id) = ngr(id) + lcm(this%mpisize, 2) &
                     - mod(ngr(id), lcm(this%mpisize, 2))
             end if
          end do
       end if

       do id = 1, 3
          do while (.not. isfactorable(ngr(id), primes))
             if (this%mpisize > 1 .and. id > 1) then
                ngr(id) = ngr(id) &
                     + lcm(this%mpisize,2)
             else
                ngr(id) = ngr(id) + 2
             end if
          end do
       end do
       boxlen = ngr * this%grid%spacing

       cuv_dimension_changed = .false.
       do id = 1, 2, 3
          ngr_size = merge(ngr(id), ngr(id) / this%mpisize, id /= 3)
          if (ubound(this%cuv, id) /= ngr_size &
               .or. ubound(this%oldcuv, id) /= ngr_size) then
             cuv_dimension_changed = .true.
             exit
          end if
       end do
       if (.not. associated(this%cuv) &
            .or. .not. associated(this%oldcuv) &
            .or. cuv_dimension_changed) then
          call reallocateBox(this, ngr, this%grid%spacing)
       end if

    else ! Fixed box size.

       do id = 2, 3
          if (this%mpirank == 0 .and. (mod(this%fixedNumGridPoints(id), this%mpisize) /= 0 &
               .or. mod(this%fixedNumGridPoints(id), 2) /= 0) ) then
             call rism_report_error( &
                  "Sorry, MPI 3D-RISM requires that fixed grid sizes be "// &
                  "divisible by two and the number of processes in the "// &
                  "y and z dimensions")
          end if
       end do
       call reallocateBox(this, this%fixedNumGridPoints, this%fixedBoxDimensionsR/this%fixedNumGridPoints)
    end if

  end subroutine resizeBox


  !> Using the current box size and resize all associated grids and
  !! variables.
  !! @param[in,out] this rism3d object.
  !! @param[in] ngr Number of grid points along each axis.
  !! @param[in] grdspc Grid spacing along each axis.
  subroutine reallocateBox(this, ngr, grdspc)
    use constants, only : pi
    use rism3d_fft_c
    use safemem
    implicit none
    type(rism3d) :: this
    integer, intent(in) :: ngr(3)
    _REAL_, intent(in) :: grdspc(3)
    integer :: i, id, irank
    _REAL_ :: memuse
    _REAL_ :: unitCellDimensions(6)

    call rism3d_fft_setgrid(this%grid, ngr, grdspc, this%solvent%numAtomTypes, this%fft_aligned)

    !TODO: MUST CHECK THIS FOR MPI AND NON_MPI CASE

    ! All grid sizes have been determined so report details about the
    ! grid at this point.
    if (this%verbose >= 1) then
       call rism_report_message("||Setting solvation box to")
       call rism_report_message("(3(a,i10))", "|grid size: ", &
            ngr(1), " X ", ngr(2), " X ", ngr(3))
       call rism_report_message("(3(a,f10.3))", "|box size [A]:  ", &
            ngr(1) * grdspc(1), " X ", ngr(2) * grdspc(2), " X ", ngr(3) * grdspc(3))
       call rism_report_message("(3(a,f10.3))", "|grid spacing [A]: ", &
            grdspc(1), " X ", grdspc(2), " X ", grdspc(3))
       if (this%periodic) then
          call rism_report_message("(3(a,f10.3))", "|internal angles [°]:  ", &
               this%grid%unitCellAngles(1) * 180 / pi, ", ", &
               this%grid%unitCellAngles(2) * 180 / pi, ", ", &
               this%grid%unitCellAngles(3) * 180 / pi)
          call rism_report_message("(a,f10.3)", "|inscribed sphere radius [A]: ", &
               this%grid%inscribedSphereRadius)
       else
          call rism_report_message("(3(a,f10.3))", "|effective buffer [A]:", &
               (ngr(1) * grdspc(1) - (maxval(this%solute%position(1,:)) &
               - minval(this%solute%position(1,:)))) / 2d0, ",  ", &
               (ngr(2) * grdspc(2) - (maxval(this%solute%position(2,:)) &
               - minval(this%solute%position(2,:)))) / 2d0, ",  ", &
               (ngr(3) * grdspc(3) - (maxval(this%solute%position(3,:)) &
               - minval(this%solute%position(3,:)))) / 2d0)
       end if

       !! $     memuse = 8d0 * (dble(product(ngr)) * this%solvent%numAtomTypes * (2 * this%NVec + 1+this%ncuvsteps) &
       !! $          +(dble(product(ngr)) + dble(product(ngr(1:2)))) * (4 + 2 * this%solvent%numAtomTypes) & ! see manual. Does not include FFTW scratch
       !! $          +product(ubound(this%solvent%xvv)) * this%mpisize & ! XVV
       !! $          +product(ubound(this%solvent%waveNumbers)) * this%mpisize & ! fourier_table
       !! $          +this%solute%numAtoms * 7*this%mpisize & ! solute
       !! $          +this%grid%waveNumberArraySize * 2*this%mpisize & ! XVVA and grid%waveNumbers
       !! $          +this%grid%totalLocalPointsK/2 * 4*this%mpisize) & ! grid%waveVectors and grid%waveVectors2
       !! $          +4 * (this%grid%totalLocalPointsK/2 * this%mpisize) ! grid%waveVectorWaveNumberMap
       !! $#ifdef MPI
       !! $     memuse = memuse +8d0 * this%grid%nktotal ! FFTW scratch
       !! $#else
       !! $     memuse = memuse +8d0 * product(ngr(2:3)) * 2 ! FFTW scratch
       !! $#endif /*MPI*/
       !! $     call rism_report_message("(a,g12.4,a)","|Projected mimimum distributed memory use: ", &
       !! $          memuse/1024d0**3, " GB" )
       !! $     memuse = memuse + 8d0 * (dble(product(ngr)) * this%solvent%numAtomTypes * this%ncuvsteps &
       !! $          +this%solute%numAtoms * 7*this%mpisize)  ! solute
       !! $
       !! $     call rism_report_message("(a,g12.4,a)","|with polar/apolar decomposition:          ", &
       !! $          memuse/1024d0 * *3," GB" )
       call flush(rism_report_getmunit())
    end if

    !
    ! 2) allocation
    !
    ! reallocate arrays that require preservation of their contents

    ! I THINK WE CAN GET RID OF THE MPI HERE

#if defined(MPI)
    this%cuvWRK => safemem_realloc(this%cuvWRK, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3), &
         this%solvent%numAtomTypes, this%NVec, .true., .true.)
    if (rism3d_solute_charged(this%solute)) then
       this%oldcuvChg => safemem_realloc(this%oldcuvChg, this%grid%localDimsR(1), this%grid%localDimsR(2), &
            this%grid%localDimsR(3), this%solvent%numAtomTypes, this%ncuvsteps, .true., .true.)
       this%oldcuv => this%oldcuvChg
    else
       this%oldcuvNoChg => safemem_realloc(this%oldcuvNoChg, this%grid%localDimsR(1), this%grid%localDimsR(2), &
            this%grid%localDimsR(3), this%solvent%numAtomTypes, this%ncuvsteps, .true., .true.)
       this%oldcuv => this%oldcuvNoChg
    end if
#else
    this%cuvWRK => safemem_realloc(this%cuvWRK, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3), &
         this%solvent%numAtomTypes, this%NVec, .true., .true.)
    if (rism3d_solute_charged(this%solute)) then
       this%oldcuvChg => safemem_realloc(this%oldcuvChg, this%grid%localDimsR(1), this%grid%localDimsR(2), &
            this%grid%localDimsR(3), this%solvent%numAtomTypes, this%ncuvsteps, .true., .true.)
       this%oldcuv => this%oldcuvChg
    else
       this%oldcuvNoChg => safemem_realloc(this%oldcuvNoChg, this%grid%localDimsR(1), this%grid%localDimsR(2), &
            this%grid%localDimsR(3), this%solvent%numAtomTypes, this%ncuvsteps, .true., .true.)
       this%oldcuv => this%oldcuvNoChg
    end if
    !! $  this%oldcuv => safemem_realloc(this%oldcuv, this%grid%localDimsR(1), this%grid%localDimsR(2), &
    !! $       this%grid%localDimsR(3), this%solvent%numAtomTypes, this%ncuvsteps, .true.)
#endif /*MPI*/
    ! reallocate arrays that do not require preservation of their contents
    this%guv => safemem_realloc(this%guv, this%grid%totalLocalPointsK, this%solvent%numAtomTypes, &
         o_preserve = .false., o_aligned = .true.)
    this%huv => safemem_realloc(this%huv, this%grid%totalLocalPointsK, this%solvent%numAtomTypes, &
         o_preserve = .false., o_aligned = .true.)
    this%cuvresWRK => safemem_realloc(this%cuvresWRK, this%grid%totalLocalPointsR, this%solvent%numAtomTypes, &
         this%NVec, .false.)
    this%xvva => safemem_realloc(this%xvva, this%grid%waveNumberArraySize * (this%solvent%numAtomTypes)**2, .false.)

    !  dac: only need three dimensions here, since the electronMap
    !  is written to disk as soon as it is calculated
    this%electronMap => safemem_realloc(this%electronMap, &
         this%grid%globalDimsR(1), this%grid%globalDimsR(2), &
         this%grid%globalDimsR(3), o_preserve = .false.)

    call rism3d_fft_destroy(this%fft)
    call rism3d_fft_new(this%fft, &
         this%fftw_planner, this%fftw_localtrans, this%fft_aligned, &
         this%grid, &
         this%guv, this%huv)

    ! updated pointers
    call mdiis_resize(this%mdiis_o, this%cuvWRK, this%cuvresWRK, &
         this%grid%totalLocalPointsR * this%solvent%numAtomTypes, this%nvec)
    this%cuv => this%cuvWRK(:, :, :, :, mdiis_getWorkVector(this%mdiis_o))
    this%cuvres => this%cuvresWRK(:, :, mdiis_getWorkVector(this%mdiis_o))

    !
    ! 3) the remaining variables are handled by rism_setup_wavevector, interpolateSolventSusceptibility
    !
    call interpolateSolventSusceptibility(this, this%solvent%xvv, this%xvva)
  end subroutine reallocateBox


  !> Using the current box size and resize all the dT associated grids
  !! and variables.
  !! IN:
  !!   this :: rism3d object
  subroutine reallocateBox_dT(this)
    use rism3d_fft_c
    use safemem
    implicit none
    type(rism3d), intent(inout) :: this

    this%cuvWRK_dT => safemem_realloc(this%cuvWRK_dT, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3), &
         this%solvent%numAtomTypes, this%NVec, .true.)
    this%cuvresWRK_dT => safemem_realloc(this%cuvresWRK_dT, this%grid%totalLocalPointsR, &
         this%solvent%numAtomTypes, this%NVec, .false.)
    this%cuv_dT => this%cuvWRK_dT(:, :, :, :, 1)
    this%cuvres_dT => this%cuvresWRK_dT(:, :, 1)

    this%xvva_dT => safemem_realloc(this%xvva_dT, this%grid%waveNumberArraySize* (this%solvent%numAtomTypes) **2, .false.)

    this%cuvk => safemem_realloc(this%cuvk, this%grid%totalLocalPointsK, this%solvent%numAtomTypes, &
         o_preserve = .false., o_aligned = .true.)
    this%guv_dT => safemem_realloc(this%guv_dT, this%grid%totalLocalPointsK, this%solvent%numAtomTypes, &
         o_preserve = .false., o_aligned = .true.)
    this%huv_dT => safemem_realloc(this%huv_dT, this%grid%totalLocalPointsK, this%solvent%numAtomTypes, &
         o_preserve = .false., o_aligned = .true.)

    call rism3d_fft_destroy(this%fft_dT)
    call rism3d_fft_destroy(this%fft_cuv)
    call rism3d_fft_new(this%fft_dT, &
         this%fftw_planner, this%fftw_localtrans, this%fft_aligned, &
         this%grid, &
         this%guv_dT, this%huv_dT)
    call rism3d_fft_new(this%fft_cuv, &
         this%fftw_planner, this%fftw_localtrans, this%fft_aligned, &
         this%grid, &
         this%cuvk, this%cuvk)

    call interpolateSolventSusceptibility(this, this%solvent%xvv_dT, this%xvva_dT)
  end subroutine reallocateBox_dT


  !> Prints the maximum amount of memory allocated at any one time so
  !! far in the run.
  subroutine rism3d_max_memory(this)
    use safemem
    implicit none
#ifdef MPI
    include "mpif.h"
#endif /*MPI*/
    type(rism3d) :: this
    integer * 8 :: memstats(10), tmemstats(10)
    integer :: err, irank, outunit
    outunit = rism_report_getMUnit()
    memstats = memStatus()
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
    if (this%mpirank == 0) then
       call MPI_REDUCE(MPI_IN_PLACE, memstats, ubound(memstats, 1), MPI_INTEGER8, &
            MPI_SUM, 0,this%mpicomm, err)
    else
       call MPI_REDUCE(memstats, memstats, ubound(memstats, 1), MPI_INTEGER8, &
            MPI_SUM, 0,this%mpicomm, err)
    end if
#  else /*USE_MPI_IN_PLACE*/
    call MPI_REDUCE(memstats, tmemstats, ubound(memstats, 1), MPI_INTEGER8, &
         MPI_SUM, 0,this%mpicomm, err)
    memstats = tmemstats
#  endif /*USE_MPI_IN_PLACE*/
    if (err/=0) call rism_report_warn("RISM_MAX_MEMORY: MPI_REDUCE failed.")
#endif
    if (this%mpirank == 0) then
       write(outunit, '(a)')
       write(outunit, '(a)') "|3D-RISM memory allocation summary"
       write(outunit, '(a)') "|Type          Maximum"
       write(outunit, '(a,f12.5,a)') "|Integer  ", &
            dble(memstats(1))/BYTES_PER_GB, " GB"
       write(outunit, '(a,f12.5,a)') "|Real     ", &
            dble(memstats(2))/BYTES_PER_GB, " GB"
       write(outunit, '(a,f12.5,a)') "|Logical  ", &
            dble(memstats(3))/BYTES_PER_GB, " GB"
       write(outunit, '(a,f12.5,a)') "|Character", &
            dble(memstats(4))/BYTES_PER_GB, " GB"
       write(outunit, '(a)') "|------------------------"
       write(outunit, '(a,f12.5,a)') "|Total    ", &
            dble(memstats(5))/BYTES_PER_GB, " GB"
    end if
  end subroutine rism3d_max_memory


  !> Interpolate the solvent-solvent susceptibility, solved on the
  !! 1D-RISM grid, to the 3D-RISM grid.
  !! @param[in,out] this rism3d object.
  !! @param[in] xvv 1D-RISM Xvv or Xvv_dT data.
  !! @param[out] xvva Interpolated result.
  subroutine interpolateSolventSusceptibility(this, xvv, xvva)
    use rism_util, only : polynomialInterpolation
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: xvv(:, :, :)
    _REAL_, intent(out) :: xvva(:)
    integer :: iwn, igk, igk1, iv1, iv2
    _REAL_ :: err
    !> Maximum number of points to interpolate.
    integer :: maxPointsToInterp
    parameter (maxPointsToInterp = 5)

#ifdef RISM_DEBUG
    write(0,*) "INTERPOLATESOLVENTSUSCEPTIBILITY", &
      sum(this%solvent%waveNumbers), sum(this%grid%waveNumbers), &
      this%grid%waveNumberArraySize; call flush(0)
#endif /*RISM_DEBUG*/

    ! Checking R-grid size.
    if (this%grid%waveNumbers(this%grid%waveNumberArraySize) > this%solvent%waveNumbers(this%solvent%numRDFpoints)) then
       call rism_report_error('(a,1pe16.8,a,1pe16.8)', &
            'DISTVV: bulk solvent Kmax=', this%solvent%waveNumbers(this%solvent%numRDFpoints), &
            'insufficient for 3D-grid Kmax=', this%grid%waveNumbers(this%grid%waveNumberArraySize))
    else if (maxPointsToInterp > this%solvent%numRDFpoints) then
       call rism_report_error('(a,i7,a,i7)', &
            'DISTVV: bulk solvent grid size Nr=', this%solvent%numRDFpoints, &
            'insufficient for interpolation maxPointsToInterp=', maxPointsToInterp)
    end if

    ! Interpolate Xvv(k) on 1D-RISM grid to 3D-RISM grid.
    ! Interpolation is performed about the 1D-RISM wave number just
    ! larger than the current 3D-RISM wave number, with a point range
    ! of +/- maxPointsToInterp / 2.
    do iwn = 1, this%grid%waveNumberArraySize
       ! Find the smallest 1D-RISM wave number just larger than the
       ! current 3D-RISM wave number.  The range starts at
       ! this point - maxPointsToInterp/2 and goes up to + maxPointsToInterp/2.
       do igk = 1, this%solvent%numRDFpoints - maxPointsToInterp + 1
          igk1 = igk
          if (this%solvent%waveNumbers(igk1 + maxPointsToInterp/2) > this%grid%waveNumbers(iwn)) then
             exit
          end if
       end do
       ! Interpolate from 1D-RISM to 3D-RISM grids about the midpoint
       ! +/- maxPointsToInterp/2.
       do iv2 = 1, this%solvent%numAtomTypes
          do iv1 = 1, this%solvent%numAtomTypes
             call polynomialInterpolation( &
                  this%solvent%waveNumbers(igk1:igk1 + maxPointsToInterp), &
                  xvv(igk1:igk1 + maxPointsToInterp, iv1, iv2), maxPointsToInterp, &
                  this%grid%waveNumbers(iwn), &
                  xvva(iwn + (iv1 - 1) * this%grid%waveNumberArraySize &
                  + (iv2 - 1) * this%grid%waveNumberArraySize * this%solvent%numAtomTypes), &
                  err)
          end do
       end do
    end do
  end subroutine interpolateSolventSusceptibility


  !> Center the solute in the solvent box.
  !! @param[in,out] this rism3d object.
  !! FIXME: This needs to be thoroughly inspected prior to merge.
  subroutine centerSolute(this)
    use rism_util, only : findCenterOfMass, translate
    ! following seems not be used (anywhere):
    ! use rism3d_grid_c, only : fractionalCoordFromCartesianCoord
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: gridpoints(3)
    _REAL_ :: mass(this%solute%numAtoms)
    _REAL_ :: projection(3)
    _REAL_ :: unitCellCenter(3)
    integer :: id
    ! If centering is not equal to 0 we want to move the solute to the
    ! center of the solvent box. However, if centering < 0 we only
    ! figure out the displacement required the _first_ time we see the
    ! solute.  Thus, for centering <= 0, the solute's CM can move
    ! relative to the grid.
    if (mod(abs(this%centering), 2) == 1) then
       mass = this%solute%mass
    else if (mod(abs(this%centering), 2) == 0) then
       mass = 1
    end if
    if (this%centering > 0 .or. (this%centering < 0 .and. this%nsolution == 0)) then
       call findCenterOfMass(this%solute%position, this%solute%centerOfMass, &
            mass, this%solute%numAtoms)
    end if
    if (this%centering >= 3 .and. this%centering <= 4) then
       ! Move the center to the nearest multiple of the grid spacing.
       ! 1. Project center of mass vector onto the unit cell axis.
       do id = 1, 3
          projection(id) = dot_product(this%solute%centerOfMass, &
               this%grid%unitCellUnitVectorsR(id, :))
       end do
       ! 2. Convert to grid point space, round to nearest grid point.
       projection = anint(projection / this%grid%spacing)
       ! 3. Project back to orthonormal Cartesian space.
       this%solute%centerOfMass = 0
       do id = 1, 3
          this%solute%centerOfMass = this%solute%centerOfMass &
               + projection(id) * this%grid%voxelVectorsR(id, :)
       end do
    end if
    if (this%centering /= 0) then
       unitCellCenter = 0
       do id = 1, 3
          unitCellCenter = unitCellCenter + this%grid%unitCellVectorsR(id,:)
       end do
       unitCellCenter = unitCellCenter / 2
       call translate(this%solute%position, this%solute%numAtoms, &
            unitCellCenter - this%solute%centerOfMass)
    end if
  end subroutine centerSolute


  !!!!
  !! Subroutines to find the iterative 3D-RISM solution.
  !!!!


  !> Main driver for the 3D-RISM solver.
  !! Makes an initial guess of the direct correlation function and
  !! then solve the RISM and closure relations until either the
  !! solution converges or the maximum of steps is reached.
  !! @param[in,out] this rism3d object.
  !! @param[in] ksave Save itermediate results every ksave interations
  !!  (0 means no saves).
  !! @param[in] kshow Print parameter for relaxation steps every kshow
  !!  iteration (0 means no saves).
  !! @param[in] maxSteps Maximum number of rism relaxation steps.
  !! @param[in] tolerance Tolerance in.
  subroutine solve3DRISM(this, ksave, kshow, maxSteps, tolerance)
    use mdiis_c
    use rism3d_restart
    implicit none
#include "def_time.h"
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: ksave, kshow, maxSteps
    _REAL_, intent(in) :: tolerance
    character(72) :: cuvsav = 'rism.csv', guvfile
    integer :: guv_local = 77

    integer :: iatv, igx, igy, igz
    logical :: found, converged = .false.
    integer :: ig, igr, iv, iv2, istep
    _REAL_ :: residual = 0

    ! Absolute first time in solve3DRISM.
    logical, save :: first = .true.

    ! MPI rank counter.
    integer :: irank
    ! iostat
    integer :: stat
    integer :: ientry, nentry
    ! MPI error.
    logical :: ierr

    call rism_timer_start(this%solve3DRISMTimer)

    if (this%verbose >= 1 .and. size(this%closureList) > 1) &
         call rism_report_message("|Using "// &
         trim(rism3d_closure_type(this%closure))//" closure")

    ! Make initial guess for DCF.
    this%cuvres = 0
    !! FIXME: delete this...
#ifdef MPI
#else
    if (this%mpirank == 0) then
       inquire (file = cuvsav, exist = found)
       if (found .and. first .and. ksave /= 0) then
#ifdef RISM_DEBUG
          write(0,*)'reading saved Cuv file:  ', cuvsav
#endif /*RISM_DEBUG*/
          call readRestartFile(cuvsav, this%cuv(:, :, :, :), &
               this%grid%totalLocalPointsR, this%solvent%numAtomTypes)
       else
#endif /*MPI*/
          !! END FIXME

          if (this%nsolution == 0 .or. this%ncuvsteps == 0) then
             ! Default initial guess for DCF is zero everywhere.
             this%cuv(:,:,:,:) = 0.d0
             ! Add long-range part.  This long-range part is
             ! subtracted in next routine, so the initial guess
             ! for the DCF remains zero.
             if (this%solute%charged .and. .not. this%periodic) then
                do iv = 1, this%solvent%numAtomTypes
                   do igz = 1, this%grid%localDimsR(3)
                      do igy = 1, this%grid%localDimsR(2)
                         do igx = 1, this%grid%localDimsR(1)
                            ig = igx + (igy-1) * this%grid%localDimsR(1) + &
                                 (igz-1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                            this%cuv(igx, igy, igz, iv) = this%cuv(igx, igy, igz, iv) &
                                 + this%solvent%charge(iv) * this%potential%dcfLongRangeAsympR(ig)
                         end do
                      end do
                   end do
                end do
             end if
          end if
#ifdef MPI
#else
       end if
    end if
#endif /*MPI*/

    ! Solve 3D-RISM for the current solute-solvent system.
    !
    ! This is done by first using the DCF guess above, the calculated
    ! potentials, and the 1D-RISM solvent-solvent susceptibility to
    ! solve the 3D-RISM equation for the TCF. These values are then
    ! used in the bridge function to obtain a new solution for TCF.
    ! The difference between the bridge function TCF and the 3D-RISM
    ! TCF is the residual, which is used to create a new guess for the
    ! DCF.
    !
    ! The process repeats until the residual is minimized to a desired
    ! amount, at which point the solution is considered converged.
    !
    ! MDIIS is used to increase the rate of convergence by combining
    ! knowledge of previous guesses of the DCF with the current TCF
    ! residual to give an improved guess for the DCF.

    call mdiis_reset(this%mdiis_o)
    this%cuv => this%cuvWRK(:, :, :, :, mdiis_getWorkVector(this%mdiis_o))
    this%cuvres => this%cuvresWRK(:, :, mdiis_getWorkVector(this%mdiis_o))

    call timer_start(TIME_R1RISM)
    do istep = 1, maxSteps
       !  -----------------------------------------------------------------
       ! One iteration of 3D-RISM and closure relation, advancing w/ mdiis.
       !  -----------------------------------------------------------------
       call single3DRISMsolution(this, residual, converged, tolerance, &
           this%periodic)

       ! Showing selected and last relaxation steps.
       if (kshow /= 0 .and. this%mpirank == 0 .and. this%verbose >= 2) then
          if (converged .or. mod(istep, kshow) == 0 .or. &
               ksave > 0 .and. mod(istep, max(ksave, 1)) == 0) then
             call rism_report_message('(a,i5,5x,a,1pg10.3,5x,a,i3)', &
                  ' Step=', istep, 'Resid=', residual, 'IS=', &
                  getCurrentNVec(this%mdiis_o))
             call rism_report_flush()
          end if
       end if

       ! Exiting relaxation loop on convergence.
       if (converged) exit
    end do

    if (.not. converged) then
       call rism_report_error('(a,i5)', &
          'RXRISM: reached limit # of relaxation steps: ', maxSteps)
    end if
    call timer_stop(TIME_R1RISM)

    first = .false.
    if (this%mpirank == 0 .and. this%verbose >= 1) then
       call rism_report_message('(a,i5,a)', &
           "|RXRISM converged in ", istep, " steps")
    end if
    call rism_timer_stop(this%solve3DRISMTimer)
    return
  end subroutine solve3DRISM

  !> One relaxation step for the UV-RISM equation with the HNC closure,
  !! Guv(r) = exp(-Uuv(r) + Tuv(r) - DelHv0) + DelHv0
  !! Cuv(r) = Guv(r) - 1 - Tvv(r)
  !! Huv(k) = Cuv(k) * (Wvv(k) + Density * Hvv(k))
  !! TuvRes(r) = Huv(r) - Guv(r) - 1
  !! @param[in,out] this A rism3d object.
  !! @param[in,out] residual ???
  !! @param[in,out] converged Returns true if the solution has converged.
  !! @param[in] tolerance Target residual tolerance for convergence.
  subroutine single3DRISMsolution(this, residual, converged, tolerance, &
      periodic )

    use rism3d_fft_c
    implicit none
#include "def_time.h"
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d), intent(inout) :: this
    logical, intent(inout) :: converged
    _REAL_, intent(inout) :: residual
    _REAL_, intent(in) :: tolerance
    logical, intent(in) :: periodic
    integer :: iis
    _REAL_ :: earg, tuv0, tvvr
    integer :: istep

    integer ::  ig1, iga, iv, iv1, iv2, igx, igy, igz, igk
#ifdef FFW_THREADS
    integer :: nthreads, totthreads
    integer, external :: OMP_get_max_threads, OMP_get_num_threads
    logical, external :: OMP_get_dynamic, OMP_get_nested
#endif
    integer :: ierr, irank
    call rism_timer_start(this%single3DRISMsolutionTimer)


    ! --------------------------------------------------------------
    ! Subtract short-range part from Cuv(r) (if not periodic);
    ! Cuv(r) is then loaded into the guv array.
    ! --------------------------------------------------------------

#if defined(MPI)
    do iv = 1, this%solvent%numAtomTypes
       do igz = 1, this%grid%localDimsR(3)
          do igy = 1, this%grid%localDimsR(2)
             if (this%solute%charged .and. .not.periodic) then
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy - 1) * this%grid%localDimsR(1) &
                       + (igz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                   igk = igx + (igy - 1) * (this%grid%localDimsR(1) + 2) &
                        + (igz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
                   this%guv(igk, iv) = this%cuv(igx, igy, igz, iv) &
                        - this%solvent%charge(iv) * this%potential%dcfLongRangeAsympR(ig1)
                end do
             else
                do igx = 1, this%grid%localDimsR(1)
                   igk = igx + (igy-1) * (this%grid%localDimsR(1) + 2) &
                        + (igz-1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
                   this%guv(igk, iv) = this%cuv(igx, igy, igz, iv)
                end do
             end if
             ! Zero out extra space.
             igk = this%grid%localDimsR(1) + 1 + (igy-1) * (this%grid%localDimsR(1) + 2) &
                  + (igz-1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
             this%guv(igk:igk + 1, iv) = 0.d0
          end do
       end do
    end do
#else
    do iv = 1, this%solvent%numAtomTypes
       if (this%solute%charged .and. .not.periodic) then
          do igz = 1, this%grid%localDimsR(3)
             do igy = 1, this%grid%localDimsR(2)
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) &
                             + (igz-1) * this%grid%localDimsR(2) &
                                       * this%grid%localDimsR(1)
                   this%guv(ig1, iv) = this%cuv(igx, igy, igz, iv) &
                        - this%solvent%charge(iv) * this%potential%dcfLongRangeAsympR(ig1)
                end do
             end do
          end do
       else
          do igz = 1, this%grid%localDimsR(3)
             do igy = 1, this%grid%localDimsR(2)
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) &
                             + (igz-1) * this%grid%localDimsR(2) &
                                       * this%grid%localDimsR(1)
                   this%guv(ig1, iv) = this%cuv(igx, igy, igz, iv)
                end do
             end do
          end do
       end if
       ! Zero out extra space.
       this%guv(this%grid%totalLocalPointsR + 1:this%grid%totalLocalPointsK, iv) = 0.d0
    end do
#endif /*defined(MPI)*/

    ! --------------------------------------------------------------
    ! [Short-range part of] Cuv(r) FFT>K.
    ! --------------------------------------------------------------
    call timer_start(TIME_RISMFFT)
    call rism_timer_start(this%fftTimer)
#if defined(MPI)
    call  rism3d_fft_fwd(this%fft, this%guv)
    this%guv(2:this%grid%totalLocalPointsK:2, :) = &
         -this%guv(2:this%grid%totalLocalPointsK:2, :)
#else
    call  rism3d_fft_fwd(this%fft, this%guv)
#endif /*defined(MPI)*/
    call rism_timer_stop(this%fftTimer)
    call timer_stop(TIME_RISMFFT)

    ! --------------------------------------------------------------
    ! Add long-range part to Cuv(k) in K-space.
    ! --------------------------------------------------------------
    if (this%solute%charged .and. .not.periodic) then
       do iv = 1, this%solvent%numAtomTypes
          do ig1 = 1, this%grid%totalLocalPointsK
             this%guv(ig1, iv) = this%guv(ig1, iv) - this%solvent%charge(iv) &
                * this%potential%dcfLongRangeAsympK(ig1)
          end do
       end do
    end if

    ! --------------------------------------------------------------
    ! Huv(k) by RISM.
    ! --------------------------------------------------------------
    do iv1 = 1, this%solvent%numAtomTypes
       do ig1 = 1, this%grid%totalLocalPointsK
          this%huv(ig1, iv1) = 0d0
          iga = this%grid%waveVectorWaveNumberMap((ig1 + 1) / 2)
          do iv2 = 1, this%solvent%numAtomTypes
             this%huv(ig1, iv1) = this%huv(ig1, iv1) &
                  + this%guv(ig1, iv2) * this%xvva(iga + (iv2 - 1) &
                     * this%grid%waveNumberArraySize &
                  + (iv1 - 1) * this%grid%waveNumberArraySize &
                     * this%solvent%numAtomTypes)
          end do
       end do
    end do

    if( .not.periodic ) then
       ! --------------------------------------------------------------
       ! Add long-range part of Huv(k) at k = 0, which was estimated by
       ! long-range part of Cuv(k) at k = 0.
       ! --------------------------------------------------------------
       if (this%mpirank == 0) then
          do ig1 = 1, 2
             do iv = 1, this%solvent%numAtomTypes
                this%huv(ig1, iv) = this%huv(ig1, iv) + this%potential%huvk0(ig1, iv)
             end do
          end do
       end if

       ! --------------------------------------------------------------
       ! Subtract long-range part from huv in K-space.
       ! --------------------------------------------------------------
       if (this%solvent%ionic) then
#if defined(MPI)
          do iv = 1, this%solvent%numAtomTypes
             if (this%mpirank == 0) then
                do ig1 = 3, this%grid%totalLocalPointsK
                   this%huv(ig1, iv) = this%huv(ig1, iv) &
                        + 1d0/this%solvent%dielconst * this%solvent%charge_sp(iv) * this%potential%tcfLongRangeAsympK(ig1)
                end do
             else
                do ig1 = 1, this%grid%totalLocalPointsK
                   this%huv(ig1, iv) = this%huv(ig1, iv) &
                     + 1d0/this%solvent%dielconst &
                        * this%solvent%charge_sp(iv) &
                        * this%potential%tcfLongRangeAsympK(ig1)
                end do
             end if
          end do
#else
          do iv = 1, this%solvent%numAtomTypes
             do ig1 = 3, this%grid%totalLocalPointsK
                this%huv(ig1, iv) = this%huv(ig1, iv) &
                  + 1d0/this%solvent%dielconst * this%solvent%charge_sp(iv) &
                      * this%potential%tcfLongRangeAsympK(ig1)
             end do
          end do
#endif /*defined(MPI)*/
       end if
    end if ! periodic

    ! --------------------------------------------------------------
    ! Short-range part of Huv(k) FFT>R.
    ! --------------------------------------------------------------
    call timer_start(TIME_RISMFFT)
    call rism_timer_start(this%fftTimer)
#if defined(MPI)
    this%huv(2:this%grid%totalLocalPointsK:2, :) = &
         -this%huv(2:this%grid%totalLocalPointsK:2, :)
    call  rism3d_fft_bwd(this%fft, this%huv)
#else
    call  rism3d_fft_bwd(this%fft, this%huv)
#endif /*defined(MPI)*/
    call rism_timer_stop(this%fftTimer)
    call timer_stop(TIME_RISMFFT)

    ! --------------------------------------------------------------
    ! Add long-range part to huv in R-space.
    ! --------------------------------------------------------------
    if (this%solvent%ionic .and. .not.periodic) then
       do iv = 1, this%solvent%numAtomTypes
          do igz = 1, this%grid%localDimsR(3)
             do igy = 1, this%grid%localDimsR(2)
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) + &
                        (igz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#if defined(MPI)
                   igk = igx + (igy - 1) * (this%grid%localDimsR(1) + 2) &
                        + (igz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
                   this%huv(igk, iv) = this%huv(igk, iv) &
                        + this%solvent%charge_sp(iv) * this%potential%tcfLongRangeAsympR(ig1)
#else
                   this%huv(ig1, iv) = this%huv(ig1, iv) &
                        + this%solvent%charge_sp(iv) * this%potential%tcfLongRangeAsympR(ig1)
#endif /*defined(MPI)*/
                end do
             end do
          end do
       end do
    end if

    ! --------------------------------------------------------------
    ! Solve the closure for the RDF.
    ! --------------------------------------------------------------
    call rism3d_closure_guv(this%closure, this%guv, this%huv, this%cuv)

    ! --------------------------------------------------------------
    ! Calculate TCF residual for use in estimating DCF residual.
    ! --------------------------------------------------------------
    this%cuvres(:, :) = 0
    do iv = 1, this%solvent%numAtomTypes
       do igz = 1, this%grid%localDimsR(3)
          do igy = 1, this%grid%localDimsR(2)
             do igx = 1, this%grid%localDimsR(1)
                ig1 = igx + (igy - 1) * this%grid%localDimsR(1) + &
                     (igz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#if defined(MPI)
                igk = igx + (igy - 1) * (this%grid%localDimsR(1) + 2) &
                     + (igz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
#else
                igk = ig1
#endif /*defined(MPI)*/
                this%cuvres(ig1, iv) = this%guv(igk, iv) - 1d0 - this%huv(igk, iv)
             end do
          end do
       end do
    end do

    ! --------------------------------------------------------------
    ! MDIIS
    ! --------------------------------------------------------------
    call timer_start(TIME_MDIIS)
    call mdiis_advance(this%mdiis_o, residual, converged, tolerance)
    this%cuv => this%cuvWRK(:, :, :, :, mdiis_getWorkVector(this%mdiis_o))
    this%cuvres => this%cuvresWRK(:, :, mdiis_getWorkVector(this%mdiis_o))
    call timer_stop(TIME_MDIIS)
    call rism_timer_stop(this%single3DRISMsolutionTimer)

  end subroutine single3DRISMsolution


  subroutine solve3DRISM_dT(this, kshow, maxSteps, tolerance)

  !> Main driver for the 3D-RISM ΔT solver.
  !! IN:
  !!   this :: rism3d object
  !!   kshow  :: print parameter for relaxation steps every kshow iterations
  !!   maxSteps :: maximum number of rism relaxation steps
  !!   mdiis_method :: MDIIS implementation to use

    use rism3d_fft_c
    use mdiis_c
    use rism3d_restart
    implicit none
#include "def_time.h"
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: kshow, maxSteps
    _REAL_, intent(in) :: tolerance
    character(72) ::  cuvsav = 'rism.csv', guvfile
    integer :: guv_local = 77
    !    character*(*)  cuvsav
    ! mdiis :: MDIIS solver object
    type(mdiis), save :: mdiis_o

    integer :: iatv, igx, igy, igz, igk
    logical ::  found, converged = .false.
    integer ::  ig, iv, iv2, istep
    _REAL_ ::  residual = 0

    !_REAL_ ::  delhv0(this%solvent%numAtomTypes)

    logical, save :: first = .true.

    ! irank  :: mpi rank counter
    ! stat :: iostat
    integer :: irank, stat, ientry, nentry
    ! ierr :: mpi error
    logical :: ierr

    integer :: ig1
    call rism_timer_start(this%solve3DRISM_dTTimer)

#ifdef MPI
    call mdiis_new_mpi(mdiis_o, this%mdiis_method, &
         this%deloz, tolerance, &
         this%MDIIS_restart, this%mpirank, this%mpisize, this%mpicomm)
#else
    call mdiis_new(mdiis_o, this%mdiis_method, &
         this%deloz, tolerance, this%MDIIS_restart)
#endif /*MPI*/
    ! call mdiis_setTimerParent(mdiis_o, this%single3DRISM_dTsolutiontimer)
    call mdiis_setData(mdiis_o, this%cuvWRK_dT, this%cuvresWRK_dT, &
         this%grid%totalLocalPointsR * this%solvent%numAtomTypes, this%nvec)
    this%cuv_dT => this%cuvWRK_dT(:, :, :, :, mdiis_getWorkVector(mdiis_o))
    this%cuvres_dT => this%cuvresWRK_dT(:, :, mdiis_getWorkVector(mdiis_o))

    this%cuvres_dT = 0
    !! FIX - delete this...

#ifdef MPI
#else
    if (this%mpirank == 0) then
       !         inquire (file = cuvsav, exist = found)
       !         if (found .and. first .and. ksave/=0) then
       !            call  readRestartFile (cuvsav, this%cuv(:, :, :, :), this%grid%totalLocalPointsR, &
       !                 this%solvent%numAtomTypes)
       !         else
#endif /*MPI*/
       !! ENDFIX
       !..... initial Cuv(r)
       if (first .or. this%ncuvsteps == 0) then
          !       if (.true.) then
          this%cuv_dT = 0
          !..... add long-range part,
          !..... because this long-range part is subtracted in next routine..
          if (this%solute%charged) then
             do iv = 1, this%solvent%numAtomTypes
                do igz = 1, this%grid%localDimsR(3)
                   do igy = 1, this%grid%localDimsR(2)
                      do igx = 1, this%grid%localDimsR(1)
                         ig = igx + (igy-1) * this%grid%localDimsR(1) + &
                              (igz-1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                         this%cuv_dT(igx, igy, igz, iv) = -1.0d0 * this%solvent%charge(iv) &
                              *this%potential%dcfLongRangeAsympR(ig)
                      end do
                   end do
                end do
             end do
          end if
       end if
#ifdef MPI
#else
       !         end if
    end if

#endif /*MPI*/


    !----- Cuv(r) > k ----
    !.....subtract short-range part from Cuv(r)
    !.....short-range part of Cuv(r) is loaded in cuvres array
#if defined(MPI)
    do iv = 1, this%solvent%numAtomTypes
       do igz = 1, this%grid%localDimsR(3)
          do igy = 1, this%grid%localDimsR(2)
             if (this%solute%charged) then
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) &
                        + (igz-1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                   igk = igx + (igy-1) * (this%grid%localDimsR(1) + 2) &
                        + (igz-1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
                   this%cuvk(igk, iv) = this%cuv(igx, igy, igz, iv) &
                        - this%solvent%charge(iv) * this%potential%dcfLongRangeAsympR(ig1)
                end do
             else
                do igx = 1, this%grid%localDimsR(1)
                   igk = igx + (igy - 1) * (this%grid%localDimsR(1) + 2) &
                        + (igz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
                   this%cuvk(igk, iv) = this%cuv(igx, igy, igz, iv)
                end do
             end if
             igk = this%grid%localDimsR(1) + 1 + (igy - 1) * (this%grid%localDimsR(1) + 2) &
                  + (igz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
             this%cuvk(igk:igk +1, iv) =0
          end do
       end do
    end do
#else
    if (this%solute%charged) then
       do iv = 1, this%solvent%numAtomTypes
          do igz = 1, this%grid%localDimsR(3)
             do igy = 1, this%grid%localDimsR(2)
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) &
                        + (igz-1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                   this%cuvk(ig1, iv) = this%cuv(igx, igy, igz, iv) &
                        - this%solvent%charge(iv) * this%potential%dcfLongRangeAsympR(ig1)
                end do
             end do
          end do
       end do
    else
       do iv = 1, this%solvent%numAtomTypes
          do igz = 1, this%grid%localDimsR(3)
             do igy = 1, this%grid%localDimsR(2)
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) &
                        + (igz-1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                   this%cuvk(ig1, iv) = this%cuv(igx, igy, igz, iv)
                end do
             end do
          end do
       end do
    end if
#endif /*defined(MPI)*/

    !.....short-range part of Cuv(r) FFT>K
    call timer_start(TIME_RISMFFT)
    call rism_timer_start(this%fft_dTTimer)
    !! $    do iv = 1, this%solvent%numAtomTypes
#if defined(MPI)
    call  rism3d_fft_fwd(this%fft_cuv, this%cuvk)
    this%cuvk(2:this%grid%totalLocalPointsK:2, :) = &
         -this%cuvk(2:this%grid%totalLocalPointsK:2, :)
#else
    call  rism3d_fft_fwd(this%fft_cuv, this%cuvk)
#endif /*defined(MPI)*/
    !! $    end do
    call rism_timer_stop(this%fft_dTTimer)
    call timer_stop(TIME_RISMFFT)
    !.....add long-range part to Cuv(k) in K-space
    if (this%solute%charged) then
       do iv = 1, this%solvent%numAtomTypes
          do ig1 = 1, this%grid%totalLocalPointsK
             this%cuvk(ig1, iv) = this%cuvk(ig1, iv) -  this%solvent%charge(iv) * this%potential%dcfLongRangeAsympK(ig1)
          end do
       end do
    end if


    !.......................... relaxing UV RISM  ..........................
#ifdef RISM_DEBUG
    write(0, *)'relaxing 3D uv RISM DT:'
    call flush(0)
#endif

    do istep = 1, maxSteps

       !................... one relaxation step of UV RISM ....................
       call timer_start(TIME_R1RISM)
       call single3DRISMsolution_dT(this, residual, converged, mdiis_o, tolerance)
       call timer_stop(TIME_R1RISM)


       !............. showing selected and last relaxation steps ..............
       if (kshow /= 0 .and. this%mpirank == 0 .and. this%verbose >= 2) then
          if (converged .or. mod(istep, kshow) == 0) then
             call rism_report_message('(a, i5, 5x, a,1pg10.3, 5x, a,i3)', &
                  ' Step=', istep, 'Resid=', residual, 'IS=', getCurrentNVec(mdiis_o))
             call rism_report_flush()
          end if
       end if

       !! FIX _ DELETE this
       !.............. saving selected and last relaxation steps ..............
#ifdef MPI
#else
       !         if (ksave /= 0 .and. first) then
       !            if (converged .or. ksave > 0 .and. mod(istep, ksave) == 0) then
       !               call  writeRestartFile (cuvsav, this%cuv(:, :, :, :), this%grid%totalLocalPointsR, this%solvent%numAtomTypes)
       !            end if
       !         end if
#endif /*MPI*/
       !! endfix
       !............... exiting relaxation loop on convergence ................
       if (converged) exit
    end do

    if (.not. converged) then
       call rism_report_error('(a,i5)','RXRISMDT: reached limit # of relaxation steps: ', maxSteps)
    end if
    first = .false.
    if (this%mpirank == 0 .and. this%verbose >= 1) then
       call rism_report_message('(a,i5,a)', "|RXRISMDT converged in ", istep)!, " steps")
    end if
    call mdiis_destroy(mdiis_o)
    this%cuv_dT => this%cuvWRK_dT(:, :, :, :, 1)
    this%cuvres_dT => this%cuvresWRK_dT(:, :, 1)
    call rism_timer_stop(this%solve3DRISM_dTTimer)
  end subroutine solve3DRISM_dT


  !! One relaxation step for the UV-RISM equation with the HNC
  !! closure,
  !! Guv(r) = exp(-this%potential%uuv(r) + Tuv(r) - DelHv0) + DelHv0
  !! Cuv(r) = Guv(r) - 1 - Tvv(r)
  !! Huv(k) = Cuv(k) * (Wvv(k) + Density * Hvv(k))
  !! TuvRes(r) = Huv(r) - Guv(r) - 1
  !! IN:
  !!  this :: rism3d object
  !!  residual ::
  !!  converged ::
  !!  mdiis  :: MDIIS object to accelerate convergence
  subroutine single3DRISMsolution_dT(this, residual, converged, mdiis_o, tolerance)
    use rism3d_fft_c
    use mdiis_c
    implicit none
#include "def_time.h"
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d) :: this
    logical ::  converged

    type(mdiis)::mdiis_o
    integer :: iis
    _REAL_, intent(inout) ::  residual
    _REAL_, intent(in) :: tolerance
    _REAL_ :: earg, tuv0, tvvr
    integer :: istep

    integer ::  ig1, iga, iv, iv1, iv2, igx, igy, igz, igk
#ifdef FFW_THREADS
    integer :: nthreads, totthreads
    integer, external :: OMP_get_max_threads, OMP_get_num_threads
    logical, external :: OMP_get_dynamic, OMP_get_nested
#endif

#ifdef RISM_DEBUG
    write(0,*)"R1RISMDT"
    call flush(0)
#endif
    call rism_timer_start(this%single3DRISM_dTsolutionTimer)

    ! Subtract short-range part from direct correlation function c(r).
    ! Short-range part of c(r) is loaded in guv array.
#if defined(MPI)
    do iv = 1, this%solvent%numAtomTypes
       do igz = 1, this%grid%localDimsR(3)
          do igy = 1, this%grid%localDimsR(2)
             if (this%solute%charged) then
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) &
                        + (igz-1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                   igk = igx + (igy-1) * (this%grid%localDimsR(1) + 2) &
                        + (igz-1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
                   this%guv_dT(igk, iv) = this%cuv_dT(igx, igy, igz, iv) &
                        + this%solvent%charge(iv) * this%potential%dcfLongRangeAsympR(ig1)
                end do
             else
                do igx = 1, this%grid%localDimsR(1)
                   igk = igx + (igy-1) * (this%grid%localDimsR(1) + 2) &
                        + (igz-1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
                   this%guv_dT(igk, iv) = this%cuv_dT(igx, igy, igz, iv)
                end do
             end if
             igk = this%grid%localDimsR(1) + 1 + (igy-1) * (this%grid%localDimsR(1) + 2) &
                  + (igz-1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
             this%guv_dT(igk:igk +1, iv) =0d0
          end do
       end do
    end do
#else
    if (this%solute%charged) then
       do iv = 1, this%solvent%numAtomTypes
          do igz = 1, this%grid%localDimsR(3)
             do igy = 1, this%grid%localDimsR(2)
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) &
                        + (igz-1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                   this%guv_dT(ig1, iv) = this%cuv_dT(igx, igy, igz, iv) &
                        + this%solvent%charge(iv) * this%potential%dcfLongRangeAsympR(ig1)
                end do
             end do
          end do
       end do
    else
       do iv = 1, this%solvent%numAtomTypes
          do igz = 1, this%grid%localDimsR(3)
             do igy = 1, this%grid%localDimsR(2)
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) &
                        + (igz-1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                   this%guv_dT(ig1, iv) = this%cuv_dT(igx, igy, igz, iv)
                end do
             end do
          end do
       end do
    end if
#endif /*defined(MPI)*/
    ! Short-range part of c(r) FFT>K.
    call timer_start(TIME_RISMFFT)
    call rism_timer_start(this%fft_dTTimer)
#if defined(MPI)
    call  rism3d_fft_fwd(this%fft_dT, this%guv_dT)
    this%guv_dT(2:this%grid%totalLocalPointsK:2, :) = &
         -this%guv_dT(2:this%grid%totalLocalPointsK:2, :)
#else
    call  rism3d_fft_fwd(this%fft_dT, this%guv_dT)
#endif /*defined(MPI)*/
    call rism_timer_stop(this%fft_dTTimer)
    call timer_stop(TIME_RISMFFT)
    ! Add long-range part to c(k) in K-space.
    if (this%solute%charged) then
       do iv = 1, this%solvent%numAtomTypes
          do ig1 = 1, this%grid%totalLocalPointsK
             this%guv_dT(ig1, iv) = this%guv_dT(ig1, iv) +  this%solvent%charge(iv) * this%potential%dcfLongRangeAsympK(ig1)
          end do
       end do
    end if
    ! h(k) by RISM.
    do iv1 = 1, this%solvent%numAtomTypes
       do ig1 = 1, this%grid%totalLocalPointsK
          this%huv_dT(ig1, iv1) = 0d0
          iga = this%grid%waveVectorWaveNumberMap((ig1 + 1)/2)
          do iv2 = 1, this%solvent%numAtomTypes
             this%huv_dT(ig1, iv1) = this%huv_dT(ig1, iv1) + &
                  this%cuvk(ig1, iv2) * this%xvva_dT(iga + (iv2-1) * this%grid%waveNumberArraySize + &
                  (iv1 - 1) * this%grid%waveNumberArraySize * this%solvent%numAtomTypes)+ &
                  this%guv_dT(ig1, iv2) * this%xvva(iga + (iv2-1) * this%grid%waveNumberArraySize + &
                  (iv1 - 1) * this%grid%waveNumberArraySize * this%solvent%numAtomTypes)
          end do
       end do
    end do

    ! Add long-range part of h(k) at k = 0
    ! which was estimated by long-range part of c(k) at k = 0.
    if (this%mpirank == 0) then
       do ig1 = 1, 2
          do iv = 1, this%solvent%numAtomTypes
             this%huv_dT(ig1, iv) = this%huv_dT(ig1, iv) + this%potential%huvk0_dT(ig1, iv)
          end do
       end do
    end if

    ! Subtract long-range part from h(k) in k-space.
    if (this%solvent%ionic) then
#if defined(MPI)
       do iv = 1, this%solvent%numAtomTypes
          if (this%mpirank == 0) then
             do ig1 = 3, this%grid%totalLocalPointsK
                this%huv_dT(ig1, iv) = this%huv_dT(ig1, iv) &
                     - 1d0/this%solvent%dielconst * this%solvent%charge_sp(iv) * this%potential%tcfLongRangeAsympK(ig1)
             end do
          else
             do ig1 = 1, this%grid%totalLocalPointsK
                this%huv_dT(ig1, iv) = this%huv_dT(ig1, iv) &
                     - 1d0/this%solvent%dielconst * this%solvent%charge_sp(iv) * this%potential%tcfLongRangeAsympK(ig1)
             end do
          end if
       end do
#else
       do iv = 1, this%solvent%numAtomTypes
          do ig1 = 3, this%grid%totalLocalPointsK
             this%huv_dT(ig1, iv) = this%huv_dT(ig1, iv) &
                  - 1d0/this%solvent%dielconst * this%solvent%charge_sp(iv) * this%potential%tcfLongRangeAsympK(ig1)
          end do
       end do
#endif /*defined(MPI)*/
    end if

    ! Short-range part of h(k) FFT>R.
    call timer_start(TIME_RISMFFT)
    call rism_timer_start(this%fft_dTTimer)
#if defined(MPI)
    this%huv_dT(2:this%grid%totalLocalPointsK:2, :)= &
         -this%huv_dT(2:this%grid%totalLocalPointsK:2, :)
    call  rism3d_fft_bwd(this%fft_dT, this%huv_dT)
#else
    call  rism3d_fft_bwd(this%fft_dT, this%huv_dT)
#endif /*defined(MPI)*/
    call rism_timer_stop(this%fft_dTTimer)
    call timer_stop(TIME_RISMFFT)

    if (this%solvent%ionic) then
       !.....add long-range part to huv in R-space
       do iv = 1, this%solvent%numAtomTypes
          do igz = 1, this%grid%localDimsR(3)
             do igy = 1, this%grid%localDimsR(2)
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) + &
                        (igz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#if defined(MPI)
                   igk = igx + (igy - 1) * (this%grid%localDimsR(1) + 2) &
                        + (igz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
                   this%huv_dT(igk, iv) = this%huv_dT(igk, iv) &
                        -  this%solvent%charge_sp(iv) * this%potential%tcfLongRangeAsympR(ig1)
#else
                   this%huv_dT(ig1, iv) = this%huv_dT(ig1, iv) &
                        -  this%solvent%charge_sp(iv) * this%potential%tcfLongRangeAsympR(ig1)
#endif /*defined(MPI)*/
                end do
             end do
          end do
       end do
    end if
    !    this%cuvres_dT(:, :) = 0

    call rism3d_closure_guv_dT(this%closure, this%guv_dT, this%huv_dT, this%cuv_dT, &
         this%guv, this%huv, this%cuv)
    do iv = 1, this%solvent%numAtomTypes
       do igz = 1, this%grid%localDimsR(3)
          do igy = 1, this%grid%localDimsR(2)
             do igx = 1, this%grid%localDimsR(1)
                ig1 = igx + (igy-1) * this%grid%localDimsR(1) + &
                     (igz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#if defined(MPI)
                igk = igx + (igy - 1) * (this%grid%localDimsR(1) + 2) &
                     + (igz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
#else
                igk = ig1
#endif /*defined(MPI)*/
                this%cuvres_dT(ig1, iv) = this%guv_dT(igk, iv) - this%huv_dT(igk, iv)
             end do
          end do
       end do
    end do
    call timer_start(TIME_MDIIS)
    call mdiis_advance(mdiis_o, residual, converged, tolerance)
    this%cuv_dT => this%cuvWRK_dT(:, :, :, :, mdiis_getWorkVector(mdiis_o))
    this%cuvres_dT => this%cuvresWRK_dT(:, :, mdiis_getWorkVector(mdiis_o))
    call timer_stop(TIME_MDIIS)
    call rism_timer_stop(this%single3DRISM_dTsolutionTimer)

  end subroutine single3DRISMsolution_dT


  ! PROPAGATE PREVIOUS SOLUTIONS

  !> Calculates a new initial guess for CUV based on the final solutions
  !! from previous timesteps.  The maximum number of previous time
  !! steps to use is provided by the user in ncuvsteps.  However, if
  !! there are not enough previous timesteps only nsolution previous
  !! timesteps will be used.
  !!
  !! See section 2.3.1 and eqs. 8-13 of doi:10.1021/ct900460m for
  !! details.
  !! @param[in] this rism3d object.
  subroutine guessDCF(this)
    implicit none
    type(rism3d) :: this
    integer :: iv, n
    n = this%grid%totalLocalPointsR * this%solvent%numAtomTypes
    if (this%ncuvsteps >= 5 .and. this%nsolution >= 5) then
       call dscal(n, 5d0, this%cuv(:, :, :, :), 1)
       call daxpy(n, -10d0, this%oldcuv(:, :, :, :, 2), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, 10d0, this%oldcuv(:, :, :, :, 3), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, -5d0, this%oldcuv(:, :, :, :, 4), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, 1d0, this%oldcuv(:, :, :, :, 5), 1, this%cuv(:, :, :, :), 1)
    else if (this%ncuvsteps >= 4 .and. this%nsolution >= 4) then
       call dscal(n, 4d0, this%cuv(:, :, :, :), 1)
       call daxpy(n, -6d0, this%oldcuv(:, :, :, :, 2), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, 4d0, this%oldcuv(:, :, :, :, 3), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, -1d0, this%oldcuv(:, :, :, :, 4), 1, this%cuv(:, :, :, :), 1)
    else if (this%ncuvsteps >= 3 .and. this%nsolution >= 3) then
       call dscal(n, 3d0, this%cuv(:, :, :, :), 1)
       call daxpy(n, -1d0, this%oldcuv(:, :, :, :, 2), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, 1d0, this%oldcuv(:, :, :, :, 3), 1, this%cuv(:, :, :, :), 1)
    else if (this%ncuvsteps >= 2 .and. this%nsolution >= 2) then
       call dscal(n, 2d0, this%cuv(:, :, :, :), 1)
       call daxpy(n, -1d0, this%oldcuv(:, :, :, :, 2), 1, this%cuv(:, :, :, :), 1)
    else if (this%ncuvsteps == 0) then
       this%cuv(:, :, :, :) = 0
    end if
  end subroutine guessDCF


  !> Updates the values in the this%oldcuv queue.  The oldest value (the
  !! ncuvstep index) is pushed out, the remainder of the data is
  !! shifted and the newest solution is placed in the first index.
  !! @param[in] this rism3d object.
  subroutine updateDCFguessHistory(this)
    implicit none
    type(rism3d) :: this
    integer :: istep, iv
#ifdef RISM_DEBUG
    write(0, *) "DCF_GUESS_HISTORY_UPDATE"; call flush(0)
#endif /*RISM_DEBUG*/
    if (this%ncuvsteps == 0) return
    do istep = min(this%ncuvsteps, this%nsolution), 2, -1
       call dcopy(this%grid%totalLocalPointsR * this%solvent%numAtomTypes, this%oldcuv(:, :, :, :, istep-1), 1, &
            this%oldcuv(:, :, :, :, istep), 1)
    end do
    call dcopy(this%grid%totalLocalPointsR * this%solvent%numAtomTypes, this%cuv(:, :, :, :), 1, &
         this%oldcuv(:, :, :, :, 1), 1)
  end subroutine updateDCFguessHistory

  !> Convert the user supplied a,b,a1,b1 Universal Correction
  !! coefficients to a0,b0,a1,b1 coefficients used in
  !! rism3d_closure_c.
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   coeff: coefficients for correction.  For the original correction
  !!          (a, b) = coeff(1:2).
  !!          Extra coefficients are a1 and b1 from Johnson et al. 2016.
  !! OUT:
  !!   array of length 4 containing a0,b0,a1,b1
  function UC_temperature_coeff(this,coeff) result(tcoeff)
    implicit none
    type(rism3d) :: this
    _REAL_, intent(in) :: coeff(:)
    _REAL_ :: tcoeff(4)
    ! a0 & a1
    if (size(coeff) > 2) then
       tcoeff(1) = coeff(1) - coeff(3)*this%solvent%temperature
       tcoeff(3) = coeff(3)
    else
       tcoeff(1) = coeff(1)
       tcoeff(3) = 0d0
    end if
    ! b0 & b1
    if (size(coeff) > 3) then
       tcoeff(2) = coeff(2) - coeff(4)*this%solvent%temperature
       tcoeff(4) = coeff(4)
    else
       tcoeff(2) = coeff(2)
       tcoeff(4) = 0d0
    end if
  end function UC_temperature_coeff

  !> Create an electron density map from a 3D solute-solvent RDF by
  !! smearing the 3D RDF with a 1D solvent electron denisty RDF.
  !! TODO: Currently only water oxygen is supported for electron
  !! smearing.
  !! @param[in] this rism3d object.
  !! @param[in] electronRDF Solvent 1D electron density map.
  !! @param[in] totalSolventElectrons Total electrons in solvent. Ex.: for water, Z = 10.
  !! @param[out] electronMap Resulting smeared electron density map.
  subroutine createElectronDensityMap(this, iv, electronRDF, &
       electronRDFGridSpacing, totalSolventElectrons, density, electronMap)
    implicit none
    type(rism3d), intent(in) :: this
    integer, intent(in) :: iv
    _REAL_, intent(in) :: electronRDF(:)
    _REAL_, intent(in) :: electronRDFGridSpacing
    integer, intent(in) :: totalSolventElectrons
    _REAL_, intent(in) :: density
    _REAL_, intent(out) :: electronMap(:,:,:)

    integer :: numSmearGridPoints(3)
    integer :: igx, igy, igz, igk, igxCenter, igyCenter, igzCenter
    integer :: igCenter
    _REAL_ :: rx(3), ry(3), rz(3), rxCenter(3), ryCenter(3), rzCenter(3)
    _REAL_ :: distance
    _REAL_ :: numElectronsAtGridCenter
    integer :: numUnitCells(3)

    _REAL_ :: electronRDFSum
    integer :: centerGridIndex(3), smearGridIndex(3)
    _REAL_ :: centerGridPoint(3), smearGridPoint(3)
    _REAL_ :: distanceVector(3)
    integer :: rdfIndex

    !TODO: Make this user defined?
    _REAL_, parameter :: cutoff = 2.0
    ! Value near zero.
    _REAL_, parameter :: nearZero = 1E-10

    electronMap(:,:,:) = 0

    ! Number of grid points to smear in each direction.
    numSmearGridPoints = floor(cutoff / this%grid%spacing)

    ! Determine how many electrons to place in the center grid point.
    do igz = -numSmearGridPoints(3), numSmearGridPoints(3)
       rz = igz * this%grid%voxelVectorsR(3, :)
       do igy = -numSmearGridPoints(2), numSmearGridPoints(2)
          ry = igy * this%grid%voxelVectorsR(2, :)
          do igx = -numSmearGridPoints(1), numSmearGridPoints(1)
             rx = igx * this%grid%voxelVectorsR(1, :)

             if (igx == 0 .and. igy == 0 .and. igz == 0) cycle

             distance = sqrt(dot_product((/ rx, ry, rz /), (/ rx, ry, rz /)))
             rdfIndex = ceiling(distance / electronRDFGridSpacing)
             ! If the rdfIndex falls outside the range then assume
             ! there are no electron outside, which is reasonable for
             ! an RDF with a 4 Angstrom maximum distance.
             if (rdfIndex <= size(electronRDF)) then
                electronRDFSum = electronRDFSum + &
                     electronRDF(rdfIndex) * this%grid%voxelVolume
             end if
          end do
       end do
    end do
    ! Sum contains all contribution from grid points within the cutoff
    ! EXCEPT the center grid point.
    numElectronsAtGridCenter = totalSolventElectrons - electronRDFSum

#ifdef _OPENMP_
! #pragma omp parallel for schedule(dynamic, 10) shared (dx, conc, result, stepx, stepy, stepz, center_grid)
!$omp parallel do private(local_equal), shared(this, density, electronMap, numSmearGridPoints, numElectronsAtGridCenter)
#endif
    do igzCenter = 0, this%grid%globalDimsR(3) - 1
       rzCenter = igzCenter * this%grid%voxelVectorsR(3, :)
       do igyCenter = 0, this%grid%globalDimsR(2) - 1
          ryCenter = igyCenter * this%grid%voxelVectorsR(2, :)
          do igxCenter = 0, this%grid%globalDimsR(1) - 1
             rxCenter = igxCenter * this%grid%voxelVectorsR(1, :)

             igCenter = 1 + igxCenter + igyCenter * this%grid%globalDimsR(1) + &
                  igzCenter * this%grid%globalDimsR(2) * this%grid%globalDimsR(1)

             centerGridIndex = (/ igxCenter, igyCenter, igzCenter /)
             centerGridPoint = rxCenter + ryCenter + rzCenter

             if (this%guv(igCenter, iv) > nearZero) then
                ! Center grid.
!$omp critical
                electronMap(igxCenter + 1, igyCenter + 1, igzCenter + 1) = &
                     electronMap(igxCenter + 1, igyCenter + 1, igzCenter + 1) + &
                     numElectronsAtGridCenter * &
                     this%guv(igCenter, iv) * density
!$omp end critical
                do igz = igzCenter - numSmearGridPoints(3), igzCenter + numSmearGridPoints(3)
                   rz = igz * this%grid%voxelVectorsR(3, :)
                   do igy = igyCenter - numSmearGridPoints(2), igyCenter + numSmearGridPoints(2)
                      ry = igy * this%grid%voxelVectorsR(2, :)
                      do igx = igxCenter - numSmearGridPoints(1), igxCenter + numSmearGridPoints(1)
                         rx = igx * this%grid%voxelVectorsR(1, :)

                         smearGridIndex = (/ igx, igy, igz /)
                         smearGridPoint = rx + ry + rz

                         distanceVector = smearGridPoint - centerGridPoint
                         distance = sqrt(dot_product(distanceVector, distanceVector))
                         if (distance .le. cutoff) then
                            ! Outside of box?
                            if (any(smearGridIndex < 0) .or. any(smearGridIndex >= this%grid%globalDimsR)) then
                               if (this%periodic) then
                                  ! Jump to grid image inside box.
                                  numUnitCells = floor(real(smearGridIndex) / this%grid%globalDimsR)
                                  smearGridIndex = smearGridIndex - numUnitCells * this%grid%globalDimsR
                               else
                                  cycle
                               end if
                            ! Skip center grid point.
                            else if (all(smearGridIndex == centerGridIndex)) then
                               cycle
                            end if

                            rdfIndex = ceiling(distance / electronRDFGridSpacing)
                            if (rdfIndex .le. size(electronRDF)) then
!$omp critical
                               smearGridIndex = smearGridIndex + 1
                               electronMap(smearGridIndex(1), smearGridIndex(2), smearGridIndex(3)) = &
                                    electronMap(smearGridIndex(1), smearGridIndex(2), smearGridIndex(3)) &
                                    + electronRDF(rdfIndex) * this%grid%voxelVolume &
                                    * this%guv(igCenter, iv) * density
!$omp end critical
                            end if
                         end if
                      end do
                   end do
                end do
             end if
          end do
       end do
    end do
#ifdef _OPENMP_
!$omp end parallel do
#endif
  end subroutine createElectronDensityMap


end module rism3d_c
