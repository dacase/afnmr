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

!> Object for solute information to be used in 3D-RISM.
module rism3d_solute_c
  use safemem
  use rism_report_c
  implicit none
  
  type rism3d_solute
     !> Number of solute atoms.
     integer :: numAtoms = 0
     ! Mass of each solute atom. [au]
     _REAL_, pointer :: mass(:) => NULL()
     ! Partial charge of each solute atom. [sqrt(kT A)]
     _REAL_, pointer :: charge(:) => NULL()
     !> Holds the orignal values in charge() when _unsetcharges()
     !! is used.  Otherwise it is not used. [sqrt(kT A)]
     _REAL_, pointer :: origCharge(:) => NULL()
     !> Cartesian coordinates of each solute atom (dim, atom). [A]
     _REAL_, pointer :: position(:,:) => NULL()
     !> Lennard-Jones parameter sigma* or r_min. [A]
     _REAL_, pointer :: ljSigma(:) => NULL()
     !> Lennard-Jones parameter epsilon. [kT]
     !! U_LJ = eps * ((ljSigma / r)**12 - 2 * (ljSigma / r)**6)
     _REAL_, pointer :: ljEpsilon(:) => NULL()

     !> False if all partial charges are zero.
     logical :: charged
     !> Sum of solute atom partial charges. [sqrt(kT A)]
     _REAL_ :: totalCharge
     !> Center of mass (or geometry) of the entire solute. [A]
     _REAL_ :: centerOfMass(3)
  end type rism3d_solute

  interface rism3d_solute_new
     module procedure rism3d_solute_new_internal, rism3d_solute_new_sander, &
          rism3d_solute_new_parm7
  end interface rism3d_solute_new

  interface rism3d_solute_setCoord
     module procedure rism3d_solute_setCoord_ratu, rism3d_solute_setCoord_crd
  end interface rism3d_solute_setCoord

  private :: allocateSolute, setCharge

contains


  !> Constructor using internal data representation.
  !! Used during cloning of solute objects.
  !! @param[in] this Solute object.
  !! @param[in] natom Number of atoms.
  !! @param[in] mass Mass of each atom. [au]
  !! @param[in] charge Charge of each atom. [sqrt(kT A)]
  !! @param[in] ratu Atom coordinates. [A]
  !! @param[in] ljSigma r_min / 2. [A]
  !! @param[in] eps LJ epsilon. [kT]
  !! @param[in] o_mpicomm (optional) MPI communicator.
  subroutine rism3d_solute_new_internal (this, natom, mass, charge, ratu, ljSigma, eps, o_mpicomm)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    type(rism3d_solute), intent(inout) :: this
    integer, intent(in) :: natom
    _REAL_, intent(in) :: mass(natom), charge(natom), ratu(3, natom), ljSigma(natom), eps(natom)
  
    integer, optional, intent(in) :: o_mpicomm
    integer :: mpicomm, mpirank, mpisize, err
    mpicomm = 0
    mpirank = 0
    mpisize = 1
#ifdef MPI
    if (present(o_mpicomm)) then
       mpicomm = o_mpicomm
       if (mpicomm == MPI_COMM_NULL) &
            call rism_report_error("RISM3D_SOLUTE: received NULL MPI communicator")
       call mpi_comm_rank(mpicomm, mpirank, err)
       if (err /=0) call rism_report_error&
            ("(a, i8)", "RISM3D SOLU: could not get MPI rank for communicator ", mpicomm)
       call mpi_comm_size(mpicomm, mpisize, err)
       if (err /=0) call rism_report_error&
            ("(a, i8)", "RISM3D SOLU: could not get MPI size for communicator ", mpisize)
    end if
#endif
    if (mpirank == 0) then
       this%numAtoms = natom
       call allocateSolute(this)
       call setCharge(this, charge)
       this%mass = mass
       this%ljSigma = ljSigma
       this%ljEpsilon = eps
  
       call rism3d_solute_setCoord(this, ratu)
    end if
#ifdef MPI
    if (mpisize /=1) &
         call rism3d_solute_mpi_clone(this, mpirank, mpicomm)
#endif /*MPI*/
  end subroutine rism3d_solute_new_internal


  !> Constructor for dealing with Amber's Sander internal molecular
  !! representation.  Currently used for both NAB and Sander.
  !! @param[in, out] this Solute object.
  !! @param[in] natom Number of atoms.
  !! @param[in] numTypes Number of atom types.
  !! @param[in] atomTypeIndex
  !! @param[in] nonbondedParmIndex
  !! @param[in] charge Charge of each atom. [e*18.2223]
  !! @param[in] ljA LJ A parameter. [kcal/mol*A^12]
  !! @param[in] ljB LJ B parameter. [kcal/mol*A^6]
  !! @param[in] mass Mass of each atom. [au]
  !! @param[in] temperature Temperature of the solvent.
  !! @param[in] gridSpacing Number of spaces between grid points in each dimension.
  !! @param[in] o_mpicomm (optional) MPI communicator.
  subroutine rism3d_solute_new_sander (this, numAtoms, numTypes, atomTypeIndex, &
       nonbondedParmIndex, charge, ljA, ljB, mass, temperature, &
       o_mpicomm)
    use constants, only: KB, COULOMB_CONST_E, BOLTZMANN, &
         AVOGADRO, AMBER_ELECTROSTATIC, PI
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    type(rism3d_solute), intent(inout) :: this
    integer, intent(in) :: numAtoms, atomTypeIndex(numAtoms), nonbondedParmIndex(numTypes**2), numTypes
    _REAL_, intent(in) :: charge(numAtoms), ljA(numTypes * (numTypes + 1) / 2), &
         ljB(numTypes * (numTypes + 1) / 2), mass(numAtoms), temperature
    integer, optional, intent(in) :: o_mpicomm

    integer :: id, iu, iv, ir    

    integer :: mpicomm, mpirank, mpisize, err
    mpicomm = 0
    mpirank = 0
    mpisize = 1
#ifdef MPI
    if (present(o_mpicomm)) then
       mpicomm = o_mpicomm
       if (mpicomm == MPI_COMM_NULL) &
            call rism_report_error("RISM3D_SOLUTE: received NULL MPI communicator")
       call mpi_comm_rank(mpicomm, mpirank, err)
       if (err /=0) call rism_report_error&
            ("(a, i8)", "RISM3D SOLU: could not get MPI rank for communicator ", mpicomm)
       call mpi_comm_size(mpicomm, mpisize, err)
       if (err /=0) call rism_report_error&
            ("(a, i8)", "RISM3D SOLU: could not get MPI size for communicator ", mpisize)
    end if
#endif
    
    if (mpirank == 0) then

       this%numAtoms = numAtoms
       call allocateSolute(this)
       ! Unit conversion from [e*18.2223] to [sqrt(kT A)].
       call setCharge(this, charge(1:numAtoms) / AMBER_ELECTROSTATIC &
            * sqrt(COULOMB_CONST_E / (KB * temperature)))
       this%mass = mass
       do iu = 1, numAtoms
          id = nonbondedParmIndex(numTypes * (atomTypeIndex(iu) - 1) + atomTypeIndex(iu))
          !! FIXME: Is this the proper way of handling no LJ?  No.  The
          !! best is to find the volume of the enclosing atom and
          !! calculating a radius for the enclosed atom that coincides
          !! with this and a well depth of 1/10 the parent atom.
          if (ljB(id) == 0d0) then
             this%ljSigma(iu) = 0.7d0
          else
             this%ljSigma(iu) = (2 * ljA(id) / ljB(id))**(1.d0 / 6.d0) / 2d0
          end if
       end do
       do iu = 1, numAtoms
          id = nonbondedParmIndex(numTypes * (atomTypeIndex(iu) - 1) + atomTypeIndex(iu))
          if (ljA(id) == 0d0) then
             this%ljEpsilon(iu) = 1d-2
          else
             this%ljEpsilon(iu) = ljB(id)**2 / (4 * ljA(id))
             this%ljEpsilon(iu) = this%ljEpsilon(iu)
          end if
          this%ljEpsilon(iu) = this%ljEpsilon(iu) / (KB * temperature)
       end do       

    end if
#ifdef MPI
    if (mpisize /= 1) &
         call rism3d_solute_mpi_clone(this, mpirank, mpicomm)
#endif /*MPI*/
  end subroutine rism3d_solute_new_sander


  !> Constructor for dealing with Amber parm7 and crd7 files.
  !! @param[in, out] this Solute object.
  !! @param[in] parm Parameter file name.
  !! @param[in] temperature Temperature of the solvent. [K]
  !! @param[in] o_mpicomm (optional) MPI communicator.
  subroutine rism3d_solute_new_parm7(this, parm, temperature, o_mpicomm)
    use rism_util, only : freeUnit
    use rism_parm
    implicit none
    type(rism3d_solute), intent(inout) :: this
    character(len=*), intent(in) :: parm
    _REAL_, intent(in) :: temperature
    integer, optional, intent(in) :: o_mpicomm

    character(len=80) fmt
    character(len=80) fmtin, ifmt, afmt, rfmt, type, line
    integer :: i, iok, unit
    integer :: numAtoms, numTypes, nbonh, mbona, ntheth, mtheta, nphih, mphia, &
         nhparm, nparm, nnb, nres, nbona, ntheta, nphia, &
         numbnd, numang, nptra, natyp, nphb, ifpert, nbper, ngper, ndper, &
         mbper, mgper, mdper, ifbox, nmxrs, ifcap, numextra, ncopy
    integer, pointer :: atomTypeIndex(:)=>NULL(), nonbondedParmIndex(:)=>NULL()
    _REAL_, pointer :: charge(:)=>NULL(), ljA(:)=>NULL(), &
         ljB(:)=>NULL(), mass(:)=>NULL()
    integer :: mpicomm
    mpicomm = 0
    if (present(o_mpicomm)) mpicomm = o_mpicomm

    !this liberally uses code from sander's rdparm.f
    unit = freeUnit()
    open (unit, file=parm, status='old')

    ifmt = '(12I6)'
    afmt = '(20A4)'
    rfmt = '(5E16.8)'
    call rism_parm_nxtsec_init()

    fmtin = ifmt
    type = 'POINTERS'
    call rism_parm_nxtsec(unit,  6,  0, fmtin,  type,  fmt,  iok)
    read(unit, fmt) numAtoms, numTypes, nbonh, mbona, ntheth, mtheta, nphih, mphia, &
         nhparm, nparm, nnb, nres, nbona, ntheta, nphia, &
         numbnd, numang, nptra, natyp, nphb, ifpert, nbper, ngper, ndper, &
         mbper, mgper, mdper, ifbox, nmxrs, ifcap, numextra, ncopy

    !allocate some temporary storage.  This uses more memory than we
    !would like but we can use the Sander constructor directly this
    !way
    atomTypeIndex => safemem_realloc(atomTypeIndex, numAtoms, .false.)
    nonbondedParmIndex => safemem_realloc(nonbondedParmIndex, numTypes**2, .false.)
    charge => safemem_realloc(charge, numAtoms, .false.)
    ljA => safemem_realloc(ljA, numTypes * (numTypes + 1) / 2, .false.)
    ljB => safemem_realloc(ljB, numTypes * (numTypes + 1) / 2, .false.)
    mass => safemem_realloc(mass, numAtoms, .false.)

    fmtin = ifmt
    type = 'ATOM_TYPE_INDEX'
    call rism_parm_nxtsec(unit,  6,  0, fmtin,  type,  fmt,  iok)
    read(unit, fmt) atomTypeIndex

    fmtin = ifmt
    type = 'NONBONDED_PARM_INDEX'
    call rism_parm_nxtsec(unit,  6,  0, fmtin,  type,  fmt,  iok)
    read(unit, fmt) nonbondedParmIndex

    fmtin = rfmt
    type = 'CHARGE'
    call rism_parm_nxtsec(unit,  6,  0, fmtin,  type,  fmt,  iok)
    read(unit, fmt) charge

    fmtin = rfmt
    type = 'LENNARD_JONES_ACOEF'
    call rism_parm_nxtsec(unit,  6,  0, fmtin,  type,  fmt,  iok)
    read(unit, fmt) ljA

    fmtin = rfmt
    type = 'LENNARD_JONES_BCOEF'
    call rism_parm_nxtsec(unit,  6,  0, fmtin,  type,  fmt,  iok)
    read(unit, fmt) ljB

    fmtin = rfmt
    type = 'MASS'
    call rism_parm_nxtsec(unit,  6,  0, fmtin,  type,  fmt,  iok)
    read(unit, fmt) mass

    !FIXME: Shouldn't this be where the box dimensions come from?
    ! fmtin = rfmt
    ! type = 'BOX_DIMENSIONS'
    ! call rism_parm_nxtsec(unit,  6,  0, fmtin,  type,  fmt,  iok)
    ! read(unit, fmt) unitCellDimensions

    close(unit)
    call rism3d_solute_new_sander (this, numAtoms, numTypes, atomTypeIndex, nonbondedParmIndex, charge, ljA, ljB, &
         mass, temperature, mpicomm)

    if (safemem_dealloc(atomTypeIndex) /= 0) call rism_report_error("RISM3D_SOLUTE: Failed to deallocate ATOMTYPEINDEX")
    if (safemem_dealloc(nonbondedParmIndex) /= 0) call rism_report_error("RISM3D_SOLUTE: Failed to deallocate NONBONDEDPARMINDEX")
    if (safemem_dealloc(charge) /= 0) call rism_report_error("RISM3D_SOLUTE: Failed to deallocate CHARGE")
    if (safemem_dealloc(ljA) /= 0) call rism_report_error("RISM3D_SOLUTE: Failed to deallocate LJA")
    if (safemem_dealloc(ljB) /= 0) call rism_report_error("RISM3D_SOLUTE: Failed to deallocate LJB")
    if (safemem_dealloc(mass) /= 0) call rism_report_error("RISM3D_SOLUTE: Failed to deallocate MASS")

  end subroutine rism3d_solute_new_parm7


  !> Clone constructor using internal data representation.
  !! @param[in] this Object to be copied.
  !! @param[in, out] clone Clone of the object.  No memory space is
  !!   shared.
  subroutine rism3d_solute_clone (this, clone)
    implicit none
    type(rism3d_solute), intent(inout) :: clone
    type(rism3d_solute), intent(in) :: this

    call rism3d_solute_new(clone, this%numAtoms, this%mass, this%charge, this%position, &
         this%ljSigma, this%ljEpsilon)
  end subroutine rism3d_solute_clone


#ifdef MPI
  !> Allocates memory on non-master nodes and then distributes information out
  !! from the master.  It is assumed that the object on the master
  !! already exists.
  !! @param[in, out] this Solute object.
  !! @param[in] rank MPI rank.
  !! @param[in] comm MPI communicator.
  subroutine rism3d_solute_mpi_clone(this, rank, comm)
    implicit none
    type(rism3d_solute), intent(inout) :: this
    integer, intent(in) :: rank, comm
    integer :: err
    include 'mpif.h'
    !first distribute the pieces of information needed to allocate memory
    call mpi_bcast(this%numAtoms, 1,mpi_integer, 0,comm, err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLUTE: could not broadcast NUMATOMS")

    !non-master processes should now allocate memory
    if (rank /= 0) then
       call allocateSolute(this)
    end if

    !now distribute the arrays to the non-master processes
    call mpi_bcast(this%charge, size(this%charge), mpi_double_precision, 0,comm, err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLUTE: could not broadcast CHARGE")
    call mpi_bcast(this%mass, size(this%mass), mpi_double_precision, 0,comm, err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLUTE: could not broadcast MASS")
    call mpi_bcast(this%position, size(this%position, 1)*size(this%position, 2), mpi_double_precision, 0,comm, err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLUTE: could not broadcast RATU")
    call mpi_bcast(this%ljSigma, size(this%ljSigma), mpi_double_precision, 0,comm, err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLUTE: could not broadcast SIG_S")
    call mpi_bcast(this%ljEpsilon, size(this%ljEpsilon), mpi_double_precision, 0,comm, err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLUTE: could not broadcast EPS")

    if (rank /=0) call setCharge(this, this%charge)
  end subroutine rism3d_solute_mpi_clone
#endif /*MPI*/


  !> Set the solute coordinates
  !! @param[in, out] this Solute object.
  !! @param[in] ratu Coordinates.
  subroutine rism3d_solute_setCoord_ratu(this, ratu)
    implicit none
    type(rism3d_solute), intent(inout) :: this
    _REAL_, intent(in) :: ratu(:, :)
    if (ubound(ratu, 1) /= 3 .and. ubound(ratu, 2) /= this%numAtoms) then
       call rism_report_error("solute size and coordinate array do not match")
       stop
    end if
    this%position = ratu
  end subroutine rism3d_solute_setCoord_ratu


  !> Set the solute coordinates.
  !! @param[in, out] this Solute object.
  !! @param[in] crd Coordinate file. Amber crd/rst or PDB 3.2 file.
  subroutine rism3d_solute_setCoord_crd(this, crd)
    use rism_util, only : freeUnit
    implicit none
    type(rism3d_solute), intent(inout) :: this
    character(len=*), intent(in) :: crd
    character(len=128) :: test
    integer :: unit, numAtoms, iostat
    unit = freeUnit()
    open (unit, file = crd, status = 'old')

    ! First try Amber format.
    ! Read off title and discard.
    read(unit, *)
    ! Read in the number of atoms, we can discard the time if it is there.
    read(unit, *) test
    read(test, *, iostat = iostat) numAtoms
    if (numAtoms == this%numAtoms .and. iostat == 0) then
       read(unit, '(6F12.7)') this%position
    else
       ! Try PDB format.
       numAtoms = 0
       rewind(unit)
       iostat = 0
       do
          read(unit, '(a)', iostat = iostat) test
          if (iostat < 0) exit
          if (iostat > 0) call rism_report_error("RISM3D_SOLUTE: failed on reading PDB")
          if (test(1:6) .ne. "ATOM  " .and. test(1:6) .ne. "HETATM") cycle
          numAtoms = numAtoms + 1
          if (numAtoms > this%numAtoms) then
             call rism_report_error('(a, i6, a)', &
                  "RISM3D_SOLUTE: more PDB atoms than present in parameter file (", &
                  this%numAtoms, ")")
          end if
          read(test(31:54), '(3(f8.3))') this%position(:, numAtoms)
       end do
       if (numAtoms /= this%numAtoms) &
            call rism_report_error("RISM3D_SOLUTE: failed reading coordinate file")
    end if
    close(unit)
  end subroutine rism3d_solute_setCoord_crd


  !> Indicates if the partial charges are on or not.
  !! @param[in, out] this Solute object.
  !! @return True if partial charges are set, false if they are no.
  function rism3d_solute_charged(this) result(charged)
    implicit none
    type(rism3d_solute), intent(inout) :: this
    logical :: charged
    charged=.not.associated(this%origCharge)
  end function rism3d_solute_charged


  !> Sets all solute partial charges to zero.
  !! @param[in, out] this Solute object.
  subroutine rism3d_solute_unsetCharges(this)
    implicit none
    type(rism3d_solute), intent(inout) :: this
    this%origCharge => safemem_realloc(this%origCharge, ubound(this%charge, 1))
    this%origCharge = this%charge
    call setCharge(this, (/0d0/))
  end subroutine rism3d_solute_unsetCharges


  !> Sets all solute partial charges to to their original
  !! values. (Undoes rism3d_solute_unsetCharge().)
  !! @param[in, out] this Solute object.
  subroutine rism3d_solute_resetCharges(this)
    implicit none
    type(rism3d_solute), intent(inout) :: this
    call setCharge(this, this%origCharge)
    if (safemem_dealloc(this%origCharge)/=0) &
         call rism_report_error("unable to deallocate 'origCharge'")
  end subroutine rism3d_solute_resetCharges


  !> Destructor.
  !! @param[in] this Solute object.
  subroutine rism3d_solute_destroy(this)
    implicit none
    type(rism3d_solute) :: this
    if (safemem_dealloc(this%mass) /=0 ) &
         call rism_report_error("RISM3D_SOLUTE: could not deallocate MASS")
    if (safemem_dealloc(this%charge) /=0 ) &
         call rism_report_error("RISM3D_SOLUTE: could not deallocate CHARGE")
    if (safemem_dealloc(this%origCharge) /=0 ) &
         call rism_report_error("RISM3D_SOLUTE: could not deallocate ORIGCHARGE")
    if (safemem_dealloc(this%position) /=0 ) &
         call rism_report_error("RISM3D_SOLUTE: could not deallocate RATU")
    if (safemem_dealloc(this%ljSigma) /=0 ) &
         call rism_report_error("RISM3D_SOLUTE: could not deallocate SIG_s")
    if (safemem_dealloc(this%ljEpsilon) /=0 ) &
         call rism_report_error("RISM3D_SOLUTE: could not deallocate EPS")
  end subroutine rism3d_solute_destroy


  !> Allocate major per-atom solute properties.
  !! @param[in, out] this Solute object.
  subroutine allocateSolute(this)
    implicit none
    type(rism3d_solute), intent(inout) :: this
    this%mass => safemem_realloc(this%mass, this%numAtoms, .false.)
    this%charge => safemem_realloc(this%charge, this%numAtoms, .false.)
    this%ljSigma => safemem_realloc(this%ljSigma, this%numAtoms, .false.)
    this%ljEpsilon => safemem_realloc(this%ljEpsilon, this%numAtoms, .false.)
    this%position => safemem_realloc(this%position, 3, this%numAtoms, .false.)
    !ensure that ratu is initialized in case we have to do an MPI broadcast
    !before getting coordinates
    this%position = 0d0
  end subroutine allocateSolute


  !> Private subroutine to take care of book keeping when setting charges.
  !! @param[in, out] this Solute object.
  !! @param[in] charge Array of charges, one per charge site or a single
  !!   value may be passed which will be applied to all sites.
  subroutine setCharge(this, charge)
    implicit none
    type(rism3d_solute), intent(inout) :: this
    _REAL_, intent(in) :: charge(:)
    if (ubound(charge, 1) /= ubound(this%charge, 1) .and. &
         ubound(charge, 1) /= 1) then
       call rism_report_error("Bad solute charge set")
    end if
    if (ubound(charge, 1) /= 1) then
       this%charge = charge
    else
       this%charge = charge(1)
    end if
    this%totalCharge = sum(this%charge)
    if (sum(abs(this%charge)) > 0d0) then
       this%charged = .true.
    else
       this%charged = .false.
    end if
  end subroutine setCharge
end module rism3d_solute_c
