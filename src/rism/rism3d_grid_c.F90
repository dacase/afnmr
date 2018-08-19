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

!> Grid class for 3D-RISM.  Handles the size and shape of the grid and
!! manages precalculated tables of wavevectors.
!!
!! The class is not MPI aware but explicitly supports grid
!! decomposition. The size of the global grid is stored as *Global*
!! and *Local* is the local grid.
module rism3d_grid_c
  use safemem
  implicit none

  type rism3d_grid
     !> Grid spacing for all of the grids. [A]
     _REAL_ :: spacing(3) = 0.5d0
     !> Volume of a grid cell (the smallest volumetric grid element). [A^3]
     !TODO: Voxel = 3D grid point, so this is incorrect terminology.
     _REAL_ :: voxelVolume
     !> Box size for 3D-RISM.  For PBC calculations, these should
     !! generally be equal. [A]
     _REAL_ :: boxLength(3)
     !> Volume of the box we are doing the calculation. [A^3]
     _REAL_ :: boxVolume

     !> Number of global r-space grid spaces in each dimension.
     integer :: globalDimsR(3)
     !> Number of global k-space grid spaces in each dimension.
     integer :: globalDimsK(3)
     !> Number of MPI local r-space grid spaces in each dimension.
     integer :: localDimsR(3)
     !> Number of MPI local k-space grid spaces in each dimension.
     integer :: localDimsK(3)
     !> r-space grid offset from (0, 0, 0) for MPI calculations.
     integer :: offsetR(3)
     !> k-space grid offset from (0, 0, 0) for MPI calculations.
     integer :: offsetK(3)

     !> Total number of MPI local r-space grid points.
     integer :: totalLocalPointsR
     !> Total number of MPI local k-space grid points.
     integer :: totalLocalPointsK
     !> Total number of global r-space grid points.
     integer :: totalGlobalPointsR

     !> Number of indices in array of wave vector magnitudes.
     integer :: waveNumberArraySize
     !> Absolute value of the wave vector, |k|, also known as the wave
     !! number (this%waveNumberArraySize long).
     _REAL_, pointer :: waveNumbers(:) => NULL()
     !> Wave vector at every k-space grid point (vector, k-space grid points).
     _REAL_, pointer :: waveVectors(:,:) => NULL()
     !> Self dot product of the wave vector at every k-space grid point.
     _REAL_, pointer :: waveVectors2(:) => NULL()
     !> Map wave vector index to wave number index in ascending wave
     !! number value.
     integer, pointer :: waveVectorWaveNumberMap(:) => NULL()

     !> True if grid has a defined periodic unit cell.
     logical :: periodic

     !TODO: Swap indices to be (vector, axis) since this is used
     ! throughout rest of codebase.
     !> Unit cell side lengths. [A]
     _REAL_ :: unitCellLengths(3)
     !> Unit cell interior angles. [radians]
     _REAL_ :: unitCellAngles(3)
     !> Unit cell coordinate vectors of length unity, in Cartesian
     !! coordinates of the box grid.  These are derived from the unit
     !! cell interior angles (axis, vector). [A]
     _REAL_ :: unitCellUnitVectorsR(3, 3)
     ! _REAL_ :: unitCellUnitVectorsK(3, 3)
     !> Unit cell coordinate vectors with lengths of the respective
     !! unit cell side, in Cartesian coordinates of the box grid.
     !! These are derived from the unit cell side lengths and interior
     !! angles (axis, vector). [A]
     _REAL_ :: unitCellVectorsR(3, 3)
     _REAL_ :: unitCellVectorsK(3, 3)
     !> Unit cell coordinate vectors with length of the grid spacing
     !! in each respective dimension. These are derived from the grid
     !! spacing and unit cell interior angles (axis, vector). [A]
     _REAL_ :: voxelVectorsR(3, 3)
     _REAL_ :: voxelVectorsK(3, 3)

     ! !> Unit cell lattice vectors in frequency space.
     ! !! TODO: This should be equivalent to unitCellVectorsK, no?
     ! _REAL_ :: unitCellRecipVectors(3, 3)

     _REAL_ :: inscribedSphereRadius

     !> MPI rank number.
     integer :: mpirank
     !> Number of MPI processes.
     integer :: mpisize
     !> MPI communicator.
     integer :: mpicomm
  end type rism3d_grid

  private setup_wavevector


contains

  !> Constructor for a new rism3d_grid.
  !! @param[in,out] this Grid object.
  !! @param[in] mpicomm MPI communicator.
  subroutine rism3d_grid_new(this, mpicomm)
    implicit none
    type(rism3d_grid), intent(inout) :: this
    integer, optional, intent(in) :: mpicomm
    if (present(mpicomm)) then
       call rism3d_grid_setmpi(this, mpicomm)
    else
       call rism3d_grid_setmpi(this, 0)
    end if
  end subroutine rism3d_grid_new

  
  !> Sets the MPI communicator.  If compiled without MPI, rank = 0 and
  !! size = 1.
  !! @param[in] this Grid object.
  !! @param[in] comm MPI communicator.
  subroutine rism3d_grid_setmpi(this, comm)
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif
    type(rism3d_grid), intent(inout) :: this
    integer, intent(in) :: comm
    integer :: err
    this%mpicomm = comm
    this%mpirank = 0
    this%mpisize = 1
#if defined(MPI)
    if (this%mpicomm == MPI_COMM_NULL)&
         call rism_report_error("RISM3D_GRID: received NULL MPI communicator")
    call mpi_comm_rank(this%mpicomm, this%mpirank, err)
    if (err /= 0) call rism_report_error&
         ("(a, i8)","RISM3D GRID: could not get MPI rank for communicator ", comm)
    call mpi_comm_size(this%mpicomm, this%mpisize, err)
    if (err /= 0) call rism_report_error&
         ("(a, i8)","RISM3D GRID: could not get MPI size for communicator ", comm)
#endif /*defined(MPI)*/

  end subroutine rism3d_grid_setmpi


  !> Sets the grid spacing for future calculations.
  !! @param[in,out] this Grid object.
  !! @param[in] grdpsc Grid spacing in each dimension. [A]
  subroutine rism3d_grid_setSpacing(this, grdspc)
    implicit none
    type(rism3d_grid), intent(inout) :: this
    _REAL_, intent(in) :: grdspc(3)

    this%spacing = grdspc
  end subroutine rism3d_grid_setSpacing


  !> Constructor for a new rism3d_grid.
  !! @param[in,out] this Grid object.
  !! @param[in] gridSpacing Space between consecutive grid points in each dimension. [A]
  !! @param[in] globalDimsR Number of global r-space grid points in each dimension.
  !! @param[in] globalDimsK Number of global k-space grid points in each dimension.
  !! @param[in] localDimsR Number of r-space grid points in each dimension.
  !! @param[in] localDimsK Number of k-space grid points in each dimension.
  !! @param[in] offsetR R-space offset from (0, 0, 0).
  !! @param[in] offsetK K-space offset from (0, 0, 0).
  !! @param[in] unitCellDimensions Unit cell dimensions lengths and beta angle.
  subroutine rism3d_grid_resize(this, gridSpacing, globalDimsR, &
       globalDimsK, localDimsR, localDimsK, offsetR, offsetK)
    use rism_util, only: cross, magnitude
    implicit none
    type(rism3d_grid), intent(inout) :: this
    integer, intent(in) :: globalDimsR(3), globalDimsK(3)
    integer, intent(in) :: localDimsR(3), localDimsK(3)
    integer, intent(in) :: offsetR(3), offsetK(3)
    _REAL_, optional, intent(in) :: gridSpacing(3)

    _REAL_ :: unitCellDimensions(6)

    integer id, idd
    _REAL_ :: width_a, width_b, width_c

    this%spacing = gridSpacing

    this%globalDimsR = globalDimsR
    this%globalDimsK = globalDimsK
    this%localDimsR = localDimsR
    this%localDimsK = localDimsK
    this%offsetR = offsetR
    this%offsetK = offsetK

    this%totalLocalPointsR = product(this%localDimsR)
    this%totalLocalPointsK = product(this%localDimsK)
    this%totalGlobalPointsR = product(this%globalDimsR)

    do id = 1, 3
       this%voxelVectorsR(id, :) = this%spacing(id) * this%unitCellUnitVectorsR(id, :)
       do idd = 1, 3
          if (this%voxelVectorsR(id, idd) /= 0) then
             this%voxelVectorsK(id, idd) = 1 / this%voxelVectorsR(id, idd)
          else
             this%voxelVectorsK(id, idd) = 0
          end if
       end do
    end do

    ! this%voxelVectorsK = this%voxelVectorsR

    this%boxLength = this%globalDimsR * this%spacing

    !TODO: Temporary workaround since boxLength is not defined for
    ! buffer + grid space combo until here.
    if (.not. this%periodic) then
       unitCellDimensions = (/ &
            this%boxLength(1), &
            this%boxLength(2), &
            this%boxLength(3), &
            90d0, 90d0, 90d0 /)
       call rism3d_grid_setUnitCellDimensions(this, unitCellDimensions, this%periodic)
    end if
    
    if (this%periodic) then
       ! Area of base times height
       ! = cross product of base with the dot product of height.
       this%boxVolume = dot_product( &
            & cross((this%unitCellVectorsR(1,:)), (this%unitCellVectorsR(2,:))), &
            & this%unitCellVectorsR(3,:))

       this%voxelVolume = dot_product( &
            & cross((this%voxelVectorsR(1,:)), (this%voxelVectorsR(2,:))), &
            & this%voxelVectorsR(3,:))
    else
       this%boxVolume = product(this%boxLength)
       this%voxelVolume = product(this%spacing)
    end if

    width_a = abs(dot_product(this%unitCellVectorsR(1,:), &
         cross(this%unitCellVectorsR(2,:), this%unitCellVectorsR(3,:)))) &
         / abs(magnitude(cross(this%unitCellVectorsR(2,:), this%unitCellVectorsR(3,:))))
    width_b = abs(dot_product(this%unitCellVectorsR(2,:), &
         cross(this%unitCellVectorsR(3,:), this%unitCellVectorsR(1,:)))) &
         / abs(magnitude(cross(this%unitCellVectorsR(3,:), this%unitCellVectorsR(1,:))))
    width_c = abs(dot_product(this%unitCellVectorsR(3,:), &
         cross(this%unitCellVectorsR(1,:), this%unitCellVectorsR(2,:)))) &
         / abs(magnitude(cross(this%unitCellVectorsR(1,:), this%unitCellVectorsR(2,:))))
    this%inscribedSphereRadius = 0.5 * min(width_a, width_b, width_c)

    this%waveVectors => safemem_realloc(this%waveVectors, 3, this%totalLocalPointsK / 2, .false.)
    this%waveVectors2 => safemem_realloc(this%waveVectors2, this%totalLocalPointsK / 2, .false.)
    this%waveVectorWaveNumberMap => safemem_realloc(this%waveVectorWaveNumberMap, this%totalLocalPointsK / 2, .false.)

    call setup_wavevector(this)
  end subroutine rism3d_grid_resize

  
  !> Set the solute unit cell dimensions.
  !! @param[in, out] this Grid object.
  !! @param[in] unitCellDimensions
  !!   Unit cell dimensions in order:
  !!       a, b, c, alpha, beta, gamma.
  !!   Distances are in Angstroms and angles are in degrees.
  !! @param[in] periodic
  !!   Determine if the 'unit cell' represents an isolated
  !!   box or a periodic unit.
  !! isolated simulation box
  subroutine rism3d_grid_setUnitCellDimensions(this, unitCellDimensions, periodic)
    use constants, only: PI
    use rism_util, only: cross ! rotationMatrixFromEulerAngles
    implicit none
    type(rism3d_grid), intent(inout) :: this
    _REAL_, intent(in) :: unitCellDimensions(6)
    logical, intent(in) :: periodic

    _REAL_ :: u12(3), u23(3), u31(3), unitCellVolume

    integer :: id, ie

    this%periodic = periodic

    this%unitCellLengths = unitCellDimensions(:3)
    this%unitCellAngles = unitCellDimensions(4:) * PI / 180

    ! Code from ew_box.F90 in sander; each vector is of unit length.
    this%unitCellUnitVectorsR(1, :) = (/ 1d0, 0d0, 0d0 /)
    this%unitCellUnitVectorsR(2, :) = (/ cos(this%unitCellAngles(3)), &
         sin(this%unitCellAngles(3)), 0d0 /)
    this%unitCellUnitVectorsR(3, 1) = cos(this%unitCellAngles(2))
    this%unitCellUnitVectorsR(3, 2) =  &
         (cos(this%unitCellAngles(1)) - this%unitCellUnitVectorsR(3, 1) * &
         this%unitCellUnitVectorsR(2, 1)) / this%unitCellUnitVectorsR(2, 2)
    this%unitCellUnitVectorsR(3, 3) = &
         sqrt(1d0 - this%unitCellUnitVectorsR(3, 1)**2 &
         - this%unitCellUnitVectorsR(3, 2)**2)

    ! Avoid errors introduced by Fortran trigonometric functions which
    ! sometimes produce near-zero values when the exact answer is
    ! zero.
    do ie = 1, 3
       do id = 1, 3
          if (abs(this%unitCellUnitVectorsR(id, ie)) < 1E-10) then
             this%unitCellUnitVectorsR(id, ie) = 0
          end if
       end do
    end do

    this%unitCellVectorsR(1, :) = this%unitCellLengths(1) * this%unitCellUnitVectorsR(1, :)
    this%unitCellVectorsR(2, :) = this%unitCellLengths(2) * this%unitCellUnitVectorsR(2, :)
    this%unitCellVectorsR(3, :) = this%unitCellLengths(3) * this%unitCellUnitVectorsR(3, :)

    this%voxelVectorsR(1, :) = this%spacing(1) * this%unitCellUnitVectorsR(1, :)
    this%voxelVectorsR(2, :) = this%spacing(2) * this%unitCellUnitVectorsR(2, :)
    this%voxelVectorsR(3, :) = this%spacing(3) * this%unitCellUnitVectorsR(3, :)

    ! Get the reciprocal basis vectors (from same sander code as above).
    ! See The Basics of Crystallography Diffraction by Christopher
    ! Hammond, p. 100, ISBN 0198559453.
    u23 = cross(this%unitCellVectorsR(2, :), this%unitCellVectorsR(3, :))
    u31 = cross(this%unitCellVectorsR(3, :), this%unitCellVectorsR(1, :))
    u12 = cross(this%unitCellVectorsR(1, :), this%unitCellVectorsR(2, :))
    unitCellVolume = dot_product(this%unitCellVectorsR(1, :), u23)
    this%unitCellVectorsK(1, :) = u23(:) / unitCellVolume
    this%unitCellVectorsK(2, :) = u31(:) / unitCellVolume
    this%unitCellVectorsK(3, :) = u12(:) / unitCellVolume

#if 0
    if( this%periodic ) then
       write(6,'(a)') 'Unit cell vectors:'
       write(6,'(3f12.5)') this%unitCellVectorsR(1,:)
       write(6,'(3f12.5)') this%unitCellVectorsR(2,:)
       write(6,'(3f12.5)') this%unitCellVectorsR(3,:)
       write(6,'(a)') 'Reciprocal cell vectors:'
       write(6,'(3f12.5)') this%unitCellVectorsK(1,:)
       write(6,'(3f12.5)') this%unitCellVectorsK(2,:)
       write(6,'(3f12.5)') this%unitCellVectorsK(3,:)
    endif
#endif

  end subroutine rism3d_grid_setUnitCellDimensions

  
  !> Print the properties of a set of unit cell axis vectors,
  !! including the vectors themselves, their magnitudes, and angles.
  !! Currently only used for debugging.
  !! @param[in] vectors Set of three unit cell axis vectors.
  !! @param[in] vectorsLabel Text briefly describing the unit cell
  !!                         vectors.
  subroutine printUnitCellVectorProperties(vectors, vectorsLabel)
    use constants, only: PI
    use rism_util, only: magnitude
    implicit none
    _REAL_, intent(in) :: vectors(3, 3)
    character(len=*), intent(in) :: vectorsLabel
    integer :: id, ie
    
    print *, vectorsLabel, " vectors: ", vectors(1, :)
    print *, vectors(2, :)
    print *, vectors(3, :)

    print *, vectorsLabel, " vector magnitudes: "
    do id = 1, 3
       print *, magnitude(vectors(id, :))
    end do

    print *, vectorsLabel, " vector angles: "
    do id = 1, 2
       do ie = 2, 3
          print *, id, ie, acos(dot_product(vectors(id, :), vectors(ie, :)) &
               / (magnitude(vectors(id, :)) * magnitude(vectors(ie, :)))) &
               * 180 / PI
       end do
    end do
  end subroutine printUnitCellVectorProperties
  

  !> Destroys a rism3d_grid object.
  !! @param[in,out] this Grid object.
  subroutine rism3d_grid_destroy(this)
    implicit none
    type(rism3d_grid), intent(inout) :: this

    this%mpirank = 0
    this%mpisize = 0
    this%mpicomm = 0
    if (safemem_dealloc(this%waveNumbers) /= 0) then
       call rism_report_error("RISM3D_GRID: failed to deallocate ga")
    end if
    if (safemem_dealloc(this%waveVectors) /= 0) then
       call rism_report_error("RISM3D_GRID: failed to deallocate gv")
    end if
    if (safemem_dealloc(this%waveVectors2) /= 0) then
       call rism_report_error("RISM3D_GRID: failed to deallocate g2")
    end if
    if (safemem_dealloc(this%waveVectorWaveNumberMap) /= 0) then
       call rism_report_error("RISM3D_GRID: failed to deallocate indga")
    end if
  end subroutine rism3d_grid_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Does first round of set up after a new solvation box size has been set.
  !! Specifically the code:
  !!  - Calculates the wave vector and wave vector self-product.
  !!  - Calculates, sorts and indexes the wave number.
  !!  - Allocates memory for the wave number array.
  !! @param[in,out] this rism3d_grid object.
  subroutine setup_wavevector(this)
    use constants, only : PI
    use safemem
    use rism_util, only: indexArray, checksum
    implicit none
    type(rism3d_grid), intent(inout) :: this

    integer :: gDimX, gDimY, gDimZ
    integer :: lDimX, lDimY, lDimZ
    integer :: lgx, lgy, lgz
    integer :: igx, igy, igz
    integer :: igk, igs
    _REAL_ ::  waveX, waveY, waveZ
    _REAL_ :: waveVector2, largestWaveVector2
    integer, pointer :: waveVectorToWaveVector2Map(:) => NULL()

    integer :: ierr

    waveVectorToWaveVector2Map => safemem_realloc(waveVectorToWaveVector2Map, &
         this%totalLocalPointsK / 2, .false.)

    ! Getting and checking box arrays sizes.
    gDimX = this%globalDimsR(1)
    gDimY = this%globalDimsR(2)
    gDimZ = this%globalDimsR(3)
    ! gDimX = this%globalDimsK(1)
    ! gDimY = this%globalDimsK(2)
    ! gDimZ = this%globalDimsK(3)
    ! print *, 'globalDimsR', this%globalDimsR
    ! print *, 'globalDimsK', this%globalDimsK

#if defined(MPI)
    ! Initializing wave vector for real FFT in wrap-around order.
    do igz = 0, this%localDimsK(3) - 1
       do igy = 0, this%localDimsK(2) - 1
          do igx = 0, this%localDimsK(1) / 2 - 1
             igk = 1 + igx + &
                  igy * (this%localDimsK(1) / 2) + &
                  igz * this%localDimsK(2) * (this%localDimsK(1) / 2)

             lgx = mod(igx + &
                  (gDimX / 2 - 1), gDimX) - (gDimX / 2 - 1)
             lgy = mod(igz + this%offsetK(3) + &
                  (gDimY / 2 - 1), gDimY) - (gDimY / 2 - 1)
             lgz = mod(igy + &
                  (gDimZ / 2 - 1), gDimZ) - (gDimZ / 2 - 1)

             ! Calculating wave vector and its self-product.
             if (.not. this%periodic) then
                waveX = 2d0 * PI / this%boxLength(1) * lgx
                waveY = 2d0 * PI / this%boxLength(2) * lgy
                waveZ = 2d0 * PI / this%boxLength(3) * lgz
                this%waveVectors(1, igk) = waveX
                this%waveVectors(2, igk) = waveY
                this%waveVectors(3, igk) = waveZ
                this%waveVectors2(igk) = waveX**2 + waveY**2 + waveZ**2
             else
                this%waveVectors(:,igk) = &
                     2d0 * PI * (this%unitCellVectorsK(1,:) * lgx &
                     + this%unitCellVectorsK(2,:) * lgy &
                     + this%unitCellVectorsK(3,:) * lgz)
                this%waveVectors2(igk) = dot_product(this%waveVectors(:,igk), &
                     this%waveVectors(:,igk))
             end if
          end do
       end do
    end do

#else
    
    ! Initializing wave vector for real-space FFT in wrap-around order 
    ! (as in Numerical Recipes).
    do igz = 0, this%globalDimsR(3) - 1
       ! Center the indices about zero, so indices less than
       ! half the max dimension are positive increasing (1 to
       ! N/2) and indices over half are negative increasing
       ! (-N/2 to -1).
       ! This allows cancellation of the complex parts in a Fourier series.
       lgz = mod(igz + (this%globalDimsR(3) / 2 - 1), &
            this%globalDimsR(3)) - (this%globalDimsR(3) / 2 - 1)
       do igy = 0, this%globalDimsR(2) - 1
          lgy = mod(igy + (this%globalDimsR(2) / 2 - 1), &
               this%globalDimsR(2)) - (this%globalDimsR(2) / 2 - 1)
          do igx = 0, this%globalDimsR(1) / 2 - 1
             !TODO: igx only goes up to half grid points, so this
             ! wrapping is unnecessary.
             lgx = mod(igx + (this%globalDimsR(1) / 2 - 1), &
                  this%globalDimsR(1)) - (this%globalDimsR(1) / 2 - 1)
             
             ! Over the range igk = 1, gDimX / 2 * gDimY * gDimZ.
             igk = 1 + igx + igy * this%globalDimsR(1) / 2 + igz * this%globalDimsR(2) * this%globalDimsR(1) / 2
             
             ! Calculate wave vector and its self-product.
             if (.not. this%periodic) then
                !   dac note: next three lines look like orthogonal code(?)
                waveX = 2d0 * PI / this%boxLength(1) * lgx
                waveY = 2d0 * PI / this%boxLength(2) * lgy
                waveZ = 2d0 * PI / this%boxLength(3) * lgz
                this%waveVectors(1, igk) = waveX
                this%waveVectors(2, igk) = waveY
                this%waveVectors(3, igk) = waveZ
                this%waveVectors2(igk) = waveX**2 + waveY**2 + waveZ**2
                ! write(0,'(a,2i5,3f12.5)'), 'orthog:  ', igk, lgx, this%waveVectors(:,igk)
             else
                this%waveVectors(:,igk) = &
                     2d0 * PI * (this%unitCellVectorsK(1,:) * lgx &
                     + this%unitCellVectorsK(2,:) * lgy &
                     + this%unitCellVectorsK(3,:) * lgz)
                this%waveVectors2(igk) = dot_product(this%waveVectors(:,igk), &
                     this%waveVectors(:,igk))
                ! write(0,'(a,2i5,3f12.5)'), 'northog: ', igk, lgx, this%waveVectors(:,igk)
             end if
          end do
       end do
    end do
    ! Initializing Nyquist frequency of wave vector for real-space FFT.
    do igz = 0, this%globalDimsR(3) - 1
       do igy = 0, this%globalDimsR(2) - 1
          igx = this%globalDimsR(1) / 2

          ! Over the range igk = 1 + gDimX / 2 * gDimY * gDimZ,
          ! gDimY * gDimZ + gDimx / 2 * gDimY * gDimZ.
          igk = 1 + igy + igz * this%globalDimsR(2) &
               + this%globalDimsR(1) / 2 * this%globalDimsR(2) * this%globalDimsR(3)
          
          lgx = mod(igx + (this%globalDimsR(1) / 2 - 1), this%globalDimsR(1)) - (this%globalDimsR(1) / 2 - 1)
          lgy = mod(igy + (this%globalDimsR(2) / 2 - 1), this%globalDimsR(2)) - (this%globalDimsR(2) / 2 - 1)
          lgz = mod(igz + (this%globalDimsR(3) / 2 - 1), this%globalDimsR(3)) - (this%globalDimsR(3) / 2 - 1)

          ! Calculate wave vector and its self-product.
          if (.not. this%periodic) then
             waveX = 2d0 * PI / this%boxLength(1) * lgx
             waveY = 2d0 * PI / this%boxLength(2) * lgy
             waveZ = 2d0 * PI / this%boxLength(3) * lgz
             this%waveVectors(1, igk) = waveX
             this%waveVectors(2, igk) = waveY
             this%waveVectors(3, igk) = waveZ
             this%waveVectors2(igk) = waveX**2 + waveY**2 + waveZ**2
             ! write(0,'(a,2i5,3f12.5)'), 'orthog:  ', igk, lgx, this%waveVectors(:,igk)
          else
             this%waveVectors(:,igk) = &
                  2d0 * PI * this%unitCellVectorsK(1, :) * lgx &
                  + 2d0 * PI * this%unitCellVectorsK(2, :) * lgy &
                  + 2d0 * PI * this%unitCellVectorsK(3, :) * lgz
             this%waveVectors2(igk) = &
                  dot_product(this%waveVectors(:, igk), this%waveVectors(:, igk))
             ! write(0,'(a,2i5,3f12.5)'), 'northog: ', igk, lgx, this%waveVectors(:,igk)
          end if
       end do
    end do
    
#endif /*defined(MPI)*/

    ! Ceate an index array that maps the box wave vector index to an
    ! index in the wave vector self-product.
    call indexArray(this%waveVectors2, waveVectorToWaveVector2Map, this%totalLocalPointsK / 2)

    ! Create an array of unique, increasing wave vector magnitudes by
    ! taking the square root of the wave vector self-products.
    this%waveNumbers => safemem_realloc(this%waveNumbers, 1, .false.)
    largestWaveVector2 = -1d0
    this%waveNumberArraySize = 0
    do igk = 1, this%totalLocalPointsK / 2
       igs = waveVectorToWaveVector2Map(igk)
       waveVector2 = this%waveVectors2(igs)
       if (waveVector2 > largestWaveVector2)  then
          this%waveNumberArraySize = this%waveNumberArraySize + 1
          ! Allocate more memory than immediately needed to avoid
          ! frequent memory allocations.  We will reduce the size
          ! later if it is too big.
          if (this%waveNumberArraySize > size(this%waveNumbers)) then
             this%waveNumbers => safemem_realloc( &
                  this%waveNumbers, 4 * this%waveNumberArraySize / 3, .true.)
          endif
          this%waveNumbers(this%waveNumberArraySize) = sqrt(waveVector2)
          largestWaveVector2 = waveVector2
       endif
       this%waveVectorWaveNumberMap(igs) = this%waveNumberArraySize
    end do

    ! Reduce allocated memory to what is required.
    this%waveNumbers => safemem_realloc(this%waveNumbers, this%waveNumberArraySize, .true.)
    ierr = safemem_dealloc(waveVectorToWaveVector2Map)
  end subroutine setup_wavevector

  
  !> Convert Cartesian coordinates to fractional coordinates
  !! (or, optionally, skew coordinates) of the unit cell.
  !! Source: Wikipedia article on fractional coordinates, accessed 2014-12-06.
  !!TODO: This should use projection, not explicit tranform matrix as
  !! currently done.
  function fractionalCoordFromCartesianCoord(this, cartesianCoord) result(fractionalCoord)
    implicit none
    type(rism3d_grid), intent(in) :: this
    _REAL_, intent(in) :: cartesianCoord(3)
    _REAL_ :: fractionalCoord(3)
    _REAL_ :: unitVolume
    
    unitVolume = sqrt(1 - &
         cos(this%unitCellAngles(1))**2 - cos(this%unitCellAngles(2))**2 - cos(this%unitCellAngles(3))**3 &
         + 2 * cos(this%unitCellAngles(1)) * cos(this%unitCellAngles(2)) * cos(this%unitCellAngles(3)))

    fractionalCoord(1) = cartesianCoord(1) / this%unitCellLengths(1) &
         - cartesianCoord(2) * cos(this%unitCellAngles(3)) / (this%unitCellLengths(1) * sin(this%unitCellAngles(3))) &
         + cartesianCoord(3) * (cos(this%unitCellAngles(1)) * cos(this%unitCellAngles(3)) - cos(this%unitCellAngles(2))) &
         / (this%unitCellLengths(1) * unitVolume * sin(this%unitCellAngles(3)))
    fractionalCoord(2) = cartesianCoord(2) / (this%unitCellLengths(2) * sin(this%unitCellAngles(3))) &
         + cartesianCoord(3) * (cos(this%unitCellAngles(2)) * cos(this%unitCellAngles(3)) - cos(this%unitCellAngles(1))) &
         / (this%unitCellLengths(2) * unitVolume * sin(this%unitCellAngles(3)))
    fractionalCoord(3) = cartesianCoord(3) * sin(this%unitCellAngles(3)) / (this%unitCellLengths(3) * unitVolume)
  end function fractionalCoordFromCartesianCoord
  
end module rism3d_grid_c
