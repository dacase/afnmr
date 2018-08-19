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

!> Minimal support for binary format CCP4 volumetric data output.
!! The format is described here: http://www.ccp4.ac.uk/html/maplib.html#description
!! It involves a structured file header of 256 longwords, then
!! symmetry information, then the map stored a 3-dimensional array.
!! The header itself is broken into 56 longwords followed by ten
!! 80-character text labels.
!! Supports triclinic unit cells.
!! TODO: There is parallel support for writing but this is done through
!! communicating with the master node so it is very slow.
module rism3d_ccp4
  use safemem
  use rism_report_c
  implicit none

contains

  !> Write volumetric data to a file in CCP4 format.  When writing in
  !! parallel, each process must call this function with its local
  !! data. Data transfer is handled internally.  We assume
  !! decomposition in the z-axis.
  !! @param[in] file File name to write to.
  !! @param[in] data Data to write in a n(1)*n(2)*n(3) linear array.
  !! @param[in] grid Grid object.
  !! @param[in] solute Solute object.
  !! @param[in] o_rank (optional) MPI process rank.
  !! @param[in] o_nproc (optional) MPI number of processes.
  !! @param[in] o_comm (optional) MPI communicator.
  subroutine rism3d_ccp4_map_write (file, data, grid, solute, o_rank, o_nproc, o_comm)
    use constants, only: PI
    use rism_util, only : freeUnit, rmExPrec
    use rism3d_grid_c
    use rism3d_solute_c
    implicit none
#if defined(MPI)
    include 'mpif.h' 
#endif /*defined(MPI)*/
    character(len=*), intent(in) :: file
    type(rism3d_grid), intent(in) :: grid
    _REAL_, target, intent(in) :: data(grid%localDimsR(1), grid%localDimsR(2), grid%localDimsR(3)) !, centerOfMass(3)
    type(rism3d_solute), intent(in) :: solute
    integer, optional :: o_rank, o_nproc, o_comm
    
    integer :: rank = 0, nproc = 1, comm = 0
    integer :: i,j,k, irank, err, count, icount
    integer, parameter :: dataperline = 3
    integer :: unit, iostat

#ifdef MPI
    _REAL_, pointer :: z_data(:) => NULL()
    integer, pointer :: nz_offset(:) => NULL(), nz_local(:) => NULL()
#endif /*MPI*/

    integer :: id
    _REAL_ :: minValue, maxValue, meanValue, rmsd!, totalValue
    logical, parameter :: bigEndian = ichar(transfer(1,'a')) == 0
    ! Up to 80-character long label describing file origin.
    character(len=*), parameter :: amberLabel = 'Amber 3D-RISM CCP4 map volumetric data. Format revision A.'
    
#ifdef RISM_DEBUG
    write(0, *) "writeCCP4", rank
    call flush(6)
#endif /*RISM_DEBUG*/
    
    unit = freeUnit()
    if (present(o_rank)) rank = o_rank
    if (present(o_nproc)) nproc = o_nproc
    if (present(o_comm)) comm = o_comm
! #ifdef MPI
!     nz_offset => safemem_realloc(nz_offset, nproc, .false.)
!     nz_local => safemem_realloc(nz_local, nproc, .false.)
! #endif /*MPI*/
    if (rank == 0) then
       ! Unfortunately gfortran does not support form='BINARY' for
       ! compiler-independent binary output, but that probably won't
       ! be an issue for ccp4 readers.
       ! if (gfortran) then
       open(unit=unit, file=file, iostat=iostat, access='stream', status='replace', form='unformatted')
       ! else ! Intel Fortran
       !     open(unit=unit, file=file, iostat=iostat, access='sequential', status='replace', form='binary')
       ! end if
       if (iostat /= 0) then
          call rism_report_error("opening "//trim(file))
       end if
    end if
    if (rank == 0) then
       ! Write header.

       ! Number of columns, rows, and sections (fastest to slowest changing).
       write(unit) int(grid%localDimsR, 4)
       ! Since values are stored as reals, mode == 2.
       write(unit) int(2, 4)
       ! There is no offset for column, row, or section.
       write(unit) int((/ 0, 0, 0 /), 4)
       ! Number of intervals along X, Y, Z.
       write(unit) int(grid%localDimsR, 4)
       ! Cell dimensions (Angstroms).
       write(unit) real(grid%boxLength, 4)
       ! Cell angles (degrees).
       write(unit) real(grid%unitCellAngles * 180 / PI, 4)
       ! Map column, rows, sects to X, Y, Z (1, 2, 3).
       write(unit) int((/ 1, 2, 3 /), 4)
       ! Minimum, maximum, and mean density values.
       !TODO: This should be MPI-ified.
       minValue = minval(data)
       maxValue = maxval(data)
       meanValue = sum(data) / size(data)
       ! minValue = HUGE(1d0)
       ! do k = 1, grid%localDimsR(3)
       !    do j = 1, grid%localDimsR(2)
       !       do i = 1, grid%localDimsR(1)
       !          if (data(i, j, k) < minValue) minValue = data(i, j, k)
       !          if (data(i, j, k) > maxValue) maxValue = data(i, j, k)
       !          ! totalValue = totalValue + data(i, j, k)
       !          rmsd = rmsd + (meanValue - data(i, j, k))**2
       !       end do
       !    end do
       ! end do
       rmsd = sqrt(sum((meanValue - data)**2) / size(data))
       ! rmsd = sqrt(rmsd / size(data))
       ! meanValue = totalValue / size(data)
       write(unit) real(minValue, 4), real(maxValue, 4), real(meanValue, 4)
       ! Space group number.  We assume P 1.
       write(unit) int(1, 4)
       ! Number of bytes used for storing symmetry operators.
       ! In our case, none.
       write(unit) int(0, 4)
       ! Flag for skew transformation, where 0 = none, 1 = foll.  Skew
       ! transformation is from standard orthogonal coordinate frame
       ! (as used for atoms) to orthogonal map frame, as:
       !          Xo(map) = S * (Xo(atoms) - t).
       write(unit) int(0, 4)
       ! Skew matrix S (in order S11, S12, S13, S21, etc.) if above
       ! flag is non-zero.
       write(unit) int((/ 0, 0, 0 /), 4)
       write(unit) int((/ 0, 0, 0 /), 4)
       write(unit) int((/ 0, 0, 0 /), 4)
       ! Skew translation t if above flag is non-zero.
       write(unit) int((/ 0, 0, 0 /), 4)
       ! Future use (zero default).
       do id = 1, 15
          write(unit) int(0, 4)
       end do
       ! Character string 'MAP ' to identify file type.
       write(unit) 'MAP '
       ! Machine stamp indicating endianness.
       if (bigEndian) then
          write(unit) real(z'11110000', 4)
       else
          write(unit) real(z'00004144', 4)
       end if
       ! RMS deviation of map from mean density.
       write(unit) real(rmsd, 4)
       ! Number of labels being used.
       write(unit) int(1, 4)
       ! Ten 80-character labels.
       write(unit) amberLabel
       do id = 1, (9 * 80) + (80 - len(amberLabel))
          write(unit) int(0, 1)
       end do

       ! Symmetry records would go here, but we are not using any.
       
       ! Write volumetric data.
       !TODO: This needs to be MPI-ified
       write(unit) real(data, 4)
    end if
!     count = 0
! #if defined(MPI)
!     z_data => safemem_realloc(z_data, grid%globalDimsR(3), .false.)
!     call mpi_gather(grid%localDimsR(3), 1, mpi_integer, nz_local, 1, mpi_integer, 0, comm, err)
!     if (err /= 0) call rism_report_error&
!          ("RISM3D_CCP4: could not gather N")
!     nz_offset(1) = 0
!     do i = 2, nproc
!        nz_offset(i) = sum(nz_local(1:i-1))
!     end do
! #endif /*defined(MPI)*/
!     do i = 1, grid%localDimsR(1)
!        do j = 1, grid%localDimsR(2)
! ! #if defined(MPI)
! !           call mpi_gatherv(data(i, j, :), grid%localDimsR(3), mpi_double_precision, &
! !                z_data, nz_local, nz_offset, mpi_double_precision, &
! !                0, comm, err)
! !           if (err /= 0) call rism_report_error&
! !                ("RISM3D_CCP4: could not gather DATA")
! ! #endif /*defined(MPI)*/
!           if (rank == 0) then
!              do k = 1, grid%globalDimsR(3)
! ! #if defined(MPI)
! !                 write(unit, "(E16.5E3)", advance = "no") rmExPrec(z_data(k))
! ! #else
! !                 write(unit, "(E16.5E3)", advance = "no") rmExPrec(data(i, j, k))
! ! #endif /*defined(MPI)*/
!                 count = count + 1
!                 if (mod(count, dataperline) == 0) then
!                    ! write(unit, "()")
!                    count = 0;
!                 end if
!              end do
!           end if
!        end do
!     end do
!     if (rank == 0) then
!        if (mod(count, dataperline) /= 0) then
!           ! write(unit, "()")
!        end if
!        ! write(unit,"(a)") 'object "Untitled" class field'
!        close(unit)
!     end if
! #ifdef MPI
!     if (safemem_dealloc(nz_offset) /= 0) call rism_report_error("WRITEDX: nz_offset deallocation failed")
!     if (safemem_dealloc(nz_local) /= 0) call rism_report_error("WRITEDX: nz_local deallocation failed")
!     if (safemem_dealloc(z_data) /= 0) call rism_report_error("WRITEDX: z_data deallocation failed")
! #endif /*MPI*/
  end subroutine rism3d_ccp4_map_write
  
end module rism3d_ccp4
