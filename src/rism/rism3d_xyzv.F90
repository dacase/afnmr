!<compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2011-2012 by
! Tyler Luchko and David A. Case.
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

!> Very basic text-based volumetric data format.  Each line is a
!! data/grid point and consists of x, y and z coordinates followed by
!! the value:
!! 
!! [x] [y] [z] [v]
!! 
!! All number are white space separated and there is no order to the
!! grid points.
module rism3d_xyzv
  use safemem
  use rism_report_c
  implicit none

contains

  !> Write volumetric data to a xyzv file.  When writing in parallel,
  !! each process must call this function with its local data. Data
  !! transfer is handled internally.
  !! @param[in] file File name to write to.
  !! @param[in] data Data to write in a n(1)*n(2)*n(3) linear array.
  !! @param[in] grid Grid object.
  !! @param[in] solute Solute object.
  !! @param[in] o_rank (optional) MPI process rank.
  !! @param[in] o_nproc (optional) MPI number of processes.
  !! @param[in] o_comm (optional) MPI communicator.
  subroutine rism3d_xyzv_write (file, data, grid, solute, o_rank, o_nproc, o_comm) 
    use rism_util, only : freeUnit
    use rism3d_grid_c
    use rism3d_solute_c
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    character(len=*), intent(in) :: file
    type(rism3d_grid), intent(in) :: grid
    _REAL_, target, intent(in) :: data(grid%localDimsR(1), grid%localDimsR(2), grid%localDimsR(3))
    type(rism3d_solute), intent(in) :: solute
    integer, optional :: o_rank, o_nproc, o_comm

    _REAL_ :: origin(3)
    
    integer :: rank = 0, nproc = 1, comm = 0
    integer :: i, j, k, irank, kOffset
    integer :: unit, iostat
    _REAL_ :: x, y, z
    _REAL_, pointer :: dataWRK(:,:,:) => NULL()

    _REAL_ :: rx(3), ry(3), rz(3), gridPoint(3)
    
#if defined(MPI)
    integer nWRK(3), status(MPI_STATUS_SIZE), err
#endif /*defined(MPI)*/

#ifdef RISM_DEBUG
    write(0,*) "xyzv_write", rank
    call flush(0)
#endif /*RISM_DEBUG*/
    
    unit = freeUnit()
    if (present(o_rank)) rank = o_rank
    if (present(o_nproc)) nproc = o_nproc
    if (present(o_comm)) comm = o_comm
    if (rank == 0) then
       open(unit = unit, file = file, iostat = iostat, status = 'REPLACE', form = 'FORMATTED')
       if (iostat /= 0) then
          call rism_report_error("could not open " // trim(file))
       end if
    end if
    kOffset = 0
    
#if defined(MPI)
    ! Data order doesn't matter for this file type. Each process will
    ! in turn send data to the master node that will then write to
    ! file. The master node will not need to transfer data, so this is
    ! a special case.
    do irank = 0, nproc - 1
       if (irank == 0) then
          dataWRK => data
       else if (rank == irank .or. rank == 0) then
          ! Get array size from process.
          if (rank == irank) &
               call mpi_send(grid%localDimsR, 3, MPI_INTEGER, 0, 0, comm, err)
          if (rank == 0) &
               call mpi_recv(nWRK, 3, MPI_INTEGER, irank, 0, comm, status, err)
          if (err /= 0) call rism_report_error &
               ("RISM3D_XYZV: could not send array size")
          ! Allocate memory.
          if (rank == 0) &
               dataWRK => safemem_realloc(dataWRK, nWRK(1), nWRK(2), nWRK(3), &
               .false.)
          kOffset = kOffset + nWRK(3)
          ! Transfer data from process.
          if (rank == irank) &
               call mpi_send(data, grid%totalLocalPointsR, MPI_DOUBLE_PRECISION, 0, 1, comm, err)
          if (rank == 0) &
               call mpi_recv(dataWRK, product(nWRK), MPI_DOUBLE_PRECISION, irank, 1, comm, status, err)
          if (err /= 0) call rism_report_error &
               ("RISM3D_XYZV: could not send array")
       end if
#else 
       dataWRK => data
#endif /*defined(MPI)*/

       if (rank == 0) then
          if (grid%periodic) then
             do k = 1, ubound(dataWRK, 3)
                rz = (k - 1) * grid%voxelVectorsR(3, :)
                do j = 1, ubound(dataWRK, 2)
                   ry = (j - 1) * grid%voxelVectorsR(2, :)
                   do i = 1, ubound(dataWRK, 1)
                      rx = (i - 1) * grid%voxelVectorsR(1, :)
                      gridPoint = rx + ry + rz
                      write(unit, "(1p, 3(E16.8E3, 1x), E16.8E3)") &
                           gridPoint, dataWRK(i, j, k)
                   end do
                end do
             end do
          else
             origin = solute%centerOfMass - grid%boxLength / 2
             do k = 1, ubound(dataWRK, 3)
                z = grid%spacing(3) * (k - 1 + kOffset) + origin(3)
                do j = 1, ubound(dataWRK, 2)
                   y = grid%spacing(2) * (j - 1) + origin(2)
                   do i = 1, ubound(dataWRK, 1)
                      x = grid%spacing(1) * (i - 1) + origin(1)
                      write(unit, "(1p, 3(E16.8E3, 1x), E16.8E3)") &
                           x, y, z, dataWRK(i, j, k)
                   end do
                end do
             end do
          end if
       end if

#if defined(MPI)
       if (irank == 0) &
            nullify(dataWRK)
    end do
    if (safemem_dealloc(dataWRK) /= 0) &
         call rism_report_error("Could not deallocate dataWRK")
#endif /*defined(MPI)*/
    
    close(unit)
  end subroutine rism3d_xyzv_write

end module rism3d_xyzv
