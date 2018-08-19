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

!> Closure super class for 3D-RISM.  Closure sub-classes (i.e.,
!! actual closure implementations) are registered here.  Subroutine
!! calls then call the appropriate subroutine of the subclass
!! interface.
!!
!! This is an explicit implementation of class inheritance. See
!! V. K. Decyk, C. D. Norton, B. K. Szymanski.  How to express C++
!! concepts in Fortran 90. Scientific Programming. 6, 363-390 (1997).
!!
!! Some closure independent properties are calculated within this
!! class. Uvv is the site-site potential.
!!
!! In gerneral, this class is MPI aware only through the rism3d_grid
!! class (it knows about the total system size and its own piece of
!! it).  It does not know about processes and does not perform
!! reductions (all thermodynamic quantities are distributed and each
!! process has only the contribution of the local slab).
module rism3d_closure_c
  use rism3d_potential_c
  use rism3d_grid_c
  use rism3d_kh_c
  use rism3d_hnc_c
  use rism3d_psen_c
  use rism_report_c
  use safemem
  implicit none

  type rism3d_closure
     !> Currenly selected closure.
     character(len=4) :: type
     !> Kovalenko-Hirata closure.
     type(rism3d_kh), pointer :: kh => NULL()
     !> Hypernetted chain equation closure.
     type(rism3d_hnc), pointer :: hnc => NULL()
     !> Partial series expansion of order n (PSE-n) closure.
     type(rism3d_psen), pointer :: psen => NULL()
     !> Electric potential object.
     type(rism3d_potential), pointer :: potential => NULL()
     !> Box grid, stored in potential object.
     type(rism3d_grid), pointer :: grid => NULL()
     !> Solvent, stored in potential object.
     type(rism3d_solvent), pointer :: solvent => NULL()
     !> Solute, stored in potential object.
     type(rism3d_solute), pointer :: solute => NULL()
  end type rism3d_closure

contains

  !> Creates a new closure object of the requested type.
  !! @param[in,out] this the closure object
  !! @param[in] pot rism3d_potential object.  Must be initialized.
  !! @param[in] type one of 'KH', 'HNC', 'PSEn', where 'n' is the
  !!          order of the PSE-n closure
  subroutine rism3d_closure_new(this,type,pot)
    use rism_util, only : caseup
    implicit none
    type(rism3d_closure), intent(inout) :: this
    type(rism3d_potential), target, intent(in)  :: pot
    character(len=*), intent(in) :: type
    integer :: order, iostat
    this%potential => pot
    this%grid => this%potential%grid
    this%solvent => this%potential%solvent
    this%solute => this%potential%solute
    this%type = trim(type)
    call caseup(this%type)
    if (this%type .eq. "KH") then
       allocate(this%kh)
       call rism3d_kh_new(this%kh,this%potential)
    else if (index(this%type,"PSE") == 1) then
       read(this%type(4:),*, iostat=iostat) order
       if (iostat/=0)&
          call rism_report_error("'"//trim(this%type)//"' not a valid closure")
       allocate(this%psen)
       call rism3d_psen_new(this%psen,this%potential,order)
    else if (trim(this%type) .eq. "HNC") then
       allocate(this%hnc)
       call rism3d_hnc_new(this%hnc,this%potential)
    else
       call rism_report_error("'"//trim(this%type)//"' not a valid closure")
    end if
  end subroutine rism3d_closure_new


  !> Check if we can perform a temperature derivative calculation
  !! (i.e. does the closure support temperature derivatives).
  !! @param[in] this Closure object.
  !! @return True if we can, false if we can't.
  function rism3d_closure_canCalc_DT(this) result(can_dT)
    implicit none
    type(rism3d_closure), intent(in) :: this
    logical :: can_dT
    can_dT=.false.
    if (associated(this%kh)) then
       can_dT=.true.
    else if (associated(this%psen)) then
       can_dT=.true.
    else if (associated(this%hnc)) then
       can_dT=.true.
    end if
  end function rism3d_closure_canCalc_DT


  !> Returns a identifier string for the closure type.
  !! @param[in] this The closure object.
  function rism3d_closure_type(this) result(type)
    implicit none
    type(rism3d_closure), intent(in) :: this
    character(len=4) :: type
    type=this%type
  end function rism3d_closure_type


  !> Calculates Guv from Uuv, Huv, and Cuv using the associated
  !! closure.
  !! @param[in] this The closure object.
  !! @param[in] guv Site-site pair correlation function.
  !! @param[in] huv Site-site total correlation function.
  !! @param[in] cuv Site-site direct correlation function.
  subroutine rism3d_closure_guv(this,guv, huv, cuv)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(out) :: guv(:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    if (associated(this%kh)) then
       call rism3d_kh_guv(this%kh,guv,huv,cuv)
    else if (associated(this%psen)) then
       call rism3d_psen_guv(this%psen,guv,huv,cuv)
    else if (associated(this%hnc)) then
       call rism3d_hnc_guv(this%hnc,guv,huv,cuv)
    end if
  end subroutine rism3d_closure_guv


  !> Calculates Guv temperature derivative (T*d/dT) from Uuv, Huv_dT,
  !! Cuv_dT, Huv and Cuv using the associated closure
  !! @param[in,out] this The closure object.
  !! @param[out] guv_dT (T*d/dT) site-site pair correlation function.
  !! @param[in] huv_dT (T*d/dT) site-site total correlation function.
  !! @param[in] cuv_dT (T*d/dT) site-site direct correlation function.
  !! @param[in] guv Site-site pair correlation function.
  !! @param[in] huv Site-site total correlation function.
  !! @param[in] cuv Site-site direct correlation function.
  subroutine rism3d_closure_guv_dt(this, guv_dT, huv_dT, cuv_dT, guv, huv, cuv)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(out) :: guv_dT(:,:)
    _REAL_, intent(in) :: huv_dT(:,:), cuv_dT(:,:,:,:), guv(:,:), huv(:,:),&
         cuv(:,:,:,:)
    if (associated(this%kh)) then
       call rism3d_kh_guv_dt(this%kh,guv_dT,huv_dT,cuv_dT,guv)
    else if (associated(this%psen)) then
       call rism3d_psen_guv_dT(this%psen,guv_dT,huv_dT,cuv_dT,guv, huv, cuv)
    else if (associated(this%hnc)) then
       call rism3d_hnc_guv_dT(this%hnc,guv_dT,huv_dT,cuv_dT,guv)
    end if
  end subroutine rism3d_closure_guv_dt


  !> Calculates the excess chemical potential in kT for each site.
  !! @param[in,out] this The closure object.
  !! @param[in] huv Site-site total correlation function.
  !! @param[in] cuv Site-site direct correlation function.
  function rism3d_closure_excessChemicalPotential(this, huv, cuv) &
           result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:), cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential(this%solvent%numAtomTypes)
    if (associated(this%kh)) then
       excessChemicalPotential = rism3d_kh_excessChemicalPotential(this%kh,huv,cuv)
    else if (associated(this%psen)) then
       excessChemicalPotential = rism3d_psen_excessChemicalPotential(this%psen,huv,cuv)
    else if (associated(this%hnc)) then
       excessChemicalPotential = rism3d_hnc_excessChemicalPotential(this%hnc,huv,cuv)
    end if
  end function rism3d_closure_excessChemicalPotential


  !!Returns a 3D map of the excess chemical potential in kT.
  !!Integrating this map gives the total excess chemical potential as
  !!returned from rism3d_closure_excessChemicalPotential(,.false.).  Memory is allocated into a
  !!pointer and must be freed by the calling function. For MPI, only
  !!the grid points local to this process are allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!   a 3D-grid of the excess chemical potential contributions
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_closure_excessChemicalPotential_tot_map(this, huv, cuv) &
          result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_,pointer :: excessChemicalPotential(:,:,:)
    integer :: ix, iy, iz, iv
    nullify(excessChemicalPotential)
    excessChemicalPotential => safemem_realloc(excessChemicalPotential, &
         this%grid%localDimsR(1), &
         this%grid%localDimsR(2), &
         this%grid%localDimsR(3))
    excessChemicalPotential = 0d0
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
             excessChemicalPotential(ix,iy,iz) = sum(&
                  rism3d_closure_excessChemicalPotential_IJK&
                  (this,huv,cuv,(/ix,iy,iz/)))
          end do
       end do
    end do
  end function rism3d_closure_excessChemicalPotential_tot_map

  !!Returns a 3D map the excess chemical potential in kT for a specific grid point.
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!   ijk  : 3d-grid index
  !!OUT:
  !!   an array of excess chemical potentials at i,j,k. One for each solvent site.
  function rism3d_closure_excessChemicalPotential_IJK(this, huv, cuv,ijk) &
          result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential(this%potential%solvent%numAtomTypes)
    integer, intent(in) :: ijk(3)
    integer :: iv
    
    if (associated(this%kh)) then
       do iv=1,this%potential%solvent%numAtomTypes
          excessChemicalPotential(iv) = &
               rism3d_kh_excessChemicalPotential_IJK&
               (this%kh,huv,cuv,ijk,iv)
       end do
    else if (associated(this%psen)) then
       do iv=1,this%potential%solvent%numAtomTypes
          excessChemicalPotential(iv) = &
               rism3d_psen_excessChemicalPotential_IJK&
               (this%psen,huv,cuv,ijk,iv)
       end do
    else if (associated(this%hnc)) then
       do iv=1,this%potential%solvent%numAtomTypes
          excessChemicalPotential(iv) = &
               rism3d_hnc_excessChemicalPotential_IJK&
               (this%hnc,huv,cuv,ijk,iv)
       end do
    end if
  end function rism3d_closure_excessChemicalPotential_IJK


  !!Returns a 3D map of the excess chemical potential in kT with the
  !!Gaussian fluctuation correction.  Integrating this map gives the
  !!total excess chemical potential as returned from
  !!rism3d_closure_excessChemicalPotential(,.false.).  Memory is allocated into a
  !!pointer and must be freed by the calling function. For MPI, only
  !!the grid points local to this process are allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!   a 3D-grid of the excess chemical potential contributions
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_closure_excessChemicalPotentialGF_tot_map(this, huv, cuv) &
           result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_,pointer :: excessChemicalPotential(:,:,:)
    integer :: ix, iy, iz, iv
    nullify(excessChemicalPotential)
    excessChemicalPotential => safemem_realloc(excessChemicalPotential, &
         this%grid%localDimsR(1), &
         this%grid%localDimsR(2), &
         this%grid%localDimsR(3))
    excessChemicalPotential=0d0
    do iv=1,this%potential%solvent%numAtomTypes
       do iz = 1, this%grid%localDimsR(3)
          do iy = 1, this%grid%localDimsR(2)
             do ix = 1, this%grid%localDimsR(1)
                excessChemicalPotential(ix,iy,iz) = excessChemicalPotential(ix,iy,iz)+&
                     rism3d_closure_excessChemicalPotentialGF_ijk&
                     (this,huv,cuv,(/ix,iy,iz/),iv)
             end do
          end do
       end do
    end do
  end function rism3d_closure_excessChemicalPotentialGF_tot_map

  !!Returns a 3D map of the excess chemical potential in kT from the
  !!PC+/3D-RISM Correction.  Integrating this map gives the total excess
  !!chemical potential as returned from
  !!rism3d_closure_excessChemicalPotential(,.false.).  Memory is allocated into a
  !!pointer and must be freed by the calling function. For MPI, only
  !!the grid points local to this process are allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!   a 3D-grid of the excess chemical potential contributions
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_closure_excessChemicalPotentialPCPLUS_tot_map(this, huv, cuv)&
           result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_,pointer :: excessChemicalPotential(:,:,:)
    integer :: ix, iy, iz, iv
    nullify(excessChemicalPotential)
    excessChemicalPotential => safemem_realloc(excessChemicalPotential, &
         this%grid%localDimsR(1), &
         this%grid%localDimsR(2), &
         this%grid%localDimsR(3))
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
             excessChemicalPotential(ix,iy,iz) = &
                  rism3d_closure_excessChemicalPotentialPCPLUS_ijk&
                  (this,huv,cuv,(/ix,iy,iz/))
          end do
       end do
    end do
  end function rism3d_closure_excessChemicalPotentialPCPLUS_tot_map
  
  !!Returns a 3D map of the excess chemical potential in kT from the
  !!Palmer correction.  Integrating this map gives the total excess
  !!chemical potential as returned from
  !!rism3d_closure_excessChemicalPotential(,.false.).  Memory is allocated into a
  !!pointer and must be freed by the calling function. For MPI, only
  !!the grid points local to this process are allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!   coeff: coefficients for correction.  For the original correction
  !!          (a,b)=coeff(1:2). Extra coefficients are ignored.
  !!OUT:
  !!   a 3D-grid of the excess chemical potential contributions
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_closure_excessChemicalPotentialUC_tot_map(this, huv, cuv,coeff)&
       result(excessChemicalPotential)
    use rism_util, only : caseup
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_, intent(in) :: coeff(:)
    _REAL_,pointer :: excessChemicalPotential(:,:,:)
    integer :: ix, iy, iz, iv
    nullify(excessChemicalPotential)
    excessChemicalPotential => safemem_realloc(excessChemicalPotential, &
         this%grid%localDimsR(1), &
         this%grid%localDimsR(2), &
         this%grid%localDimsR(3))
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
             excessChemicalPotential(ix,iy,iz) = &
                  rism3d_closure_excessChemicalPotentialUC_ijk&
                  (this,huv,cuv,coeff,(/ix,iy,iz/))
          end do
       end do
    end do

  end function rism3d_closure_excessChemicalPotentialUC_tot_map


  !!Returns a 3D map of the excess chemical potential in kT.
  !!Integrating this map gives the excess chemical potential
  !!contributions from each site as returned from
  !!rism3d_closure_excessChemicalPotential(,.false.).  Memory is allocated into a
  !!pointer and must be freed by the calling function. For MPI, only
  !!the grid points local to this process are allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!   a 3D-grid of the excess chemical potential contributions (nx,ny,nz,nsite)
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_closure_excessChemicalPotential_site_map(this, huv, cuv) &
       result(excessChemicalPotential)
    use rism_util, only : caseup
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_,pointer :: excessChemicalPotential(:,:,:,:)
    integer :: ix, iy, iz, iv
    nullify(excessChemicalPotential)
    excessChemicalPotential => safemem_realloc(excessChemicalPotential, &
         this%grid%localDimsR(1), &
         this%grid%localDimsR(2), &
         this%grid%localDimsR(3), &
         this%solvent%numAtomTypes)
    excessChemicalPotential=0d0
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
             excessChemicalPotential(ix,iy,iz,:) = &
                  rism3d_closure_excessChemicalPotential_IJK&
                  (this,huv,cuv,(/ix,iy,iz/))
          end do
       end do
    end do
  end function rism3d_closure_excessChemicalPotential_site_map


  !> Returns a 3D map of the excess chemical potential in kT from the
  !! GF approximation.  Integrating this map gives the excess chemical
  !! potential contributions from each site as returned from
  !! rism3d_closure_excessChemicalPotential(,.false.).  Memory is allocated into a
  !! pointer and must be freed by the calling function. For MPI, only
  !! the grid points local to this process are allocated and
  !! calculated.
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!   a 3D-grid of the excess chemical potential contributions (nx,ny,nz,nsite)
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_closure_excessChemicalPotentialGF_site_map(this, huv, cuv) &
       result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_,pointer :: excessChemicalPotential(:,:,:,:)
    integer :: ix, iy, iz, iv
    nullify(excessChemicalPotential)
    excessChemicalPotential => safemem_realloc(excessChemicalPotential, &
         this%grid%localDimsR(1), &
         this%grid%localDimsR(2), &
         this%grid%localDimsR(3), &
         this%solvent%numAtomTypes)
    excessChemicalPotential=0d0
    do iv=1,this%potential%solvent%numAtomTypes
       do iz = 1, this%grid%localDimsR(3)
          do iy = 1, this%grid%localDimsR(2)
             do ix = 1, this%grid%localDimsR(1)
                excessChemicalPotential(ix,iy,iz,iv) = &
                     rism3d_closure_excessChemicalPotentialGF_ijk&
                     (this,huv,cuv,(/ix,iy,iz/),iv)
             end do
          end do
       end do
    end do
  end function rism3d_closure_excessChemicalPotentialGF_site_map


  !!Calculates the excess chemical potential w/ asymptotic correction in kT for
  !!each site
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  function rism3d_closure_aexcessChemicalPotential(this, huv, cuv) result(excessChemicalPotential)
    use rism_util, only : gaussquad_legendre, checksum
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential(this%solvent%numAtomTypes)
    if (associated(this%kh)) then
       excessChemicalPotential = rism3d_kh_aexcessChemicalPotential(this%kh,huv,cuv)
    else if (associated(this%psen)) then
       excessChemicalPotential = rism3d_psen_aexcessChemicalPotential(this%psen,huv,cuv)
    else if (associated(this%hnc)) then
       excessChemicalPotential = rism3d_hnc_aexcessChemicalPotential(this%hnc,huv,cuv)
    end if
  end function rism3d_closure_aexcessChemicalPotential


  !!Calculates the excess chemical potential in kT for each site using the
  !!Gaussian fluctuation expression.  This is closure independent.
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  function rism3d_closure_excessChemicalPotentialGF(this, huv, cuv) result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential(this%solvent%numAtomTypes)
    _REAL_ :: tuv
    integer :: ix, iy, iz, iv, igk
    excessChemicalPotential = 0.d0
    do iv=1,this%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
#ifdef MPI
                igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) &
                     + (iz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
#else
                igk = ix + (iy - 1) * this%grid%localDimsR(1) &
                     + (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#endif
                tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                excessChemicalPotential(iv) = excessChemicalPotential(iv) - &
                     (1.d0+0.5d0*huv(igk,iv))*cuv(ix,iy,iz,iv)
             end do
          end do
       end do
    end do
    excessChemicalPotential = this%solvent%density*excessChemicalPotential*this%grid%voxelVolume
  end function rism3d_closure_excessChemicalPotentialGF

  !!Calculates the excess chemical potential in kT for each site using
  !!the PC+/3D-RISM Correction. without long
  !!range corrections. See rism3d_closure_aexcessChemicalPotentialPCPLUS().
  !!
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!    The total excess chemical potential.  It is not site decomposible.
  function rism3d_closure_excessChemicalPotentialPCPLUS(this, huv, cuv) &
       result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential
    integer :: ix, iy,iz
    excessChemicalPotential = 0d0
      excessChemicalPotential = &
           sum(rism3d_closure_excessChemicalPotential(this,huv,cuv)) &
           - 0.5d0 * rism3d_closure_partialMolarVolume(this,cuv) * &
           (1d0/this%solvent%xikt + sum(this%solvent%density_sp))
  end function rism3d_closure_excessChemicalPotentialPCPLUS

  !!Calculates the excess chemical potential in kT for each site at the
  !!requested grid point using the Gaussian fluctuation expression.
  !!This is closure independent.
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!   ijk  : 3d-grid index
  !!   iv   : solvent site
  function rism3d_closure_excessChemicalPotentialGF_IJK(this, huv, cuv,ijk,iv) &
           result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    integer, intent(in) :: ijk(3), iv
    _REAL_ :: excessChemicalPotential
    _REAL_ :: tuv
    integer :: ix, iy, iz, ig, igk
    ix=ijk(1)
    iy=ijk(2)
    iz=ijk(3)
    ig = ix + (iy-1)*this%grid%localDimsR(1) &
            + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#ifdef MPI
    igk = ix + (iy-1)*(this%grid%localDimsR(1)+2) &
             + (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
    igk = ix + (iy-1)*this%grid%localDimsR(1) &
             + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#endif
    tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
    excessChemicalPotential = - (1.d0+0.5d0*huv(igk,iv)) &
         * cuv(ix,iy,iz,iv)*this%potential%solvent%density(iv)
  end function rism3d_closure_excessChemicalPotentialGF_IJK

  !!Calculates the excess chemical potential in kT for each site using
  !!the PC+/3D-RISM Correction expression at the
  !!requested grid point.  See rism3d_closure_aexcessChemicalPotentialPCPLUS().
  !!
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!    The total excess chemical potential.  It is not site decomposible.
  function rism3d_closure_excessChemicalPotentialPCPLUS_ijk(this, huv, cuv, ijk) &
       result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    integer, intent(in) :: ijk(3)
    _REAL_ :: excessChemicalPotential
    integer :: iv 
    
    excessChemicalPotential = sum( &
            rism3d_closure_excessChemicalPotential_IJK(this,huv,cuv,ijk))

    excessChemicalPotential = excessChemicalPotential  &
         - 0.5*rism3d_closure_partialMolarVolume_IJK(this,cuv,ijk) &
         * (1d0/this%solvent%xikt + sum(this%solvent%density_sp))
  end function rism3d_closure_excessChemicalPotentialPCPLUS_ijk

  !!Calculates the excess chemical potential in kT for each site using the
  !!Gaussian fluctuation expression.  This is closure independent.
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  function rism3d_closure_aexcessChemicalPotentialGF(this, huv, cuv) &
           result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential(this%solvent%numAtomTypes)
    _REAL_ :: excessChemicalPotentiallr(this%solvent%numAtomTypes)
    _REAL_ :: excessChemicalPotentialh2lr(this%solvent%numAtomTypes)
    _REAL_ :: excessChemicalPotentialhclr(this%solvent%numAtomTypes)
    _REAL_ :: tuv, cuvlr, huvlr
    integer :: ix, iy, iz, iv, ig, igk

    if (.not.this%solute%charged) then
       excessChemicalPotential = &
           rism3d_closure_excessChemicalPotentialGF(this,huv,cuv)
       return
    end if

    !
    !Closure independent, long-range part
    !
    call rism3d_potential_int_h2_hc(this%potential, &
         excessChemicalPotentialh2lr,excessChemicalPotentialhclr)
    !
    !Short-range part
    !
    huvlr=0d0
    excessChemicalPotential = 0.d0
    do iv=1,this%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
                ig = ix + (iy-1)*this%grid%localDimsR(1) &
                        + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#ifdef MPI
                igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) &
                     + (iz - 1) * this%grid%localDimsR(2) &
                                * (this%grid%localDimsR(1) + 2)
#else
                igk = ix + (iy - 1) * this%grid%localDimsR(1) &
                     + (iz - 1) * this%grid%localDimsR(2) &
                                * this%grid%localDimsR(1)
#endif
                tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                if (this%solvent%ionic) then
                   huvlr = this%solvent%charge_sp(iv) &
                           * this%potential%tcfLongRangeAsympR(ig)
                end if
                cuvlr = this%solvent%charge(iv) &
                        * this%potential%dcfLongRangeAsympR(ig)
                excessChemicalPotential(iv) = excessChemicalPotential(iv) &
                     - (1.d0+0.5d0*huv(igk,iv))*cuv(ix,iy,iz,iv) &
                     + 0.5d0*huvlr*cuvlr
             end do
          end do
       end do
       excessChemicalPotential(iv) =  &
           excessChemicalPotential(iv)*this%grid%voxelVolume
    end do
    excessChemicalPotential = this%solvent%density &
          * (excessChemicalPotential-excessChemicalPotentialhclr/2d0)
  end function rism3d_closure_aexcessChemicalPotentialGF

  !!Calculates the excess chemical potential in kT for each site using
  !!the PC+/3D-RISM Correction expression of Sergiievskyi et al. and
  !!Misin et al. with long range corrections. This is closure
  !!independent.
  !!
  !!We compute the expression as in Johnson et al.,
  !!
  !! \Delta G=\Delta\mu^RISM - 1/2*kT * v * (1/xikt + \rho_Tot)
  !!
  !!where \mu^RISM is from the closure, v is the PMV and xikt is the
  !!isothermal compressibility
  !!
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!    The total excess chemical potential.  It is not site decomposible.
  function rism3d_closure_aexcessChemicalPotentialPCPLUS(this, huv, cuv) &
           result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential
    excessChemicalPotential = &
         sum(rism3d_closure_aexcessChemicalPotential(this,huv,cuv)) &
         - 0.5d0 * rism3d_closure_partialMolarVolume(this,cuv) &
         * (1d0/this%solvent%xikt + sum(this%solvent%density_sp))
  end function rism3d_closure_aexcessChemicalPotentialPCPLUS

  !! Calculates the solvation energy, dE = dmu + dTS, from the
  !! temperature derivative (T*d/dT) w/o asymptotic correction in kT
  !! for each site.
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv     : site-site total correlation function
  !!   cuv     : site-site direct correlation function
  !! @return Excess chemical potential with long-range correction.
  function rism3d_closure_solvationEnergy(this, huv_dT, cuv_dT, huv, cuv) &
           result(solvationEnergy)
    use rism_util, only : gaussquad_legendre, checksum
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: solvationEnergy(this%solvent%numAtomTypes)

    solvationEnergy = huge(1d0)
    if (.not.rism3d_solvent_canCalc_dT(this%solvent)) then
       return
    else if (associated(this%kh)) then
       solvationEnergy = rism3d_kh_solvationEnergy(this%kh,huv_dT,cuv_dT,huv,cuv)
    else if (associated(this%psen)) then
       solvationEnergy = rism3d_psen_solvationEnergy(this%psen,huv_dT,cuv_dT,huv,cuv)
    else if (associated(this%hnc)) then
       solvationEnergy = rism3d_hnc_solvationEnergy(this%hnc,huv_dT,cuv_dT,huv,cuv)
    end if
  end function rism3d_closure_solvationEnergy


  !!Returns a 3D map of the solvation energy in kT.
  !!Integrating this map gives the total solvation energy as
  !!returned from rism3d_closure_excessChemicalPotential(,.false.).  Memory is allocated into a
  !!pointer and must be freed by the calling function. For MPI, only
  !!the grid points local to this process are allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   cuv  : site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!OUT:
  !!   a 3D-grid of the solvation energy contributions
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid

  function rism3d_closure_solvationEnergy_tot_map(this, huv_dT, cuv_dT, huv, cuv)&
       result(solvationEnergy)
    use rism_util, only : caseup
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_,pointer :: solvationEnergy(:,:,:)
    integer :: ix, iy, iz, iv
    nullify(solvationEnergy)
    solvationEnergy => safemem_realloc(solvationEnergy, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3))
    solvationEnergy=0d0
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
             solvationEnergy(ix,iy,iz) = sum(&
                  rism3d_closure_solvationEnergy_ijk(&
                  this, huv_dT, cuv_dT,huv,cuv,(/ix,iy,iz/)))
          end do
       end do
    end do
  end function rism3d_closure_solvationEnergy_tot_map


  !!Returns the solvation energy in kT for all solvent sites at a
  !!specify grid point.
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   cuv  : site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!   ijk     : 3d-grid index
  !!OUT:
  !!   an array of solvation energies at i,j,k. One for each solvent site.
  function rism3d_closure_solvationEnergy_ijk(this, huv_dT, cuv_dT, huv, cuv, ijk)&
       result(solvationEnergy)
    use rism_util, only : caseup
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    integer, intent(in) :: ijk(3)
    _REAL_ :: solvationEnergy(this%potential%solvent%numAtomTypes)
    integer :: iv
    if (associated(this%kh)) then
       do iv=1, this%solvent%numAtomTypes
          solvationEnergy(iv)=&
               rism3d_kh_solvationEnergy_ijk&
               (this%kh,huv_dT,cuv_dT,huv,cuv,ijk,iv)
       end do
    else if (associated(this%psen)) then
       do iv=1, this%solvent%numAtomTypes
          solvationEnergy(iv) = &
               rism3d_psen_solvationEnergy_ijk&
               (this%psen,huv_dT,cuv_dT,huv,cuv,ijk,iv)
       end do
    else if (associated(this%hnc)) then
       do iv=1, this%solvent%numAtomTypes
          solvationEnergy(iv) = &
               rism3d_hnc_solvationEnergy_ijk&
               (this%hnc,huv_dT,cuv_dT,huv,cuv,ijk,iv)
       end do
    end if
  end function rism3d_closure_solvationEnergy_ijk

  !!Returns a 3D map of the solvation energy in kT with the Gaussian
  !!fluctuation correction.  Integrating this map gives the total
  !!solvation energy as returned from rism3d_closure_excessChemicalPotential(,.false.).
  !!Memory is allocated into a pointer and must be freed by the calling
  !!function. For MPI, only the grid points local to this process are
  !!allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   cuv  : site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!OUT:
  !!   a 3D-grid of the solvation energy contributions
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid

  function rism3d_closure_solvationEnergyGF_tot_map(this, huv_dT, cuv_dT, huv, cuv) &
       result(solvationEnergy)
    use rism_util, only : caseup
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_,pointer :: solvationEnergy(:,:,:)
    integer :: ix, iy, iz, iv
    nullify(solvationEnergy)
    solvationEnergy => safemem_realloc(solvationEnergy, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3))
    solvationEnergy=0d0
    do iv=1,this%potential%solvent%numAtomTypes
       do iz = 1, this%grid%localDimsR(3)
          do iy = 1, this%grid%localDimsR(2)
             do ix = 1, this%grid%localDimsR(1)
                solvationEnergy(ix,iy,iz) = solvationEnergy(ix,iy,iz)+&
                     rism3d_closure_solvationEnergyGF_ijk&
                     (this,huv_dT,cuv_dT,huv,cuv,(/ix,iy,iz/),iv)
             end do
          end do
       end do
    end do
  end function rism3d_closure_solvationEnergyGF_tot_map

  !!Returns a 3D map of the solvation energy in kT with the 
  !!PC+/3D-RISM Correction.  Integrating this map gives the total solvation
  !!energy as returned from rism3d_closure_excessChemicalPotential(,.false.).  Memory is
  !!allocated into a pointer and must be freed by the calling
  !!function. For MPI, only the grid points local to this process are
  !!allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   cuv  : site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!OUT:
  !!   a 3D-grid of the solvation energy contributions
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_closure_solvationEnergyPCPLUS_tot_map(this, huv_dT, cuv_dT, huv, cuv) &
       result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_,pointer :: solvationEnergy(:,:,:)
    integer :: ix, iy, iz, iv
    nullify(solvationEnergy)
    solvationEnergy => safemem_realloc(solvationEnergy, this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3))
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
             solvationEnergy(ix,iy,iz) = rism3d_closure_solvationEnergyPCPLUS_ijk&
                  (this,huv_dT,cuv_dT,huv,cuv,(/ix,iy,iz/))
          end do
       end do
    end do
  end function rism3d_closure_solvationEnergyPCPLUS_tot_map

  !!Returns a 3D map of the solvation energy in kT with the Palmer et
  !!al. Universal Correction.  Integrating this map gives the total
  !!solvation energy as returned from rism3d_closure_excessChemicalPotential(,.false.).
  !!Memory is allocated into a pointer and must be freed by the calling
  !!function. For MPI, only the grid points local to this process are
  !!allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   cuv  : site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!   coeff: coefficients for correction.  For the original correction
  !!          (a,b)=coeff(1:2). Extra coefficients are ignored.
  !!OUT:
  !!   a 3D-grid of the solvation energy contributions
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid

  function rism3d_closure_solvationEnergyUC_tot_map(this, huv_dT, cuv_dT, huv, cuv,&
       coeff) result(solvationEnergy)
    use rism_util, only : caseup
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_, intent(in) :: coeff(:)
    _REAL_,pointer :: solvationEnergy(:,:,:)
    integer :: ix, iy, iz, iv
    nullify(solvationEnergy)
    solvationEnergy => safemem_realloc(solvationEnergy, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3))
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
             solvationEnergy(ix,iy,iz) = rism3d_closure_solvationEnergyUC_ijk&
                  (this,huv_dT,cuv_dT,huv,cuv,coeff,(/ix,iy,iz/))
          end do
       end do
    end do
  end function rism3d_closure_solvationEnergyUC_tot_map


  !!Returns a 3D map of the solvation energy in kT.
  !!Integrating this map gives the site solvation energy contributions as
  !!returned from rism3d_closure_excessChemicalPotential(,.false.).  Memory is allocated into a
  !!pointer and must be freed by the calling function. For MPI, only
  !!the grid points local to this process are allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   cuv  : site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!   a 3D-grid of the solvation energy contributions (nx,ny,nz,nsite)
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid

  function rism3d_closure_solvationEnergy_site_map(this, huv_dT, cuv_dT, huv, cuv) result(solvationEnergy)
    use rism_util, only : caseup
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_,pointer :: solvationEnergy(:,:,:,:)
    integer :: ix, iy, iz, iv
    nullify(solvationEnergy)
    solvationEnergy => safemem_realloc(solvationEnergy, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3),&
         this%solvent%numAtomTypes)
    solvationEnergy=0d0
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
             solvationEnergy(ix,iy,iz,:) = &
                  rism3d_closure_solvationEnergy_ijk&
                  (this,huv_dT,cuv_dT,huv,cuv,(/ix,iy,iz/))
          end do
       end do
    end do
  end function rism3d_closure_solvationEnergy_site_map


  !!Returns a 3D map of the solvation energy in kT from the GF
  !!approximation.  Integrating this map gives the site solvation
  !!energy contributions as returned from
  !!rism3d_closure_excessChemicalPotential(,.false.).  Memory is allocated into a
  !!pointer and must be freed by the calling function. For MPI, only
  !!the grid points local to this process are allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   cuv  : site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!   a 3D-grid of the solvation energy contributions (nx,ny,nz,nsite)
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid

  function rism3d_closure_solvationEnergyGF_site_map(this, huv_dT, cuv_dT, huv, cuv) &
       result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_,pointer :: solvationEnergy(:,:,:,:)
    integer :: ix, iy, iz, iv
    nullify(solvationEnergy)
    solvationEnergy => safemem_realloc(solvationEnergy, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3),&
         this%solvent%numAtomTypes)
    solvationEnergy=0d0
    do iv=1,this%potential%solvent%numAtomTypes
       do iz = 1, this%grid%localDimsR(3)
          do iy = 1, this%grid%localDimsR(2)
             do ix = 1, this%grid%localDimsR(1)
                solvationEnergy(ix,iy,iz,iv) = &
                     rism3d_closure_solvationEnergyGF_ijk&
                     (this,huv_dT,cuv_dT,huv,cuv,(/ix,iy,iz/),iv)
             end do
          end do
       end do
    end do
  end function rism3d_closure_solvationEnergyGF_site_map


  !!Calculates the solvation energy, dE = dmu + dTS, from the
  !!temperature derivative (T*d/dT) w/ asymptotic correction in kT for
  !!each site
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv     : site-site total correlation function
  !!   cuv     : site-site direct correlation function
  !!OUT:
  !!    excess chemical potential with long-range correction

  function rism3d_closure_asolvationEnergy(this, huv_dT, cuv_dT, huv, cuv) result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: solvationEnergy(this%solvent%numAtomTypes)

    solvationEnergy = huge(1d0)
    if (.not.rism3d_solvent_canCalc_dT(this%solvent)) then
       return
    else if (associated(this%kh)) then
       solvationEnergy = rism3d_kh_asolvationEnergy(this%kh,huv_dT,cuv_dT,huv,cuv)
    else if (associated(this%psen)) then
       solvationEnergy = rism3d_psen_asolvationEnergy(this%psen,huv_dT,cuv_dT,huv,cuv)
    else if (associated(this%hnc)) then
       solvationEnergy = rism3d_hnc_asolvationEnergy(this%hnc,huv_dT,cuv_dT,huv,cuv)
    end if
  end function rism3d_closure_asolvationEnergy


  !!Calculates the solvation energy, dE = dmu + dTS, for each site
  !!using the Gaussian fluctuation expression temperature derivative.
  !!
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv     : site-site total correlation function
  !!   cuv     : site-site direct correlation function
  !!OUT:
  !!    excess chemical potential with long-range correction

  function rism3d_closure_solvationEnergyGF(this, huv_dT, cuv_dT, huv, cuv) result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: solvationEnergy(this%solvent%numAtomTypes)
    integer :: ix, iy, iz, iv, ig, igk

    solvationEnergy = huge(1d0)
    if (.not.rism3d_solvent_canCalc_dT(this%solvent))&
         return
    solvationEnergy =0d0
    do iv=1,this%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
#ifdef MPI
                igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) &
                     + (iz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
#else
                igk = ix + (iy - 1) * this%grid%localDimsR(1) &
                     + (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#endif
                  solvationEnergy(iv) = solvationEnergy(iv) + cuv_dT(ix,iy,iz,iv) + &
                       0.5d0 * (huv(igk,iv)*cuv_dT(ix,iy,iz,iv) + &
                       huv_dT(igk,iv)*cuv(ix,iy,iz,iv))
             end do
          end do
       end do
       solvationEnergy(iv) =  this%potential%solvent%density(iv)*&
            solvationEnergy(iv)*this%grid%voxelVolume
    end do
  end function rism3d_closure_solvationEnergyGF

  !!Calculates the excess solvation energy in kT using the PC+/3D-RISM 
  !!Correction expression with long range
  !!corrections.  See rism3d_closure_asolvationEnergyPCPLUS().
  !!
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!    The total solvation energy.  It is not site decomposible.
  function rism3d_closure_solvationEnergyPCPLUS(this, huv_dT, cuv_dT, huv, cuv) result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: solvationEnergy
    integer :: ix, iy,iz
    if(this%solvent%xikt_dT == huge(1d0)) then
       solvationEnergy = huge(1d0)
       return
    end if
    solvationEnergy = 0d0
      solvationEnergy = sum(rism3d_closure_solvationEnergy(this,huv_dT,cuv_dT,huv,cuv)) &
           + 0.5d0 * rism3d_closure_partialMolarVolume_dT(this,cuv_dT,cuv) &
           * (1d0/this%solvent%xikt + sum(this%solvent%density_sp)) &
           - 0.5d0 * rism3d_closure_partialMolarVolume(this,cuv) &
           * (1d0/this%solvent%xikt)**2 *  (this%solvent%xikt+ this%solvent%xikt_dT)
  end function rism3d_closure_solvationEnergyPCPLUS

  !!Calculates the solvation energy, dE = dmu + dTS, for each site at
  !!the requested grid point using the Gaussian fluctuation expression
  !!temperature derivative.
  !!IN:
  !!   this    : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv     : site-site total correlation function
  !!   cuv     : site-site direct correlation function
  !!   ijk     : 3d-grid index
  !!   iv      : solvent site
  !!OUT:
  !!    excess chemical potential with long-range correction

  function rism3d_closure_solvationEnergyGF_IJK(this, huv_dT, cuv_dT, huv, cuv,ijk,iv) result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    integer, intent(in) :: ijk(3), iv
    _REAL_ :: solvationEnergy
    integer :: ix, iy, iz, ig, igk
    ix=ijk(1)
    iy=ijk(2)
    iz=ijk(3)
#ifdef MPI
    igk = ix + (iy-1)*(this%grid%localDimsR(1)+2) + (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
    igk = ix + (iy-1)*this%grid%localDimsR(1) + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#endif
    solvationEnergy = cuv_dT(ix,iy,iz,iv) + &
         0.5d0 * (huv(igk,iv)*cuv_dT(ix,iy,iz,iv) + &
         huv_dT(igk,iv)*cuv(ix,iy,iz,iv))
    solvationEnergy= solvationEnergy*this%potential%solvent%density(iv)
  end function rism3d_closure_solvationEnergyGF_IJK

  !!Calculates the solvation energy in kT using the PC+/3D-RISM
  !!Correction expression at the requested grid point.
  !!
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!   ijk     : 3d-grid index
  !!OUT:
  !!    The total solvation energy.  It is not site decomposible.
  function rism3d_closure_solvationEnergyPCPLUS_IJK(this, huv_dT, cuv_dT, huv, cuv,&
       ijk) result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    integer, intent(in) :: ijk(3)
    _REAL_ :: solvationEnergy
    integer :: iv
    solvationEnergy = sum(&
            rism3d_closure_solvationEnergy_IJK(this,huv_dT,cuv_dT,huv,cuv,ijk))

    solvationEnergy = solvationEnergy &
         + 0.5d0 * rism3d_closure_partialMolarVolume_dT_IJK(this,cuv_dT,cuv,ijk) &
         * (1d0/this%solvent%xikt + sum(this%solvent%density_sp)) &
         - 0.5d0 * rism3d_closure_partialMolarVolume_IJK(this,cuv,ijk) &
         * (1d0/this%solvent%xikt)**2 *  (this%solvent%xikt+ this%solvent%xikt_dT)
  end function rism3d_closure_solvationEnergyPCPLUS_IJK

  !!Calculates the solvation energy, dE = dmu + dTS, for each site
  !!using the Gaussian fluctuation expression temperature derivative.
  !!
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv     : site-site total correlation function
  !!   cuv     : site-site direct correlation function
  !!OUT:
  !!    excess chemical potential with long-range correction
  function rism3d_closure_asolvationEnergyGF(this, huv_dT, cuv_dT, huv, cuv) result(solvationEnergy)
    use rism_util, only : gaussquad_legendre, checksum
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: solvationEnergy(this%solvent%numAtomTypes)
    _REAL_ :: solvationEnergyh2lr(this%solvent%numAtomTypes),solvationEnergyhclr(this%solvent%numAtomTypes)
    _REAL_ ::cuvlr, huvlr
    integer :: ix, iy, iz, iv, ig, igk

    solvationEnergy = huge(1d0)
    if (.not.rism3d_solvent_canCalc_dT(this%solvent))&
         return
    if (.not.this%solute%charged) then
       solvationEnergy = rism3d_closure_solvationEnergyGF(this, huv_dT, cuv_dT,huv,cuv)
       return
    end if

    ! Closure independent, long-range part.
    call rism3d_potential_int_h2_hc(this%potential,solvationEnergyh2lr,solvationEnergyhclr)

    huvlr=0d0
    solvationEnergy = 0d0
    do iv=1,this%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
                ig = ix + (iy - 1) * this%grid%localDimsR(1) &
                        + (iz - 1) * this%grid%localDimsR(2) &
                                   * this%grid%localDimsR(1)
#ifdef MPI
                igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                                    * (this%grid%localDimsR(1) + 2)
#else
                igk = ix + (iy - 1) * this%grid%localDimsR(1) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                                    * this%grid%localDimsR(1)
#endif

                if (this%solvent%ionic) then
                   huvlr = this%solvent%charge_sp(iv)*this%potential%tcfLongRangeAsympR(ig)
                end if
                cuvlr = this%solvent%charge(iv)*this%potential%dcfLongRangeAsympR(ig)
                solvationEnergy(iv) = solvationEnergy(iv) + cuv_dT(ix,iy,iz,iv) &
                     + 0.5d0 * (huv(igk,iv)*cuv_dT(ix,iy,iz,iv) &
                     + huv_dT(igk,iv)*cuv(ix,iy,iz,iv))&
                     - huvlr*cuvlr
             end do
          end do
       end do
    end do
    solvationEnergy =  solvationEnergy*this%grid%voxelVolume
    solvationEnergy =  (solvationEnergy+solvationEnergyhclr)*this%solvent%density
  end function rism3d_closure_asolvationEnergyGF

  !!Calculates the excess solvation energy in kT using the PC+/3D-RISM 
  !!Correction expression with long range
  !!corrections.  See rism3d_closure_aexcessChemicalPotentialPCPLUS().
  !!
  !!We compute the expression as in Johnson et al.,
  !!
  !!  \Delta\epsilon=\Delta\epsilon^RISM + 1/2*kT*v_dT*(1/(xikt)+\rho_Tot)
  !!                 - 1/2*kT*v(1/(xikt))^{2} * (xikt + xikt_dT)
  !!
  !!where \mu^RISM is from the closure, v_dT is the PMV temperature
  !!derivative and xikt_dT is the isothermal compressibility
  !!temperature derivative.
  !!
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!    The total excess chemical potential.  It is not site decomposible.
  function rism3d_closure_asolvationEnergyPCPLUS(this, huv_dT, cuv_dT, huv, cuv) result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: solvationEnergy

    if(this%solvent%xikt_dT == huge(1d0)) then
       solvationEnergy = huge(1d0)
       return
    end if
    solvationEnergy = sum(rism3d_closure_asolvationEnergy(this,huv_dT,cuv_dT,huv,cuv)) &
         + 0.5d0 * rism3d_closure_partialMolarVolume_dT(this,cuv_dT,cuv)&
         * (1d0/this%solvent%xikt + sum(this%solvent%density_sp)) &
         - 0.5d0 * rism3d_closure_partialMolarVolume(this,cuv)&
         * (1d0/this%solvent%xikt)**2 * (this%solvent%xikt+this%solvent%xikt_dT)

  end function rism3d_closure_asolvationEnergyPCPLUS

  !!Calculates the excess chemical potential in kT using the
  !!temperature dependent Universal Correction with out long range
  !!asymptotics. See rism3d_closure_aexcessChemicalPotentialUC() for
  !!more details.
  !!
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!   coeff: coefficients for correction: (a0,b0,a1,b1). Extra
  !!          coefficients are ignored.
  !!OUT:
  !!    The total excess chemical potential.  It is not site decomposible.
  function rism3d_closure_excessChemicalPotentialUC(this, huv, cuv, coeff) result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:),coeff(:)
    _REAL_ :: excessChemicalPotential
    integer :: ix, iy,iz
    excessChemicalPotential = 0d0
    excessChemicalPotential = sum(rism3d_closure_excessChemicalPotential(this,huv,cuv)) &
         + (coeff(1) + coeff(3)*this%solvent%temperature)&
         *rism3d_closure_partialMolarVolume(this,cuv)
    if (this%grid%mpirank==0)&
         excessChemicalPotential = excessChemicalPotential + (coeff(2) + coeff(4)*this%solvent%temperature)
  end function rism3d_closure_excessChemicalPotentialUC


  !!Calculates the excess chemical potential in kT for the temperature
  !!dependent Universal Correction expression at the requested grid
  !!point. See rism3d_closure_aexcessChemicalPotentialUC() for more
  !!details.
  !!
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!   coeff: coefficients for correction: (a0,b0,a1,b1). Extra
  !!          coefficients are ignored.
  !!   ijk     : 3d-grid index
  !!OUT:
  !!    The total excess chemical potential.  It is not site decomposible.
  function rism3d_closure_excessChemicalPotentialUC_ijk(this, huv, cuv, coeff, ijk) result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:),coeff(:)
    integer, intent(in) :: ijk(3)
    _REAL_ :: excessChemicalPotential
    integer :: iv

    excessChemicalPotential = 0d0
    ! \Delta G
    excessChemicalPotential = sum(&
         rism3d_closure_excessChemicalPotential_ijk(this,huv,cuv,ijk))
    
    ! (a0+a1*T)*V
    excessChemicalPotential = excessChemicalPotential  &
         + (coeff(1) + coeff(3)*this%solvent%temperature)&
         * rism3d_closure_partialMolarVolume_IJK(this,cuv,ijk)

    ! b0 + b1*T
    excessChemicalPotential = excessChemicalPotential &
         + (coeff(2) + coeff(4)*this%solvent%temperature)&
         / this%grid%voxelVolume / this%grid%totalGlobalPointsR
  end function rism3d_closure_excessChemicalPotentialUC_ijk


  !! Calculates the excess chemical potential in kT for each site using
  !! the temperature dependent Universal Correction with long
  !! range corrections.
  !!
  !! Uses the total molecular density of the solvent.  For pure water,
  !! this gives the original Palmer correction.  There is no definition
  !! for mixed solvents but this seems reasonable.
  !!
  !! David S Palmer, Andrey I Frolov, Ekaterina L Ratkova and Maxim V Fedorov
  !! J. Phys.: Condens. Matter 22 (2010) 492101
  !! doi:10.1088/0953-8984/22/49/492101
  !!
  !! Temperature dependence from Johnson et al.
  !!
  !! The correction has the form
  !!
  !! \Delta G^UC = \Delta G + (a0+a1*T)*V + b0+b1*T
  !!
  !! To recover the original correction, let a1=b1=0.
  !!
  !!IN:
  !!   this : the closure object
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!   coeff: coefficients for correction: (a0,b0,a1,b1). Extra
  !!          coefficients are ignored.
  !!OUT:
  !!    The total excess chemical potential.  It is not site decomposible.
  function rism3d_closure_aexcessChemicalPotentialUC(this, huv, cuv, coeff) result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:),coeff(:)
    _REAL_ :: excessChemicalPotential

    excessChemicalPotential = sum(rism3d_closure_aexcessChemicalPotential(this,huv,cuv)) &
         + (coeff(1) + coeff(3)*this%solvent%temperature)&
         *rism3d_closure_partialMolarVolume(this,cuv)

    if (this%grid%mpirank==0)&
         excessChemicalPotential = excessChemicalPotential + (coeff(2) + coeff(4)*this%solvent%temperature)
  end function rism3d_closure_aexcessChemicalPotentialUC


  !!Calculates the solvation energy in kT using the temperature
  !!dependent Universal Correction expression (UCT) without long range
  !!corrections. See rism3d_closure_asolvationEnergyUC() for more
  !!details.
  !!
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!   coeff: coefficients for correction: (a0,b0,a1,b1). Extra
  !!          coefficients are ignored.
  !!OUT:
  !!    The total solvation energy.  It is not site decomposible.
  function rism3d_closure_solvationEnergyUC(this, huv_dT, cuv_dT, huv, cuv, coeff) result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:),coeff(:)
    _REAL_ :: solvationEnergy
    integer :: ix, iy,iz
    if (this%solvent%xikt_dT == huge(1d0)) then
       solvationEnergy = huge(1d0)
       return
    end if
    solvationEnergy = 0d0
    ! \Delta E + a*V - (a0-a1*T)*V_dT
    solvationEnergy = sum(rism3d_closure_solvationEnergy(this,huv_dT,cuv_dT,huv,cuv)) &
         + coeff(1)&
         * rism3d_closure_partialMolarVolume(this,cuv)&
         - (coeff(1) + coeff(3)*this%solvent%temperature)&
         * rism3d_closure_partialMolarVolume_dT(this,cuv_dT,cuv)
    ! b 
      if (this%grid%mpirank==0)&
           solvationEnergy = solvationEnergy + coeff(2)
  end function rism3d_closure_solvationEnergyUC


  !! Calculates the solvation energy in kT using the temperature
  !! dependent Universal Correction expression (UCT) at the requested
  !! grid point. See rism3d_closure_asolvationEnergyUC() for more
  !! details.
  !!
  !! Uses the total molecular density of the solvent.  For pure water,
  !! this gives the original Palmer correction.  There is no definition
  !! for mixed solvents but this seems reasonable.
  !!
  !! Original correction:
  !! David S Palmer, Andrey I Frolov, Ekaterina L Ratkova and Maxim V Fedorov
  !! J. Phys.: Condens. Matter 22 (2010) 492101
  !! doi:10.1088/0953-8984/22/49/492101
  !!
  !! Temperature dependent version Johnson et al.
  !!
  !! The correction has the form
  !!
  !! \Delta E^UCT = \Delta E + (a0)*V - (a0-a1*T)*V_dT + b0
  !!
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!   coeff: coefficients for correction: (a0,b0,a1,b1). Extra
  !!          coefficients are ignored.
  !!   ijk     : 3d-grid index
  !!OUT:
  !!    The total solvation energy.  It is not site decomposible.
  function rism3d_closure_solvationEnergyUC_IJK(this, huv_dT, cuv_dT, huv, cuv, coeff,&
       ijk) result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:),coeff(:)
    integer, intent(in) :: ijk(3)
    _REAL_ :: solvationEnergy
    integer :: iv

    ! uncorrected \Delta E
    solvationEnergy = sum(rism3d_closure_solvationEnergy_IJK(&
         this,huv_dT,cuv_dT,huv,cuv,ijk))
    
    ! a0*V - (a0+a1*T)*V_dT

    solvationEnergy = solvationEnergy &
         + coeff(1) &
         * rism3d_closure_partialMolarVolume_IJK(this,cuv,ijk)&
         - (coeff(1) + coeff(3)*this%solvent%temperature)&
         * rism3d_closure_partialMolarVolume_dT_IJK(this,cuv_dT,cuv,ijk)
         
    ! b0
    solvationEnergy = solvationEnergy &
         + coeff(2)&
         /this%grid%totalGlobalPointsR/this%grid%voxelVolume
  end function rism3d_closure_solvationEnergyUC_IJK


  !! Calculates the solvation energy in kT using the temperature dependent Universal
  !! Correction expression (UCT) with long range asymptotic correction. 
  !!
  !! Uses the total molecular density of the solvent.  For pure water,
  !! this gives the original Palmer correction.  There is no definition
  !! for mixed solvents but this seems reasonable.
  !!
  !! Original correction:
  !! David S Palmer, Andrey I Frolov, Ekaterina L Ratkova and Maxim V Fedorov
  !! J. Phys.: Condens. Matter 22 (2010) 492101
  !! doi:10.1088/0953-8984/22/49/492101
  !!
  !! Temperature dependent version Johnson et al.
  !!
  !! The correction has the form
  !!
  !! \Delta E^UCT = \Delta E + (a0)*V - (a0-a1*T)*V_dT + b0
  !!
  !!IN:
  !!   this : the closure object
  !!   huv_dT  : (T*d/dT) site-site total correlation function
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   huv  : site-site total correlation function
  !!   cuv  : site-site direct correlation function
  !!   coeff: coefficients for correction: (a0,b0,a1,b1). Extra
  !!          coefficients are ignored.
  !!   ijk     : 3d-grid index
  !!OUT:
  !!    The total solvation energy.  It is not site decomposible.
  function rism3d_closure_asolvationEnergyUC(this, huv_dT, cuv_dT, huv, cuv, coeff) result(solvationEnergy)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:),coeff(:)
    _REAL_ :: solvationEnergy

    if (this%solvent%xikt_dT == huge(1d0)) then
       solvationEnergy = huge(1d0)
       return
    end if

    ! \Delta E + a*V - (a0-a1*T)*V_dT
    solvationEnergy = sum(rism3d_closure_asolvationEnergy(this,huv_dT,cuv_dT,huv,cuv)) &
         + coeff(1)&
         * rism3d_closure_partialMolarVolume(this,cuv)&
         - (coeff(1) + coeff(3)*this%solvent%temperature)&
         * rism3d_closure_partialMolarVolume_dT(this,cuv_dT,cuv)

    ! b 
    if (this%grid%mpirank==0)&
         solvationEnergy = solvationEnergy + coeff(2)
  end function rism3d_closure_asolvationEnergyUC


  !!Calculate the total solvation interaction energy: de = density sum g*u for
  !!each solvent site.  I.e., the direct intection potential energy of
  !!solute and solvent and not the total solvation energy (see solvationEnergy).
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    guv  :: site-site pair distribution function
  !!OUT:
  !!    the contribution of each solvent site to the total
  !!    solute-solvent potential energy [kT]
  function rism3d_closure_solvPotEne (this,guv) result(ene)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    integer ::  igk, iv, ix,iy,iz
    _REAL_ ::  ene(this%solvent%numAtomTypes)
    ene = 0.d0
#ifdef RISM_DEBUG
    write(6,*)"EXENER"
#endif /*RISM_DEBUG*/

    !!FIX - can we avoid the loops and use BLAS array operations?
    do iv=1,this%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
                ene(iv) = ene(iv) &
                   + rism3d_closure_solvPotEne_ijk (this,guv,(/ix,iy,iz/),iv)
             end do
          end do
       end do
    end do
    !!endfix
    ene = ene*this%grid%voxelVolume
  end function rism3d_closure_solvPotEne


  !!Calculate the total solvation interaction energy: de = density sum g*u
  !!for the request grid point and solvent site.  I.e., the direct
  !!intection potential energy of solute and solvent and not the total
  !!solvation energy (see solvationEnergy).
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    guv  :: site-site pair distribution function
  !!    ijk  : 3d-grid index
  !!OUT:
  !!    the contribution of each solvent site to the total
  !!    solute-solvent potential energy [kT]
  function rism3d_closure_solvPotEne_ijk (this,guv,ijk,iv) result(ene)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    integer, intent(in) :: ijk(3), iv
    integer ::  igk, ix,iy,iz
    _REAL_ ::  ene
    ix=ijk(1)
    iy=ijk(2)
    iz=ijk(3)
#ifdef MPI
    igk = ix + (iy-1)*(this%grid%localDimsR(1)+2) &
             + (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
    igk = ix + (iy-1)*this%grid%localDimsR(1) &
             + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#endif
    ene = guv(igk,iv) * this%potential%uuv(ix,iy,iz,iv) &
                      * this%potential%solvent%density(iv)
  end function rism3d_closure_solvPotEne_ijk


  !!Returns a 3D map of the solvent-solute potential energy in kT.
  !!Integrating this map gives the total solvent-solute potential
  !!energy as returned from rism3d_closure_solvPotEne(,.false.).
  !!Memory is allocated into a pointer and must be freed by the calling
  !!function. For MPI, only the grid points local to this process are
  !!allocated and calculated.
  !!IN:
  !!    this :: rism3d object with computed solution
  !!    guv  :: site-site pair distribution function
  !!OUT:
  !!   a 3D-grid of the solvent-solute energy contributions
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_closure_solvPotEne_tot_map (this,guv) result(solvPotEne)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    integer ::  igk, iv, ix,iy,iz,iatu
    _REAL_,pointer ::  solvPotEne(:,:,:), ulj, uc

    nullify(solvPotEne)
    solvPotEne => safemem_realloc(solvPotEne, &
         this%grid%localDimsR(1), &
         this%grid%localDimsR(2), &
         this%grid%localDimsR(3))
    solvPotEne=0d0
    do iv=1,this%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
                solvPotEne(ix,iy,iz) = solvPotEne(ix,iy,iz) +&
                     rism3d_closure_solvPotEne_ijk(this,guv,(/ix,iy,iz/),iv)
             end do
          end do
       end do
    end do
  end function rism3d_closure_solvPotEne_tot_map


  !!Returns a 3D map of the solvent-solute potential energy in kT.
  !!Integrating this map gives the site solvent-solute potential energy
  !!contribution as returned from rism3d_closure_solvPotEne(,.false.).
  !!Memory is allocated into a pointer and must be freed by the calling
  !!function. For MPI, only the grid points local to this process are
  !!allocated and calculated.
  !!IN:
  !!    this :: rism3d object with computed solution
  !!    guv  :: site-site pair distribution function
  !!OUT:
  !!   a 3D-grid of the solvent-solute energy contributions (nx,ny,nz,nsite)
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_closure_solvPotEne_site_map (this,guv) result(solvPotEne)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    integer ::  igk, iv, ix,iy,iz,iatu
    _REAL_,pointer ::  solvPotEne(:,:,:,:), ulj, uc

    nullify(solvPotEne)
    solvPotEne=>safemem_realloc(solvPotEne, &
         this%grid%localDimsR(1), &
         this%grid%localDimsR(2), &
         this%grid%localDimsR(3), &
         this%solvent%numAtomTypes)
    solvPotEne=0d0
    do iv=1,this%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
                solvPotEne(ix,iy,iz,iv) = solvPotEne(ix,iy,iz,iv) +&
                     rism3d_closure_solvPotEne_ijk(this,guv,(/ix,iy,iz/),iv)
             end do
          end do
       end do
    end do
  end function rism3d_closure_solvPotEne_site_map


  !> Calculates the partial molar volume.
  !! @param[in] this rism3d object with calculated solution.
  !! @param[in] cuv Site-site direct correlation function.
  !! @return The calculated partial molar volume. [A^3]
  function rism3d_closure_partialMolarVolume(this, cuv) result(partialMolarVolume)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: cuv(:,:,:,:)
    _REAL_ :: partialMolarVolume
    _REAL_ :: rCuv(this%solvent%numAtomTypes),rCuv_dT(this%solvent%numAtomTypes)
    integer ::  ix, iy, iz, iv, ierr
    partialMolarVolume = 0d0
    do iv = 1, this%solvent%numAtomTypes
       rcuv(iv) = sum(cuv(:,:,:,iv))
    end do
    partialMolarVolume = -this%solvent%xikt &
         * sum(rcuv* this%solvent%density)&
         *this%grid%voxelVolume
    if (this%grid%mpirank == 0) then
       partialMolarVolume = partialMolarVolume +this%solvent%xikt
    end if
  end function rism3d_closure_partialMolarVolume


  !> Returns a 3D map of the partial molar volume integrand.
  !! Integrating this map gives the partial molar volume as returned
  !! from rism3d_closure_partialMolarVolume(). Memory is allocated
  !! into a pointer and must be freed by the calling function. For
  !! MPI, only the grid points local to this process are allocated and
  !! calculated.
  !! @param[in] this The closure object.
  !! @param[in] cuv Solute-solvent direct correlation function.
  !! @param[out] partialMolarVolume
  !!  A 3D-grid of the PARTIALMOLARVOLUME integrand.
  !! @sideeffect
  !!  Memory is allocated for the grid.
  function rism3d_closure_partialMolarVolume_tot_map(this, cuv) result(partialMolarVolume)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: cuv(:,:,:,:)
    _REAL_,pointer :: partialMolarVolume(:,:,:)
    integer :: ix, iy, iz, iv
    nullify(partialMolarVolume)
    partialMolarVolume => safemem_realloc(partialMolarVolume, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3))
    partialMolarVolume = 0d0
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
             partialMolarVolume(ix,iy,iz) = partialMolarVolume(ix,iy,iz) + &
                  rism3d_closure_partialMolarVolume_ijk &
                  (this, cuv, (/ix,iy,iz/))
          end do
       end do
    end do
  end function rism3d_closure_partialMolarVolume_tot_map


  !> Calculates the partial molar volume at the requested grid point.
  !! @param[in] this rism3d Object with calculated solution.
  !! @param[in] cuv Site-site direct correlation function.
  !! @param[in] ijk 3D-grid index.
  !! @return The calculated partial molar volume. [unitless here]
  function rism3d_closure_partialMolarVolume_IJK(this, cuv, ijk) result(partialMolarVolume)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: cuv(:,:,:,:)
    integer, intent(in) :: ijk(3)
    _REAL_ :: partialMolarVolume
    integer :: iv, ierr
    partialMolarVolume = this%solvent%xikt &
          * (1d0 / this%grid%totalGlobalPointsR / this%grid%voxelVolume &
         - sum(cuv(ijk(1), ijk(2), ijk(3), :) * this%solvent%density))
  end function rism3d_closure_partialMolarVolume_IJK


  !> Calculates the partial molar volume temperature derivative (T*d/dT).
  !! @param[in] this rism3d object with calculated solution.
  !! @param[in] cuv_dT Site-site (T*d/dT) direct correlation function.
  !! @param[in] cuv Site-site direct correlation function.
  !! @return The calculated partial molar volume. [A^3]
  function rism3d_closure_partialMolarVolume_dT(this, cuv_dT, cuv) result(partialMolarVolume_dT)
    use rism_util, only: checksum
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: cuv(:,:,:,:),cuv_dT(:,:,:,:)
    _REAL_ :: partialMolarVolume_dT
    _REAL_ :: rCuv(this%solvent%numAtomTypes),rCuv_dT(this%solvent%numAtomTypes)
    integer :: iv
    if (this%solvent%xikt_dT == huge(1d0)) then
       partialMolarVolume_dT = huge(1d0)
       return
    end if

    do iv = 1, this%solvent%numAtomTypes
       rcuv(iv) = sum(cuv(:,:,:,iv))
       rcuv_dT(iv) = sum(cuv_dT(:,:,:,iv))
    end do
    partialMolarVolume_dT = rism3d_closure_partialMolarVolume(this, cuv) &
         -(this%solvent%xikt_dT &
         * sum(rcuv* this%solvent%density) &
         + this%solvent%xikt&
         *sum(rcuv_dT * this%solvent%density))&
         *this%grid%voxelVolume
    if (this%grid%mpirank == 0) then
       partialMolarVolume_dT = partialMolarVolume_dT +this%solvent%xikt_dT
    end if
  end function rism3d_closure_partialMolarVolume_dT


  !!Returns a 3D map of the partial molar volume temperature derivative
  !!(T*d/dT) integrand.  Integrating this map gives the partial molar
  !!volume temperature derivative as returned from
  !!rism3d_closure_partialMolarVolume_dT().  Memory is allocated into a pointer and
  !!must be freed by the calling function. For MPI, only the grid
  !!points local to this process are allocated and calculated.
  !!IN:
  !!   this : the closure object
  !!   cuv_dT  : (T*d/dT) site-site direct correlation function
  !!   cuv  : site-site direct correlation function
  !!OUT:
  !!   a 3D-grid of the PARTIALMOLARVOLUME integrand
  !!SIDEEFFECTS:
  !!   memory is allocated for the grid
  function rism3d_closure_PARTIALMOLARVOLUME_dT_tot_map(this, cuv_dT, cuv) result(partialMolarVolume_dT)
    use rism_util, only : caseup
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: cuv_dT(:,:,:,:)
    _REAL_, intent(in) :: cuv(:,:,:,:)
    _REAL_,pointer :: partialMolarVolume_dT(:,:,:)
    integer :: ix, iy, iz, iv
    nullify(partialMolarVolume_dT)
    partialMolarVolume_dT => safemem_realloc(partialMolarVolume_dT, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3))
    partialMolarVolume_dT=0d0
    do iz = 1, this%grid%localDimsR(3)
       do iy = 1, this%grid%localDimsR(2)
          do ix = 1, this%grid%localDimsR(1)
             partialMolarVolume_dT(ix,iy,iz) = partialMolarVolume_dT(ix,iy,iz)+&
                  rism3d_closure_partialMolarVolume_dT_ijk&
                  (this, cuv_dT,cuv,(/ix,iy,iz/))
          end do
       end do
    end do
  end function rism3d_closure_partialMolarVolume_dT_tot_map

  
  !> Calculates the partial molar volume temperature derivative
  !! (T*d/dT) at the requested grid point.
  !! @param[in] this rism3d object with calculated solution.
  !! @param[in] cuv_dT Site-site (T*d/dT) direct correlation function.
  !! @param[in] cuv Site-site direct correlation function.
  !! @param[in] ijk 3D-grid index.
  !! @return The calculated partial molar volume. [unitless here]
  function rism3d_closure_partialMolarVolume_dT_IJK(this, cuv_dT, cuv, ijk) result(partialMolarVolume_dT)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: cuv(:,:,:,:), cuv_dT(:,:,:,:)
    integer, intent(in) :: ijk(3)
    _REAL_ :: partialMolarVolume_dT
    integer :: ix, iy, iz, iv, ierr
    ix = ijk(1)
    iy = ijk(2)
    iz = ijk(3)
    partialMolarVolume_dT = rism3d_closure_partialMolarVolume_IJK(this, cuv, ijk) &
         + this%solvent%xikt_dT &
         * (1d0 / this%grid%totalGlobalPointsR / this%grid%voxelVolume &
         - sum(cuv(ix, iy, iz, :) * this%solvent%density)) &
         - sum(cuv_dT(ix, iy, iz, :) * this%solvent%density) * this%solvent%xikt
  end function rism3d_closure_partialMolarVolume_dT_IJK

  
  !> Calculate the excess number of particles about the solute
  !! compared to the bulk solvation.  No attempt is made to account
  !! for excluded volume.  This is also-known-as a molar preferential
  !! interaction parameter.
  !! @param[in] this rism3d object with computed solution.
  !! @param[in] guv site-site pair distribution function.
  !! @return The number of excess particles for each solvent site.
  function rism3d_closure_excessParticles (this, guv) result(num)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    _REAL_ :: num(this%solvent%numAtomTypes)
    num = this%solvent%density * rism3d_closure_kirkwoodBuff(this, guv)
  end function rism3d_closure_excessParticles

  
  !> Calculate the excess number of particles about the solute
  !! compared to the bulk solvation corrected with the long-range
  !! asymptotics.  No attempt is made to account for excluded volume.
  !! This is also-known-as as a molar preferential interaction
  !! parameter.
  !! @param[in] this rism3d object with computed solution.
  !! @param[in,out] guv Site-site pair distribution function.
  !! @return The number of excess particles for each solvent site.
  function rism3d_closure_aexcessParticles (this, guv) result(num)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    _REAL_ :: num(this%solvent%numAtomTypes)
    num = this%solvent%density * rism3d_closure_akirkwoodBuff(this, guv)
  end function rism3d_closure_aexcessParticles


  !!Calculate the temperature derivative (T*d/dT) of the excess number
  !!of particles about the solute compared to the bulk solvation.  No
  !!attempt is made to account for excluded volume.  This is
  !!also-known-as as a molar preferential interaction parameter
  !!temperature derivative.
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    guv_dT  :: (T*d/dT) site-site pair distribution function
  !!OUT:
  !!    The temperature derivative of the number of excess particles
  !!    for each solvent site.
  function rism3d_closure_excessParticles_dT (this, guv_dT) result(num)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv_dT(:,:)
    _REAL_ :: num(this%solvent%numAtomTypes)
    num = this%solvent%density * rism3d_closure_kirkwoodBuff_dT(this, guv_dT)
  end function rism3d_closure_excessParticles_dT


  !> Calculate the temperature derivative (T*d/dT) of the excess number of
  !! particles about the solute compared to the bulk solvation corrected
  !! with the longrange asymptotics.  No attempt is made to account for
  !! excluded volume.  This is also-known-as as a molar preferential
  !! interaction parameter temperature derivative.
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    guv_dT  :: (T*d/dT) site-site pair distribution function
  !!OUT:
  !!    The temperature derivative of the number of excess particles
  !!    for each solvent site.
  function rism3d_closure_aexcessParticles_dT (this, guv_dT) result(num)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv_dT(:,:)
    _REAL_ :: num(this%solvent%numAtomTypes)
    num = this%solvent%density * rism3d_closure_akirkwoodBuff_dT(this, guv_dT)
  end function rism3d_closure_aexcessParticles_dT


  !> Calculate the Kirkwood-Buff integral for the solute. This is the
  !! all space integral of h_{uv}.
  !!
  !! J. G. Kirkwood; F. P. Buff. J. Chem. Phys. 1951, 19, 774-777
  !! @param[in] this rism3d object with computed solution.
  !! @param[in] guv Site-site pair distribution function.
  !! @return Kirkwood-Buff integral for each solvent site.
  function rism3d_closure_kirkwoodBuff (this, guv) result(H)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    _REAL_ :: H(this%solvent%numAtomTypes)
    integer :: iv, ix, iy, iz, igk
    
    H = 0
    do iv = 1, this%solvent%numAtomTypes
       do iz = 1, this%grid%localDimsR(3)
          do iy = 1, this%grid%localDimsR(2)
             do ix = 1, this%grid%localDimsR(1)
#ifdef MPI
                igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                                    * (this%grid%localDimsR(1) + 2)
#else
                igk = ix + (iy - 1) * this%grid%localDimsR(1) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                         * this%grid%localDimsR(1)
#endif
                H(iv) = H(iv) + (guv(igk, iv) - 1d0)
             end do
          end do
       end do
    end do
    H = H * this%grid%voxelVolume
  end function rism3d_closure_kirkwoodBuff


  !> Calculate the Kirkwood-Buff integral for the solute w/ long-range
  !! asymptotic correction. This is the all space integral of h_{uv}.
  !!
  !! J. G. Kirkwood; F. P. Buff. J. Chem. Phys. 1951, 19, 774-777
  !! @param[in] this rism3d object with computed solution.
  !! @param[in] guv Site-site pair distribution function.
  !! @return Kirkwood-Buff integral for each solvent site.
  function rism3d_closure_akirkwoodBuff (this,guv) result(H)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    _REAL_ :: H(this%solvent%numAtomTypes)
    integer :: iv, ix, iy, iz, ig, igk

    if (.not. this%solvent%ionic) then
       H = rism3d_closure_kirkwoodbuff(this, guv)
       return
    end if
    H = 0
    do iv = 1, this%solvent%numAtomTypes
       do iz = 1, this%grid%localDimsR(3)
          do iy = 1, this%grid%localDimsR(2)
             do ix = 1, this%grid%localDimsR(1)
                ig = ix + (iy - 1) * this%grid%localDimsR(1) &
                     + (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#ifdef MPI
                igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                                    * (this%grid%localDimsR(1) + 2)
#else
                igk = ix + (iy - 1) * this%grid%localDimsR(1) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                                    * this%grid%localDimsR(1)
#endif
                H(iv) = H(iv) + (guv(igk, iv) - 1d0) &
                              - this%potential%solvent%charge_sp(iv) &
                              * this%potential%tcfLongRangeAsympR(ig)
             end do
          end do
       end do
    end do
    H = H * this%grid%voxelVolume
    ! Use a smear of 0 to ensure electroneutrality.
    if (sum(this%grid%offsetR) == 0) &
         H = H + rism3d_potential_int_h(this%potential)
  end function rism3d_closure_akirkwoodBuff


  !!Calculates the Kirkwood-Buff integral temperature derivative for the
  !!solute. This is the all space integral of huv_dT.
  !!
  !!J. G. Kirkwood; F. P. Buff. J. Chem. Phys. 1951, 19, 774-777
  !! IN:
  !!    this    :: rism3d object with computed solution
  !!    guv_dT  :: (T*d/dT) site-site pair distribution function
  !!OUT:
  !!    Kirkwood-Buff integeral temperature derivative for each solvent site
  function rism3d_closure_kirkwoodBuff_dT (this,guv_dT) result(G_dT)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv_dT(:,:)
    _REAL_ ::  G_dT(this%solvent%numAtomTypes)
    integer :: iv, ix, iy, iz, igk
    G_dT=0
    do iv=1,this%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
#ifdef MPI
                igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                                    * (this%grid%localDimsR(1) + 2)
#else
                igk = ix + (iy - 1) * this%grid%localDimsR(1) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                                    * this%grid%localDimsR(1)
#endif
                G_dT(iv) = G_dT(iv) + (guv_dT(igk,iv))
             end do
          end do
       end do
    end do
    G_dT = G_dT*this%grid%voxelVolume
  end function rism3d_closure_kirkwoodBuff_dT


  !!Calculate the Kirkwood-Buff integral temperature derivative for the
  !!solute w/ long-range asymptotic correction. This is the all space
  !!integral of huv_dT.
  !!
  !!J. G. Kirkwood; F. P. Buff. J. Chem. Phys. 1951, 19, 774-777
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    guv_dT  :: (T*d/dT) site-site pair distribution function
  !!OUT:
  !!    Kirkwood-Buff integral temperature derivative for each solvent site.
  function rism3d_closure_akirkwoodBuff_dT (this,guv_dT) result(G_dT)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv_dT(:,:)
    _REAL_ ::  G_dT(this%solvent%numAtomTypes)
    integer :: iv, ix, iy, iz, ig, igk

    if (.not.this%solvent%ionic .or. .not.this%solute%charged) then
       G_dT=rism3d_closure_kirkwoodbuff_dT(this,guv_dT)
       return
    end if
    G_dT=0
    do iv=1,this%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
                ig = ix + (iy-1)*this%grid%localDimsR(1) + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#ifdef MPI
                igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                                    * (this%grid%localDimsR(1) + 2)
#else
                igk = ix + (iy - 1) * this%grid%localDimsR(1) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                                    * this%grid%localDimsR(1)
#endif
                G_dT(iv) = G_dT(iv) + (guv_dT(igk,iv)) &
                     - this%potential%solvent%charge_sp(iv) &
                       * this%potential%tcfLongRangeAsympR(ig)
             end do
          end do
       end do
    end do
    G_dT = G_dT*this%grid%voxelVolume
    if (sum(this%grid%offsetR) == 0)&
         G_dT=G_dT+rism3d_potential_int_h(this%potential)
  end function rism3d_closure_akirkwoodBuff_dT


  !!Calculates the direct correlation function integral for the solute. This is the
  !!all space integral of cuv.
  !!
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    cuv  :: site-site pair direct correlation function
  !!OUT:
  !!    DCF integeral for each solvent site
  function rism3d_closure_DCFintegral (this,cuv) result(C)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: cuv(:,:,:,:)
    _REAL_ ::  C(this%solvent%numAtomTypes)
    integer :: iv, ix, iy, iz, igk
    do iv=1,this%solvent%numAtomTypes
       C(iv) = sum(cuv(:,:,:,iv))
    end do
    C = C*this%grid%voxelVolume
  end function rism3d_closure_DCFintegral


  !!Calculates the direct correlation function integral for the solute
  !!w/ long-range asymptotic correction. This is the all space integral
  !!of cuv.
  !!
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    cuv  :: site-site pair direct correlation function
  !!OUT:
  !!    DCF integeral for each solvent site

!!$    function rism3d_closure_aDCFintegral (this,cuv) result(C)
!!$      implicit none
!!$      type(rism3d_closure) :: this
!!$      _REAL_, intent(in) :: cuv(:,:)
!!$      _REAL_ ::  C(this%solvent%numAtomTypes)
!!$      integer :: iv, ix, iy, iz, ig, igk
!!$      !can we use asymphk(0) here?
!!$      do iv=1,this%solvent%numAtomTypes
!!$         C(iv) = sum(cuv(:,:,:,iv)) -sum(this%potential%dcfLongRangeAsympR)*this%solvent%charge(iv)
!!$      end do
!!$      C = C*this%grid%voxelVolume
!!$      if (sum(this%grid%offsetR) == 0)&
!!$           C=C+rism3d_potential_int_c(this%potential)
!!$    end function rism3d_closure_aDCFintegral


  !!Calculates the direct correlation function temperature derivative
  !!integral for the solute. This is the all space integral of cuv_dT.
  !!
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    cuv_dT :: (T*d/dT) site-site pair direct correlation function
  !!OUT:
  !!    DCF temperature derivative integeral for each solvent site

  function rism3d_closure_DCFintegral_dT (this,cuv_dT) result(C_dT)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: cuv_dT(:,:,:,:)
    _REAL_ ::  C_dT(this%solvent%numAtomTypes)
    integer :: iv, ix, iy, iz, igk
    do iv=1,this%solvent%numAtomTypes
       C_dT(iv) = sum(cuv_dT(:,:,:,iv))
    end do
    C_dT = C_dT*this%grid%voxelVolume
  end function rism3d_closure_DCFintegral_dT


  !> Calculates the forces on the solute contributed by the solvent
  !! according to 3D-RISM. In fact, this subroutine calls the
  !! appropriate subroutines to calculate this.
  subroutine rism3d_closure_force(this, ff, guv, periodicPotential)
    implicit none
    type(rism3d_closure):: this !< Closure object with computed solution.
    _REAL_, intent(out) :: ff(3,this%solute%numAtoms) !< 3D-RISM forces [kT/A].
    _REAL_, intent(in) :: guv(:,:) !< Site-site pair distribution function.
    character(len=*), intent(in) :: periodicPotential !< Label for
    ! periodic potential, else empty string.
    integer :: atom
    integer :: atomRange
#ifdef RISM_DEBUG
    write(6,*) "RISM_FF"
    call flush(6)
#endif /*RISM_DEBUG*/
    ff = 0
    if (periodicPotential == '') then
       if (this%potential%cutoff2 > sum(this%grid%boxLength**2)) then
          call force_brute(this, ff, guv)
       else
          call lennardJonesForce(this, ff, guv, .false.)
          call coulombicForce(this, ff, guv)
       endif
    else
       call lennardJonesForcePeriodic(this%potential, ff, guv)
       if (periodicPotential == 'pme') then
          call particleMeshEwaldForce(this%potential, ff, guv)
       else
          call ewaldSumForce(this%potential, ff, guv)
       end if
    end if

  end subroutine rism3d_closure_force


  !!Frees memory and resets object state
  !!IN:
  !!   this : the closure object

  subroutine rism3d_closure_destroy(this)
    use safemem
    implicit none
    type(rism3d_closure), intent(inout) :: this
    if (associated(this%kh)) then
       call rism3d_kh_destroy(this%kh)
       deallocate(this%kh)
    end if
    if (associated(this%psen)) then
       call rism3d_psen_destroy(this%psen)
       deallocate(this%psen)
    end if
    if (associated(this%hnc)) then
       call rism3d_hnc_destroy(this%hnc)
       deallocate(this%hnc)
    end if
    nullify(this%potential)
    nullify(this%grid)
    nullify(this%solute)
    nullify(this%solvent)
  end subroutine rism3d_closure_destroy


  !!                         PRIVATE

  !> Lennard-Jones contributions to the total mean solvation forces.
  !! Forces remain distributed over the nodes for MPI calculations.
  subroutine lennardJonesForce(this, ff, guv, periodic)
    use rism3d_potential_c, only : minimumImage
    implicit none
    type(rism3d_closure), intent(in):: this !< Closure object with computed solution.
    _REAL_, intent(inout) :: ff(3, this%solute%numAtoms) !< Force array [kT/A].
    _REAL_, intent(in) :: guv(:, :) !< Site-site pair correlation function.
    !TODO: Maybe move this to closure object?
    logical, intent(in) :: periodic
    
    integer :: igx, igy, igz, igA(3), ig, ig2, igk, iu, iv, id, idim1, idim2
    _REAL_ :: rx, ry, rz, dX(3), dz2, dyz2, rs2, rs6i, r5, r4, r3, r2, ra
    _REAL_ :: offset
    _REAL_, parameter :: rcor = 0.002d0
    _REAL_, parameter :: rcor2 = rcor**2

    _REAL_ :: dUlj_dr,dU_dr(3)
    ! Linear spacing of the grid.
    _REAL_ :: ff_temp(3, this%solute%numAtoms)
    ! Number of gridpoints in each direction to use closest point to
    ! solute atom.
    integer :: grdpnts(3), cp(3), first(3), last(3), iu2, id2, ierr

#ifdef RISM_DEBUG
    write(6,*) "FORCE_LJ"; call flush(6)
#endif /*RISM_DEBUG*/
    
    ff = 0d0

    offset = this%grid%spacing(3) * this%grid%offsetR(3)
    
    ! Calculate the number of grid point necessary to cover this
    ! cutoff range. In case of a large cutoff, we need to protect
    ! against integer overflow.
    grdpnts = nint(min(dble(huge(1) - 1), &
         sqrt(this%potential%cutoff2) / this%grid%spacing)) + 1
    do iu = 1, this%solute%numAtoms
       do id = 1, 3
          cp(id) = nint(this%solute%position(id, iu) / this%grid%spacing(id))
          first(id) = max(0, cp(id) - grdpnts(id))
          ! Note: we have to protect against cp(id) + grdpnts(id) overflowing.
          last(id) = min(this%grid%globalDimsR(id) - 1, &
               cp(id) + min(huge(1) - cp(id), grdpnts(id)))
       end do
       cp(3) = cp(3) - this%grid%offsetR(3)
       first(3) = max(0, first(3) - this%grid%offsetR(3))
       last(3) = min(this%grid%localDimsR(3) - 1, last(3) - this%grid%offsetR(3))
       do igz = first(3), last(3)
          rz = igz * this%grid%spacing(3) + offset
          dX(3) = this%solute%position(3, iu) - rz
          dz2 = dX(3) * dX(3)
          do igy = first(2), last(2)
             ry = igy * this%grid%spacing(2)
             dX(2) = this%solute%position(2, iu) - ry
             dyz2 = dX(2) * dX(2) + dz2
             do igx = first(1), last(1)
                rx = igx * this%grid%spacing(1)
#ifdef MPI
                ig = 1 + igx + igy * (this%grid%globalDimsR(1) + 2) &
                       + igz * this%grid%globalDimsR(2) &
                             * (this%grid%globalDimsR(1) + 2)
#else
                !FIXME: Force functions seem to use an odd mix of
                ! local and global grid dimensions.
                ig = 1 + igx + igy * this%grid%globalDimsR(1) &
                       + igz * this%grid%globalDimsR(2) &
                             * this%grid%globalDimsR(1)
#endif

                dX(1) = this%solute%position(1,iu) - rx
                r2 = dX(1) * dX(1) + dyz2
                if (r2 < rcor2) r2 = rcor2
                if (r2 < this%potential%cutoff2) then
                   do iv = 1, this%solvent%numAtomTypes
                      rs2 = r2 / this%potential%ljSigmaUV(iu, iv)**2
                      rs6i = 1d0 / rs2**3
                      dUlj_dr = 12.d0 * this%potential%ljEpsilonuv(iu, iv) &
                           * rs6i * (rs6i - 1.d0) &
                           / r2 * this%solvent%density(iv) * guv(ig, iv)
                      do id = 1, 3
                         ff(id, iu) = ff(id, iu) + dUlj_dr * dX(id)
                      end do
                   end do
                end if
             end do
          end do
       end do
    end do
    ff = ff * this%grid%voxelVolume
  end subroutine lennardJonesForce


  !> Tabulate the solute-solvent 12-6 Lennard-Jones force in the
  !! box subject to the minimum image convention.
  subroutine lennardJonesForcePeriodic(this, ff, guv)
    use constants, only : PI
    use rism_util, only : checksum
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_potential), intent(in) :: this !> potential object.
    _REAL_, intent(inout) :: ff(3, this%solute%numAtoms) !> Force on each atom.
    _REAL_, intent(in) :: guv(:, :) !< Site-site pair correlation function.
    
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
    ! Minimum distance to prevent divide by zero.
    _REAL_, parameter :: minDistance = 0.002d0
    _REAL_, parameter :: minDistanceSquared = minDistance**2
    _REAL_ :: solutePosition(3)
    ! Lennard-Jones force between two particles.
    _REAL_ :: dUlj_dr
    
    ! Offset used for dividing the box along the z-axis
    ! for MPI.
    !TODO:JAJ Need to adjust for MPI z-slab decomposition.
    offset = this%grid%spacing(3) * this%grid%offsetR(3)

    ff = 0d0

    do iu = 1, this%solute%numAtoms
       ! Calculate the solute-solvent Lennard-Jones potential at each
       ! grid point in the box.
       do igz = 1, this%grid%localDimsR(3)
          rz(:) = (igz - 1) * this%grid%voxelVectorsR(3,:)
          rz(:) = rz(:) + offset
          do igy = 1, this%grid%localDimsR(2)
             ry(:) = (igy - 1) * this%grid%voxelVectorsR(2,:)
             do igx = 1, this%grid%localDimsR(1)
                rx(:) = (igx - 1) * this%grid%voxelVectorsR(1,:)
#ifdef MPI
                ig = 1 + igx + igy * (this%grid%globalDimsR(1) + 2) &
                       + igz * this%grid%globalDimsR(2) &
                             * (this%grid%globalDimsR(1) + 2)
#else
                ig = 1 + (igx - 1) + (igy - 1) * this%grid%globalDimsR(1) &
                       + (igz - 1) * this%grid%globalDimsR(2) &
                                   * this%grid%globalDimsR(1)
#endif

                gridPoint(:) = rx(:) + ry(:) + rz(:)
                
                ! Vector from solute atom to grid point, in real space.
                solutePosition(:) = gridPoint(:) - this%solute%position(:, iu)
                
                ! Adjust distance to be to the closest solute atom
                ! image, which may be in an adjacent cell.
                solutePosition = minimumImage(this, solutePosition)
                
                soluteDistanceSquared = dot_product(solutePosition, solutePosition)
                if (soluteDistanceSquared < minDistanceSquared) &
                     soluteDistanceSquared = minDistanceSquared
                if (soluteDistanceSquared < this%cutoff2) then
                   dUlj_dr = 0
                   do iv = 1, this%solvent%numAtomTypes
                      ljBaseTerm = soluteDistanceSquared / this%ljSigmaUV(iu, iv)**2
                      ljBaseTerm = 1d0 / ljBaseTerm**3
                      dUlj_dr = dUlj_dr + this%ljEpsilonUV(iu, iv) &
                           * ljBaseTerm * (ljBaseTerm - 1.d0) &
                            * this%solvent%density(iv) * guv(ig, iv)
                   end do
                   ff(:, iu) = ff(:, iu) &
                        - dUlj_dr * solutePosition(:) / soluteDistanceSquared
                end if
             end do
          end do
       end do
    end do
    ff = ff * 12d0 * this%grid%voxelVolume
  end subroutine lennardJonesForcePeriodic

  
  !> Electrostatic contributions to the total mean solvation forces
  !! Forces remain distributed over the nodes for MPI calculations.
  subroutine coulombicForce(this, ff, guv)
    implicit none
    type(rism3d_closure), intent(in) :: this !< Closure object with
                                             !! computed solution.
    _REAL_, intent(inout) :: ff(3, this%solute%numAtoms) !< Force array [kT/A].
    _REAL_, intent(in) :: guv(:,:) !< Site-site pair correlation function.
    
    _REAL_ :: frc(3)
    integer ::  igx, igy, igz, igA(3), ig, igk, iu, iv, id, idim1, idim2
    _REAL_ ::  rx, ry, rz, dX(3), dz2, dyz2, r2, ra
    _REAL_, parameter :: minimumDistance = 0.002d0
    _REAL_, parameter :: minimumDistance2 = minimumDistance**2

    _REAL_ ::  cff
    _REAL_ ::  dUc_dr, dU_dr(3)
    _REAL_ :: qutemp(this%solute%numAtoms), qvtemp(this%solvent%numAtomTypes)
    !linear spacing of the grid
    _REAL_ :: offset, ff_temp(3, this%solute%numAtoms)
    !smstart :: initial grid point where the sparse and fine grids have the same value
    integer :: smstart
    integer :: ierr, irank
    !number of gridpoints in each direction to use
    !closest point to solute atom
    integer :: grdpnts(3), cp(3), first(3), last(3), iu2, id2

    smstart = mod(this%grid%offsetR(3), 2)
    qutemp(1:this%solute%numAtoms) = this%solute%charge(1:this%solute%numAtoms)
    qvtemp(1:this%solvent%numAtomTypes) = &
        this%solvent%charge(1:this%solvent%numAtomTypes)
    offset = this%grid%spacing(3) * this%grid%offsetR(3)

    ! Calculate the number of grid point necessary to cover this
    ! cutoff range. In case of a large cutoff, we need to protect
    ! against integer overflow.
    ! Make it an odd number.
    grdpnts = nint(min(dble(huge(1) - 1), &
         sqrt(this%potential%cutoff2) / this%grid%spacing)) + 1
    do iu = 1, this%solute%numAtoms
       frc = 0
       do id = 1, 3
          cp(id) = nint(this%solute%position(id, iu) / this%grid%spacing(id))
          first(id) = max(0, cp(id) - grdpnts(id))
          ! Start on an even number.
          first(id) = first(id) + mod(first(id), 2)
          ! End on an odd number.
          last(id) = min(this%grid%globalDimsR(id) - 1, &
               first(id) + 2 * grdpnts(id) - 1)
       end do
       cp(3) = cp(3) - this%grid%offsetR(3)
       first(3) = max(0, first(3) - this%grid%offsetR(3))
       last(3) = min(this%grid%localDimsR(3)-1, last(3)-this%grid%offsetR(3))
       do igz = first(3), last(3)
          rz = igz * this%grid%spacing(3) + offset
          dX(3) =  this%solute%position(3, iu) - rz
          dz2 = dX(3) * dX(3)
          do igy = first(2), last(2)
             ry = igy * this%grid%spacing(2)
             dX(2) =  this%solute%position(2, iu) - ry
             dyz2 = dX(2) * dX(2) + dz2
             do igx = first(1), last(1)

#ifdef MPI
                ig = 1 + igx + igy * (this%grid%globalDimsR(1) + 2) &
                       + igz * this%grid%globalDimsR(2) &
                             * (this%grid%globalDimsR(1) + 2)
#else
                ig = 1 + igx + igy * this%grid%globalDimsR(1) &
                       + igz * this%grid%globalDimsR(2) &
                             * this%grid%globalDimsR(1)
#endif
                rx = igx * this%grid%spacing(1)

                dX(1) =  this%solute%position(1, iu) - rx
                r2 = dX(1) * dX(1) + dyz2
                if (r2 < minimumDistance2) r2 = minimumDistance2
                ra = sqrt(r2)
                do iv = 1, this%solvent%numAtomTypes
                   dUc_dr  = qutemp(iu) * qvtemp(iv) / r2 / ra &
                        * this%solvent%density(iv) * guv(ig, iv)
                   do id = 1, 3
                      frc(id) = frc(id) + dUc_dr * dX(id)
                   end do
                end do
             end do
          end do
       end do
       
       ff(:,iu) = ff(:,iu) + frc * this%grid%voxelVolume
       frc = 0
       
       ! Long range coarse grid.
       do igz = smstart, this%grid%localDimsR(3) - 1, 2
          rz = igz * this%grid%spacing(3) + offset
          dX(3) =  this%solute%position(3, iu) - rz
          dz2 = dX(3) * dX(3)
          do igy = 0, this%grid%localDimsR(2) - 1, 2
             ry = igy * this%grid%spacing(2)
             dX(2) =  this%solute%position(2, iu) - ry
             dyz2 = dX(2) * dX(2) + dz2
             do igx = 0, this%grid%localDimsR(1) - 1, 2
                if ((igx >= first(1) .and. igx <= last(1)) .and. &
                     (igy >= first(2) .and. igy <= last(2)) .and. &
                     (igz >= first(3) .and. igz <= last(3))) cycle
#ifdef MPI
                ig = 1 + igx + igy * (this%grid%localDimsR(1) + 2) &
                       + igz * this%grid%localDimsR(2) &
                             * (this%grid%localDimsR(1) + 2)
#else
                ig = 1 + igx + igy * this%grid%localDimsR(1) &
                       + igz * this%grid%localDimsR(2) &
                             * this%grid%localDimsR(1)
#endif
                rx = igx*this%grid%spacing(1)

                dX(1) =  this%solute%position(1, iu) - rx
                r2 = dX(1) * dX(1) + dyz2
                ra = sqrt(r2)
                do iv = 1, this%solvent%numAtomTypes
                   dUc_dr  = -qutemp(iu) * qvtemp(iv) / r2 / ra &
                        * this%solvent%density(iv) * guv(ig, iv)
                   do id = 1, 3
                      frc(id) = &
                           frc(id) - dUc_dr * dX(id)
                   end do
                end do
             end do
          end do
       end do
       ff(:, iu) = ff(:, iu) + frc * this%grid%voxelVolume * 8d0
    end do
  end subroutine coulombicForce

  
  !> Short-range portion of the Ewald sum electric force exerted on
  !! solute by the solvent.
  !! Note: never called in current code
  !! Reference: 
  !! Deserno M and Holm C. How To Mesh Up Ewald Sums I. Eq. 14.
  !! J. Chem. Phys., Vol. 109, No. 18, 8 November 1998.
  subroutine ewaldSumShortRangeForce(this, ff, guv)
    use constants, only: pi
    use rism3d_opendx, only : rism3d_opendx_write
    implicit none
    type(rism3d_potential), intent(inout) :: this !< potential object.
    _REAL_, intent(inout) :: ff(3, this%solute%numAtoms) !< Force array [kT/A].
    _REAL_, intent(in) :: guv(:,:) !< Site-site pair correlation function.

    ! Minimum distance to prevent division by zero.
    _REAL_, parameter :: minDistance = 0.002
    _REAL_, parameter :: minDistance2 = minDistance**2
    
    integer :: igx, igy, igz
    integer :: ig
    integer :: iu, iv
    integer :: id

    _REAL_ :: rx(3), ry(3), rz(3)
    _REAL_ :: gridPoint(3), solutePosition(3)
    _REAL_ :: soluteDistance, soluteDistance2
    _REAL_ :: dUc_dr, frc(3)

    ! For MPI spatial decomposition.
    ! Distance in z-axis due to spatial decomp
    _REAL_ :: offset

    _REAL_ :: smear2

    _REAL_, external :: erfc

    smear2 = this%solvent%smear**2
    
    ! Calculate the number of grid point necessary to cover this
    ! cutoff range.  In case of a large cutoff, we need to protect
    ! against integer overflow.
    offset = this%grid%spacing(3) * (this%grid%offsetR(3))
    
    ! Calculate short-range term of Ewald sum and combine with
    ! previously calculated long-range term.
    do iu = 1, this%solute%numAtoms
       frc = 0
       do igz = 0, this%grid%globalDimsR(3) - 1
          rz = igz * this%grid%voxelVectorsR(3, :)
          do igy = 0, this%grid%globalDimsR(2) - 1
             ry = igy * this%grid%voxelVectorsR(2, :)
             do igx = 0, this%grid%globalDimsR(1) - 1
                ig = 1 + igx + igy * this%grid%globalDimsR(1) + &
                     igz * this%grid%globalDimsR(2) * this%grid%globalDimsR(1)
                rx = igx * this%grid%voxelVectorsR(1, :)

                gridPoint = rx + ry + rz
                
                ! Solute atom position relative to a grid point.
                solutePosition = gridPoint - this%solute%position(:, iu)

                ! Apply minimum image convention.
                solutePosition = minimumImage(this, solutePosition)

                ! Distance from solute atom to grid point.
                soluteDistance2 = dot_product(solutePosition, solutePosition)

                ! Short-range term of Ewald sum.
                if (soluteDistance2 < this%cutoff2) then
                   if (soluteDistance2 < minDistance2) &
                        soluteDistance2 = minDistance2
                   soluteDistance = sqrt(soluteDistance2)
                   
                   dUc_dr = 0
                   do iv = 1, this%solvent%numAtomTypes
                      dUc_dr = dUc_dr + this%solvent%density(iv) &
                           * guv(ig, iv) * this%solvent%charge(iv)
                   end do
                   dUc_dr = dUc_dr * &
                        (2 / (this%solvent%smear * sqrt(pi)) &
                        * exp(-soluteDistance2 / smear2) &
                        + erfc(soluteDistance / this%solvent%smear) &
                        / soluteDistance)
                   
                   frc = frc + dUc_dr * soluteDistance / soluteDistance2
                end if
             end do
          end do
       end do
       ff(:, iu) = ff(:, iu) + frc * this%grid%voxelVolume * this%solute%charge(iu)
    end do
  end subroutine ewaldSumShortRangeForce


  !> 
  subroutine ewaldSumForce(this, ff, guv)
    use constants, only : PI, FOURPI
#ifdef __PGI
    use ieee_arithmetic, only : ieee_is_nan
#endif
    implicit none
    type(rism3d_potential), intent(inout) :: this !< potential object.
    _REAL_, intent(inout) :: ff(3, this%solute%numAtoms) !< Force array [kT/A].
    _REAL_, intent(in) :: guv(:,:) !< Site-site pair correlation function.

    integer :: ig, igx, igy, igz, igk
    integer :: iu, iv

    _REAL_ :: rx(3), ry(3), rz(3)
    _REAL_ :: gridPoint(3), solutePosition(3)
    _REAL_ :: soluteDistance, soluteDistance2
    _REAL_ :: dUc_dk, kforce(3), kforceSum(3), phase
    _REAL_ :: smear2_4
    _REAL_ :: totalSolventCharge

    smear2_4 = this%solvent%smear**2 / 4d0
    
    
    ! if (this%grid%offsetK(3) == 0) then
    !    this%dcfLongRangeAsympK(1) = 0d0
    !    this%dcfLongRangeAsympK(2) = 0d0
    !    if (this%solvent%ionic) then
    !       this%tcfLongRangeAsympK(1) = 0d0
    !       this%tcfLongRangeAsympK(2) = 0d0
    !    end if
    !    ig0 = 2
    ! else
    !    ig0 = 1
    ! end if

    do iu = 1, this%solute%numAtoms
       kforceSum = 0d0
       do igz = 0, this%grid%globalDimsR(3) - 1
          rz = igz * this%grid%voxelVectorsK(3, :)
          do igy = 0, this%grid%globalDimsR(2) - 1
             ry = igy * this%grid%voxelVectorsK(2, :)
             do igx = 0, this%grid%globalDimsR(1) - 1
                rx = igx * this%grid%voxelVectorsK(1, :)
                ig = 1 + igx + igy * this%grid%globalDimsR(1) + &
                     igz * this%grid%globalDimsR(2) * this%grid%globalDimsR(1)
                ! k = 0 is skipped.
                ! if (ig == 1) cycle
                ! if (ig < 3) cycle

                gridPoint = rx + ry + rz
                solutePosition = gridPoint - this%solute%position(:, iu)
                ! solutePosition = minimumImage(this, solutePosition)

                ! dUc_dk = 0
                kforce = 0d0

                ! k = 0 is skipped.
                do igk = 2, this%grid%totalLocalPointsK / 2
                   dUc_dk = FOURPI * exp(-smear2_4 * this%grid%waveVectors2(igk))
                   dUc_dk = dUc_dk / this%grid%waveVectors2(igk)
                   phase = dot_product(this%grid%waveVectors(:, igk), solutePosition)
                   dUc_dk = dUc_dk * sin(phase)
                   kforce = kforce + dUc_dk * this%grid%waveVectors(:, igk)
                end do

                totalSolventCharge = 0

                ! Calculate total solvent charge at grid point.
                do iv = 1, this%solvent%numAtomTypes
                   totalSolventCharge = totalSolventCharge &
                        + this%solvent%density(iv) &
                        * guv(ig, iv) * this%solvent%charge(iv)
                end do
                
                kforceSum = kforceSum + kforce * totalSolventCharge
                
#ifdef __PGI
                if (ieee_is_nan(this%grid%waveVectors2(ig))) &
#else
                if (isnan(this%grid%waveVectors2(ig))) &
#endif
                     print *, ig, this%grid%totalLocalPointsK, this%grid%waveVectors2(ig)
             end do
          end do
       end do
       kforceSum = kforceSum * this%grid%voxelVolume * this%solute%charge(iu) &
            / this%grid%boxVolume
       ! Convert force from reciprocal space to real space.
       kforceSum = 1 / kforceSum
       ff(:, iu) = kforceSum
    end do

    !TODO: Add in short-range Ewald sum by calling ewaldSumShortRangeForce().
  end subroutine ewaldSumForce


  subroutine particleMeshEwaldForce(this, ff, guv)
    use bspline
    use constants, only : pi, fourpi
    use FFTW3
    use rism_util, only : r2c_pointer
    implicit none
    type(rism3d_potential), intent(inout) :: this !< potential object.
    _REAL_, intent(inout) :: ff(3, this%solute%numAtoms) !< Force array [kT/A].
    _REAL_, intent(in) :: guv(:,:) !< Site-site pair correlation function.

    ! Order of b-spline interpolation.
    integer, parameter :: splineOrder = 4

    integer :: i, j, k
    integer :: ix, iy, iz, ixyz
    integer :: igx, igy, igz
    integer :: ig, ig_fort, ig_cpp, igk, igk_fort, igk_cpp
    integer :: id
    integer :: iu, iv
    integer :: gridPointsYZ

    integer :: ierr
    
    integer :: gridPoints(splineOrder, 3)

    _REAL_ :: byz, bxyz
    _REAL_ :: weights(splineOrder, 3)
    _REAL_ :: weightDerivs(splineOrder, 3)
    _REAL_, pointer :: bsplineFourierCoeffX(:) => NULL()
    _REAL_, pointer :: bsplineFourierCoeffY(:) => NULL()
    _REAL_, pointer :: bsplineFourierCoeffZ(:) => NULL()
    _REAL_ :: kernel(this%grid%localDimsR(1) &
         * this%grid%localDimsR(2) * (this%grid%localDimsR(3) / 2 + 1))
    _REAL_ :: k2
    _REAL_, pointer :: kxi(:,:) => NULL()
    _REAL_, pointer :: kyi(:,:) => NULL()
    _REAL_, pointer :: kzi(:,:) => NULL()
    _REAL_, pointer :: elecPotentialR(:) => NULL()
    _REAL_, pointer :: elecPotentialK(:) => NULL()

    _REAL_ :: potGrad(3), weightYZ(3)

    type(C_PTR) ::  planfwd = C_NULL_PTR, planbwd = C_NULL_PTR
    _REAL_, pointer :: inr(:), outr(:)
    complex(kind(1d0)), pointer :: outk(:), ink(:)

    _REAL_ :: waveVector(3)

    _REAL_ :: reciprocalPos(3)

    ! Allocate heap memory.
    kxi => safemem_realloc(kxi,  this%grid%localDimsR(1), 3, .false.)
    kyi => safemem_realloc(kyi,  this%grid%localDimsR(2), 3, .false.)
    kzi => safemem_realloc(kzi,  this%grid%localDimsR(3) / 2 + 1, 3, .false.)

    elecPotentialR => safemem_realloc(elecPotentialR, this%grid%localDimsR(1) &
         * this%grid%localDimsR(2) * (this%grid%localDimsR(3) + 2), &
         o_preserve = .false., o_aligned = .true.)
    elecPotentialK => safemem_realloc(elecPotentialK, this%grid%localDimsR(1) &
         * this%grid%localDimsR(2) * (this%grid%localDimsR(3) + 2), &
         o_preserve = .false., o_aligned = .true.)

    bsplineFourierCoeffX => safemem_realloc( &
         bsplineFourierCoeffX, this%grid%localDimsR(1), .false.)
    bsplineFourierCoeffY => safemem_realloc( &
         bsplineFourierCoeffY, this%grid%localDimsR(2), .false.)
    bsplineFourierCoeffZ => safemem_realloc( &
         bsplineFourierCoeffZ, this%grid%localDimsR(3) / 2 + 1, .false.)

    ! Angular wave vectors.
    do igk = 0, this%grid%localDimsR(1) - 1
       kxi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(1,:) &
            * merge(igk, -(this%grid%localDimsR(1) - igk), &
            igk < this%grid%localDimsR(1) / 2 + 1)
    end do
    do igk = 0, this%grid%localDimsR(2) - 1
       kyi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(2,:) &
            * merge(igk, -(this%grid%localDimsR(2) - igk), &
            igk < this%grid%localDimsR(2) / 2 + 1)
    end do
    do igk = 0, this%grid%localDimsR(3) / 2 + 1 - 1
       kzi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(3,:) * igk
    end do

    ! Compute the discrete Fourier transform coefficients of the b-spline.
    call cardinal_bspline(merge(0d0, 0.5d0, mod(splineOrder, 2) == 0), &
         splineOrder, weights(:,1))
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%localDimsR(1), splineOrder, &
         weights(:,1), bsplineFourierCoeffX, .false.)
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%localDimsR(2), splineOrder, &
         weights(:,1), bsplineFourierCoeffY, .false.)
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%localDimsR(3), splineOrder, &
         weights(:,1), bsplineFourierCoeffZ, .true.)

    ! Reciprocal space kernel with b-spline discrete Fourier transform
    ! correction.
    ixyz = 1
    do ix = 0, this%grid%localDimsR(1) - 1
       do iy = 0, this%grid%localDimsR(2) - 1
          byz = bsplineFourierCoeffX(ix + 1) * bsplineFourierCoeffY(iy + 1)
          do iz = 0, this%grid%localDimsR(3) / 2 + 1 - 1
             bxyz = bsplineFourierCoeffZ(iz + 1) * byz
             waveVector = kxi(ix+1,:) + kyi(iy+1,:) + kzi(iz+1,:)
             k2 = dot_product(waveVector, waveVector)
             kernel(ixyz) = (4 * pi / k2) / bxyz
             ! Note that bxyz is not squared here like it would normally be.
             ixyz = ixyz + 1
          end do
       end do
    end do
    ! Remove the k = 0 term (tinfoil boundary conditions).
    kernel(1) = 0d0

    ! Initialize grids.
    elecPotentialR = 0
    elecPotentialK = 0

    ! Calculate the solvent charge at each grid point.
    !  dac: stored in elecPotentialR(:)
    ! Since solvent density was already solved on a grid by the
    ! 3D-RISM, no need for smearing.
    do iz = 0, this%grid%localDimsR(3) - 1
       do iy = 0, this%grid%localDimsR(2) - 1
          do ix = 0, this%grid%localDimsR(1) - 1
             ig_fort = 1 + ix + (iy + iz * this%grid%localDimsR(2)) &
                  * this%grid%localDimsR(1)
             ig_cpp = 1 + iz + (iy + ix * this%grid%localDimsR(2)) &
                  * this%grid%localDimsR(3)
             do iv = 1, this%solvent%numAtomTypes
                elecPotentialR(ig_cpp) = elecPotentialR(ig_cpp) &
                     + guv(ig_fort, iv) * this%solvent%density(iv) &
                     * this%grid%voxelVolume * this%solvent%charge(iv)
             end do
          end do
       end do
    end do

    ! Allocate the FFT.
    inr => elecPotentialR
    outk => r2c_pointer(elecPotentialK)
    ink => r2c_pointer(elecPotentialK)
    outr => elecPotentialR

    planfwd = fftw_plan_dft_r2c_3d( &
         this%grid%globalDimsR(1), &
         this%grid%globalDimsR(2), &
         this%grid%globalDimsR(3), &
         inr, outk, FFTW_ESTIMATE)

    planbwd = fftw_plan_dft_c2r_3d( &
         this%grid%globalDimsR(1), &
         this%grid%globalDimsR(2), &
         this%grid%globalDimsR(3), &
         ink, outr, FFTW_ESTIMATE)

    ! Compute the reciprocal space charge distribution.
    call fftw_execute_dft_r2c(planfwd, inr, outk)
    
    ! In reciprocal space, convolute the charge density with the
    ! Green function.
    do igk = 1, this%grid%localDimsR(1) * this%grid%localDimsR(2) &
         * (this%grid%localDimsR(3) / 2 + 1)
       outk(igk) = outk(igk) * kernel(igk)
    end do

    call fftw_execute_dft_c2r(planbwd, ink, outr)
    elecPotentialR = elecPotentialR / this%grid%boxVolume
    
    ! Interpolate the gradient of the electric potential from the grid
    ! onto each solute atom.
    ! b-spline gradient + interpolation.

    do iu = 1, this%solute%numAtoms
       ! Convert Cartesian position to reciprocal space by projecting
       ! to reciprocal unit cell vectors.
       do id = 1, 3
          reciprocalPos(id) = dot_product(this%solute%position(:,iu), &
               this%grid%unitCellVectorsK(id, :))
       end do

       ! Set the box length to unity since reciprocal space already
       ! divides positions by the box length.
       !TODO: These could be saved from potential calculation.
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(1), 1d0, this%grid%localDimsR(1), &
            splineOrder, gridPoints(:,1), weights(:,1), weightDerivs(:,1))
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(2), 1d0, this%grid%localDimsR(2), &
            splineOrder, gridPoints(:,2), weights(:,2), weightDerivs(:,2))
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(3), 1d0, this%grid%localDimsR(3), &
            splineOrder, gridPoints(:,3), weights(:,3), weightDerivs(:,3))

       potGrad = 0
       !TODO: Switch index loop order.
       do i = 1, splineOrder
          do j = 1, splineOrder
             gridPointsYZ = (gridPoints(j,2) + gridPoints(i,1) &
                  * this%grid%localDimsR(2)) * this%grid%localDimsR(3)
             weightYZ(1) = weightDerivs(i,1) * weights(j,2)
             weightYZ(2) = weights(i,1) * weightDerivs(j,2)
             weightYZ(3) = weights(i,1) * weights(j,2)
             
             do k = 1, splineOrder
                potGrad(1) = potGrad(1) + weights(k,3) * weightYZ(1) &
                     * elecPotentialR(1 + gridPoints(k,3) + gridPointsYZ)
                potGrad(2) = potGrad(2) + weights(k,3) * weightYZ(2) &
                     * elecPotentialR(1 + gridPoints(k,3) + gridPointsYZ)
                potGrad(3) = potGrad(3) + weightDerivs(k,3) * weightYZ(3) &
                     * elecPotentialR(1 + gridPoints(k,3) + gridPointsYZ)
             end do
          end do
       end do
       ! Apply chain rule from converting b-spline derivatives from
       ! reciprocal space to real space positions.
       potGrad = this%grid%localDimsR * this%solute%charge(iu) * potGrad
       do id = 1, 3
          ff(id, iu) = ff(id, iu) - &
               dot_product(this%grid%unitCellVectorsK(:,id), potGrad)
       end do
    end do

    ! No short-range correction is needed since the solvent is treated
    ! as a smooth distribution on the grid, not as point charges that
    ! need to be smeared out by Gaussians.

    ierr = safemem_dealloc(bsplineFourierCoeffX)
    ierr = safemem_dealloc(bsplineFourierCoeffY)
    ierr = safemem_dealloc(bsplineFourierCoeffZ)
    ierr = safemem_dealloc(kxi)
    ierr = safemem_dealloc(kyi)
    ierr = safemem_dealloc(kzi)
    ierr = safemem_dealloc(elecPotentialR,o_aligned=.true.)
    ierr = safemem_dealloc(elecPotentialK,o_aligned=.true.)

  end subroutine particleMeshEwaldForce
  
  
  !> Calculates both Lennard-Jones and electrostatic solvation forces
  !! over the entire grid without cutoffs.
  subroutine force_brute(this, ff, guv)
    implicit none
    type(rism3d_closure), intent(in) :: this !< Closure object.
    _REAL_, intent(inout) :: ff(3,this%solute%numAtoms) !< Array for forces [kT/A].
    _REAL_, intent(in) :: guv(:,:) !< Site-site pair correlation function.

    integer ::  igx, igy, igz, ig, igoffset, iu, iv, id, idim1, idim2
    _REAL_ ::  rx, ry, rz, dX(3), rs2, rs6i, r5, r4, r3, r2, ra
    _REAL_ :: ffcoef
    _REAL_ ::  rcor, rcor2
    parameter (rcor = 0.002d0, rcor2 = rcor**2)

    _REAL_ ::  sum, dUlj_dr, dUc_dr, dU_dr(3)
    _REAL_ :: qutemp(this%solute%numAtoms), qvtemp(this%solvent%numAtomTypes)
    _REAL_ :: offset
    
    offset = this%grid%spacing(3) * (this%grid%offsetR(3))

    ff = 0d0
    qutemp(1:this%solute%numAtoms) = &
        this%solute%charge(1:this%solute%numAtoms)! * AMBER_ELECTROSTATIC
    qvtemp(1:this%solvent%numAtomTypes) = &
        this%solvent%charge(1:this%solvent%numAtomTypes)! * AMBER_ELECTROSTATIC

    do igz = 0, this%grid%localDimsR(3) - 1
       rz = igz * this%grid%spacing(3) + offset
       do igy = 0, this%grid%localDimsR(2) - 1
          ry = igy * this%grid%spacing(2)
          do igx = 0, this%grid%localDimsR(1) - 1
#ifdef MPI
             ig = 1 + igx + igy * (this%grid%globalDimsR(1) + 2) &
                    + igz * this%grid%globalDimsR(2) &
                          * (this%grid%globalDimsR(1) + 2)
#else
             ig = 1 + igx + igy * this%grid%globalDimsR(1) &
                    + igz * this%grid%globalDimsR(2) &
                          * this%grid%globalDimsR(1)
#endif
             rx = igx * this%grid%spacing(1)

             do iu = 1, this%solute%numAtoms
                dX(1) = this%solute%position(1, iu) - rx
                dX(2) = this%solute%position(2, iu) - ry
                dX(3) = this%solute%position(3, iu) - rz

                r2 = dX(1) * dX(1) + dX(2) * dX(2) + dX(3) * dX(3)
                ra = sqrt(r2)
                r3 = ra * r2
                r4 = r2 * r2
                r5 = ra**5

                if (ra < rcor) then
                   cycle
                end if
                do iv = 1, this%solvent%numAtomTypes
                   rs2 = r2 / this%potential%ljSigmaUV(iu, iv)**2
                   rs6i = 1.d0 / rs2**3
                   dUlj_dr = -12.d0 * this%potential%ljEpsilonuv(iu, iv) &
                        * rs6i * (rs6i - 1.d0) / ra
                   dUc_dr  = -qutemp(iu) * qvtemp(iv) / r2
                   do id = 1, 3
                      dU_dR(id) = (dUlj_dr + dUc_dr) * dX(id) / ra
                      ff(id, iu) = ff(id, iu) &
                          - dU_dR(id) * this%solvent%density(iv) * guv(ig, iv)
                   end do
                end do
             end do
          end do
       end do
    end do
    ff = ff * this%grid%voxelVolume
    return
  end subroutine force_brute
end module rism3d_closure_c
