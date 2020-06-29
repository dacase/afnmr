!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2011-2012 by 
!Andriy Kovalenko, Tyler Luchko and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
!ibid. 112:10391-10417 (2000).   
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Partial series expansion of order-n (PSE-n) closure class for 3D-RISM.
!!!Stefan M. Kast and Thomas Kloss. J. Chem. Phys. 129, 236101 (2008)
!!!
!!!Interpolates between the Kovalenko-Hirata (KH) and hypernetted chain equation
!!!closures (HNC).  Order-1 gives KH and order-infinity gives HNC.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module rism3d_psen_c
    use rism3d_potential_c
    use rism3d_grid_c
    !the PSEN type
    type rism3d_psen
       type(rism3d_potential),pointer :: pot => NULL()
       !grid : points to grid in potential object
       type(rism3d_grid),pointer :: grid => NULL()
       !order    : order of the series expansion.  >= 1
       !order1   : order+1
       integer :: order=0, order1
       !order1fac: (order+1)!
       _REAL_ :: order1fac
    end type rism3d_psen

    public rism3d_psen_new, rism3d_psen_destroy!, rism3d_psen_guv
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the PSEN closure
!!!IN:
!!!   this : PSEN object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_psen_new(this,pot, order)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      type(rism3d_potential), target, intent(in) :: pot
      integer, intent(in) :: order
      integer :: i, orderfac
      this%pot => pot
      this%grid => this%pot%grid
      if(order <1)then
         call rism_report_error('(a,i4)',"PSE-n closure: order must be > 0:",order)
      end if
      this%order = order
      this%order1 = order+1
      this%order1fac = 1
      orderfac=1

      !calculate the factorial of (order+1) for future reference.  Check that we 
      !are not losing precision or overflowing the precision
      do i=2,this%order1
         orderfac = this%order1fac
         this%order1fac = this%order1fac*i
         !keep 16 significant digits
         if( abs((this%order1fac/i)/orderfac -1d0) >1d-16 )then
            call rism_report_error('(a,i4)',&
                 "PSE-n: numerical precision exceeded.  Use a smaller order. Try ",i-1)
         end if
      end do
    end subroutine rism3d_psen_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv from Uuv, Huv, and Cuv using the PSEN closure
!!!IN:
!!!   this : the PSEN closure object
!!!   guv  : site-site pair correlation function
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_psen_guv(this,guv, huv, cuv)
      implicit none
      type(rism3d_psen), intent(in) :: this
      _REAL_, intent(out) :: guv(:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      integer :: i, iv, ir, ix, iy, iz, ig, ig1
      _REAL_ :: tuv, orderfac

      do iv = 1,this%pot%solvent%numAtomTypes
         do iz = 1, this%grid%localDimsR(3)
            do iy = 1, this%grid%localDimsR(2)
               do ix = 1, this%grid%localDimsR(1)
                  ig1 = ix + (iy - 1) * this%grid%localDimsR(1) + &
                       (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#if defined(MPI)
                  ig = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) &
                       + (iz - 1) * (this%grid%localDimsR(1) + 2) * this%grid%localDimsR(2)
#else
                  ig = ix + (iy - 1) * this%grid%localDimsR(1) + &
                       (iz - 1) * this%grid%localDimsR(1) * this%grid%localDimsR(2)
#endif /*defined(MPI)*/
                  tuv = -this%pot%uuv(ix,iy,iz,iv) + huv(ig,iv) - cuv(ix,iy,iz,iv)
                  if(tuv >= 0d0)then
                     guv(ig,iv) = 1d0
                     orderfac = 1d0
                     do i=1,this%order
                        orderfac = orderfac*i
                        guv(ig,iv) = guv(ig,iv) + (tuv**i)/orderfac
                     end do
                  else
                     guv(ig,iv) = exp(tuv)
                  endif
               end do
            end do
         end do
      end do
    end subroutine rism3d_psen_guv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv_dT from Uuv, Huv_dT, Cuv_dT and Guv using the PSEN closure
!!!IN:
!!!   this    : the PSEN closure object
!!!   guv_dT  : temperature derivative site-site pair correlation function
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   guv     : site-site pair correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_psen_guv_dT(this,guv_dT, huv_dT, cuv_dT, guv, huv, cuv)
      implicit none
      type(rism3d_psen), intent(in) :: this
      _REAL_, intent(out) :: guv_dT(:,:)
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:), guv(:,:), huv(:,:),&
           cuv(:,:,:,:)
      integer :: i, iv, ir, ix, iy, iz, ig, ig1
      !tuv :: t* = beta*u + h - c
      !tuv :: dT t* = beta*u + dT h - dT c
      _REAL_ ::  tuv, tuv_dT
      _REAL_ :: orderfac

      do iv = 1,this%pot%solvent%numAtomTypes
         do iz = 1, this%grid%localDimsR(3)
            do iy = 1, this%grid%localDimsR(2)
               do ix = 1, this%grid%localDimsR(1)
                  ig1 = ix + (iy - 1) * this%grid%localDimsR(1) + &
                       (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#if defined(MPI)
                  ig = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) + &
                       (iz - 1) * (this%grid%localDimsR(1) + 2) * this%grid%localDimsR(2)
#else
                  ig = ix + (iy - 1) * this%grid%localDimsR(1) + &
                       (iz - 1) * this%grid%localDimsR(1) * this%grid%localDimsR(2)
#endif /*defined(MPI)*/

                  tuv_dT = this%pot%uuv(ix,iy,iz,iv) + huv_dT(ig,iv) - cuv_dT(ix,iy,iz,iv)
                  if(guv(ig,iv) < 1d0)then
                     guv_dT(ig,iv) = guv(ig,iv)*tuv_dT
                  else
                     tuv = -this%pot%uuv(ix,iy,iz,iv) + huv(ig,iv) - cuv(ix,iy,iz,iv)
                     guv_dT(ig,iv)=1d0
                     orderfac = 1d0
                     do i=1,this%order-1
                        orderfac = orderfac*i
                        guv_dT(ig,iv) = guv_dT(ig,iv) + (tuv**i)/orderfac
                     end do
                     guv_dT(ig,iv) = guv_dT(ig,iv)*tuv_dT
                  endif
               end do
            end do
         end do
      end do
    end subroutine rism3d_psen_guv_dT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_psen_excessChemicalPotential(this, huv, cuv) result(excessChemicalPotential)
      implicit none
      type(rism3d_psen), intent(in) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: excessChemicalPotential(this%pot%solvent%numAtomTypes)
      _REAL_ :: tuv, tsuv
      integer :: ix, iy, iz, iv, igk
      excessChemicalPotential = 0.d0
      do iv=1,this%pot%solvent%numAtomTypes
         do iz=1,this%grid%localDimsR(3)
            do iy=1,this%grid%localDimsR(2)
               do ix=1,this%grid%localDimsR(1)
#if defined(MPI)
                  igk = ix + (iy-1)*(this%grid%localDimsR(1)+2) + (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
                  igk = ix + (iy-1)*this%grid%localDimsR(1) + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#endif /*defined(MPI)*/
                  tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                  excessChemicalPotential(iv) = excessChemicalPotential(iv) + 0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv)
                  if (huv(igk,iv) > 0d0)  then
                     tsuv = tuv - this%pot%uuv(ix,iy,iz,iv)
                     excessChemicalPotential(iv) = excessChemicalPotential(iv) - tsuv**(this%order1)/this%order1fac
                  endif
               end do
            end do
         end do
         excessChemicalPotential(iv) =  this%pot%solvent%density(iv)*&
              excessChemicalPotential(iv)*this%grid%voxelVolume
      enddo
    end function rism3d_psen_excessChemicalPotential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT at the requested grid point
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!   ijk  : 3d-grid index
!!!OUT:
!!!   The contribution to the excess chemical potential from grid point i, j, k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_psen_excessChemicalPotential_IJK(this, huv, cuv,ijk,iv) result(excessChemicalPotential)
      implicit none
      type(rism3d_psen), intent(in) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      integer, intent(in) :: ijk(3), iv
      _REAL_ :: excessChemicalPotential
      _REAL_ :: tuv, tsuv
      integer :: ix, iy, iz, igk
      ix=ijk(1)
      iy=ijk(2)
      iz=ijk(3)
#if defined(MPI)
      igk = ix + (iy-1)*(this%grid%localDimsR(1)+2) + (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
      igk = ix + (iy-1)*this%grid%localDimsR(1) + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#endif /*defined(MPI)*/
      tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
      excessChemicalPotential = 0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv)
      if (huv(igk,iv) > 0d0)  then
         tsuv = tuv - this%pot%uuv(ix,iy,iz,iv)
         excessChemicalPotential = excessChemicalPotential - tsuv**(this%order1)/this%order1fac
      endif
      excessChemicalPotential = excessChemicalPotential *this%pot%solvent%density(iv)
    end function rism3d_psen_excessChemicalPotential_IJK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the short range part of the excess chemical potential w/ 
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_psen_aexcessChemicalPotential(this, huv, cuv) result(excessChemicalPotential)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: excessChemicalPotential(this%pot%solvent%numAtomTypes)
      !excessChemicalPotentialh2lr :: contribution from tcfLongRangeAsympR**2
      !excessChemicalPotentialhclr :: contribution from tcfLongRangeAsympR*dcfLongRangeAsympR
      !excessChemicalPotentialhnlr :: contribution from tcfLongRangeAsympR**n      
      _REAL_ :: excessChemicalPotentialh2lr(this%pot%solvent%numAtomTypes), &
           excessChemicalPotentialhclr(this%pot%solvent%numAtomTypes), &
           excessChemicalPotentialhnlr(this%pot%solvent%numAtomTypes)
      _REAL_ :: tuv, tsuv, huvlr, cuvlr, tuvlr
      integer :: ix, iy, iz, iv, ig, igk

      if(.not.this%pot%solute%charged)then
         excessChemicalPotential = rism3d_psen_excessChemicalPotential(this,huv,cuv)
         return
      end if


      !
      !long range part
      !
      call rism3d_potential_int_h2_hc(this%pot,excessChemicalPotentialh2lr,excessChemicalPotentialhclr)
      excessChemicalPotentialhnlr = rism3d_potential_int_hn(this%pot,this%order1)
      
      do iv=1,this%pot%solvent%numAtomTypes
         if(this%pot%solute%totalCharge*this%pot%solvent%charge_sp(iv) > 0.d0) &
              excessChemicalPotentialhnlr(iv) = 0
      end do

      !
      !short range part
      !
      excessChemicalPotential = 0.d0
      huvlr=0d0
      do iv=1,this%pot%solvent%numAtomTypes
         do iz=1,this%grid%localDimsR(3)
            do iy=1,this%grid%localDimsR(2)
               do ix=1,this%grid%localDimsR(1)
                  ig = ix + (iy - 1) * this%grid%localDimsR(1) + &
                       (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#if defined(MPI)
                  igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) + &
                       (iz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
#else
                  igk = ix + (iy - 1) * this%grid%localDimsR(1) + &
                       (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#endif /*defined(MPI)*/
                  tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                  if(this%pot%solvent%ionic)&
                       huvlr = this%pot%solvent%charge_sp(iv)*this%pot%tcfLongRangeAsympR(ig)
                  cuvlr = this%pot%solvent%charge(iv)*this%pot%dcfLongRangeAsympR(ig)
                  tuvlr = huvlr - cuvlr
                  excessChemicalPotential(iv) = excessChemicalPotential(iv) + 0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv)
!!$                  excessChemicalPotential(iv) = excessChemicalPotential(iv) - (0.5d0*huvlr*tuvlr - cuvlr)
                  excessChemicalPotential(iv) = excessChemicalPotential(iv) - (0.5d0*huvlr*tuvlr)
                  if (huv(igk,iv) > 0d0)  then
                     tsuv = tuv - this%pot%uuv(ix,iy,iz,iv)
                     excessChemicalPotential(iv) = excessChemicalPotential(iv) - tsuv**(this%order1)/this%order1fac
                  endif

                  if (this%pot%solute%totalCharge*this%pot%solvent%charge_sp(iv) <= 0.d0) then
                     !recall that cuvlr = - uuvlr so 
                     !(huvlr - cuvlr - uuvlr)**n = huvlr**n
                     excessChemicalPotential(iv) = excessChemicalPotential(iv)  + huvlr**(this%order1)/this%order1fac
                  endif
               end do
            end do
         end do
         excessChemicalPotential(iv) =  this%pot%solvent%density(iv) &
              * (excessChemicalPotential(iv) * this%grid%voxelVolume &
              +  (excessChemicalPotentialh2lr(iv) - excessChemicalPotentialhclr(iv)) / 2d0 &
              -  excessChemicalPotentialhnlr(iv) / this%order1fac)
      enddo

    end function rism3d_psen_aexcessChemicalPotential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the temperature derivative w/o
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this    : the PSEN closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!OUT:
!!!    solvation energy without long-range correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_psen_solvationEnergy(this, huv_dT, cuv_dT, huv, cuv) result(solvationEnergy)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: solvationEnergy(this%pot%solvent%numAtomTypes)
      !tuv :: t* = beta*u + h - c
      !tuv :: dT t* = beta*u + dT h - dT c
      _REAL_ :: tuv, tuv_dT
      !orderfac :: order-n factorial
      _REAL_ orderfac
      integer :: ix, iy, iz, iv, igk
      orderfac = this%order1fac/this%order1
      solvationEnergy = 0.d0
      do iv=1,this%pot%solvent%numAtomTypes
         do iz=1,this%grid%localDimsR(3)
            do iy=1,this%grid%localDimsR(2)
               do ix=1,this%grid%localDimsR(1)
#if defined(MPI)
                  igk = ix + (iy-1)*(this%grid%localDimsR(1)+2)&
                       + (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
                  igk = ix + (iy-1)*this%grid%localDimsR(1) &
                       + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#endif /*defined(MPI)*/
                  solvationEnergy(iv) = solvationEnergy(iv) + huv(igk,iv)*huv_dT(igk,iv) &
                       - (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
                       - 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
                  if (huv(igk,iv) > 0d0)  then
                     tuv = -this%pot%uuv(ix,iy,iz,iv) + huv(igk,iv) - cuv(ix,iy,iz,iv)
                     tuv_dT = this%pot%uuv(ix,iy,iz,iv) + huv_dT(igk,iv) - cuv_dT(ix,iy,iz,iv)
                     solvationEnergy(iv) = solvationEnergy(iv) - tuv**this%order/orderfac * tuv_dT
                  endif
               end do
            end do
         end do
         solvationEnergy(iv) =  -1d0*this%pot%solvent%density(iv)*&
              solvationEnergy(iv)*this%grid%voxelVolume
      enddo
    end function rism3d_psen_solvationEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the
!!!temperature derivative w/o asymptotic correction in kT for the grid
!!!point.
!!!IN:
!!!   this    : the PSEN closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!   ijk  : 3d-grid index
!!!OUT:
!!!    Contribution to the solvation energy from grid point i, j, k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_psen_solvationEnergy_IJK(this, huv_dT, cuv_dT, huv, cuv, ijk, iv) result(solvationEnergy)
      implicit none
      type(rism3d_psen), intent(in) :: this
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      integer, intent(in) :: ijk(3), iv
      _REAL_ :: solvationEnergy
      !tuv :: t* = beta*u + h - c
      !tuv :: dT t* = beta*u + dT h - dT c
      _REAL_ :: tuv, tuv_dT
      !orderfac :: order-n factorial
      _REAL_ orderfac      
      integer :: ix, iy, iz,igk

      !use the precomputed (n+1)! to get n!
      orderfac = this%order1fac/this%order1

      ix=ijk(1)
      iy=ijk(2)
      iz=ijk(3)
#if defined(MPI)
      igk = ix + (iy-1)*(this%grid%localDimsR(1)+2) + (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
      igk = ix + (iy-1)*this%grid%localDimsR(1) + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#endif /*defined(MPI)*/
      solvationEnergy = -huv(igk,iv)*huv_dT(igk,iv) &
           + (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
           + 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
      if (huv(igk,iv) > 0d0)  then
         tuv = -this%pot%uuv(ix,iy,iz,iv) + huv(igk,iv) - cuv(ix,iy,iz,iv)
         tuv_dT = this%pot%uuv(ix,iy,iz,iv) + huv_dT(igk,iv) - cuv_dT(ix,iy,iz,iv)
         solvationEnergy = solvationEnergy + tuv**this%order/orderfac * tuv_dT
      endif
      solvationEnergy=solvationEnergy*this%pot%solvent%density(iv)
    end function rism3d_psen_solvationEnergy_IJK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the temperature derivative w/
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this    : the PSEN closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!OUT:
!!!    excess chemical potential with long-range correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_psen_asolvationEnergy(this, huv_dT, cuv_dT, huv, cuv) result(solvationEnergy)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: solvationEnergy(this%pot%solvent%numAtomTypes)
      _REAL_ :: solvationEnergyh2lr(this%pot%solvent%numAtomTypes),solvationEnergyhclr(this%pot%solvent%numAtomTypes),&
           solvationEnergyhnlr(this%pot%solvent%numAtomTypes)
      !tuv :: t* = beta*u + h - c
      !tuv :: dT t* = beta*u + dT h - dT c
      _REAL_ :: tuv, tuv_dT, huvlr, cuvlr
      !orderfac :: order-n factorial
      _REAL_ orderfac
      integer :: ix, iy, iz, iv, ig, igk

      if(.not.this%pot%solute%charged)then
         solvationEnergy = rism3d_psen_solvationEnergy(this,huv_dT,cuv_dT,huv,cuv)
         return
      end if

      !
      !long range part
      !
      call rism3d_potential_int_h2_hc(this%pot,solvationEnergyh2lr,solvationEnergyhclr)
      solvationEnergyhnlr = rism3d_potential_int_hn(this%pot,this%order1)

      do iv=1,this%pot%solvent%numAtomTypes
         if(this%pot%solute%totalCharge*this%pot%solvent%charge_sp(iv) <= 0.d0) &
              solvationEnergyhnlr(iv)=0
      end do

      !
      !short range part
      !
      solvationEnergy = 0.d0
      huvlr=0d0
      !use the precomputed (n+1)! to get n!
      orderfac = this%order1fac/this%order1
      do iv=1,this%pot%solvent%numAtomTypes
         do iz=1,this%grid%localDimsR(3)
            do iy=1,this%grid%localDimsR(2)
               do ix=1,this%grid%localDimsR(1)
                  ig = ix + (iy - 1) * this%grid%localDimsR(1) + &
                       (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#if defined(MPI)
                  igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) + &
                       (iz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
#else
                  igk = ix + (iy - 1) * this%grid%localDimsR(1) + &
                       (iz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#endif /*defined(MPI)*/
                  if(this%pot%solvent%ionic)&
                       huvlr = this%pot%solvent%charge_sp(iv)*this%pot%tcfLongRangeAsympR(ig)
                  cuvlr = this%pot%solvent%charge(iv)*this%pot%dcfLongRangeAsympR(ig)
                  solvationEnergy(iv) = solvationEnergy(iv) + huv(igk,iv)*huv_dT(igk,iv) &
                       - (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
                       - 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
                  if (huv(igk,iv) > 0d0)  then
                     tuv = -this%pot%uuv(ix,iy,iz,iv) + huv(igk,iv) - cuv(ix,iy,iz,iv)
                     tuv_dT = this%pot%uuv(ix,iy,iz,iv) + huv_dT(igk,iv) - cuv_dT(ix,iy,iz,iv)
                     solvationEnergy(iv) = solvationEnergy(iv) - tuv**this%order/orderfac * tuv_dT
                  endif

                  solvationEnergy(iv) = solvationEnergy(iv) + huvlr*huvlr - huvlr*cuvlr
                  if (this%pot%solute%totalCharge*this%pot%solvent%charge_sp(iv) <= 0.d0) then
                     !recall that cuvlr = - uuvlr so 
                     !(huvlr - cuvlr - uuvlr)**n = huvlr**n
                     solvationEnergy(iv) = solvationEnergy(iv) - huvlr**this%order1/orderfac
                  endif
               end do
            end do
         end do
         solvationEnergy(iv) =  -1d0*(this%pot%solvent%density(iv)*(solvationEnergy(iv)*this%grid%voxelVolume &
              -  (solvationEnergyh2lr(iv) - solvationEnergyhclr(iv) - solvationEnergyhnlr(iv)/orderfac)))
      enddo

    end function rism3d_psen_asolvationEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the PSEN closure
!!!IN:
!!!   this : PSEN object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_psen_destroy(this)
      implicit none
      type(rism3d_psen), intent(inout) :: this
      nullify(this%pot)
      nullify(this%grid)
    end subroutine rism3d_psen_destroy
  end module rism3d_psen_c
