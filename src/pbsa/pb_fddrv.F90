! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"
#include "timer.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver for various finite-difference algorithms
subroutine pb_fddrv( npbstep,npbgrid,nstlim,atmfirst,atmlast,npbopt,solvopt,level,nfocus,bcopt,&
                     h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                     xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                     maxitn,itn,fmiccg,fmiccg2,accept,laccept,wsor,lwsor,inorm,norm,&
                     pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,osmene,&
                     ngrdcrg,grdcrg,qgrdcrg,&
                     gcrd,acrg,&
                     nbnd,iepsav,insas,epsx,epsy,epsz,&
                     chgrd,saltgrd,ioncrg,phi,&
                     bv,cphi,xs&
                   )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Solving A * x + f(x) = b, where A is dielectric (and salt if linear) map,
   ! f(x) is the nonlienar salt map, b is the charge map, x is phi map.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! R. Qi
   use pbtimer_module
   implicit none

#include "pb_constants.h"

   ! all the driver variables are shared among all "contained" routines, so are
   ! not redeclared in the containted routines anymore, except there is a need
   ! to remap their dimensions and to copy to other variables.
    
   ! passed variables

   integer npbstep, npbgrid, nstlim, atmfirst, atmlast
   integer npbopt, solvopt, level, nfocus, bcopt
   _REAL_ h, savh(MAXLEVEL), gox, goy, goz, savgox(MAXLEVEL), savgoy(MAXLEVEL), savgoz(MAXLEVEL)
   integer xm, ym, zm, xmym, xmymzm, savxm(MAXLEVEL), savym(MAXLEVEL), savzm(MAXLEVEL)
   integer maxitn, itn
   _REAL_ fmiccg, fmiccg2, accept, laccept, wsor, lwsor, inorm, norm
   _REAL_ pbkappa, pbkb, pbtemp, ivalence, istrng, eps0, epsin, epsout
   _REAL_ ionene, osmene
   integer ngrdcrg,grdcrg(3,*)
   _REAL_ qgrdcrg(*)
   _REAL_  gcrd(3,atmfirst:atmlast), acrg(atmfirst:atmlast)
   integer nbnd, iepsav(4,xmymzm)
   integer insas(*)
   _REAL_ epsx(xmymzm+ym*zm), epsy(xmymzm+xm*zm), epsz(xmymzm+xm*ym), chgrd(xmymzm)
   _REAL_ saltgrd(xmymzm), ioncrg(xmymzm)
   _REAL_ phi(xmymzm)
   _REAL_ bv(xmymzm), cphi(xmymzm), xs(xmymzm+2*xmym)
    
   ! local variables
    
   _REAL_ rh
   integer i, j, k
 
   rh = ONE/h

   ! prepare the b vector

   ! part a. place induced polarization charges or atomic charges on the grid points and
   ! store in the b vector for the singularity-free PB, whether outside is total field or
   ! reaction field.

   call pbtimer_start(PBTIME_PBBUILDSYS)
   cphi = ZERO
   if ( bcopt == 2 .or. bcopt > 5 .and. bcopt < 10 ) then
      bv = ZERO
      call pb_dbcgrd( chgrd,insas,epsx,epsy,epsz,bv,cphi )
   else
      bv = chgrd*rh
   end if

   ! part b. set the space boundary condition, note except the first level, focusing is
   ! needed. the boundary conditions (potentials) are converted into charges on boundary
   ! grid points and saved in the b vector

   call pb_bndcnd( level,bcopt,chgrd,phi,epsx,epsy,epsz,bv )
   call pbtimer_stop(PBTIME_PBBUILDSYS)

   ! enter the core iteration routine
   ! npbstep and npbgrid should be replaced with pbgrid or pbinitial logical variable

   call pbtimer_start(PBTIME_PBSOLV) ! Timing
   if ( npbopt == 0) then
      call solve_lpb( xm,ym,zm,xmym,xmymzm,maxitn,bcopt, &
                      fmiccg,fmiccg2,accept,pbkappa,epsout,h,wsor, &
                      bv,saltgrd,xs )
   else
      call solve_npb( xm,ym,zm,xmym,xmymzm,itn,maxitn, &
                      npbstep,npbgrid, &
                      inorm,norm,wsor,lwsor, &
                      saltgrd,bv &
                    )
   end if
   call pbtimer_stop(PBTIME_PBSOLV) ! Timing

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Dielectric Boundary
subroutine pb_dbcgrd( chgrd,insas,epsx,epsy,epsz,bv,cphi )
   
   implicit none

   ! passed variables

   _REAL_  chgrd(xm,ym,zm)
   integer insas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_  epsx(0:xm,1:ym,1:zm)
   _REAL_  epsy(1:xm,0:ym,1:zm)
   _REAL_  epsz(1:xm,1:ym,0:zm)
   _REAL_  bv(xm,ym,zm)
   _REAL_  cphi(xm,ym,zm)

   ! local variables

   integer i, j, k
   integer i0, j0, k0
   integer ip
   _REAL_ tmp, tmp0, epst
   
   ! collecting grid charges into working arrays
 
   ngrdcrg = 0
   do k = 1, zm; do j = 1, ym; do i = 1, xm
      if ( chgrd(i,j,k) == ZERO ) cycle
      ngrdcrg = ngrdcrg + 1

      grdcrg(1,ngrdcrg) = i
      grdcrg(2,ngrdcrg) = j
      grdcrg(3,ngrdcrg) = k
      qgrdcrg(ngrdcrg) = chgrd(i,j,k)
   end do; end do; end do

   tmp = ZERO; tmp0 = ZERO
   do ip = 1, nbnd
      i0 = iepsav(1,ip); j0 = iepsav(2,ip); k0 = iepsav(3,ip)
      epst = ZERO         
   
      i = i0 - 1; epst = epst + epsx(i ,j0,k0)
      if ( insas(i ,j0,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i ,j0,k0) == ZERO ) then
            call get_coulpot(i ,j0,k0,tmp)
            cphi(i ,j0,k0) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsx(i ,j0,k0)*cphi(i ,j0,k0)
      end if

      i = i0 + 1; epst = epst + epsx(i0,j0,k0)
      if ( insas(i ,j0,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i ,j0,k0) == ZERO ) then
            call get_coulpot(i ,j0,k0,tmp)
            cphi(i ,j0,k0) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsx(i0,j0,k0)*cphi(i ,j0,k0)
      end if
   
      j = j0 - 1; epst = epst + epsy(i0,j ,k0)
      if ( insas(i0,j ,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j ,k0) == ZERO ) then
            call get_coulpot(i0,j ,k0,tmp); 
            cphi(i0,j ,k0) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsy(i0,j ,k0)*cphi(i0,j ,k0)
      end if

      j = j0 + 1; epst = epst + epsy(i0,j0,k0)
      if ( insas(i0,j ,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j ,k0) == ZERO ) then
            call get_coulpot(i0,j ,k0,tmp)
            cphi(i0,j ,k0) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsy(i0,j0,k0)*cphi(i0,j ,k0)
      end if

      k = k0 - 1; epst = epst + epsz(i0,j0,k )
      if ( insas(i0,j0,k ) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j0,k ) == ZERO ) then
            call get_coulpot(i0,j0,k ,tmp)
            cphi(i0,j0,k ) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsz(i0,j0,k )*cphi(i0,j0,k )
      end if

      k = k0 + 1; epst = epst + epsz(i0,j0,k0)
      if ( insas(i0,j0,k ) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j0,k ) == ZERO ) then
            call get_coulpot(i0,j0,k ,tmp)
            cphi(i0,j0,k ) = tmp/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) + epsz(i0,j0,k0)*cphi(i0,j0,k )
      end if

      if ( insas(i0,j0,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j0,k0) == ZERO ) then
            call get_coulpot(i0,j0,k0,tmp0)
            cphi(i0,j0,k0) = tmp0/epsin
         end if
         bv(i0,j0,k0) = bv(i0,j0,k0) - epst*cphi(i0,j0,k0)
      end if

      bv(i0,j0,k0) = bv(i0,j0,k0) + chgrd(i0,j0,k0)*rh
   end do

end subroutine pb_dbcgrd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ get coulomb potential
subroutine get_coulpot(i,j,k,pot)

   implicit none

   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green

   integer i,j,k
   _REAL_ pot

   integer iatm
   integer itmp,jtmp,ktmp
   integer idx,idy,idz
   _REAL_ qtmp,rinv

   pot = ZERO
   do iatm = 1, ngrdcrg
      itmp = grdcrg(1,iatm); jtmp = grdcrg(2,iatm); ktmp = grdcrg(3,iatm)
      qtmp = qgrdcrg(iatm)
           
      idx = abs(i-itmp); idy = abs(j-jtmp); idz = abs(k-ktmp)
      if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
         rinv = green(idx,idy,idz)
         pot = pot + qtmp*rinv
      else
         rinv = ONE/sqrt(dble(idx**2 + idy**2 + idz**2))
         pot = pot + qtmp*rinv
      end if
   end do  !  iatm = 1, ngrdcrg
   pot = pot*INV_FOURPI*rh

end subroutine get_coulpot
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign Debye-Huckel potential for the boundary charge grid
subroutine pb_bndcnd( level,bcopt,chgrd,phi,epsx,epsy,epsz,bv )
   
   implicit none

   ! Common variables
    
   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green
     
   ! Passed variables

   integer level, bcopt
   _REAL_ chgrd(xm,ym,zm)
   _REAL_ bv(xm,ym,zm)
   _REAL_ phi(xm,ym,zm)
   _REAL_ epsx(0:xm,1:ym,1:zm)
   _REAL_ epsy(1:xm,0:ym,1:zm)
   _REAL_ epsz(1:xm,1:ym,0:zm)
    
   ! Local variables
    
   integer i, j, k, iatm
   integer xmtmp, ymtmp, zmtmp, ix, iy, iz, itmp, jtmp, ktmp, idx, idy, idz
   _REAL_ htmp, goxtmp, goytmp, goztmp
   _REAL_ qtmp, xtmp, ytmp, ztmp, dx, dy, dz
   _REAL_ x, y, z
   _REAL_ xi, yi, zi, aa, bb, cc, aa1, bb1, cc1
   _REAL_ r, rinv
   _REAL_ factor
   
   ! part a: level = 1 cases ::::: 
   ! bcopt = 1
   ! conductor boundary: use zero potential for boundary grid
   ! the boundary will be all solvent, so do nothing below ...
   ! to solve the singular PB
    
   if ( level == 1 .and. bcopt == 1 ) then
    
   ! bcopt = 2
   ! conductor boundary: use zero potential for boundary grid
   ! the boundary will be all solvent, so do nothing below ...
   ! to solve the singularity-free PB
    
   else if ( level == 1 .and. bcopt == 2 ) then
    
   ! bcopt = 3
   ! sum of residue dipolar debye-huckel contribution. the boundary will be all solvent.
    
   else if ( level == 1 .and. bcopt == 3 ) then
      write(6, *) "PB bomb in pb_bndcnd(): residue dipolar BC not supported"
      call mexit(6, 1)
    
   ! bcopt = 4 or 5 or 6 or 7
   ! to solve the singular (4 and 5) or singularity-free (6 and 7) PB
   ! sum of atom charge contributions (4 or 6).
   ! sum of grid charge contributions (5 or 7).
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. ( bcopt == 4 .or. bcopt == 6 ) ) then

      do iatm = atmfirst, atmlast
         ! h unit, which is just what is needed in bv as it is in e/h
         xtmp = gcrd(1,iatm); ytmp = gcrd(2,iatm); ztmp = gcrd(3,iatm)
         qtmp = INV_FOURPI*acrg(iatm)

         ! k=0 and k=zm+1 faces

         do j = 1, ym; do i = 1, xm
            dx = i-xtmp; dy = j-ytmp; dz = ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h*r))*qtmp/r

            dz = zm+1-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h*r))*qtmp/r
         end do; end do

         ! j=0 and ym+1 faces

         do k = 1, zm; do i = 1, xm
            dx = i-xtmp; dy  = ytmp; dz  = k-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h*r))*qtmp/r

            dy = ym+1-ytmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h*r))*qtmp/r
         end do; end do

         ! i=0 and i=xm+1 faces

         do k = 1, zm; do j = 1, ym
            dx = xtmp; dy = j-ytmp; dz = k-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h*r))*qtmp/r

            dx = xm+1-xtmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h*r))*qtmp/r
         end do; end do
      end do  !  iatm = atmfirst, atmlast
       
   else if ( level == 1 .and. ( bcopt == 5 .or. bcopt == 7 ) ) then

      factor = INV_FOURPI*rh
      do itmp = 1, xm; do jtmp = 1, ym; do ktmp = 1, zm
         qtmp = chgrd(itmp,jtmp,ktmp)
         if ( qtmp == ZERO ) cycle
         qtmp = factor*qtmp
          
         ! k=0 and k=zm+1 faces

         do j = 1, ym; do i = 1, xm
            idx = abs(i-itmp); idy = abs(j-jtmp); idz = ktmp
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h*r))*qtmp/r
            end if
             
            idz = abs(zm+1-ktmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do
           
         ! j=0 and ym+1 faces
           
         do k = 1, zm; do i = 1, xm
            idx = abs(i-itmp); idy  = jtmp; idz  = abs(k-ktmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
              
            idy = abs(ym+1-jtmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do
      
         ! i=0 and i=xm+1 faces
      
         do k = 1, zm; do j = 1, ym
            idx = itmp; idy = abs(j-jtmp); idz = abs(k-ktmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
            
            idx = abs(xm+1-itmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do

      end do; end do; end do  !  itmp = 1, xm; jtmp = 1, ym; ktmp = 1, zm

   ! bcopt = 8 or 9
   ! to solve the singularity-free PB.
   ! sum of grid charge reaction field contribution (9). the boundary will be all solvent.
   ! sum of atom charge contributions (8) is silently treated as sum of grid charge
   ! contributions
    
   else if ( level == 1 .and. ( bcopt == 8 .or. bcopt == 9 ) ) then
       
      factor = INV_FOURPI * rh * epsout * (ONE/epsout - ONE/epsin) 
      do itmp = 1, xm; do jtmp = 1, ym; do ktmp = 1, zm
         qtmp = chgrd(itmp,jtmp,ktmp)
         if ( qtmp == ZERO ) cycle
         qtmp = factor*qtmp
          
         ! k=0 and k=zm+1 faces

         do j = 1, ym; do i = 1, xm
            idx = abs(i-itmp); idy = abs(j-jtmp); idz = ktmp
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,j,1 ) = bv(i,j,1 ) + exp(pbkappa*(-h*r))*qtmp/r
            end if
             
            idz = abs(zm+1-ktmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,j,zm) = bv(i,j,zm) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do
           
         ! j=0 and ym+1 faces
           
         do k = 1, zm; do i = 1, xm
            idx = abs(i-itmp); idy  = jtmp; idz  = abs(k-ktmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,1 ,k) = bv(i,1 ,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
              
            idy = abs(ym+1-jtmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,ym,k) = bv(i,ym,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do
      
         ! i=0 and i=xm+1 faces
      
         do k = 1, zm; do j = 1, ym
            idx = itmp; idy = abs(j-jtmp); idz = abs(k-ktmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(1 ,j,k) = bv(1 ,j,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
            
            idx = abs(xm+1-itmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h/rinv))*qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(xm,j,k) = bv(xm,j,k) + exp(pbkappa*(-h*r))*qtmp/r
            end if
         end do; end do

      end do; end do; end do  !  itmp = 1, xm; jtmp = 1, ym; ktmp = 1, zm

   ! bcopt = 10
   ! periodic boundary
   ! nothing needs to be done here just like the conductor case, bcopt=1 
   ! the linear system has been changed to impose the PBC

   else if ( level == 1 .and. bcopt == 10 ) then

   ! part b. level > 1 cases
   ! electrostatic focusing
    
   else if ( level > 1 ) then
      xmtmp  = savxm(level-1) ; ymtmp  = savym(level-1) ; zmtmp  = savzm(level-1)
      htmp   = savh(level-1)
      goxtmp = savgox(level-1); goytmp = savgoy(level-1); goztmp = savgoz(level-1)
       
      ! k=0 and k=zm+1 faces
       
      do j = 1, ym; do i = 1, xm
          
         x  = gox + h*i        ; y  = goy + h*j        ; z  = goz
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa  = xi - dble( ix ); bb  = yi - dble( iy ); cc  = zi - dble( iz )
         aa1 = ONE - aa       ; bb1 = ONE - bb       ; cc1 = ONE - cc
         bv(i,j,1 ) = bv(i,j,1 ) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         z  = goz + h*(zm+1)
         zi = (z - goztmp)/htmp
         iz = int( zi )
         cc  = zi - dble( iz )
         cc1 = ONE - cc
         bv(i,j,zm) = bv(i,j,zm) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
   
      ! j=0 and j=ym+1 faces
   
      do k = 1, zm; do i = 1, xm
          
         x  = gox + h*i        ; y  = goy              ; z  = goz + h*k
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa = xi - dble( ix ); bb = yi - dble( iy ); cc = zi - dble( iz )
         aa1 = ONE - aa      ; bb1 = ONE - bb      ; cc1 = ONE - cc
         bv(i,1 ,k) = bv(i,1 ,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         y  = goy + h*(ym+1)
         yi = (y - goytmp)/htmp
         iy = int( yi )
         bb  = yi - dble( iy )
         bb1 = ONE - bb
         bv(i,ym,k) = bv(i,ym,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
        
      ! i=0 and i=xm+1 faces
        
      do k = 1, zm; do j = 1, ym
          
         x  = gox              ; y  = goy + h*j        ; z  = goz + h*k
         xi = (x - goxtmp)/htmp; yi = (y - goytmp)/htmp; zi = (z - goztmp)/htmp
         ix = int( xi )        ; iy = int( yi )        ; iz = int( zi )
         aa  = xi - dble( ix ); bb  = yi - dble( iy ); cc  = zi - dble( iz )
         aa1 = ONE - aa       ; bb1 = ONE - bb       ; cc1 = ONE - cc
         bv(1 ,j,k) = bv(1 ,j,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
         x  = gox + h * (xm+1)
         xi = (x - goxtmp)/htmp
         ix = int( xi )
         aa  = xi - dble( ix )
         aa1 = ONE - aa
         bv(xm,j,k) = bv(xm,j,k) + epsout*phintp( xmtmp, ymtmp, zmtmp, ix, iy, iz, aa, bb, cc, aa1, bb1, cc1 )
          
      end do; end do
       
   else
       
      ! unknown combinations
       
      write(6, *) 'PB bomb in pb_bndcnd(): unknown BC option'
      call mexit(6, 1)
   end if 

end subroutine pb_bndcnd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ phi interpretation, xm, ym, zm are for the previous phi map
_REAL_ function phintp( xmtmp,ymtmp,zmtmp,ix,iy,iz,aa,bb,cc,aa1,bb1,cc1 )
  
   implicit none
  
   ! Passed variables
    
   integer, intent(in) :: xmtmp, ymtmp, zmtmp, ix, iy, iz
   _REAL_, intent(in) :: aa, bb, cc, aa1, bb1, cc1
    
   ! Local Variables
    
   _REAL_ bb1cc1, bb_cc1, bb1cc, bb_cc
    
   ! determine the position of the point w.r.t. the map
    
   bb1cc1 = bb1*cc1; bb_cc1 = bb *cc1; bb1cc  = bb1*cc ; bb_cc  = bb *cc

   ! triliner interpolation
    
   phintp = aa1*bb1cc1*phi( ix   + xmtmp*( iy-1 + ymtmp*( iz-1 ) ) ) + &
            aa *bb1cc1*phi( ix+1 + xmtmp*( iy-1 + ymtmp*( iz-1 ) ) ) + &
            aa1*bb_cc1*phi( ix   + xmtmp*( iy   + ymtmp*( iz-1 ) ) ) + &
            aa *bb_cc1*phi( ix+1 + xmtmp*( iy   + ymtmp*( iz-1 ) ) ) + &
            aa1*bb1cc *phi( ix   + xmtmp*( iy-1 + ymtmp*( iz   ) ) ) + &
            aa *bb1cc *phi( ix+1 + xmtmp*( iy-1 + ymtmp*( iz   ) ) ) + &
            aa1*bb_cc *phi( ix   + xmtmp*( iy   + ymtmp*( iz   ) ) ) + &
            aa *bb_cc *phi( ix+1 + xmtmp*( iy   + ymtmp*( iz   ) ) )


end function phintp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ driver of linear finite-difference algorithms
subroutine solve_lpb( nx,ny,nz,nxny,nxnynz,p_maxitn,p_bcopt, &
                      p_fmiccg,p_fmiccg2,p_accept,p_pbkappa,p_epsout,p_h,p_wsor, &
                      p_bv,p_iv,p_xs )

   use pb_lsolver
   implicit none

   integer nx,ny,nz,nxny,nxnynz,p_maxitn,p_bcopt
   _REAL_ p_fmiccg,p_fmiccg2,p_accept,p_epsout,p_pbkappa,p_h,p_wsor
   _REAL_ p_bv(1:nxnynz),p_iv(1:nxnynz)
   _REAL_ p_xs(*)

#ifdef CUDA
   ! For single precision
   real , dimension(1:nxnynz) :: s_phi
   real , dimension(1:nxnynz) :: s_p_xs
   ! For returning functions
   integer get_itn
   real get_inorm, get_norm ! C function returning a value should be called as F fucntion, otherwise as F subroutine
   external get_itn, get_inorm, get_norm

   s_phi = phi(1:nxnynz)
   s_p_xs = p_xs(1:nxnynz)
 
   if ( solvopt == 2 .or. solvopt == 4 ) then
      call init_param_c(nx, ny, nz, p_maxitn, p_bcopt, &
              real(p_accept), real(p_pbkappa), real(p_epsout), real(p_h), real(p_wsor))
      call allocate_array_cuda(solvopt)
      call init_array_cuda(solvopt, real(epsx), real(epsy), real(epsz), real(p_bv), real(p_iv), real(p_xs(1:nxnynz)))
   else if ( solvopt == 3) then
      call init_param( nx,ny,nz,nxny,nxnynz,&
              p_maxitn,p_bcopt,p_fmiccg,p_fmiccg2,p_accept,p_pbkappa,p_epsout,p_h,p_wsor)
      call allocate_array(solvopt)
      call init_array(solvopt,epsx,epsy,epsz,p_bv,p_iv,p_xs)
   else
      write(6, *) 'PB Bomb in pb_fddrv(): unknown solver'
      call mexit(6, 1)
   end if
#else
   call init_param( nx,ny,nz,nxny,nxnynz,&
           p_maxitn,p_bcopt,p_fmiccg,p_fmiccg2,p_accept,p_pbkappa,p_epsout,p_h,p_wsor)
   call allocate_array(solvopt)
   call init_array(solvopt,epsx,epsy,epsz,p_bv,p_iv,p_xs)
#endif

#ifdef CUDA
   select case ( solvopt )
   case (2)
      ! p_xs index from 1, not 1-nxny, since only 1:l_xmymzm of passed xs will be used
      call pb_mg_cuda(s_phi, s_p_xs(1))
   case (3)
      if ( l_bcopt == 10 ) then
      call pb_pcg(phi,p_xs)
      else
      call pb_cg(phi,p_xs)
      end if
   case (4)
      ! p_xs index from 1, not 1-nxny, since only 1:l_xmymzm of passed xs will be used
      call pb_sor_cuda(s_phi, s_p_xs(1))
   case default
      write(6, *) 'PB Bomb in pb_fddrv(): unknown solver'
      call mexit(6, 1)
   end select

   ! return back
   if ( solvopt == 2 .or. solvopt == 4 ) then
      phi(1:nxnynz) = s_phi(1:nxnynz)
      p_xs(1:nxnynz) = s_p_xs(1:nxnynz)
      itn = get_itn()
      inorm = get_inorm()
      norm = get_norm()

      call deallocate_array_cuda()
   else if ( solvopt == 3) then
      itn = l_itn
      inorm = l_inorm
      norm = l_norm

      call deallocate_array(solvopt)
   else
      write(6, *) 'PB Bomb in pb_fddrv(): unknown solver'
      call mexit(6, 1)
   end if


#else
   select case ( solvopt )
   case (1)
      if ( l_bcopt == 10 ) then
         call pb_piccg(phi,p_xs(1:(nx+2)*(ny+2)*(nz+2)))
      else
         call pb_iccg(phi,p_xs)
      end if
   case (2)
      call pb_mg(phi,p_xs)
   case (3)
      if ( l_bcopt == 10 ) then
         call pb_pcg(phi,p_xs)
      else
         call pb_cg(phi,p_xs)
      end if
   case (4)
      if ( l_bcopt == 10 ) then
         call pb_psor(phi,p_xs)
      else
         call pb_sor(phi,p_xs)
      end if
   case default
      write(6, *) 'PB Bomb in pb_fddrv(): unknown solver'
      call mexit(6, 1)
   end select

   ! return back
   itn = l_itn
   inorm = l_inorm
   norm = l_norm

   call deallocate_array(solvopt)
#endif

end subroutine solve_lpb
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ driver of nonlinear finite-difference algorithms
subroutine solve_npb(nx,ny,nz,nxny,nxnynz,p_itn,p_maxitn, &
                     p_npbstep,p_npbgrid, &
                     p_inorm,p_norm,p_wsor,p_lwsor, &
                     p_iv,p_bv &
                    )

   use pb_nlsolver

   implicit none

   integer nx,ny,nz,nxny,nxnynz
   integer p_itn,p_maxitn,p_npbstep,p_npbgrid
   _REAL_ p_inorm,p_norm,p_wsor,p_lwsor
   _REAL_ p_iv(nxnynz),p_bv(nxnynz)

   npbstep = p_npbstep

   call init_param(nx,ny,nz,nxny,nxnynz,p_maxitn,p_npbstep,p_npbgrid, &
                   ivalence,h,pbkb,pbtemp,istrng,p_wsor,p_lwsor)
   call allocate_array(solvopt)
   call init_array(xs,epsx,epsy,epsz,p_iv,p_bv, &
                   solvopt,npbopt,h,epsout,eps0,pbtemp,pbkappa,nbnd,iepsav )
   
   select case ( solvopt )
   case (1)
      call pb_nticcg( phi, xs, p_bv(1), accept, npbopt, nbnd, iepsav )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case (2) 
      call pb_nmg( phi, xs, p_bv(1), epsout, accept, npbopt, nbnd, iepsav )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case (3)
      if ( npbopt == 1 .and. ( npbstep == 0 .or. mod(npbstep+1,npbgrid) == 0 )) then
         npbopt = 0
         call pb_ncg( phi, xs, laccept, npbopt )
         bv(1:xmymzm) = p_bv(1:xmymzm) 
         npbopt = 1
      endif
      call pb_ncg( phi, xs, accept, npbopt )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case (4)
      if ( npbopt == 1 .and. ( npbstep == 0 .or. mod(npbstep+1,npbgrid) == 0 ) ) then
         npbopt = 0
         call pb_nsor( phi, xs, lwsor, laccept, npbopt )
         npbopt = 1
      endif
      call pb_nsor( phi, xs, wsor, accept, npbopt )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case (5)
      if ( npbopt == 1 .and. ( npbstep == 0 .or. mod(npbstep+1,npbgrid) == 0 ) ) then
         npbopt = 0
         call pb_asor( phi, xs, lwsor, laccept, npbopt )
         npbopt = 1
      endif
      call pb_asor( phi, xs, wsor, accept, npbopt )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case (6)
      if ( npbopt == 1 .and. ( npbstep == 0 .or. mod(npbstep+1,npbgrid) == 0 ) ) then
         npbopt = 0
         call pb_dsor( phi, xs, lwsor, laccept, npbopt )
         npbopt = 1
      endif
      call pb_dsor( phi, xs, wsor, accept, npbopt )
      phi(1:xmymzm) = xs(1+xmym:xmymzm+xmym)
   case default
      write(6, *) 'PB bomb in pb_fddrv(): unknown solver'
      call mexit(6, 1)
   end select

   p_itn = itn
   p_inorm = inorm
   p_norm = norm

   !  calculate ionic energy

   if ( npbopt /= 0 ) then
      call pb_ionene(nfocus,level,savgox,savgoy,savgoz, &
         savxm,savym,savzm,savh,phi,ioncrg,ionene,osmene)
   end if

   call deallocate_array(solvopt)

end subroutine solve_npb

end subroutine pb_fddrv
