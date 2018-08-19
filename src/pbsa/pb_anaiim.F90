! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"
#include "timer.h"

module ana_iim

#include "pb_constants.h"

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ IIM algoritm interface routine
subroutine pb_anaiim( npbstep,npbgrid,nstlim,atmfirst,atmlast,npbopt,solvopt,level,nfocus,bcopt,&
                      natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                      xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                      maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
                      pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
                      gcrd,acrg,&
                      nbnd,iepsav,insas,lvlset,&
                      chgrd,saltgrd,phi,&
                      bv,cphi,xs&
                    )

   use pbtimer_module
   implicit none

   ! all the driver variables are shared among all "contained" routines, so are
   ! not redeclared in the containted routines anymore, except there is a need
   ! to remap their dimensions and to copy to other variables.

   ! passed variables

   integer npbstep, npbgrid, nstlim, atmfirst, atmlast,natom
   integer npbopt, solvopt, level, nfocus, bcopt
   _REAL_  h, savh(nfocus), gox, goy, goz, savgox(nfocus), savgoy(nfocus), savgoz(nfocus)
   integer xm, ym, zm, xmym, xmymzm, savxm(nfocus), savym(nfocus), savzm(nfocus)
   integer maxitn, itn
   _REAL_  fmiccg, accept, laccept, wsor, lwsor, inorm, norm
   _REAL_  pbkappa, pbkb, pbtemp, ivalence, istrng, eps0, epsin, epsout, ionene
   _REAL_  gcrd(3,natom), acrg(natom)
   integer nbnd, iepsav(4,xmymzm)

   integer insas(xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8)
   _REAL_  lvlset(xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8)

   _REAL_  chgrd(xmymzm), saltgrd(xmymzm), phi(xmymzm)
   _REAL_  bv(xmymzm), cphi(xmymzm), xs(xmymzm+2*xmym)
   integer ii

   ! local varialbes

   _REAL_ rh

   rh = ONE/h

   ! for singular PB equation bv() is initialized as the grid charges

   call pbtimer_start(PBTIME_PBBUILDSYS)
   if ( bcopt /= 2 .and. bcopt < 6 ) then
      bv = chgrd
   end if

   ! put effective charges at the space boundary grid points into bv()
   ! to take care of the boundary conditions

   call pb_bndcnd( bv(1), chgrd(1) )
   call pbtimer_stop(PBTIME_PBBUILDSYS)

   call iim(gox,goy,goz,xm,ym,zm,lvlset,insas,nbnd,iepsav, &
            epsin/eps0,epsout/eps0, &
            bv(1),phi(1),xs(1),h,atmfirst,atmlast,bcopt,&
            maxitn,itn,accept,norm,inorm)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign Debye-Huckel potential for the boundary charge grid
subroutine pb_bndcnd( bv, chgrd )

   implicit none
   
   ! Common variables
    
   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green
     
   ! Passed variables
    
   _REAL_ bv(xm,ym,zm), chgrd(xm,ym,zm)
    
   ! Local variables
    
   integer i, j, k, iatm, ii
   integer xmtmp, ymtmp, zmtmp, ix, iy, iz, itmp, jtmp, ktmp, idx, idy, idz
   _REAL_ htmp, goxtmp, goytmp, goztmp
   _REAL_ qtmp, factor
   _REAL_ x, y, z, dx, dy, dz, xtmp, ytmp, ztmp
   _REAL_ xi, yi, zi, aa, bb, cc, aa1, bb1, cc1
   _REAL_ r, rinv

   ! part a: level = 1 cases ::::: 
   ! bcopt = 1
   ! zero potential in the singular PB.
   ! the boundary will be all solvent

   if ( level == 1 .and. bcopt == 1 ) then
    
   ! bcopt = 2
   ! zero potential in the singularty-free PB.
   ! the boundary will be all solvent
    
   else if ( level == 1 .and. bcopt == 2 ) then
    
   ! bcopt = 3
   ! sum of residue dipolar debye-huckel potentials in the singular PB.
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. bcopt == 3 ) then
      write(6, *) "PB bomb in pb_iimdrv(): residue dipolar BC not supported"
      call mexit(6, 1)
    
   ! bcopt = 4 .or. bcopt = 6
   ! sum of atom charge deby-huckel potentials in the singular (4) or singularity-free (6) PB.
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. ( bcopt == 4 .or. bcopt == 6 ) ) then

      do iatm = atmfirst, atmlast
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
    
   ! bcopt = 5 .or. bcopt = 7
   ! sum of grid charge debye-huckel potentials in the singular (5) or singularity-free (7) PB.
   ! the boundary will be all solvent.

   else if ( level == 1 .and. (bcopt == 5 .or. bcopt == 7) ) then
       
      do itmp = 1, xm; do jtmp = 1, ym; do ktmp = 1, zm
         qtmp = chgrd(itmp,jtmp,ktmp)
         if ( qtmp == ZERO ) cycle

         qtmp = INV_FOURPI*qtmp 

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
       
   ! bcopt = 8
   ! sum of atom charge reaction field potentials in the singularity-free PB.
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. bcopt == 8 ) then
       
      factor = INV_FOURPI*epsout*(ONE/epsout - ONE/epsin)
      do iatm = atmfirst, atmlast
         xtmp = gcrd(1,iatm); ytmp = gcrd(2,iatm); ztmp = gcrd(3,iatm)
         qtmp = factor*acrg(iatm)
          
         ! k=0 and k=zm+1 faces
          
         do j = 1, ym; do i = 1, xm
            dx = i-xtmp; dy = j-ytmp; dz = ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,j,1 ) = bv(i,j,1 ) + qtmp/r
             
            dz = zm+1-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,j,zm) = bv(i,j,zm) + qtmp/r
         end do; end do
          
         ! j=0 and ym+1 faces
          
         do k = 1, zm; do i = 1, xm
            dx = i-xtmp; dy  = ytmp; dz  = k-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,1 ,k) = bv(i,1 ,k) + qtmp/r
             
            dy = ym+1-ytmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(i,ym,k) = bv(i,ym,k) + qtmp/r
         end do; end do
          
         ! i=0 and i=xm+1 faces
          
         do k = 1, zm; do j = 1, ym
            dx = xtmp; dy = j-ytmp; dz = k-ztmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(1 ,j,k) = bv(1 ,j,k) + qtmp/r
             
            dx = xm+1-xtmp
            r = sqrt(dx**2 + dy**2 + dz**2)
            bv(xm,j,k) = bv(xm,j,k) + qtmp/r
         end do; end do
      end do  !  iatm = atmfirst, atmlast
    
   ! bcopt = 9
   ! sum of grid charge reaction field potentials in the singularity-free PB.
   ! the boundary will be all solvent.
    
   else if ( level == 1 .and. bcopt == 9 ) then
       
      factor = INV_FOURPI*epsout*(ONE/epsout - ONE/epsin)
      do itmp = 1, xm; do jtmp = 1, ym; do ktmp = 1, zm
         qtmp = chgrd(itmp,jtmp,ktmp)
         if ( qtmp == ZERO ) cycle

         qtmp = factor*qtmp
          
         ! k=0 and k=zm+1 faces
          
         do j = 1, ym; do i = 1, xm
            idx = abs(i-itmp); idy = abs(j-jtmp); idz = ktmp
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,j,1 ) = bv(i,j,1 ) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,j,1 ) = bv(i,j,1 ) + qtmp/r
            end if
             
            idz = abs(zm+1-ktmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,j,zm) = bv(i,j,zm) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,j,zm) = bv(i,j,zm) + qtmp/r
            end if
         end do; end do
           
         ! j=0 and ym+1 faces
           
         do k = 1, zm; do i = 1, xm
            idx = abs(i-itmp); idy  = jtmp; idz  = abs(k-ktmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,1 ,k) = bv(i,1 ,k) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,1 ,k) = bv(i,1 ,k) + qtmp/r
            end if
              
            idy = abs(ym+1-jtmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(i,ym,k) = bv(i,ym,k) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(i,ym,k) = bv(i,ym,k) + qtmp/r
            end if
         end do; end do
      
         ! i=0 and i=xm+1 faces
      
         do k = 1, zm; do j = 1, ym
            idx = itmp; idy = abs(j-jtmp); idz = abs(k-ktmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(1 ,j,k) = bv(1 ,j,k) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(1 ,j,k) = bv(1 ,j,k) + qtmp/r
            end if
            
            idx = abs(xm+1-itmp)
            if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
               rinv = green(idx,idy,idz)
               bv(xm,j,k) = bv(xm,j,k) + qtmp*rinv
            else
               r = sqrt(dble(idx**2 + idy**2 + idz**2))
               bv(xm,j,k) = bv(xm,j,k) + qtmp/r
            end if
         end do; end do

      end do; end do; end do  !  itmp = 1, xm; jtmp = 1, ym; ktmp = 1, zm

   ! part b: level > 1 case 
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
       
      ! unknown bcopt
       
      write(6, *) 'PB bomb in pb_iimdrv(): unknown BC option', bcopt
      call mexit(6, 1)
   end if 
 
end subroutine pb_bndcnd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ phi interpretation, xm, ym, zm are for the previous phi map
_REAL_ function phintp( xmtmp,ymtmp,zmtmp,ix,iy,iz,aa,bb,cc,aa1,bb1,cc1 )
    
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

end subroutine pb_anaiim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ This is the real IIM driver 
subroutine iim(gox,goy,goz,l,m,n,lvlset,insas,nbnd,iepsav,epsin,epsout,bv,u,u0,h,&
               atmfirst,atmlast,bcopt,maxitn,itn,accept,norm,inorm)

!    xs <- gox, ys <- goy, zs <- goz
!    l  <- xm , m  <- ym , n  <- zm
!    bi <- epsin, bo <- epsout
!    u  <- phi, u0 <- xs

   use pbtimer_module
   use density_surface, only : index, index2, x, y, z, cirreg
   implicit none

   ! passed variables
      
   _REAL_  gox,goy,goz,h,epsout,epsin
   integer l,m,n,nbnd,atmfirst,atmlast,bcopt
   integer maxitn,itn
   _REAL_ accept,norm,inorm
   integer insas(0:l+1,0:m+1,0:n+1),iepsav(1:4,1:l*m*n)
   _REAL_  bv(l,m,n),u0(l,m,n),u(l,m,n)
   _REAL_  lvlset(0:l+1,0:m+1,0:n+1)

   ! local variables

   _REAL_, parameter :: eps0 = 8.8542D-12 / (1.6022D-19)**2 /  &
                               (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer, parameter :: nq=27
   integer i,j,k,IFAIL,IFAIL2,ii,ir,nc,nz_num
   _REAL_ bi,bo
   _REAL_ beta_max,rhs
   _REAL_ coe2(27)
   _REAL_ err
   
   _REAL_, allocatable :: wp (:), qp (:)
   _REAL_, allocatable :: qyp(:), qzp(:)
   _REAL_, allocatable :: wyp(:), wzp(:)
   _REAL_, allocatable :: wyyp(:), wzzp(:), wyzp(:)

   _REAL_, allocatable :: c(:,:,:,:), c2(:, :)
   _REAL_, allocatable :: f(:,:,:)

   ! interfacing variables 

   bi = epsin; bo = epsout

   ! allocating working arrays

   call pbtimer_start(PBTIME_PBBUILDSYS)
   allocate (wp(nbnd), qp(nbnd))  ! jump conditions at irregular points
   allocate (qyp(nbnd),qzp(nbnd)) ! tangential derivatives of field jump conditions
   allocate (wyp(nbnd),wzp(nbnd)) ! tangential derivatives of potential jump conditions
   allocate (wyyp(nbnd),wzzp(nbnd),wyzp(nbnd)) ! second derivatives of potential jump conditions
   allocate (c(l,m,n,7), c2(nbnd, 27)) ! A matrix coefficients
   allocate (f(l,m,n)) ! rhs B vector

   ! the surface related infrastructure data from density_surface module is used

   ! setting up jump conditions and their derivatives at the projection points

   call jumps(l,m,n,nbnd,bcopt,bi,bo,atmfirst,atmlast,cirreg,&
               wp,wyp,wzp,wyyp,wzzp,wyzp,qp,qyp,qzp)

   ! setting up linear system coefficient matrix

   beta_max = max(bi,bo)
   nz_num = 0
   do k = 1, n
   do j = 1, m
   do i = 1, l
      if ( index(i,j,k) == 1 ) then
         do ii = 2, 7
            c(i,j,k,ii) = epsin/h/h
         end do
         if ( i == 1 ) c(i,j,k,2) = 0.d0
         if ( i == l ) c(i,j,k,3) = 0.d0
         if ( j == 1 ) c(i,j,k,4) = 0.d0
         if ( j == m ) c(i,j,k,5) = 0.d0
         if ( k == 1 ) c(i,j,k,6) = 0.d0
         if ( k == n ) c(i,j,k,7) = 0.d0
         c(i,j,k,1) = 6.d0*epsin/h/h
         f(i,j,k)   = bv(i,j,k)/h/h/h*FOURPI
         do ii = 1, 7
            if ( abs(c(i,j,k,ii)) > 1.d-10 ) then
               nz_num=nz_num+1
            end if
         end do
      else if (index(i,j,k) == 5 ) then
         do ii = 2, 7
            c(i,j,k,ii) = epsout/h/h
         end do
         if ( i == 1 ) c(i,j,k,2) = 0.d0
         if ( i == l ) c(i,j,k,3) = 0.d0
         if ( j == 1 ) c(i,j,k,4) = 0.d0
         if ( j == m ) c(i,j,k,5) = 0.d0
         if ( k == 1 ) c(i,j,k,6) = 0.d0
         if ( k == n ) c(i,j,k,7) = 0.d0
         c(i,j,k,1) = 6.d0*epsout/h/h
         f(i,j,k) = bv(i,j,k)/h/h/h*FOURPI
         do ii = 1, 7
            if ( c(i,j,k,ii ) > 1.d-10 ) then
               nz_num=nz_num+1
            end if
         end do
      else

         call irre31(l,m,n,h,h,h,h,IFAIL, &
                     i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,lvlset,index, &
                     nq,nbnd,index2,cirreg,coe2,rhs, &
                     wp,wyp,wzp,wyyp,wzzp,wyzp,qp,qyp,qzp)
         if (IFAIL.gt.10) &
         call irre32(l,m,n,h,h,h,h,IFAIL2, &
                     i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,lvlset,index, &
                     nq,nbnd,index2,cirreg,coe2,rhs, &
                     wp,wyp,wzp,wyyp,wzzp,wyzp,qp,qyp,qzp)
         ir = index2(i,j,k)
         do ii = 1, 27
            c2(ir,ii) = coe2(ii)
         end do
         f(i,j,k) = rhs
         do ii = 1, 27
            if ( abs(c2(ir,ii)) > 1.d-10 ) then
               nz_num=nz_num+1
            end if
         end do
      endif
   end do
   end do
   end do
   ! cleaning up working arrays
   ! TODO: for future interface field interpolation, the jump conditions
   ! should be saved
   
   deallocate (wp, qp)  ! jump conditions at irregular points
   deallocate (qyp, qzp) ! tangential derivatives of field jump conditions
   deallocate (wyp, wzp) ! tangential derivatives of potential jump conditions
   deallocate (wyyp, wzzp, wyzp) ! second derivatives of potential jump conditions
   call pbtimer_stop(PBTIME_PBBUILDSYS)

   ! entering the linear system solver
   ! Only allows ILU-preconditioned BiCG

   call pbtimer_start(PBTIME_PBSOLV)
   call bicg(l,m,n,nbnd,nz_num,c,c2,index,index2,f,u,u0,maxitn,itn,accept,err)
   call pbtimer_stop(PBTIME_PBSOLV)

   ! return error related data

   inorm = sum(abs(f))
   norm = inorm*err

   ! converting back to the PBSA unit for potential
 
   do k =1, n
   do j =1, m
   do i =1, l
      u(i,j,k) = u(i,j,k) * INV_FOURPI / eps0
   end do
   end do
   end do
   
   deallocate (c, c2)
   deallocate (f)

end subroutine iim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine irre31 here]
subroutine irre31(l,m,n,h,hx,hy,hz,ifail,i0,j0,k0,info,bmax,b_in,b_out,x,y,z,phi,index, &
              nq,nirreg,index2,cirreg,xx,rhs, &
              wp,wyp,wzp,wyyp,wzzp,wyzp,qp,qyp,qzp)

   ! set up the linear equation for an irregular point (i0,j0,k0)
   ! determin the coefficients in matrix A and the rhs to satisfy
   ! the specified jump conditions

   implicit none

   ! parameter settings 

   integer im,ime,immax,ns,in,inmax,imnn,lwar,liwar
   parameter(im=36, ime=10, immax=im, ns=27, in=ns, &
         inmax=in, imnn=im+2*in, &
         lwar=3*inmax*inmax/2+10*inmax+2*immax+1, &
         liwar=in)
 
   ! external functions
  
   !_REAL_, external :: fmin,fb_in,fb_out,fk_in,fk_out,ff_in,ff_out

   ! interface variables
   
   integer ifail,i0,j0,k0,info,nq,nirreg
   integer index(l,m,n),index2(l,m,n)
   _REAL_  bmax,b_out,b_in
   _REAL_  rhs
   _REAL_  x(0:l+1), y(0:m+1), z(0:n+1)
   _REAL_  phi(0:l+1,0:m+1, 0:n+1)
   _REAL_  cirreg(15, nirreg)
   _REAL_  xx(in)
   _REAL_  wp(nirreg), &
           wyp(nirreg), wzp(nirreg), wyyp(nirreg), wzzp(nirreg), wyzp(nirreg), &
           qp(nirreg), qyp(nirreg), qzp(nirreg)

   ! local variables
   
   integer i,j,nc,kr,kctr,ndis,l,m,n,k,ii,ir,iout,iprint,ji
   integer iwar(liwar) 
   
   _REAL_  wy,wz,wyy,wzz,wyz,q,qy,qz,fkk_in,fkk_out &
          ,hxx,hyy,hzz,h1,eps,rho,rho1,bl2jmp,bl3jmp &
          ,hxyz,tmp,tmp1,fkk,corr &
          ,fkjmp,fjmp,uump,tmp2,temp,tmp3,x1,y1,z1  &
          ,xyy,xzz,xyz,h,hx,hy,hz,hmax,ujmp
   _REAL_  t(3,3),a(20)
   _REAL_  xloc(3),tmpx(3)
   _REAL_  bl_out(3),bl_in(3)
   _REAL_  cc(inmax,inmax), ca(immax,inmax), caa(immax,inmax)
   _REAL_  cd(inmax),cb(immax),xl(in),xu(in),uu(imnn)
   _REAL_  war(lwar)
   _REAL_  tmpcoe(7)

   ! retrieving info for the irregular point

   ir = index2(i0, j0, k0)
 
   x1     = cirreg( 1, ir)
   y1     = cirreg( 2, ir)
   z1     = cirreg( 3, ir)

   xyy    = cirreg( 4, ir)
   xzz    = cirreg( 5, ir)
   xyz    = cirreg( 6, ir)

   t(1,1) = cirreg( 7, ir)
   t(2,1) = cirreg( 8, ir)
   t(3,1) = cirreg( 9, ir)
   t(1,2) = cirreg(10, ir)
   t(2,2) = cirreg(11, ir)
   t(3,2) = cirreg(12, ir)
   t(1,3) = cirreg(13, ir)
   t(2,3) = cirreg(14, ir)
   t(3,3) = cirreg(15, ir)

   ! this is the jump of the rhs ...

   !call fjmps(x1,y1,z1, fjmp)
   fjmp = ZERO

   ! this is the jump of the linear term ...

   !call fkjmps(x1,y1,z1, fkjmp)
   fkjmp = ZERO

   ! retrieve the jump in u, its first derivatives,
   ! and second derivatives

   ujmp = wp(ir)
   wy   = wyp(ir)
   wz   = wzp(ir)
   wyy  = wyyp(ir)
   wzz  = wzzp(ir)
   wyz  = wyzp(ir)

   ! retrieve the jump in beta u_n, its first derivative,
   ! and second derivatives

   q =  qp(ir)
   qy = qyp(ir)
   qz = qzp(ir)
 
   ! compute the derivatives of beta of inside and out
   ! and transform them into local coordinate system

   !call betas(hx,hy,hz,x1,y1,z1,t,bl_in,bl_out,b_in,b_out)
   bl_in = ZERO
   bl_out = ZERO

   ! the linear term inside and outside

   fkk_in = ZERO !fk_in(x1,y1,z1)
   fkk_out = ZERO !fk_out(x1,y1,z1)

   hxx=hx*hx
   hyy=hy*hy
   hzz=hz*hz
   h1 = min(hxx,hyy,hzz)

   eps=0.1e-11

   ! set up ca() and cb() for the nearest grid points
   ! Jun: what are these?
   ! ca: distance/distance products
   ! cb: beta/beta derivatives

   do i=1,immax
      do j=1,inmax
         ca(i,j) = 0.0
      end do
      cb(i) = 0.0
   end do

   nc = 0
   kr = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1
            kr = kr +1
            if (i == i0 .and. j == j0 .and. k == k0) then
               kr = kr -1
               kctr = nc
            else
               ca(kr+ime,nc) = 1.0d0
               cb(kr+ime) = -1.0d-20
            end if
            xl(nc) = 0.0
            xu(nc) = 2.0*bmax/h1
         end do
      end do
   end do

   xl(kctr)=-12.0*bmax/h1
   xu(kctr)=0.0

   rho = b_in/b_out
   rho1 = 1 - rho
   bl2jmp = bl_out(2)-bl_in(2)
   bl3jmp = bl_out(3)-bl_in(3)
   temp = fkjmp/b_out

   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1

            ! ----------- local coordinates transformation

            tmpx(1) = x(i) - x1
            tmpx(2) = y(j) - y1
            tmpx(3) = z(k) - z1
            call matvec(3,3,t, tmpx, xloc)

            ! ----------- form the coefficients of equality constraints

            if (index(i,j,k) <= 3) then
               ca(1,nc)  = 1.0
               ca(2,nc)  = xloc(1)
               ca(3,nc)  = xloc(2)
               ca(4,nc)  = xloc(3)
               ca(5,nc)  = 0.5*xloc(1)*xloc(1)
               ca(6,nc)  = 0.5*xloc(2)*xloc(2)
               ca(7,nc)  = 0.5*xloc(3)*xloc(3)
               ca(8,nc)  = xloc(1)*xloc(2)
               ca(9,nc)  = xloc(1)*xloc(3)
               ca(10,nc) = xloc(2)*xloc(3)
            else
               ca(1,nc)  = 1.0-0.5*temp*xloc(1)*xloc(1)
               ca(2,nc)  = rho*xloc(1) &
                     - 0.5*xloc(1)*xloc(1)*(rho1*(xyy+xzz) &
                     - bl_in(1)/b_out + rho*bl_out(1)/b_out) &
                     + 0.5*xloc(2)*xloc(2)*xyy*rho1 &
                     + 0.5*xloc(3)*xloc(3)*xzz*rho1 &
                     + xloc(1)*xloc(2)*(bl_in(2)/b_out &
                     - rho*bl_out(2)/b_out) &
                     + xloc(1)*xloc(3)*(bl_in(3)/b_out &
                     - rho*bl_out(3)/b_out) &
                     + xloc(2)*xloc(3)*xyz*rho1
               ca(3,nc)  = - 0.5*xloc(1)*xloc(1)*bl2jmp/b_out &
                     + xloc(1)*xloc(2)*xyy*rho1 &
                     + xloc(1)*xloc(3)*xyz*rho1 &
                     + xloc(2)
               ca(4,nc)  = - 0.5*xloc(1)*xloc(1)*bl3jmp/b_out &
                     + xloc(1)*xloc(3)*xzz*rho1 &
                     + xloc(1)*xloc(2)*xyz*rho1 &
                     + xloc(3)
               ca(5,nc)  = 0.5*xloc(1)*xloc(1)*rho
               ca(6,nc)  = 0.5*xloc(2)*xloc(2) &
                     + 0.5*(rho-1)*xloc(1)*xloc(1)
               ca(7,nc)  = 0.5*xloc(3)*xloc(3) &
                     + 0.5*(rho-1)*xloc(1)*xloc(1)
               ca(8,nc)  = xloc(1)*xloc(2)*rho
               ca(9,nc)  = xloc(1)*xloc(3)*rho
               ca(10,nc) = xloc(2)*xloc(3)
            end if  ! (index(i,j,k) <= 3)
         end do  !  k=k0-1,k0+1
      end do  !  j=j0-1,j0+1
   end do  !  i=i0-1,i0+1
 
   if (info > 3) then
      cb(1) = fkjmp
   end if

   cb(2)  = -bl_in(1)
   cb(3)  = -bl_in(2)
   cb(4)  = -bl_in(3)
   cb(5)  = -b_in
   cb(6)  = -b_in
   cb(7)  = -b_in

   do i=1,inmax
      do j=1,inmax
         cc(i,j) = 0.0
      end do
      cc(i,i) = 1.0
   end do
   hxyz=hx*hx+hy*hy+hz*hz
   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1
            ndis=abs(i-i0)+abs(j-j0)+abs(k-k0)
            if (ndis <= 1) then
               if (phi(i,j,k) <= 0.0) then
                  cd(nc)=-3.0*b_in/hxyz!fb_in(b_in,b_out,x(i),y(j),z(k))/hxyz
               else
                  cd(nc)=-3.0*b_out/hxyz!fb_out(b_in,b_out,x(i),y(j),z(k))/hxyz
               end if
            else
               cd(nc)=0.0
            end if

            if (ndis == 0) cd(nc)=-6.0*cd(nc)
            
         end do
      end do
   end do

   if (info > 3) then
      call cpymat(immax,inmax,ca,caa)
   end if

   iprint = 0
   iwar(1) = 1

   ! ----- solve quadratic programming
         
   call ql0001(im,ime,immax,in,inmax,imnn,cc,cd,ca,cb,xl,xu, &
      xx,uu,iout,ifail,iprint,war,lwar,iwar,liwar,eps)
 
   tmp=0.0
   do ii=1,10
      tmp1=0.0
      do ji=1,inmax
         tmp1=tmp1+ca(ii,ji)*xx(ji)
      end do
      if(abs(tmp1+cb(ii)) > tmp) tmp = abs(tmp1+cb(ii))
   end do

   if(ifail > 10) then
      write(*,*) "Warning!", i0,j0,k0, "Irre31 Inconsistent Constraints"
      return
      do i=1, in
         xx(i) = 0.0
      end do
      call regula(l,m,n,hx,hy,hz,b_in,b_out,i0,j0,k0,info,x,y,z, tmpcoe,rhs)
      xx(5)  = tmpcoe(2)
      xx(23) = tmpcoe(3)
      xx(11) = tmpcoe(4)
      xx(17) = tmpcoe(5)
      xx(13) = tmpcoe(6)
      xx(15) = tmpcoe(7)
      xx(14) = tmpcoe(1)
      return
   end if

   do i=1,20
      a(i) = 0.0
   end do

   !About the local coordinates:note that in the test point i dir and ksi are opposite

   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1

            nc = nc + 1
            tmpx(1) = x(i) - x1
            tmpx(2) = y(j) - y1
            tmpx(3) = z(k) - z1
            call matvec(3,3,t, tmpx, xloc)

            if (index(i,j,k) <= 3) then

               a(1) =a(1) +xx(nc)
               a(3) =a(3) +xx(nc)*xloc(1)
               a(5) =a(5) +xx(nc)*xloc(2)
               a(7) =a(7) +xx(nc)*xloc(3)
               a(9) =a(9) +xx(nc)*xloc(1)*xloc(1)*0.5
               a(11)=a(11)+xx(nc)*xloc(2)*xloc(2)*0.5
               a(13)=a(13)+xx(nc)*xloc(3)*xloc(3)*0.5
               a(15)=a(15)+xx(nc)*xloc(1)*xloc(2)
               a(17)=a(17)+xx(nc)*xloc(1)*xloc(3)
               a(19)=a(19)+xx(nc)*xloc(2)*xloc(3)
            else
               a(2) =a(2) +xx(nc)
               a(4) =a(4) +xx(nc)*xloc(1)
               a(6) =a(6) +xx(nc)*xloc(2)
               a(8) =a(8) +xx(nc)*xloc(3)
               a(10)=a(10)+xx(nc)*xloc(1)*xloc(1)*0.5
               a(12)=a(12)+xx(nc)*xloc(2)*xloc(2)*0.5
               a(14)=a(14)+xx(nc)*xloc(3)*xloc(3)*0.5
               a(16)=a(16)+xx(nc)*xloc(1)*xloc(2)
               a(18)=a(18)+xx(nc)*xloc(1)*xloc(3)
               a(20)=a(20)+xx(nc)*xloc(2)*xloc(3)
            end if
         end do
      end do  !  j=j0-1,j0+1
   end do  !  i=i0-1,i0+1
   tmp  = 1.0/b_out
   tmp1 =   a(4) + a(10)*(xyy+xzz-bl_out(1)*tmp) &
         - a(12)*xyy - a(14)*xzz - a(16)*bl_out(2)*tmp &
         - a(18)*bl_out(3)*tmp - a(20)*xyz
   tmp2 =   a(6) - a(10)*bl_out(2)*tmp &
         + a(16)*xyy + a(18)*xyz
   tmp3 =   a(8) - a(10)*bl_out(3)*tmp &
         + a(16)*xyz + a(18)*xzz
   corr =   a(10)*((fjmp-fkk_out*ujmp)*tmp-wyy-wzz) &
         + a(12)*wyy + a(14)*wzz + a(16)*qy*tmp &
         + a(18)*qz*tmp + a(20)*wyz + a(2)*ujmp &
         + tmp*tmp1*q + tmp2*wy + tmp3*wz

   if (info > 3) then
      corr = corr + fkk_out*ujmp - fjmp
   end if

   if (info <= 3) then
      fkk = ZERO!fk_in(x(i0),y(j0),z(k0))
      rhs = ZERO!ff_in(x(i0),y(j0),z(k0))
   else
      fkk = ZERO!fk_out(x(i0),y(j0),z(k0))
      rhs = ZERO!ff_out(x(i0),y(j0),z(k0))
   end if

   xx(kctr) = xx(kctr)+fkk
   rhs = rhs + corr
  
end subroutine irre31 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine irre31 here]
subroutine irre32(l,m,n,h,hx,hy,hz,ifail, &
              i0,j0,k0,info,bmax,b_in,b_out,x,y,z,phi,index, &
              nq,nirreg,index2,cirreg,xx,rhs, &
              wp,wyp,wzp,wyyp,wzzp,wyzp,qp,qyp,qzp)

   ! for an irregular point (i0,j0,k0), find corresponding coefficients

   implicit none

   integer,parameter :: im=4, ime=4, immax=im, ns=27, in=ns, &
         inmax=in, imnn=im+2*in, &
         lwar=3*inmax*inmax/2+10*inmax+2*immax+1, &
         liwar=in

   ! Passed variables

   integer l,m,n,ifail,i0,j0,k0,info,nq,nirreg
   integer index(l,m,n),index2(l,m,n)
   _REAL_  h,hx,hy,hz,bmax,b_in,b_out,rhs
   _REAL_  x(0:l+1), y(0:m+1), z(0:n+1)
   _REAL_  phi(0:l+1,0:m+1, 0:n+1)
   _REAL_  cirreg(15, nirreg)
   _REAL_  wp(nirreg), &
           wyp(nirreg), wzp(nirreg), wyyp(nirreg), wzzp(nirreg), wyzp(nirreg), &
           qp(nirreg), qyp(nirreg), qzp(nirreg)


   ! Local variables

   integer  nc,kr,kctr,ndis,ji,i,j,k,ir,iprint,ii,iout,ijk
   _REAL_   xyy,xzz,xyz,x1,y1,z1,ujmp,wy,wz,wyy,wzz,wyz,q,qy,qz,&
            tmp,temp,tmp1,tmp2,tmp3,rho,&
            rho1,hxx,hyy,hzz,hxyz,fmin,h1,bl2jmp,bl3jmp,&
            fjmp,corr,eps,fkk_in,fkk_out,fkk,fkjmp

   _REAL_   t(3,3),a(20)
   _REAL_   xloc(3),tmpx(3)
   _REAL_   bl_out(3),bl_in(3)
   _REAL_   cc(inmax,inmax), ca(immax,inmax), caa(immax,inmax)
   _REAL_   cd(inmax),cb(immax),xx(in),xl(in),xu(in),uu(imnn)
   _REAL_   war(lwar),iwar(liwar)
   _REAL_   tmpcoe(7)

   ir = index2(i0, j0, k0)

   x1     = cirreg( 1, ir)
   y1     = cirreg( 2, ir)
   z1     = cirreg( 3, ir)

   xyy    = cirreg( 4, ir)
   xzz    = cirreg( 5, ir)
   xyz    = cirreg( 6, ir)

   t(1,1) = cirreg( 7, ir)
   t(2,1) = cirreg( 8, ir)
   t(3,1) = cirreg( 9, ir)
   t(1,2) = cirreg(10, ir)
   t(2,2) = cirreg(11, ir)
   t(3,2) = cirreg(12, ir)
   t(1,3) = cirreg(13, ir)
   t(2,3) = cirreg(14, ir)
   t(3,3) = cirreg(15, ir)

   ! ----- find jump conditions

!   call fjmps(x1,y1,z1, fjmp)
   fjmp=ZERO
!   call fkjmps(x1,y1,z1, fkjmp)
   fkjmp=ZERO

   ujmp = wp(ir)
   wy   = wyp(ir)
   wz   = wzp(ir)
   wyy  = wyyp(ir)
   wzz  = wzzp(ir)
   wyz  = wyzp(ir)

   q =  qp(ir)
   qy = qyp(ir)
   qz = qzp(ir)


!   call betas(hx,hy,hz,x1,y1,z1,t,bl_in,bl_out,b_in,b_out)
   bl_in=ZERO
   bl_out=ZERO

   fkk_in = 0.0d0!fk_in(x1,y1,z1), the kappa term inside
   fkk_out = 0.0d0!fk_out(x1,y1,z1), the kappa term outside

   hxx=hx*hx
   hyy=hy*hy
   hzz=hz*hz
   h1=fmin(hxx,fmin(hyy,hzz))

   eps=0.1e-11

   do i=1,immax
      do j=1,inmax
         ca(i,j) = 0.0
      end do
      cb(i) = 0.0
   end do

   nc = 0
   kr = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1
            kr = kr +1
            if (i == i0 .and. j == j0 .and. k == k0) then
               !               kr = kr -1
               kctr = nc
               !               ca(kr+ime,nc) = -1.0d0
               !               cb(kr+ime) = -1.0d-20
            else
               !               ca(kr+ime,nc) = 1.0d0
               !               cb(kr+ime) = -100.0d0
            end if
            xl(nc) = 0.0
            xu(nc) = 2.0*bmax/h1
         end do
      end do
   end do

   xl(kctr)=-12.0*bmax/h1
   xu(kctr)=0.0

   rho = b_in/b_out
   rho1 = 1 - rho
   bl2jmp = bl_out(2)-bl_in(2)
   bl3jmp = bl_out(3)-bl_in(3)
   temp = fkjmp/b_out

   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1

            ! ----------- local coordinates transformation

            tmpx(1) = x(i) - x1
            tmpx(2) = y(j) - y1
            tmpx(3) = z(k) - z1
            call matvec(3,3,t, tmpx, xloc)

            ! ----------- form the coefficients of equality constraints

            if (index(i,j,k) <= 3) then
               ca(1,nc)  = 1.0
               ca(2,nc)  = xloc(1)
               ca(3,nc)  = xloc(2)
               ca(4,nc)  = xloc(3)
            else
               ca(1,nc)  = 1.0-0.5*temp*xloc(1)*xloc(1)
               ca(2,nc)  = rho*xloc(1) &
                     - 0.5*xloc(1)*xloc(1)*(rho1*(xyy+xzz) &
                     - bl_in(1)/b_out + rho*bl_out(1)/b_out) &
                     + 0.5*xloc(2)*xloc(2)*xyy*rho1 &
                     + 0.5*xloc(3)*xloc(3)*xzz*rho1 &
                     + xloc(1)*xloc(2)*(bl_in(2)/b_out &
                     - rho*bl_out(2)/b_out) &
                     + xloc(1)*xloc(3)*(bl_in(3)/b_out &
                     - rho*bl_out(3)/b_out) &
                     + xloc(2)*xloc(3)*xyz*rho1
               ca(3,nc)  = - 0.5*xloc(1)*xloc(1)*bl2jmp/b_out &
                     + xloc(1)*xloc(2)*xyy*rho1 &
                     + xloc(1)*xloc(3)*xyz*rho1 &
                     + xloc(2)
               ca(4,nc)  = - 0.5*xloc(1)*xloc(1)*bl3jmp/b_out &
                     + xloc(1)*xloc(3)*xzz*rho1 &
                     + xloc(1)*xloc(2)*xyz*rho1 &
                     + xloc(3)
            end if
         end do  !  k=k0-1,k0+1
      end do  !  j=j0-1,j0+1
   end do  !  i=i0-1,i0+1

   if (info > 3) then
      cb(1) = fkjmp
   end if

   cb(2)  = -bl_in(1)
   cb(3)  = -bl_in(2)
   cb(4)  = -bl_in(3)

   do i=1,inmax
      do j=1,inmax
         cc(i,j) = 0.0
      end do
      cc(i,i) = 1.0
   end do

   hxyz=hx*hx+hy*hy+hz*hz
   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1
            nc = nc + 1
            ndis=abs(i-i0)+abs(j-j0)+abs(k-k0)
            if (ndis <= 1) then
               if (phi(i,j,k) <= 0.0) then
                  cd(nc)=-3.0*b_in!fb_in(b_in,b_out,x(i),y(j),z(k))/hxyz
               else
                  cd(nc)=-3.0*b_out!fb_out(b_in,b_out,x(i),y(j),z(k))/hxyz
               end if
            else
               cd(nc)=0.0
            end if

            if (ndis == 0) cd(nc)=-6.0*cd(nc)
         end do
      end do
   end do

   if (info > 3) then
      call cpymat(immax,inmax,ca,caa)
   end if

   iprint = 1
   iwar(1) = 1

   ! ----- solve quadratic programming

   call ql0001(im,ime,immax,in,inmax,imnn,cc,cd,ca,cb,xl,xu, &
         xx,uu,iout,ifail,iprint,war,lwar,iwar,liwar,eps)

   tmp=0.0
   do ii=1,4
      tmp1=0.0
      do ji=1,inmax
         tmp1=tmp1+ca(ii,ji)*xx(ji)
      end do
      if(abs(tmp1+cb(ii)) > tmp) tmp = abs(tmp1+cb(ii))
   end do

   if(ifail > 10) then
      write(*,*) "irre32-Warning!",i0,j0,k0,"Inconsistent Constraints"
      do i=1, in
         xx(i) = 0.0
      end do
      call regula(l,m,n,hx,hy,hz,b_in,b_out,i0,j0,k0,info,x,y,z, tmpcoe,rhs)
      xx(5)  = tmpcoe(2)
      xx(23) = tmpcoe(3)
      xx(11) = tmpcoe(4)
      xx(17) = tmpcoe(5)
      xx(13) = tmpcoe(6)
      xx(15) = tmpcoe(7)
      xx(14) = tmpcoe(1)
      return
   end if

   do i=1,20
      a(i) = 0.0
   end do

   nc = 0
   do i=i0-1,i0+1
      do j=j0-1,j0+1
         do k=k0-1,k0+1

            nc = nc + 1
            tmpx(1) = x(i) - x1
            tmpx(2) = y(j) - y1
            tmpx(3) = z(k) - z1
            call matvec(3,3,t, tmpx, xloc)

            if (index(i,j,k) <= 3) then
               a(1) =a(1) +xx(nc)
               a(3) =a(3) +xx(nc)*xloc(1)
               a(5) =a(5) +xx(nc)*xloc(2)
               a(7) =a(7) +xx(nc)*xloc(3)
               a(9) =a(9) +xx(nc)*xloc(1)*xloc(1)*0.5
               a(11)=a(11)+xx(nc)*xloc(2)*xloc(2)*0.5
               a(13)=a(13)+xx(nc)*xloc(3)*xloc(3)*0.5
               a(15)=a(15)+xx(nc)*xloc(1)*xloc(2)
               a(17)=a(17)+xx(nc)*xloc(1)*xloc(3)
               a(19)=a(19)+xx(nc)*xloc(2)*xloc(3)
            else
               a(2) =a(2) +xx(nc)
               a(4) =a(4) +xx(nc)*xloc(1)
               a(6) =a(6) +xx(nc)*xloc(2)
               a(8) =a(8) +xx(nc)*xloc(3)
               a(10)=a(10)+xx(nc)*xloc(1)*xloc(1)*0.5
               a(12)=a(12)+xx(nc)*xloc(2)*xloc(2)*0.5
               a(14)=a(14)+xx(nc)*xloc(3)*xloc(3)*0.5
               a(16)=a(16)+xx(nc)*xloc(1)*xloc(2)
               a(18)=a(18)+xx(nc)*xloc(1)*xloc(3)
               a(20)=a(20)+xx(nc)*xloc(2)*xloc(3)
            end if
         end do
      end do  !  j=j0-1,j0+1
   end do  !  i=i0-1,i0+1

   tmp  = 1.0/b_out
   tmp1 =   a(4) + a(10)*(xyy+xzz-bl_out(1)*tmp) &
         - a(12)*xyy - a(14)*xzz - a(16)*bl_out(2)*tmp &
         - a(18)*bl_out(3)*tmp - a(20)*xyz
   tmp2 =   a(6) - a(10)*bl_out(2)*tmp &
         + a(16)*xyy + a(18)*xyz
   tmp3 =   a(8) - a(10)*bl_out(3)*tmp &
         + a(16)*xyz + a(18)*xzz

   corr =   a(10)*((fjmp-fkk_out*ujmp)*tmp-wyy-wzz) &
         + a(12)*wyy + a(14)*wzz + a(16)*qy*tmp &
         + a(18)*qz*tmp + a(20)*wyz + a(2)*ujmp &
         + tmp*tmp1*q + tmp2*wy + tmp3*wz

   !.......if the grid point is on (+) side

   if (info > 3) then
      corr = corr + fkk_out*ujmp - fjmp
   end if

   if (info <= 3) then
      fkk = 0.0d0!fk_in(x(i0),y(j0),z(k0))
      rhs = 0.0d0!ff_in(x(i0),y(j0),z(k0))
   else
      fkk = 0.0d0!fk_out(x(i0),y(j0),z(k0))
      rhs = 0.0d0!ff_out(x(i0),y(j0),z(k0))
   end if

   xx(kctr) = xx(kctr)+fkk
   rhs = rhs + corr

end subroutine irre32
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ local copy of matvec
subroutine matvec(m,n,a,x,y)

   implicit none
   integer m,n,i,j
   _REAL_  a(m,n), x(n),y(m)

   do i=1,m
      y(i) = 0.0
      do j=1,n
         y(i) = y(i) + a(i,j)*x(j)
      end do
   end do

end subroutine matvec
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ driver to set up the jump conditions and their derivatives
subroutine jumps(l,m,n,nirreg,bcopt,bi,bo,atmfirst,atmlast,cirreg,&
                  wp,wyp,wzp,wyyp,wzzp,wyzp,qp,qyp,qzp)

   implicit none

   ! passed variables

   integer l, m, n, nirreg
   integer bcopt, atmfirst, atmlast
   _REAL_ cirreg(15, nirreg)
   _REAL_ t(3,3), c(3)
   _REAL_ bi,bo
   _REAL_ wp(nirreg), &
          wyp(nirreg), wzp(nirreg), wyyp(nirreg), wzzp(nirreg), wyzp(nirreg), &
          qp(nirreg), qyp(nirreg), qzp(nirreg)

   ! local variables

   integer ir
   _REAL_ xx, yy, zz
   _REAL_, allocatable :: lacrd(:,:)

   wp  = ZERO
   wyp = ZERO
   wzp = ZERO
   wyyp= ZERO
   wzzp= ZERO
   wyzp= ZERO
   qp  = ZERO
   qyp = ZERO
   qzp = ZERO

   allocate(lacrd(4,atmlast-atmfirst+1))

   do ir = 1, nirreg
      xx     = cirreg( 1, ir)
      yy     = cirreg( 2, ir)
      zz     = cirreg( 3, ir)
      c(1)   = cirreg( 4, ir)
      c(2)   = cirreg( 5, ir)
      c(3)   = cirreg( 6, ir)
      t(1,1) = cirreg( 7, ir)
      t(1,2) = cirreg( 8, ir)
      t(1,3) = cirreg( 9, ir)
      t(2,1) = cirreg(10, ir)
      t(2,2) = cirreg(11, ir)
      t(2,3) = cirreg(12, ir)
      t(3,1) = cirreg(13, ir)
      t(3,2) = cirreg(14, ir)
      t(3,3) = cirreg(15, ir)

      ! cirreg(4:6,ir) is the curvature for point ir
      ! the first two are the values of w at projection point, the next two are the
      ! first derivatives along the two tangential directions, the next three are
      ! the second derivatives wyy, wzz, wyz; for q the first two are the values of q at
      ! projection point, the next two are the first derivatives of q along the two
      ! tangential directions

      call problem(bcopt,ir,atmfirst,atmlast,bi,bo,xx,yy,zz,c,t)
   end do

   deallocate(lacrd)

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The IIM method is capable of solving a generic PB equation.
! bcopt = 4 or 5: inside is total field and outside is total field.
! bcopt = 2 or 6 or 7: inside is reaction field and outside is total field.
! bcopt = 8 or 9: inside is reaction field and outside is reaction field.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Jump conditions of potential (phi) at projection point (x,y,z)
! bcopt = 4 or 5: it is zero.
! bcopt = 2 or 6 or 7: it is coulomb potential.
! bcopt = 8 or 9: it is zero.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Jump conditions of field (b*E_n) at projection point (x,y,z)
! bcopt = 4 or 5: it is zero. 
! bcopt = 2 or 6 or 7: it is  bi  * coulomb field.
! bcopt = 8 or 9: it is (bi - bo) * coulomb field.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ The jump conditions and their derivatives at projection point (x,y,z)
subroutine problem(bcopt,l,atmfirst,atmlast,bi,bo,x,y,z,curv,matrix)

   use poisson_boltzmann, only : acrg, acrd
   implicit none

   integer bcopt,l,atmfirst,atmlast
   _REAL_ :: bi,bo
   _REAL_ x,y,z
   _REAL_, dimension(3) :: curv
   _REAL_, dimension(3,3) :: matrix

   integer k
   _REAL_ :: dinv2, dinv3, dinv5
   _REAL_ :: lacrg
   _REAL_ :: lacrg_time_dinv5, lacrg_time_dinv3

   ! the first three are atom position in local coordinate, the last one is the
   ! reverse distance between the atom and the local origin

   if ( bcopt == 4 .or. bcopt == 5 ) return

   ! compute needed displacement vector in the local frame
   ! will test the use of cutoff later

   do k = atmfirst, atmlast
      lacrd(1,k)=(acrd(1,k)-x)*matrix(1,1)+&
                 (acrd(2,k)-y)*matrix(2,1)+&
                 (acrd(3,k)-z)*matrix(3,1)
      lacrd(2,k)=(acrd(1,k)-x)*matrix(1,2)+&
                 (acrd(2,k)-y)*matrix(2,2)+&
                 (acrd(3,k)-z)*matrix(3,2)
      lacrd(3,k)=(acrd(1,k)-x)*matrix(1,3)+&
                 (acrd(2,k)-y)*matrix(2,3)+&
                 (acrd(3,k)-z)*matrix(3,3)
      lacrd(4,k)=ONE/sqrt(lacrd(1,k)**2+lacrd(2,k)**2+lacrd(3,k)**2)
   end do

   ! compute the jump conditionsi and their derivatives using Coulomb's law

   do k = atmfirst, atmlast
      dinv2=lacrd(4,k)**2
      dinv3=dinv2*lacrd(4,k)
      dinv5=dinv3*dinv2
      lacrg=acrg(k)
      lacrg_time_dinv3=lacrg*dinv3
      lacrg_time_dinv5=lacrg*dinv5

      wp (l) = wp(l) + lacrg*lacrd(4,k)

      wyp(l) = wyp(l) + lacrg_time_dinv3*lacrd(2,k)
      wzp(l) = wzp(l) + lacrg_time_dinv3*lacrd(3,k)

      wyyp(l)= wyyp(l) + THREE*lacrg_time_dinv5*lacrd(2,k)**2-&
                               lacrg_time_dinv3+&
                               lacrg_time_dinv3*lacrd(1,k)*curv(1)
      wzzp(l)= wzzp(l) + THREE*lacrg_time_dinv5*lacrd(3,k)**2-&
                               lacrg_time_dinv3+&
                               lacrg_time_dinv3*lacrd(1,k)*curv(2)
      wyzp(l)= wyzp(l) + THREE*lacrg_time_dinv5*lacrd(2,k)*lacrd(3,k)+&
                               lacrg_time_dinv3*lacrd(1,k)*curv(3)

      qp (l) = qp(l)   +       lacrg_time_dinv3*lacrd(1,k)

      qyp(l) = qyp(l)  + THREE*lacrg_time_dinv5*lacrd(1,k)*lacrd(2,k)
      qzp(l) = qzp(l)  + THREE*lacrg_time_dinv5*lacrd(1,k)*lacrd(3,k)
   end do

   if ( bcopt ==2 .or. bcopt == 6 .or. bcopt == 7 ) then
      qp (l) =  bi      *qp (l)
      qyp(l) =  bi      *qyp(l)
      qzp(l) =  bi      *qzp(l)
   end if
   if ( bcopt == 8 .or. bcopt == 9 ) then
      qp (l) = (bi - bo)*qp (l)
      qyp(l) = (bi - bo)*qyp(l)
      qzp(l) = (bi - bo)*qzp(l)

      wp (l) = ZERO
      wyp(l) = ZERO
      wzp(l) = ZERO
      wyyp(l) = ZERO
      wzzp(l) = ZERO
      wyzp(l) = ZERO
   end if

end subroutine problem

end subroutine jumps
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine indexg here]
subroutine ana_irreg_grd(natm,gcrd,radih,dprobh,xm,ym,zm,nbnd,iepsav,goxh,goyh,gozh,h)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! set up irregular grid point index and data structure for later use:
! 1) projection point on the interface
! 2) normal vector
! 3) transformation matrix to local frame
! 4) curvature
! analytical version
! 
! input:
! gcrd(3,natm) atomic coordinates in the FD grid frame.
! radih(natm)  atomic radii in the FD grid unit.
! dprobh       solvent probe in the FD grid unit.
!
! xm,ym,zm     FD grid dimensions
!
! nbnd         number of irregular grid points
! bndatmptr(nbnd+1)
!              for finding to which atom the irregular point belongs
! bndatmlst(maxlist)
!              atom neighbor list
!
! iepsav(4,nbnd)
!              i,j,k of irregular grid point, and the atom owner
!
! output:
! cirreg(15,nbnd)
!              coords, curvatures, transformation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use density_surface, only: cirreg, bndatmlst, bndatmptr
   implicit none

   integer :: natm, nbnd
   _REAL_ :: dprobh
   integer :: xm, ym, zm
   _REAL_, dimension(3,natm) :: gcrd
   _REAL_, dimension(natm) :: radih
   integer, dimension(4,nbnd) :: iepsav
   _REAL_ :: goxh, goyh, gozh, h

   _REAL_ :: factor
   _REAL_, parameter :: small = 1.0d-6

   integer :: i, j, k, l

   _REAL_ :: x, y, z ! coordinates in the grid frame
   _REAL_ :: nextlvlset, middlelvlset
   _REAL_ :: x_proj, y_proj, z_proj, x_proj_end, y_proj_end, z_proj_end
   _REAL_, dimension(3) :: proj_point
   _REAL_, dimension(3) :: deriv, deriv_square
   _REAL_, dimension(3) :: gradient
   _REAL_ :: deriv_inv, deriv2_inv, deriv3_inv, deriv2
   _REAL_, dimension(3,3) :: trans_matrix
   _REAL_, dimension(3) :: bnd_deriv

   ! conversion factor from grid-based distance to 2*dprobe based distance

   factor = HALF/dprobh

   ! a few words on the notation
   ! x means the normal direction, y and z arethe  tangential directions
   ! the first one is the second-order tangential derivative along the first
   ! direction(use yy as the indicator), the second one is that along the second
   ! direction(use zz as the indicator), the last one is the mixed derivative
   ! of both of the direction(use yz as the indicator)

   ! For printing the interface
   !character(10) str
   !open(unit=58,file='interface.dot')
   !write (58, *) nbnd
   !write (58, '("DOTS")')

   do l = 1, nbnd
      i = iepsav(1, l); j = iepsav(2, l); k = iepsav(3, l)
      x = dble(i); y = dble(j); z = dble(k)

      ! finding gradient vector at irregular grid point
      ! note the gradient here can't be zero

      gradient(1:3) = find_gradient(x,y,z,l)
      if ( max(abs(gradient(1)),abs(gradient(2)),abs(gradient(3))) < small ) then
         print *, 'Warning: There is no gradient for irregular point at', x,y,z
         call mexit(6,1)
      end if

      ! finding root/proj point using Newton's Method 
      
      call find_projection(x,y,z,l,proj_point)
      x_proj=proj_point(1)
      y_proj=proj_point(2)
      z_proj=proj_point(3) 
  
      ! compute the 1st derivative at this point
      ! the minus sign here is to be consistent with the convention in the numerical code

      deriv(1:3) = -find_gradient(x_proj,y_proj,z_proj,l)

      ! set up global to local transformation matrix

      call transformation(deriv,trans_matrix)

      ! comptue curvature at the projection point

      call curvature(x_proj,y_proj,z_proj,l,trans_matrix,bnd_deriv)

      ! time to save to the global array
      ! convert projection coord to the angstrom unit

      cirreg(1, l) = (x_proj+goxh)*h
      cirreg(2, l) = (y_proj+goyh)*h
      cirreg(3, l) = (z_proj+gozh)*h
      !write(str,'(i10)') l
      !str = adjustl(str)
      !write(58,'(a,3(f20.15),2x)') "H"//str, cirreg(1:3,l)

      ! convert the curvature to the angstrom unit
      ! most curvatures are positive because gradients point inward towards the atoms

      cirreg(4:6, l) = bnd_deriv(1:3)/h

      ! here is the matrix

      cirreg(7, l) = trans_matrix(1, 1)
      cirreg(8, l) = trans_matrix(1, 2)
      cirreg(9, l) = trans_matrix(1, 3)
      cirreg(10, l) = trans_matrix(2, 1)
      cirreg(11, l) = trans_matrix(2, 2)
      cirreg(12, l) = trans_matrix(2, 3)
      cirreg(13, l) = trans_matrix(3, 1)
      cirreg(14, l) = trans_matrix(3, 2)
      cirreg(15, l) = trans_matrix(3, 3)

   end do
   !close(58)

contains

subroutine find_projection(x,y,z,l,proj_point)

   _REAL_ :: x, y, z
   integer :: l
   _REAL_, dimension(3) :: proj_point

   x_proj = x
   y_proj = y
   z_proj = z
   nextlvlset = find_lvlset(x_proj,y_proj,z_proj,l)
   if ( nextlvlset*&
        find_lvlset(x-sign(ONE,nextlvlset)*gradient(1),&
                    y-sign(ONE,nextlvlset)*gradient(2),&
                    z-sign(ONE,nextlvlset)*gradient(3),l)&
        <=0 ) then
      do while ( abs(nextlvlset) > small )
         x_proj = x_proj - nextlvlset*gradient(1)
         y_proj = y_proj - nextlvlset*gradient(2)
         z_proj = z_proj - nextlvlset*gradient(3)
         middlelvlset = nextlvlset
         nextlvlset = find_lvlset(x_proj,y_proj,z_proj,l)
         if ( middlelvlset*nextlvlset<=0 ) then
            x_proj = x_proj + HALF*middlelvlset*gradient(1)
            y_proj = y_proj + HALF*middlelvlset*gradient(2)
            z_proj = z_proj + HALF*middlelvlset*gradient(3)
            nextlvlset = find_lvlset(x_proj,y_proj,z_proj,l)
         end if
      end do
   else ! bisection method
      if      ( nextlvlset*find_lvlset(x+ONE,y,z,l) <= ZERO ) then
         x_proj_end = x+ONE
         y_proj_end = y
         z_proj_end = z
      else if ( nextlvlset*find_lvlset(x-ONE,y,z,l) <= ZERO ) then
         x_proj_end = x-ONE
         y_proj_end = y
         z_proj_end = z
      else if ( nextlvlset*find_lvlset(x,y+ONE,z,l) <= ZERO ) then
         x_proj_end = x
         y_proj_end = y+ONE
         z_proj_end = z
      else if ( nextlvlset*find_lvlset(x,y-ONE,z,l) <= ZERO ) then
         x_proj_end = x
         y_proj_end = y-ONE
         z_proj_end = z
      else if ( nextlvlset*find_lvlset(x,y,z+ONE,l) <= ZERO ) then
         x_proj_end = x
         y_proj_end = y
         z_proj_end = z+ONE
      else
         x_proj_end = x
         y_proj_end = y
         z_proj_end = z-ONE
         if ( nextlvlset*find_lvlset(x,y,z-ONE,l) > ZERO ) then
            write(6,*) ' PBSA Bomb in ana_irreg_grid/find_projection(): '
            write(6,*) ' Failed to find projection point ',&
            l,i,j,k,x,y,z,nextlvlset,&
            find_lvlset(x+1,y,z,l),find_lvlset(x-1,y,z,l),&
            find_lvlset(x,y+1,z,l),find_lvlset(x,y-1,z,l),&
            find_lvlset(x,y,z+1,l),find_lvlset(x,y,z-1,l)
            call mexit(6,1)
         end if
      end if
      middlelvlset = find_lvlset(HALF*(x_proj+x_proj_end),HALF*(y_proj+y_proj_end),HALF*(z_proj+z_proj_end),l)
      do while ( abs(middlelvlset) > small )
         if ( nextlvlset*middlelvlset <= 0 ) then
            x_proj_end = HALF*(x_proj+x_proj_end)
            y_proj_end = HALF*(y_proj+y_proj_end)
            z_proj_end = HALF*(z_proj+z_proj_end)
         else
            x_proj     = HALF*(x_proj+x_proj_end)
            y_proj     = HALF*(y_proj+y_proj_end)
            z_proj     = HALF*(z_proj+z_proj_end)
            nextlvlset = middlelvlset
         end if
         middlelvlset = find_lvlset(HALF*(x_proj+x_proj_end),HALF*(y_proj+y_proj_end),HALF*(z_proj+z_proj_end),l)
      end do
      x_proj = HALF*(x_proj+x_proj_end)
      y_proj = HALF*(y_proj+y_proj_end)
      z_proj = HALF*(z_proj+z_proj_end)
   end if
   proj_point(1)=x_proj
   proj_point(2)=y_proj
   proj_point(3)=z_proj

end subroutine find_projection

subroutine transformation(deriv,trans_matrix)

   _REAL_, dimension(3) :: deriv
   _REAL_, dimension(3,3) :: trans_matrix

   deriv_square(1:3) = deriv(1:3)**2
   deriv_inv = ONE/sqrt(sum(deriv_square(1:3)))
   if ( deriv_square(3) < deriv_square(2) ) then
      deriv2 = sqrt(sum(deriv_square(1:2)))
      deriv2_inv = ONE/deriv2
      deriv3_inv = ONE/sqrt(&
           deriv_square(1)*deriv_square(3) +&
           deriv_square(2)*deriv_square(3) +&
          (deriv_square(1)+deriv_square(2))**2)

      trans_matrix(1,1) =  deriv(1)*deriv_inv
      trans_matrix(2,1) =  deriv(2)*deriv_inv
      trans_matrix(3,1) =  deriv(3)*deriv_inv
      trans_matrix(1,2) =  deriv(2)*deriv2_inv
      trans_matrix(2,2) = -deriv(1)*deriv2_inv
      trans_matrix(3,2) =  ZERO
      trans_matrix(1,3) =  deriv(1)*deriv(3)*deriv3_inv
      trans_matrix(2,3) =  deriv(2)*deriv(3)*deriv3_inv
      trans_matrix(3,3) =(-deriv_square(1)-deriv_square(2))*deriv3_inv
   else
      deriv2=sqrt(deriv_square(1)+deriv_square(3))
      deriv2_inv=ONE/deriv2
      deriv3_inv=ONE/sqrt(&
           deriv_square(1)*deriv_square(2) +&
           deriv_square(2)*deriv_square(3) +&
          (deriv_square(1)+deriv_square(3))**2)

      trans_matrix(1,1) =  deriv(1)*deriv_inv
      trans_matrix(2,1) =  deriv(2)*deriv_inv
      trans_matrix(3,1) =  deriv(3)*deriv_inv
      trans_matrix(1,2) =  deriv(3)*deriv2_inv
      trans_matrix(2,2) =  ZERO
      trans_matrix(3,2) = -deriv(1)*deriv2_inv
      trans_matrix(1,3) = -deriv(1)*deriv(2)*deriv3_inv
      trans_matrix(2,3) = (deriv_square(1)+deriv_square(3))*deriv3_inv
      trans_matrix(3,3) = -deriv(3)*deriv(2)*deriv3_inv
   end if
 
end subroutine transformation

subroutine curvature(x,y,z,l,matrix,bnd_deriv)
   implicit none

   _REAL_ :: x, y, z
   _REAL_, dimension(3,3) :: matrix
   _REAL_, dimension(3) :: bnd_deriv
   integer :: l

   _REAL_ :: factor
   _REAL_ :: dx, dy, dz, ldx, ldy, ldz
   _REAL_ :: disth, dist
   _REAL_, dimension(3) :: lvlset_deriv_norm
   integer :: atmfirst, atmlast
   integer :: jp, jatm

   factor = HALF/dprobh

   bnd_deriv(1:3) = ZERO
   lvlset_deriv_norm(1:3) = ZERO
   atmfirst=bndatmptr(l-1)+1
   atmlast=bndatmptr(l)
   do jp=atmfirst,atmlast
      jatm=bndatmlst(jp)
      dx = gcrd(1,jatm) - x
      dy = gcrd(2,jatm) - y
      dz = gcrd(3,jatm) - z
      ! ldx represents the normal direction, the rest
      ! two are tangential directions
      ldx = dx*matrix(1,1) + dy*matrix(2,1) + dz*matrix(3,1)
      ldy = dx*matrix(1,2) + dy*matrix(2,2) + dz*matrix(3,2)
      ldz = dx*matrix(1,3) + dy*matrix(2,3) + dz*matrix(3,3)
      disth = sqrt(ldx**2+ldy**2+ldz**2)
      dist = (disth-radih(jatm))*factor
      lvlset_deriv_norm(1:3) = lvlset_deriv_norm(1:3)+&
         density_1deriv_atom(dist,-ldx/disth*factor,ZERO,ZERO)
      bnd_deriv(1:3) = bnd_deriv(1:3)+&
         density_2deriv_atom(dist,disth,-ldy,-ldz)
   end do
   bnd_deriv(1) = -bnd_deriv(1)/lvlset_deriv_norm(1)
   bnd_deriv(2) = -bnd_deriv(2)/lvlset_deriv_norm(1)
   bnd_deriv(3) = -bnd_deriv(3)/lvlset_deriv_norm(1)

end subroutine curvature

function find_lvlset(x,y,z,l)
   implicit none
   
   _REAL_ :: find_lvlset
   _REAL_ :: x, y, z
   integer :: l

   integer :: jp, jatm
   integer :: atmfirst, atmlast 
   _REAL_ :: dx, dy, dz, disth
   
   find_lvlset = -ONE
   atmfirst = bndatmptr(l-1)+1
   atmlast = bndatmptr(l)
   do jp = atmfirst, atmlast
      jatm = bndatmlst(jp)
      dx = x-gcrd(1,jatm)
      dy = y-gcrd(2,jatm)
      dz = z-gcrd(3,jatm)
      disth = sqrt(dx**2+dy**2+dz**2)
      find_lvlset = find_lvlset + density_atom((disth-radih(jatm))*factor)
   end do

end function find_lvlset

function find_gradient(x,y,z,l)
   implicit none
   
   _REAL_, dimension(3) :: find_gradient
   _REAL_ :: x,y,z
   integer :: l

   integer :: jp, jatm
   _REAL_ :: dx, dy, dz
   _REAL_ :: disth, dist_inv, dist
   integer :: atmfirst, atmlast
   _REAL_ :: norm
 
   find_gradient(1:3) = ZERO
   atmfirst = bndatmptr(l-1)+1
   atmlast = bndatmptr(l)
   do jp = atmfirst, atmlast
      jatm = bndatmlst(jp) 
      dx = x-gcrd(1,jatm)
      dy = y-gcrd(2,jatm)
      dz = z-gcrd(3,jatm)
      disth = sqrt(dx**2+dy**2+dz**2)
      dist_inv = factor/disth
      dx = dx*dist_inv
      dy = dy*dist_inv
      dz = dz*dist_inv
      dist = (disth-radih(jatm))*factor
      find_gradient(1:3) = find_gradient(1:3) + density_1deriv_atom(dist,dx,dy,dz)
   end do

   ! normalize the vector

   norm = ONE/sqrt(sum(find_gradient(1:3)**2))
   find_gradient(1:3) = find_gradient(1:3)*norm

end function find_gradient

function density_atom(dist)
   implicit none

   _REAL_ :: density_atom
   _REAL_ :: dist

   integer m
   _REAL_, parameter :: dash(6) = (/0.00d0,0.20d0,0.40d0,0.60d0,0.80d0,1.00d0 /)
   _REAL_, parameter :: spcoef(5,4) = &
       reshape((/1.000000  ,0.2100000  ,0.1500000  ,0.0500000  ,0.010000   ,  &
                -4.527143  ,-2.067608  ,0.0475730  ,-0.522686  ,-0.056828  ,  &
                -3.640532  ,15.938209  ,-5.362303  ,2.5110050  ,-0.181716  ,  &
                32.631235  ,-35.500854 ,13.122180  ,-4.487867  ,1.079289/), (/5,4/))

   ! All dist and components are in the unit of solvent probe diameter

   density_atom = ZERO
   if ( dist > ONE ) then
   else if ( dist <= ZERO ) then
      density_atom = 1.0d0 - 4.527143d0*dist
   else
      do m = 1, 5
         if ( dist > dash(m) .and. dist <= dash(m+1) ) then
            density_atom = spcoef(m,1)                   + &
                           spcoef(m,2)*(dist-dash(m))    + &
                           spcoef(m,3)*(dist-dash(m))**2 + &
                           spcoef(m,4)*(dist-dash(m))**3
         end if
      end do
   end if

end function density_atom

function density_1deriv_atom(dist,dx,dy,dz)
   implicit none

   _REAL_, dimension(3) :: density_1deriv_atom
   _REAL_ :: dist, dx, dy, dz

   integer m
   _REAL_, parameter :: dash(6) = (/0.00d0,0.20d0,0.40d0,0.60d0,0.80d0,1.00d0 /)
   _REAL_, parameter :: spcoef(5,4) = &
       reshape((/1.000000  ,0.2100000  ,0.1500000  ,0.0500000  ,0.010000   ,  &
                -4.527143  ,-2.067608  ,0.0475730  ,-0.522686  ,-0.056828  ,  &
                -3.640532  ,15.938209  ,-5.362303  ,2.5110050  ,-0.181716  ,  &
                32.631235  ,-35.500854 ,13.122180  ,-4.487867  ,1.079289/), (/5,4/))
    
   ! All dist and components are in the unit of solvent probe diameter

   density_1deriv_atom(1:3) = ZERO
   if ( dist > ONE ) then
   else if ( dist <= ZERO ) then
      density_1deriv_atom(1) = -4.527143d0 * dx
      density_1deriv_atom(2) = -4.527143d0 * dy
      density_1deriv_atom(3) = -4.527143d0 * dz
   else
      do m = 1, 5
         if ( dist > dash(m) .and. dist <= dash(m+1) ) then
            density_1deriv_atom(1) = spcoef(m,2)                  *dx + &
                               TWO  *spcoef(m,3)*(dist-dash(m))   *dx + &
                               THREE*spcoef(m,4)*(dist-dash(m))**2*dx
            density_1deriv_atom(2) = spcoef(m,2)                  *dy + &
                               TWO  *spcoef(m,3)*(dist-dash(m))   *dy + &
                               THREE*spcoef(m,4)*(dist-dash(m))**2*dy
            density_1deriv_atom(3) = spcoef(m,2)                  *dz + &
                               TWO  *spcoef(m,3)*(dist-dash(m))   *dz + &
                               THREE*spcoef(m,4)*(dist-dash(m))**2*dz
         end if
      end do
   end if

end function density_1deriv_atom

function density_2deriv_atom(dist,disth,dx1,dx2)
   implicit none

   _REAL_,dimension(3) :: density_2deriv_atom
   _REAL_ :: dist, disth, dx1, dx2

   integer m
   _REAL_, parameter :: dash(6) = (/0.00d0,0.20d0,0.40d0,0.60d0,0.80d0,1.00d0 /)
   _REAL_, parameter :: spcoef(5,4) = &
       reshape((/1.000000  ,0.2100000  ,0.1500000  ,0.0500000  ,0.010000   ,  &
                -4.527143  ,-2.067608  ,0.0475730  ,-0.522686  ,-0.056828  ,  &
                -3.640532  ,15.938209  ,-5.362303  ,2.5110050  ,-0.181716  ,  &
                32.631235  ,-35.500854 ,13.122180  ,-4.487867  ,1.079289/), (/5,4/))

   ! All dist and components are in the unit of solvent probe diameter
   ! unfortunately the real distance (in h) has to be used as well here.

   _REAL_ :: factor, factor2
   _REAL_ :: disth_inv, disth_inv2, disth_inv3

   factor = HALF/dprobh
   factor2 = factor**2
   disth_inv = ONE/disth
   disth_inv2 = disth_inv**2
   disth_inv3 = disth_inv2*disth_inv

   density_2deriv_atom(1:3) = ZERO
   if      ( dist > ONE ) then
   else if ( dist <= ZERO ) then
      density_2deriv_atom(1) = 4.527143d0*dx1*dx1*factor*disth_inv3-&
                               4.527143d0*factor*disth_inv
      density_2deriv_atom(2) = 4.527143d0*dx2*dx2*factor*disth_inv3-&
                               4.527143d0*factor*disth_inv
      density_2deriv_atom(3) = 4.527143d0*dx1*dx2*factor*disth_inv3
   else
      do m = 1, 5
         if ( dist > dash(m) .and. dist <= dash(m+1) ) then
            density_2deriv_atom(1) = ( spcoef(m,2)+&
                                 TWO  *spcoef(m,3)*(dist-dash(m))+&
                                 THREE*spcoef(m,4)*(dist-dash(m))**2 )*&
                                               factor *disth_inv-&
                                     ( spcoef(m,2)+&
                                 TWO  *spcoef(m,3)*(dist-dash(m))+&
                                 THREE*spcoef(m,4)*(dist-dash(m))**2 )*&
                                       dx1*dx1*factor *disth_inv3+&
                               ( TWO  *spcoef(m,3)+&
                                 SIX  *spcoef(m,4)*(dist-dash(m)) )*&
                                       dx1*dx1*factor2*disth_inv2

            density_2deriv_atom(2) = ( spcoef(m,2)+&
                                 TWO  *spcoef(m,3)*(dist-dash(m))+&
                                 THREE*spcoef(m,4)*(dist-dash(m))**2 )*&
                                               factor *disth_inv-&
                                     ( spcoef(m,2)+&
                                 TWO  *spcoef(m,3)*(dist-dash(m))+&
                                 THREE*spcoef(m,4)*(dist-dash(m))**2 )*&
                                       dx2*dx2*factor *disth_inv3+&
                               ( TWO  *spcoef(m,3)+&
                                 SIX  *spcoef(m,4)*(dist-dash(m)) )*&
                                       dx2*dx2*factor2*disth_inv2

            density_2deriv_atom(3) =-( spcoef(m,2)+&
                                 TWO  *spcoef(m,3)*(dist-dash(m))+&
                                 THREE*spcoef(m,4)*(dist-dash(m))**2 )*&
                                       dx1*dx2*factor *disth_inv3+&
                               ( TWO  *spcoef(m,3)+&
                                 SIX  *spcoef(m,4)*(dist-dash(m)) )*&
                                       dx1*dx2*factor2*disth_inv2
         end if
      end do
   end if

end function density_2deriv_atom

end subroutine ana_irreg_grd

end module ana_iim
