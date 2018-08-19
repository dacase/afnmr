! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"

module iim_util

#  include "pb_constants.h"

   _REAL_, allocatable :: wp(:),qp(:)
   _REAL_, allocatable :: qyp(:),qzp(:)
   _REAL_, allocatable :: wyp(:),wzp(:)
   _REAL_, allocatable :: wyyp(:),wzzp(:),wyzp(:)

contains

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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ IIM algoritm interface routine
subroutine pb_iimdrv( npbstep,npbgrid,nstlim,atmfirst,atmlast,npbopt,solvopt,level,nfocus,bcopt,&
                      natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                      xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                      maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
                      pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
                      gcrd,acrg,&
                      nbnd,iepsav,insas,lvlset,&
                      chgrd,saltgrd,phi,&
                      bv,cphi,xs&
                    )

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

   ! for singularity-free PB equation we need to compute cphi() at the
   ! dielectric interface grid points

   if ( bcopt > 5 .and. bcopt < 10 ) then
      cphi(1:xmymzm) = ZERO
      call pb_dbcgrd( cphi(1), insas )
   end if

   ! for singular PB equation bv() is initialized as the grid charges

   if ( bcopt < 6 ) then
      bv = chgrd
   end if

   ! now we put effective charges at the space boundary grid points into bv()
   ! to take care of the boundary conditions

   call pb_bndcnd( bv(1), chgrd(1) )

   call iim(gox,goy,goz,xm,ym,zm,lvlset,insas,nbnd,iepsav, &
            epsin/eps0,epsout/eps0, &
            bv(1),phi(1),xs(1),h,atmfirst,atmlast,bcopt,&
            maxitn,itn,accept,norm,inorm)

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_dbcgrd( cphi, insas )

   _REAL_ cphi(xm,ym,zm)
   integer insas(0:xm+1,0:ym+1,0:zm+1)

   integer i, j, k
   integer i0, j0, k0
   integer ip
   _REAL_ tmp

   tmp = ZERO

   do ip = 1, nbnd
      i0 = iepsav(1,ip); j0 = iepsav(2,ip); k0 = iepsav(3,ip)

      i = i0 - 1
      if ( insas(i ,j0,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i ,j0,k0) == ZERO ) then
            call get_coulpot(i ,j0,k0,tmp)
            cphi(i ,j0,k0) = tmp/epsin
         end if
      end if

      i = i0 + 1
      if ( insas(i ,j0,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i ,j0,k0) == ZERO ) then
            call get_coulpot(i ,j0,k0,tmp)
            cphi(i ,j0,k0) = tmp/epsin
         end if
      end if

      j = j0 - 1
      if ( insas(i0,j ,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j ,k0) == ZERO ) then
            call get_coulpot(i0,j ,k0,tmp);
            cphi(i0,j ,k0) = tmp/epsin
         end if
      end if

      j = j0 + 1
      if ( insas(i0,j ,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j ,k0) == ZERO ) then
            call get_coulpot(i0,j ,k0,tmp)
            cphi(i0,j ,k0) = tmp/epsin
         end if
      end if

      k = k0 - 1
      if ( insas(i0,j0,k ) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j0,k ) == ZERO ) then
            call get_coulpot(i0,j0,k ,tmp)
            cphi(i0,j0,k ) = tmp/epsin
         end if
      end if

      k = k0 + 1
      if ( insas(i0,j0,k ) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j0,k ) == ZERO ) then
            call get_coulpot(i0,j0,k ,tmp)
            cphi(i0,j0,k ) = tmp/epsin
         end if
      end if

      if ( insas(i0,j0,k0) > 0 .or. bcopt > 7 ) then
         if ( cphi(i0,j0,k0) == ZERO ) then
            call get_coulpot(i0,j0,k0,tmp)
            cphi(i0,j0,k0) = tmp/epsin
         end if
      end if

   end do

end subroutine pb_dbcgrd
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_coulpot(i,j,k,pot)

   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green

   integer i,j,k
   _REAL_ pot

   integer iatm
   integer itmp,jtmp,ktmp
   integer idx,idy,idz

   _REAL_ factor,qtmp,rinv,xtmp,ytmp,ztmp,dx,dy,dz
   _REAL_ a,a1,b,b1,c,c1

   factor = ONE/(FOURPI)/h

   pot = ZERO
   do iatm = atmfirst, atmlast
      xtmp = gcrd(1,iatm); ytmp = gcrd(2,iatm); ztmp = gcrd(3,iatm)
      qtmp = factor*acrg(iatm)

      dx = abs(i-xtmp); dy = abs(j-ytmp); dz = abs(k-ztmp)
      if (dx < 40.d0 .and. dy < 40.d0 .and. dz < 40.d0) then
         idx = floor(dx); idy = floor(dy); idz = floor(dz)
         a=dx-idx;b=dy-idy;c=dz-idz
         a1 = 1 - a; b1 = 1 - b; c1 = 1 - c
         rinv = a1*b1*c1*green(idx  ,idy  ,idz  ) &
               +a *b1*c1*green(idx+1,idy  ,idz  ) &
               +a1*b *c1*green(idx  ,idy+1,idz  ) &
               +a *b *c1*green(idx+1,idy+1,idz  ) &
               +a1*b1*c *green(idx  ,idy  ,idz+1) &
               +a *b1*c *green(idx+1,idy  ,idz+1) &
               +a1*b *c *green(idx  ,idy+1,idz+1) &
               +a *b *c *green(idx+1,idy+1,idz+1)
      else
         rinv = ONE/sqrt(dble(dx**2 + dy**2 + dz**2))
      end if
      pot = pot + qtmp*rinv
   end do  !  iatm = atmfirst, atmlast

end subroutine get_coulpot
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign Debye-Huckel potential for the boundary charge grid
subroutine pb_bndcnd( bv, chgrd )

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
    bv=0.0d0 !XP: This can make imin=6 work with ipb4
   if ( level == 1 .and. bcopt == 1 ) then
      write(6, *) "PB bomb in pb_bndcnd(): zero BC not supported"
      call mexit(6, 1)

   ! bcopt = 2
   ! molecule dipolar debye-huckel potential in the singular PB.
   ! the boundary will be all solvent.

   else if ( level == 1 .and. bcopt == 2 ) then
      write(6, *) "PB bomb in pb_iimdrv(): molecular dipolar BC not supported"
      call mexit(6, 1)

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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ This is the real IIM driver
subroutine iim(gox,goy,goz,l,m,n,lvlset,insas,nbnd,iepsav,epsin,epsout,bv,u,u0,h,&
               atmfirst,atmlast,bcopt,maxitn,itn,accept,norm,inorm)

!    xs <- gox, ys <- goy, zs <- goz
!    l  <- xm , m  <- ym , n  <- zm
!    bi <- epsin, bo <- epsout
!    u  <- phi, u0 <- xs

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
   integer nirreg
   integer i,j,k,IFAIL,IFAIL2,ii,ir,nc,nz_num
   _REAL_ xs,ys,zs,xf,yf,zf,hx,hy,hz,hmax,bi,bo
   _REAL_ beta_max,rhs,err
   _REAL_ coe1(7), coe2(27)

   integer it,jt,kt,iit,it1,jt1,kt1,iit1
   _REAL_,allocatable :: q0p(:),w0p(:)
   _REAL_,allocatable :: wcoe(:,:),wxcoe(:, :),wycoe(:, :)
   _REAL_,allocatable :: wzcoe(:,:),wxxcoe(:, :),wyycoe(:, :)
   _REAL_,allocatable :: wzzcoe(:, :),wxycoe(:, :)
   _REAL_,allocatable :: wxzcoe(:, :),wyzcoe(:, :)
   _REAL_,allocatable :: sss1(:),sss2(:)
   _REAL_,allocatable :: c(:,:,:,:), c2(:, :)
   _REAL_,allocatable :: f(:,:,:)

   ! interfacing variables

   nirreg = nbnd
   xs = gox; ys = goy; zs = goz
   hx = h; hy = h; hz = h; hmax = h
   bi = epsin; bo = epsout

   ! allocating working arrays

   allocate (wp(nbnd), qp(nbnd))  ! jump conditions at irregular points
   allocate (q0p(nbnd),qyp(nbnd),qzp(nbnd)) ! tangential derivatives of field jump conditions
   allocate (w0p(nbnd),wyp(nbnd),wzp(nbnd)) ! tangential derivatives of potential jump conditions
   allocate (wyyp(nbnd),wzzp(nbnd),wyzp(nbnd)) ! second derivatives of potential jump conditions
   allocate (wcoe(nbnd,nq),wxcoe(nbnd, nq),wycoe(nbnd, nq)) ! coefficients obtained by SVD
   allocate (wzcoe(nbnd,nq),wxxcoe(nbnd, nq),wyycoe(nbnd, nq))
   allocate (wzzcoe(nbnd, nq),wxycoe(nbnd, nq))
   allocate (wxzcoe(nbnd, nq),wyzcoe(nbnd, nq))
   allocate (sss1(nbnd),sss2(nbnd)) !dummy arrays
   allocate (c(l,m,n,7), c2(nbnd, 27))
   allocate (f(l,m,n))

   ! the surface related infrastructure data from density_surface module is used

   ! setting up the jump conditions at the projection points of the irreular (boundary) grid points

   call prodis(l,m,n,nirreg,bcopt,bi,bo,atmfirst,atmlast,nbnd,cirreg)

   ! calculate the first and second derivatives of the jump conditions in the
   ! surface tangential directions in the local coordinate system
   ! step 1: this is the first derivatives

   call coed20(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,nbnd,nq,index2,cirreg, &
               wcoe,wxcoe,wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)
   call qint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg, &
             wcoe,wycoe,wzcoe,q0p)
   call wint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,wp, &
             wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,wyp,wzp,wyyp,wzzp,wyzp)

   ! step 2: this is the second derivatives

   call coed6(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,nbnd,nq,index2,cirreg, &
              wcoe,wxcoe,wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)
   call wint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,wp, &
             wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,sss1,sss2,wyyp,wzzp,wyzp)

   ! setting up linear system coefficient matrix

   beta_max = max(bi,bo)
   nz_num = 0
   do k=1, n
   do j=1, m
   do i=1, l
      if ( index(i,j,k).eq.1 ) then
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
         ii = i+(j-1)*l+(k-1)*l*m
         f(i,j,k)   = bv(i,j,k)/h/h/h*FOURPI
         do ii = 1 , 7
            if ( abs(c(i,j,k,ii)) > 1.d-10 ) nz_num=nz_num+1
         end do
      else if (index(i,j,k).eq.5 ) then
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
         ii = i+(j-1)*l+(k-1)*l*m
         f(i,j,k) = bv(i,j,k)/h/h/h*FOURPI
         do ii = 1, 7
            if ( c(i,j,k,ii ) > 1.d-10 ) nz_num=nz_num+1
         end do
      else
         ii = i+(j-1)*l+(k-1)*l*m
         call irre31(l,m,n,h,hx,hy,hz,IFAIL, &
                     i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,lvlset,index, &
                     nq,nbnd,index2,cirreg,coe2,rhs)

         if (IFAIL.gt.10) then
            call irre32(l,m,n,h,hx,hy,hz,IFAIL2, &
                        i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,lvlset,index, &
                        q0p,w0p, &
                        nq,nbnd,index2,cirreg,coe2,rhs)
         end if

         ir = index2(i,j,k)

         do nc = 1, 27
            c2(ir,nc) = coe2(nc)
         end do
         f(i,j,k) = rhs
         do ii = 1, 27
            if ( abs(c2(ir,ii)) > 1.d-10 ) nz_num=nz_num+1
         end do
      endif
   end do
   end do
   end do

   ! cleaning up working arrays
   ! not the jump conditions which will be used by other routines

   deallocate (q0p,w0p)
   deallocate (wcoe,wxcoe,wycoe) ! coefficients obtained by SVD
   deallocate (wzcoe,wxxcoe,wyycoe)
   deallocate (wzzcoe,wxycoe)
   deallocate (wxzcoe,wyzcoe)
   deallocate (sss1,sss2) !dummy arrays

   ! entering the linear system solver

   !if ( solvopt == 1 ) then
   !   ! algebraic multigrid
   !   ! call amg(l,m,n,nbnd,c,c2,index,index2,f,u,u0,accept) ! should be modified
   !   write(6,'(a)') "AMG Driver is retired."; call mexit(6,1)
   !else if ( solvopt == 2 ) then
   !   ! ILU preconditioned GMRES
   !   call gmres(l,m,n,nbnd,nz_num,c,c2,index,index2,f,u,u0,accept)
   !else if ( solvopt == 3 ) then
   ! Only allows the ILU preconditioned BiCG
   call bicg(l,m,n,nbnd,nz_num,c,c2,index,index2,f,u,u0,maxitn,itn,accept,err)
   !else
   !   write(6,'(a)') "Unsupported Solver option"; call mexit(6,1)
   !end if

   ! returning error related data

   inorm = sum(abs(f))
   norm = inorm*err

   ! converting back to the PBSA unit for potential

   do k=1, n
   do j=1, m
   do i=1, l
      u(i,j,k) = u(i,j,k) * INV_FOURPI / eps0
   end do
   end do
   end do

   deallocate (c, c2)
   deallocate (f)

end subroutine iim

end subroutine pb_iimdrv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine prodis here]
subroutine prodis(l,m,n,nirreg,bcopt,bi,bo,atmfirst,atmlast,maxirr,cirreg)

   implicit none

   ! passed variables

   integer l, m, n, nirreg
   integer bcopt, atmfirst, atmlast, maxirr
   _REAL_ cirreg(15, maxirr)
   _REAL_ t(3,3)
   _REAL_ bi,bo

   ! local variables

   integer ir
   _REAL_ xx, yy, zz

   ! external functions

   do ir = 1, nirreg
      xx     = cirreg( 1, ir)
      yy     = cirreg( 2, ir)
      zz     = cirreg( 3, ir)
      t(1,1) = cirreg( 7, ir)
      t(1,2) = cirreg( 8, ir)
      t(1,3) = cirreg( 9, ir)
      t(2,1) = cirreg(10, ir)
      t(2,2) = cirreg(11, ir)
      t(2,3) = cirreg(12, ir)
      t(3,1) = cirreg(13, ir)
      t(3,2) = cirreg(14, ir)
      t(3,3) = cirreg(15, ir)
      wp(ir) = fw(bcopt,bi,bo,atmfirst,atmlast,xx,yy,zz)
      qp(ir) = fq(bcopt,bi,bo,atmfirst,atmlast,xx,yy,zz,t)
   end do

end subroutine prodis
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine qint here]
subroutine qint(l,m,n,h,hx,hy,hz,xs,ys,zs,maxirr,nq,index,index2,cirreg, &
              wcoe,wycoe,wzcoe,q0p)

   implicit none

   ! Passed variables

   integer l,m,n,maxirr,nq
   _REAL_  h,hx,hy,hz,xs,ys,zs
   integer index(l,m,n), index2(l,m,n)
   _REAL_  cirreg(15, maxirr)
   _REAL_  wcoe(maxirr,nq),wycoe(maxirr, nq),wzcoe(maxirr, nq)
   _REAL_  q0p(maxirr)

   !Local variables

   integer i,j,k,ir
   _REAL_  q,qy,qz

   do i=1, l
   do j=1, m
   do k=1, n
      if (index(i,j,k) > 1 .and. index(i,j,k) < 5) then
         ir = index2(i, j, k)
         call qint2(l,m,n,h,hx,hy,hz,xs,ys,zs,ir,maxirr,nq,index2,cirreg,qp,&
            wcoe,wycoe,wzcoe,q,qy,qz)
         q0p(ir) = q
         qyp(ir) = qy
         qzp(ir) = qz
      end if
   end do
   end do
   end do

end subroutine qint
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine qint2 here]
subroutine qint2(l,m,n,h,hx,hy,hz,xs,ys,zs,ir, maxirr,nq,index2,cirreg,qp,&
              wcoe,wycoe,wzcoe,q,qy,qz)

   implicit none
   integer,parameter :: nqq = 27

   ! Passed variables

   integer l,m,n,ir,maxirr,nq
   _REAL_  xs,ys,zs,h,hx,hy,hz
   integer index2(l,m,n)
   _REAL_  cirreg(15, maxirr)
   _REAL_  wcoe(maxirr,nq),wycoe(maxirr, nq),wzcoe(maxirr, nq)
   _REAL_  qp(maxirr),t(3,3)
   _REAL_  q,qy,qz
   _REAL_  itemp(nqq)

   !Local variables

   integer if0,i0,j0,k0,nsub,i1,j1,k1,i2,j2,k2,idis,ii,ip,nn2
   _REAL_  hmax,alf,x1,y1,z1,xyy,xzz,xyz,x2,y2,z2,dis,x3,y3,z3

   hmax=h
   if (nq /= nqq) then
      write(*,*) "   Please change nqq in qint2() s.t. nqq=nq in main()!"
      stop
   end if

   alf = 10.1*hmax
   if0 = int(alf/hmax) + 1

   q  = 0.0
   qy = 0.0
   qz = 0.0

   x1     = cirreg( 1, ir)
   y1     = cirreg( 2, ir)
   z1     = cirreg( 3, ir)

   xyy    = cirreg( 4, ir)
   xzz    = cirreg( 5, ir)
   xyz    = cirreg( 6, ir)

   t(1,1) = cirreg( 7, ir)
   t(1,2) = cirreg( 8, ir)
   t(1,3) = cirreg( 9, ir)
   t(2,1) = cirreg(10, ir)
   t(2,2) = cirreg(11, ir)
   t(2,3) = cirreg(12, ir)
   t(3,1) = cirreg(13, ir)
   t(3,2) = cirreg(14, ir)
   t(3,3) = cirreg(15, ir)

   i0 = nint((x1-xs)/hx)
   j0 = nint((y1-ys)/hy)
   k0 = nint((z1-zs)/hz)

   nsub = 0
   do i2 = 0, if0          !# Starting from the center
      do i1 = i0-i2, i0+i2
      do j1 = j0-i2, j0+i2
      do k1 = k0-i2, k0+i2
         if (i1 < 1 .or. j1 < 1 .or. k1 < 1 .or. i1 > l .or. j1 > m .or. k1 > n) goto 555

         idis = abs(i1-i0) + abs(j1-j0) + abs(k1-k0)
         nn2 = index2(i1,j1,k1)

         if (idis == i2 .and. nn2 > 0) then
            x2 = cirreg(1, nn2)
            y2 = cirreg(2, nn2)
            z2 = cirreg(3, nn2)

            do ii = 1, nsub
               ip = itemp(ii)
               x3 = cirreg(1, ip)
               y3 = cirreg(2, ip)
               z3 = cirreg(3, ip)
               call distan(x2,y2,z2,x3,y3,z3, dis)
               if (dis < 0.1*hmax) then
                  goto 555
               end if
            end do

            call distan(x1,y1,z1,x2,y2,z2, dis)
            if ( dis < alf ) then
               nsub = nsub + 1
               itemp(nsub) = nn2
               ip = index2(i1,j1,k1)
               q  = q + wcoe(ir, nsub)*qp(ip)
               qy = qy + wycoe(ir, nsub)*qp(ip)
               qz = qz + wzcoe(ir, nsub)*qp(ip)
               if ( nsub .eq. nq) goto 577
            end if
         end if

         555 continue
      end do  ! k1 = k0-i2, k0+i2
      end do  ! j1 = j0-i2, j0+i2
      end do  ! i1 = i0-i2, i0+i2
   end do  ! i2=0, if0
   577 continue

end subroutine qint2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine wint here]
subroutine wint(l,m,n,h,hx,hy,hz,xs,ys,zs,maxirr,nq,index,index2,cirreg,wp,&
              wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,wyp,wzp,wyyp,wzzp,wyzp)

   implicit none

   ! Passed variables

   integer l,m,n,maxirr,nq
   _REAL_  xs,ys,zs,hx,hy,hz,h
   integer index(l,m,n), index2(l,m,n)
   _REAL_  cirreg(15, maxirr)
   _REAL_  wcoe(maxirr,nq),wycoe(maxirr, nq),wzcoe(maxirr, nq)
   _REAL_  wyycoe(maxirr, nq)
   _REAL_  wzzcoe(maxirr, nq), wyzcoe(maxirr, nq)
   _REAL_  wp(maxirr)
   _REAL_  w0p(maxirr),wyp(maxirr),wzp(maxirr)
   _REAL_  wyyp(maxirr),wzzp(maxirr),wyzp(maxirr)

   ! Local variables

   integer i,j,k,ir
   _REAL_  w0,wy,wz,wyy,wzz,wyz

   do i=1,l
   do j=1,m
   do k=1,n
      if (index(i,j,k) > 1 .and. index(i,j,k) < 5) then
         ir = index2(i, j, k)
         call wint2(l,m,n,h,hx,hy,hz,xs,ys,zs,ir,maxirr,nq,index2,cirreg,wp,&
            wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0,wy,wz,wyy,wzz,wyz)
         w0p(ir)  = w0
         wyp(ir)  = wy
         wzp(ir)  = wz
         wyyp(ir) = wyy
         wzzp(ir) = wzz
         wyzp(ir) = wyz
      end if
   end do
   end do
   end do

end subroutine wint
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine wint2 here]
subroutine wint2(l,m,n,h,hx,hy,hz,xs,ys,zs,ir, maxirr,nq,index2,cirreg, wp, &
      wcoe,wycoe,wzcoe, wyycoe,wzzcoe,wyzcoe, &
      w0,wy,wz,wyy,wzz,wyz)

   implicit none

   ! Passed variables

   integer, parameter :: nqq = 27
   integer l,m,n,ir,maxirr,nq
   _REAL_  xs,ys,zs,h,hx,hy,hz
   integer index2(l,m,n)
   _REAL_  cirreg(15, maxirr)
   _REAL_  wcoe(maxirr,nq),wycoe(maxirr, nq),wzcoe(maxirr, nq)
   _REAL_  wyycoe(maxirr, nq)
   _REAL_  wzzcoe(maxirr, nq), wyzcoe(maxirr, nq)
   _REAL_  wp(maxirr)
   _REAL_  w0,wy,wz,wyy,wzz,wyz

   ! Local variables

   integer if0,i0,j0,k0,nsub,i1,j1,k1,i2,j2,k2,ip,ii,idis,nn2
   _REAL_  hmax,alf,x1,y1,z1,x2,y2,z2,x3,y3,z3,dis,xyy,xzz,xyz
   _REAL_  t(3,3)
   integer itemp(nqq)

   hmax=h
   if (nq /= nqq) then
      write(*,*) "   Please change nqq in wint2() s.t. nqq=nq in main()!"
      stop
   end if

   alf = 10.1*hmax
   if0 = int(alf/hmax) + 1

   w0  = 0.0
   wy  = 0.0
   wz  = 0.0
   wyy = 0.0
   wzz = 0.0
   wyz = 0.0

   x1     = cirreg( 1, ir)
   y1     = cirreg( 2, ir)
   z1     = cirreg( 3, ir)

   xyy    = cirreg( 4, ir)
   xzz    = cirreg( 5, ir)
   xyz    = cirreg( 6, ir)

   t(1,1) = cirreg( 7, ir)
   t(1,2) = cirreg( 8, ir)
   t(1,3) = cirreg( 9, ir)
   t(2,1) = cirreg(10, ir)
   t(2,2) = cirreg(11, ir)
   t(2,3) = cirreg(12, ir)
   t(3,1) = cirreg(13, ir)
   t(3,2) = cirreg(14, ir)
   t(3,3) = cirreg(15, ir)

   i0 = nint((x1-xs)/hx)
   j0 = nint((y1-ys)/hy)
   k0 = nint((z1-zs)/hz)

   nsub = 0
   do i2 = 0, if0          !# Starting from the center
      do i1 = i0-i2, i0+i2
      do j1 = j0-i2, j0+i2
      do k1 = k0-i2, k0+i2
         if (i1 < 1 .or. j1 < 1 .or. k1 < 1 .or. i1 > l .or. j1 > m .or. k1 > n) goto 555

         idis = abs(i1-i0) + abs(j1-j0) + abs(k1-k0)
         nn2 = index2(i1,j1,k1)

         if (idis == i2 .and. nn2 > 0) then
            x2 = cirreg(1, nn2)
            y2 = cirreg(2, nn2)
            z2 = cirreg(3, nn2)

            do ii = 1, nsub
               ip = itemp(ii)
               x3 = cirreg(1, ip)
               y3 = cirreg(2, ip)
               z3 = cirreg(3, ip)
               call distan(x2,y2,z2,x3,y3,z3, dis)
               if (dis < 0.1*hmax) then
                  goto 555
               end if
            end do

            call distan(x1,y1,z1,x2,y2,z2, dis)
            if ( dis < alf ) then
               nsub = nsub + 1
               itemp(nsub) = nn2
               ip = index2(i1,j1,k1)
               w0  = w0 + wcoe(ir, nsub)*wp(ip)
               wy = wy + wycoe(ir, nsub)*wp(ip)
               wz = wz + wzcoe(ir, nsub)*wp(ip)
               wyy = wyy + wyycoe(ir, nsub)*wp(ip)
               wzz = wzz + wzzcoe(ir, nsub)*wp(ip)
               wyz = wyz + wyzcoe(ir, nsub)*wp(ip)
               if ( nsub .eq. nq) goto 577
            end if
         end if
         555 continue
      end do  ! k1 = k0-i2, k0+i2
      end do  ! j1 = j0-i2, j0+i2
      end do  ! i1 = i0-i2, i0+i2
   end do  ! i2=0, if0
   577 continue

end subroutine wint2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine coed20 here]
subroutine coed20(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,maxirr,nq,index2,cirreg, &
              wcoe, wxcoe, wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)

   implicit none

   ! Passed variables

   integer l,m,n,nirreg,maxirr,nq
   integer index2(l,m,n)
   _REAL_  xs,ys,zs,h,hx,hy,hz
   _REAL_  cirreg(15, maxirr)
   _REAL_  wcoe(maxirr,nq),wxcoe(maxirr, nq),wycoe(maxirr, nq)
   _REAL_  wzcoe(maxirr,nq),wxxcoe(maxirr, nq),wyycoe(maxirr, nq)
   _REAL_  wzzcoe(maxirr, nq),wxycoe(maxirr, nq)
   _REAL_  wxzcoe(maxirr, nq),wyzcoe(maxirr, nq)

   ! Local variables

   integer,parameter :: nqq = 27
   integer i0,j0,k0,ir,isvd,i,j,k,i1,j1,k1,ip,job,ii,&
           ij,jj,nsub,nn2,i2,if0,idis,inf
   _REAL_  alf,x1,y1,z1,x2,y2,z2,x3,y3,z3,xyy,xzz,xyz,dis,hmax
   _REAL_  w1(10, nqq), w4(10)
   _REAL_  sd(11), uw(10,10), v(nqq, nqq)
   _REAL_  ew(nqq)
   _REAL_  ew1(nqq), ew2(nqq), ew3(nqq)
   _REAL_  ew4(nqq), ew5(nqq), ew6(nqq)
   _REAL_  ew7(nqq), ew8(nqq), ew9(nqq), ew10(nqq)
   _REAL_  w31(nqq), w32(nqq), w33(nqq)
   _REAL_  w34(nqq), w35(nqq), w36(nqq)
   _REAL_  w37(nqq), w38(nqq), w39(nqq), w310(nqq)
   _REAL_  sd2(10), work(1000)
   _REAL_  t(3,3)
   _REAL_  tempx(3),tempy(3)
   integer itemp(nqq)

   hmax=h
   isvd = 1

   if (nq /= nqq) then
      write(*,*) "   Please change nqq in coedef() s.t. nqq=nq in main()!"
      stop
   end if

   alf = 18.1*hmax
   if0 = int(alf/hmax) + 1

   do i=1, nirreg
      do j=1, nq
         wcoe(i,j)   = 0.0
         wxcoe(i,j)  = 0.0
         wycoe(i,j)  = 0.0
         wzcoe(i,j)  = 0.0
         wxxcoe(i,j) = 0.0
         wyycoe(i,j) = 0.0
         wzzcoe(i,j) = 0.0
         wxycoe(i,j) = 0.0
         wxzcoe(i,j) = 0.0
         wyzcoe(i,j) = 0.0
      end do
   end do

   ! ----- for each projection point

   do ir = 1, nirreg

      x1     = cirreg( 1, ir)
      y1     = cirreg( 2, ir)
      z1     = cirreg( 3, ir)

      xyy    = cirreg( 4, ir)
      xzz    = cirreg( 5, ir)
      xyz    = cirreg( 6, ir)

      t(1,1) = cirreg( 7, ir)
      t(1,2) = cirreg( 8, ir)
      t(1,3) = cirreg( 9, ir)
      t(2,1) = cirreg(10, ir)
      t(2,2) = cirreg(11, ir)
      t(2,3) = cirreg(12, ir)
      t(3,1) = cirreg(13, ir)
      t(3,2) = cirreg(14, ir)
      t(3,3) = cirreg(15, ir)

      i0 = nint((x1-xs)/hx)
      j0 = nint((y1-ys)/hy)
      k0 = nint((z1-zs)/hz)

      ! -------- initialize

      do jj=1, nq
         do ii=1, nq
            v(ii, jj) = 0.0      !# the V matrix in SVD: U*S*V
         end do
         do ii=1, 20
            w1(ii, jj) = 0.0     !# the coefficient matrix to be SVDed
         end do
      end do

      ! -------- Generate the matrix

      nsub = 0
      do i2=0, if0          !# Starting from the center
         do i1 = i0-i2, i0+i2
         do j1 = j0-i2, j0+i2
         do k1 = k0-i2, k0+i2
            if (i1 < 1 .or. j1 < 1 .or. k1 < 1 .or. i1 > l .or. j1 > m .or. k1 > n) goto 55

            idis = abs(i1-i0) + abs(j1-j0) + abs(k1-k0)
            nn2 = index2(i1,j1,k1)
            if (idis == i2 .and. nn2 > 0) then

               x2 = cirreg(1, nn2)
               y2 = cirreg(2, nn2)
               z2 = cirreg(3, nn2)
               do ii = 1, nsub
                  ip = itemp(ii)
                  x3 = cirreg(1, ip)
                  y3 = cirreg(2, ip)
                  z3 = cirreg(3, ip)
                  call distan(x2,y2,z2,x3,y3,z3, dis)
                  if (dis < 0.1*hmax) then
                     goto 55
                  end if
               end do

               call distan(x1,y1,z1,x2,y2,z2, dis)
               if ( dis < alf ) then
                  nsub = nsub + 1
                  itemp(nsub) = nn2
                  tempx(1) = x2 - x1
                  tempx(2) = y2 - y1
                  tempx(3) = z2 - z1
                  call matvec(3,3,t, tempx, tempy)
                  w1(1, nsub)  = 1.0
                  w1(2, nsub)  = tempy(1)
                  w1(3, nsub)  = tempy(2)
                  w1(4, nsub)  = tempy(3)
                  w1(5, nsub)  = 0.5*tempy(1) * tempy(1)
                  w1(6, nsub)  = 0.5*tempy(2) * tempy(2)
                  w1(7, nsub)  = 0.5*tempy(3) * tempy(3)
                  w1(8, nsub)  = tempy(1) * tempy(2)
                  w1(9, nsub)  = tempy(1) * tempy(3)
                  w1(10, nsub) = tempy(2) * tempy(3)
                  if ( nsub .eq. nq) goto 77
               end if
            end if  ! (idis == i2 .and. nn2 > 0)
            55 continue
         end do  ! k1 = k0-i2, k0+i2
         end do  ! j1 = j0-i2, j0+i2
         end do ! i1 = i0-i2, i0+i2
      end do ! i2=0, if0
      77 continue

      if (nsub < 12) then
         write(*,*) "   nsub = ",nsub,"alf is too small in coede20f()!"
         stop
      end if

      ! -------- Call least square routine

      if (isvd == 1) then
         job = 11
         call dsvdc(w1,10,10,nsub,sd,ew,uw,10,v,nq,w4,job,inf)
      else
         call dgesvd('A','A',10,nsub,w1,10,sd2,uw,10,v,nq, &
               work, 1000, inf)

         do ij = 1, 10
            sd(ij) = sd2(ij)
         end do
      end if

      if (inf /= 0) then
         write(*,*) inf, " - ssvdc() or sgesvd() failed in coedef!"
         stop
      end if

      do i1=1, 10
         if (abs(sd(i1)) > 1.0e-14 ) then
            ew1(i1)=uw(1, i1)/sd(i1)
            ew2(i1)=uw(2, i1)/sd(i1)
            ew3(i1)=uw(3, i1)/sd(i1)
            ew4(i1)=uw(4, i1)/sd(i1)
            ew5(i1)=uw(5, i1)/sd(i1)
            ew6(i1)=uw(6, i1)/sd(i1)
            ew7(i1)=uw(7, i1)/sd(i1)
            ew8(i1)=uw(8, i1)/sd(i1)
            ew9(i1)=uw(9, i1)/sd(i1)
            ew10(i1)=uw(10, i1)/sd(i1)
         else
            ew1(i1)= 0.0
            ew2(i1)= 0.0
            ew3(i1)= 0.0
            ew4(i1)= 0.0
            ew5(i1)= 0.0
            ew6(i1)= 0.0
            ew7(i1)= 0.0
            ew8(i1)= 0.0
            ew9(i1)= 0.0
            ew10(i1)= 0.0
         end if
      end do

      ! -------- w31(i),w32(i),......,w36(i) are the solutions

      do i1=1, nsub
         w31(i1) = 0.0
         w32(i1) = 0.0
         w33(i1) = 0.0
         w34(i1) = 0.0
         w35(i1) = 0.0
         w36(i1) = 0.0
         w37(i1) = 0.0
         w38(i1) = 0.0
         w39(i1) = 0.0
         w310(i1) = 0.0

         if (isvd == 1) then
            do j1=1,10
               w31(i1) = w31(i1) + v(i1,j1)*ew1(j1)
               w32(i1) = w32(i1) + v(i1,j1)*ew2(j1)
               w33(i1) = w33(i1) + v(i1,j1)*ew3(j1)
               w34(i1) = w34(i1) + v(i1,j1)*ew4(j1)
               w35(i1) = w35(i1) + v(i1,j1)*ew5(j1)
               w36(i1) = w36(i1) + v(i1,j1)*ew6(j1)
               w37(i1) = w37(i1) + v(i1,j1)*ew7(j1)
               w38(i1) = w38(i1) + v(i1,j1)*ew8(j1)
               w39(i1) = w39(i1) + v(i1,j1)*ew9(j1)
               w310(i1) = w310(i1) + v(i1,j1)*ew10(j1)
            end do
         else
            do j1=1,10
               w31(i1) = w31(i1) + v(j1,i1)*ew1(j1)
               w32(i1) = w32(i1) + v(j1,i1)*ew2(j1)
               w33(i1) = w33(i1) + v(j1,i1)*ew3(j1)
               w34(i1) = w34(i1) + v(j1,i1)*ew4(j1)
               w35(i1) = w35(i1) + v(j1,i1)*ew5(j1)
               w36(i1) = w36(i1) + v(j1,i1)*ew6(j1)
               w37(i1) = w37(i1) + v(j1,i1)*ew7(j1)
               w38(i1) = w38(i1) + v(j1,i1)*ew8(j1)
               w39(i1) = w39(i1) + v(j1,i1)*ew9(j1)
               w310(i1) = w310(i1) + v(j1,i1)*ew10(j1)
            end do
         end if

         wcoe(ir,i1)   = w31(i1)
         wxcoe(ir,i1)  = w32(i1)
         wycoe(ir,i1)  = w33(i1)
         wzcoe(ir,i1)  = w34(i1)
         wxxcoe(ir,i1) = w35(i1)
         wyycoe(ir,i1) = w36(i1)
         wzzcoe(ir,i1) = w37(i1)
         wxycoe(ir,i1) = w38(i1)
         wxzcoe(ir,i1) = w39(i1)
         wyzcoe(ir,i1) = w310(i1)
      end do  !  i1=1, nsub

   end do  ! ir = 1, nirreg

end subroutine coed20
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine coed6 here]
subroutine coed6(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,maxirr,nq,index2,cirreg, &
              wcoe, wxcoe, wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)

   implicit none

   !Passed variables
   integer l, m, n, nirreg,maxirr,nq
   integer index2(l,m,n)
   _REAL_  xs,ys,zs,h,hx,hy,hz
   _REAL_  cirreg(15, maxirr)
   _REAL_  wcoe(maxirr,nq),wxcoe(maxirr, nq),wycoe(maxirr, nq)
   _REAL_  wzcoe(maxirr,nq),wxxcoe(maxirr, nq),wyycoe(maxirr, nq)
   _REAL_  wzzcoe(maxirr, nq),wxycoe(maxirr, nq)
   _REAL_  wxzcoe(maxirr, nq),wyzcoe(maxirr, nq)

   ! Local variables

   integer,parameter :: nqq = 27
   integer isvd,alf,i0,j0,k0,ip,i,j,k,i1,j1,k1,i2,j2,k2,&
           inf,ij,ii,jj,job,nn2,nsub,ir,if0
   _REAL_  hmax,x1,y1,z1,x2,y2,z2,x3,y3,z3,dis,&
           idis,xyy,xyz,xzz
   _REAL_  w1(6, nqq), w4(6)
   _REAL_  sd(7), uw(6,6), v(nqq, nqq)
   _REAL_  ew(nqq)
   _REAL_  ew1(nqq), ew2(nqq), ew3(nqq)
   _REAL_  ew4(nqq), ew5(nqq), ew6(nqq)
   _REAL_  w31(nqq), w32(nqq), w33(nqq)
   _REAL_  w34(nqq), w35(nqq), w36(nqq)
   _REAL_  sd2(6), work(1000)
   _REAL_  t(3,3)
   _REAL_  tempx(3),tempy(3)
   integer itemp(nqq)

   ! ----- set SVD subroutine: isvd=1  call ssvdc()
   !                           isvd=2  call sgecvd()

   isvd = 1
   hmax=h

   if (nq /= nqq) then
      write(*,*) "   Please change nqq in coedef() s.t. nqq=nq in main()!"
      stop
   end if
   alf = 18.1*hmax
   if0 = int(alf/hmax) + 1

   do i=1, nirreg
      do j=1, nq
         wcoe(i,j)   = 0.0
         wycoe(i,j)  = 0.0
         wzcoe(i,j)  = 0.0
         wyycoe(i,j) = 0.0
         wzzcoe(i,j) = 0.0
         wyzcoe(i,j) = 0.0
      end do
   end do

   ! ----- for each projection point

   do ir = 1, nirreg

      x1     = cirreg( 1, ir)
      y1     = cirreg( 2, ir)
      z1     = cirreg( 3, ir)

      xyy    = cirreg( 4, ir)
      xzz    = cirreg( 5, ir)
      xyz    = cirreg( 6, ir)

      t(1,1) = cirreg( 7, ir)
      t(1,2) = cirreg( 8, ir)
      t(1,3) = cirreg( 9, ir)
      t(2,1) = cirreg(10, ir)
      t(2,2) = cirreg(11, ir)
      t(2,3) = cirreg(12, ir)
      t(3,1) = cirreg(13, ir)
      t(3,2) = cirreg(14, ir)
      t(3,3) = cirreg(15, ir)

      i0 = nint((x1-xs)/hx)
      j0 = nint((y1-ys)/hy)
      k0 = nint((z1-zs)/hz)

      ! -------- initialize

      do jj=1, nq
         do ii=1, nq
            v(ii, jj) = 0.0      !# the V matrix in SVD: U*S*V
         end do
         do ii=1, 6
            w1(ii, jj) = 0.0     !# the coefficient matrix to be SVDed
         end do
      end do

      ! -------- Generate the matrix

      nsub = 0
      do i2=0, if0          !# Starting from the center
         do i1 = i0-i2, i0+i2
         do j1 = j0-i2, j0+i2
         do k1 = k0-i2, k0+i2
            if (i1 < 1 .or. j1 < 1 .or. k1 < 1 .or. i1 > l .or. j1 > m .or. k1 > n) goto 55

            idis = abs(i1-i0) + abs(j1-j0) + abs(k1-k0)
            nn2 = index2(i1,j1,k1)
            if (idis == i2 .and. nn2 > 0) then

               x2 = cirreg(1, nn2)
               y2 = cirreg(2, nn2)
               z2 = cirreg(3, nn2)
               do ii = 1, nsub
                  ip = itemp(ii)
                  x3 = cirreg(1, ip)
                  y3 = cirreg(2, ip)
                  z3 = cirreg(3, ip)
                  call distan(x2,y2,z2,x3,y3,z3, dis)
                  if (dis < 0.1*hmax) then
                     goto 55
                  end if
               end do

               call distan(x1,y1,z1,x2,y2,z2, dis)
               if ( dis < alf ) then
                  nsub = nsub + 1
                  itemp(nsub) = nn2
                  tempx(1) = x2 - x1
                  tempx(2) = y2 - y1
                  tempx(3) = z2 - z1
                  call matvec(3,3,t, tempx, tempy)
                  w1(1, nsub)  = 1.0
                  w1(2, nsub)  = tempy(2)
                  w1(3, nsub)  = tempy(3)
                  w1(4, nsub)  = 0.5*tempy(2) * tempy(2)
                  w1(5, nsub)  = 0.5*tempy(3) * tempy(3)
                  w1(6, nsub)  = tempy(2) * tempy(3)
                  if ( nsub .eq. nq) goto 77
               end if
            end if
         55 continue
         end do  ! k1 = k0-i2, k0+i2
         end do  ! j1 = j0-i2, j0+i2
         end do  ! i1 = i0-i2, i0+i2
      end do ! i2=0, if0
   77 continue

      if (nsub < 12) then
         write(*,*) "   alf is too small in coede6f()!"
         stop
      end if

      ! -------- Call least square routine

      if (isvd == 1) then
         job = 11
         call dsvdc(w1,6,6,nsub,sd,ew,uw,6,v,nq,w4,job,inf)
      else
         call dgesvd('A','A',6,nsub,w1,6,sd2,uw,6,v,nq, &
               work, 1000, inf)
         do ij = 1, 6
            sd(ij) = sd2(ij)
         end do
      end if

      if (inf /= 0) then
         write(*,*) inf, " - ssvdc() or sgesvd() failed in coedef!"
         stop
      end if

      do i1=1, 6
         if (abs(sd(i1)) > 1.0e-14) then
            ew1(i1)=uw(1, i1)/sd(i1)
            ew2(i1)=uw(2, i1)/sd(i1)
            ew3(i1)=uw(3, i1)/sd(i1)
            ew4(i1)=uw(4, i1)/sd(i1)
            ew5(i1)=uw(5, i1)/sd(i1)
            ew6(i1)=uw(6, i1)/sd(i1)
         else
            ew1(i1)= 0.0
            ew2(i1)= 0.0
            ew3(i1)= 0.0
            ew4(i1)= 0.0
            ew5(i1)= 0.0
            ew6(i1)= 0.0
         end if
      end do

      ! -------- w31(i),w32(i),......,w36(i) are the solutions

      do i1=1, nsub
         w31(i1) = 0.0
         w32(i1) = 0.0
         w33(i1) = 0.0
         w34(i1) = 0.0
         w35(i1) = 0.0
         w36(i1) = 0.0

         if (isvd == 1) then
            do j1=1,6
               w31(i1) = w31(i1) + v(i1,j1)*ew1(j1)
               w32(i1) = w32(i1) + v(i1,j1)*ew2(j1)
               w33(i1) = w33(i1) + v(i1,j1)*ew3(j1)
               w34(i1) = w34(i1) + v(i1,j1)*ew4(j1)
               w35(i1) = w35(i1) + v(i1,j1)*ew5(j1)
               w36(i1) = w36(i1) + v(i1,j1)*ew6(j1)
            end do
         else
            do j1=1,6
               w31(i1) = w31(i1) + v(j1,i1)*ew1(j1)
               w32(i1) = w32(i1) + v(j1,i1)*ew2(j1)
               w33(i1) = w33(i1) + v(j1,i1)*ew3(j1)
               w34(i1) = w34(i1) + v(j1,i1)*ew4(j1)
               w35(i1) = w35(i1) + v(j1,i1)*ew5(j1)
               w36(i1) = w36(i1) + v(j1,i1)*ew6(j1)
            end do
         end if

         wcoe(ir,i1)   = w31(i1)
         wycoe(ir,i1)  = w32(i1)
         wzcoe(ir,i1)  = w33(i1)
         wyycoe(ir,i1) = w34(i1)
         wzzcoe(ir,i1) = w35(i1)
         wyzcoe(ir,i1) = w36(i1)
      end do  !  i1=1, nsub
   end do  ! ir = 1, nirreg

end subroutine coed6

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The IIM method solves the singular free PB equation.
! bcopt = 4 or 5: inside is total field and outside is total field.
! bcopt = 6 or 7: inside is reaction field and outside is total field.
! bcopt = 8 or 9: inside is reaction field and outside is reaction field.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! the jump condition of potential at projection point (x,y,z)
! bcopt = 4 or 5: it is zero.
! bcopt = 6 or 7: it is the coulomb potential.
! bcopt = 8 or 9: it is zero.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ the jump condition of potential at projection point (x,y,z)
function fw(bcopt,bi,bo,atmfirst,atmlast,x,y,z)

   use poisson_boltzmann, only : acrg,acrd
   implicit none

   _REAL_ :: bi,bo

   integer bcopt,atmfirst,atmlast
   _REAL_ fw
   _REAL_ x,y,z

   integer k
   _REAL_ r(1:3),d
   _REAL_ phi0

   fw = 0.d0
   if ( bcopt == 6 .or. bcopt == 7 ) then
      do k = atmfirst, atmlast
         r(1) = x - acrd(1,k)
         r(2) = y - acrd(2,k)
         r(3) = z - acrd(3,k)
         d = sqrt(sum(r*r))
         fw = fw + acrg(k) / d
      end do
   end if

end function fw
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! the jump condition of field (b*E_n) at projection point (x,y,z)
! bcopt = 4 or 5: it is zero.
! bcopt = 6 or 7: it is  bi       * coulomb field.
! bcopt = 8 or 9: it is (bi - bo) * coulomb field.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ the jump condition of field (b*E_n) at projection point (x,y,z)
function fq(bcopt,bi,bo,atmfirst,atmlast,x,y,z,t)

   use poisson_boltzmann, only : acrg,acrd
   implicit none

   _REAL_ :: bi,bo

   integer bcopt, atmfirst, atmlast
   _REAL_ fq
   _REAL_ x,y,z,t(3,3)

   integer k
   _REAL_ r(1:3),rn(1:3),d,d2,f(3),ft(3)

   fq = 0.d0
   if ( bcopt > 5 .and. bcopt < 10 )  then
      do k = atmfirst, atmlast
         r(1) = x - acrd(1,k)
         r(2) = y - acrd(2,k)
         r(3) = z - acrd(3,k)
         d2 = sum(r*r)
         d = sqrt(d2)
         r = r / d

         f = - acrg(k) / d2 * r
         call matvec(3,3,t,f,ft)
         fq = fq + ft(1)
      end do
      if ( bcopt == 6 .or. bcopt == 7 ) fq =  bi      *fq
      if ( bcopt == 8 .or. bcopt == 9 ) fq = (bi - bo)*fq
   end if

end function fq


end module iim_util
