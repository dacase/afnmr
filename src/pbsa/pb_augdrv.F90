! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"

module iimaug

#  include "pb_constants.h"

   integer l_xm, l_ym, l_zm, l_xmym, l_xmymzm
   _REAL_ l_fmiccg,l_wsor
   _REAL_ l_norm, l_inorm, l_accept, l_epsout, l_pbkappa, l_h
   integer l_itn,l_maxitn
   integer mg_nlevel,ncyc_before,ncyc_after

   _REAL_,allocatable :: l_am1(:), l_am2(:), l_am3(:), l_ad(:)
   _REAL_,allocatable :: l_bv(:)
   _REAL_,allocatable :: l_pv(:), l_tv(:), l_zv(:), l_rd(:)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ AUG algoritm interface routine
subroutine pb_augdrv( npbstep,npbgrid,nstlim,atmfirst,atmlast,npbopt,solvopt,level,nfocus,bcopt,&
                      natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                      xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                      maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
                      pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
                      gcrd,acrg,&
                      nbnd,iepsav,insas,lvlset,&
                      chgrd,saltgrd,phi,&
                      bv,cphi,xs,&
                      pbverbose, &
                      augtoltype, augctf, augtol &
                    )
   use iim_util
   implicit none

   ! all the driver variables are shared among all "contained" routines, so are
   ! not redeclared in the containted routines anymore, except there is a need
   ! to remap their dimensions and to copy to other variables.

   ! passed variables

   logical pbverbose
   integer npbstep, npbgrid, nstlim, atmfirst, atmlast
   integer npbopt, solvopt, level, nfocus, bcopt ! level, nfocus presumably is always 1
   integer natom
   integer maxitn, itn
   integer xm, ym, zm, xmym, xmymzm
   integer savxm(nfocus), savym(nfocus), savzm(nfocus)
   integer nbnd, iepsav(4,xmymzm)
   _REAL_  fmiccg, accept, laccept, wsor, lwsor, inorm, norm
   _REAL_  pbkappa, pbkb, pbtemp, ivalence, istrng, eps0, epsin, epsout, ionene
   _REAL_  h, savh(nfocus), gox, goy, goz, savgox(nfocus), savgoy(nfocus), savgoz(nfocus)
   _REAL_  gcrd(3,natom), acrg(natom)

   integer insas(xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8)
   _REAL_  lvlset(xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8)

   _REAL_  chgrd(xmymzm), saltgrd(xmymzm), phi(xmymzm)
   _REAL_  bv(xmymzm), cphi(xmymzm), xs(xmymzm+2*xmym)
   integer ii

   integer augtoltype !WMBS 0 for absolute 1 for relative
   _REAL_  augctf !WMBS cluster radius factor <= 0 -> h*h,
                  !     > 0 -> h*augctf
   _REAL_  augtol !WMBS tolerance for gmres routine

   ! local varialbes

   _REAL_ rh

   ! for cluster cutoff

   rh = ONE/h

   ! for singularity-free PB equation we need to compute cphi() at the
   ! dielectric interface grid points

   if ( bcopt > 5 .and. bcopt < 10 ) then
      cphi(1:xmymzm) = ZERO
      call pb_dbcgrd( cphi(1), insas )
   end if

   ! for singular PB equation bv() is initialized as the grid charges

   if ( bcopt < 6 .or. bcopt > 9 ) then
      bv = chgrd
   end if

   ! now we put effective charges at the space boundary grid points into bv()
   ! to take care of the boundary conditions

   call pb_bndcnd( bv(1), chgrd(1) )
   call aug(gox,goy,goz,xm,ym,zm,lvlset,insas,nbnd,iepsav, &
            epsin/eps0,epsout/eps0, &
            bv(1),phi(1),accept,h,atmfirst,atmlast,bcopt,solvopt,&
            augtoltype,augctf,augtol)

contains

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
   if (.not. ( (bcopt == 10) .and. (solvopt == 7 &
        .or. solvopt == 3))) then
           bv=0.0d0 !XP: This can make imin=6 work with ipb4 5
        !WMBS: This will cause problems for fft or pcg solver
   end if

   if ( level == 1 .and. bcopt == 10 ) then
      if (.not. (solvopt == 7 .or. solvopt == 3)) then
        write(6, '(a)') "PB bomb in pb_bndcnd(): zero BC only allowed for fft or pcg"
        call mexit(6, 1)
      end if

   ! bcopt = 2
   ! molecule dipolar debye-huckel potential in the singular PB.
   ! the boundary will be all solvent.

   else if ( level == 1 .and. bcopt == 2 ) then
      write(6, '(a)') "PB bomb in pb_iimdrv(): molecular dipolar BC not supported"
      call mexit(6, 1)

   ! bcopt = 3
   ! sum of residue dipolar debye-huckel potentials in the singular PB.
   ! the boundary will be all solvent.

   else if ( level == 1 .and. bcopt == 3 ) then
      write(6, '(a)') "PB bomb in pb_iimdrv(): residue dipolar BC not supported"
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

      write(6, '(a,i6)') 'PB bomb in pb_iimdrv(): unknown BC option', bcopt
      call mexit(6, 1)
   end if


end subroutine pb_bndcnd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ This is the AUG driver
subroutine aug(gox,goy,goz,l,m,n,lvlset,insas,nbnd,iepsav,epsin,epsout,bv,u,accept,h,&
               atmfirst,atmlast,bcopt,solvopt,augtoltype,augctf,augtol)
   ! the AUG method used to solve the PBE with different delectric constants of
   ! two regions separated by the interface

   !    xs <- gox, ys <- goy, zs <- goz
   !    l  <- xm , m  <- ym , n  <- zm
   !    bi <- epsin, bo <- epsout
   !    u  <- phi,u0 <- xs

   use density_surface, only : index, index2, x, y, z, cirreg
   use iim_util
   implicit none

   ! passed variables

   !  go[x-z]                starting crds of the box r
   !  l,m,n                  dimension of the 3-D box
   !  lvlset                 level set function
   !  insas                  index of grids to differentiate inner and outer
   !                         sider, the defination varies to the sasopt
   !  nbnd                   number of irregular grid points
   !  iepsav                 array storing the grid crds of irregular ptns, it
   !                         originally also store the number of atom the grid
   !                         point belongs to, yet it does not apply the AUG
   !  epsin,epsout           dielectric constants of inside and outside
   !  bv                     initial rhs of the linear eqns of the whole grids
   !  u                      potential of whole grids
   !  accept                 convergence criteria of the linear solver
   !  h                      grid spacing
   !  atmfirst,atmlast       the first and last atm index

   _REAL_  gox,goy,goz,h,accept,epsout,epsin,tol,dist
   integer l,m,n,nbnd,atmfirst,atmlast,bcopt,solvopt,iter
   integer imax,mm
   integer insas(0:l+1,0:m+1,0:n+1),iepsav(1:4,1:l*m*n)
   _REAL_  bv(l,m,n),u(l,m,n)
   _REAL_  lvlset(0:l+1,0:m+1,0:n+1)

   integer augtoltype ! 0 absolute, 1 relative
   _REAL_  augctf !cluster radius factor
   _REAL_  augtol !gmres tolerance

   _REAL_ xs,ys,zs,xf,yf,zf,hx,hy,hz,hmax,bi,bo,error
   integer  nirreg,nind3

   ! local variables
   ! x(),y(),z()             coordinates of grid points
   ! cirreg cirreg1          geometrical parameters at irregular points
   ! index()                 flag of all grid points with 1-5
   ! index2()                number index for irregular grid points
   ! index3()                number index for inner side irregular points
   ! wp,qp,wp1,qp1           [u]and[betaUn]for irregular and inner irregular points
   ! unj(),unjf()            [Un] for irregular points
   ! bf(), fvec()            rhs of Ag=b,output of matrix-vector multiply
   ! w0p,wyp,wzp,w0p,        tangential derivatives of jump conditions
   ! wyp,wzp
   ! wcoe,wxcoe,wycoe,wzcoe  coefficients obtained by SVD
   ! wxxcoe,wyycoe wzzcoe
   ! wxycoe,wxzcoe,wyzcoe
   ! c,c2                    correction terms
   ! sss1,sss2               dummy arrays

   _REAL_  ctf
   integer ctn,icn
   integer it,jt,kt,iit,it1
   integer i1,j1,k1,flag,i2,j2,k2
   _REAL_, parameter :: eps0 = 8.8542D-12 / (1.6022D-19)**2 /  &
                               (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   integer, parameter :: nq=27
   integer i,j,k,IFAIL,IFAIL2,ii,ir,nc,nz_num
   _REAL_ rhs
   _REAL_, allocatable :: cirreg1(:,:),ct(:,:)! geometrical parameters at irregular points
   _REAL_,allocatable :: unj(:),unjf(:),bf(:),fvec(:)
   _REAL_,allocatable :: q0p(:)
   _REAL_,allocatable :: w0p(:)
   _REAL_,allocatable :: wcoe(:,:),wxcoe(:, :),wycoe(:, :)
   _REAL_,allocatable :: wzcoe(:,:),wxxcoe(:, :),wyycoe(:, :)
   _REAL_,allocatable :: wzzcoe(:, :),wxycoe(:, :)
   _REAL_,allocatable :: wxzcoe(:, :),wyzcoe(:, :)
   _REAL_,allocatable :: sss1(:),sss2(:)
   _REAL_,allocatable :: c(:,:,:,:), c2(:, :)
   _REAL_,allocatable :: f(:,:,:)
   character(10) str
   integer icall

   character*11 potentialdxfilename
   character*9  potentialdxdataname
   integer      potentialdxfilenum

   logical use_average_tolerance

   if (augtoltype == 1) then
      use_average_tolerance = .true.
   else
      use_average_tolerance = .false.
   endif

   potentialdxfilename= "pbsa_phi.dx"
   potentialdxdataname= "Potential"
   potentialdxfilenum= 90885

   ! interface variables

   nirreg = nbnd
   xs = gox; ys = goy; zs = goz
   hx = h; hy = h; hz = h; hmax = h

   bi = epsin; bo = epsout

   ! allocating working arrays

   allocate(cirreg1(1:15, 1:nbnd),ct(1:3, 1:nbnd))
   allocate(wp(nbnd),qp(nbnd))
   allocate(unj(nbnd),unjf(nbnd),bf(nbnd),fvec(nbnd))
   allocate(q0p(nbnd),qyp(nbnd),qzp(nbnd))
   allocate(w0p(nbnd),wyp(nbnd),wzp(nbnd))
   allocate(wyyp(nbnd),wzzp(nbnd),wyzp(nbnd))
   allocate(wcoe(nbnd,nq),wxcoe(nbnd, nq),wycoe(nbnd, nq))
   allocate(wzcoe(nbnd,nq),wxxcoe(nbnd, nq),wyycoe(nbnd, nq))
   allocate(wzzcoe(nbnd, nq),wxycoe(nbnd, nq))
   allocate(wxzcoe(nbnd, nq),wyzcoe(nbnd, nq))
   allocate(sss1(nbnd),sss2(nbnd))
   allocate(c(l,m,n,7), c2(nbnd, 27))
   allocate(f(l,m,n))

   ! XP: cluster irregular grid points to get the new index2 and cirreg

   if ( augctf > 0 .and. augctf <= 1) then
      ctf = h*augctf
   else
      ctf=h*h             ! cut off
   endif

   ctn=0;index2=0         ! initialization
   do ii=1,nbnd
      if( ii==1) then       !1st one
         ct(1:3,1)=cirreg(1:3,ii)
         ctn=ctn+1
         cirreg1(1:15,ctn)=cirreg(1:15,ii)
         i1=iepsav(1,ii);j1=iepsav(2,ii);k1=iepsav(3,ii)
         index2(i1,j1,k1)=ctn
      else
         do icn=1,ctn         !not the 1st one, compare to the nearest current point
            dist=sqrt((ct(1,icn)-cirreg(1,ii))**2+(ct(2,icn)-cirreg(2,ii))**2+(ct(3,icn)-cirreg(3,ii))**2)
            if(dist <= ctf) then ! fall into a current one
               i1=iepsav(1,ii);j1=iepsav(2,ii);k1=iepsav(3,ii)
               index2(i1,j1,k1)=icn
               exit
            end if
            if(icn==ctn) then    ! not fall into a current one
               ctn=ctn+1
               ct(1:3,ctn)=cirreg(1:3,ii)
               cirreg1(1:15,ctn)=cirreg(1:15,ii)
               i1=iepsav(1,ii);j1=iepsav(2,ii);k1=iepsav(3,ii)
               index2(i1,j1,k1)=ctn
            end if
         end do
      end if
   end do

   nirreg=ctn ! ctn the after-cluster dimension

   ! setting up mm for gmres by experience

   if ( ctn < 1000 ) then
      mm= ceiling(real(ctn)/16.0d0)
   else
      mm= ceiling(real(ctn)/60.0d0)
   end if

   ! setting up the jump conditions at the projection points of the irreular (boundary) grid points

   wp=0.0d0;qp=0.0d0  ! Initialization

   call prodis(l,m,n,nirreg,bcopt,bi,bo,atmfirst,atmlast,nbnd,cirreg1)!,wp[u] qp [beta_un]

   ! prepare to call matvec the first time to compute F2
   ! set g=0(unjf)  and also b(bf)=0 to get F2 as in IIM_augment

   do i=1, ctn !
      unj(i)=0.0d0
      unjf(i)=0.0d0
      bf(i)=0.0d0
   end do

   icall = 1

   ! 1st time to call matvec to get the bf as in the book ch 6.1
   ! To solve (T-EA^-1 B)G=F2,  matvec(0)=(T-EA^-1 B)0-F2=-F2

   call matvec3(l,m,n,h,nbnd,ctn,x,y,z,xs,ys,zs,lvlset,u,f,&
           cirreg1,&
           unj,unjf,fvec,bf,epsin,epsout,index,index2,&
           q0p,w0p,&
           wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
           wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
           accept,bv,icall,solvopt,norm,inorm)

   do i=1,  ctn
      bf(i)= -fvec(i)
      unjf(i)=qp(i)
   end do

   imax= 3*max(l,m,n)! ;stop                             !set values for imax

   if (augtol > 0 .and. augtol < 1) then
      tol = augtol
   else
      tol = 1.0d-5
   endif

   if (use_average_tolerance) then
      tol = tol*ctn
   end if

   unj=0.0d0; unjf=0.0d0;u=0.0d0;f=0.0d0  !Initialization..

   call gmresx(mm,l,m,n,nbnd,ctn,imax,h,lvlset,x,y,z,xs,ys,zs,&
           index,index2,&
           cirreg1,&
           unj,u,f,unjf,bf,&
           tol,iter,error,epsout,epsin,&
           q0p,w0p,&
           wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
           wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
           accept,bv,solvopt,pbverbose,norm,inorm)

   ! converting back to the PBSA unit for potential

   do k=1, n
   do j=1, m
   do i=1, l
      u(i,j,k) = u(i,j,k) * INV_FOURPI / eps0
   end do
   end do
   end do

   deallocate (cirreg1,ct)
   deallocate (wp,qp)
   deallocate (unj,unjf,bf,fvec)
   deallocate (q0p,qyp,qzp)
   deallocate (w0p,wyp,wzp)
   deallocate (wyyp,wzzp,wyzp)
   deallocate (wcoe,wxcoe,wycoe)
   deallocate (wzcoe,wxxcoe,wyycoe)
   deallocate (wzzcoe,wxycoe)
   deallocate (wxzcoe,wyzcoe)
   deallocate (sss1,sss2)
   deallocate (c, c2)
   deallocate (f)

end subroutine aug
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

end subroutine pb_augdrv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Caculate A*X with the Un inputted, which is [beta*un] + bf  ]
subroutine matvec3(xm,ym,zm,h,nbnd,nind3,x,y,z,xs,ys,zs,phi,u,f,&
             cirreg,&
             unj,unjf,fvec,bf,epsin,epsout,index,index2,&
             q0p,w0p,&
             wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe,&
             wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2,&
             accept,bv,icall,solvopt,norm,inorm )

   use iim_util

   ! input unj, output fvec
   ! nind3 is the  reduced number of  all the irregular points
   ! note unjf and  bf are passed in u&u0 are not
   ! note this is two sides aug approaches.

   implicit none
   integer, parameter :: nq=27
   logical alive
   integer :: status = 0

   ! passed variables

   ! nbnd,nind3              the initial dimension; reduced dimension of aug
   ! x(),y(),z()             coordinates of grid points
   ! xs,ys,zs                staring crds of the box
   ! phi,u                   lvlset function; potential of the whole grids
   ! cirreg                  geometrical parameters at irregular points
   ! index()                 flag of all grid points with 1-5
   ! index2()                number index for irregular grid points
   ! wp,qp                   [u]and[betaUn]for the interface
   ! unj(),unjf()            [Un] for irregular points, unjf is the input
   ! bf(), fvec()            rhs of Ag=b; output of matrix-vector multiply
   ! w0p,wyp,wzp,w0p,        tangential derivatives of jump conditions
   ! wyp,wzp
   ! wcoe,wxcoe,wycoe,wzcoe  coefficients obtained by SVD
   ! wxxcoe,wyycoe wzzcoe
   ! wxycoe,wxzcoe,wyzcoe
   ! bv                      initial rhs of the linear eqns of the whole grids
   ! solvopt,icall           solvopt for linear solvers above, debuging variable
   ! c,c2                    linear coefficients
   ! accept                  convergence criteria of the linear solvers
   ! sss1,sss2               dummy arrays

   ! local variables

   ! xp,yp,zp                projection crds on interface of the irregular pnts
   ! eps[x-z]                the coefs of grid pnts, which is 1 for AUG method
   ! f(i,j,k)                rhs of the eqns for the whole grid for AUG method
   ! xyy,xzz,xyz             derivatives of the geometry on projection pnts
   ! t                       trans matrix between local crds and grid crds
   ! unin                    derivative of u from inner side on proj pnts

   integer xm,ym,zm,nbnd,nind3,n3,n1,icall
   integer index(1:xm,1:ym,1:zm),index2(1:xm,1:ym,1:zm)
   _REAL_ xyz,xyy,xzz,gox,goy,goz,h,accept,epsout,epsin
   _REAL_  q0p(nbnd),w0p(nbnd)

   _REAL_  wcoe(nbnd,nq),wxcoe(nbnd, nq),wycoe(nbnd,nq),wzcoe(nbnd,nq),&
           wxxcoe(nbnd, nq),wyycoe(nbnd, nq),wzzcoe(nbnd, nq),&
           wxycoe(nbnd, nq),wxzcoe(nbnd, nq),wyzcoe(nbnd,nq)
   _REAL_  sss1(nbnd),sss2(nbnd)
   _REAL_  c(xm,ym,zm,7),c2(nbnd, 27)
   _REAL_  t(3,3)
   _REAL_  bv(xm,ym,zm)
   _REAL_  x(0:xm+1),y(0:ym+1),z(0:zm+1)
   _REAL_  phi(0:xm+1,0:ym+1,0:zm+1)
   _REAL_  u(1:xm,1:ym,1:zm),f(xm,ym,zm)
   _REAL_  unj(1:nbnd),cirreg(1:15,1:nbnd),&
           unjf(1:nbnd),fvec(1:nbnd),bf(1:nbnd)

   ! variables to use mg

   _REAL_  inorm,norm,dummy
   integer itn , solvopt
   _REAL_  epsx,epsy,epsz
   _REAL_  iv(1:xm*ym*zm)
   _REAL_  xso(xm*ym*zm+2*xm*ym)

   _REAL_ xs,ys,zs,xf,yf,zf,hx,hy,hz,hmax,bi,bo
   _REAL_ uu,dudx,dudy,dudz,unin,unout
   integer l, m, n, nirreg

   !sgni,j,k is not used for now... XL

   integer :: sgni=0
   integer :: sgnj=0
   integer :: sgnk=0
   integer i0,j0,k0,nn1

   !WMBS - for use in neutralizing charge on f for periodic solvers

   integer cntirreg !count of the number of irregular grid nodes
   _REAL_ fsum !holds the average charge per irregular grid node to
               !be neutralized

   integer i,j,k,IFAIL,IFAIL2,ii,ir,nc,nz_num
   _REAL_ beta_max,rhs
   _REAL_ coe2(27)

   _REAL_ xp,yp,zp

   _REAL_ qpsav(nbnd)

   l = xm; m = ym; n = zm; nirreg = nind3
   hx=h;hy=h;hz=h

   epsx=1.0d0;epsy=1.0d0;epsz=1.0d0
   bi=epsin;bo=epsout

   cntirreg = 0

   do i = 1, nind3
      unj(i) = unjf(i)
      qpsav(i)=qp(i)
      qp(i)=unjf(i)
   end do

   ! calculate the first and second derivatives of the jump conditions in the
   ! surface tangential directions in the local coordinate system
   ! step 1: this is the first derivatives of w g for  irregular points only

   call coed20(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,nbnd,nq,index2,cirreg, &
           wcoe,wxcoe,wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)

   call qint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg, &
           wcoe,wycoe,wzcoe,q0p)

   call wint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,wp, &
           wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,wyp,wzp,wyyp,wzzp,wyzp)


   ! step 2: this is the second derivatives ( when secodary order, should be used )

   call coed6(l,m,n,h,hx,hy,hz,xs,ys,zs,nirreg,nbnd,nq,index2,cirreg, &
           wcoe,wxcoe,wycoe,wzcoe,wxxcoe,wyycoe,wzzcoe,wxycoe,wxzcoe,wyzcoe)

   call wint(l,m,n,h,hx,hy,hz,xs,ys,zs,nbnd,nq,index,index2,cirreg,wp, &
           wcoe,wycoe,wzcoe,wyycoe,wzzcoe,wyzcoe,w0p,sss1,sss2,wyyp,wzzp,wyzp)

   ! setting up linear system coefficient matrix
   ! the current version uses irregular points on both sides

   beta_max=1.0
   nz_num = 0

   do k = 1, n
   do j = 1, m
   do i = 1, l

      ! inside regular points

      if ( index(i,j,k) == 1 ) then

         ! set up the 7-band laplassian operator in uniform beta_max

         do ii = 2, 7
            c(i,j,k,ii) =1.0d0/h/h             !epsin/h/h
         end do
         if ( i == 1 ) c(i,j,k,2) = 0.d0
         if ( i == l ) c(i,j,k,3) = 0.d0
         if ( j == 1 ) c(i,j,k,4) = 0.d0
         if ( j == m ) c(i,j,k,5) = 0.d0
         if ( k == 1 ) c(i,j,k,6) = 0.d0
         if ( k == n ) c(i,j,k,7) = 0.d0
         c(i,j,k,1) = 6.0d0/h/h ! epsin/h/h  XP: this is still for general IIM, epsin /= epsout
         f(i,j,k) = bv(i,j,k)/h/h/h/epsin*FOURPI
         do ii = 1 , 7
            if ( abs(c(i,j,k,ii)) > 1.d-10 ) nz_num = nz_num + 1
         end do

      ! outside regular points

      else if (index(i,j,k) == 5 ) then

         do ii = 2, 7
            c(i,j,k,ii) = 1.0d0/h/h          !epsout/h/h
         end do
         if ( i == 1 ) c(i,j,k,2) = 0.d0
         if ( i == l ) c(i,j,k,3) = 0.d0
         if ( j == 1 ) c(i,j,k,4) = 0.d0
         if ( j == m ) c(i,j,k,5) = 0.d0
         if ( k == 1 ) c(i,j,k,6) = 0.d0
         if ( k == n ) c(i,j,k,7) = 0.d0
         c(i,j,k,1) = 6.0d0/h/h ! epsout/h/h ! XP: this is still for general IIM, epsin /= epsout
         f(i,j,k) = bv(i,j,k)/h/h/h/epsout*FOURPI
         do ii = 1, 7
            if ( c(i,j,k,ii ) > 1.d-10 ) nz_num=nz_num+1
         end do

      ! irregular points

      else
         cntirreg = cntirreg+1  !WS
         coe2=0.0d0

         ! XP: Note the interface has changed. b_in, b_out, and bi bo are all passed in.

         bo=bi
         call irre31(l,m,n,h,hx,hy,hz,IFAIL, &
                 i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,phi,index, &
                 nq,nbnd,index2,cirreg,coe2,rhs)

         if (IFAIL.gt.10) call irre32(l,m,n,h,hx,hy,hz,IFAIL2, &
                 i,j,k,index(i,j,k),beta_max,bi,bo,x,y,z,phi,index, &
                 q0p,w0p, &
                 nq,nbnd,index2,cirreg,coe2,rhs)

         ir = index2(i,j,k)
         bo=epsout
         f(i,j,k) = -rhs            !*(-6.0d0)/coe2(14)

         ! Out of the 27 neighbors, 7 is nonzero and contribute

         do nc = 1, 27
            c2(ir,nc) = coe2(nc)
         end do
         do ii = 1, 27
            if ( abs(c2(ir,ii)) > 1.d-10 ) nz_num = nz_num + 1
         end do
      end if
   end do
   end do
   end do

   f=f*h*h
   dummy=0.0d0;iv=0.0d0; !use iccg
   call init_param(l,m,n,l*m,l*m*n,10000,dummy,accept, 0.0d0,1.0d0,h,1.9d0)
   call allocate_array(solvopt)
   xso(:) = 0.d0
   call init_array(solvopt,epsx,epsy,epsz,f,iv,xso)

   if (solvopt /= 7 ) then !.and. solvopt/=5) then
      call pb_iccg(u,xso)
   end if

   itn = l_itn
   inorm = l_inorm
   norm = l_norm

   call deallocate_array(solvopt)

   n3=10! second order interplation
   do k=1,n
   do j=1,m
   do i=1,l
      nn1=index2(i,j,k)

      if (nn1 .gt. 0) then

         ! projection crds and geometrical info

         xp     = cirreg( 1, nn1)
         yp     = cirreg( 2, nn1)
         zp     = cirreg( 3, nn1)

         xyy    = cirreg( 4, nn1)
         xzz    = cirreg( 5, nn1)
         xyz    = cirreg( 6, nn1)

         t(1,1) = cirreg( 7, nn1)
         t(1,2) = cirreg( 8, nn1)
         t(1,3) = cirreg( 9, nn1)
         t(2,1) = cirreg(10, nn1)
         t(2,2) = cirreg(11, nn1)
         t(2,3) = cirreg(12, nn1)
         t(3,1) = cirreg(13, nn1)
         t(3,2) = cirreg(14, nn1)
         t(3,3) = cirreg(15, nn1)

         !the nearest grid point

         i0 = nint((xp - x(0))/h)
         j0 = nint((yp - y(0))/h)
         k0 = nint((zp - z(0))/h)

         bo=bi
         call twosided(l,m,n,bi,bo,n3,i0,j0,k0,&
                 nbnd,xp,yp,zp,uu,dudx,dudy,dudz,&
                 xyy,xzz,xyz,t,u,phi,x(0),y(0),z(0),h,&
                 nn1)
         bo=epsout

         ! call interp output Un^-, dudx in the local crds

         unin = dudx
         unout = unin + unj(nn1)

         ! Preconditioning, which results in smaller residual
         ! unout = (qp(nn1)-unj(nn1))/79.0d0 , as we are using 1:80

         fvec(nn1)=epsout*unout-epsin*unin-qpsav(nn1)
         fvec(nn1)=fvec(nn1)+bf(nn1)
      end if
   end do
   end do
   end do

   do i = 1, nind3
      qp(i) = qpsav(i)
   end do

end subroutine matvec3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ the core iterative driver for augmented solver
subroutine gmresx(mm,l,m,n,nbnd,ctn,imax,h,phi,x,y,z,xs,ys,zs,&
              index,index2,&
              cirreg,&
              unj,u,f,x0,bf,&
              tol,iter,error,bout,bin,&
              q0p,w0p,&
              wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
              wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
              accept,bv,solvopt,pbverbose,norm,inorm)

   ! The following process is the GMRES(mm) iterative method for  the algebraic
   ! system  A*x = b  with A being any non-symmetric matrix. The original code
   ! is under FIIM LEVEL from Zhilin and is transformed form 2D to 3D now. nbnd
   ! is number of all irregular points,while ctn is the reduced one after clustering now.
   !
   ! The algorithm for GMRES is described as below:
   ! To solve A*X=f; x0 is the guess vector for solution r0=f-Ax0 and v1=r0/||r0||
   ! The key aim is to use Schmidt orthogonalization method to orthogonalize the
   ! Krylov subspace {v1,Av1,A^2V1,...,A^(k-1)v1}
   ! Iterate: for j=1,2,...,m  ;
   ! Vj+1=Avj-sum((Avj*vi)*vi,i=1,2,..,j, vi is normalized.
   ! vj+1=Vj+1/||Vj+1||  Got the nomalized and orthogonalized basis.
   ! solution x=x0+Vm*ym, please refer to related material to know the math of
   ! implicit residual.
   ! And A*x is done my subroutine matvec.
   !
   ! Implicit residual, methods to calculat the residual without needing to
   ! calculate ym at the time:
   ! AVk=Vk+1 H'  res=||f-A(x0+vmym||=||r0-Vk+1 H'ym||
   ! => res=||Vk+1(beta*e1-H'ym)||=||beta*e1-H'ym||
   ! Please be advised that only the last element of the last row of H' is nonzero,
   ! which makes that when Q applys, res=||Q(beta*e1-H'ym)||=||gk-R*ym|| last row
   ! of R is zero, and the minimum norm would be the others are zero, thus the last
   !  element of gk is the residual, the implicit one.

   use iim_util
   implicit none
   integer, parameter :: nq=27

   !  Passed variables
   !  pbverbose              logical flag variable, .true. to print detail
   !  mm	             INTEGER number of iterations between restarts
   !  l,m,n                  dimension of the 3-D box
   !  [x-z]s,x-z             starting crds and 3-D crds array of the box
   !  nbnd,ctn               original dimension of irregular ptn; clustered one
   !  imax                   max times of restarts of the GMRES algorithm
   !  h                      grid spacing
   !  phi                    level set function
   !  index                  the index labeling for all grid points
   !  index2                 index of number for irregular points
   !  cirreg                 geometrical parameters at irregular points
   !  unj                    jump of 1st derivative of potential on interface
   !  wp,qp                  [u]and[betaUn]for the interface
   !  x0	             REAL initial guess vector
   !  bf f                   REAL right hand side vector & variable to store
   !  tol                    determine resolution of the solution
   !  bout,bin               dielectric coefficients of outside and inside
   !  w0p,wyp,wzp,w0p,       tangential derivatives of jump conditions
   !  wyp,wzp
   !  wcoe,wxcoe,wycoe,wzcoe coefficients obtained by SVD
   !  wxxcoe,wyycoe wzzcoe
   !  wxycoe,wxzcoe,wyzcoe
   !  accept                 convergence criteria of the linear solvers
   !  bv                     initial rhs of the linear eqns of the whole grids
   !  solvopt                solvopt for linear solver of whole grid potential
   !
   ! Output:
   !  u                      potential of whole grids
   !  x0                     The solution vector
   !  iter          	     The number of iterations
   !  error	             The 2-norm of the residue of the solution

   logical pbverbose
   _REAL_  accept,tol,bin,bout,h, ssum,nsum,xs,ys,zs
   integer mm,m,n,l,n2,nbnd,imax,solvopt,ctn
   _REAL_  q0p(nbnd),w0p(nbnd)

   _REAL_  wcoe(nbnd,nq),wxcoe(nbnd, nq),wycoe(nbnd,nq),wzcoe(nbnd,nq),&
           wxxcoe(nbnd, nq),wyycoe(nbnd, nq),wzzcoe(nbnd, nq),&
           wxycoe(nbnd, nq),wxzcoe(nbnd, nq),wyzcoe(nbnd,nq)
   _REAL_  sss1(nbnd),sss2(nbnd)
   _REAL_  c(l,m,n,7),c2(nbnd, 27)
   _REAL_  x0(nbnd),bf(nbnd),r(nbnd),s(nbnd), &
           x11(nbnd),vk(nbnd),w(nbnd),hj(mm,2),yy(nbnd),&
           x110(nbnd),w0(nbnd),yy0(nbnd)
   _REAL_  x(0:l+1),y(0:m+1),z(0:n+1)
   _REAL_  f(l,m,n),u(1:l,1:m,1:n),phi(0:l+1,0:m+1,0:n+1),bv(l,m,n)
   _REAL_  unj(nbnd)
   _REAL_  cirreg(1:15,1:nbnd)
   integer index(1:l,1:m,1:n),index2(1:l,1:m,1:n)
   _REAL_  norm, inorm

   ! Local variables
   !
   !  hg                     stores the Hessenberg matrix.
   !  v                      stores the modified Gram-Schmidt orthonormal
   !		             which forms the Krylov subspace.
   !  ser1,ser2              implicit residual; explicit residual
   !  RSD,MAXD,ABSD          all kinds of criteria for convergence,debuging use
   !
   ! n2<-nind3  index<-index2    index2<-index3

   integer idimf,infon,i1,j2,i2,k,k1,k2,mi,n1,j,i,iter,j1,icall
   _REAL_  drl,ser0,a,b,d,dr1,dsqrt,&
           hit,dt,st,ser1,serr,error,ser2,RSD0,MAXD0,ABSD0,RSD,MAXD,ABSD
   _REAL_,allocatable :: hg(:,:),v(:,:)

   allocate (hg(nbnd,mm),v(nbnd,mm))
   write(6,'(a,i8)') 'The gmres dimension is ',ctn
   iter = 0

   ! Loop 2, this is the restart of the algorthm with a guess vector x0

   do j=1,imax

      ! Initializaiton for variables used each start of the algorithm

      ser0 = 0.0;hg=0.0d0;v=0.0d0;w0=0.0d0
      icall = 20  ! Flag for debuging

      ! first time to call matvec to perform r0=Ax0-b, get residual r0 and calculate v1

      call matvec3(l,m,n,h,nbnd,ctn,x,y,z,xs,ys,zs,phi,u,f,&
              cirreg,&
              unj,x0,yy0,bf,bin,bout,index,index2,&
              q0p,w0p,&
              wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
              wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
              accept,bv,icall,solvopt,norm,inorm )

      do i=1,ctn
         yy(i)=yy0(i)
      end do

      call residX(ctn,yy,bf,r)       !chng resid to residX  note we use 2 sides

      ! Initial values of RSD MAXD and ABSD values, used for comparative convergence
      ! criteria

      dr1 = dsqrt(dot_product(r(1:ctn),r(1:ctn)))

      ! if x0 satisfy the eqn..  already get the desired potential u, return

      if(dr1 <= 1e-12) return

      ! Get v1
      do i=1,ctn
         v(i,1) = r(i)/dr1
         s(i) = 0.0d0
      end do

      s(1) = dr1

      !--Loop 1-- For i=1,2, ..., mm -----calculate w v and matrix h-------------------------

      do i=1,mm-1
         iter = iter + 1

         do i1=1,ctn
            vk(i1) = v(i1,i)
         end do

         ! do Avi, calculate the next basis vector of the Krylov subspace

         call matvec3(l,m,n,h,nbnd,ctn,x,y,z,xs,ys,zs,phi,u,f,&
                 cirreg,&
                 unj,vk,w0,bf,bin,bout,index,index2,&
                 q0p,w0p,&
                 wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
                 wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
                 accept,bv,icall,solvopt,norm,inorm)             !input vk,output w

         do k=1,ctn
            w(k)=w0(k)
         end do

         ! Calculate the projs of the new  basis in direction of former basises

         do k=1,i

            do i1=1,ctn
               vk(i1) = v(i1,k)
            end do

            hg(k,i) = dot_product(w(1:ctn),vk(1:ctn)) ! basis in vk direction

            ! orthogonalizing..

            do i1=1,ctn
               w(i1) = w(i1) - hg(k,i)*vk(i1)
            end do

         end do

         ! norm of the next orthogonalized basis, the i+1 th one

         hg(i+1,i) = dsqrt(dot_product(w(1:ctn),w(1:ctn)))

         ! If hg == 0, then we got the complete basis of the space.

         if( abs(hg(i+1,i)) <= 1e-14) then
            infon = 1
            write(6,'(a)') 'info=1, Got complete basis of space.'
         else
            infon = 0

            ! the nomalized new vector of the i+1 th.

            do i1=1,ctn
               v(i1,i+1) = w(i1)/hg(i+1,i)
            end do
         end if

         !----------------------------------------------------------------
         ! Convert matrix h to get the residual implicitly
         !
         !----- Apply J_1, j_2, ..., J_{i-1} on (hg_{1,i}, ..., hg_{i+1,i} ---------
         !   Suppose J_i =
         !                  | I 0 |    p = | cos(\alf)  sin(\alf) |
         !                  | 0 P |            | -sin(\alf) cos(\alf) |
         !
         !       cos(\alf) = hg(i,i)/sqrt(hg(i,i)^2 + hg(i+1,i)^2)
         !       sin(\alf) = hg(i+1,i)/sqrt(hg(i,i)^2 + hg(i+1,i)^2)
         !
         !------- Form J_i so that the (i+1)th component of J_i*h(:,i) is zero.

         do k=1,i-1
            hit = hj(k,1)* hg(k,i) + hj(k,2)*hg(k+1,i)
            hg(k+1,i) = -hj(k,2)*hg(k,i) + hj(k,1)*hg(k+1,i)
            hg(k,i) = hit
         end do

         if(infon == 0) then
            dt = dsqrt(hg(i,i)*hg(i,i)+hg(i+1,i)*hg(i+1,i))
            hj(i,1) = hg(i,i)/dt
            hj(i,2) = hg(i+1,i)/dt
            st = hj(i,1)*s(i) + hj(i,2)*s(i+1)
            s(i+1) = -hj(i,2)*s(i) + hj(i,1)*s(i+1)
            s(i) = st
            hit = hj(i,1)* hg(i,i) + hj(i,2)*hg(i+1,i)
            hg(i+1,i) = -hj(i,2)*hg(i,i) + hj(i,1)*hg(i+1,i) ! found by XP
            hg(i,i) = hit
         end if

         ser1 = abs(s(i+1))

         ! begine explicit test by XP

         mi=i
         yy(mi) = s(mi)/(hg(mi,mi))!+1.0d-14)

         do k=mi-1,1,-1
            yy(k) = s(k)
            do j1 = k+1,mi
               yy(k) = yy(k) - hg(k,j1)*yy(j1)
            end do
            yy(k) = yy(k)/hg(k,k)   ! The coefficients for each basises.
         end do

         do i2=1,ctn
            x11(i2) = x0(i2)

            do k=1,mi
               x11(i2) = x11(i2) + yy(k)*v(i2,k)
            end do
         end do

         do i2=1,ctn
            x110(i2)=x11(i2)
         end do

         ! calculate the real residual, explicit residual

         call matvec3(l,m,n,h,nbnd,ctn,x,y,z,xs,ys,zs,phi,u,f,&
                 cirreg,&
                 unj,x110,w,bf,bin,bout,index,index2,&
                 q0p,w0p,&
                 wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
                 wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
                 accept,bv,icall,solvopt,norm,inorm)       !intput x110 out put w

         do i2=1,ctn
            x11(i2)=x110(i2)
         end do

         call residX(ctn,w,bf,r)
         ser2 = dsqrt(dot_product(r(1:ctn),r(1:ctn)))

         ! Judge residual if OK goto calculate it explicitly if not go on loop1
         !XP: Here we can try to use the explicit residual as the criteria

         !Choices for criteria for convergence.

         !if(ser1 <= tol .or. serr<=tol ) then  Notes: serr= ser1/(last ser1),
         if (ser2 <= tol  ) then
            serr = tol - 1.0e-15
            mi = i
            goto 100
         end if

         ! restart if the explicit and implicit residuals are far different

         if(ser1/ser2 > 10.0d0 .or. ser1/ser2 < 1.0e-1) then
            mi=i
            goto 100
         end if
         ser0 = ser1
      end do
      mi = mm - 1

      !end loop 1 for i=1,2,...mm-----------------------------------

      ! Update(x_i,i) Calculate y x and call matvec to calculate the residual

 100  yy(mi) = s(mi)/(hg(mi,mi)+1.0d-14) !100    yy(mi) = s(mi)/hg(mi,mi)

      do k=mi-1,1,-1
         yy(k) = s(k)
         do j1 = k+1,mi
            yy(k) = yy(k) - hg(k,j1)*yy(j1)
         end do
         yy(k) = yy(k)/hg(k,k)
      end do

      do i=1, ctn
         x11(i) = x0(i)
         do k=1,mi
            x11(i) = x11(i) + yy(k)*v(i,k)
         end do
      end do

      do i=1,ctn
         x110(i)=x11(i)
      end do
      icall = 100
      call matvec3(l,m,n,h,nbnd,ctn,x,y,z,xs,ys,zs,phi,u,f,&
              cirreg,&
              unj,x110,w,bf,bin,bout,index,index2,&
              q0p,w0p,&
              wcoe, wxcoe,wycoe, wzcoe,wxxcoe,wyycoe,wzzcoe ,&
              wxycoe,wxzcoe,wyzcoe,sss1,sss2,c,c2 ,&
              accept,bv,icall,solvopt,norm,inorm)       !intput x110 out put w

      do i=1,ctn
         x11(i)=x110(i)
      end do

      ! explicit residual

      call residX(ctn,w,bf,r)
      s(mi+1) = dsqrt(dot_product(r(1:ctn),r(1:ctn)))

      ! update x0 for potential restart.

      do k=1,ctn
         x0(k) = x11(k)
      end do

      !-----Judge if residual OK reture, or back to loop 2 j=1 imax------------

      if( abs(s(mi+1)) < tol ) then
         deallocate(hg,v)
         return
      end if

      !avoid dead lock, when convergence can't be minimized, just restart.

      if (j > 4) then
         deallocate(hg,v)
         return
      end if

      error = s(mi+1)
   end do  ! j=1,imax
   deallocate(hg,v)

end subroutine gmresx

!*******************************************************************************

subroutine residX(n,x,b,r)

   !+ [The process is to perform y=Ax.]

   implicit none

   integer i,j,n
   _REAL_  x(n), b(n), r(n)

   do i=1,n
      r(i) = b(i) - x(i)
   end do
   return
end subroutine residX

function RSD(n,x)

   implicit none

   integer n,i
   _REAL_ x(n),RSD

   RSD=0.0d0
   do i=1,n
   RSD=RSD+x(i)*x(i)
   end do
   RSD=sqrt(RSD)
end function RSD

function MAXD(n,x)

  implicit none
  integer n,i
  _REAL_  x(n),MAXD
  MAXD=0.0d0
  do i=1,n
  if( MAXD < abs(x(i))) MAXD=x(i)
  end do
end function MAXD

function ABSD(n,x)

  implicit none
  integer n,i
  _REAL_  x(n),ABSD

  ABSD=0.0d0
  do i=1,n
  ABSD=ABSD+abs(x(i))
  end do
end function ABSD

!===========================================================================

subroutine init_param(nx,ny,nz,nxny,nxnynz,p_maxitn,p_fmiccg,p_accept,p_pbkappa,p_epsout,p_h,p_wsor)

   implicit none

   integer nx,ny,nz,nxny,nxnynz,p_maxitn
   _REAL_ p_fmiccg,p_accept,p_epsout,p_pbkappa,p_wsor,p_h

   l_xm = nx
   l_ym = ny
   l_zm = nz
   l_xmym = nxny
   l_xmymzm = nxnynz
   l_maxitn = p_maxitn
   l_accept = p_accept

!  ICCG
   l_fmiccg = p_fmiccg
!  MG
   mg_nlevel = 4
   ncyc_before = 10
   ncyc_after = 10
   l_pbkappa = p_pbkappa
   l_epsout = p_epsout
   l_h       = p_h
!  SOR
   l_wsor = p_wsor

end subroutine

!===========================================================================

subroutine allocate_array(solvopt)

   implicit none
   integer solvopt

   integer l,m,n,i

      allocate(l_ad(1:l_xmymzm+l_xmym),l_am1(1-l_xmym:l_xmymzm+l_xmym))
      allocate(l_am2(1-l_xmym:l_xmymzm+l_xmym), l_am3(1-l_xmym:l_xmymzm+l_xmym))
      allocate(l_rd(1-l_xmym:l_xmymzm),l_bv(1-l_xmym:l_xmymzm))
      allocate(l_tv(1-l_xmym:l_xmymzm+l_xmym),l_zv(1-l_xmym:l_xmymzm+l_xmym))
      allocate(l_pv(1-l_xmym:l_xmymzm+l_xmym))

end subroutine allocate_array

!===========================================================================

subroutine deallocate_array(solvopt)

   implicit none
   integer solvopt

      deallocate(l_ad,l_am1)
      deallocate(l_am2,l_am3)
      deallocate(l_rd,l_bv)
      deallocate(l_tv,l_zv)
      deallocate(l_pv)

end subroutine deallocate_array

!==============================================================================

!===========================================================================

subroutine init_array( solvopt,epsx,epsy,epsz,p_bv,p_iv,p_xs )

   implicit none

   integer solvopt
   _REAL_ epsx,epsy,epsz
   _REAL_ p_bv(1:l_xmymzm),p_iv(1:l_xmymzm)
   _REAL_ p_xs(1-l_xmym:l_xmymzm+l_xmym)

   integer lxmym,l,m,n,i
   _REAL_,allocatable :: lepsx(:), lepsy(:), lepsz(:)
   _REAL_ lfactor

      call feedepsintoam(l_xm, l_ym, l_zm, l_am1(1:l_xmymzm),  &
                         l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                                              epsx, epsy, epsz )
      l_am1(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      l_am2(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      l_am3(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO

      call feedepsintoad(l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz)

      l_am1(1-l_xmym:0) = ZERO
      l_am2(1-l_xmym:0) = ZERO
      l_am3(1-l_xmym:0) = ZERO
      call pb_setupper(l_am1(1),l_am2(1),l_am3(1))
      l_ad (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      l_bv (1-l_xmym:0) = ZERO
      l_bv (1:l_xmymzm) = p_bv(1:l_xmymzm)
      l_pv (1-l_xmym:0) = ZERO
      l_pv (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
      l_tv = ZERO
      l_zv = ZERO
      l_rd = ONE
contains

subroutine feedepsintoam(xm, ym, zm, am1, am2, am3, eps1, eps2, eps3)
    implicit none
    integer xm, ym, zm
    _REAL_ am1(1:xm,1:ym,1:zm)
    _REAL_ am2(1:xm,1:ym,1:zm)
    _REAL_ am3(1:xm,1:ym,1:zm)
    _REAL_ eps1
    _REAL_ eps2
    _REAL_ eps3
    am1(1:xm,1:ym,1:zm) = eps1
    am2(1:xm,1:ym,1:zm) = eps2
    am3(1:xm,1:ym,1:zm) = eps3
end subroutine feedepsintoam

subroutine feedepsintoad(xm, ym, zm, ad, eps1, eps2, eps3)
    implicit none
    integer xm, ym, zm
    _REAL_ ad(1:xm,1:ym,1:zm)
    _REAL_ eps1
    _REAL_ eps2
    _REAL_ eps3
    ad(1:xm,1:ym,1:zm) =                    eps1
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps1
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3
end subroutine feedepsintoad

end subroutine init_array

!===========================================================================

subroutine pb_iccg(phi,xs)

!  use poisson_boltzmann, only: level
   implicit none

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

! Local variables

   logical uconvg
   integer i, j
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

   ! initialization

   do i = 1, l_xmymzm
      l_rd(i) = ONE/( l_ad(i) - &
         l_am1(i-1   )*(       l_am1(i-1   )+l_fmiccg*l_am2(i-1   )+l_fmiccg*l_am3(i-1   ))*l_rd(i-1   ) - &
         l_am2(i-l_xm  )*(l_fmiccg*l_am1(i-l_xm  )+       l_am2(i-l_xm  )+l_fmiccg*l_am3(i-l_xm  ))*l_rd(i-l_xm  ) - &
         l_am3(i-l_xmym)*(l_fmiccg*l_am1(i-l_xmym)+l_fmiccg*l_am2(i-l_xmym)+       l_am3(i-l_xmym))*l_rd(i-l_xmym) )
   end do

   do i = 1, l_xmymzm
      l_ad(i) = l_ad(i)*l_rd(i)
      l_rd(i) = sqrt(l_rd(i))
      l_bv(i) = l_bv(i)*l_rd(i)
      l_am1(i-1   ) = l_am1(i-1   )*l_rd(i)*l_rd(i-1   )
      l_am2(i-l_xm  ) = l_am2(i-l_xm  )*l_rd(i)*l_rd(i-l_xm  )
      l_am3(i-l_xmym) = l_am3(i-l_xmym)*l_rd(i)*l_rd(i-l_xmym)
   end do

   l_inorm = ZERO
   do i = 1, l_xmymzm
      l_ad(i) = l_ad(i) - TWO

      l_bv(i) = l_bv(i) + l_am1(i-1   )*l_bv(i-1   ) &
                    + l_am2(i-l_xm  )*l_bv(i-l_xm  ) &
                    + l_am3(i-l_xmym)*l_bv(i-l_xmym)
      l_inorm = l_inorm + abs(l_bv(i))
   end do

   do i = l_xmymzm, 1, -1
      l_tv(i) = xs(i) + l_am1(i     )*l_tv(i+1     ) &
                      + l_am2(i     )*l_tv(i+l_xm  ) &
                      + l_am3(i     )*l_tv(i+l_xmym)
   end do
   do i = 1, l_xmymzm
      l_zv(i) = xs(i) + l_ad (i       )*l_tv(i       ) &
                      + l_am1(i-1     )*l_zv(i-1     ) &
                      + l_am2(i-l_xm  )*l_zv(i-l_xm  ) &
                      + l_am3(i-l_xmym)*l_zv(i-l_xmym)
   end do
   bdotb1 = ZERO
   do i = l_xmymzm, 1, -1
      l_zv(i) = l_zv(i) + l_tv(i)
      l_bv(i) = l_bv(i) - l_zv(i)

      ! iteration 0.

      bdotb1 = bdotb1 + l_bv(i)*l_bv(i)
      l_pv(i)  = l_bv(i)

      ! first step of the matrix vector multiplication, see below

      l_tv(i) = l_pv(i) + l_am1(i     )*l_tv(i+1     ) &
                        + l_am2(i     )*l_tv(i+l_xm  ) &
                        + l_am3(i     )*l_tv(i+l_xmym)
   end do

   l_itn = 0
   uconvg = .true.

   ! the main loop of iccg solver

   do while ( uconvg )

      ! second and third steps of the matrix vector multiplication

      pdotz = ZERO
      do i = 1, l_xmymzm+l_xmym
         l_zv(i) = l_pv(i) + l_ad (i     )*l_tv(i     ) &
                           + l_am1(i-1   )*l_zv(i-1   ) &
                           + l_am2(i-l_xm  )*l_zv(i-l_xm  ) &
                           + l_am3(i-l_xmym)*l_zv(i-l_xmym)

         j = i - l_xmym
         l_zv(j) = l_zv(j) + l_tv(j)

         pdotz = pdotz + l_pv(j)*l_zv(j)
      end do
      alpha = bdotb1/pdotz

      l_norm = ZERO
      bdotb2 = ZERO
      l_itn = l_itn + 1
      do i = 1, l_xmymzm
         xs(i) = xs(i) + alpha*l_pv(i)
         l_bv(i) = l_bv(i) - alpha*l_zv(i)
         l_norm  = l_norm  + abs(l_bv(i))

         bdotb2= bdotb2+ l_bv(i)*l_bv(i)
      end do

      ! check convergence

      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then

         uconvg = .false.
         if ( l_itn >= l_maxitn ) then
            write(6, *) 'PBSA Warning in pb_miccg(): CG l_maxitn exceeded!'
         end if

      else

         beta = bdotb2/bdotb1
         bdotb1 = bdotb2

         ! first step of the matrix vector multiplication

         do i = l_xmymzm, 1, -1
            l_pv(i) = l_bv(i) + beta*l_pv(i)

            l_tv(i) = l_pv(i) + l_am1(i)*l_tv(i+1   ) &
                          + l_am2(i)*l_tv(i+l_xm  ) &
                          + l_am3(i)*l_tv(i+l_xmym)
         end do
      end if
   end do  !  while ( uconvg ), end of the main iccg loop

   ! back scaling of the solution

   do i = l_xmymzm, 1, -1
      l_tv(i)  = xs(i) + l_am1(i)*l_tv(i+1     ) &
                       + l_am2(i)*l_tv(i+l_xm  ) &
                       + l_am3(i)*l_tv(i+l_xmym)

      phi(i) = l_tv(i)*l_rd(i)
   end do

end subroutine pb_iccg

subroutine pb_cg(phi,xs)

    implicit none

! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

! Local variables

   logical uconvg
   integer i, j
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

   l_inorm = sum(abs(l_bv(1:l_xmymzm)))

!  compute b - A * x(0) and save it in r(0)
!  p(0) = r(0)
!
!  iteration 0:
!  compute <r(0),r(0)>
!
   l_itn = 0
   bdotb1 = ZERO
   do i = 1,l_xmymzm
      l_bv(i)  = l_bv(i)  + l_am3(i-l_xmym)*xs(i-l_xmym) &
                      + l_am2(i-l_xm  )*xs(i-l_xm  ) &
                      + l_am1(i-1   )*xs(i-1   ) &
                      - l_ad(i)      *xs(i     ) &
                      + l_am1(i     )*xs(i+1   ) &
                      + l_am2(i     )*xs(i+l_xm  ) &
                      + l_am3(i     )*xs(i+l_xmym)
      bdotb1 = bdotb1 + l_bv(i)*l_bv(i)
      l_pv(i)  = l_bv(i)
   end do
!
! the main loop of the CG solver
!
   uconvg = .true.
   do while ( uconvg )
!
! iteration i:
!
!   compute Ap(i) = A * p(i)
!   compute alpha(i) = <r(i),r(i)>/<p(i),Ap(i)>
!
      pdotz = ZERO
      do i = 1,l_xmymzm
         l_zv(i) = - l_am3(i-l_xmym)*l_pv(i-l_xmym) &
                 - l_am2(i-l_xm  )*l_pv(i-l_xm  ) &
                 - l_am1(i-1   )*l_pv(i-1   ) &
                 - l_am1(i     )*l_pv(i+1   ) &
                 - l_am2(i     )*l_pv(i+l_xm  ) &
                 - l_am3(i     )*l_pv(i+l_xmym) &
                 + l_ad(i)      *l_pv(i)
         pdotz = pdotz + l_pv(i)*l_zv(i)
      end do
!
! iteration i+1:
!
      l_itn = l_itn + 1
!
!   update x(i+1) = x(i) + alpha(i) p(i)
!          r(i+1) = r(i) - alpha(i) Ap(i)
!
      alpha  = bdotb1/pdotz
      l_norm   = ZERO
      bdotb2 = ZERO
      do i = 1,l_xmymzm
         xs(i) = xs(i) + alpha*l_pv(i)
         l_bv(i) = l_bv(i) - alpha*l_zv(i)
         l_norm  = l_norm  +  abs(l_bv(i))
!
!   compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part one
!
         bdotb2 = bdotb2 + l_bv(i)*l_bv(i)
      end do
!     write(6, *)  'itn & norm ',l_itn, l_norm
!
!   check convergence
!
      if ( l_itn .ge. l_maxitn .or. l_norm .le. l_accept*l_inorm ) then

         uconvg = .false.
         if ( l_itn .ge. l_maxitn ) then
            write(6, *) 'PBSA WARNING: CG l_maxitn exceeded!'
         endif

      else
!
!   compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part two
!
         beta   = bdotb2/bdotb1
         bdotb1 = bdotb2
!
!   update p(i+1) = r(i+1) + beta(i) p(i)
!
         l_pv(1:l_xmymzm) = l_bv(1:l_xmymzm) + beta*l_pv(1:l_xmymzm)
      endif
   enddo
!
! end of the main CG loop
!
!
   phi(1:l_xmymzm) = xs(1:l_xmymzm)


end subroutine pb_cg

subroutine pb_setupper( l_am1, l_am2, l_am3 )

   implicit none

   _REAL_ l_am1(l_xm,l_ym,l_zm), l_am2(l_xm,l_ym,l_zm), l_am3(l_xm,l_ym,l_zm)

   integer i,j,k

   do j = 1, l_ym; do k = 1, l_zm
      l_am1(l_xm,j,k) = ZERO
   end do; end do
   do i = 1, l_xm; do k = 1, l_zm
      l_am2(i,l_ym,k) = ZERO
   end do; end do
   do i = 1, l_xm; do j = 1, l_ym
      l_am3(i,j,l_zm) = ZERO
   end do; end do

end subroutine pb_setupper

end module iimaug

