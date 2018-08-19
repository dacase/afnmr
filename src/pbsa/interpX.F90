#include "../include/dprec.fh"
#if !defined SANDER && !defined LIBPBSA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ twosided interpolation
subroutine twosided(l,m,n,bi,bo,n3,i0,j0,k0,nbnd,xp,yp,zp,uu,dudx,dudy,dudz, &
                   xyy,xzz,xyz,t,u,phi,xs,ys,zs,h,ir)

   use iim_util
   implicit none

   integer, parameter :: nq=27

   ! passed variables

   ! l,m,n                   box dimension
   ! h                       grid spacing
   ! n3                      when n3=fitorder 10 is 2nd order and 4 is 1st order
   ! i0,j0,k0                projection pnts in grid crds
   ! nbnd                    initial dimension without clustering
   ! xp,yp,zp                projection crds
   ! uu,dudx,dudy,dudz       potential and derivatives of the pnt on interface
   !                         uu is potential, dudx is the normal direct field,
   !                         dudy and dudz are tangential fields.
   ! phi,u                   lvlset function; potential of whole grids
   ! wp,qp                   [u]and[betaUn]for the interface
   ! wp,wyp,wzp,wyyp,        jump conditions and derivatives of [u]
   ! wyzp,wzzp
   ! qp,qyp,qzp              jump conditions and derivatives of [beta*u]
   ! xyy,xzz,xyz             derivative of the geometry on the pnt on interface
   ! t                       crds trans matrix

   integer ir,l,m,n,nbnd
   _REAL_ phi(0:l+1,0:m+1, 0:n+1)
   _REAL_ u(1:l,1:m, 1:n)
   _REAL_ xp,yp,zp
   _REAL_ xs,ys,zs,h
   integer i0,j0,k0
   _REAL_ t(3,3)
   _REAL_ bi,bo
   _REAL_ xyy,xzz,xyz

   ! local variables

   ! xloc                    local crds distance
   ! rho,rho1                ratio of the dielectric coefficients

   integer isvd,ij,i01,j01,k01,i02,j02,k02
   _REAL_  sd2(n3),work(1000)
   _REAL_ hx,hy,hz,hmax
   _REAL_ x(0:l+1), y(0:m+1), z(0:n+1)
   _REAL_ uu,dudx,dudy,dudz,x1,y1,z1
   _REAL_ ds
   _REAL_ xloc(3),tmpx(3)
   _REAL_ bl_out(3),bl_in(3)
   _REAL_ fjmp,fkjmp,ujmp
   _REAL_ rho,rho1
   _REAL_ bl2jmp,bl3jmp,temp
   integer nc,i,j,k,job,n3,nsub,info
   _REAL_ wy,wz,wyy,wyz,wzz,q,qy,qz
   _REAL_ fbi,fbo

   integer, parameter :: nn = 100
   _REAL_,dimension(1:n3,1:nn) :: ca
   _REAL_,dimension(1:nn) :: b,bc
   _REAL_,dimension(1:nn) :: ew,w4
   _REAL_,dimension(1:nn) :: w3ux,w3uy,w3uz,w3u
   _REAL_,dimension(1:nn) :: ewux,ewuy,ewuz,ewu
   _REAL_,dimension(1:n3,1:n3) :: uw
   _REAL_,dimension(1:nn,1:nn) :: vl
   _REAL_,dimension(1:n3+1) :: sd
   _REAL_ dist,dist0

   ! initialization

   hx = h; hy = h; hz = h; hmax = h

   do i = 0, l+1
      x(i) = xs + i*hx
   end do
   do j = 0, m+1
      y(j) = ys + j*hy
   end do
   do k = 0, n+1
      z(k) = zs + k*hz
   end do

   x1=xp ! (xp,yp,zp) is already in the lab frame
   y1=yp
   z1=zp

   fjmp = ZERO
   fkjmp = ZERO

   ujmp = wp(ir)
   wy   = wyp(ir)
   wz   = wzp(ir)
   wyy  = wyyp(ir)
   wzz  = wzzp(ir)
   wyz  = wyzp(ir)

   q  =  qp(ir)
   qy = qyp(ir)
   qz = qzp(ir)

   bl_out= ZERO ; bl_in= ZERO
   rho = bi/bo
   rho1 = 1 - rho
   bl2jmp = bl_out(2)-bl_in(2)
   bl3jmp = bl_out(3)-bl_in(3)
   temp = fkjmp/bo

   ! This step is to make sure the center point is the nearest irregular point.
   ! If i0, j0, k0 did these before call twosided, we can skip this step
   ! i0 = nint((xp - x(0))/h)
   ! j0 = nint((yp - y(0))/h)
   ! k0 = nint((zp - z(0))/h)

   dist0=999.0d0
   i02 = 0; j02 = 0; k02 = 0
   do i = -1, 1
   do j = -1, 1
   do k = -1, 1
      i01 = i0 + i
      j01 = j0 + j
      k01 = k0 + k
      dist = (xp-x(i01))**2+(yp-y(j01))**2+(zp-z(k01))**2
      if ( dist < dist0 ) then
         i02 = i01
         j02 = j01
         k02 = k01
         dist0 = dist
       end if
   end do
   end do
   end do

   ! set up the linear systems for all 27 grid points in the cubic searching box

   nc = 0
   do i = i02-1, i02+1
   do j = j02-1, j02+1
   do k = k02-1, k02+1
      tmpx(1) = x(i) - x1
      tmpx(2) = y(j) - y1
      tmpx(3) = z(k) - z1
      nc = nc + 1
      ds = ONE ! this is to weight each grid differently if needed.
      call matvec(3,3,t,tmpx,xloc)
      if (phi(i,j,k) <= ZERO ) then
         ca(1,nc)  = ONE*ds
         ca(2,nc)  = xloc(1)*ds
         ca(3,nc)  = xloc(2)*ds
         ca(4,nc)  = xloc(3)*ds
         if ( n3 == 10 ) then
            ca(5,nc)  = HALF*xloc(1)*xloc(1)*ds
            ca(6,nc)  = HALF*xloc(2)*xloc(2)*ds
            ca(7,nc)  = HALF*xloc(3)*xloc(3)*ds
            ca(8,nc)  = xloc(1)*xloc(2)*ds
            ca(9,nc)  = xloc(1)*xloc(3)*ds
            ca(10,nc) = xloc(2)*xloc(3)*ds
         end if
         b(nc)  = u(i,j,k)*ds
      else if (phi(i,j,k) > ZERO ) then
         if ( n3 == 4 ) then
            ca(1,nc) = ONE*ds
            ca(2,nc) = rho*xloc(1)
            ca(3,nc) = xloc(2)
            ca(4,nc) = xloc(3)
            b(nc)    = u(i,j,k)*ds - ujmp*ds &
                       - xloc(1)*q/bo*ds &
                       - xloc(2)*wy*ds &
                       - xloc(3)*wz*ds
         end if
         if (n3 == 10) then
            ca(1,nc)  = ONE*ds-HALF*temp*xloc(1)*xloc(1)
            ca(2,nc)  = rho*xloc(1) &
                        - HALF*xloc(1)*xloc(1)*(rho1*(xyy+xzz) &
                        - bl_in(1)/bo + rho*bl_out(1)/bo) &
                        + HALF*xloc(2)*xloc(2)*xyy*rho1 &
                        + HALF*xloc(3)*xloc(3)*xzz*rho1 &
                        + xloc(1)*xloc(2)*(bl_in(2)/bo &
                        - rho*bl_out(2)/bo) &
                        + xloc(1)*xloc(3)*(bl_in(3)/bo &
                        - rho*bl_out(3)/bo) &
                        + xloc(2)*xloc(3)*xyz*rho1
            ca(3,nc)  = - HALF*xloc(1)*xloc(1)*bl2jmp/bo &
                        + xloc(1)*xloc(2)*xyy*rho1 &
                        + xloc(1)*xloc(3)*xyz*rho1 &
                        + xloc(2)
            ca(4,nc)  = - HALF*xloc(1)*xloc(1)*bl3jmp/bo &
                        + xloc(1)*xloc(3)*xzz*rho1 &
                        + xloc(1)*xloc(2)*xyz*rho1 &
                        + xloc(3)
            ca(5,nc)  = HALF*xloc(1)*xloc(1)*rho*ds
            ca(6,nc)  = HALF*xloc(2)*xloc(2)*ds &
                        + HALF*(rho-1)*xloc(1)*xloc(1)*ds
            ca(7,nc)  = HALF*xloc(3)*xloc(3)*ds &
                        + HALF*(rho-1)*xloc(1)*xloc(1)*ds
            ca(8,nc)  = xloc(1)*xloc(2)*rho*ds
            ca(9,nc)  = xloc(1)*xloc(3)*rho*ds
            ca(10,nc) = xloc(2)*xloc(3)*ds

            b(nc)    = u(i,j,k)*ds - ujmp*ds &
                       - xloc(1)*q/bo*ds &
                       - xloc(2)*wy*ds &
                       - xloc(3)*wz*ds &
                       - xloc(2)*xloc(3)*(-q/bo*xyz+wyz) &
                       - HALF*xloc(2)*xloc(2)*(-q/bo*xyy+wyy) &
                       - HALF*xloc(3)*xloc(3)*(-q/bo*xzz+wzz) &
                       - xloc(1)*xloc(2)*(wy*xyy+wz*xyz+qy/bo) &
                       - xloc(1)*xloc(3)*(wy*xyz+wz*xzz+qz/bo)&
                       - HALF*xloc(1)*xloc(1)*(q/bo*(xyy+xzz)-wyy-wzz)
         end if
      end if
   end do
   end do
   end do

   ! now send to SVD

   job = 11
   nsub = nc
   isvd = 2
   if ( isvd == 1 ) then
      call dsvdc(ca,n3,n3,nsub,sd,ew,uw,n3,vl,nn,w4,job,info)
   else
      call dgesvd('A','A',n3,nsub,ca,n3,sd2,uw,n3,vl,nn,work,1000,info)
      do ij = 1, n3
         sd(ij) = sd2(ij)
      end do
   end if

   ! calculate E_in using the returned least-squared coefficients

   do i = 1, n3
      if ( sd(i) > 1.0d-12 ) then
         ewu(i)  = uw(1,i)/sd(i)
         ewux(i) = uw(2,i)/sd(i)
         ewuy(i) = uw(3,i)/sd(i)
         ewuz(i) = uw(4,i)/sd(i)
      else
         ewu(i)  = ZERO
         ewux(i) = ZERO
         ewuy(i) = ZERO
         ewuz(i) = ZERO
      end if
   end do

   if ( isvd == 1 ) then
      do i = 1, nsub
         w3u(i)  = ZERO
         w3ux(i) = ZERO
         w3uy(i) = ZERO
         w3uz(i) = ZERO
         do j = 1, n3
            w3u(i)  = w3u(i)  + vl(i,j)*ewu(j)
            w3ux(i) = w3ux(i) + vl(i,j)*ewux(j)
            w3uy(i) = w3uy(i) + vl(i,j)*ewuy(j)
            w3uz(i) = w3uz(i) + vl(i,j)*ewuz(j)
         end do
      end do
   else
      do i = 1, nsub
         w3u(i)  = ZERO
         w3ux(i) = ZERO
         w3uy(i) = ZERO
         w3uz(i) = ZERO
         do j = 1, n3
            w3u(i)  = w3u(i)  + vl(j,i)*ewu(j)
            w3ux(i) = w3ux(i) + vl(j,i)*ewux(j)
            w3uy(i) = w3uy(i) + vl(j,i)*ewuy(j)
            w3uz(i) = w3uz(i) + vl(j,i)*ewuz(j)
         end do
      end do
   end if

   uu   = ZERO
   dudx = ZERO
   dudy = ZERO
   dudz = ZERO

   do i = 1, nsub
      uu   = uu   +  w3u(i) *b(i)
      dudx = dudx +  w3ux(i)*b(i)
      dudy = dudy +  w3uy(i)*b(i)
      dudz = dudz +  w3uz(i)*b(i)
   end do

end subroutine  twosided
#endif /*ndef SANDER or LIBPBSA*/
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine onesided(l,m,n,n3,n2,xp,yp,zp,up,dudx,dudy,dudz,u,phi,h)

   implicit none

#  include "pb_constants.h"

! passed variables:
! phip: flag to choose inside or outside points
!       -1, inside
!       +1, outside
! n2: no. of grid points to be used for least square fitting
! n3: no. of unknowns to be fitted
!        4, first-order linear fitting
!       10, second-roder quadratic fitting
! (Mengjuei strongly opposes variable name phi.)

   integer l,m,n
   integer n2,n3
   _REAL_  xp,yp,zp
   _REAL_  u(1:l,1:m,1:n), phi(0:l+1,0:m+1,0:n+1)
   _REAL_  up,dudx,dudy,dudz,h

   ! local variable

   integer job,info,isvd
   integer ix0,iy0,iz0,ix1,iy1,iz1,ix,iy,iz
   integer nsub
   integer i1, i2, j2, k2, i, j, k, ij
   integer, parameter :: nn = 100
   _REAL_  w1(1:n3,1:nn)
   _REAL_  b(1:nn)
   _REAL_  sd(1:n3+1)
   _REAL_  ew(1:nn),w4(1:nn)
   _REAL_  w3u(1:nn),w3ux(1:nn),w3uy(1:nn),w3uz(1:nn)
   _REAL_  ewu(1:nn),ewux(1:n3),ewuy(1:n3),ewuz(1:n3)
   _REAL_  uw(1:n3,1:n3)
   _REAL_  vl(1:nn,1:nn)
   _REAL_  dist, dist0
   _REAL_  dx,dy,dz
   _REAL_  sd2(n3),work(1000)

   ix = 0; iy = 0; iz = 0
   isvd = 2

   ! select the grid (ix,iy,iz) closest to the surface charge (xp,yp,zp)

   ix0 = nint(xp)
   iy0 = nint(yp)
   iz0 = nint(zp)
   dist0 = 9999.0d0
   do i = -1, 1
   do j = -1, 1
   do k = -1, 1
      ix1 = ix0 + i
      iy1 = iy0 + j
      iz1 = iz0 + k
      if ( phi(ix1,iy1,iz1) > ZERO ) cycle ! only choose interior points
      dist = (xp-ix1)**2+(yp-iy1)**2+(zp-iz1)**2
      if ( dist < dist0 ) then
         ix = ix1
         iy = iy1
         iz = iz1
         dist0 = dist
      end if
   end do
   end do
   end do

   ! select the closest inside grid points to interplate E_in
   ! and construcut the matrix for SVD

   nsub = 0
   i1 = 0
   do while ( nsub < n2 )
      do i2 = -i1, i1
      do j2 = -i1, i1
      do k2 = -i1, i1
         i = i2 + ix
         j = j2 + iy
         k = k2 + iz
         if ( i > l .or. i < 1 ) cycle
         if ( j > m .or. j < 1 ) cycle
         if ( k > n .or. k < 1 ) cycle
         if ( i2*i2 + j2*j2 + k2*k2 == i1 ) then
            if ( phi(i,j,k) <= ZERO ) then
               nsub = nsub + 1
               dx = (i - xp)*h
               dy = (j - yp)*h
               dz = (k - zp)*h
               w1(1,nsub) = ONE
               w1(2,nsub) = dx
               w1(3,nsub) = dy
               w1(4,nsub) = dz
               if ( n3 == 10 ) then
                  w1(5,nsub) = HALF*dx*dx
                  w1(6,nsub) = HALF*dy*dy
                  w1(7,nsub) = HALF*dz*dz
                  w1(8,nsub) = dx*dy
                  w1(9,nsub) = dx*dz
                  w1(10,nsub) = dy*dz
               end if
               b(nsub) = u(i,j,k)
            end if
         end if
      end do
      end do
      end do
      i1 = i1 + 1
   end do

   job = 11
   if ( isvd == 1 ) then
      call dsvdc(w1,n3,n3,nsub,sd,ew,uw,n3,vl,nn,w4,job,info)
   else
      call dgesvd('A','A',n3,nsub,w1,n3,sd2,uw,n3,vl,nn,work,1000,info)
      do ij = 1, n3
         sd(ij) = sd2(ij)
      end do
   end if

   ! calculate E_in using the returned least-squared coefficients

   do i = 1, n3
      if ( sd(i) > 1.0d-12 ) then
         ewu(i)  = uw(1,i)/sd(i)
         ewux(i) = uw(2,i)/sd(i)
         ewuy(i) = uw(3,i)/sd(i)
         ewuz(i) = uw(4,i)/sd(i)
      else
         ewu(i)  = ZERO
         ewux(i) = ZERO
         ewuy(i) = ZERO
         ewuz(i) = ZERO
      endif
   enddo

   if ( isvd == 1 ) then
      do i = 1, nsub
         w3u(i)  = ZERO
         w3ux(i) = ZERO
         w3uy(i) = ZERO
         w3uz(i) = ZERO
         do j = 1, n3
            w3u(i)  = w3u(i)  + vl(i,j)*ewu(j)
            w3ux(i) = w3ux(i) + vl(i,j)*ewux(j)
            w3uy(i) = w3uy(i) + vl(i,j)*ewuy(j)
            w3uz(i) = w3uz(i) + vl(i,j)*ewuz(j)
         enddo
      enddo
   else
      do i = 1, nsub
         w3u(i)  = ZERO
         w3ux(i) = ZERO
         w3uy(i) = ZERO
         w3uz(i) = ZERO
         do j = 1, n3
            w3u(i)  = w3u(i)  + vl(j,i)*ewu(j)
            w3ux(i) = w3ux(i) + vl(j,i)*ewux(j)
            w3uy(i) = w3uy(i) + vl(j,i)*ewuy(j)
            w3uz(i) = w3uz(i) + vl(j,i)*ewuz(j)
         enddo
      enddo
   end if

   up   = ZERO
   dudx = ZERO
   dudy = ZERO
   dudz = ZERO
   do i = 1, nsub
      up   = up   +  w3u(i) *b(i)
      dudx = dudx +  w3ux(i)*b(i)
      dudy = dudy +  w3uy(i)*b(i)
      dudz = dudz +  w3uz(i)*b(i)
   end do

end subroutine onesided
