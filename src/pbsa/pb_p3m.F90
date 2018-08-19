#include "../include/dprec.fh"
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)

module pb_p3m

   use density_surface, only : tempcnt, bndatmptr, bndatmlst, cirreg

#  include "pb_constants.h"

   integer ncrg
   integer, allocatable :: tmpind(:,:,:)
   _REAL_, allocatable :: grdcrd(:,:)
   _REAL_, allocatable :: grdcrg(:)
   _REAL_, allocatable :: du(:,:)
   _REAL_, allocatable :: aphi(:,:,:,:)
   integer, allocatable :: ndenatm(:)

   contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ initializaiton of p3m routines
subroutine p3m_init( p3mopt,nbnd,natom,xm,ym,zm,iepsav,acrd,acrg,icrd,gcrg,&
              iar1pb,iprshrt,phi,cphi,chgrd,repsin )

   implicit none

   ! p3mopt controls what lists are generated, i.e. PM correction or
   ! PP direct addback.

   integer p3mopt
   integer nbnd, natom
   integer xm, ym, zm
   integer iepsav(4,xm*ym*zm)
   _REAL_ acrd(3,natom), acrg(natom)
   integer icrd(3,natom)
   _REAL_ gcrg(8,natom)
   integer iar1pb(4,0:natom)
   integer iprshrt(*)
   _REAL_ phi(xm,ym,zm), cphi(xm,ym,zm), chgrd(xm,ym,zm)
   _REAL_ repsin

   integer i, j, k, i0, j0, k0, l, cnt
   integer atmfirst, atmlast, ier
   integer iatm, jatm, jp, ibnd, jbnd
   _REAL_ srfcrg

   ! option 1: generate lists for FD/PM correction

   if ( p3mopt == 1 ) then

      ! allocation and preparation
      ! the following should be done sparsely during MD simulations
      ! so they are only called at most once per MD step.

      ! grdcrg is used for both PM correction and PP addback
      ! grdcrd is used for PP addback.
      ! iepsav is used for PM correction and borrowed to store the
      ! grid indices, so not reallocated

      if ( allocated(tmpind) ) then
         deallocate(tmpind, stat = ier); REQUIRE(ier==0)
      end if
      allocate(tmpind(3,8,natom), stat = ier); REQUIRE(ier==0)
      if ( allocated(du) ) then
         deallocate(du, stat = ier); REQUIRE(ier==0)
      end if
      allocate(du(3,natom), stat = ier); REQUIRE(ier==0)
      if ( allocated(grdcrg) ) then
         deallocate(grdcrg, stat = ier); REQUIRE(ier==0)
      end if
      allocate(grdcrg(nbnd+8*natom), stat = ier); REQUIRE(ier==0)
      if ( allocated(grdcrd) ) then
         deallocate(grdcrd, stat = ier); REQUIRE(ier==0)
      end if
      allocate(grdcrd(3,nbnd+8*natom), stat = ier); REQUIRE(ier==0)
      if ( allocated(aphi) ) then
         deallocate(aphi, stat = ier); REQUIRE(ier==0)
      end if
      allocate (aphi(-1:1,-1:1,-1:1,natom), stat = ier); REQUIRE(ier==0)
      if ( allocated(ndenatm) ) then
         deallocate(ndenatm, stat = ier); REQUIRE(ier==0)
      end if
      allocate(ndenatm(natom), stat = ier); REQUIRE(ier==0)

      ! the charge list includes both boundary charges and 8 corners
      ! of each atom in PM correction

      ncrg = nbnd + 8*natom

      ! compute and pack the boundary charges into the global list
      ! iepsav is used so grid indices are already there.

      srfcrg = ZERO
      call get_charge_pol(nbnd,xm,ym,zm,phi,cphi,chgrd,grdcrg,srfcrg)

      ! pack atom grid charges and indices into the lists

      cnt = nbnd
      do iatm = 1, natom

         ! stack integer indices for the 8 corners into tmp array

         i0 = icrd(1,iatm); j0 = icrd(2,iatm); k0 = icrd(3,iatm)
         l = 0
         do k = 0, 1
         do j = 0, 1
         do i = 0, 1
            l = l + 1
            tmpind(1,l,iatm) = i0 + i
            tmpind(2,l,iatm) = j0 + j
            tmpind(3,l,iatm) = k0 + k
         end do
         end do
         end do

         ! pack atom portions into the global lists
         ! these charges should be scaled by relative epsin

         do l = 1, 8
            cnt = cnt + 1
            grdcrg(cnt) = gcrg(l,iatm)*repsin
            iepsav(1:3,cnt) = tmpind(1:3,l,iatm)
         end do

      end do

      ! prepare to pack nonbonded list
      ! the following code is revised from density.bndatm()

      ! part a: count how many neighboring atoms are there within cutfd

      ndenatm = 1 ! count self always for self energy correction
      do iatm = 1, natom
         atmfirst = iar1pb(4, iatm-1) + 1
         atmlast  = iar1pb(2, iatm)
         do jp = atmfirst, atmlast
            jatm = iprshrt(jp)
            ndenatm(iatm) = ndenatm(iatm) + 1
            ndenatm(jatm) = ndenatm(jatm) + 1
         end do
      end do

      ! part b: set up nonbonded list pointers
      ! bndatmptr, bndatmlst and tempcnt should have been allocated to be large
      ! enough.

      cnt = nbnd
      do iatm = 1, natom
         do l = 1, 8
            cnt = cnt + 1
            bndatmptr(cnt) = bndatmptr(cnt-1) + ndenatm(iatm)
         end do
      end do

      ! part c: pack the nonbonded list for the grid charges

      tempcnt(nbnd+1:nbnd+8*natom) = 0
      do iatm = 1, natom
         ! beginning of the iatom's grid charges
         ibnd = nbnd + 8*(iatm - 1)

         ! first pack iatom's neighbors
         atmfirst = iar1pb(4, iatm-1) + 1
         atmlast  = iar1pb(2, iatm)

         do jp = atmfirst, atmlast
            jatm = iprshrt(jp)
            ! beginning of the jatom's grid charges
            jbnd = nbnd + 8*(jatm - 1)
            do l = 1, 8
               tempcnt(ibnd+l) = tempcnt(ibnd+l) + 1
               bndatmlst(bndatmptr(ibnd+l-1)+tempcnt(ibnd+l)) = jatm
               tempcnt(jbnd+l) = tempcnt(jbnd+l) + 1
               bndatmlst(bndatmptr(jbnd+l-1)+tempcnt(jbnd+l)) = iatm
            end do
         end do

         ! second pack iatom itself
         do l = 1, 8
            tempcnt(ibnd+l) = tempcnt(ibnd+l) + 1
            bndatmlst(bndatmptr(ibnd+l-1)+tempcnt(ibnd+l)) = iatm
         end do
      end do

      return

   end if

   ! option 2: generate lists for PP addback

   if ( p3mopt == 2 ) then

      ! the list for direct PP portion is shorter since
      ! we are using pairwise atoms

      ncrg = nbnd + natom

      ! pack the boundary portion into the global lists
      ! the charge is the same so no need to repack
      ! the coordinates are shifted so need to repack

      grdcrd(1:3,1:nbnd) = cirreg(1:3,1:nbnd)

      ! pack the atom portion into the global lists
      ! the charge is now on the atom, so need to repack
      ! the coordinates of the atom

      cnt = nbnd
      do iatm = 1, natom
         cnt = cnt + 1
         grdcrg(cnt) = acrg(iatm)*repsin
         grdcrd(1:3,cnt) = acrd(1:3,iatm)
      end do

      ! prepare to pack nonbonded list
      ! the following code is revised from density.bndatm()

      ! part a: count how many neighboring atoms are there
      ! this is shorter for correction because there are
      ! no excluded atoms any more.

      ndenatm(1:natom) = 0
      do iatm = 1, natom
         atmfirst = iar1pb(1, iatm) + 1
         atmlast  = iar1pb(2, iatm)
         do jp = atmfirst, atmlast
            jatm = iprshrt(jp)
            ndenatm(iatm) = ndenatm(iatm) + 1
            ndenatm(jatm) = ndenatm(jatm) + 1
         end do
      end do

      ! part b: set up nonbonded list pointers

      cnt = nbnd
      do iatm = 1, natom
         cnt = cnt + 1
         bndatmptr(cnt) = bndatmptr(cnt-1) + ndenatm(iatm)
      end do

      ! part c: pack the nonbonded list for the pairwise atoms
      ! in direct pairwise summation loops

      tempcnt(nbnd+1:nbnd+natom) = 0
      do iatm = 1, natom
         ! beginning of the iatom's grid charges
         ibnd = nbnd + iatm

         ! first pack iatom's neighbors
         atmfirst = iar1pb(1, iatm) + 1
         atmlast  = iar1pb(2, iatm)
         do jp = atmfirst, atmlast
            jatm = iprshrt(jp)
            ! beginning of the jatom's grid charges
            jbnd = nbnd + jatm
            tempcnt(ibnd) = tempcnt(ibnd) + 1
            bndatmlst(bndatmptr(ibnd-1)+tempcnt(ibnd)) = jatm
            tempcnt(jbnd) = tempcnt(jbnd) + 1
            bndatmlst(bndatmptr(jbnd-1)+tempcnt(jbnd)) = iatm
         end do
      end do

      return

   end if

   if ( p3mopt /= 1 .and. p3mopt /= 2 ) then
      write(6,'(a)') 'PB Bomb in p3m_init(): Unknow p3mopt'
      call mexit(6,1)
   end if

end subroutine p3m_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ correction of short-range FD potential
subroutine p3m_fdpotential( nbnd,natom,iepsav,icrd,h,eps0 )

   implicit none

   ! Common variables green's function values

   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green

   ! Passed variables

   ! ncrg = nbnd + 8*natom to correct short p-m pairs
   ! iepsav and grdcrg contain both dielectric and atomic grid charges
   ! bndatm arrays contains pair lists for both types of grid charges
   ! here atom grid charges always appear in sets of 8.

   integer nbnd, natom
   integer iepsav(4,nbnd+8*natom)
   integer icrd(3,natom)
   _REAL_ h, eps0

   ! Local variables

   integer iatm, ip, jp
   integer atmfirst, atmlast
   !_REAL_ xi(3), xii(3), dx(3)
   _REAL_ charge, factor
   integer ix(3), iix(3), dij(3)
   integer dijx, dijxm, dijxp, dijy, dijym, dijyp, dijz, dijzm, dijzp

   ! zero it out first

   aphi = ZERO
   factor = INV_FOURPI/(h*eps0)

   ! loop over all grid charges

   do ip = 1, ncrg

      ! collect this charge's data

      iix(1:3) = iepsav(1:3, ip)
      !xii(1:3) = grdcrd(1:3, ip)
      charge = grdcrg(ip)

      ! loop over all closeby atoms of this grid point

      atmfirst = bndatmptr(ip-1) + 1
      atmlast = bndatmptr(ip)
      do jp = atmfirst, atmlast

         ! collect this atom's data

         iatm = bndatmlst(jp)
         ix(1:3) = icrd(1:3,iatm)
         !xi(1:3) = grdcrd(1:3,nbnd+iatm)

         ! the cutoff is needed for now because the density function uses a
         ! larger cutoff. We can do a two-tier list as the iar1pb() array

         ! this has to go
         !dx(1:3) = xii(1:3) - xi(1:3)
         !if ( dx(1)**2 + dx(2)**2 + dx(3)**2 > cutfd ) cycle

         ! all pairwise grid distances to be considered

         dij(1:3) = iix(1:3) - ix(1:3)

         dijxm = abs(dij(1)+1); dijym = abs(dij(2)+1); dijzm = abs(dij(3)+1)
         dijxp = abs(dij(1)-1); dijyp = abs(dij(2)-1); dijzp = abs(dij(3)-1)
         dijx  = abs(dij(1)  ); dijy  = abs(dij(2)  ); dijz  = abs(dij(3)  )

         ! correct all FD potentials to be used for energies and forces for this atom
         ! the points used are in the cube of k = -1, +1; j = -1, +1; i = -1, +1,
         ! centered at [ix, iy, iz]

         aphi(-1,-1,-1,iatm) = aphi(-1,-1,-1,iatm) + charge*green(dijxm,dijym,dijzm)
         aphi( 0,-1,-1,iatm) = aphi( 0,-1,-1,iatm) + charge*green(dijx ,dijym,dijzm)
         aphi(+1,-1,-1,iatm) = aphi(+1,-1,-1,iatm) + charge*green(dijxp,dijym,dijzm)
         aphi(-1, 0,-1,iatm) = aphi(-1, 0,-1,iatm) + charge*green(dijxm,dijy ,dijzm)
         aphi( 0, 0,-1,iatm) = aphi( 0, 0,-1,iatm) + charge*green(dijx ,dijy ,dijzm)
         aphi(+1, 0,-1,iatm) = aphi(+1, 0,-1,iatm) + charge*green(dijxp,dijy ,dijzm)
         aphi(-1,+1,-1,iatm) = aphi(-1,+1,-1,iatm) + charge*green(dijxm,dijyp,dijzm)
         aphi( 0,+1,-1,iatm) = aphi( 0,+1,-1,iatm) + charge*green(dijx ,dijyp,dijzm)
         aphi(+1,+1,-1,iatm) = aphi(+1,+1,-1,iatm) + charge*green(dijxp,dijyp,dijzm)
         aphi(-1,-1, 0,iatm) = aphi(-1,-1, 0,iatm) + charge*green(dijxm,dijym,dijz )
         aphi( 0,-1, 0,iatm) = aphi( 0,-1, 0,iatm) + charge*green(dijx ,dijym,dijz )
         aphi(+1,-1, 0,iatm) = aphi(+1,-1, 0,iatm) + charge*green(dijxp,dijym,dijz )
         aphi(-1, 0, 0,iatm) = aphi(-1, 0, 0,iatm) + charge*green(dijxm,dijy ,dijz )
         aphi( 0, 0, 0,iatm) = aphi( 0, 0, 0,iatm) + charge*green(dijx ,dijy ,dijz )
         aphi(+1, 0, 0,iatm) = aphi(+1, 0, 0,iatm) + charge*green(dijxp,dijy ,dijz )
         aphi(-1,+1, 0,iatm) = aphi(-1,+1, 0,iatm) + charge*green(dijxm,dijyp,dijz )
         aphi( 0,+1, 0,iatm) = aphi( 0,+1, 0,iatm) + charge*green(dijx ,dijyp,dijz )
         aphi(+1,+1, 0,iatm) = aphi(+1,+1, 0,iatm) + charge*green(dijxp,dijyp,dijz )
         aphi(-1,-1,+1,iatm) = aphi(-1,-1,+1,iatm) + charge*green(dijxm,dijym,dijzp)
         aphi( 0,-1,+1,iatm) = aphi( 0,-1,+1,iatm) + charge*green(dijx ,dijym,dijzp)
         aphi(+1,-1,+1,iatm) = aphi(+1,-1,+1,iatm) + charge*green(dijxp,dijym,dijzp)
         aphi(-1, 0,+1,iatm) = aphi(-1, 0,+1,iatm) + charge*green(dijxm,dijy ,dijzp)
         aphi( 0, 0,+1,iatm) = aphi( 0, 0,+1,iatm) + charge*green(dijx ,dijy ,dijzp)
         aphi(+1, 0,+1,iatm) = aphi(+1, 0,+1,iatm) + charge*green(dijxp,dijy ,dijzp)
         aphi(-1,+1,+1,iatm) = aphi(-1,+1,+1,iatm) + charge*green(dijxm,dijyp,dijzp)
         aphi( 0,+1,+1,iatm) = aphi( 0,+1,+1,iatm) + charge*green(dijx ,dijyp,dijzp)
         aphi(+1,+1,+1,iatm) = aphi(+1,+1,+1,iatm) + charge*green(dijxp,dijyp,dijzp)

      end do  !  jp = atmfirst, atmlast

   end do  !  ip = 1, ncrg

   aphi = aphi*factor
!write(502,'(3e30.12)') aphi

end subroutine p3m_fdpotential
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ qefrc that uses one-sided interpolation for energy and force
subroutine p3m_qefrc(natom,xm,ym,zm,icrd,acrd,acrg,grdnrg,pbfrc,phi)

   use density_surface, only : index, x, y, z
   implicit none

   ! Passed variables

   integer natom
   integer xm, ym, zm
   integer icrd(3,natom)
   _REAL_ acrd(3,natom), acrg(natom)
   _REAL_ grdnrg
   _REAL_ pbfrc(3,natom)
   _REAL_ phi(xm,ym,zm)

   ! Local variables

   integer iatm
   integer i0, j0, k0
   _REAL_ xp, yp, zp, charge
   _REAL_ up, dudx(3)

   grdnrg = ZERO
   pbfrc = ZERO
   do iatm = 1, natom
      i0 = icrd(1,iatm); j0 = icrd(2,iatm); k0 = icrd(3,iatm)
      xp = acrd(1,iatm); yp = acrd(2,iatm); zp = acrd(3,iatm)
      charge = acrg(iatm)

      call p3m_onesided(iatm,i0,j0,k0,xp,yp,zp,up,dudx)

      grdnrg = grdnrg + charge*up

      pbfrc(1:3,iatm) = pbfrc(1:3,iatm) - charge*dudx(1:3)
   end do
!write(602,'(3e30.12)') pbfrc

   grdnrg = grdnrg*HALF

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!+ field interpolation using grid potentials from inside of the interface only
subroutine p3m_onesided(iatm,i0,j0,k0,xp,yp,zp,up,dudx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate up and dudx/dudy/dudz for an atom using inside potentials only.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! passed variables:
   ! iatm:     atom number
   ! i0,j0,k0: grid node assigned to the atom
   ! xp,yp,zp: grid coordinates of the atom
   ! up,dudx,  up is exterpolated potential dudx is derivative

   integer iatm,i0,j0,k0
   _REAL_ xp,yp,zp
   _REAL_ up,dudx(3)

   ! local variable

   integer, parameter :: norder = 10
   integer, parameter :: nn = 100

   integer info
   integer nsub
   integer i, j, k
   _REAL_ dx(3)
   _REAL_ w1(1:norder,1:nn)
   _REAL_ b(1:nn)

   _REAL_ sd(1:norder), uw(1:norder,1:norder)
   _REAL_ vl(1:nn,1:nn)
   _REAL_ work(1000)

   _REAL_ ewu(1:norder),ewux(1:norder),ewuy(1:norder),ewuz(1:norder)
   _REAL_ w3u(1:nn),w3ux(1:nn),w3uy(1:nn),w3uz(1:nn)

   ! select inside grid points to construcut the matrix for SVD

   nsub = 0
   do i = i0-1, i0+1
   do j = j0-1, j0+1
   do k = k0-1, k0+1
      if ( index(i,j,k) .ge. 4 ) cycle
      nsub = nsub + 1
      dx(1) = x(i) - xp
      dx(2) = y(j) - yp
      dx(3) = z(k) - zp

      w1(1,nsub) = 1.0d0
      w1(2,nsub) = dx(1)
      w1(3,nsub) = dx(2)
      w1(4,nsub) = dx(3)
      if ( norder == 10 ) then
         w1(5 ,nsub) = 0.5d0*dx(1)*dx(1)
         w1(6 ,nsub) = 0.5d0*dx(2)*dx(2)
         w1(7 ,nsub) = 0.5d0*dx(3)*dx(3)
         w1(8 ,nsub) = dx(1)*dx(2)
         w1(9 ,nsub) = dx(1)*dx(3)
         w1(10,nsub) = dx(2)*dx(3)
      end if
      b(nsub) = phi(i,j,k) - aphi(i-i0,j-j0,k-k0,iatm)
   end do
   end do
   end do

   ! This is to comptue the matrices of eqn (17) of Wang et al. second IIM paper
   ! sd() is the sigma values, uw() is the u matrix, vl() is the v* matrix.

   call dgesvd('A','A',norder,nsub,w1,norder,sd,uw,norder,vl,nn,work,1000,info)
   if ( info /= 0 ) then
      write(6,'(a,2i6)') 'PB Bomb in p3m_qefrc(): SVD failed on atom', iatm, info
      call mexit(6,1)
   end if

   ! calculate E_in which means we only need the first four values of the beta() vector
   ! compute the beta vector, eqn(18) of Wang et al. IIM second paper.

   ! this multiplies u with sigma inverse.

   do i = 1, norder
      if ( abs(sd(i)) > 1.0d-9 ) then
         ewu(i)  = uw(1,i)/sd(i)
         ewux(i) = uw(2,i)/sd(i)
         ewuy(i) = uw(3,i)/sd(i)
         ewuz(i) = uw(4,i)/sd(i)
      else
         ewu(i)  = 0.0d0
         ewux(i) = 0.0d0
         ewuy(i) = 0.0d0
         ewuz(i) = 0.0d0
      end if
   end do

   ! this multiplies with v* matrix.

   do i = 1, nsub
      w3u(i)  = 0.0d0
      w3ux(i) = 0.0d0
      w3uy(i) = 0.0d0
      w3uz(i) = 0.0d0
      do j = 1, norder
         w3u(i)  = w3u(i)  + vl(j,i)*ewu(j)
         w3ux(i) = w3ux(i) + vl(j,i)*ewux(j)
         w3uy(i) = w3uy(i) + vl(j,i)*ewuy(j)
         w3uz(i) = w3uz(i) + vl(j,i)*ewuz(j)
      end do
   end do

   ! this multiplies with the f() vector, the right hand side.

   up   = 0.0d0
   dudx = 0.0d0
   do i = 1, nsub
      up      = up      + w3u(i) *b(i)
      dudx(1) = dudx(1) + w3ux(i)*b(i)
      dudx(2) = dudx(2) + w3uy(i)*b(i)
      dudx(3) = dudx(3) + w3uz(i)*b(i)
   end do

end subroutine p3m_onesided

end subroutine p3m_qefrc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dbfrc that uses a simple method to interpolate interface field
subroutine p3m_dbfrc( natom,nbnd,xm,ym,zm,h,dprob,radi,pbfrc,iepsav,&
              u,phi,eps0,epsx,epsy,epsz )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This dbfrc uses a extremely simple method to calculate D_n for an irregular grid
! The method is rooted in retrieving the correct field information from the
! finite-volume method with the weighted-harmonic average, i.e. first-order IIM,
! to treat the dielectric jump.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use density_surface, only : x, y, z, density_deriv_atom
   implicit none

   ! passed variables

   integer natom, nbnd
   integer xm, ym, zm
   _REAL_ h, dprob
   _REAL_ radi(natom)
   _REAL_ pbfrc(3,natom)
   integer iepsav(4,nbnd)
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm)
   _REAL_ eps0, epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)

   ! local variables

   integer ip, jp
   integer i, j, k
   integer iatm
   integer atmfirst, atmlast
   _REAL_ xi(3), dist, dist0, charge, rh
   _REAL_ dx(3), df(3), n(3)
   _REAL_ dn, dbf, dum(3), coef
   _REAL_ factor, weight, weight2
   _REAL_ dusum(3)

   ! initialization

!Changhao: comparing with analytical test case
!pbfrc = ZERO

   rh = -ONE/h/eps0 ! the minus sign is to convert gradient to field
   factor = ONE/(TWO*dprob)

   do ip = 1, nbnd

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip)

      xi(1:3) = grdcrd(1:3,ip)
      charge = grdcrg(ip)
!write(4602,'(a,4f10.5)') '  H', xi(1:3), sqrt(xi(1)**2+xi(2)**2+xi(3)**2)

      ! compute density partial derivatives and gradient

      du(1:3, 1:natom) = ZERO
      dusum(1:3) = ZERO
      atmfirst = bndatmptr(ip-1)+1
      atmlast = bndatmptr(ip)
      do jp = atmfirst, atmlast
         iatm = bndatmlst(jp)

         dx(1:3) = xi(1:3) - grdcrd(1:3,nbnd+iatm)
         dist0 = sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         dist = (dist0 - radi(iatm))*factor
         dx(1:3) = dx(1:3)*factor
         if ( radi(iatm) == ZERO .or. dist > ONE ) cycle

         du(1,iatm) = density_deriv_atom(dist, dist0, dx(1))
         du(2,iatm) = density_deriv_atom(dist, dist0, dx(2))
         du(3,iatm) = density_deriv_atom(dist, dist0, dx(3))
         dusum(1:3) = dusum(1:3) + du(1:3,iatm)
      end do

      ! normal/gradient unit vector

      weight2 = ONE/(dusum(1)**2+dusum(2)**2+dusum(3)**2)
      weight = sqrt(weight2)
      n(1:3) = dusum(1:3)*weight

      ! compute D vector
      ! D component is nonzero only in a direction along which there is a fractional edge
      ! the division by h is done below, i.e. *rh

      if      ( xi(1) .gt. x(i) ) then
         df(1) = (phi(i+1,j,k)-phi(i,j,k))*epsx(i,j,k)
      else if ( xi(1) .lt. x(i) ) then
         df(1) = (phi(i,j,k)-phi(i-1,j,k))*epsx(i-1,j,k)
      else
         df(1) = ZERO
      end if
      if      ( xi(2) .gt. y(j) ) then
         df(2) = (phi(i,j+1,k)-phi(i,j,k))*epsy(i,j,k)
      else if ( xi(2) .lt. y(j) ) then
         df(2) = (phi(i,j,k)-phi(i,j-1,k))*epsy(i,j-1,k)
      else
         df(2) = ZERO
      end if
      if      ( xi(3) .gt. z(k) ) then
         df(3) = (phi(i,j,k+1)-phi(i,j,k))*epsz(i,j,k)
      else if ( xi(3) .lt. z(k) ) then
         df(3) = (phi(i,j,k)-phi(i,j,k-1))*epsz(i,j,k-1)
      else
         df(3) = ZERO
      end if

      ! Dn and DBF

      dn = dot_product(df(1:3), n(1:3))
      dbf = HALF*charge*rh*dn
      dum(1:3) = dbf*n(1:3)

! Changhao: comparing with analytical test case
!
!if ( xi(1)>0.0d0 .and. xi(2)>0.0d0 .and. xi(3)>0.0d0 ) then
!   write(311,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)>0.0d0 .and. xi(3)<0.0d0 ) then
!   write(312,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)<0.0d0 .and. xi(3)>0.0d0 ) then
!   write(313,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)<0.0d0 .and. xi(3)<0.0d0 ) then
!   write(314,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)>0.0d0 .and. xi(3)>0.0d0 ) then
!   write(315,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)>0.0d0 .and. xi(3)<0.0d0 ) then
!   write(316,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)<0.0d0 .and. xi(3)>0.0d0 ) then
!   write(317,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)<0.0d0 .and. xi(3)<0.0d0 ) then
!   write(318,*) dum(1),dum(2),dum(3)
!end if
!
!if ( xi(1)>=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(411,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(412,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(413,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(414,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(415,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(416,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(417,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(418,*) dum(1),dum(2),dum(3)
!end if
!
!if ( xi(1)==0.0d0 .and. xi(2)==0.0d0 .and. xi(3)>0.0d0 ) then
!   write(511,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)==0.0d0 .and. xi(2)==0.0d0 .and. xi(3)<0.0d0 ) then
!   write(512,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)==0.0d0 .and. xi(2)>0.0d0 .and. xi(3)==0.0d0 ) then
!   write(513,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)==0.0d0 .and. xi(2)<0.0d0 .and. xi(3)==0.0d0 ) then
!   write(514,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)==0.0d0 .and. xi(3)==0.0d0 ) then
!   write(515,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)==0.0d0 .and. xi(3)==0.0d0 ) then
!   write(516,*) dum(1),dum(2),dum(3)
!end if

      ! partition the force element to atoms

      do jp = atmfirst, atmlast
         iatm = bndatmlst(jp)
         coef = dot_product(du(1:3,iatm), dusum(1:3))*weight2
         pbfrc(1:3,iatm) = pbfrc(1:3,iatm) + coef*dum(1:3)
      end do

   end do

!do iatm = 1, natom
!write(4671,'(3f15.6)') pbfrc(1:3, iatm )*FOURPI*epsin*AMBER_ELECTROSTATIC2
!enddo
!write(6,'(a,3f15.6)') 'dbfrc', pbfrc(1:3, 1)*FOURPI*epsin*AMBER_ELECTROSTATIC2

end subroutine p3m_dbfrc
#if !defined SANDER && !defined LIBPBSA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dbfrc1d that uses a simple method to interpolate interface field
subroutine p3m_dbfrc1d( natom,nbnd,xm,ym,zm,h,gox,goy,goz,dprob,radi,acrg,pbfrc,iepsav,&
              u,phi,epsin )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This dbfrc uses a extremely simple method to calculate D_n for an irregular grid
! by extrapolate reaction filed in each dimension from inside, so bcopt=6 should
! be used, but WHA or IIM should both be fine (ipb=2 or 4).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use density_surface, only : x, y, z, density_deriv_atom
   implicit none

   ! passed variables

   integer natom, nbnd
   integer xm, ym, zm
   _REAL_ h, dprob
   _REAL_ gox, goy, goz
   _REAL_ radi(natom), acrg(natom)
   _REAL_ pbfrc(3,natom)
   integer iepsav(4,nbnd)
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm)
   _REAL_ epsin

   ! local variables

   integer ip, jp
   integer i, j, k
   integer iatm
   integer atmfirst, atmlast
   _REAL_ xi(3), dist, dist0, charge
   _REAL_ dcc, coulomb(3)
   _REAL_ dx(3), df(3), n(3)
   _REAL_ dn, dbf, dum(3), coef
   _REAL_ factor, weight, weight2
   _REAL_ dusum(3)

   ! initialization

!Changhao: comparing with analytical test case
!pbfrc = ZERO

   factor = ONE/(TWO*dprob)
   do ip = 1, nbnd

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip)

      xi(1:3) = grdcrd(1:3,ip)
      charge = grdcrg(ip)

      ! compute density partial derivatives and gradient

      du(1:3, 1:natom) = ZERO
      dusum(1:3) = ZERO
      atmfirst = bndatmptr(ip-1)+1
      atmlast = bndatmptr(ip)
      do jp = atmfirst, atmlast
         iatm = bndatmlst(jp)

         dx(1:3) = xi(1:3) - grdcrd(1:3,nbnd+iatm)
         dist0 = sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         dist = (dist0 - radi(iatm))*factor
         dx(1:3) = dx(1:3)*factor
         if ( radi(iatm) == ZERO .or. dist > ONE ) cycle

         du(1,iatm) = density_deriv_atom(dist, dist0, dx(1))
         du(2,iatm) = density_deriv_atom(dist, dist0, dx(2))
         du(3,iatm) = density_deriv_atom(dist, dist0, dx(3))
         dusum(1:3) = dusum(1:3) + du(1:3,iatm)
      end do

      ! normal/gradient unit vector

      weight2 = ONE/(dusum(1)**2+dusum(2)**2+dusum(3)**2)
      weight = sqrt(weight2)
      n(1:3) = dusum(1:3)*weight

      ! compute D vector
      ! first compute reaction field

      call onedimen(xi(1),xi(2),xi(3),df(1),df(2),df(3))
      df = -df

      ! second compute coulomb field

      coulomb(1:3) = ZERO
      do iatm = 1, natom
         dx(1:3) = xi(1:3) - grdcrd(1:3,nbnd+iatm)
         dist0 = ONE/(dx(1)**2 + dx(2)**2 + dx(3)**2); dist = sqrt(dist0)

         dcc = acrg(iatm)*dist0*dist

         coulomb(1) = coulomb(1) + dx(1)*dcc
         coulomb(2) = coulomb(2) + dx(2)*dcc
         coulomb(3) = coulomb(3) + dx(3)*dcc
      end do
      coulomb = coulomb/(FOURPI*epsin)

      ! third this is D

      df(1:3) = df(1:3) + coulomb(1:3)

      ! Dn and DBF

      dn = dot_product(df(1:3), n(1:3))
      dbf = HALF*charge*dn
      dum(1:3) = dbf*n(1:3)

!Changhao: comparing with analytical test case
!
!if ( xi(1)>0.0d0 .and. xi(2)>0.0d0 .and. xi(3)>0.0d0 ) then
!   write(311,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)>0.0d0 .and. xi(3)<0.0d0 ) then
!   write(312,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)<0.0d0 .and. xi(3)>0.0d0 ) then
!   write(313,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)<0.0d0 .and. xi(3)<0.0d0 ) then
!   write(314,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)>0.0d0 .and. xi(3)>0.0d0 ) then
!   write(315,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)>0.0d0 .and. xi(3)<0.0d0 ) then
!   write(316,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)<0.0d0 .and. xi(3)>0.0d0 ) then
!   write(317,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)<0.0d0 .and. xi(3)<0.0d0 ) then
!   write(318,*) dum(1),dum(2),dum(3)
!end if
!
!if ( xi(1)>=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(411,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(412,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(413,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(414,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(415,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(416,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(417,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(418,*) dum(1),dum(2),dum(3)
!end if
!
!if ( xi(1)==0.0d0 .and. xi(2)==0.0d0 .and. xi(3)>0.0d0 ) then
!   write(511,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)==0.0d0 .and. xi(2)==0.0d0 .and. xi(3)<0.0d0 ) then
!   write(512,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)==0.0d0 .and. xi(2)>0.0d0 .and. xi(3)==0.0d0 ) then
!   write(513,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)==0.0d0 .and. xi(2)<0.0d0 .and. xi(3)==0.0d0 ) then
!   write(514,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)==0.0d0 .and. xi(3)==0.0d0 ) then
!   write(515,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)==0.0d0 .and. xi(3)==0.0d0 ) then
!   write(516,*) dum(1),dum(2),dum(3)
!end if

      ! partition the force element to atoms

      do jp = atmfirst, atmlast
         iatm = bndatmlst(jp)
         coef = dot_product(du(1:3,iatm), dusum(1:3))*weight2
         pbfrc(1:3,iatm) = pbfrc(1:3,iatm) + coef*dum(1:3)
      end do

   end do

!do iatm = 1, natom
!write(4671,'(3f15.6)') pbfrc(1:3, iatm)*FOURPI*epsin*AMBER_ELECTROSTATIC2
!enddo
!write(6,'(a,3f15.6)') 'dbfrc', pbfrc(1:3, 1)*FOURPI*epsin*AMBER_ELECTROSTATIC2

contains

subroutine onedimen(xp,yp,zp,dudx,dudy,dudz)

   implicit none

   _REAL_ xp,yp,zp
   _REAL_ dudx,dudy,dudz

   ! local variable

   integer ix0,iy0,iz0,ix1,iy1,iz1,ix,iy,iz
   integer nsub
   integer i2, j2, k2, i, j, k
   _REAL_  w1(1:3,1:4)
   _REAL_  b(1:3,1:4)

   _REAL_ dist, dist0
   _REAL_ dx,dy,dz

   ! select the inner grid point (ix,iy,iz) closest to the surface charge (xp,yp,zp)
   ! this is not necessarily the original grid point to be projected.

   ix = 0; iy = 0; iz = 0
   ix0 = nint((xp-gox)/h)
   iy0 = nint((yp-goy)/h)
   iz0 = nint((zp-goz)/h)
   dist0 = 9999.0d0
   do i = -1, 1
   do j = -1, 1
   do k = -1, 1
      ix1 = ix0 + i
      iy1 = iy0 + j
      iz1 = iz0 + k
      if ( u(ix1,iy1,iz1) > ZERO ) cycle
      dist = ((xp-gox)/h-ix1)**2+((yp-goy)/h-iy1)**2+((zp-goz)/h-iz1)**2
      if ( dist < dist0 ) then
         ix = ix1
         iy = iy1
         iz = iz1
         dist0 = dist
      end if
   end do
   end do
   end do
   if ( ix-2 < 1 .or. ix+2 > xm ) then
      write(6,'(a)') 'PB Bomb in onedimen(): project point too close to the grid edge'
      call mexit(6,1)
   end if
   if ( iy-2 < 1 .or. iy+2 > ym ) then
      write(6,'(a)') 'PB Bomb in onedimen(): project point too close to the grid edge'
      call mexit(6,1)
   end if
   if ( iz-2 < 1 .or. iz+2 > zm ) then
      write(6,'(a)') 'PB Bomb in onedimen(): project point too close to the grid edge'
      call mexit(6,1)
   end if

   ! select the closest innner grid points to interplate E_in
   ! No need for SVD any more because it is a deterministic linear system

   ! working on the x-dimension

   nsub = 0
   if ( u(ix+1,iy,iz)<=0.0d0 .and. u(ix-1,iy,iz)<=0.0d0 ) then
      do i2 = -1, 1
         i = i2+ix
         j = iy
         k = iz
         nsub = nsub + 1
         dx = (i - (xp-gox)/h)*h ! this is i0
         w1(1,nsub) = dx
         b(1,nsub) = phi(i,j,k)
      end do
   else if ( u(ix-1,iy,iz)<=0.0d0 ) then
      do i2 = -2, 0
         i = i2 + ix
         j = iy
         k = iz
         if ( u(i,j,k) <= ZERO ) then
         nsub = nsub + 1
         dx = (i - (xp-gox)/h)*h
         w1(1,nsub) = dx
         b(1,nsub) = phi(i,j,k)
         end if
      end do
   else if ( u(ix+1,iy,iz)<=0.0d0 ) then
      do i2 = 0, 2
         i = i2 + ix
         j = iy
         k = iz
         if ( u(i,j,k) <= ZERO ) then ! phip = -1 for interior points
         nsub = nsub + 1
         dx = (i - (xp-gox)/h)*h
         w1(1,nsub) = dx
         b(1,nsub) = phi(i,j,k)
         end if
      end do
   else if ( u(ix+1,iy,iz)>0.0d0 .and. u(ix-1,iy,iz)>0.0d0 ) then
      nsub = 1
   end if

   if ( nsub > 2 )   then
      dudx = ((w1(1,1)*b(1,1)+w1(1,2)*b(1,2)+w1(1,3)*b(1,3))*3*(w1(1,1)**4+w1(1,2)**4+w1(1,3)**4)-&
         (b(1,1)+b(1,2)+b(1,3))*(w1(1,1)+w1(1,2)+w1(1,3))*(w1(1,1)**4+w1(1,2)**4+w1(1,3)**4)+&
         (b(1,1)+b(1,2)+b(1,3))*(w1(1,1)**2+w1(1,2)**2+w1(1,3)**2)*(w1(1,1)**3+w1(1,2)**3+w1(1,3)**3)-&
         (w1(1,1)**2*b(1,1)+w1(1,2)**2*b(1,2)+w1(1,3)**2*b(1,3))*3*(w1(1,1)**3+w1(1,2)**3+w1(1,3)**3)-&
         (w1(1,1)*b(1,1)+w1(1,2)*b(1,2)+w1(1,3)*b(1,3))*(w1(1,1)**2+w1(1,2)**2+w1(1,3)**2)**2+&
         (w1(1,1)**2*b(1,1)+w1(1,2)**2*b(1,2)+w1(1,3)**2*b(1,3))*(w1(1,1)+w1(1,2)+w1(1,3))*(w1(1,1)**2+w1(1,2)**2+w1(1,3)**2))/&
         (3*(w1(1,1)**2+w1(1,2)**2+w1(1,3)**2)*(w1(1,1)**4+w1(1,2)**4+w1(1,3)**4)-&
         (w1(1,1)+w1(1,2)+w1(1,3))**2*(w1(1,1)**4+w1(1,2)**4+w1(1,3)**4)-&
         3*(w1(1,1)**3+w1(1,2)**3+w1(1,3)**3)**2+&
         2*(w1(1,1)+w1(1,2)+w1(1,3))*(w1(1,1)**2+w1(1,2)**2+w1(1,3)**2)*(w1(1,1)**3+w1(1,2)**3+w1(1,3)**3)-&
         (w1(1,1)**2+w1(1,2)**2+w1(1,3)**2)**3)
   else
      dudx = ZERO
   end if

   ! working on the y-dimension

   nsub = 0
   if ( u(ix,iy-1,iz)<=0.0d0 .and. u(ix,iy+1,iz)<=0.0d0 ) then
      do j2 = -1, 1
         i = ix
         j = j2+iy
         k = iz
         nsub = nsub + 1
         dy = (j - (yp-goy)/h)*h
         w1(2,nsub) = dy
         b(2,nsub) = phi(i,j,k)
      end do
   else if ( u(ix,iy-1,iz)<=0.0d0 ) then
      do j2 = -2, 0
         i = ix
         j = j2+iy
         k = iz
         if ( u(i,j,k) <= ZERO ) then ! phip = -1 for interior points
         nsub = nsub + 1
         dy = (j - (yp-goy)/h)*h
         w1(2,nsub) = dy
         b(2,nsub) = phi(i,j,k)
         end if
      end do
   else if ( u(ix,iy+1,iz)<=0.0d0 ) then
      do j2 = 0, 2
         i = ix
         j = j2+iy
         k = iz
         if ( u(i,j,k) <= ZERO ) then ! phip = -1 for interior points
         nsub = nsub + 1
         dy = (j - (yp-goy)/h)*h
         w1(2,nsub) = dy
         b(2,nsub) = phi(i,j,k)
         end if
      end do
   else if ( u(ix,iy+1,iz)>0.0d0 .and. u(ix,iy-1,iz)>0.0d0 ) then
      nsub = 1
   end if

   if ( nsub > 2 ) then
      dudy = ((w1(2,1)*b(2,1)+w1(2,2)*b(2,2)+w1(2,3)*b(2,3))*3*(w1(2,1)**4+w1(2,2)**4+w1(2,3)**4)-&
         (b(2,1)+b(2,2)+b(2,3))*(w1(2,1)+w1(2,2)+w1(2,3))*(w1(2,1)**4+w1(2,2)**4+w1(2,3)**4)+&
         (b(2,1)+b(2,2)+b(2,3))*(w1(2,1)**2+w1(2,2)**2+w1(2,3)**2)*(w1(2,1)**3+w1(2,2)**3+w1(2,3)**3)-&
         (w1(2,1)**2*b(2,1)+w1(2,2)**2*b(2,2)+w1(2,3)**2*b(2,3))*3*(w1(2,1)**3+w1(2,2)**3+w1(2,3)**3)-&
         (w1(2,1)*b(2,1)+w1(2,2)*b(2,2)+w1(2,3)*b(2,3))*(w1(2,1)**2+w1(2,2)**2+w1(2,3)**2)**2+&
         (w1(2,1)**2*b(2,1)+w1(2,2)**2*b(2,2)+w1(2,3)**2*b(2,3))*(w1(2,1)+w1(2,2)+w1(2,3))*(w1(2,1)**2+w1(2,2)**2+w1(2,3)**2))/&
         (3*(w1(2,1)**2+w1(2,2)**2+w1(2,3)**2)*(w1(2,1)**4+w1(2,2)**4+w1(2,3)**4)-&
         (w1(2,1)+w1(2,2)+w1(2,3))**2*(w1(2,1)**4+w1(2,2)**4+w1(2,3)**4)-&
         3*(w1(2,1)**3+w1(2,2)**3+w1(2,3)**3)**2+&
         2*(w1(2,1)+w1(2,2)+w1(2,3))*(w1(2,1)**2+w1(2,2)**2+w1(2,3)**2)*(w1(2,1)**3+w1(2,2)**3+w1(2,3)**3)-&
         (w1(2,1)**2+w1(2,2)**2+w1(2,3)**2)**3)
   else
      dudy = 0
   end if

   ! working on the z-dimension

   nsub = 0
   if ( u(ix,iy,iz-1)<=0.0d0 .and. u(ix,iy,iz+1)<=0.0d0 ) then
      do k2 = -1, 1
         i = ix
         j = iy
         k = k2+iz
         nsub = nsub+1
         dz = (k - (zp-goz)/h)*h
         w1(3,nsub) = dz
         b(3,nsub) = phi(i,j,k)
      end do
   else if ( u(ix,iy,iz-1)<=0.0d0 ) then
      do k2 = -2, 0
         i = ix
         j = iy
         k = k2+iz
         if ( u(i,j,k) <= ZERO ) then ! phip = -1 for interior points
         nsub = nsub + 1
         dz = (k - (zp-goz)/h)*h
         w1(3,nsub) = dz
         b(3,nsub) = phi(i,j,k)
         end if
      end do
   else if ( u(ix,iy,iz+1)<=0.0d0 ) then
      do k2 = 0, 2
         i = ix
         j = iy
         k = k2+iz
         if ( u(i,j,k) <= ZERO ) then
         nsub = nsub + 1
         dz = (k - (zp-goz)/h)*h
         w1(3,nsub) = dz
         b(3,nsub) = phi(i,j,k)
         end if
      end do
   else if ( u(ix,iy,iz+1)>0.0d0 .and. u(ix,iy,iz-1)>0.0d0 ) then
      nsub = 1
   end if

   if ( nsub > 2 ) then
      dudz = ((w1(3,1)*b(3,1)+w1(3,2)*b(3,2)+w1(3,3)*b(3,3))*3*(w1(3,1)**4+w1(3,2)**4+w1(3,3)**4)-&
         (b(3,1)+b(3,2)+b(3,3))*(w1(3,1)+w1(3,2)+w1(3,3))*(w1(3,1)**4+w1(3,2)**4+w1(3,3)**4)+&
         (b(3,1)+b(3,2)+b(3,3))*(w1(3,1)**2+w1(3,2)**2+w1(3,3)**2)*(w1(3,1)**3+w1(3,2)**3+w1(3,3)**3)-&
         (w1(3,1)**2*b(3,1)+w1(3,2)**2*b(3,2)+w1(3,3)**2*b(3,3))*3*(w1(3,1)**3+w1(3,2)**3+w1(3,3)**3)-&
         (w1(3,1)*b(3,1)+w1(3,2)*b(3,2)+w1(3,3)*b(3,3))*(w1(3,1)**2+w1(3,2)**2+w1(3,3)**2)**2+&
         (w1(3,1)**2*b(3,1)+w1(3,2)**2*b(3,2)+w1(3,3)**2*b(3,3))*(w1(3,1)+w1(3,2)+w1(3,3))*(w1(3,1)**2+w1(3,2)**2+w1(3,3)**2))/&
         (3*(w1(3,1)**2+w1(3,2)**2+w1(3,3)**2)*(w1(3,1)**4+w1(3,2)**4+w1(3,3)**4)-&
         (w1(3,1)+w1(3,2)+w1(3,3))**2*(w1(3,1)**4+w1(3,2)**4+w1(3,3)**4)-&
         3*(w1(3,1)**3+w1(3,2)**3+w1(3,3)**3)**2+&
         2*(w1(3,1)+w1(3,2)+w1(3,3))*(w1(3,1)**2+w1(3,2)**2+w1(3,3)**2)*(w1(3,1)**3+w1(3,2)**3+w1(3,3)**3)-&
         (w1(3,1)**2+w1(3,2)**2+w1(3,3)**2)**3)
   else
      dudz = 0
   end if

end subroutine onedimen

end subroutine p3m_dbfrc1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dbfrc1side that uses the standard 1-sided method to interpolate interface field
subroutine p3m_dbfrc1side( natom,nbnd,xm,ym,zm,h,gox,goy,goz,dprob,radi,acrg,pbfrc,iepsav,&
              u,phi,epsin )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This dbfrc uses the standard 1 sided method to calculate D_n for an irregular grid
! bcopt=6 should be used, and IIM would be the most consistent (ipb=4).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use density_surface, only : x, y, z, density_deriv_atom
   implicit none

   ! passed variables

   integer natom, nbnd
   integer xm, ym, zm
   _REAL_ h, dprob
   _REAL_ gox, goy, goz
   _REAL_ radi(natom), acrg(natom)
   _REAL_ pbfrc(3,natom)
   integer iepsav(4,nbnd)
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm)
   _REAL_ epsin

   ! local variables

   integer ip, jp
   integer i, j, k
   integer iatm
   integer atmfirst, atmlast
   _REAL_ xi(3), dist, dist0, charge
   _REAL_ dcc, coulomb(3)
   _REAL_ dx(3), df(3), n(3)
   _REAL_ dn, dbf, dum(3), coef
   _REAL_ factor, weight, weight2
   _REAL_ dusum(3)
   _REAL_ xp,yp,zp,up
   _REAL_ du(1:3,natom)

   ! initialization

!Changhao: comparing with analytical test case
!pbfrc = ZERO

   factor = ONE/(TWO*dprob)
   do ip = 1, nbnd
   du = ZERO
      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip)

      xi(1:3) = grdcrd(1:3,ip)
      charge = grdcrg(ip)

      ! compute density partial derivatives and gradient

      dusum(1:3) = ZERO
      du(1:3, 1:natom) = ZERO
      atmfirst = bndatmptr(ip-1)+1
      atmlast = bndatmptr(ip)
      do jp = atmfirst, atmlast
         iatm = bndatmlst(jp)

         dx(1:3) = xi(1:3) - grdcrd(1:3,nbnd+iatm)
         dist0 = sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         dist = (dist0 - radi(iatm))*factor
         dx(1:3) = dx(1:3)*factor
         if ( radi(iatm) == ZERO .or. dist > ONE ) cycle

         du(1,iatm) = density_deriv_atom(dist, dist0, dx(1))
         du(2,iatm) = density_deriv_atom(dist, dist0, dx(2))
         du(3,iatm) = density_deriv_atom(dist, dist0, dx(3))
         dusum(1:3) = dusum(1:3) + du(1:3,iatm)
      end do

      ! normal/gradient unit vector

      weight2 = ONE/(dusum(1)**2+dusum(2)**2+dusum(3)**2)
      weight = sqrt(weight2)
      n(1:3) = dusum(1:3)*weight

      ! compute D vector
      ! first compute reaction field

      xp = (xi(1) - gox)/h
      yp = (xi(2) - goy)/h
      zp = (xi(3) - goz)/h

      call onesided(xm,ym,zm,10,27,xp,yp,zp,up,df(1),df(2),df(3),phi,u,h)
      df = -df

      ! second compute coulomb field

      coulomb(1:3) = ZERO
      do iatm = 1, natom
         dx(1:3) = xi(1:3) - grdcrd(1:3,nbnd+iatm)
         dist0 = ONE/(dx(1)**2 + dx(2)**2 + dx(3)**2); dist = sqrt(dist0)

         dcc = acrg(iatm)*dist0*dist

         coulomb(1) = coulomb(1) + dx(1)*dcc
         coulomb(2) = coulomb(2) + dx(2)*dcc
         coulomb(3) = coulomb(3) + dx(3)*dcc
      end do
      coulomb = coulomb/(FOURPI*epsin)

      ! third this is D

      df(1:3) = df(1:3) + coulomb(1:3)

      ! Dn and DBF

      dn = dot_product(df(1:3), n(1:3))
      dbf = HALF*charge*dn
      dum(1:3) = dbf*n(1:3)

      ! partition the force element to atoms

      do jp = atmfirst, atmlast
         iatm = bndatmlst(jp)
         dx(1:3) = xi(1:3) - grdcrd(1:3,nbnd+iatm)
         dist0 = sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         dist = (dist0 - radi(iatm))*factor
         dx(1:3) = dx(1:3)*factor
         if ( radi(iatm) == ZERO .or. dist > ONE ) cycle

         coef = dot_product(du(1:3,iatm), dusum(1:3))*weight2
         pbfrc(1:3,iatm) = pbfrc(1:3,iatm) + coef*dum(1:3)
      end do

   end do

!do i = 1, natom
!   write(4671,*) pbfrc(1:3,i)*FOURPI*epsin*AMBER_ELECTROSTATIC2
!end do
!write(6,'(a,3f15.6)') 'dbfrc', pbfrc(1:3, 1)*FOURPI*epsin*AMBER_ELECTROSTATIC2

end subroutine p3m_dbfrc1side
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dbfrc2side uses a standard 2sided method to interpolate interface field
subroutine p3m_dbfrc2side( natom,nbnd,xm,ym,zm,h,gox,goy,goz,dprob,radi,acrg,pbfrc,iepsav,&
              u,phi,epsin,eps0 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This dbfrc uses the standard 2 sided method to calculate D_n for an irregular grid
! bcopt=6 should be used, and IIM would be the most consistent (ipb=4).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use density_surface, only : x, y, z, density_deriv_atom
   use iim_util
   implicit none

   ! passed variables

   integer natom, nbnd
   integer xm, ym, zm
   _REAL_ h, dprob
   _REAL_ gox, goy, goz
   _REAL_ radi(natom), acrg(natom)
   _REAL_ pbfrc(3,natom)
   integer iepsav(4,nbnd)
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm)
   _REAL_ epsin
   _REAL_ eps0

   ! local variables

   integer ip, jp
   integer i, j, k, i0, j0, k0
   integer iatm
   integer atmfirst, atmlast
   _REAL_ xi(3), dist, dist0, charge
   _REAL_ dcc, coulomb(3)
   _REAL_ dx(3), df(3), n(3)
   _REAL_ dn, dbf, dum(3), coef
   _REAL_ factor, weight, weight2
   _REAL_ dusum(3)
   _REAL_ ui(1:xm,1:ym,1:zm), up
   _REAL_ t(3,3), xyy, xzz, xyz

   ! initialization

! Changhao: comparing with analytical test case
pbfrc = ZERO

   ! the following convert potential back to eletron/A unit to be consistent
   ! with IIM utility routines that are deployed here.

   do k = 1,zm
   do j = 1,ym
   do i = 1,xm
      ui(i,j,k) = phi(i,j,k)/INV_FOURPI*eps0
   end do
   end do
   end do

   factor = ONE/(TWO*dprob)
   do ip = 1, nbnd

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip)

      xi(1:3) = grdcrd(1:3,ip)
      charge = grdcrg(ip)

      ! compute density partial derivatives and gradient

      du(1:3, 1:natom) = ZERO
      dusum(1:3) = ZERO
      atmfirst = bndatmptr(ip-1)+1
      atmlast = bndatmptr(ip)
      do jp = atmfirst, atmlast
         iatm = bndatmlst(jp)

         dx(1:3) = xi(1:3) - grdcrd(1:3,nbnd+iatm)
         dist0 = sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         dist = (dist0 - radi(iatm))*factor
         dx(1:3) = dx(1:3)*factor
         if ( radi(iatm) == ZERO .or. dist > ONE ) cycle

         du(1,iatm) = density_deriv_atom(dist, dist0, dx(1))
         du(2,iatm) = density_deriv_atom(dist, dist0, dx(2))
         du(3,iatm) = density_deriv_atom(dist, dist0, dx(3))
         dusum(1:3) = dusum(1:3) + du(1:3,iatm)
      end do

      ! normal/gradient unit vector

      weight2 = ONE/(dusum(1)**2+dusum(2)**2+dusum(3)**2)
      weight = sqrt(weight2)
      n(1:3) = dusum(1:3)*weight

      ! compute D vector
      ! first compute reaction field

      i0 = nint((cirreg(1,ip) - gox)/h)
      j0 = nint((cirreg(2,ip) - goy)/h)
      k0 = nint((cirreg(3,ip) - goz)/h)

      xyy = cirreg(4,ip)
      xzz = cirreg(5,ip)
      xyz = cirreg(6,ip)

      t(1,1) = cirreg( 7,ip)
      t(1,2) = cirreg( 8,ip)
      t(1,3) = cirreg( 9,ip)
      t(2,1) = cirreg(10,ip)
      t(2,2) = cirreg(11,ip)
      t(2,3) = cirreg(12,ip)
      t(3,1) = cirreg(13,ip)
      t(3,2) = cirreg(14,ip)
      t(3,3) = cirreg(15,ip)

      call twosided(xm,ym,zm,1.0d0,80.0d0,10,i0,j0,k0,nbnd,&
         cirreg(1,ip),cirreg(2,ip),cirreg(3,ip),up,df(1),df(2),df(3),&
         xyy,xzz,xyz,t,ui,u,gox,goy,goz,h,ip)
      df = -df ! this is reaction field in the local coordinate frame

      ! second compute coulomb field

      coulomb(1:3) = ZERO
      do iatm = 1, natom
         dx(1:3) = xi(1:3) - grdcrd(1:3,nbnd+iatm)
         dist0 = ONE/(dx(1)**2 + dx(2)**2 + dx(3)**2); dist = sqrt(dist0)

         dcc = acrg(iatm)*dist0*dist

         coulomb(1) = coulomb(1) + dx(1)*dcc
         coulomb(2) = coulomb(2) + dx(2)*dcc
         coulomb(3) = coulomb(3) + dx(3)*dcc
      end do

      ! Dn and DBF, back to the pbsa internal unit

      dn = (df(1) + dot_product(coulomb(1:3), n(1:3)))/(FOURPI*epsin)
      dbf = HALF*charge*dn
      dum(1:3) = dbf*n(1:3)

! Changhao: comparing with analytical test case
!
!if ( xi(1)>0.0d0 .and. xi(2)>0.0d0 .and. xi(3)>0.0d0 ) then
!   write(311,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)>0.0d0 .and. xi(3)<0.0d0 ) then
!   write(312,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)<0.0d0 .and. xi(3)>0.0d0 ) then
!   write(313,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)<0.0d0 .and. xi(3)<0.0d0 ) then
!   write(314,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)>0.0d0 .and. xi(3)>0.0d0 ) then
!   write(315,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)>0.0d0 .and. xi(3)<0.0d0 ) then
!   write(316,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)<0.0d0 .and. xi(3)>0.0d0 ) then
!   write(317,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)<0.0d0 .and. xi(3)<0.0d0 ) then
!   write(318,*) dum(1),dum(2),dum(3)
!end if
!
!if ( xi(1)>=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(411,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(412,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(413,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(414,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(415,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)>=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(416,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)>=0.0d0 ) then
!   write(417,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<=0.0d0 .and. xi(2)<=0.0d0 .and. xi(3)<=0.0d0 ) then
!   write(418,*) dum(1),dum(2),dum(3)
!end if
!
!if ( xi(1)==0.0d0 .and. xi(2)==0.0d0 .and. xi(3)>0.0d0 ) then
!   write(511,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)==0.0d0 .and. xi(2)==0.0d0 .and. xi(3)<0.0d0 ) then
!   write(512,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)==0.0d0 .and. xi(2)>0.0d0 .and. xi(3)==0.0d0 ) then
!   write(513,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)==0.0d0 .and. xi(2)<0.0d0 .and. xi(3)==0.0d0 ) then
!   write(514,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)>0.0d0 .and. xi(2)==0.0d0 .and. xi(3)==0.0d0 ) then
!   write(515,*) dum(1),dum(2),dum(3)
!end if
!if ( xi(1)<0.0d0 .and. xi(2)==0.0d0 .and. xi(3)==0.0d0 ) then
!   write(516,*) dum(1),dum(2),dum(3)
!end if

      ! partition the force element to atoms

      do jp = atmfirst, atmlast
         iatm = bndatmlst(jp)
         coef = dot_product(du(1:3,iatm), dusum(1:3))*weight2
         pbfrc(1:3,iatm) = pbfrc(1:3,iatm) + coef*dum(1:3)
      end do

   end do

!write(6,'(a,3f15.6)') 'dbfrc', pbfrc(1:3, 1)*FOURPI*eps0*AMBER_ELECTROSTATIC2

end subroutine p3m_dbfrc2side
#endif /*ndef SANDER or LIBPBSA*/
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Add back short-range electrostatic energy and forces
subroutine p3m_direct( nbnd,natom,acrg,eelrf,pbfrc )

   implicit none

   ! Passed variables

   integer nbnd, natom
   _REAL_ acrg(natom)
   _REAL_ eelrf, pbfrc(3,natom)

   ! Local variables

   integer jatm, ip, jp
   integer atmfirst, atmlast
   _REAL_ xi(3), dx(3)
   _REAL_ dirnrg, dirfrc(3,natom)
   _REAL_ charge, dist2, d2inv, df2
   _REAL_ factor

   ! the ncrg = nbnd + natom here due to the p-p addback.

   dirnrg = ZERO
   dirfrc = ZERO
   do ip = 1, ncrg
      xi(1:3) = grdcrd(1:3,ip)
      charge = grdcrg(ip)

      ! loop over all closeby atoms of this charge

      atmfirst = bndatmptr(ip-1) + 1
      atmlast = bndatmptr(ip)
      do jp = atmfirst, atmlast
         jatm = bndatmlst(jp)

         dx(1:3) = xi(1:3) - grdcrd(1:3,nbnd+jatm)
         dist2 = dx(1)**2 + dx(2)**2 + dx(3)**2
         ! this has to go
         !if ( dist2 > cutfd ) cycle

         d2inv = ONE/dist2
         df2 = charge*acrg(jatm)*sqrt(d2inv)
         dirnrg = dirnrg + df2
         df2 = df2*d2inv
         dirfrc(1:3,jatm) = dirfrc(1:3,jatm) - df2*dx(1:3)
      end do  !  jp = atmfirst, atmlast

   end do  !  ip = 1, ncrg

   dirnrg = dirnrg*AMBER_ELECTROSTATIC2*HALF; eelrf = eelrf + dirnrg
   dirfrc = dirfrc*AMBER_ELECTROSTATIC2; pbfrc = pbfrc + dirfrc

end subroutine p3m_direct
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ pairwise nonbonded forces
subroutine p3m_vdw( natom,iprshrt,iar1pb,cn1pb,cn2pb,x,f,enb)

   implicit none

   ! Passed variables

   integer natom, iprshrt(*), iar1pb(4,0:natom)
   _REAL_ cn1pb(*), cn2pb(*)
   _REAL_ x(3,natom)
   _REAL_ f(3,natom)
   _REAL_ enb

   ! Local variables

   integer i, j, jp, ilast, jfirst, jlast
   _REAL_ xi, yi, zi
   _REAL_ dumx, dumy, dumz
   _REAL_ dx, dy, dz, d2inv, r6
   _REAL_ f2, f1, df, fw1, fw2, fw3

   ! initialization

   enb = ZERO
   ilast = natom - 1
   do i = 1, ilast
      xi = x(1,i); yi = x(2,i); zi = x(3,i)
      dumx = ZERO; dumy = ZERO; dumz = ZERO

      jfirst = iar1pb(1, i) + 1; jlast = iar1pb(4, i)
      do jp = jfirst, jlast
         j = iprshrt(jp)
         dx = xi - x(1,j); dy = yi - x(2,j); dz = zi - x(3,j)
         d2inv = ONE/(dx**2+dy**2+dz**2)
         r6 = d2inv**3
         f2 = cn2pb(jp)*r6
         f1 = cn1pb(jp)*(r6*r6)
         enb = enb + (f2-f1)
         df = SIX*((f2-f1)-f1)*d2inv

         fw1 = dx*df; fw2 = dy*df; fw3 = dz*df
         dumx = dumx + fw1
         dumy = dumy + fw2
         dumz = dumz + fw3
         f(1,j) = f(1,j) + fw1
         f(2,j) = f(2,j) + fw2
         f(3,j) = f(3,j) + fw3
      end do

      f(1,i) = f(1,i) - dumx
      f(2,i) = f(2,i) - dumy
      f(3,i) = f(3,i) - dumz
   end do  !  i = 1, ilast

   enb = -enb

end subroutine p3m_vdw

end module pb_p3m
