! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!+ pack all apparent and real charges on the grid into atom-based lists for P3M
subroutine pb_crgview( verbose,eneout,natom,f,epsx,epsy,epsz,insas,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd, dotarc, arcatm, narcdot, ntri, triatm, triopt
   use poisson_boltzmann, only : pos_crg, ipos_crg, surf_crg, crg_num, &
      iepsav, acrd, nbnd, h, gox, goy, goz, xm, ym, zm

   implicit none

#  include "pb_constants.h"

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ f(3,natom)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)
   _REAL_ insas(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)

   ! local variables

   integer, parameter :: n_point = 27
   integer i, j, k, iatm, jatm, matm, natm, iarc, ip
   _REAL_, parameter :: small_dn = 1.d-3 * AMBER_ELECTROSTATIC
   _REAL_ g(3), x(3), dx(3), crd(3), rn(3), sgn, rdx, rsphere, crg
   _REAL_ srfcrg, coulomb(3)
   _REAL_ fx0, fy0, fz0, rh
   _REAL_ up, dudxi0, dudyi0, dudzi0
   _REAL_ dn, dt2, dum(3)
   _REAL_ mvec(3), nvec(3), mdotn, mxnv(3), rmdotn2, fdotm, fdotn
   _REAL_ dfm, dfn, dum_norm(3), dum_tang(3), dumnorm
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)

   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_, allocatable :: pol_charge(:)
   integer   atmflag, crgflag

   integer reentopt, m, n, mp, np, patm, qatm
  ! _REAL_, allocatable :: a(:,:), u(:,:), w(:), v(:,:), b(:), t(:)
   _REAL_ d1, d2, wmax, thresh, TOL
   _REAL_ pvec(3), qvec(3), mdist, ndist, pdist, qdist
   _REAL_ scalar
   mp = 3; np = 3;
   TOL = 1.d-5

   allocate (pol_charge(1:nbnd))

   ! initialization

   rh = ONE/h
   !  a = ZERO; u = ZERO; w = ZERO; v = ZERO; b = ZERO; t = ZERO
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO
   x = ZERO; rsphere = ZERO; sgn = ONE; coulomb = ZERO
   crg_num = ZERO

   ! compute polarization charges on the boundary grid points and
   ! report total in srfcrg
   ! pol_charge is in kcal/mol

   srfcrg = ZERO
   call get_charge_pol( nbnd,xm,ym,zm,phi,cphi,chgrd,pol_charge,srfcrg )

   ! main double loops over polarization charges and atom charges

   do ip = 1, nbnd

      ! collect boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
      crg = pol_charge(ip)

      ! project the boundary grid point onto the molecular surface, crd() is the
      ! projected coord, x() is the atom/probe center coord, and fx/y/z0 is the grid
      ! version of crd()

      if ( iatm == 0 ) then
         write(6,'(a)') 'PB Bomb in pb_crgview(): can not find projection atom/probe'
         call mexit(6, 1)
      end if

      if ( abs(insas(i,j,k)) == TWO ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         sgn = ONE

         atmflag = iatm
         crg_num(atmflag) = crg_num(atmflag) + 1
         crgflag = crg_num(atmflag)
         surf_crg(crgflag,atmflag) = pol_charge(ip)

      else if ( abs(insas(i,j,k)) == ONE ) then
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         sgn = -ONE

         if ( triopt == 1 ) then
            if ( -iatm > narcdot-ntri ) then
               matm = triatm(1,ntri-iatm-narcdot)
               natm = triatm(2,ntri-iatm-narcdot)
               patm = triatm(3,ntri-iatm-narcdot)
            else
               iarc = dotarc(-iatm)
               matm = arcatm(1,iarc)
               natm = arcatm(2,iarc)
            end if
         else
            iarc = dotarc(-iatm)
            matm = arcatm(1,iarc)
            natm = arcatm(2,iarc)
         end if
         atmflag = matm
         crg_num(atmflag) = crg_num(atmflag) + 1
         crgflag = crg_num(atmflag)
         surf_crg(crgflag,atmflag) = pol_charge(ip)

      end if

      if ( crgflag > MAXSURFC ) then
         write(6,'(a,2i6)') 'PB Bomb in pb_crgview: Maximum stored surface charges reached for atom', &
            atmflag, crgflag
         call mexit(6,1)
      end if

      dx = g - x
      rdx = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
      crd = x + dx*rdx*rsphere
!     crd = g ! if using the option of not shifting the grid charges
      pos_crg(1,crgflag,atmflag) = crd(1)
      pos_crg(2,crgflag,atmflag) = crd(2)
      pos_crg(3,crgflag,atmflag) = crd(3)
      ipos_crg(1,crgflag,atmflag) = i
      ipos_crg(2,crgflag,atmflag) = j
      ipos_crg(3,crgflag,atmflag) = k

   end do

   deallocate (pol_charge)

!do i = 1, natom
!do j = 1, crg_num(i)
!write(401,'(4f15.6)')ipos_crg(1:3,j,i), surf_crg(j,i)
!write(401,'(2i10)') i, crg_num(i)
!end do
!end do

end subroutine pb_crgview
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ FD coulombic energy and forces.
subroutine pb_fdreaction( natom, atmfirst, atmlast, idecomp, grdreac, outflag, pbfrc )

   use poisson_boltzmann, only : pos_crg, ipos_crg, surf_crg, crg_num, icrd, gcrg, iar1pb, &
      iprshrt, h, epsin, frcfac

   implicit none

   ! Common variables

   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green

#  include "pb_constants.h"

   ! Passed variables

   integer natom, atmfirst, atmlast, idecomp
   integer outflag(natom)
   _REAL_ grdreac
   _REAL_ pbfrc(3, natom)

   ! Local variables

   integer iatm, jatm, jp, ilast, jfirst, jlast
   integer dijx, dijx0, dijx1, dijx2, dijx3, dijy, dijy0, dijy1, dijy2, dijy3, dijz, dijz0, dijz1, dijz2, dijz3
   integer ix, iy, iz, jx, jy, jz
   integer iix, iiy, iiz, jjx, jjy, jjz, i, j
   _REAL_ gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ gcj1, gcj2, gcj3, gcj4, gcj5, gcj6, gcj7, gcj8
   _REAL_ gcij(27)
   _REAL_ factor, factor1, decfac
   _REAL_ ffx, ffy, ffz
   _REAL_ dumx, dumy, dumz
   _REAL_ frc(3, natom)
   _REAL_ grdreactmp

   integer atmflag, crgflag
   integer isrf, jsrf, icrg, jcrg

   ! begin code

   factor  = ONE/( FOURPI*epsin*h )
   factor1 = HALF/( FOURPI*epsin*h*h )
   decfac  = -factor*frcfac

   grdreac = ZERO
   frc = ZERO

   ! first loop over iatm's grid charges with iatm's surface charges

   do iatm = atmfirst, atmlast
!     if ( ligand .and. realflag(iatm) == 0 ) then
!        cycle
!     else if ( outflag(iatm) == 1 ) then
!        cycle
!     end if
      ix = icrd(1,iatm); iy = icrd(2,iatm); iz = icrd(3,iatm)

      ! loop over iatm's grid charges over iatm's surface charges

      grdreactmp = ZERO
      do isrf = 1, crg_num(iatm)
            jx = ipos_crg(1, isrf, iatm)
            jy = ipos_crg(2, isrf, iatm)
            jz = ipos_crg(3, isrf, iatm)
            dijx  =       jx-ix; dijy  =       jy-iy; dijz  =       jz-iz
            dijx1 = abs(dijx-1); dijy1 = abs(dijy-1); dijz1 = abs(dijz-1)
            dijx  = abs(dijx  ); dijy  = abs(dijy  ); dijz  = abs(dijz  )
            grdreactmp = grdreactmp + &
            surf_crg(isrf,iatm)*( gcrg(1,iatm)*l_green(dijx ,dijy ,dijz ) + gcrg(2,iatm)*l_green(dijx1,dijy ,dijz ) + &
                                  gcrg(3,iatm)*l_green(dijx ,dijy1,dijz ) + gcrg(4,iatm)*l_green(dijx1,dijy1,dijz ) + &
                                  gcrg(5,iatm)*l_green(dijx ,dijy ,dijz1) + gcrg(6,iatm)*l_green(dijx1,dijy ,dijz1) + &
                                  gcrg(7,iatm)*l_green(dijx ,dijy1,dijz1) + gcrg(8,iatm)*l_green(dijx1,dijy1,dijz1) )

      end do
         grdreac = grdreac + grdreactmp

   end do

   ! second loop over iatm's grid charges with iatm's surface charges

   do iatm = atmfirst, atmlast
!     if ( ligand .and. realflag(iatm) == 0 ) then
!        cycle
!     else if ( outflag(iatm) == 1 ) then
!        cycle
!     end if
      ix = icrd(1,iatm); iy = icrd(2,iatm); iz = icrd(3,iatm)

      jfirst = iar1pb(4, iatm-1) + 1
      jlast  = iar1pb(2, iatm)
      dumx = ZERO; dumy = ZERO; dumz = ZERO
      do jp = jfirst, jlast
         jatm = iprshrt(jp)
         if( jatm > atmlast ) cycle
         jx = icrd(1,jatm); jy = icrd(2,jatm); jz = icrd(3,jatm)

         ! loop over jatm's grid charges over iatm's surface charges

         grdreactmp = ZERO
         do isrf = 1, crg_num(iatm)
               iix = ipos_crg(1, isrf, iatm)
               iiy = ipos_crg(2, isrf, iatm)
               iiz = ipos_crg(3, isrf, iatm)
            dijx  =       iix-jx; dijy  =       iiy-jy; dijz  =       iiz-jz
            dijx1 = abs(dijx-1); dijy1 = abs(dijy-1); dijz1 = abs(dijz-1)
            dijx  = abs(dijx  ); dijy  = abs(dijy  ); dijz  = abs(dijz  )
            grdreactmp = grdreactmp + &
            surf_crg(isrf,iatm)*( gcrg(1,jatm)*l_green(dijx ,dijy ,dijz ) + gcrg(2,jatm)*l_green(dijx1,dijy ,dijz ) + &
                                  gcrg(3,jatm)*l_green(dijx ,dijy1,dijz ) + gcrg(4,jatm)*l_green(dijx1,dijy1,dijz ) + &
                                  gcrg(5,jatm)*l_green(dijx ,dijy ,dijz1) + gcrg(6,jatm)*l_green(dijx1,dijy ,dijz1) + &
                                  gcrg(7,jatm)*l_green(dijx ,dijy1,dijz1) + gcrg(8,jatm)*l_green(dijx1,dijy1,dijz1) )

         end do
         grdreac = grdreac + grdreactmp

         ! loop over iatm's grid charges over jatm's surface charges

         grdreactmp = ZERO
         do jsrf = 1, crg_num(jatm)
               jjx = ipos_crg(1, jsrf, jatm)
               jjy = ipos_crg(2, jsrf, jatm)
               jjz = ipos_crg(3, jsrf, jatm)
            dijx  =       jjx-ix; dijy  =       jjy-iy; dijz  =       jjz-iz
            dijx1 = abs(dijx-1); dijy1 = abs(dijy-1); dijz1 = abs(dijz-1)
            dijx  = abs(dijx  ); dijy  = abs(dijy  ); dijz  = abs(dijz  )
            grdreactmp = grdreactmp + &
            surf_crg(jsrf,jatm)*( gcrg(1,iatm)*l_green(dijx ,dijy ,dijz ) + gcrg(2,iatm)*l_green(dijx1,dijy ,dijz ) + &
                                  gcrg(3,iatm)*l_green(dijx ,dijy1,dijz ) + gcrg(4,iatm)*l_green(dijx1,dijy1,dijz ) + &
                                  gcrg(5,iatm)*l_green(dijx ,dijy ,dijz1) + gcrg(6,iatm)*l_green(dijx1,dijy ,dijz1) + &
                                  gcrg(7,iatm)*l_green(dijx ,dijy1,dijz1) + gcrg(8,iatm)*l_green(dijx1,dijy1,dijz1) )

         end do
         grdreac = grdreac + grdreactmp

      end do  !  jp = jfirst, jlast

   end do  !  iatm = 1, ilast

   grdreac = grdreac*AMBER_ELECTROSTATIC*HALF/h
!   print *, grdreac, 'FD reaction'
contains

function l_green (i,j,k)

   implicit none
   integer i,j,k
   _REAL_ l_green

   if ( i <= 40  .and. j <= 40 .and. k <= 40 ) then
      l_green = green(i,j,k)
   else
      l_green = ONE/sqrt(dble(i*i+j*j+k*k))
   end if

end function l_green

end subroutine pb_fdreaction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Direct coulombic energy and forces.
subroutine pb_direct_reaction( natom, atmfirst, atmlast, idecomp, drcreac, outflag, pbfrc )

   use poisson_boltzmann, only : pos_crg, ipos_crg, surf_crg, crg_num, acrd, acrg, iar1pb, &
      iprshrt, h, frcfac, epsin

   implicit none

#  include "pb_constants.h"

   ! Passed variables

   integer natom, atmfirst, atmlast, idecomp
   integer outflag(natom)
   _REAL_ drcreac
   _REAL_ pbfrc(3, natom)

   ! Local variables

   integer iatm, jatm, jp, ilast, jfirst, jlast
   integer dijx, dijx0, dijx1, dijx2, dijx3, dijy, dijy0, dijy1, dijy2, dijy3, dijz, dijz0, dijz1, dijz2, dijz3
   integer ix, iy, iz, jx, jy, jz
   integer iix, iiy, iiz, jjx, jjy, jjz, i, j
   _REAL_ gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ gcj1, gcj2, gcj3, gcj4, gcj5, gcj6, gcj7, gcj8
   _REAL_ gcij(27)
   _REAL_ factor, factor1, decfac
   _REAL_ ffx, ffy, ffz
   _REAL_ dumx, dumy, dumz
   _REAL_ frc(3, natom)
   _REAL_ drcreactmp
   _REAL_ xx_i, yy_i, zz_i, xx_j, yy_j, zz_j, dx(1:3)

   integer atmflag, crgflag, ncount
   integer isrf, jsrf, icrg, jcrg

   ! begin code

   factor  = ONE/( FOURPI*epsin*h )
   factor1 = HALF/( FOURPI*epsin*h*h )
   decfac  = -factor*frcfac

   drcreac = ZERO
   frc = ZERO

   ! first loop over iatm's grid charges with iatm's surface charges

   do iatm = atmfirst, atmlast
!     if ( ligand .and. realflag(iatm) == 0 ) then
!        cycle
!     else if ( outflag(iatm) == 1 ) then
!        cycle
!     end if
      xx_i = acrd(1,iatm); yy_i = acrd(2,iatm); zz_i = acrd(3,iatm)

      ! loop over iatm's grid charges over iatm's surface charges

      drcreactmp = ZERO
      do isrf = 1, crg_num(iatm)
            dx(1) = pos_crg(1, isrf, iatm) - xx_i
            dx(2) = pos_crg(2, isrf, iatm) - yy_i
            dx(3) = pos_crg(3, isrf, iatm) - zz_i
            drcreactmp = drcreactmp + surf_crg(isrf,iatm)*acrg(iatm)/sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
      end do
         drcreac = drcreac + drcreactmp

   end do

   ! second loop over iatm's grid charges with iatm's surface charges

   do iatm = atmfirst, atmlast
!     if ( ligand .and. realflag(iatm) == 0 ) then
!        cycle
!     else if ( outflag(iatm) == 1 ) then
!        cycle
!     end if
      xx_i = acrd(1,iatm); yy_i = acrd(2,iatm); zz_i = acrd(3,iatm)

      jfirst = iar1pb(4, iatm-1) + 1
      jlast  = iar1pb(2, iatm)
      dumx = ZERO; dumy = ZERO; dumz = ZERO
      do jp = jfirst, jlast
         jatm = iprshrt(jp)
         if( jatm > atmlast ) cycle
         xx_j = acrd(1,jatm); yy_j = acrd(2,jatm); zz_j = acrd(3,jatm)

         ! loop over jatm's grid charges over iatm's surface charges

         drcreactmp = ZERO
         do isrf = 1, crg_num(iatm)
               dx(1) = pos_crg(1, isrf, iatm) - xx_j
               dx(2) = pos_crg(2, isrf, iatm) - yy_j
               dx(3) = pos_crg(3, isrf, iatm) - zz_j
            drcreactmp = drcreactmp + surf_crg(isrf,iatm)*acrg(jatm)/sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         end do
         drcreac = drcreac + drcreactmp

         ! loop over iatm's grid charges over jatm's surface charges

         drcreactmp = ZERO
         do jsrf = 1, crg_num(jatm)
               dx(1) = pos_crg(1, jsrf, jatm) - xx_i
               dx(2) = pos_crg(2, jsrf, jatm) - yy_i
               dx(3) = pos_crg(3, jsrf, jatm) - zz_i
            drcreactmp = drcreactmp + surf_crg(jsrf,jatm)*acrg(iatm)/sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         end do
         drcreac = drcreac + drcreactmp

      end do  !  jp = jfirst, jlast

   end do  !  iatm = 1, ilast

   drcreac = drcreac*AMBER_ELECTROSTATIC*HALF
!  print *, drcreac, 'direct reaction'
end subroutine pb_direct_reaction
