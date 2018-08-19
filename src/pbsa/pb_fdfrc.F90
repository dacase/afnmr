! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"
#include "timer.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ driver for FDPB forces and energy
!     call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,inatm,npdec,idecomp,irespw, &
!             ipres,jgroup,ibgwat,ibgion,pbfrc,enb,eelrffd,npbstep,npbgrid,nstlim)
subroutine pb_fdfrc( pbverbose,pbprint,pbgrid,ifcap,ipb,inp,imin,natom,atmlast,npdec,idecomp,irespw, &
              ipres,jgroup,ibgwat,ibgion,pbfrc,enb,eelrf,npbstep,npbgrid,nstlim )

   use poisson_boltzmann
   use solvent_accessibility, only : dprob,iprob,radi,radip3,nzratm, &
       narcdot,maxarc,marc,m2narc,fstarc,arcatm,arccrd,savarc,dotarc
	   
#if !defined SANDER && !defined LIBPBSA
   use iim_util, only : pb_iimdrv
   use ana_iim, only : pb_anaiim
   use iimaug, only : pb_augdrv
#endif /*#if !defined SANDER && !defined LIBPBSA*/
   
   use pb_p3m, only : p3m_init, p3m_fdpotential, p3m_qefrc, p3m_direct, p3m_dbfrc, p3m_vdw
   use pbtimer_module
   implicit none

#  include "flocntrl.h"
#  include "parallel.h"

   ! passed variables

   logical pbverbose, pbprint, pbgrid
   integer ifcap, ipb, inp, imin, natom, atmlast, npdec, idecomp, ibgwat, ibgion
   integer npbstep, npbgrid, nstlim
   integer irespw(*), ipres(*), jgroup(*)
   _REAL_ pbfrc(3,natom)
   _REAL_ enb, eelrf

   ! local variables

   integer atmfirst
   integer iatm, lastatm, mpdec, i
   integer cnt, j, k
   _REAL_ eelself, eelcoul, rh, fcrd(3,atmlast)
   _REAL_ aa, bb, cc, aa1, bb1, cc1
   _REAL_ bb1cc1, bb_cc1, bb1cc, bb_cc
   _REAL_ acrgtmp, factor
   _REAL_ grdreac, drcreac
   _REAL_ ionene, osmene

   ! map output options and files

   character(len=16) phifilename
   character(len=32) phidataname
   integer phifilenum
   integer ii,wloops

   if ( phiform == 0 .OR. phiform == 1 ) then
      phifilename="pbsa_phi.phi"
   else
      phifilename="pbsa_phi.dx "
   end if
   phidataname="Electrostatic Potential"
   phifilenum=64

   if ( do_pbfd == 0 ) return

   ionene = ZERO; osmene = ZERO

   mpdec = 1
   if ( idecomp > 2 ) mpdec = npdec

   !-- PB pairwise decomp <<1,mpdec>>
   ! atmfirst is supposedly to be passed

   atmfirst = 1
   do m = 1, mpdec

      ! do fdpb calculations upto nfocus

      eelself = ZERO; eelcoul = ZERO

      if ( ligand ) then
         icrd = 0
         phi = ZERO; bv = ZERO; xs = ZERO!; fcrd = ZERO; chgrd = ZERO
      end if

      xsoffset = 1
      do level = 1, nfocus

         call pbtimer_start(PBTIME_PBBUILDSYS)

         ! retrieving saved grid data into working variables

         bcopt = savbcopt(level)
         xm = savxm(level); ym = savym(level); zm = savzm(level)
         xmym = xm*ym; xmymzm = xmym*zm
         h = savh(level)
         gox = savgox(level); goy = savgoy(level); goz = savgoz(level)

         ! grid-unit version of the atom data array at the working level

         rh = ONE/h
         if ( ifcap == 1 ) then
            gcrd(1,0) = (acrd(1,0) - gox)*rh
            gcrd(2,0) = (acrd(2,0) - goy)*rh
            gcrd(3,0) = (acrd(3,0) - goz)*rh
         end if
         do iatm = atmfirst, atmlast
            gcrd(1,iatm) = (acrd(1,iatm) - gox)*rh
            gcrd(2,iatm) = (acrd(2,iatm) - goy)*rh
            gcrd(3,iatm) = (acrd(3,iatm) - goz)*rh
         end do

         do iatm = atmfirst, atmlast
            if ( ligand ) then
               if ( level > 1 .and. realflag(iatm) == 0 ) cycle
            endif
            !if ( level > 1 .and. ligand     .and. realflag(iatm) == 0 ) cycle
            icrd(1,iatm) = floor(gcrd(1,iatm))
            icrd(2,iatm) = floor(gcrd(2,iatm))
            icrd(3,iatm) = floor(gcrd(3,iatm))
            fcrd(1,iatm) = dble(icrd(1,iatm))
            fcrd(2,iatm) = dble(icrd(2,iatm))
            fcrd(3,iatm) = dble(icrd(3,iatm))
            if ( icrd(1,iatm) > xm-1 .or. icrd(2,iatm) > ym-1 .or. icrd(3,iatm) > zm-1 ) then
               write(6,'(3f12.4)') acrd(1:3,iatm)
               write(6,'(a,3i6)') "pb_fdfrc(): Atom out of focusing box",icrd(1:3,iatm)
               call mexit(6,1)
            end if
         end do

         do iatm = atmfirst, atmlast
            if ( ligand ) then
               if ( level > 1 .and. realflag(iatm) == 0 ) cycle
            end if
            !if ( level > 1 .and. ligand     .and. realflag(iatm) == 0 ) cycle
            aa = gcrd(1,iatm) - fcrd(1,iatm)
            bb = gcrd(2,iatm) - fcrd(2,iatm)
            cc = gcrd(3,iatm) - fcrd(3,iatm)
            bb1 = ONE - bb; cc1 = ONE - cc

            !-- PB decomp

            if (idecomp < 3 ) then
               acrgtmp = acrg(iatm)
            else if (iatm >= ipres(irespw(m)) .and. iatm < ipres(irespw(m)+1)) then
               acrgtmp = acrg(iatm)
            else
               acrgtmp = ZERO
            end if

            aa  = acrgtmp*aa; aa1 = acrgtmp - aa
            bb1cc1 = bb1*cc1; bb_cc1 = bb *cc1
            bb1cc  = bb1*cc ; bb_cc  = bb *cc
            if ( (ifcap == 2 .or. ifcap == 5) .and. outflag(iatm) == 1 ) then
               gcrg(1,iatm) = ZERO; gcrg(2,iatm) = ZERO
               gcrg(3,iatm) = ZERO; gcrg(4,iatm) = ZERO
               gcrg(5,iatm) = ZERO; gcrg(6,iatm) = ZERO
               gcrg(7,iatm) = ZERO; gcrg(8,iatm) = ZERO
            else
               gcrg(1,iatm) = aa1*bb1cc1; gcrg(2,iatm) = aa *bb1cc1
               gcrg(3,iatm) = aa1*bb_cc1; gcrg(4,iatm) = aa *bb_cc1
               gcrg(5,iatm) = aa1*bb1cc ; gcrg(6,iatm) = aa *bb1cc
               gcrg(7,iatm) = aa1*bb_cc ; gcrg(8,iatm) = aa *bb_cc
            end if
         end do

         ! part I. set up dielectric map
         ! when ifcap == 0, do dielectric map assignment everystep
         ! when ifcap /= 0, do dielectric map assignment once only when grid is set up

         if ( ifcap /= 1 ) then

            ! part I.a.
            ! for membrane simulations, set up membrane zmin and zmax

            if ( membraneopt > 0 ) then
               if ( abs(mctrdz) < 0.0001d0 ) then
               ! The default option is just use the box center as the membrane
               ! center
                  mzmin = (czbox(1) + mctrdz - mthick*HALF - goz)/h
                  mzmax = (czbox(1) + mctrdz + mthick*HALF - goz)/h
               else
               ! Otherwise, use the read in center as the membrane center
                  mzmin = mctrdz - mthick*HALF
                  mzmax = mctrdz + mthick*HALF
                  if ( mzmin < zmin .or. mzmax > zmax ) then
                     write(6,'(a)') 'PBSA BOMB: in pb_fdfrc()'
                     write(6,'(a)') 'Membrane slab outside protein: '
                     write(6,'(a,2f10.3)') 'Membrane mzmin/mzmax', mzmin, mzmax
                     write(6,'(a,2f10.3)') 'Protein zmin/zmax', zmin, zmax
                     call mexit(6,1)
                  end if
                  mzmin = (mzmin - goz)/h
                  mzmax = (mzmax - goz)/h
               end if
            end if

            ! part I.b.
            ! when solving systems with salt, set up stern layer map, on the
            ! grid points

            if ( istrng /= ZERO ) call pb_ionmap( pbverbose,ifcap,natom,&
               membraneopt,mzmin,mzmax,iprob,h,gox,goy,goz,xm,ym,zm,xmymzm,&
               outflag,gcrd,radi,atmsas,insas,saltgrd)

            ! part I.c.
            ! here comes the dielectric map on the grid edges: x, y, and z

            call pb_exmol( pbverbose,ifcap,mpdec,ipb,inp,savbcopt,saopt,sasopt,natom,&
               smoothopt,epsin,epsout,epsmem,membraneopt,mzmin,mzmax,&
               h,gox,goy,goz,xm,ym,zm,xmymzm,&
               level,nfocus,&
               nbnd,nbndx,nbndy,nbndz,&
               outflag,gcrd,acrd,&
               atmsas,insas,lvlset,mlvlset,zv,epsx,epsy,epsz,&
               iepsav,iepsavx,iepsavy,iepsavz,fedgex,fedgey,fedgez)
         elseif ( pbgrid ) then
            call pb_exmol_cap( pbverbose,ifcap )
         end if

         ! part II. set up grid charges and store in chgrd, which is in the unit
         ! of eletron

         chgrd = ZERO
         call pb_crggrd( natom,1,atmlast,xm,ym,zm,icrd,gcrg,chgrd(1) )

         call pbtimer_stop(PBTIME_PBBUILDSYS)

         ! for surface area calculations ... we are done

         if ( saopt < 0 ) then
            if ( level == nfocus ) then
               return
            else
               cycle
            end if
         end if

         ! part III, call linear system solution drivers

         if ( ipb >= 1 .and. ipb <=3 ) then
            if ( solvopt /= 7 ) then
               call pb_fddrv( npbstep,npbgrid,nstlim,1,natom,npbopt,solvopt,level,nfocus,bcopt,&
                  h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                  xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                  maxitn,itn,fmiccg,fmiccg2,accept,laccept,wsor,lwsor,inorm,norm,&
                  pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,osmene,&
                  ngrdcrg,grdcrg,qgrdcrg,&
                  gcrd,acrg,&
                  nbnd,iepsav,insas,epsx,epsy,epsz,&
                  chgrd,saltgrd,ioncrg,phi,&
                  bv,cphi,xs(xsoffset) )
            end if
            if ( solvopt == 7 ) then
               write(6,'(a)') "PB_FFT is no longer supported in any release of PBSA"
               call mexit(6,1)
            end if
         else if ( ipb == 4 ) then
            if ( ligand ) then
               write(6,'(a)') "IIM doesn't work with ligandmask."; call mexit(6,1)
            end if
#if !defined SANDER && !defined LIBPBSA
            call pb_iimdrv( npbstep,npbgrid,nstlim,1,natom,npbopt,solvopt,level,nfocus,bcopt,&
               natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
               xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
               maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
               pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
               gcrd,acrg,&
               nbnd,iepsav,insas,lvlset,&
               chgrd,saltgrd,phi,&
               bv,cphi,xs(xsoffset) )
#else
            write(6,'(a)') "PB_IIM is not supported in SANDER or LIBPBSA";call mexit(6,1)
#endif /*#if !defined SANDER && !defined LIBPBSA*/
         else if ( ipb == 5 ) then
            if ( ligand ) then
               write(6,'(a)') "AUG doesn't work with ligandmask."; call mexit(6,1)
            end if
#if !defined SANDER && !defined LIBPBSA
            call pb_augdrv( npbstep,npbgrid,nstlim,1,natom,npbopt,solvopt,level,nfocus,bcopt,&
               natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
               xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
               maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
               pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
               gcrd,acrg,&
               nbnd,iepsav,insas,lvlset,&
               chgrd,saltgrd,phi,&
               bv,cphi,xs(xsoffset),&
               pbverbose, &
               augtoltype, augctf, augtol )
#else
            write(6,'(a)') "PB_AUG is not supported in SANDER or LIBPBSA";call mexit(6,1)
#endif /*#if !defined SANDER && !defined LIBPBSA*/
         else if ( ipb == 6 ) then
            if ( ligand ) then
               write(6,'(a)') "ANA-IIM doesn't work with ligandmask."; call mexit(6,1)
            end if
#if !defined SANDER && !defined LIBPBSA
            call pb_anaiim( npbstep,npbgrid,nstlim,1,natom,npbopt,solvopt,level,nfocus,bcopt,&
               natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
               xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
               maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
               pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
               gcrd,acrg,&
               nbnd,iepsav,insas,lvlset,&
               chgrd,saltgrd,phi,&
               bv,cphi,xs(xsoffset) )
#else
            write(6,'(a)') "ANAIIM is not supported in SANDER or LIBPBSA";call mexit(6,1)
#endif /*#if !defined SANDER && !defined LIBPBSA*/
         else
            write(6,'(a)') "PBSA Bomb: pb_fdfrc: unknown ipb";call mexit(6,1)
         end if

         ! part IV,
         ! do some output if requested
         ! first, print a summary when the grid is set up

         if ( pbverbose .and. pbprint ) then
            call pb_print( ifcap, ipb, natom )
            write(6, '(a,I6     )') '   Iterations required        :', itn
            write(6, '(a,F21.10 )') '   Norm of the constant vector:', inorm
            write(6, '(a,F21.13 )') '   Norm of the residual vector:', norm
            write(6, '(a,ES24.16)') '   Convergence achieved       :', norm/inorm
            write(6, '()')
         end if  !  pbverbose .and. pbprint

         ! second, output phi map

         if ( outphi .and. level == nfocus ) then
            if ( phiform == 0 ) then ! write delphi format phi
               write(6,'(a)') 'writing potential map in delphi format'
               open(64,file=phifilename,form="unformatted")
               write(64) ' start of phimap    '
               write(64) ' potential', ' ------ AMBER PBSA ------ phimap in kT/e (0.593kcal/mol-e)  '
               write(64) real((frcfac/0.593d0)*phi(1:xmymzm))
               write(64) ' end of phimap  '
               write(64) real(1.0d0/h),real(cxbox(level)),real(cybox(level)),real(czbox(level)),xm
               close(64)
            else if ( phiform == 1 ) then ! write amber format phi
               write(6,'(a)') 'writing potential map in amber format'
               open(64,file=phifilename,form="formatted")
               write(64,'(a)') '# the following data is provided:'
               write(64,'(a)') '# h, gox, goy, goz'
               write(64,'(a)') '# xm, ym, zm'
               write(64,'(a)') '# phi(1:xmymzm) in kcal/mol-e'
               write(64,'(a)') '# mapping between (i,j,k) and phi index:'
               write(64,'(a)') '# i + xm * ( j-1 + ym * ( k-1 ) )'
               write(64,'(a)') '# grid coordinates: xg = gox + h*i; '
               write(64,'(a)') '# yg = goy + h*j; zg = goz + h*k'
               write(64,'(ES16.8E2 ,ES16.8E2, ES16.8E2, ES16.8E2)') h, gox, goy, goz
               write(64,'(I8 , I8, I8)') xm, ym, zm
               wloops = floor(xmymzm/6d0)
               do ii = 0, wloops-1
                  write(64,FMT='(ES18.8E3, ES18.8E3, ES18.8E3, ES18.8E3, ES18.8E3, ES18.8E3)') frcfac*phi((ii*6+1):(ii*6+6))
               end do
               write(64,FMT='(ES18.8E3, ES18.8E3, ES18.8E3, ES18.8E3, ES18.8E3, ES18.8E3)') frcfac*phi((6*wloops+1):xmymzm)
               close(64)
            else if ( phiform == 2 ) then ! write dx format phi
               write(6,'(a)') 'writing potential map in dx format'
               call gen_dx_file(xm,ym,zm,h,gox,goy,goz,frcfac*phi(1:xmymzm),phifilename,phifilenum,phidataname)
            else
               write(6,'(a)') 'PBSA Warning, unrecognizable phimap output format.'
               write(6,'(a)') 'No phimap will be written out.'
            end if
         end if  ! outphi .and. level == nfocus

         xsoffset = xsoffset + savxmymzm(level) + 2*savxmym(level)
      end do  !  level = 1, nfocus

      ! part V, finally we are ready to interpolate energy and forces.

      ! option 1: compute fd energy and force by the qE option
      ! note that self forces are zero
      ! delete fd coulomb energy and forces for all close pairs
      ! dbf is computed by Gilson et al

      if ( eneopt == 1 ) then

         ! compute total qE energy and forces and delete self energy, note that self forces are zero

         if ( intopt == 1 ) call pb_qefrc( natom,1,atmlast,idecomp,xm,ym,zm,eelrf,eelself,pbfrc,phi )
         zv = -dble(insas) ! pseudo signed distance function
         if ( intopt == 2 ) then
            if ( bcopt < 6 .or. bcopt > 9 ) then
               write(6,'(a)') ' PBSA Bomb in pb_force(): intop = 2 does not support bcopt < 6'
               call mexit(6, 1)
            end if
            call pb_qefrc2( natom,1,atmlast,xm,ym,zm,eelrf,pbfrc,zv(1),phi )
         end if

         if ( pbverbose .and. pbprint .and. npbopt > 0 ) then
            write(6,'(a,e20.9)') ' full nonlinear PB FDPB energy in kcal/mol', frcfac*(eelrf+ionene)
         end if

         ! delete fd grid energy and forces for all close pairs
         ! when bcopt >= 6, we only have reaction field energy in eelrf
         ! no need to any corrections, coulombic energy should come from
         ! pb_direct ...

         if ( bcopt == 1 .or. bcopt == 4 .or. bcopt == 5 .or. bcopt == 10 ) then
            call pb_fdcoulomb( natom, 1, atmlast, idecomp, eelcoul, outflag, pbfrc )
            eelrf = eelrf - eelself - eelcoul
         end if

         if ( pbverbose .and. pbprint ) then
            write(6,'(a,2e15.4)') 'final eelrf/ionene in kcal/mol', frcfac*eelrf, frcfac*ionene
            write(6,'(a,e15.4)') 'final eelself in kcal/mol', frcfac*eelself
            write(6,'(a,e15.4)') 'final eelcoul in kcal/mol', frcfac*eelcoul
         end if

         ! add ion contributions for nonlinear PBE ...

         if ( npbopt /= 0 ) then
            eelrf = eelrf + ionene
         end if

         ! returning total ES energy after converting to kcal/mol

         eelrf = frcfac*eelrf

         ! the final piece of force calculation, dbf

         if ( frcopt == 1 ) then ! only force option supported for eneopt==1
            if ( pbverbose .and. pbprint ) then ! force dumping statements
               open (unit = 102, file = 'force.dat')
               write(102,'(a)') ' :::: Atomic QE forces ::::'
               do iatm = 1, natom
                  write(102,'(3e20.6)') pbfrc(1:3,iatm)*frcfac
               end do
            end if
            ! Atomic DB forces dumping is inside
            call pb_dbfrc_fld(pbverbose,pbprint,natom,pbfrc,epsx,epsy,epsz,phi,cphi)
            pbfrc = frcfac*pbfrc
         else if ( frcopt == 0 ) then
            pbfrc = ZERO
         else
            write(6,'(a)') 'PBSA Bomb in pb_force(): unsupported frcopt for eneopt=1'
            call mexit(6, 1)
         end if

      ! option 2: compute fdfrc by the charge option
      ! dbf is computed by Cai and Ye et al

      else if ( eneopt == 2 .and. epsin /= epsout ) then
         if (ifcap == 2 .or. ifcap == 5) then ! Only consider protein atoms for calc of total qE energy
            if (ibgion /= 0) then
               lastatm = ipres(ibgion) - 1
            else if (ibgwat /= 0) then
               lastatm = ipres(ibgwat) - 1
            else
               lastatm = atmlast
            end if
         else ! Consider all atoms
            lastatm = atmlast
         end if

         ! the following routines compute both eelrf and forces

         if ( frcopt == 2 ) then
            zv = -dble(insas) ! pseudo signed distance function
            call pb_dbfrc_crg ( pbverbose,pbprint,natom,eelrf,pbfrc, &
               epsx,epsy,epsz,zv,phi,chgrd,cphi )
         else if ( frcopt == 3 ) then
            zv = -dble(insas) ! pseudo signed distance function
            call pb_dbfrc_fld2( pbverbose,pbprint,natom,eelrf,pbfrc, &
               epsx,epsy,epsz,zv,phi,chgrd,cphi )
            !call pb_dbfrc_fld3( pbverbose,pbprint,natom,eelrf,pbfrc, &
            !   epsx,epsy,epsz,zv,phi,chgrd,cphi )
         else if ( frcopt == 4 ) then
            zv = -dble(insas) ! pseudo signed distance function
            call pb_dbfrc_crg2( pbverbose,pbprint,natom,eelrf,pbfrc, &
               epsx,epsy,epsz,zv,phi,chgrd,cphi )
         else ! the default option is not to compute forces
            if ( npbopt == 0 ) then
            call pb_lpbene( pbverbose,pbprint,natom,lastatm,ifcap,&
               npdec,idecomp,m,irespw,ipres,eelrf,insas,phi,chgrd,cphi )
            else
            call pb_npbene(nbnd,xm,ym,zm,chgrd,insas,phi,ioncrg,osmene)
            end if
         end if

      ! option 3: compute fd energy and force by the P3M option
      ! note that self forces are zero
      ! delete fd coulomb energy and forces for all close pairs
      ! delete fd reaction energy and forces for all close pairs

      else if ( eneopt == 3 ) then
         if (ifcap == 2 .or. ifcap == 5) then
            ! Only consider protein atoms for calc of total qE energy
            if (ibgion /= 0) then
               lastatm = ipres(ibgion) - 1
            else if (ibgwat /= 0) then
               lastatm = ipres(ibgwat) - 1
            else
               lastatm = atmlast
            end if
         else
            ! Consider all atoms
            lastatm = atmlast
         end if

         ! compute total qE energy and forces and delete self energy, note that self forces are zero

         if ( intopt == 1 ) call pb_qefrc( natom,1,atmlast,idecomp,xm,ym,zm,eelrf,eelself,pbfrc,phi )
         zv = -dble(insas) ! pseudo signed distance function
         if ( intopt == 2 ) then
            if ( bcopt < 6 .or. bcopt > 9 ) then
               write(6,'(a)') 'PBSA Bomb in pb_force(): intop = 2 does not support bcopt < 6'
               call mexit(6, 1)
            end if
            call pb_qefrc2( natom,1,atmlast,xm,ym,zm,eelrf,pbfrc,zv(1),phi )
         end if

         ! delete fd coulomb energy and forces for all close pairs
         ! when bcopt >= 6, we only have reaction field energy in eelrf
         ! no need to any corrections, coulombic energy should come from
         ! pb_direct ...

         if ( bcopt < 6 .and. bcopt /= 2 ) then
            call pb_fdcoulomb( natom, 1, atmlast, idecomp, eelcoul, outflag, pbfrc )
            eelrf = eelrf - eelself - eelcoul
         end if

         grdreac = ZERO
         drcreac = ZERO
         allocate (pos_crg(3,MAXSURFC,natom))
         allocate (ipos_crg(3,MAXSURFC,natom))
         allocate (surf_crg(MAXSURFC,natom))
         allocate (crg_num(1:natom))
         if ( epsin /= epsout ) then
            call pb_crgview( pbverbose,pbprint,natom,pbfrc,epsx,epsy,epsz,zv(1),phi,chgrd,cphi )
            call pb_fdreaction( natom, 1, atmlast, idecomp, grdreac, outflag, pbfrc )
            call pb_direct_reaction( natom, 1, atmlast, idecomp, drcreac, outflag, pbfrc )
         end if
         deallocate (pos_crg)
         deallocate (ipos_crg)
         deallocate (surf_crg)
         deallocate (crg_num)

         ! add ion contributions for nonlinear PBE ...

         if ( npbopt /= 0 ) then
            eelrf = eelrf + ionene
         end if

         eelrf = frcfac*eelrf
         if ( pbverbose .and. pbprint ) write(6,'(a,f20.10)') 'Standard reaction field energy:', eelrf
         eelrf = eelrf - grdreac + drcreac
         if ( pbverbose .and. pbprint ) write(6,'(a,f20.10)') 'P3M reaction field energy', eelrf

         if ( frcopt == 2 ) then
            zv = -dble(insas) ! pseudo signed distance function
            call pb_dbfrc_crg ( pbverbose,pbprint,natom,eelrf,pbfrc, &
               epsx,epsy,epsz,zv(1),phi,chgrd,cphi )
         else
            write(6,'(a)') 'PBSA Bomb in pb_force(): unsupported frcopt for eneopt=3'
            call mexit(6, 1)
         end if

      ! option 4: compute fd energy and force by the new P3M option (aiming for CUDA)

      else if ( eneopt == 4 ) then

         ! part a: PM calculations/corrections

         call p3m_init( 1,nbnd,natom,xm,ym,zm,iepsav,acrd,acrg,icrd,gcrg,&
            iar1pb,iprshrt,phi,cphi,chgrd,eps0/epsin )
         call p3m_fdpotential( nbnd,natom,iepsav,icrd,h,eps0 )
         call p3m_qefrc( natom,xm,ym,zm,icrd,acrd,acrg,eelrf,pbfrc,phi )

         !         Add FDPB ionene from nonlinear PBE

         if ( npbopt /= 0 ) then
            eelrf = eelrf + ionene
         end if

         ! Need to update the list

         call p3m_init( 2,nbnd,natom,xm,ym,zm,iepsav,acrd,acrg,icrd,gcrg,&
            iar1pb,iprshrt,phi,cphi,chgrd,eps0/epsin )

         ! part b: DBF calculations, there is no PP correction for this term

         call p3m_dbfrc( natom,nbnd,xm,ym,zm,h,dprob,radip3,pbfrc,iepsav,&
            lvlset,phi,eps0,epsx,epsy,epsz )
         !call p3m_dbfrc1d( natom,nbnd,xm,ym,zm,h,gox,goy,goz,dprob,radi,acrg,pbfrc,iepsav,&
         !  lvlset,phi,epsin,frcfac )
         !call p3m_dbfrc1side( natom,nbnd,xm,ym,zm,h,gox,goy,goz,dprob,radi,acrg,pbfrc,iepsav,&
         !   lvlset,phi,epsin,frcfac )
         !call p3m_dbfrc2side( natom,nbnd,xm,ym,zm,h,gox,goy,goz,dprob,radi,acrg,pbfrc,iepsav,&
         !  lvlset,phi,epsin,eps0,frcfac )

         ! Now converting to Amber unit, i.e. kcal/mol-A.

         eelrf = eelrf*frcfac
         pbfrc = pbfrc*frcfac

         ! The following are in the Amber unit

         ! part c: PP calculations, direct electrostatics interactions

         call p3m_direct( nbnd,natom,acrg,eelrf,pbfrc )

         ! part d: PP calcualtions, direct VDW interactions

         call p3m_vdw( natom,iprshrt,iar1pb,cn1pb,cn2pb,acrd(1,1),pbfrc,enb )

         if ( pbverbose .and. pbprint ) then ! dumping forces here
            open (unit = 102, file = 'force.dat')
            write(102,'(a)') ' :::: Atomic total forces ::::'
            do iatm = 1, natom
               write(102,'(3e20.6)') pbfrc(1:3,iatm)
            end do
            write(102,'(a)') ' :::: Average net forces ::::'
            do i=1,3
               write(102,'(3e20.6)') sum(pbfrc(i,1:natom))/dble(natom)
            enddo
            close (102)
         end if

      ! option 5: probably a uniform dielectric run for reference state or development

      else
         eelrf = ZERO
         pbfrc = ZERO
      end if

   end do !-- PB pairwise decomp <<1,mpdec>>

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ total finite difference es energy and forces for intopt = 1
subroutine pb_qefrc( natom,atmfirst,atmlast,idecomp,xm,ym,zm,grdnrg,grdself,pbfrc,phi )

   use decomp, only: decpair
   implicit none

   ! Passed variables

   integer natom, atmfirst, atmlast, idecomp
   integer xm, ym, zm
   _REAL_ grdnrg, grdself
   _REAL_ pbfrc(3,natom)
   _REAL_ phi(xm,ym,zm)

   ! Local variables

   integer iatm
   integer i, j, k
   _REAL_ :: g000, g100, g110, g111
   _REAL_ :: gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ :: grdnrgtmp, grdselftmp, hpsnrv

   ! begin code

   g000 = INV_FOURPI*l_green(0,0,0)
   g100 = INV_FOURPI*l_green(1,0,0)
   g110 = INV_FOURPI*l_green(1,1,0)
   g111 = INV_FOURPI*l_green(1,1,1)

   hpsnrv = ONE/(h*epsin)

   grdnrg = ZERO
   grdself = ZERO

   ! split each atoms charge over the eight surrounding
   ! grid points according to the trilinear weighting
   ! function and add up each of the contributions.

   do iatm = atmfirst, atmlast
      if ( ligand ) then
         if ( liveflag(iatm) == 0 ) cycle
      else if ( outflag(iatm) == 1 ) then
         cycle
      end if
      i = icrd(1,iatm); j = icrd(2,iatm); k = icrd(3,iatm)
      gci1 = gcrg(1,iatm); gci2 = gcrg(2,iatm)
      gci3 = gcrg(3,iatm); gci4 = gcrg(4,iatm)
      gci5 = gcrg(5,iatm); gci6 = gcrg(6,iatm)
      gci7 = gcrg(7,iatm); gci8 = gcrg(8,iatm)

      grdnrgtmp = &
         gci1*phi(i  ,j  ,k  ) + gci2*phi(i+1,j  ,k  ) + &
         gci3*phi(i  ,j+1,k  ) + gci4*phi(i+1,j+1,k  ) + &
         gci5*phi(i  ,j  ,k+1) + gci6*phi(i+1,j  ,k+1) + &
         gci7*phi(i  ,j+1,k+1) + gci8*phi(i+1,j+1,k+1)

      grdnrg = grdnrg + grdnrgtmp

      grdselftmp = &
         g000 * (gci1*gci1 + gci2*gci2 + gci3*gci3 + gci4*gci4 + &
                 gci5*gci5 + gci6*gci6 + gci7*gci7 + gci8*gci8 )*HALF + &
         g100 * (gci1*gci2 + gci1*gci3 + gci1*gci5 + gci2*gci4 + &
                 gci2*gci6 + gci4*gci3 + gci4*gci8 + gci3*gci7 + &
                 gci5*gci6 + gci5*gci7 + gci6*gci8 + gci8*gci7 ) + &
         g110 * (gci1*gci4 + gci1*gci6 + gci1*gci7 + gci2*gci3 + &
                 gci2*gci5 + gci2*gci8 + gci4*gci6 + gci4*gci7 + &
                 gci3*gci5 + gci3*gci8 + gci5*gci8 + gci6*gci7 ) + &
         g111 * (gci1*gci8 + gci2*gci7 + gci4*gci5 + gci6*gci3 )

      grdself = grdself + grdselftmp

      !-- PB decomp
      if(idecomp == 1 .or. idecomp == 2) then
         grdnrgtmp = frcfac*(HALF*grdnrgtmp - grdselftmp*hpsnrv)
         call decpair(1,iatm,iatm,grdnrgtmp)
      end if

      if (frcopt /= 0) then
         pbfrc(1,iatm) = &
         gci1*ex(i  ,j  ,k  ,xm,ym,zm,h,phi) + gci2*ex(i+1,j  ,k  ,xm,ym,zm,h,phi) + &
         gci3*ex(i  ,j+1,k  ,xm,ym,zm,h,phi) + gci4*ex(i+1,j+1,k  ,xm,ym,zm,h,phi) + &
         gci5*ex(i  ,j  ,k+1,xm,ym,zm,h,phi) + gci6*ex(i+1,j  ,k+1,xm,ym,zm,h,phi) + &
         gci7*ex(i  ,j+1,k+1,xm,ym,zm,h,phi) + gci8*ex(i+1,j+1,k+1,xm,ym,zm,h,phi)
         pbfrc(2,iatm) = &
         gci1*ey(i  ,j  ,k  ,xm,ym,zm,h,phi) + gci2*ey(i+1,j  ,k  ,xm,ym,zm,h,phi) + &
         gci3*ey(i  ,j+1,k  ,xm,ym,zm,h,phi) + gci4*ey(i+1,j+1,k  ,xm,ym,zm,h,phi) + &
         gci5*ey(i  ,j  ,k+1,xm,ym,zm,h,phi) + gci6*ey(i+1,j  ,k+1,xm,ym,zm,h,phi) + &
         gci7*ey(i  ,j+1,k+1,xm,ym,zm,h,phi) + gci8*ey(i+1,j+1,k+1,xm,ym,zm,h,phi)
         pbfrc(3,iatm) = &
         gci1*ez(i  ,j  ,k  ,xm,ym,zm,h,phi) + gci2*ez(i+1,j  ,k  ,xm,ym,zm,h,phi) + &
         gci3*ez(i  ,j+1,k  ,xm,ym,zm,h,phi) + gci4*ez(i+1,j+1,k  ,xm,ym,zm,h,phi) + &
         gci5*ez(i  ,j  ,k+1,xm,ym,zm,h,phi) + gci6*ez(i+1,j  ,k+1,xm,ym,zm,h,phi) + &
         gci7*ez(i  ,j+1,k+1,xm,ym,zm,h,phi) + gci8*ez(i+1,j+1,k+1,xm,ym,zm,h,phi)
      end if
   end do

   grdnrg  = HALF*grdnrg
   grdself = grdself/( h*epsin )

end subroutine pb_qefrc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ e field vector x
_REAL_ function ex(i,j,k,xm,ym,zm,h,phi)

   implicit none
   integer, intent(in) :: i, j, k
   integer, intent(in) :: xm, ym, zm
   _REAL_, intent(in) :: h, phi(xm,ym,zm)

   ex = ( phi(i-1,j  ,k  ) - phi(i+1,j  ,k  ) )/(2*h)

end function ex
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ e field vector y
_REAL_ function ey(i,j,k,xm,ym,zm,h,phi)

   implicit none
   integer, intent(in) :: i, j, k
   integer, intent(in) :: xm, ym, zm
   _REAL_, intent(in) :: h, phi(xm,ym,zm)

   ey = ( phi(i  ,j-1,k  ) - phi(i  ,j+1,k  ) )/(2*h)

end function ey
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ e field vector z
_REAL_ function ez(i,j,k,xm,ym,zm,h,phi)

   implicit none
   integer, intent(in) :: i, j, k
   integer, intent(in) :: xm, ym, zm
   _REAL_, intent(in) :: h, phi(xm,ym,zm)

   ez = ( phi(i  ,j  ,k-1) - phi(i  ,j  ,k+1) )/(2*h)

end function ez
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ total finite difference es energy and forces for intopt = 2
subroutine pb_qefrc2( natom,atmfirst,atmlast,xm,ym,zm,grdnrg,pbfrc,insas,phi )

   implicit none

   ! Passed variables

   integer natom, atmfirst, atmlast
   integer xm, ym, zm
   _REAL_ grdnrg
   _REAL_ pbfrc(3, natom)
   _REAL_ insas(0:xm+1, 0:ym+1, 0:zm+1), phi(xm,ym,zm)

   ! Local variables

   integer iatm
   integer i, j, k
   _REAL_ :: gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ :: fx0, fy0, fz0, charge
   _REAL_ :: up, dudxi0, dudyi0, dudzi0
   integer, parameter :: n_point = 8

   ! begin code

   grdnrg = ZERO

   ! split each atoms charge over the eight surrounding
   ! grid points according to the trilinear weighting
   ! function and add up each of the contributions.

   do iatm = atmfirst, atmlast
      if ( ligand ) then
         if ( liveflag(iatm) == 0 ) cycle
      else if ( outflag(iatm) == 1 ) then
         cycle
      end if
      i = icrd(1,iatm); j = icrd(2,iatm); k = icrd(3,iatm)
      fx0 = gcrd(1,iatm); fy0 = gcrd(2,iatm); fz0 = gcrd(3,iatm)
      charge = acrg(iatm)

      call onesided(xm,ym,zm,4,n_point,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,insas,h)

      grdnrg = grdnrg + charge*up

      pbfrc(1,iatm) = pbfrc(1,iatm) - charge*dudxi0
      pbfrc(2,iatm) = pbfrc(2,iatm) - charge*dudyi0
      pbfrc(3,iatm) = pbfrc(3,iatm) - charge*dudzi0
   end do

   grdnrg = HALF*grdnrg

end subroutine pb_qefrc2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ FD coulombic energy and forces.
subroutine pb_fdcoulomb( natom, atmfirst, atmlast, idecomp, grdcoul, outflag, pbfrc )

   use decomp, only: decpair
   implicit none

   ! Passed variables

   integer natom, atmfirst, atmlast, idecomp
   integer outflag(natom)
   _REAL_ grdcoul
   _REAL_ pbfrc(3, natom)

   ! Local variables

   integer iatm, jatm, jp, ilast, jfirst, jlast
   integer dijx, dijx0, dijx1, dijx2, dijx3, dijy, dijy0, dijy1, dijy2, dijy3, dijz, dijz0, dijz1, dijz2, dijz3
   integer ix, iy, iz, jx, jy, jz
   _REAL_ gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ gcj1, gcj2, gcj3, gcj4, gcj5, gcj6, gcj7, gcj8
   _REAL_ gcij(27)
   _REAL_ factor, factor1, decfac
   _REAL_ ffx, ffy, ffz
   _REAL_ dumx, dumy, dumz
   _REAL_ frc(3, natom)
   _REAL_ grdcoultmp
   _REAL_ pair_correct

   ! begin code

   factor  = ONE/( FOURPI*epsin*h )
   factor1 = HALF/( FOURPI*epsin*h*h )
   decfac  = -factor*frcfac

   grdcoul = ZERO
   frc = ZERO
   do iatm = atmfirst, atmlast
      if ( ligand ) then
         if ( liveflag(iatm) == 0 ) cycle
      else if ( outflag(iatm) == 1 ) then
         cycle
      end if
      ix = icrd(1,iatm); iy = icrd(2,iatm); iz = icrd(3,iatm)

      gci1 = gcrg(1,iatm); gci2 = gcrg(2,iatm)
      gci3 = gcrg(3,iatm); gci4 = gcrg(4,iatm)
      gci5 = gcrg(5,iatm); gci6 = gcrg(6,iatm)
      gci7 = gcrg(7,iatm); gci8 = gcrg(8,iatm)

      jfirst = iar1pb(4, iatm-1) + 1
      jlast  = iar1pb(2, iatm)
      dumx = ZERO; dumy = ZERO; dumz = ZERO
      do jp = jfirst, jlast
         jatm = iprshrt(jp)
         pair_correct = ONE
         jx = icrd(1,jatm); jy = icrd(2,jatm); jz = icrd(3,jatm)

         gcj1 = gcrg(1,jatm); gcj2 = gcrg(2,jatm)
         gcj3 = gcrg(3,jatm); gcj4 = gcrg(4,jatm)
         gcj5 = gcrg(5,jatm); gcj6 = gcrg(6,jatm)
         gcj7 = gcrg(7,jatm); gcj8 = gcrg(8,jatm)

         dijx  =       ix-jx; dijy  =       iy-jy; dijz  =       iz-jz
         dijx0 = abs(dijx-2); dijy0 = abs(dijy-2); dijz0 = abs(dijz-2)
         dijx1 = abs(dijx-1); dijy1 = abs(dijy-1); dijz1 = abs(dijz-1)
         dijx2 = abs(dijx+1); dijy2 = abs(dijy+1); dijz2 = abs(dijz+1)
         dijx3 = abs(dijx+2); dijy3 = abs(dijy+2); dijz3 = abs(dijz+2)
         dijx  = abs(dijx  ); dijy  = abs(dijy  ); dijz  = abs(dijz  )

         gcij( 1) = gci1*gcj1 + gci2*gcj2 + gci3*gcj3 + gci4*gcj4 + gci5*gcj5 + gci6*gcj6 + gci7*gcj7 + gci8*gcj8
         gcij( 2) = gci1*gcj2 + gci3*gcj4 + gci5*gcj6 + gci7*gcj8
         gcij( 3) = gci1*gcj3 + gci2*gcj4 + gci5*gcj7 + gci6*gcj8
         gcij( 4) = gci1*gcj4 + gci5*gcj8
         gcij( 5) = gci1*gcj5 + gci2*gcj6 + gci3*gcj7 + gci4*gcj8
         gcij( 6) = gci1*gcj6 + gci3*gcj8
         gcij( 7) = gci1*gcj7 + gci2*gcj8
         gcij( 8) = gci1*gcj8
         gcij( 9) = gci2*gcj1 + gci4*gcj3 + gci6*gcj5 + gci8*gcj7
         gcij(10) = gci2*gcj3 + gci6*gcj7
         gcij(11) = gci2*gcj5 + gci4*gcj7
         gcij(12) = gci2*gcj7
         gcij(13) = gci3*gcj1 + gci4*gcj2 + gci7*gcj5 + gci8*gcj6
         gcij(14) = gci3*gcj2 + gci7*gcj6
         gcij(15) = gci3*gcj5 + gci4*gcj6
         gcij(16) = gci3*gcj6
         gcij(17) = gci4*gcj1 + gci8*gcj5
         gcij(18) = gci4*gcj5
         gcij(19) = gci5*gcj1 + gci6*gcj2 + gci7*gcj3 + gci8*gcj4
         gcij(20) = gci5*gcj2 + gci7*gcj4
         gcij(21) = gci5*gcj3 + gci6*gcj4
         gcij(22) = gci5*gcj4
         gcij(23) = gci6*gcj1 + gci8*gcj3
         gcij(24) = gci6*gcj3
         gcij(25) = gci7*gcj1 + gci8*gcj2
         gcij(26) = gci7*gcj2
         gcij(27) = gci8*gcj1

         grdcoultmp = dble( &
              gci1*( gcj1*l_green(dijx ,dijy ,dijz ) + gcj2*l_green(dijx1,dijy ,dijz ) + &
                     gcj3*l_green(dijx ,dijy1,dijz ) + gcj4*l_green(dijx1,dijy1,dijz ) + &
                     gcj5*l_green(dijx ,dijy ,dijz1) + gcj6*l_green(dijx1,dijy ,dijz1) + &
                     gcj7*l_green(dijx ,dijy1,dijz1) + gcj8*l_green(dijx1,dijy1,dijz1) ) + &
              gci2*( gcj1*l_green(dijx2,dijy ,dijz ) + gcj2*l_green(dijx ,dijy ,dijz ) + &
                     gcj3*l_green(dijx2,dijy1,dijz ) + gcj4*l_green(dijx ,dijy1,dijz ) + &
                     gcj5*l_green(dijx2,dijy ,dijz1) + gcj6*l_green(dijx ,dijy ,dijz1) + &
                     gcj7*l_green(dijx2,dijy1,dijz1) + gcj8*l_green(dijx ,dijy1,dijz1) ) + &
              gci3*( gcj1*l_green(dijx ,dijy2,dijz ) + gcj2*l_green(dijx1,dijy2,dijz ) + &
                     gcj3*l_green(dijx ,dijy ,dijz ) + gcj4*l_green(dijx1,dijy ,dijz ) + &
                     gcj5*l_green(dijx ,dijy2,dijz1) + gcj6*l_green(dijx1,dijy2,dijz1) + &
                     gcj7*l_green(dijx ,dijy ,dijz1) + gcj8*l_green(dijx1,dijy ,dijz1) ) + &
              gci4*( gcj1*l_green(dijx2,dijy2,dijz ) + gcj2*l_green(dijx ,dijy2,dijz ) + &
                     gcj3*l_green(dijx2,dijy ,dijz ) + gcj4*l_green(dijx ,dijy ,dijz ) + &
                     gcj5*l_green(dijx2,dijy2,dijz1) + gcj6*l_green(dijx ,dijy2,dijz1) + &
                     gcj7*l_green(dijx2,dijy ,dijz1) + gcj8*l_green(dijx ,dijy ,dijz1) ) + &
              gci5*( gcj1*l_green(dijx ,dijy ,dijz2) + gcj2*l_green(dijx1,dijy ,dijz2) + &
                     gcj3*l_green(dijx ,dijy1,dijz2) + gcj4*l_green(dijx1,dijy1,dijz2) + &
                     gcj5*l_green(dijx ,dijy ,dijz ) + gcj6*l_green(dijx1,dijy ,dijz ) + &
                     gcj7*l_green(dijx ,dijy1,dijz ) + gcj8*l_green(dijx1,dijy1,dijz ) ) + &
              gci6*( gcj1*l_green(dijx2,dijy ,dijz2) + gcj2*l_green(dijx ,dijy ,dijz2) + &
                     gcj3*l_green(dijx2,dijy1,dijz2) + gcj4*l_green(dijx ,dijy1,dijz2) + &
                     gcj5*l_green(dijx2,dijy ,dijz ) + gcj6*l_green(dijx ,dijy ,dijz ) + &
                     gcj7*l_green(dijx2,dijy1,dijz ) + gcj8*l_green(dijx ,dijy1,dijz ) ) + &
              gci7*( gcj1*l_green(dijx ,dijy2,dijz2) + gcj2*l_green(dijx1,dijy2,dijz2) + &
                     gcj3*l_green(dijx ,dijy ,dijz2) + gcj4*l_green(dijx1,dijy ,dijz2) + &
                     gcj5*l_green(dijx ,dijy2,dijz ) + gcj6*l_green(dijx1,dijy2,dijz ) + &
                     gcj7*l_green(dijx ,dijy ,dijz ) + gcj8*l_green(dijx1,dijy ,dijz ) ) + &
              gci8*( gcj1*l_green(dijx2,dijy2,dijz2) + gcj2*l_green(dijx ,dijy2,dijz2) + &
                     gcj3*l_green(dijx2,dijy ,dijz2) + gcj4*l_green(dijx ,dijy ,dijz2) + &
                     gcj5*l_green(dijx2,dijy2,dijz ) + gcj6*l_green(dijx ,dijy2,dijz ) + &
                     gcj7*l_green(dijx2,dijy ,dijz ) + gcj8*l_green(dijx ,dijy ,dijz ) ) )

         grdcoul = grdcoul + grdcoultmp*pair_correct

         !-- PB decomp
         if(idecomp == 1 .or. idecomp == 2) call decpair(1,iatm,jatm,decfac*grdcoultmp)

         ffx = dble( gcij( 1)*(l_green(dijx1,dijy ,dijz ) - l_green(dijx2,dijy ,dijz )) + &
                     gcij( 2)*(l_green(dijx0,dijy ,dijz ) - l_green(dijx ,dijy ,dijz )) + &
                     gcij( 3)*(l_green(dijx1,dijy1,dijz ) - l_green(dijx2,dijy1,dijz )) + &
                     gcij( 4)*(l_green(dijx0,dijy1,dijz ) - l_green(dijx ,dijy1,dijz )) + &
                     gcij( 5)*(l_green(dijx1,dijy ,dijz1) - l_green(dijx2,dijy ,dijz1)) + &
                     gcij( 6)*(l_green(dijx0,dijy ,dijz1) - l_green(dijx ,dijy ,dijz1)) + &
                     gcij( 7)*(l_green(dijx1,dijy1,dijz1) - l_green(dijx2,dijy1,dijz1)) + &
                     gcij( 8)*(l_green(dijx0,dijy1,dijz1) - l_green(dijx ,dijy1,dijz1)) + &
                     gcij( 9)*(l_green(dijx ,dijy ,dijz ) - l_green(dijx3,dijy ,dijz )) + &
                     gcij(10)*(l_green(dijx ,dijy1,dijz ) - l_green(dijx3,dijy1,dijz )) + &
                     gcij(11)*(l_green(dijx ,dijy ,dijz1) - l_green(dijx3,dijy ,dijz1)) + &
                     gcij(12)*(l_green(dijx ,dijy1,dijz1) - l_green(dijx3,dijy1,dijz1)) + &
                     gcij(13)*(l_green(dijx1,dijy2,dijz ) - l_green(dijx2,dijy2,dijz )) + &
                     gcij(14)*(l_green(dijx0,dijy2,dijz ) - l_green(dijx ,dijy2,dijz )) + &
                     gcij(15)*(l_green(dijx1,dijy2,dijz1) - l_green(dijx2,dijy2,dijz1)) + &
                     gcij(16)*(l_green(dijx0,dijy2,dijz1) - l_green(dijx ,dijy2,dijz1)) + &
                     gcij(17)*(l_green(dijx ,dijy2,dijz ) - l_green(dijx3,dijy2,dijz )) + &
                     gcij(18)*(l_green(dijx ,dijy2,dijz1) - l_green(dijx3,dijy2,dijz1)) + &
                     gcij(19)*(l_green(dijx1,dijy ,dijz2) - l_green(dijx2,dijy ,dijz2)) + &
                     gcij(20)*(l_green(dijx0,dijy ,dijz2) - l_green(dijx ,dijy ,dijz2)) + &
                     gcij(21)*(l_green(dijx1,dijy1,dijz2) - l_green(dijx2,dijy1,dijz2)) + &
                     gcij(22)*(l_green(dijx0,dijy1,dijz2) - l_green(dijx ,dijy1,dijz2)) + &
                     gcij(23)*(l_green(dijx ,dijy ,dijz2) - l_green(dijx3,dijy ,dijz2)) + &
                     gcij(24)*(l_green(dijx ,dijy1,dijz2) - l_green(dijx3,dijy1,dijz2)) + &
                     gcij(25)*(l_green(dijx1,dijy2,dijz2) - l_green(dijx2,dijy2,dijz2)) + &
                     gcij(26)*(l_green(dijx0,dijy2,dijz2) - l_green(dijx ,dijy2,dijz2)) + &
                     gcij(27)*(l_green(dijx ,dijy2,dijz2) - l_green(dijx3,dijy2,dijz2)) )

         ffy = dble( gcij( 1)*(l_green(dijx ,dijy1,dijz ) - l_green(dijx ,dijy2,dijz )) + &
                     gcij( 2)*(l_green(dijx1,dijy1,dijz ) - l_green(dijx1,dijy2,dijz )) + &
                     gcij( 3)*(l_green(dijx ,dijy0,dijz ) - l_green(dijx ,dijy ,dijz )) + &
                     gcij( 4)*(l_green(dijx1,dijy0,dijz ) - l_green(dijx1,dijy ,dijz )) + &
                     gcij( 5)*(l_green(dijx ,dijy1,dijz1) - l_green(dijx ,dijy2,dijz1)) + &
                     gcij( 6)*(l_green(dijx1,dijy1,dijz1) - l_green(dijx1,dijy2,dijz1)) + &
                     gcij( 7)*(l_green(dijx ,dijy0,dijz1) - l_green(dijx ,dijy ,dijz1)) + &
                     gcij( 8)*(l_green(dijx1,dijy0,dijz1) - l_green(dijx1,dijy ,dijz1)) + &
                     gcij( 9)*(l_green(dijx2,dijy1,dijz ) - l_green(dijx2,dijy2,dijz )) + &
                     gcij(10)*(l_green(dijx2,dijy0,dijz ) - l_green(dijx2,dijy ,dijz )) + &
                     gcij(11)*(l_green(dijx2,dijy1,dijz1) - l_green(dijx2,dijy2,dijz1)) + &
                     gcij(12)*(l_green(dijx2,dijy0,dijz1) - l_green(dijx2,dijy ,dijz1)) + &
                     gcij(13)*(l_green(dijx ,dijy ,dijz ) - l_green(dijx ,dijy3,dijz )) + &
                     gcij(14)*(l_green(dijx1,dijy ,dijz ) - l_green(dijx1,dijy3,dijz )) + &
                     gcij(15)*(l_green(dijx ,dijy ,dijz1) - l_green(dijx ,dijy3,dijz1)) + &
                     gcij(16)*(l_green(dijx1,dijy ,dijz1) - l_green(dijx1,dijy3,dijz1)) + &
                     gcij(17)*(l_green(dijx2,dijy ,dijz ) - l_green(dijx2,dijy3,dijz )) + &
                     gcij(18)*(l_green(dijx2,dijy ,dijz1) - l_green(dijx2,dijy3,dijz1)) + &
                     gcij(19)*(l_green(dijx ,dijy1,dijz2) - l_green(dijx ,dijy2,dijz2)) + &
                     gcij(20)*(l_green(dijx1,dijy1,dijz2) - l_green(dijx1,dijy2,dijz2)) + &
                     gcij(21)*(l_green(dijx ,dijy0,dijz2) - l_green(dijx ,dijy ,dijz2)) + &
                     gcij(22)*(l_green(dijx1,dijy0,dijz2) - l_green(dijx1,dijy ,dijz2)) + &
                     gcij(23)*(l_green(dijx2,dijy1,dijz2) - l_green(dijx2,dijy2,dijz2)) + &
                     gcij(24)*(l_green(dijx2,dijy0,dijz2) - l_green(dijx2,dijy ,dijz2)) + &
                     gcij(25)*(l_green(dijx ,dijy ,dijz2) - l_green(dijx ,dijy3,dijz2)) + &
                     gcij(26)*(l_green(dijx1,dijy ,dijz2) - l_green(dijx1,dijy3,dijz2)) + &
                     gcij(27)*(l_green(dijx2,dijy ,dijz2) - l_green(dijx2,dijy3,dijz2)) )

         ffz = dble( gcij( 1)*(l_green(dijx ,dijy ,dijz1) - l_green(dijx ,dijy ,dijz2)) + &
                     gcij( 2)*(l_green(dijx1,dijy ,dijz1) - l_green(dijx1,dijy ,dijz2)) + &
                     gcij( 3)*(l_green(dijx ,dijy1,dijz1) - l_green(dijx ,dijy1,dijz2)) + &
                     gcij( 4)*(l_green(dijx1,dijy1,dijz1) - l_green(dijx1,dijy1,dijz2)) + &
                     gcij( 5)*(l_green(dijx ,dijy ,dijz0) - l_green(dijx ,dijy ,dijz )) + &
                     gcij( 6)*(l_green(dijx1,dijy ,dijz0) - l_green(dijx1,dijy ,dijz )) + &
                     gcij( 7)*(l_green(dijx ,dijy1,dijz0) - l_green(dijx ,dijy1,dijz )) + &
                     gcij( 8)*(l_green(dijx1,dijy1,dijz0) - l_green(dijx1,dijy1,dijz )) + &
                     gcij( 9)*(l_green(dijx2,dijy ,dijz1) - l_green(dijx2,dijy ,dijz2)) + &
                     gcij(10)*(l_green(dijx2,dijy1,dijz1) - l_green(dijx2,dijy1,dijz2)) + &
                     gcij(11)*(l_green(dijx2,dijy ,dijz0) - l_green(dijx2,dijy ,dijz )) + &
                     gcij(12)*(l_green(dijx2,dijy1,dijz0) - l_green(dijx2,dijy1,dijz )) + &
                     gcij(13)*(l_green(dijx ,dijy2,dijz1) - l_green(dijx ,dijy2,dijz2)) + &
                     gcij(14)*(l_green(dijx1,dijy2,dijz1) - l_green(dijx1,dijy2,dijz2)) + &
                     gcij(15)*(l_green(dijx ,dijy2,dijz0) - l_green(dijx ,dijy2,dijz )) + &
                     gcij(16)*(l_green(dijx1,dijy2,dijz0) - l_green(dijx1,dijy2,dijz )) + &
                     gcij(17)*(l_green(dijx2,dijy2,dijz1) - l_green(dijx2,dijy2,dijz2)) + &
                     gcij(18)*(l_green(dijx2,dijy2,dijz0) - l_green(dijx2,dijy2,dijz )) + &
                     gcij(19)*(l_green(dijx ,dijy ,dijz ) - l_green(dijx ,dijy ,dijz3)) + &
                     gcij(20)*(l_green(dijx1,dijy ,dijz ) - l_green(dijx1,dijy ,dijz3)) + &
                     gcij(21)*(l_green(dijx ,dijy1,dijz ) - l_green(dijx ,dijy1,dijz3)) + &
                     gcij(22)*(l_green(dijx1,dijy1,dijz ) - l_green(dijx1,dijy1,dijz3)) + &
                     gcij(23)*(l_green(dijx2,dijy ,dijz ) - l_green(dijx2,dijy ,dijz3)) + &
                     gcij(24)*(l_green(dijx2,dijy1,dijz ) - l_green(dijx2,dijy1,dijz3)) + &
                     gcij(25)*(l_green(dijx ,dijy2,dijz ) - l_green(dijx ,dijy2,dijz3)) + &
                     gcij(26)*(l_green(dijx1,dijy2,dijz ) - l_green(dijx1,dijy2,dijz3)) + &
                     gcij(27)*(l_green(dijx2,dijy2,dijz ) - l_green(dijx2,dijy2,dijz3)) )

         dumx = dumx + ffx; dumy = dumy + ffy; dumz = dumz + ffz
         frc(1,jatm) = frc(1,jatm) - ffx
         frc(2,jatm) = frc(2,jatm) - ffy
         frc(3,jatm) = frc(3,jatm) - ffz
      end do  !  jp = jfirst, jlast

      frc(1,iatm) = frc(1,iatm) + dumx
      frc(2,iatm) = frc(2,iatm) + dumy
      frc(3,iatm) = frc(3,iatm) + dumz
   end do  !  iatm = 1, ilast

   grdcoul = factor*grdcoul
   pbfrc  = pbfrc - factor1*frc

end subroutine pb_fdcoulomb
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ wrapper of the green's function common block values
_REAL_ function l_green (i,j,k)

   implicit none
   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green
   integer i,j,k

   if ( i <= 40  .and. j <= 40 .and. k <= 40 ) then
      l_green = green(i,j,k)
   else
      l_green = ONE/sqrt(dble(i*i+j*j+k*k))
   end if

end function l_green
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ pure pairwise charge interaction energy for linear PB equations. no cutoff
subroutine pb_lpbene( verbose,eneout,natom,atmlast,ifcap,npdec,idecomp,m,irespw,ipres, &
                     eel,insas,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd
   use decomp, only: decpair
   implicit none

   ! Passed variables

   logical verbose, eneout
   integer natom,atmlast,ifcap,npdec,idecomp,m,irespw(*),ipres(*)
   integer insas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)
   _REAL_ eel

   ! Local variables

   integer  i, j, k, iatm, jatm, ip
   _REAL_, parameter :: smallcrg = 0.5d0
   _REAL_ srfcrg, factor, scalfact
   _REAL_ g(3), x(3), dx(3), crd(3), dx2, dist, rdist, acg
   _REAL_ eelrf
   _REAL_ d2inv, dinv, de
   _REAL_, allocatable :: pol_charge(:)

   ! initialization

   factor = THREE*h/(TWOPI)
   scalfact = ONE

   ! compute polarization charges on the boundary grid points and
   ! report total in srfcrg
   ! pol_charge is in au

   allocate (pol_charge(1:nbnd))
   call get_charge_pol(nbnd,xm,ym,zm,phi,cphi,chgrd,pol_charge,srfcrg)

   ! for InsightII display
   !open (unit=55, file='ms.dot')
   !write (55, '("DOTS")')

   eel = ZERO
   do ip = 1, nbnd
      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k

      ! project the surface grid point on to the molecular surface, crd() is the new coord

      if      ( iatm > 0 ) then
         ! the contact boundary grid points are projected to the atom spheres
         x(1:3) = acrd(1:3,iatm)
         dist = radi(iatm)
      else if ( iatm < 0 ) then
         ! the reentry boundary grid points are projected to the solvent spheres
         x(1:3) = arccrd(1:3,-iatm)
         dist = dprob
      else
         ! otherwise do not project. The setting should make crd = g
         x(1:3) = g(1:3)
         dist = ONE
      end if

      dx = g - x; dx2 = dx(1)**2 + dx(2)**2 + dx(3)**2
      if ( dx2 == ZERO ) then
         rdist = ONE
      else
         rdist = dist*ONE/sqrt(dx2)
      end if
      crd = x + dx*rdist

      ! for InsightII display
      !write (55,'(4(f8.3,2x))') crd(1:3), 300.

      ! retrieve compouted pol charge

      acg = pol_charge(ip)

      ! compute reaction field energy and forces due to this boundary grid point

      eelrf = ZERO
      if (idecomp < 3) then
         do jatm = 1, atmlast
            if ((ifcap == 2 .or. ifcap == 5) .and. outflag(jatm) == 1) cycle
            if ( ligand ) then
               if (liveflag(jatm) == 0) cycle
            endif

            dx(1) = crd(1) - acrd(1,jatm)
            dx(2) = crd(2) - acrd(2,jatm)
            dx(3) = crd(3) - acrd(3,jatm)
            dinv = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2); d2inv = dinv**2

            de = acg*acrg(jatm)*dinv
            eelrf = eelrf + de

            !-- PB decomp
            if(idecomp == 1 .or. idecomp == 2) call decpair(1,jatm,jatm,de*AMBER_ELECTROSTATIC2*HALF)

         end do
      !-- PB pairwise decomp
      else if (idecomp >= 3) then
         do n = 1, npdec
            do jatm = ipres(irespw(n)), ipres(irespw(n)+1)-1
               if ( ligand ) then
                  if (liveflag(jatm) == 0) cycle
               endif
               dx(1) = crd(1) - acrd(1,jatm)
               dx(2) = crd(2) - acrd(2,jatm)
               dx(3) = crd(3) - acrd(3,jatm)
               dinv = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
               de = acg*acrg(jatm)*dinv
               eelrf = eelrf + de
               call decpair(-1,ipres(irespw(m)),ipres(irespw(n)),de*AMBER_ELECTROSTATIC2*HALF)
            end do
         end do
      end if

      ! collecting energy

      eel = eel + eelrf

   end do ! end of ip = 1, nbnd

   ! for InsightII display
   !close(55)
   !stop

   if ( scalerf .and. abs(totcrg) > smallcrg ) then
      scalfact = abs( totcrg/srfcrg*(ONE/epsin - ONE/epsout)*eps0 )
      srfcrg = scalfact*srfcrg
      eel = scalfact*eel
   end if

   ! converting to kcal/mol

   eel = HALF*AMBER_ELECTROSTATIC2*eel

   if ( eneout ) then
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eel
   end if

   deallocate (pol_charge)

end subroutine pb_lpbene
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ pure pairwise charge interaction energy for nonlinear PB equations. no cutoff
subroutine pb_npbene(nbnd,xm,ym,zm,chgrd,insas,phi,ion_con,osmene)

   implicit none

   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green

   ! passed variables

   integer nbnd,xm,ym,zm
   _REAL_ chgrd(xm,ym,zm)
   integer insas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ phi(xm,ym,zm)
   _REAL_ ion_con(xm,ym,zm)
   _REAL_ epb0, osmene

   ! local variables

   integer i, j, k, ip,i0,j0,k0,idx,idy,idz,nbound
   _REAL_ total_phi(7)
   _REAL_ rdis,rdisinv,qtmp,qtmp1
   _REAL_, allocatable :: bound_crg(:)
   _REAL_, allocatable :: pol_charge(:)

   epb0 = ZERO

   ! compute conductor grid boundary charges

   nbound = 2*xm*ym+2*xm*zm+2*zm*ym
   allocate(bound_crg(nbound))
   do j = 1, ym; do i = 1, xm
      bound_crg(i+(j-1)*xm)=-h*phi(i,j,1)*eps0
      bound_crg(i+(j-1)*xm+xm*ym)=-h*phi(i,j,zm)*eps0
   end do; end do
   do k = 1, zm; do i = 1, xm;
      bound_crg(i+(k-1)*xm+2*xm*ym)=-h*phi(i,1,k)*eps0
      bound_crg(i+(k-1)*xm+2*xm*ym+xm*zm)=-h*phi(i,ym,k)*eps0
   end do; end do
   do k = 1, zm; do j = 1, ym;
      bound_crg(j+(k-1)*ym+2*xm*ym+2*xm*zm)=-h*phi(1,j,k)*eps0
      bound_crg(j+(k-1)*ym+2*xm*ym+2*xm*zm+ym*zm)=-h*phi(xm,j,k)*eps0
   end do; end do

  ! compute polarization charges

   allocate (pol_charge(1:nbnd))
   do ip = 1, nbnd

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip)

      total_phi(1) = phi(i,  j,k)
      total_phi(2) = phi(i-1,j,k)
      total_phi(3) = phi(i+1,j,k)
      total_phi(4) = phi(i,j-1,k)
      total_phi(5) = phi(i,j+1,k)
      total_phi(6) = phi(i,j,k-1)
      total_phi(7) = phi(i,j,k+1)

      pol_charge(ip) = SIX*h*(  total_phi(1)-&
           SIXTH*( total_phi(2)+total_phi(3)+&
                   total_phi(4)+total_phi(5)+&
                   total_phi(6)+total_phi(7) ) )
      if ( insas(i,j,k) >= 0 ) then
         pol_charge(ip) = pol_charge(ip)*eps0-chgrd(i,j,k)/(epsin /eps0)-ion_con(i,j,k)/(epsout/eps0)
      else
         pol_charge(ip) = pol_charge(ip)*eps0-chgrd(i,j,k)/(epsout/eps0)-ion_con(i,j,k)/(epsout/eps0)
      endif

   end do

   ! charge-charge interactions

   do k = 1, zm; do j = 1, ym; do i = 1, xm

      ! This is the fixed (atomic charges on grid nodes

      if (chgrd(i,j,k)==ZERO .and. ion_con(i,j,k)==ZERO) cycle
      qtmp1 = chgrd(i,j,k)-ion_con(i,j,k)
      if ( qtmp1 == ZERO ) cycle

      ! part a: interactions with polarization charges

      do ip = 1, nbnd
         i0 = iepsav(1,ip); j0 = iepsav(2,ip); k0 = iepsav(3,ip)
         idx = abs(i0-i); idy = abs(j0-j); idz = abs(k0-k)

         if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
            rdisinv = green(idx,idy,idz)
            epb0 = epb0 + pol_charge(ip)*qtmp1*rdisinv/h
         else
            rdis = sqrt(dble(idx**2 + idy**2 + idz**2))
            epb0 = epb0 + pol_charge(ip)*qtmp1/rdis/h
         end if
      end do

      ! part b: interactions with other atomic charges and ion charges

      do k0 = 1, zm; do j0 = 1, ym; do i0 = 1, xm
         if ( insas(i0,j0,k0) >= 0 ) then
            qtmp=chgrd(i0,j0,k0)*(eps0/ epsin)+ion_con(i0,j0,k0)*(eps0/epsout)
         else
            qtmp=chgrd(i0,j0,k0)*(eps0/epsout)+ion_con(i0,j0,k0)*(eps0/epsout)
         endif
         if ( qtmp == ZERO ) cycle
         idx = abs(i0-i); idy = abs(j0-j); idz = abs(k0-k)
         if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
            rdisinv = green(idx,idy,idz)
            epb0 = epb0 + qtmp*qtmp1*rdisinv/h
         else
            rdis = sqrt(dble(idx**2 + idy**2 + idz**2))
            epb0 = epb0 + qtmp*qtmp1/rdis/h
         end if
      end do; end do; end do

      ! part c: interactions with grid boundary charges, when bcopt = 1, i.e. conductor
      ! Note that this is the only case that we can easily compare original full
      ! pb energy formulation with the charge-view formulation.

      do j0 = 1, ym; do i0 = 1, xm

         k0 = 0; qtmp=bound_crg(i0+(j0-1)*xm)
         idx = abs(i0-i); idy=abs(j0-j); idz=abs(k0-k)
         if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
            rdisinv = green(idx,idy,idz)
            epb0 = epb0 + qtmp1*qtmp*rdisinv/h
         else
            rdis = sqrt(dble(idx**2 + idy**2 + idz**2))
            epb0 = epb0 + qtmp1*qtmp/rdis/h
         end if

         k0 = zm+1; qtmp=bound_crg(i0+(j0-1)*xm+xm*ym)
         idx=abs(i0-i); idy=abs(j0-j); idz=abs(k0-k)
         if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
            rdisinv = green(idx,idy,idz)
            epb0 = epb0 + qtmp1*qtmp*rdisinv/h
         else
            rdis = sqrt(dble(idx**2 + idy**2 + idz**2))
            epb0 = epb0 + qtmp1*qtmp/rdis/h
         end if
      end do; end do

      do k0 = 1, zm; do i0 = 1, xm

         j0 = 0; qtmp=bound_crg(i0+(k0-1)*xm+2*xm*ym)
         idx = abs(i0-i); idy=abs(j0-j); idz=abs(k0-k)
         if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
            rdisinv = green(idx,idy,idz)
            epb0 = epb0 + qtmp1*qtmp*rdisinv/h
         else
            rdis = sqrt(dble(idx**2 + idy**2 + idz**2))
            epb0 = epb0 + qtmp1*qtmp/rdis/h
         end if

         j0 = ym+1; qtmp=bound_crg(i0+(k0-1)*xm+2*xm*ym+xm*zm)
         idx = abs(i0-i); idy = abs(j0-j); idz = abs(k0-k)
         if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
            rdisinv = green(idx,idy,idz)
            epb0 = epb0 + qtmp1*qtmp*rdisinv/h
         else
            rdis = sqrt(dble(idx**2 + idy**2 + idz**2))
            epb0 = epb0 + qtmp1*qtmp/rdis/h
         end if
      end do; end do

      do k0 = 1, zm; do j0 = 1, ym

         i0 = 0; qtmp=bound_crg(j0+(k0-1)*ym+2*xm*ym+2*xm*zm)
         idx = abs(i0-i); idy = abs(j0-j); idz = abs(k0-k)
         if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
            rdisinv = green(idx,idy,idz)
            epb0 = epb0 + qtmp1*qtmp*rdisinv/h
         else
            rdis = sqrt(dble(idx**2 + idy**2 + idz**2))
            epb0 = epb0 + qtmp1*qtmp/rdis/h
         end if

         i0 = xm+1; qtmp=bound_crg(j0+(k0-1)*ym+2*xm*ym+2*xm*zm+ym*zm)
         idx = abs(i0-i); idy = abs(j0-j); idz = abs(k0-k)
         if (idx <= 40 .and. idy <= 40 .and. idz <= 40) then
            rdisinv = green(idx,idy,idz)
            epb0 = epb0 + qtmp1*qtmp*rdisinv/h
         else
            rdis = sqrt(dble(idx**2 + idy**2 + idz**2))
            epb0 = epb0 + qtmp1*qtmp/rdis/h
         end if
      end do; end do

   end do; end do; end do

   deallocate(bound_crg)
   deallocate(pol_charge)

   write(6,'(a,e20.9)') ' full nonlinear PB charge-view energy in kcal/mol',&
      (0.5d0*epb0/(FOURPI*eps0)+osmene)*frcfac

end subroutine pb_npbene
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_dbfrc_fld(verbose,eneout,natom,f,epsx,epsy,epsz,phi,cphi)

   use solvent_accessibility, only : dprob, radi, arccrd
   implicit none

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ f(3,natom)
   _REAL_ epsx(0:xm,1:ym,1:zm)
   _REAL_ epsy(1:xm,0:ym,1:zm)
   _REAL_ epsz(1:xm,1:ym,0:zm)
   _REAL_ phi(xm,ym,zm)
   _REAL_ cphi(xm,ym,zm)

   ! local variables

   integer i, j, k, iatm
   integer xp, yp, zp
   _REAL_ df, factor, factor1, sfactor
   _REAL_ lEx, lEy, lEz
   _REAL_ dx(3), dum(3)
   _REAL_ dfx, dfy, dfz
   _REAL_ adx, ady, adz
   _REAL_ x(3), crd(3)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_ rsphere, cut, r, r2, h2
   _REAL_ ref(3)

   ! initialization

   ref = ZERO

   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO

   h2 = h*h
   sfactor = -HALF*(epsout-epsin)/(epsin*epsout)

   ! contributions from epsx boundary edges

   do xp = 1, nbndx
      i = iepsavx(1,xp); j = iepsavx(2,xp); k = iepsavx(3,xp); iatm = iepsavx(4,xp)
      lEx = phi(i+1,j,k)+cphi(i+1,j,k) - phi(i,j,k)-cphi(i,j,k)
      crd(1) = gox + h*(i+fedgex(xp)); crd(2) = goy + h*j; crd(3) = goz + h*k
      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         dx = crd - x
!        write(58,'(a,3(f20.15),2x)') "H", crd(1:3)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
      else
         x(1:3) = arccrd(1:3,-iatm)
         dx = x - crd
!        write(59,'(a,3(f20.15),2x)') "H", crd(1:3)
         rsphere = dprob
      end if
      adx = abs(dx(1))

      df = sfactor*lEx*lEx*epsx(i,j,k)**2
      factor1 = df/adx
      dum = factor1*dx

      ref(1)=ref(1)+dum(1)
      ref(2)=ref(2)+dum(2)
      ref(3)=ref(3)+dum(3)

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)

   enddo

   ! contributions from epsy boundary grid edges

   do yp = 1, nbndy
      i = iepsavy(1,yp); j = iepsavy(2,yp); k = iepsavy(3,yp); iatm = iepsavy(4,yp)
      lEy = phi(i,j+1,k)+cphi(i,j+1,k) - phi(i,j,k)-cphi(i,j,k)
      crd(1) = gox + h*i; crd(2) = goy + h*(j+fedgey(yp)); crd(3) = goz + h*k

      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         dx = crd - x
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
      else
         x(1:3) = arccrd(1:3,-iatm)
         dx = x - crd
         rsphere = dprob
      end if
      ady = abs(dx(2))

      df = sfactor*lEy*lEy*epsy(i,j,k)**2
      factor1= df/ady
      dum = factor1*dx

      ref(1)=ref(1)+dum(1)
      ref(2)=ref(2)+dum(2)
      ref(3)=ref(3)+dum(3)

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)

   enddo

   ! contributions from epsz boundary grid edges

   do zp = 1, nbndz
      i = iepsavz(1,zp); j = iepsavz(2,zp); k = iepsavz(3,zp); iatm = iepsavz(4,zp)
      lEz = phi(i,j,k+1)+cphi(i,j,k+1) - phi(i,j,k)-cphi(i,j,k)
      crd(1) = gox + h*i; crd(2) = goy + h*j; crd(3) = goz + h*(k+fedgez(zp))
      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         dx = crd - x
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
      else
         x(1:3) = arccrd(1:3,-iatm)
         dx = x - crd
         rsphere = dprob
      end if
      adz = abs(dx(3))

      df = sfactor*lEz*lEz*epsz(i,j,k)**2
      factor1= df/adz
      dum = factor1*dx

      ref(1)=ref(1)+dum(1)
      ref(2)=ref(2)+dum(2)
      ref(3)=ref(3)+dum(3)

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)

   end do

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz

   if ( verbose .and. eneout ) then
      ref = 0.0d0
      write(102,*) ' :::: Atomic DB forces ::::'
      do iatm = 1, natom
         write(102,'(3e20.6)') dbx(iatm)*frcfac,dby(iatm)*frcfac,dbz(iatm)*frcfac
      end do
      ref = ref*frcfac
      write(6,'(a,3f12.4)') 'Vector DB forces is',ref(1:3)
      write(6,'(a,f12.4)') 'Total DB forces is',sqrt(ref(1)**2+ref(2)**2+ref(3)**2)
   end if

   f(1,1:natom) = f(1,1:natom) + dbx(1:natom)
   f(2,1:natom) = f(2,1:natom) + dby(1:natom)
   f(3,1:natom) = f(3,1:natom) + dbz(1:natom)

end subroutine pb_dbfrc_fld
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_dbfrc_crg( verbose,eneout,natom,eelrf,f,epsx,epsy,epsz,zv,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd
   implicit none

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ eelrf
   _REAL_ f(3,natom)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)
   _REAL_ zv(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)

   ! local variables

   integer, parameter :: n_point = 27
   integer i, j, k, iatm, ip
   _REAL_, parameter :: small_dn = 1.d-3 * AMBER_ELECTROSTATIC
   _REAL_ g(3), x(3), dx(3), crd(3), rn(3), sgn, rdx, rsphere, crg
   _REAL_ srfcrg, coulomb(3)
   _REAL_ fx0, fy0, fz0, rh
   _REAL_ up, dudxi0, dudyi0, dudzi0
   _REAL_ dn, dt2, dum(3)
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_ Ex, Ey, Ez, w1, w2, w3, cnnt, factor, reps0
   integer sgnx, sgny, sgnz
   _REAL_, allocatable :: pol_charge(:)

   ! initialization

   if ( isurfchg > 0 ) open(131,file='srfchg.pos',form="formatted")

   rh = ONE/h; reps0 = ONE/eps0
   factor = FOURPI*eps0*AMBER_ELECTROSTATIC
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO
   x = ZERO; rsphere = ZERO; sgn = ONE; coulomb = ZERO

   ! compute polarization charges on the boundary grid points and
   ! report total in srfcrg
   ! pol_charge is in au

   allocate (pol_charge(1:nbnd))
   call get_charge_pol(nbnd,xm,ym,zm,phi,cphi,chgrd,pol_charge,srfcrg)

   ! main double loops over polarization charges and atom charges

   eelrf = ZERO
   do ip = 1, nbnd

      ! collect boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
      crg = pol_charge(ip)*AMBER_ELECTROSTATIC

      ! project the boundary grid point onto the molecular surface, crd() is the
      ! projected coord, x() is the atom/probe center coord, and fx/y/z0 is the grid
      ! version of crd()

      if ( iatm == 0 ) then
         write(6,'(a)') 'PB Bomb in pb_dbfrc_crg(): can not find projection atom/probe'
         call mexit(6, 1)
      end if

      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         sgn = ONE
      else if ( iatm < 0 ) then
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         sgn = -ONE
      end if
      dx = g - x
      rdx = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
      crd = x + dx*rdx*rsphere

      !  write surface charges to output

      if (isurfchg > 0 ) then
        write(131,'(3f15.5,f10.6)') crd(1:3),crg*INV_AMBER_ELECTROSTATIC
      endif

      ! normal vector, it should be sgn*dx*rdx

      rn = sgn*dx*rdx

      ! compute reaction field energy, forces, and coulomb field
      ! between the current surface charge and natom atomic charges
      ! inner loop over atoms...

      call get_coulomb(natom,crg,crd,acrg,ax,ay,az,qex,qey,qez,eelrf,coulomb)

      ! now it is time to compute dbf
      dn = ZERO
      if ( iatm > 0 ) then

         ! interpolate E

         ! grid unit of the projected boundary point

         fx0 = (crd(1)-gox)*rh
         fy0 = (crd(2)-goy)*rh
         fz0 = (crd(3)-goz)*rh

         ! compute E on the inner side of surface position crd(1:3)

        !call gradu(xm,ym,zm,-ONE,n_point,10,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,zv)
         call onesided(xm,ym,zm,10,n_point,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,zv,h)

         dudxi0 = -dudxi0*factor
         dudyi0 = -dudyi0*factor
         dudzi0 = -dudzi0*factor

         ! add the coulomb field to get the total E of inner side
         ! converted to electric displacement

         dudxi0 = (dudxi0 + coulomb(1))*epsin*reps0
         dudyi0 = (dudyi0 + coulomb(2))*epsin*reps0
         dudzi0 = (dudzi0 + coulomb(3))*epsin*reps0

         ! compute the DBF as HALF*Q*D^2/Dn
         ! when denominator dn is tiny, use the limiting law.

         ! dt2 = dudxi0*dudxi0 + dudyi0*dudyi0 + dudzi0*dudzi0
         dn = dudxi0*rn(1)+dudyi0*rn(2)+dudzi0*rn(3)

      elseif ( iatm < 0 ) then

         ! weighted sum of D
         ! bcopt = 9 is better?

         cnnt = ZERO

         !  contact
         !  consider inside, positive rn
         !  should be 1, but inside -2, outside 2
         !  reentry
         !  consider inside, positive rn
         !  should be -1, but there is sgn in rn, and inside -1, outside 1
         sgnx = -1; sgny = -1; sgnz = -1

         Ex = ZERO; Ey = ZERO; Ez = ZERO
         w1 = ZERO; w2 = ZERO; w3 = ZERO
         if ( rn(1) < ZERO ) sgnx = -sgnx
         if ( rn(2) < ZERO ) sgny = -sgny
         if ( rn(3) < ZERO ) sgnz = -sgnz
         if ( zv(i,j,k) < ZERO ) then
            sgnx = -sgnx
            sgny = -sgny
            sgnz = -sgnz
         end if
         if ( zv(i,j,k)*zv(i+sgnx,j,k) < ZERO .and. abs(rn(1)) > 1.d-2 ) then
            !  Ex = dble(sign(1,sgnx))*(phi(i,j,k)-phi(i+sgnx,j,k))
            Ex = dble(sign(1,sgnx))*(phi(i,j,k)+cphi(i,j,k)-phi(i+sgnx,j,k)-cphi(i+sgnx,j,k))
            w1 = ONE/rn(1)
            cnnt = cnnt + ONE
         end if
         if ( zv(i,j,k)*zv(i,j+sgny,k) < ZERO .and. abs(rn(2)) > 1.d-2 ) then
            !  Ey = dble(sign(1,sgny))*(phi(i,j,k)-phi(i,j+sgny,k))
            Ey = dble(sign(1,sgny))*(phi(i,j,k)+cphi(i,j,k)-phi(i,j+sgny,k)-cphi(i,j+sgny,k))
            w2 = ONE/rn(2)
            cnnt = cnnt + ONE
         end if
         if ( zv(i,j,k)*zv(i,j,k+sgnz) < ZERO .and. abs(rn(3)) > 1.d-2 ) then
            !  Ez = dble(sign(1,sgnz))*(phi(i,j,k)-phi(i,j,k+sgnz))
            Ez = dble(sign(1,sgnz))*(phi(i,j,k)+cphi(i,j,k)-phi(i,j,k+sgnz)-cphi(i,j,k+sgnz))
            w3 = ONE/rn(3)
            cnnt = cnnt + ONE
         end if
         if ( cnnt > ZERO ) then
            if ( sgnx == 1 ) sgnx = 0
            if ( sgny == 1 ) sgny = 0
            if ( sgnz == 1 ) sgnz = 0
            Ex = Ex*epsx(i+sgnx,j,k)*reps0
            Ey = Ey*epsy(i,j+sgny,k)*reps0
            Ez = Ez*epsz(i,j,k+sgnz)*reps0
            ! dt2 = (Ex*Ex+Ey*Ey+Ez*Ez)*rh**2*factor**2
            dn = (Ex*w1+Ey*w2+Ez*w3)/cnnt*rh*factor
         else
!           write(6,'(a)') "PB warning: cnnt = 0, at", i, j, k
            continue
         end if

      end if

      dum = HALF*crg*dn*rn

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)
   end do

   deallocate (pol_charge)

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz

   if ( verbose .and. eneout ) then
      open (unit = 103, file = 'force.dat')
      write(103,*) ' :::: Atomic qE forces ::::'
      do iatm = 1, natom
         write(103,'(3e20.6)') qex(iatm),qey(iatm),qez(iatm)
      end do
      write(103,*) ' :::: Atomic contact forces ::::'
      do iatm = 1, natom
         write(103,'(3e20.6)') cnx(iatm),cny(iatm),cnz(iatm)
      end do
      write(103,*) ' :::: Atomic reentry forces ::::'
      do iatm = 1, natom
         write(103,'(3e20.6)') rnx(iatm),rny(iatm),rnz(iatm)
      end do
      write(103,*) ' :::: Atomic DB forces ::::'
      do iatm = 1, natom
         write(103,'(3e20.6)') dbx(iatm),dby(iatm),dbz(iatm)
      end do
      close(103)
   end if

   f(1,1:natom) = f(1,1:natom) + dbx(1:natom) + qex(1:natom)
   f(2,1:natom) = f(2,1:natom) + dby(1:natom) + qey(1:natom)
   f(3,1:natom) = f(3,1:natom) + dbz(1:natom) + qez(1:natom)

   eelrf = HALF*eelrf

   if ( verbose .and. eneout ) then
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eelrf
      write(6, '(1x,a,3f12.4)') 'Total solvation force',sum(dbx(1:natom))+sum(qex(1:natom)),&
         sum(dby(1:natom))+sum(qey(1:natom)), sum(dbz(1:natom))+sum(qez(1:natom))
   end if

   if( isurfchg > 0 ) close(131)

end subroutine pb_dbfrc_crg
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric boundary forces, second field-based strategy
!+ loop over boundary grid edges, three loops
subroutine pb_dbfrc_fld2( verbose,eneout,natom,eelrf,f,epsx,epsy,epsz,zv,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd
   implicit none

   ! Passed variables

   logical verbose, eneout
   integer natom
   _REAL_ zv(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)
   _REAL_ eelrf, f(3,natom)

   ! Local variables

   integer i, j, k, iatm, ip
   _REAL_ x(3), crd(3), crd1(3)
   _REAL_ fbnd, dum(3)
   _REAL_ fx0, fy0, fz0
   _REAL_ rn(1:3), rsphere, dr, rn1(1:3)
   _REAL_ ds, srfarea, srfarea1
   _REAL_ crg, h2, hh, r1, r2, r3
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_ up, dudxi0, dudyi0, dudzi0
   _REAL_ dudxo0, dudyo0, dudzo0
   _REAL_ dudni, dudno, coulomb(3)
   _REAL_ E2

   rsphere = -1d0

   ! initialization for DBF

   h2 = h*h; hh = HALF*h
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO
   coulomb = ZERO

   eelrf = ZERO
   srfarea = ZERO
   srfarea1 = ZERO
   do ip = 1, nbndx
      i = iepsavx(1,ip); j = iepsavx(2,ip); k = iepsavx(3,ip); iatm = iepsavx(4,ip)
      fx0 = i + fedgex(ip); fy0 = j; fz0 = k
      crd1(1) = gox + h*i + hh; crd1(2) = goy + h*j; crd1(3) = goz + h*k
      crd(1) = crd1(1) - hh + fedgex(ip)*h; crd(2) = crd1(2); crd(3) = crd1(3)

      if ( iatm == 0 ) then
         x(1:3) = ZERO
         write(6,'(a)') 'PBSA BOMB: can not find projection atom/probe'
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         rn = crd - x
         rn1 = crd1 - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn = x - crd
         rn1 = x - crd1
      end if

      ! surface element

      dr = abs(rn1(1))
      r2 = ONE/(rn1(1)**2+rn1(2)**2+rn1(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds = dr*rsphere**2*r3*h2
      srfarea = srfarea + ds

#     include "pb_dbfrc_fld2.h"

   end do !nbndx

   do ip = 1, nbndy
      i = iepsavy(1,ip); j = iepsavy(2,ip); k = iepsavy(3,ip); iatm = iepsavy(4,ip)
      fx0 = i; fy0 = j + fedgey(ip); fz0 = k
      crd1(1) = gox + h*i; crd1(2) = goy + h*j + hh; crd1(3) = goz + h*k
      crd(1) = crd1(1); crd(2) = crd1(2) - hh + fedgey(ip)*h; crd(3) = crd1(3)

      if ( iatm == 0 ) then
         x(1:3) = ZERO
         write(6,'(a)') 'PBSA BOMB: can not find projection atom/probe'
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         rn = crd - x
         rn1 = crd1 - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn = x - crd
         rn1 = x - crd1
      end if

      ! surface element

      dr = abs(rn1(2))
      r2 = ONE/(rn1(1)**2+rn1(2)**2+rn1(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds = dr*rsphere**2*r3*h2
      srfarea = srfarea + ds

#     include "pb_dbfrc_fld2.h"

   end do !nbndy

   do ip = 1, nbndz
      i = iepsavz(1,ip); j = iepsavz(2,ip); k = iepsavz(3,ip); iatm = iepsavz(4,ip)
      fx0 = i; fy0 = j; fz0 = k + fedgez(ip)
      crd1(1) = gox + h*i ; crd1(2) = goy + h*j; crd1(3) = goz + h*k + hh
      crd(1) = crd1(1); crd(2) = crd1(2); crd(3) = crd1(3) - hh + fedgez(ip)*h

      if ( iatm == 0 ) then
         x(1:3) = ZERO
         write(6,'(a)') 'PBSA BOMB: can not find projection atom/probe'
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         rn = crd - x
         rn1 = crd1 - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn = x - crd
         rn1 = x - crd1
      end if

      ! surface element

      dr = abs(rn1(3))
      r2 = ONE/(rn1(1)**2+rn1(2)**2+rn1(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds = dr*rsphere**2*r3*h2
      srfarea = srfarea + ds

#     include "pb_dbfrc_fld2.h"

   end do !nbndz

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz

   if ( verbose .and. eneout ) then
      open (unit = 103, file = 'force.dat')
      write(103,*) ' :::: Atomic DB forces ::::'
      do iatm = 1, natom
         write(103,'(3e20.6)') dbx(iatm),dby(iatm),dbz(iatm)
      end do
   end if

   f(1,1:natom) = f(1,1:natom) + dbx(1:natom) + qex(1:natom)
   f(2,1:natom) = f(2,1:natom) + dby(1:natom) + qey(1:natom)
   f(3,1:natom) = f(3,1:natom) + dbz(1:natom) + qez(1:natom)

   if ( verbose .and. eneout ) then
      write(6, '(1x,a,f12.4)') 'Total molecular surface', srfarea
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eelrf
      write(6, '(1x,a,3f12.4)') 'Total solvation force', sum(dbx(1:natom))+sum(qex(1:natom)),&
         sum(dby(1:natom))+sum(qey(1:natom)),sum(dbz(1:natom))+sum(qez(1:natom))
   end if

end subroutine pb_dbfrc_fld2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric boundary forces, second charge-based strategy
subroutine pb_dbfrc_crg2( verbose,eneout,natom,eelrf,f,epsx,epsy,epsz,zv,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd
   implicit none

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ eelrf
   _REAL_ f(3,natom)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)
   _REAL_ zv(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)

   ! local variables

   integer, parameter :: n_point = 7
   integer i, j, k, iatm, jatm, matm, natm, iarc, ip
   _REAL_ srfcrg, crg
   _REAL_ g(3), x(3), dx(3), crd(3), rn(3), rdx, rsphere, sgn
   _REAL_ dum(3), coulomb(3)
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_ fx0, fy0, fz0, rh
   _REAL_ up,dudxi0, dudyi0, dudzi0, dudni, E2
   _REAL_, allocatable :: pol_charge(:)

   ! mjhsieh: warning eliminator

   sgn = ZERO

   ! initialization

   rh = ONE/h
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO
   x = ZERO; rsphere = ZERO; coulomb = ZERO

   ! compute polarization charges on the boundary grid points and
   ! report total in srfcrg
   ! pol_charge is in au

   allocate (pol_charge(1:nbnd))
   call get_charge_pol(nbnd,xm,ym,zm,phi,cphi,chgrd,pol_charge,srfcrg)

   ! main double loops over polarization charges and atom charges

   eelrf = ZERO
   do ip = 1, nbnd

      ! collect boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
      crg = pol_charge(ip)*AMBER_ELECTROSTATIC

      ! project the surface grid point on to the molecular surface, crd() is the
      ! new coord, and x() is the atom/probe coord, dx() is the projection
      ! direction vector, fx/y/z0 is the grid version of crd()

      if ( iatm == 0 ) then
         write(6,'(a)') 'PB Bomb in pb_dbfrc_crg2(): can not find projection atom/probe'
         call mexit(6, 1)
      end if

      ! collect data

      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         sgn = ONE
      else if ( iatm < 0 ) then
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         sgn = -ONE
      end if

      ! project onto the sphere

      dx = g - x
      rdx = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
      crd = x + dx*rdx*rsphere

      ! normal direction vector

      rn = sgn*dx*rdx

      ! grid unit of the projected boundary point

      fx0 = (crd(1)-gox)*rh
      fy0 = (crd(2)-goy)*rh
      fz0 = (crd(3)-goz)*rh

      ! compute reaction field energy, forces, and coulomb field
      ! between the current surface charge and natom atomic charges
      ! inner loop over atoms...

      call get_coulomb(natom,crg,crd,acrg,ax,ay,az,qex,qey,qez,eelrf,coulomb)

      ! now it is time to compute dbf
      ! compute E on the inner side of surface position crd(1:3)

     !call gradu(xm,ym,zm,-ONE,n_point,4,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,zv)
      call onesided(xm,ym,zm,4,n_point,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,zv,h)

      dudxi0 = -dudxi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudyi0 = -dudyi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudzi0 = -dudzi0*FOURPI*eps0*AMBER_ELECTROSTATIC

      ! add the coulomb field to get the total E of inner side
      ! converted to electric displacement

      dudxi0 = (dudxi0 + coulomb(1))*epsin/eps0
      dudyi0 = (dudyi0 + coulomb(2))*epsin/eps0
      dudzi0 = (dudzi0 + coulomb(3))*epsin/eps0

      !dt2 = dudxi0*dudxi0 + dudyi0*dudyi0 + dudzi0*dudzi0
      dudni = dudxi0*rn(1)+dudyi0*rn(2)+dudzi0*rn(3)

      !dudno = dudni
      !dudxo0 = (dudxi0 - dudni*rn(1))/epsin*epsout + dudno*rn(1)
      !dudyo0 = (dudyi0 - dudni*rn(2))/epsin*epsout + dudno*rn(2)
      !dudzo0 = (dudzi0 - dudni*rn(3))/epsin*epsout + dudno*rn(3)

      ! part b: apply the normal field approximation
      !         or use the total field

      E2 = dudni
      !E2 = dudxi0*dudxo0 + dudyi0*dudyo0 + dudzi0*dudzo0

      ! compute the DBF as HALF*Q*Dn with the normal field approximation
      ! compute the DBF as HALF*Q*Di*Do/Dn with the total field

      dum = HALF*crg*E2*rn
      !if ( abs(dn) > small_dn ) then
      !   dum = HALF*crg/dn*E2*rn
      !else
      !   dum = HALF*crg*E2*rn
      !end if

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)

   end do

   deallocate (pol_charge)

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz

   if ( verbose .and. eneout ) then
      open (unit = 103, file = 'force.dat')
      write(103,*) ' :::: Atomic qE forces ::::'
      do iatm = 1, natom
         write(103,'(3e20.6)') qex(iatm),qey(iatm),qez(iatm)
      end do
      write(103,*) ' :::: Atomic DB forces ::::'
      do iatm = 1, natom
         write(103,'(3e20.6)') dbx(iatm),dby(iatm),dbz(iatm)
      end do
   end if

   f(1,1:natom) = f(1,1:natom) + dbx(1:natom) + qex(1:natom)
   f(2,1:natom) = f(2,1:natom) + dby(1:natom) + qey(1:natom)
   f(3,1:natom) = f(3,1:natom) + dbz(1:natom) + qez(1:natom)

   eelrf = HALF*eelrf

   if ( verbose .and. eneout ) then
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eelrf
      write(6, '(1x,a,3f12.4)') 'Total solvation force', sum(dbx(1:natom))+sum(qex(1:natom)),&
         sum(dby(1:natom))+sum(qey(1:natom)), sum(dbz(1:natom))+sum(qez(1:natom))
   end if

end subroutine pb_dbfrc_crg2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric boundary forces, third field-based strategy,
!+ loop over boundary grid points, one loop
subroutine pb_dbfrc_fld3( verbose,eneout,natom,eelrf,f,epsx,epsy,epsz,zv,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd
   implicit none

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ eelrf
   _REAL_ f(3,natom)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)
   _REAL_ zv(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)

   ! local variables

   integer, parameter :: n_point = 7
   integer i, j, k, iatm, ip
   _REAL_ g(3), x(3), crd(3), dx(3), rn(3), sgn, rdx, rsphere
   _REAL_ fx0, fy0, fz0, rh
   _REAL_ ds, srfarea, epsth, repsin, repsout, repsp(3), repsm(3)
   _REAL_ E2, fbnd, dum(3)
   _REAL_ srfcrg, crg
   _REAL_ up, dudxi0, dudyi0, dudzi0
   _REAL_ dudxo0, dudyo0, dudzo0
   _REAL_ dudni, dudno, coulomb(3)
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_, allocatable :: pol_charge(:)

   ! initialization

   rh = ONE/h
   repsin = ONE/epsin; repsout = ONE/epsout
   epsth = TWO/(repsin+repsout)
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO
   x = ZERO; sgn = ONE; rsphere = ZERO; crg = ZERO; coulomb = ZERO

   ! compute polarization charges on the boundary grid points and
   ! report total in srfcrg
   ! pol_charge is in au

   allocate (pol_charge(1:nbnd))
   call get_charge_pol(nbnd,xm,ym,zm,phi,cphi,chgrd,pol_charge,srfcrg)

   eelrf = ZERO
   srfarea = ZERO
   do ip = 1, nbnd

      ! collect boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
      crg = pol_charge(ip)*AMBER_ELECTROSTATIC

      ! project the boundary grid point onto the molecular surface, crd() is the
      ! projected coord, x() is the atom/probe center coord, and fx/y/z0 is the grid
      ! version of crd()

      if ( iatm == 0 ) then
         write(6,'(a)') 'PB Bomb in pb_dbfrc_fld3(): can not find projection atom/probe'
         call mexit(6, 1)
      end if
      if      ( iatm > 0 ) then
         ! the contact boundary grid points are projected to the atom spheres
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         sgn = ONE
      else if ( iatm < 0 ) then
         ! the reentry boundary grid points are projected to the solvent spheres
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         sgn = -ONE
      else
         ! otherwise do not project. The setting should make crd = g
         x(1:3) = g(1:3)
         rsphere = ONE
         sgn = ONE
      end if
      dx = g - x; rdx = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
      crd = x + dx*rdx*rsphere

      ! normal direction vector

      rn = sgn*dx*rdx

      ! grid unit of the projected boundary point

      fx0 = (crd(1)-gox)*rh
      fy0 = (crd(2)-goy)*rh
      fz0 = (crd(3)-goz)*rh

      ! compute the effective surface element according to
      ! Cai, et al. JCC, in prep.
      ! first manipulate dielectric constants if smoothing is used

      if ( smoothopt == 0 ) then
         ds = dx(1)/epsx(i,j,k) - dx(1)/epsx(i-1,j,k) + &
              dx(2)/epsy(i,j,k) - dx(2)/epsy(i,j-1,k) + &
              dx(3)/epsz(i,j,k) - dx(3)/epsz(i,j,k-1)
      else
         repsp = repsin
         repsm = repsin
         if ( epsx(i-1,j,k) > epsth ) repsm(1) = repsout
         if ( epsx(i  ,j,k) > epsth ) repsp(1) = repsout
         if ( epsy(i,j-1,k) > epsth ) repsm(2) = repsout
         if ( epsy(i,j ,k ) > epsth ) repsp(2) = repsout
         if ( epsz(i,j,k-1) > epsth ) repsm(3) = repsout
         if ( epsz(i,j,k  ) > epsth ) repsp(3) = repsout
         ds = dx(1)*repsp(1) - dx(1)*repsm(1) + &
              dx(2)*repsp(2) - dx(2)*repsm(2) + &
              dx(3)*repsp(3) - dx(3)*repsm(3)
      end if
      ds = sgn*ds*rdx*h*h/(repsout-repsin)

      srfarea = srfarea + ds

      ! compute reaction field energy, forces, and coulomb field
      ! between the current surface charge and all atomic charges
      ! inner loop over "natom" natoms...

      call get_coulomb(natom,crg,crd,acrg,ax,ay,az,qex,qey,qez,eelrf,coulomb)

      ! now it is time to compute dbf
      ! part a: compute E_in and E_out ...
      !    compute E on the inner side of surface position crd(1:3)

     !call gradu(xm,ym,zm,-ONE,n_point,4,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,zv)
      call onesided(xm,ym,zm,4,n_point,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,zv,h)

      dudxi0 = -dudxi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudyi0 = -dudyi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudzi0 = -dudzi0*FOURPI*eps0*AMBER_ELECTROSTATIC

      !    add the coulomb field to get the total E on the inner side
      !    convert to displacement, D, for consistency with other methods

      dudxi0 = (dudxi0 + coulomb(1))*epsin/eps0
      dudyi0 = (dudyi0 + coulomb(2))*epsin/eps0
      dudzi0 = (dudzi0 + coulomb(3))*epsin/eps0

      !    get normal field component, which is continuous
      !    get D of the outer side based on the jump conditions

      dudni = dudxi0*rn(1)+dudyi0*rn(2)+dudzi0*rn(3)
      dudno = dudni
      !dudxo0 = (dudxi0 - dudni*rn(1))/epsin*epsout + dudno*rn(1)
      !dudyo0 = (dudyi0 - dudni*rn(2))/epsin*epsout + dudno*rn(2)
      !dudzo0 = (dudzi0 - dudni*rn(3))/epsin*epsout + dudno*rn(3)

      ! part b: apply the normal field approximation
      !         or use the total field

      E2 = dudni*dudni
      !E2 = dudxi0*dudxo0 + dudyi0*dudyo0 + dudzi0*dudzo0

      ! part c: compute the surface force element

      fbnd = INV_EIGHTPI*(eps0*(epsin-epsout)/(epsin*epsout))*E2*ds
      dum(1) = fbnd*rn(1)
      dum(2) = fbnd*rn(2)
      dum(3) = fbnd*rn(3)

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)
   end do

   deallocate (pol_charge)

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz

   if ( verbose .and. eneout ) then
      open (unit = 103, file = 'force.dat')
      write(103,*) ' :::: Atomic qE forces ::::'
      do iatm = 1, natom
         write(103,'(3e20.6)') qex(iatm),qey(iatm),qez(iatm)
      end do
      write(103,*) ' :::: Atomic DB forces ::::'
      do iatm = 1, natom
         write(103,'(3e20.6)') dbx(iatm),dby(iatm),dbz(iatm)
      end do
   end if

   f(1,1:natom) = f(1,1:natom) + dbx(1:natom) + qex(1:natom)
   f(2,1:natom) = f(2,1:natom) + dby(1:natom) + qey(1:natom)
   f(3,1:natom) = f(3,1:natom) + dbz(1:natom) + qez(1:natom)

   eelrf = HALF*eelrf

   if ( verbose .and. eneout ) then
      write(6, '(1x,a,f12.4)') 'Total molecular surface', srfarea
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eelrf
      write(6, '(1x,a,3f12.4)') 'Total solvation force', sum(dbx(1:natom))+sum(qex(1:natom)),&
          sum(dby(1:natom))+sum(qey(1:natom)), sum(dbz(1:natom))+sum(qez(1:natom))
   end if

end subroutine pb_dbfrc_fld3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)

   use solvent_accessibility, only : dotarc, arcatm, narcdot, ntri, triatm, triopt
   implicit none

   integer natom, iatm
   _REAL_ dum(3), x(3)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)

   integer, parameter :: mp = 3, np = 3
   integer i, j, iarc
   integer m, n, patm, matm, natm
   _REAL_ mvec(3), nvec(3), mdotn, mxnv(3), rmdotn2, fdotm, fdotn
   _REAL_ dfm, dfn, dum_norm(3), dum_tang(3), dumnorm
   _REAL_ a(mp,np), u(mp,np), w(np), v(mp,np), b(mp), t(np)
   _REAL_ wmax, thresh, TOL
   _REAL_ pvec(3), qvec(3), mdist, ndist, pdist, qdist

   TOL = 1.d-5
   patm = -1; natm = -1; matm = -1
   a = ZERO; u = ZERO; w = ZERO; v = ZERO; b = ZERO; t = ZERO

   ! try different force decomposition scheme

#  include "pb_dbfrc.h"
!#  include "pb_dbfrc1.h"

end subroutine dbfrc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_coulomb(natom,crg,crd,ac,ax,ay,az,fx,fy,fz,eelrf,coulomb)

   implicit none

   ! passed variables

   integer natom
   _REAL_ crg, crd(1:3), coulomb(1:3), eelrf
   _REAL_ ac(natom), ax(natom), ay(natom), az(natom)
   _REAL_ fx(natom), fy(natom), fz(natom)

   ! local variables

   integer jatm
   _REAL_ dinv, d2inv, de, dx(1:3), dff, dcc

   coulomb(1:3) = ZERO
   do jatm = 1, natom
      if ( ligand ) then
         if (liveflag(jatm) == 0) cycle
      endif
      !if ( ligand .and. liveflag(jatm) == 0 ) cycle
      dx(1) = crd(1) - ax(jatm)
      dx(2) = crd(2) - ay(jatm)
      dx(3) = crd(3) - az(jatm)
      d2inv = ONE/(dx(1)**2 + dx(2)**2 + dx(3)**2); dinv = sqrt(d2inv)

      de = ac(jatm)*dinv*AMBER_ELECTROSTATIC
      dcc = de*d2inv
      dff = crg*dcc

      ! calculate reaction field energy

      eelrf = eelrf + de*crg

      ! calculate the QE force on the atom

      fx(jatm) = fx(jatm) - dx(1)*dff
      fy(jatm) = fy(jatm) - dx(2)*dff
      fz(jatm) = fz(jatm) - dx(3)*dff

      ! calculate the coulomb field on the surface charge

      coulomb(1) = coulomb(1) + dx(1)*dcc
      coulomb(2) = coulomb(2) + dx(2)*dcc
      coulomb(3) = coulomb(3) + dx(3)*dcc
   end do


end subroutine get_coulomb

end subroutine pb_fdfrc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ compute induced polarization charges
subroutine get_charge_pol( nbnd,xm,ym,zm,phi,cphi,chgrd,pol_charge,srfcrg )

   use poisson_boltzmann, only : h, eps0, epsin, iepsav
   implicit none

#  include "pb_constants.h"

   ! passed variables

   integer nbnd
   integer xm,ym,zm
   _REAL_ phi(xm,ym,zm), cphi(xm,ym,zm), chgrd(xm,ym,zm)
   _REAL_ pol_charge(nbnd), srfcrg

   ! local variables

   integer i, j, k, ip
   _REAL_ total_phi(7)

   srfcrg = ZERO
   do ip = 1, nbnd
      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip)

      total_phi(1) = phi(i,  j,k) + cphi(i,  j,k)
      total_phi(2) = phi(i-1,j,k) + cphi(i-1,j,k)
      total_phi(3) = phi(i+1,j,k) + cphi(i+1,j,k)
      total_phi(4) = phi(i,j-1,k) + cphi(i,j-1,k)
      total_phi(5) = phi(i,j+1,k) + cphi(i,j+1,k)
      total_phi(6) = phi(i,j,k-1) + cphi(i,j,k-1)
      total_phi(7) = phi(i,j,k+1) + cphi(i,j,k+1)

      pol_charge(ip) = (SIX*h*eps0)*( total_phi(1)-&
           SIXTH*( total_phi(2)+total_phi(3)+&
                   total_phi(4)+total_phi(5)+&
                   total_phi(6)+total_phi(7) ) )
      pol_charge(ip) = ( pol_charge(ip) - chgrd(i,j,k)/(epsin/eps0) )
      srfcrg = srfcrg + pol_charge(ip)
   end do

end subroutine get_charge_pol
