! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)
#include "pb_def.h"
#include "timer.h"

module poisson_boltzmann

   implicit none

#  include "pb_constants.h"

   ! PBMD parameters

   _REAL_, parameter :: pbkb   = 1.3807D-23 / 1.6606D-27 / (1.00D+12)**2 * (1.00D+10)**2
   _REAL_, parameter :: fioni  = 6.0220D+23 / 1.00D+30
   _REAL_, parameter :: fiono  = ONE / fioni
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   _REAL_, parameter :: frcfac = FOURPI*eps0*AMBER_ELECTROSTATIC2

   ! PBMD FD control variables

   logical :: outphi
   logical :: srsas
   logical :: scalerf
   logical :: outlvlset
   logical :: outmlvlset

   integer :: phiform
   integer :: saopt
   integer :: sasopt
   integer :: dbfopt
   integer :: eneopt
   integer :: npbopt
   integer :: solvopt
   integer :: frcopt
   integer :: intopt
   integer :: bcopt
   integer :: smoothopt
   integer :: fold16
   integer :: isurfchg
   integer :: membraneopt
   integer :: poretype
   integer :: augtoltype
   integer :: rxm
   integer :: rym
   integer :: rzm
   integer :: xm
   integer :: ym
   integer :: zm
   integer :: xmym
   integer :: xmymzm
   integer :: nbuffer
   integer :: level
   integer :: nfocus
   integer :: fscale
   integer :: maxitn
   integer :: itn
   integer :: m, n
   integer :: savbcopt(MAXLEVEL)
   integer :: levelblock(MAXLEVEL)
   integer :: savxm(MAXLEVEL)
   integer :: savym(MAXLEVEL)
   integer :: savzm(MAXLEVEL)
   integer :: savxo(MAXLEVEL)
   integer :: savyo(MAXLEVEL)
   integer :: savzo(MAXLEVEL)
   integer :: savxmym(MAXLEVEL)
   integer :: savxmymzm(MAXLEVEL)

   _REAL_ :: h
   _REAL_ :: rgox
   _REAL_ :: rgoy
   _REAL_ :: rgoz
   _REAL_ :: gox
   _REAL_ :: goy
   _REAL_ :: goz
   _REAL_ :: fmiccg
   _REAL_ :: fmiccg2
   _REAL_ :: accept
   _REAL_ :: laccept
   _REAL_ :: wsor
   _REAL_ :: lwsor
   _REAL_ :: norm
   _REAL_ :: inorm
   _REAL_ :: xmax
   _REAL_ :: xmin
   _REAL_ :: ymax
   _REAL_ :: ymin
   _REAL_ :: zmax
   _REAL_ :: zmin
   _REAL_ :: gxmax
   _REAL_ :: gxmin
   _REAL_ :: gymax
   _REAL_ :: gymin
   _REAL_ :: gzmax
   _REAL_ :: gzmin
   _REAL_ :: savxbox(MAXLEVEL)
   _REAL_ :: savybox(MAXLEVEL)
   _REAL_ :: savzbox(MAXLEVEL)
   _REAL_ :: cxbox(MAXLEVEL)
   _REAL_ :: cybox(MAXLEVEL)
   _REAL_ :: czbox(MAXLEVEL)
   _REAL_ :: savh(MAXLEVEL)
   _REAL_ :: savgox(MAXLEVEL)
   _REAL_ :: savgoy(MAXLEVEL)
   _REAL_ :: savgoz(MAXLEVEL)
   _REAL_ :: offx
   _REAL_ :: offy
   _REAL_ :: offz
   _REAL_ :: fillratio

   _REAL_ :: epsin
   _REAL_ :: epsout
   _REAL_ :: epsmem
   _REAL_ :: epsmemb
   _REAL_ :: pbkappa
   _REAL_ :: istrng
   _REAL_ :: ivalence
   _REAL_ :: pbtemp
   _REAL_ :: totcrg
   _REAL_ :: totcrgp
   _REAL_ :: totcrgn

   _REAL_ :: pbgamma_int
   _REAL_ :: pbgamma_ext

   _REAL_ :: mzmin
   _REAL_ :: mzmax
   _REAL_ :: mthick
   _REAL_ :: mctrdz
   _REAL_ :: poreradi

   _REAL_ :: augctf
   _REAL_ :: augtol

   ! PBMD topology information

   logical              :: nocharge
   logical              :: noradius

   integer              :: lastp
   integer              :: ngrdcrg

   integer, allocatable ::   icrd(:,:)
   integer, allocatable :: grdcrg(:,:)
   _REAL_, allocatable :: qgrdcrg(:)
   _REAL_, allocatable ::    gcrd(:,:)
   _REAL_, allocatable ::    acrd(:,:)
   _REAL_, allocatable ::    acrg(:)
   _REAL_, allocatable ::    gcrg(:,:)

   ! PBMD nblist information

   integer              :: maxnbr
   integer              :: maxnba
   _REAL_              :: cutres, cutnb, cutfd, cutsa

   integer, allocatable ::   nshrt(:)
   integer, allocatable ::     nex(:)
   integer, allocatable ::     iex(:,:)
   integer, allocatable :: iprshrt(:)
   integer, allocatable ::  iar1pb(:,:)
   _REAL_, allocatable :: cn1pb(:)
   _REAL_, allocatable :: cn2pb(:)
   _REAL_, allocatable :: cn3pb(:)

   ! PBMD cap water simulation information

   integer              :: mpopt
   integer              :: lmax
   integer              :: inatm
   integer              :: outwat
   integer              :: oution
   integer, allocatable :: outflag(:)
   integer, allocatable :: outflagorig(:)
   integer, allocatable :: mapout(:)
   integer, allocatable :: ibelly(:)
   _REAL_              :: sepbuf

   ! physical variables for energy and force calculations

   integer:: nbnd
   integer:: nbndx
   integer:: nbndy
   integer:: nbndz
   _REAL_, allocatable :: pos_crg(:,:,:)
   _REAL_, allocatable :: surf_crg(:,:)
   integer, allocatable :: ipos_crg(:,:,:)
   integer, allocatable :: crg_num(:)

   ! physical variable maps for numerical solutions

   _REAL_, allocatable ::     phi(:)
   _REAL_, allocatable ::      bv(:)
   _REAL_, allocatable ::   chgrd(:)
   _REAL_, allocatable ::    epsx(:)
   _REAL_, allocatable ::    epsy(:)
   _REAL_, allocatable ::    epsz(:)
   _REAL_, allocatable :: saltgrd(:)
   _REAL_, allocatable ::  ioncrg(:)

   ! geometry maps for dielectric interface

   integer, allocatable ::  insas(:)
   integer, allocatable :: atmsas(:)
   _REAL_, allocatable ::  lvlset(:)
   _REAL_, allocatable :: mlvlset(:)
   _REAL_, allocatable ::      zv(:)

   ! physical variable maps for FD force calculations

   _REAL_, allocatable ::     cphi(:)
   integer, allocatable ::  iepsav(:,:)
   integer, allocatable :: iepsavx(:,:)
   integer, allocatable :: iepsavz(:,:)
   integer, allocatable :: iepsavy(:,:)
   _REAL_, allocatable :: fedgex(:)
   _REAL_, allocatable :: fedgey(:)
   _REAL_, allocatable :: fedgez(:)

   ! saved phi array for pbmd

   _REAL_, allocatable :: xs(:)
   integer :: xsoffset

   ! ligand focusing options

   logical :: ligand
   character(len=256) ligandmask
   integer, allocatable :: liveflag(:)
   integer, allocatable :: realflag(:)
   integer :: ntrajmol
   _REAL_ :: buffer

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver of PBMD energy and forces
!  call pb_force( natom,nres,ntypes,npdec,ix(i02),ix(i04),ix(i06),&
!                 ix(i10),cn1,cn2,xx(l15),x,f,evdw,eelt,epol)
subroutine pb_force( natom,nres,ntypes,npdec,ipres,iac,ico,natex,cn1,cn2,cg,x,f,enb,eel,eelrf )
    
   use solvent_accessibility, only : dprob, radi, radip, radip2, radip3, nzratm, sa_init, sa_driver, sa_free, sa_free_mb

   use decomp, only : irespw, jgroup
   use pbtimer_module
   use random, only: amrset, amrand

   ! Common variables

#  include "../include/md.h"
#  include "pb_md.h"
#ifdef SANDER
#  include "../sander/box.h"
   integer, parameter :: mytaskid = 0
   integer, parameter :: numtasks = 1
#elif defined LIBPBSA
#  include "box.h"
   integer, parameter :: mytaskid = 0
   integer, parameter :: numtasks = 1
#else
#  include "box.h"
   integer, parameter :: mytaskid = 0
   integer, parameter :: numtasks = 1
#endif /* SANDER */
#  include "extra.h"

   ! Passed variables

   integer natom, nres, ntypes, npdec, ipres(*), iac(*), ico(*), natex(*)
   _REAL_ cn1(*), cn2(*), cg(natom), x(3,natom), f(3,natom)
   _REAL_ enb, eel, eelrf

   ! Local variables

   integer i, j, k
   integer iatm, proatm, atmfirst, atmlast
   integer atmind(natom)
   _REAL_ acg(natom)
   _REAL_ pbcutcap, pbxcap, pbycap, pbzcap
   _REAL_ eelrffd, eelrfmp
   _REAL_ pbfrc(3,natom)

   ! Variables initialization

   enb = ZERO; eel = ZERO; eelrf = ZERO
   eelrffd = ZERO; eelrfmp = ZERO
   pbfrc = ZERO
   atmind = 0

   acg = ZERO !making sure client(s) have clean acg

   if ( ifcap /= 0 .and. (ifcap < 3 .or. ifcap > 5) ) then
      pbcutcap = cutcap+TWO; pbxcap = xcap; pbycap = ycap; pbzcap = zcap
      radi(0) = pbcutcap; acrd(1,0) = pbxcap; acrd(2,0) = pbycap; acrd(3,0) = pbzcap
      radip3(1) = radi(0); nzratm(1) = 0
   else
      pbcutcap = ZERO; pbxcap = ZERO; pbycap = ZERO; pbzcap = ZERO
   end if

   ! part a.
   ! split atoms into internal/external and update nblist whenever a new pb grid
   ! is set up if requested

   call pbtimer_start(PBTIME_PBLIST)
   if ( pbgrid ) then
      if ( mpopt == 1 ) then
         ! multipole expansion
         call pb_atmpart(pbverbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                         inatm,outwat,oution,ipres,outflag, &
                         pbxcap,pbycap,pbzcap,pbcutcap,sepbuf,x,ifcap)
      else if ( ifcap == 2 ) then
         ! Use cutcap here, not pbcutcap because the latter is augmented by TWO
         call pb_atmpart(pbverbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                         inatm,outwat,oution,ipres,outflag, &
                         pbxcap,pbycap,pbzcap,cutcap,0.0d0,x,ifcap)
      else if ( ifcap == 5 ) then
         ! Use cutcap here, not pbcutcap because the latter is augmented by TWO
         call pb_atmpart2(pbverbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                          inatm,outwat,oution,ipres,outflag, &
                          cutcap,x)
      else if ( ligand ) then
         ! This option will be visited if it is not the first time pb_force got
         ! called, however there would be something to do here once we've done
         ! MD.
         continue
      else
         ! Multiblock and other conditions go here
         outflag = 0
      end if
   else
      if ( mpopt == 1 .or. ifcap == 2 .or. ifcap == 5 ) outflag = outflagorig
   end if

   if ( pqropt == 0 ) call pb_atmconv(mpopt,ifcap,natom,ibgwat,ienwat,ibgion,ienion,atmind,ipres,x,cg,acg)

   ! part b. Set up nblist and grid
   ! This is for the global run for the coarse grid
   ! If ligand/multiple block is used, these will be updated later in docklist

   if ( ntnba == 1 .and. max(cutnb,cutsa,cutfd) > ZERO ) &
      call pb_atmlist(pbverbose,pbprint,pqropt,&
      maxnba,natom,ntypes,iac,ico,natex,nshrt,nex,iex,iar1pb,iprshrt,&
      cutnb,cutsa,cutfd,cn1,cn2,cn1pb,cn2pb,cn3pb,cg,acrd(1,1))
   if ( ntnbr == 1 ) ntnbr = 0
   if ( ntnba == 1 ) ntnba = 0
   call pbtimer_stop(PBTIME_PBLIST)

   call pbtimer_start(PBTIME_PBSETUP)
   if ( mpopt /=2 .and. pbgrid ) then
      if ( ligand ) &
         call pb_atmpart3(pbverbose,pbprint,natom,buffer,xmin,xmax,ymin,ymax,&
              zmin,zmax,liveflag,realflag,outflag,x)
      if (ifcap == 2 .or. ifcap == 5) then
         call pb_setgrd(ipb,pbverbose,pbprint,pbinit,pbgrid,ifcap,1,inatm,pbxcap,pbycap,pbzcap,pbcutcap)
      else
         call pb_setgrd(ipb,pbverbose,pbprint,pbinit,pbgrid,ifcap,1,natom,pbxcap,pbycap,pbzcap,pbcutcap)
      end if
   end if
   call pbtimer_stop(PBTIME_PBSETUP)

   ! R. Luo
   ! Moving this to pb_exmol() in pb_bldsys.F90 for general applications

   ! part c. compute grid-independent sas calculations for dielectric assignment
   ! when ifcap /= 0, no need to comptue sas

   call pbtimer_start(PBTIME_PBSAS)
   if ( ligand .and. ( ifcap == 0 .or. ifcap == 5 ) ) then
      if( ifcap == 5 ) then
         call sa_init(pbverbose,pbprint,natom,inatm,ifcap,dprob,radi,radip,radip2,outflag)
         ! the call here requires verification if we need to take care of ligand option as well
         call sa_driver(pbverbose,pbprint,pqropt,sasopt,inp,natom,inatm,dosas,ndosas,npbstep,nsaslag,&
                        ligand, outflag,&
                        acrd(1,1),iar1pb(1,0),iprshrt,nex,iex,.false.)
      else
         call sa_init(pbverbose,pbprint,natom,natom,ifcap,dprob,radi,radip,radip2,outflag)
         call sa_driver(pbverbose,pbprint,pqropt,sasopt,inp,natom,natom,dosas,ndosas,npbstep,nsaslag,&
                        ligand, outflag,&
                        acrd(1,1),iar1pb(1,0),iprshrt,nex,iex,.false.)
      end if
   end if
   call pbtimer_stop(PBTIME_PBSAS)

   ! for focussing run, liveflag, outflag, realflag are updated
   ! atom list is updated next
   ! surface area is then updated

   if ( ligand ) then
      call pb_atmlist(pbverbose,pbprint,pqropt,maxnba,natom,ntypes,iac,ico,natex, &
              nshrt,nex,iex,iar1pb,iprshrt,cutnb,cutsa,cutfd,cn1,cn2,cn1pb, &
              cn2pb,cn3pb,cg,acrd)
      call sa_driver(pbverbose,pbprint,pqropt,sasopt,inp,natom,natom,dosas,ndosas, &
              npbstep,nsaslag, ligand                 ,outflag,&
              acrd(1,1),iar1pb(1,0),iprshrt,nex,iex,.false.)
   end if

   ! part e. add FD reaction field energy/force
   ! it is possible the call is without charges as in BAR runs, so don't do any
   ! PB in this case.

   call pbtimer_start(PBTIME_PBFDFRC)
   if ( .not. nocharge .and. .not. noradius ) then
      if ( epsout /= epsin .and. mpopt /= 2 ) then
         ! In the case of ifcap == 2,5, only map crg within cap to grid (atmlast == inatm),
         ! else map all crg (atmlast == natom)
         if( ifcap == 2 .or. ifcap == 5 ) then
            call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,inp,imin,natom,inatm,npdec,idecomp,irespw, &
                    ipres,jgroup,ibgwat,ibgion,pbfrc,enb,eelrffd,npbstep,npbgrid,nstlim)
         else
            call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,inp,imin,natom,natom,npdec,idecomp,irespw, &
                    ipres,jgroup,ibgwat,ibgion,pbfrc,enb,eelrffd,npbstep,npbgrid,nstlim)
         end if
      else if (epsout == epsin) then
         ! Apparently this is for expert use only
            call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,inp,imin,natom,natom,npdec,idecomp,irespw, &
                    ipres,jgroup,ibgwat,ibgion,pbfrc,enb,eelrffd,npbstep,npbgrid,nstlim)
      else
         write(6,'(a)') 'PB Bomb in pb_force(): Unknown FD force call option'
         call mexit(6,1)
      end if
   else
      enb = ZERO ! this is from the p3m routines, which apparently can't work without charges or radii
      eelrffd = ZERO
   end if
   call pbtimer_stop(PBTIME_PBFDFRC)

   ! clean up for sas calculations

   call pbtimer_start(PBTIME_PBSETUP)
   if ( srsas .and. (ifcap == 0 .or. ifcap == 5) ) then
      call sa_free
      if ( npdec > 1 ) dosas = ndosas ! restore dosas to default for later
   end if
   call pbtimer_stop(PBTIME_PBSETUP)

   if ( saopt < 0 ) return

   ! part f. add MP reaction field energy/forces when ifcap /= 0

   call pbtimer_start(PBTIME_PBMP)
   if ( mpopt /= 0 .and. epsout /= epsin ) then
      if ( mpopt == 1 ) then      ! multipole expansion for boundary atoms
         atmfirst = inatm + 1
         atmlast  = natom
      else if ( mpopt == 2 ) then ! multipole expansion for all atoms
         atmfirst = 1
         atmlast  = natom
      end if
      call pb_mpfrc(natom,atmfirst,atmlast,lmax,pbcutcap,pbxcap,pbycap,pbzcap,&
              epsin,epsout,acrg,acrd(1,1),pbfrc,eelrfmp)
   end if
   call pbtimer_stop(PBTIME_PBMP)

   ! part g. add direct coulombic and nonbonded forces

   call pbtimer_start(PBTIME_PBDIRECT)
   if ( eneopt == 4) then
      ! for the new p3m, pb_direct SHOULD BE DONE inside the loop
      continue
   else if ( pqropt == 0 .and. cutnb == ZERO ) then
      call pb_directnocut(natom,proatm,inatm,ipres,ibgwat,ienwat,ibgion,ienion,ntypes,eneopt,idecomp,ifcap, &
              iac,ico,nex,iex,cn1,cn2,acg,acrd(1,1),pbfrc,eel,enb)
   else
      if( ifcap == 2 .or. ifcap == 5) then
         call pb_directwtcut(natom,inatm,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1), &
                 pbfrc,eel,enb)
      else if ( pqropt == 0 ) then
         call pb_directwtcut(natom,natom,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1), &
                 pbfrc,eel,enb)
      end if
   end if
   call pbtimer_stop(PBTIME_PBDIRECT)

   ! part h. returning:
   ! i. returning energies

   eel = eel * eps0/epsin
   if ( eneopt == 1 .and. (bcopt /= 2 .and. bcopt < 6 .or. bcopt > 9) ) then
      eel = eel + eelrffd + eelrfmp
      eelrf = ZERO
   else if ( eneopt == 4 ) then
      eel = eelrffd
      eelrf = ZERO
   else
      eelrf = eelrffd + eelrfmp
   end if

   ! ii. returning forces
   ! only do this for normal amber file input

   if ( pqropt == 0 ) then
      if ( ifcap == 2 .or. ifcap == 5 ) then
         atmlast = inatm
      else
         atmlast = natom
      end if

      do iatm = 1, natom
         f(1,iatm) = f(1,iatm) + pbfrc(1,mapout(iatm))
         f(2,iatm) = f(2,iatm) + pbfrc(2,mapout(iatm))
         f(3,iatm) = f(3,iatm) + pbfrc(3,mapout(iatm))
      end do
   end if

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ convert passed coordinates and charges to the internal format
subroutine pb_atmconv( mpopt,ifcap,natom,ibgwat,ienwat,ibgion,ienion,atmind,ipres,x,cg,acg )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Passed variables

   integer mpopt, ifcap, natom, ibgwat, ienwat, ibgion, ienion
   integer atmind(natom), ipres(*)
   _REAL_ x(3,natom), cg(natom), acg(natom)

   ! Local variables

   integer i, j, ifirst, ilast, iatm, ires, num

   if ( mpopt == 1 .or. ifcap == 2 .or. ifcap == 5 ) then

      ! copy reordered coord/charge to private arrays for pb/mp or ifcap == 2,5
      ! protein atoms go into internal portion (so do IONS!)

      ifirst = 1; ilast = ipres(ibgwat)-1

      do iatm = ifirst, ilast
         acrd(1,iatm) = x(1,iatm); acrd(2,iatm) = x(2,iatm); acrd(3,iatm) = x(3,iatm)
         acrg(iatm) = cg(iatm)/18.2223d0; acg(iatm) = cg(iatm); atmind(iatm) = iatm
         mapout(iatm) = iatm
      end do

      ! water atoms go into internal/external portion

      ifirst = ipres(ibgwat); ilast = natom

      i = ifirst; j =  inatm + 1
      do iatm = ifirst, ilast
         if ( outflag(iatm) == 0 ) then
            acrd(1,i   ) = x(1,iatm); acrd(2,i   ) = x(2,iatm); acrd(3,i   ) = x(3,iatm)
            acrg(i   ) = cg(iatm)/18.2223d0; acg(i   ) = cg(iatm); atmind(i   ) = iatm
            mapout(iatm) = i
            i = i + 1
         else
            acrd(1,j   ) = x(1,iatm); acrd(2,j   ) = x(2,iatm); acrd(3,j   ) = x(3,iatm);
            acrg(j   ) = cg(iatm)/18.2223d0; acg(j   ) = cg(iatm); atmind(j   ) = iatm
            mapout(iatm) = j
            j = j + 1
         end if
      end do

      ! store original outflag array and prepare an updated one for water atoms

      outflagorig = outflag
      do iatm = ifirst, ilast
         if( iatm <= inatm ) then
            outflag(iatm) = 0
         else
            outflag(iatm) = 1
         end if
      end do

   else

      do iatm = 1, natom
         acrd(1,iatm) = x(1,iatm); acrd(2,iatm) = x(2,iatm); acrd(3,iatm) = x(3,iatm)
         acrg(iatm) = cg(iatm)/18.2223d0; acg(iatm) = cg(iatm)
         mapout(iatm) = iatm
      end do
   end if

end subroutine pb_atmconv

end subroutine pb_force

end module poisson_boltzmann
