#include "copyright.h"
#include "../include/dprec.fh"

#ifndef SANDER
#ifndef LIBPBSA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open input files and read cntrl namelist.
subroutine mdread1()
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  Meng-Juei Hsieh
   !  The Luo Research Group
   !  University of California, Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   use memory_module

   implicit none

#  include "pb_constants.h"
#  include "files.h"
#  include "box.h"
#  include "../include/md.h"

   integer     ifind
   logical     mdin_cntrl  ! true if this namelist exists in mdin
   character(len=8) date
   character(len=10) time
   _REAL_ rnum

   namelist /cntrl/ ntx,    ipb,    inp,   igb,            &
                    imin,   ntf,    ntb,   dielc,  cut,    & !compatibility
                    nsnb,   maxcyc, ntmin,                 & !compatibility
                    ivcap,  cutcap, xcap,  ycap,   zcap,   & !compatibility
                    idecomp, ig

   if (mdout /= "stdout" ) call myopen(6,mdout,owrite,'F','W')

   call myopen(5,mdin,'O','F','R')

   write(6,9308)

   call date_and_time( DATE=date, TIME=time )

   write(6,'(12(a))') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)

   if (owrite /= 'N') write(6, '(2x,a)') '[-O]verwriting output'

   ! Echo the file assignments to the user:

   write(6,9700) 'MDIN'   ,mdin(1:70)  , 'MDOUT' ,mdout(1:70) , &
                 'INPCRD' ,inpcrd(1:70), 'PARM'  ,parm(1:70)

   ! Echo the input file to the user:

   call myechoin(5,6)

   !     ----- READ DATA CHARACTERIZING THE MD-RUN -----

   read(5,'(20a4)') title

   ! read input in namelist format, first setting up defaults

   imin = 1
   ntx = 1
   ipb = 2
   inp = 2
   igb = 0

   ! the following are to be retired, but included for backward compatibility

   irest = 0
   ibelly = 0
   idecomp= 0
   ntxo = 1
   ig = 71277
!  tempi = ZERO
   ntt = 0
!  temp0 = 300.0d0

   tautp = ONE
   ntp = 0
   pres0 = ONE
   comp = 44.6d0
   taup = ONE
   npscal = 1
   nscm = 1000
   nstlim = 1
   t = ZERO
   dt = 0.001d0
   ntc = 1
   tol = 0.00001
   ntpr = 50
   ntwr = 500
   ntwx = 0
   ntwv = 0
   ntwe = 0
   ntave = 0
   ioutfm = 0
   ntr = 0
   ntrx = 1
   fcap = 1.5d0

   isftrp = 0
   rwell = ONE
   dx0 = 0.01d0
   drms = 1.0d-4
   vlimit = 20.0d0
   mxsub = 1
   ntwprt = 0
!   plevel = 1
   rgbmax = 25.d0
   iyammp = 0
   vrand=1000
   iwrap = 0
   nrespa = 1
   nrespai = 1
   irespa = 1
   gamma_ln = ZERO
   iconstreff = 0
   cut_inner = EIGHT
   icfe = 0
   clambda = ZERO
   klambda = 1
   rbornstat = 0
   lastrst = 2000000
   lastist = 2000000
   restraintmask=''
   restraint_wt = ZERO
   bellymask=''
   saltcon= ZERO
   surften= 0.005d0
   offset = 0.09d0
   intdiel= ONE
   extdiel= 78.5d0
   gbsa   = 0
   ncyc   = 10

   icnstph = 0
   solvph = SEVEN
   ntcnstph = 10

   ! check what namelists are defined

   mdin_cntrl=.false.
   mdin_pb=.false.
   call mynmlsrc('cntrl',5,ifind)
   if (ifind /= 0) mdin_cntrl=.true.
   call mynmlsrc('pb',5,ifind)
   if (ifind /= 0) mdin_pb=.true.

   ! read cntrl namelist

   rewind 5
   if ( mdin_cntrl ) then
      read(5,nml=cntrl)
   else
      write(6, '(1x,a,/)') 'Error: Could not find cntrl namelist'
      call mexit(6,1)
   end if

   ! disabled cntrl flag

   ntf   =1
   ntb   =1
   dielc  = ONE
   cut    = EIGHT
   nsnb   = 25
   !scnb   = TWO
   !scee   = 1.2d0
   maxcyc = 1
   ntmin  = 1
   ivcap  = 0
   xcap   = 0
   ycap   = 0
   zcap   = 0

   if ( ig == -1 ) then
      call wallclock(rnum)
      ig = nint(rnum)
   end if

   ! read pb namelist

   if (ipb <= 0) then
      write(6,'(a)') "Error: ipb can only be > 0 in PBSA"
      call mexit(6,1)
   else
      call pb_read
   endif

   write(6,9309)

   call printflags()

   9308 format(/10x,55('-'),/10x, &
         'PBSA VERSION 2018                       UC Irvine', &
         /10x,55('-')/)
   9309 format(/80('-')/'   1.  RESOURCE   USE: ',/80('-')/)
   9700 format(/,'File Assignments:',/,12('|',a6,': ',a,/))

end subroutine mdread1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize to defaults and print the inputable variables.
subroutine mdread2(x,ix,ih)

   use decomp, only: jgroup, index, irespw

   use memory_module

   implicit none

   _REAL_ x(*)
   integer ix(*)
   character(len=4) ih(*)
   integer inerr
   integer nmropt
   integer ntmp, ngrp

#  include "pb_constants.h"
#  include "files.h"
#  include "box.h"
#  include "../include/md.h"
#  include "pb_md.h"

   close(unit=8)

   !     ----- SET THE DEFAULT VALUES FOR SOME VARIABLES -----

   nrp = natom

   if (ifbox == 1) write(6, '(/5x,''BOX TYPE: RECTILINEAR'')')
   if (ifbox == 2) write(6, '(/5x,''BOX TYPE: TRUNCATED OCTAHEDRON'')')
   if (ifbox == 3) write(6, '(/5x,''BOX TYPE: GENERAL'')')

   nsolut =  nrp
   if ( nscm > 0 .and. ntb == 0 ) then
      ndfmin = 6   ! both translation and rotation com motion removed
      if (nsolut == 1) ndfmin = 3
      if (nsolut == 2) ndfmin = 5
   else if ( nscm > 0 ) then
      ndfmin = 3    ! just translation com will be removed
   else
      ndfmin = 0
   end if
   if (ibelly > 0) then   ! No COM Motion Removal, ever.
      nscm = 9999999
      ndfmin = 0
   end if
   if(nscm <= 0) nscm = 9999999
   init = 3
   if (irest > 0) then
!     init = 4
      write(6,'(a)') "We are sorry, but irest have to be 0"
      call mexit(6,0)
   endif
   !if (scnb == ZERO ) scnb = TWO
   if (dielc <= ZERO ) dielc = ONE
   if (tautp <= ZERO ) tautp = 0.2d0
   if (taup <= ZERO ) taup = 0.2d0

   ! RESET THE CAP IF NEEDED

   if(ivcap == 2) ifcap = 0

   ! PRINT DATA CHARACTERIZING THE RUN

   nmropt = 0
   iesp = 0
   ipol = 0
   write(6,9328)
   write(6,9008) title
   write(6,'(/a)') 'General flags:'
   write(6,'(5x,2(a,i8))') 'imin    =',imin,', nmropt  =',nmropt

   write(6,'(/a)') 'Nature and format of input:'
   write(6,'(5x,4(a,i8))') 'ntx     =',ntx,', irest   =',irest, &
         ', ntrx    =',ntrx

   write(6,'(/a)') 'Nature and format of output:'
   write(6,'(5x,4(a,i8))') 'ntxo    =',ntxo,', ntpr    =',ntpr, &
         ', ntrx    =',ntrx,', ntwr    =',ntwr
   write(6,'(5x,4(a,i8))') 'iwrap   =',iwrap,', ntwx    =',ntwx, &
         ', ntwv    =',ntwv,', ntwe    =',ntwe
   write(6,'(5x,3(a,i8),a,i7)') 'ioutfm  =',ioutfm, &
         ', ntwprt  =',ntwprt, &
         ', idecomp =',idecomp,', rbornstat=',rbornstat

   write(6,'(/a)') 'Potential function:'
   write(6,'(5x,5(a,i8))') 'ntf     =',ntf,', ntb     =',ntb, &
         ', ipb     =',ipb,', nsnb    =',nsnb
   write(6,'(5x,4(a,i8))') 'ipol    =',ipol,', gbsa    =',gbsa, &
         ', iesp    =',iesp
   write(6,'(5x,3(a,f10.5))') 'dielc   =',dielc, &
         ', cut     =',cut,', intdiel =',intdiel

   !write(6,'(5x,3(a,f10.5))') 'scnb    =',scnb, &
   !      ', scee    =',scee

   write(6,'(/a)') 'Frozen or restrained atoms:'
   write(6,'(5x,4(a,i8))') 'ibelly  =',ibelly,', ntr     =',ntr

   if( imin /= 0 ) then
      write(6,'(/a)') 'Energy minimization:'
      ! print inputable variables applicable to all minimization methods.
      write(6,'(5x,4(a,i8))') 'maxcyc  =',maxcyc,', ncyc    =',ncyc, &
            ', ntmin   =',ntmin
      write(6,'(5x,2(a,f10.5))') 'dx0     =',dx0, ', drms    =',drms

!     ! Input flag ntmin determines the method of minimization
!     select case ( ntmin )
!     case ( 0, 1, 2 )
!        ! no specific output
!     case default
!        ! invalid ntmin
!        write(6,'(a,i3)') ' ERROR: Invalid NTMIN.', ntmin
!        call mexit(6, 1)
!     end select

   else
      write(6,'(/a)') 'Molecular dynamics:'
      write(6,'(5x,4(a,i8))') 'nstlim  =',nstlim,', nscm    =',nscm, &
            ', nrespa  =',nrespa
      write(6,'(5x,3(a,f10.5))') 't       =',t, &
            ', dt      =',dt,', vlimit  =',vlimit
   end if

   cut = cut*cut
   cut_inner = cut_inner*cut_inner

   ! If user has requested Poisson-Boltzmann electrostatics, set up variables

   call pb_init(ifcap,natom,nres,ntypes,nbonh,nbona,ix(i02),ix(i04),ix(i06),ix(i08),ix(i10),&
                ix(iibh),ix(ijbh),ix(iiba),ix(ijba),ix(ibellygp),ih(m02),ih(m04),ih(m06),x(l15),x(l97))

   ! we are ready to read the data from the pqr file

#if !defined SANDER && !defined LIBPBSA
   if ( pqropt /= 0 ) then
         call myopen(PQR_UNIT,pqr,'O','F','R')
         call rdpqr2(PQR_UNIT, natom, ih(m02))
   end if
#endif /*ndef SANDER or LIBPBSA*/

   ntmp = 0; ngrp = 0
   if(idecomp > 0) then
      write(6,9428)
      call rgroup(natom,ntmp,nres,ngrp,ix(i02),ih(m02), &
            ih(m04),ih(m06),ih(m08),ix(ibellygp), &
            jgroup,index,irespw,npdec, &
            x(l60),x(lcrdr),.false.,.false.,.false.,idecomp,5,.true.)
   end if

   ! checks on not supported options

   if( ifcap /= 0 ) then
      write(6, '(a,i3)') ' ERROR: cap water not supported.', ifcap
      call mexit(6, 1)
   end if

   if( igb /= 10 .and. ipb < 0 ) then
      write(6, '(a,i3)') ' ERROR: Non PB potential not supported', igb
      call mexit(6, 1)
   end if

   ! checks on bogus data ---

   inerr = 0

   if (imin /= 0 .and. imin /= 1 .and. imin /= 6) then
      write(6,'(/2x,a,i3,a)') 'IMIN (',imin,') must be 0 or 1 or 6.'
      inerr = 1
   end if
   if (ntx < 1 .or. ntx > 7) then
      write(6,'(/2x,a,i3,a)') 'NTX (',ntx,') must be in 1..7'
      inerr = 1
   end if

   ! WARNINGS:

   if (inerr == 1) then
      write(6, '(/,a)') ' *** input error(s)'
      call mexit(6,1)
   end if

   !  DEBUG input; force checking

   call load_debug(5)

   ! Standard format statements:

   9328 format(/80('-')/,'   2.  CONTROL  DATA  FOR  THE  RUN',/80('-')/)
   9428 format(/4x,'LOADING THE DECOMP ATOMS AS GROUPS',/)
   9008 format(a80)


end subroutine mdread2
#endif /* LIBPBSA */

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit defined preprocessor names, ie, flags.
subroutine printflags()

   implicit none
   integer     max_line_length
   parameter ( max_line_length = 80 )

   character(len=max_line_length) line  ! output string of active flags
   integer n                            ! len(line)

   line = '| Flags:'
   n = 8

#ifdef ISTAR2
   call printflags2(' ISTAR2',7,n,line,.false.)
#endif
#ifdef SGIFFT
   call printflags2(' SGIFFT',7,n,line,.false.)
#endif
#ifdef noBTREE
   call printflags2(' noBTREE',8,n,line,.false.)
#endif
#ifdef NMODE
   call printflags2(' NMODE',6,n,line,.false.)
#endif
#ifdef HAS_10_12
   call printflags2(' HAS_10_12',10,n,line,.false.)
#endif
#ifdef NO_SIGN
   call printflags2(' NO_SIGN',8,n,line,.false.)
#endif
#ifdef CHARMM
   call printflags2(' CHARMM',7,n,line,.false.)
#endif
#ifdef DNA_SHIFT
   call printflags2(' DNA_SHIFT',10,n,line,.false.)
#endif
#ifdef CHARGE_MIXING
   call printflags2(' CHARGE_MIXING',14,n,line,.false.)
#endif
#ifdef NO_EWGRPS
   call printflags2(' NO_EWGRPS',10,n,line,.false.)
#endif

   call printflags2(' ',1,n,line,.true.)
   return
end subroutine printflags

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Primitive pre-Fortran90 implementation of printflags.
subroutine printflags2(flag,flag_len,line_len,line,last)

   implicit none
   integer     max_line_length
   parameter ( max_line_length = 80 )

   character(*) flag                ! flag name with blank prefix, intent(in)
   integer flag_len                 ! len(flag), intent(in)
   integer line_len                 ! len(line), intent(inout)
   character(len=max_line_length) line ! intent(inout)
   logical last                     ! is this the last flag ?, intent(in)

   if (line_len + flag_len > max_line_length) then
      write( 6,'(a)') line
      ! begin another line
      line = '| Flags:'
      line_len=8
   end if
   line=line(1:line_len) // flag(1:flag_len)
   line_len=line_len+flag_len
   if(last)write( 6,'(a)') line
   return
end subroutine printflags2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of a float; abort on illegal values.
subroutine float_legal_range(string,param,lo,hi)
   implicit none
   _REAL_ param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'Ewald PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',e12.5)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',e12.5,' Upper limit: ',e12.5)
   63 format(1x,'Check ew_legal.h')
   return
end subroutine float_legal_range

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer; abort on illegal values.
subroutine int_legal_range(string,param,lo,hi)
   implicit none
   integer param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',i8)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',i8,' Upper limit: ',i8)
   63 format(1x,'The limits may be adjustable; search in the .h files ')
   return
end subroutine int_legal_range

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer option; abort on illegal values.
subroutine opt_legal_range(string,param,lo,hi)
   implicit none
   integer param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'Ewald OPTION CHECKING: ')
   60 format(1x,'option ',a,' has value ',i5)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',i5,' Upper limit: ',i5)
   63 format(1x,'Check the manual')
   return
end subroutine opt_legal_range
#endif /*ifndef SANDER*/

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBSA parsing and initialization
#ifdef LIBPBSA
subroutine pb_read(smoothopt_, radiopt_, npbopt_, solvopt_, maxitn_,   &
   nbuffer_, nfocus_, fscale_, npbgrid_, dbfopt_, bcopt_, scalec_,     &
   eneopt_, frcopt_, nsnbr_, nsnba_, phiout_, phiform_, npbverb_,      &
   npopt_, &
   decompopt_, use_rmin_, use_sav_, maxsph_, maxarc_, ndofd_, ndosas_, &
   mpopt_, lmax_, maxarcdot_, pbprint_,                                &
   epsin_, epsout_, epsmem_, istrng_, pbtemp_, dprob_, iprob_,         &
   accept_, fillratio_, space_, arcres_, cutres_, cutfd_, cutnb_,      &
   sprob_, vprob_, rhow_effect_, cavity_surften_, cavity_offset_,      &
   cutsa_, fmiccg_, ivalence_, laccept_, wsor_, lwsor_, radinc_,       &
   expthresh_, offx_, offy_, offz_, sepbuf_)
#else
subroutine pb_read
#endif

   ! Module variables

   use poisson_boltzmann, only : epsin, epsout, epsmem, smoothopt,   &
         istrng, pbtemp, npbopt, solvopt, accept, maxitn,    &
         fillratio, nbuffer, nfocus, fscale, dbfopt, bcopt,  &
         buffer, xmin, xmax, ymin, ymax, zmin, zmax,         &
         eneopt, frcopt, cutres, cutfd, cutnb, phiform,      &
         cutsa, fmiccg, ivalence, laccept, wsor, lwsor,      &
         pbkappa, offx, offy, offz, sepbuf, mpopt, lmax,     &
         saopt, intopt, ligand, ligandmask, outphi, scalerf, &
         srsas, level, lastp, savxm, savym, savzm, savxmym,  &
         savxmymzm, savh, savbcopt, savgox, savgoy, savgoz,  &
         pbgamma_int, pbgamma_ext, eps0, fioni,pbkb,isurfchg,&
         sasopt,        &
         membraneopt, mthick, mctrdz, outlvlset, outmlvlset, &
         poretype,poreradi,augtoltype, augctf, augtol,     &
         fmiccg2, fold16, rxm, rym, rzm, rgox, rgoy, rgoz, epsmemb

   use dispersion_cavity
   use solvent_accessibility
#ifdef SANDER
   use file_io_dat
#endif
   implicit none

   ! Common variables

#  include "pb_def.h"
#  include "../include/md.h"
#  include "pb_md.h"
#ifndef SANDER
#  include "files.h"
#endif

   ! Local variables

   integer npbverb, l, phiout, scalec
   _REAL_ space

   ! Passed variables if it is used as a C library and it is read by the namelist
   ! if it is standalone or as a F library
   ! WMBS: Need to add lvlset writting, membrane, and aug IIM options. See below
   ! RL: only need to support membrane function and surface area in the libraries.

#ifdef LIBPBSA
   integer smoothopt_, radiopt_, npbopt_, solvopt_, maxitn_,     &
           nbuffer_, nfocus_, fscale_, npbgrid_, dbfopt_, bcopt_,&
           scalec_, eneopt_, frcopt_, nsnbr_, nsnba_, phiout_,   &
           phiform_, npbverb_, npopt_, decompopt_, use_rmin_,    &
           use_sav_, maxsph_, maxarc_, ndofd_, ndosas_, mpopt_,  &
           lmax_, maxarcdot_
   logical pbprint_
   _REAL_  epsin_, epsout_, epsmem_, istrng_, pbtemp_, dprob_, iprob_,    &
           accept_, fillratio_, space_, arcres_, cutres_, cutfd_,&
           cutnb_, sprob_, vprob_, rhow_effect_, cavity_surften_,&
           cavity_offset_, cutsa_, fmiccg_, ivalence_, laccept_, &
           wsor_, lwsor_, radinc_, expthresh_, offx_, offy_,     &
           offz_, sepbuf_

#else

   namelist /pb/ epsin, epsout, epsmem, smoothopt, istrng, pbtemp,   &
      radiopt, dprob, mprob, iprob, npbopt, solvopt, accept, maxitn, &
      fillratio, space, nbuffer, nfocus, fscale, npbgrid,      &
      arcres,dbfopt,bcopt,scalec,eneopt,frcopt,       cutfd,   &
      cutnb,        nsnba,phiout, phiform, npbverb, npopt,     &
      decompopt, use_rmin, sprob, vprob, rhow_effect, use_sav, &
      cavity_surften, cavity_offset, maxsph, maxarc,           &
      cutsa, ndofd, ndosas, fmiccg, ivalence, laccept, wsor,   &
      lwsor, pbkappa, radinc, expthresh, offx, offy, offz,     &
      sepbuf, mpopt, lmax, saopt, intopt, ligandmask, buffer,  &
      xmin, ymin, zmin, xmax, ymax, zmax, isurfchg,            &
      triopt, sasopt, maxtri, maxarcdot, radiscale, protscale, &
      radires,membraneopt,mthick,mctrdz,outlvlset,outmlvlset,  &
      poretype, poreradi, augtoltype, augtol, augctf,          &
      fmiccg2, fold16, rxm, rym, rzm, rgox, rgoy, rgoz, epsmemb
#endif

   ! initialize with default values

   ! output options

   outphi = .false.
   outlvlset = .false.
   outmlvlset = .false.
   phiout = 0
   phiform = 0
   isurfchg = 0

   ! WMBS: membrane options
   ! this needs to be updated for LIBPBSA

   membraneopt = 0
   mthick = 40.0d0
   mctrdz = 0.0d0
   poretype = 1
   poreradi = 6.0

   ! WMBS: Augmented IIM options

   augtoltype = 1
   augctf = 0.0d0
   augtol = 1.0d-5

   ! physical constants

   epsin = 1.0d0
   epsout = 80.0d0
   epsmem= 0.0d0
   epsmemb= 0.0d0
   istrng = 0.0d0
   ivalence = 1.0d0
   pbtemp = 300.0d0

   scalec = 0
   scalerf = .false.

   ! molecular surface options

   sasopt = 0
   triopt = 1
   srsas = .true.
   dprob = 1.4d0
   mprob = 2.7d0 ! optimized basd on membrane channel visualization
   iprob = 2.0d0
   smoothopt = 1
   radiopt = 1
   radiscale = 1.0d0
   protscale = 1.0d0
   radires = ''
   radinc = 0.798d0 ! optimized radinc = dprob*0.57
   expthresh = 0.08d0 ! optimized expthresh = 0.08
   arcres = 0.25d0
   maxsph  = 400
   maxarcdot= 1500
   maxarc = 512
   maxtri = 10

   ! molecular surface area calculation options
   ! this needs to be updated for LIBPBSA

   saopt = 0

   ! nonpolar solvent model options

   npopt = 2
   decompopt = 2
   use_rmin = 1
   use_sav = 1
   rhow_effect = 1.129d0
   sprob = 0.557d0
   vprob = 1.30d0
   cavity_surften = 0.0378d0
   cavity_offset = -0.5692d0

   ! finite-difference options

   fold16 = 0
   nbuffer = 0
   nfocus = 2
   fscale = 8
   level = 1
   bcopt = 5
   space = 0.5d0
   rxm = 0
   rym = 0
   rzm = 0
   rgox = 0.0d0
   rgoy = 0.0d0
   rgoz = 0.0d0
   savxm(1) = 0
   savym(1) = 0
   savzm(1) = 0
   savxmym(1) = 0
   savxmymzm(1) = 0
   savh(1) = 0
   savbcopt(1) = 0
   offx = 0.0d0
   offy = 0.0d0
   offz = 0.0d0
   fillratio= 2.0d0
   xmin = 0.0d0
   xmax = 0.0d0
   ymin = 0.0d0
   ymax = 0.0d0
   zmin = 0.0d0
   zmax = 0.0d0
   buffer = 16*0.5 ! this needs to be updated for NAB

   ! linear system solver options

   npbopt = 0
   solvopt = 1
   fmiccg = 0.90d0
   fmiccg2= 0.90d0
   accept = 1.0d-3
   laccept = 0.1d0
   wsor = 1.9d0
   lwsor = 1.95d0
   maxitn = 100

   ! i/o and md options

   pbverbose = .false.
   pbprint = .true.
   pbgrid = .true.
   pbinit = .true.
   npbverb = 0
   npbgrid = 1
   ndofd = 1
   dofd = 1
   donpsa = .true.
   ndosas = 1
   nsaslag = 100
   nsnbr = 1
   nsnba = 1
   ntnbr = 1
   ntnba = 1
   cutres = 99.9d0
   cutfd = 5.0d0
   cutnb = 0.0d0
   cutsa = 9.0d0
   lastp = 0
   pbgamma_int = 1.0
   pbgamma_ext = 65.0

   ! multipole expansion options

   sepbuf = 4.0d0
   mpopt = 0
   lmax = 80

   ! force interpolation options

   dbfopt = -1
   eneopt = -1
   frcopt = 0
   intopt = 1 ! this needs to be updated for NAB

   ! ligand focusing option

   ligandmask=''

   ! default options for sander

#ifdef SANDER
   if ( imin == 0 .or. maxcyc > 1 ) then
      space = 0.25d0
      arcres = 0.125d0
      smoothopt = 1
      nfocus = 2
      fscale = 4
      eneopt = 2
      bcopt = 6
      frcopt = 2
      if ( imin == 0 ) fmiccg = -0.30d0
   end if
#endif

   ! reading parameters

#ifdef LIBPBSA
   smoothopt=smoothopt_
   radiopt=radiopt_
   npbopt=npbopt_
   solvopt=solvopt_
   maxitn=maxitn_
   nbuffer=nbuffer_
   nfocus=nfocus_
   fscale=fscale_
   npbgrid=npbgrid_
   dbfopt=dbfopt_
   bcopt=bcopt_
   scalec=scalec_
   eneopt=eneopt_
   frcopt=frcopt_
   nsnbr=nsnbr_
   nsnba=nsnba_
   phiout=phiout_
   phiform=phiform_
   npbverb=npbverb_
   npopt=npopt_
   decompopt=decompopt_
   use_rmin=use_rmin_
   use_sav=use_sav_
   maxsph=maxsph_
   maxarc=maxarc_
   ndofd=ndofd_
   ndosas=ndosas_
   mpopt=mpopt_
   lmax=lmax_
   maxarcdot=maxarcdot_
   pbprint=pbprint_
   epsin=epsin_
   epsout=epsout_
   epsmem=epsmem_
   istrng=istrng_
   pbtemp=pbtemp_
   dprob=dprob_
   iprob=iprob_
   accept=accept_
   fillratio=fillratio_
   space=space_
   arcres=arcres_
   cutres=cutres_
   cutfd=cutfd_
   cutnb=cutnb_
   sprob=sprob_
   vprob=vprob_
   rhow_effect=rhow_effect_
   cavity_surften=cavity_surften_
   cavity_offset=cavity_offset_
   cutsa=cutsa_
   fmiccg=fmiccg_
   ivalence=ivalence_
   laccept=laccept_
   wsor=wsor_
   lwsor=lwsor_
   radinc=radinc_
   expthresh=expthresh_
   offx=offx_
   offy=offy_
   offz=offz_
   sepbuf=sepbuf_
#else
   if ( mdin_pb ) then
      rewind 5
      read(5, nml=pb)
   end if
#endif

   ! option checking and post processing of options

   ! WMBS: membrane options compatibility

   if ( membraneopt > 0 ) then
      if ( epsmem == 0.0d0 .and. epsmemb /= 0.0d0 ) epsmem = epsmemb
      if ( epsmem == 0.0d0 .and. epsmemb == 0.0d0 ) epsmem = 1.0d0
      if ( epsmem == 0.0d0 .and. epsmemb == 0.0d0 ) then
         epsmem = 1.0d0
         epsmemb = 1.0d0
      end if
      if ( epsmem < 1.0d0 .or. epsmem > epsout ) then
         write(6,'(a)') ' PB Bomb in pb_read(): membrane eps is currently limited to values between 1.0 and epsout'
         call mexit(6,1)
      end if
      if ( ipb > 2 ) then
         write(6,'(a)') ' PB Bomb in pb_read(): membrane setup is only compatible with ipb=1 or 2'
         call mexit(6,1)
      end if
      if ( sasopt < 2 .and. ipb == 2 ) then
         write(6,'(a)') ' PB Bomb in pb_read(): membrane SAS and SES are only compatible with ipb=1'
         call mexit(6,1)
      end if
      if ( nfocus /= 1 ) then
         write(6,'(a)') ' PB Bomb in pb_read(): membrane setup is only compatible with nonfocussing FDPB'
         call mexit(6,1)
      end if
      if ( bcopt /=10 ) then
         write(6,'(a)') ' PB Bomb in pb_read(): membrane setup is only compatible with the periodic boundary'
         call mexit(6,1)
      end if
      if ( solvopt < 1 .or. solvopt > 4 ) then
         write(6,'(a)') ' PB Info in pb_read(): membrane setup requires a periodic solver'
         solvopt=1
      end if
      if ( eneopt > 1 .or. frcopt > 0 ) then
         write(6,'(a)') ' PB Info in pb_read(): membrane setup is only compatible with eneopt=1 and frcopt=0'
         eneopt = 1; frcopt = 0
      end if
      if ( cutfd <= 12.0*space ) then
         write(6,'(a)') ' PB Info in pb_read(): membrane setup requires a normal cutfd, cutnb'
         cutfd = 14.0*space; cutnb = 99.0
      end if
   else
      epsmem = 1.0d0
      epsmemb = 1.0d0
   end if

   ! surface options

   if ( ipb == 1 .and. sasopt == 2 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): the smooth density surface must be built with a level set funtion (ipb > 1)'
      call mexit(6, 1)
   end if

   ! set up solvent accessible arc resolutions and limits
   ! Mengjuei: if maxarcdot is default, automatically set it up.

   if ( arcres > space * 0.5d0 ) then
      arcres = space * 0.5d0
      write(6,'(a,f12.4)') ' PBSA Warning in pb_read(): arcres is too big and is reset to', arcres
   end if
   if ( maxarcdot == 1500 ) maxarcdot = 1500*nint(0.5d0/arcres)

   ! solvent and ion probes

   if ( dprob > mprob ) then
      write(6,'(a)') ' PBSA Warning in pb_read(): membrane probe is smaller than solvent probe'
   end if
   if ( dprob > iprob ) then
      if ( sasopt == 1 ) then
         write(6,'(a)') ' PB Bomb in pb_read(): ion probe cannot be smaller than solvent prob when the SAS is used (sasopt = 1)'
         call mexit(6, 1)
      else
         write(6,'(a)') ' PBSA Warning in pb_read(): ion probe is smaller than solvent prob'
      end if
   end if
   if ( dprob <= space .and. ipb == 1 ) then
      write(6,'(a)') ' PBSA Warning in pb_read(): setting grid spacing larger than solvent probe'
      write(6,'(a)') 'may cause numerical instability if ipb=1'
   end if
   if ( dprob <= space .and. ipb == 2 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): solvent probe cannot be smaller than grid spacing if ipb=2'
      call mexit(6, 1)
   end if
   if ( dprob < 0.1d0 .and. (sasopt == 0 .or. sasopt == 2) ) then
      write(6,'(a)') ' PB Bomb in pb_read(): solvent probe cannot be too small '
      write(6,'(a,f9.3)') '   for solvent excluded surface or density surface', dprob
      call mexit(6, 1)
   end if

   ! pb options

   if ( ipb < 1 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): ipb can only be larger than or equal to 1.'
      call mexit(6, 1)
   else if ( igb == 10 .and. ipb == 0 ) then
      ipb = 2
   else if ( igb /= 0 .and. igb /= 10 ) then
      igb = 10
   end if

   if ( ipb == 4 .or. ipb == 5 .or. ipb == 6 ) then
      if ( bcopt /=2 .and. bcopt /=6 ) then
         bcopt = 2
         write(6,'(a)') ' PB Info in pb_read(): bcopt has been reset to 2 with ipb=4/5/6'
      end if
      if ( nfocus /= 1 ) then
         nfocus = 1
         write(6,'(a)') ' PB Info in pb_read(): nfocus has been reset to 1 with ipb=4/5/6'
      end if
      if ( istrng /= 0 ) then
         write(6,'(a)') ' PB Info in pb_read(): ipb=4/5/6 cannot handle the salt solution'
         call mexit(6, 1)
      end if
      if ( ligand ) then
         write(6,'(a)') "IIM/AUG/ANAIIM doesn't work with ligandmask."; call mexit(6,1)
      end if
   end if

   ! focusing options
   ! RL: I'm commenting out the nfocus > 2 check to be consistent with Amber13.

   if ( nfocus > 1 .and. fscale < 2 ) then
      write(6,'(a)') ' PB Info in pb_read(): nfocus reset to 1 due small fscale.'
      nfocus = 1
   end if

   ! nonpolar options

   if ( inp == 2 .and. pqropt /= 0 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): PQR input cannot be used with inp=2.'
      call mexit(6, 1)
   end if
   if ( inp /= npopt ) then
      npopt = inp
   end if
   if ( inp == 1 .and. use_sav == 1 ) then
      use_sav = 0
   end if
   if ( inp == 1 .and. use_rmin == 1 ) then
      use_rmin = 0
   end if
   if ( inp == 1 .and. nint(10000 * sprob) == nint( 10000 * 0.557d0 ) ) then
      write(6, '(a)') &
         ' PBSA Warning in pb_read(): sprob=0.557 is optimized for inp=2 and &
         & should not be used with inp=1. It has been reset to 1.4.'
      sprob = 1.4d0
   end if
   if ( inp == 1 .and. nint(10000 * cavity_surften) == nint(10000 * 0.0378d0) ) then
      write(6, '(a)') &
         ' PBSA Warning in pb_read(): cavity_surften=0.0378 is optimized for inp=2 &
         & and should not be used with inp=1. It has been reset to 0.005000.'
      cavity_surften =  0.005000
   end if
   if ( inp == 1 .and. nint(10000 * cavity_offset) == nint(-0.56920d0 * 10000) ) then
      write(6, '(a)') &
         ' PBSA Warning in pb_read(): cavity_offset=-0.5692 is optimized for inp=2 &
         & and should not be used with inp=1. It has been reset to 0.000.'
      cavity_offset =  0.000d0
   end if
   if ( inp == 1 .and. radiopt /= 0) then
      write(6,'(a)') ' PBSA Warning in pb_read(): radiopt should be set to 0 when inp = 1.'
      !radiopt = 0
   end if
   if ( inp == 1 ) then
      donpsa = .false.
   end if

   ! check force options

   if ( pqropt /= 0 .and. frcopt /= 0 ) then
      write(6, '(a)') &
         ' PBSA Warning in pb_read(): PQR files do not support force calculation'
      frcopt = 0
   end if
   if ( pqropt /= 0 .and. eneopt > 2 ) then
      write(6, '(a)') &
         ' PBSA Warning in pb_read(): PQR files do not support P3M energy calculation'
      eneopt = 1
   end if
   if ( eneopt == -1 ) then
      if ( dbfopt == 0 ) then
         eneopt = 1
      else if ( dbfopt == 1 .or. dbfopt == -1 ) then
         eneopt = 2
      else
         write(6,'(a)') ' PB Info in pb_read(): only dbfopt = 0 or 1 are supported'
         write(6,'(a)') ' PB Info in pb_read(): dbfopt is replaced by eneopt'
         call mexit(6, 1)
      end if
   end if
   if ( dbfopt /= -1 ) then
      write(6,'(a)') ' PB Info in pb_read(): dbfopt is ignored when eneopt is set'
   end if
   if ( eneopt < 1 .or. eneopt > 4 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): only eneopt= 1-4 are supported'
      call mexit(6, 1)
   end if
   if ( frcopt < 0 .or. frcopt > 5 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): only frcopt= 0-5 are supported'
      call mexit(6, 1)
   end if
   if ( eneopt == 1 .and. frcopt > 1 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): combination of eneopt and frcopt is unsupported'
      call mexit(6, 1)
   end if
   if ( eneopt == 2 .and. (frcopt == 1 .or. frcopt == 5) ) then
      write(6,'(a)') ' PB Bomb in pb_read(): combination of eneopt and frcopt is unsupported'
      call mexit(6, 1)
   end if
   if ( eneopt == 3 .and. frcopt /= 0 .and. frcopt /= 2 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): combination of eneopt and frcopt is unsupported'
      call mexit(6, 1)
   end if
   if ( eneopt == 4 .and. frcopt /= 0 .and. frcopt /= 5 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): combination of eneopt and frcopt is unsupported'
      call mexit(6, 1)
   end if
   if ( (frcopt >= 2 .and. frcopt <= 4) .and. bcopt /= 7 ) then
      bcopt = 7
      write(6,'(a)') ' PB Info in pb_read(): bcopt has been reset to 7 for frcopt= 2 - 4'
   end if
   if ( ( frcopt == 1 .and. eneopt == 1 ) .and. .not.(bcopt == 4 .or. bcopt == 5) ) then
      bcopt = 4
      write(6,'(a)') ' PB Info in pb_read(): bcopt has been reset to 4 for eneopt/frcopt= 1'
   end if
!   if ( frcopt > 0 .and. smoothopt == 0 ) then
!     smoothopt = 1
!     write(6,'(a)') ' PB Info in pb_read(): smoothopt has been reset to 1 for force compuation'
!   end if

   if ( ivcap /= 0 .and. (eneopt /= 1 .or. frcopt /= 0) ) then
!     write(6,'(a)') "cap water should only be run with eneopt=1 and frcopt=0"
!     call mexit(6,1)
   end if

   ! surface area calculation options

   if ( saopt /= 0 .and. abs(saopt) /= 1 .and. abs(saopt) /= 2 ) then
      saopt = 0
      write(6,'(a)') ' PB Info in pb_read(): saopt has been reset to 0 for surface area compuation'
   end if

   ! numerical solver options

   if ( solvopt == 7 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): solvopt=7 is no longer supported in 2017 and later release '
      call mexit(6, 1)
   else if ( npbopt == 0 .and. solvopt > 5 .and. (.not. solvopt == 8) ) then
      write(6,'(a)') ' PB Bomb in pb_read(): solvopt>5 cannot be used to solve linear PB equation'
      call mexit(6, 1)
   endif
   if ( solvopt == 2 .and. bcopt == 6 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): bcopt=6 cannot be used with solvopt=2'
      call mexit(6, 1)
   end if
   if ( npbopt == 1 .and. bcopt == 10 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): bcopt=10 can be used only with npbopt=0'
      call mexit(6, 1)
   end if
   if ( nfocus > 1 .and. bcopt == 10 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): bcopt=10 can be used only with nfocus=1'
      call mexit(6, 1)
   end if
   if ( solvopt < 1 .and. solvopt > 4 .and. bcopt == 10 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): bcopt=10 can be used only with solvopt=1,2,3,4'
      call mexit(6, 1)
   end if
   if ( npbopt == 1 .and. solvopt > 6 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): unsupported solvopt (>6) for e nonlinear PB equation'
      call mexit(6, 1)
   endif
   if ( npbopt == 1 .and. eneopt == 2 ) then
      write(6,'(a)') ' PB Info in pb_read(): pure charge view method uses '
      write(6,'(a)') ' FD greens function for nonlinear PB equation'
   end if
   if ( npbopt == 1 .and. frcopt > 0 .and. eneopt /= 4 ) then
      write(6,'(a)') ' PB Bomb in pb_read(): force computatoin is only supported for'
      write(6,'(a)') ' nonlinear PB equation with the new p3m method, eneopt = 4'
      call mexit(6, 1)
   end if

   ! cutoff options

   if ( cutfd > cutsa ) then
      cutsa = cutfd
      write(6,'(a)') ' PB Info in pb_read(): cutsa has been reset to be equal to cutfd'
   end if
   if ( cutnb /= 0 .and. cutfd > cutnb ) then
      cutnb = cutfd
      write(6,'(a)') ' PB Info in pb_read(): cutnb has been reset to be equal to cutfd'
   end if
   if ( max(cutnb,cutsa,cutfd) > cutres ) then
      write(6,'(a,3f12.4)') ' PB Bomb in pb_read(): cutnb/cutfd must be <= cutres', cutnb, cutfd, cutres
      call mexit(6, 1)
   end if
   if ( cutnb /= 0 .and. (bcopt == 2 .or. bcopt > 5 .and. bcopt < 10) ) then
      cutnb = 0
      write(6,'(a,f12.4,i6)') ' PB Info in pb_read(): cutnb has been reset to 0 with bcopt=2, 6-9', cutnb, bcopt
   end if
   if ( cutnb == 0 ) then
      if ( eneopt == 1 .and. (bcopt == 1 .or. bcopt == 4 .or. bcopt == 5 .or. bcopt == 10) ) then
         write(6,'(a,f12.4,2i6)') ' PB Bomb in pb_read(): cutnb=0 cannot be used with eneopt=1 ', cutnb, eneopt, bcopt
         call mexit(6, 1)
      end if
      if ( eneopt > 2 ) then
         write(6,'(a,f12.4,i6)') ' PB Bomb in pb_read(): cutnb=0 cannot be used with eneopt=3,4', cutnb, eneopt
         call mexit(6, 1)
      end if
   end if

   ! let's set up some initial variables

   if ( npbverb > 0 ) pbverbose = .true.

   cutfd = cutfd**2
   cutsa = cutsa**2
   cutnb = cutnb**2
   cutres = cutres**2

   dofd  = ndofd
   dosas  = ndosas

   epsin  = epsin*eps0
   epsout = epsout*eps0
   epsmem= epsmem*eps0
   istrng = fioni * istrng
   pbkappa  = SQRT( 2.0d0 * istrng / (epsout * pbkb * pbtemp) )

   ! set buffer zone between the fine FD grid boundary and the solute surface:

   if ( nbuffer == 0 ) then
      if ( istrng == 0.0d0 ) then
         nbuffer = int(2.0d0*dprob/space)+1
      else
         if ( dprob >= iprob ) then
            nbuffer = int(2.0d0*dprob/space)+1
         else
            nbuffer = int(2.0d0*iprob/space)+1
         end if
      end if
      if ( nfocus > 1 ) then
         if ( nbuffer >= fscale ) then
            nbuffer = 2*nbuffer+1
         else
            nbuffer = 2*fscale+1
         end if
      end if
   end if

   ! set flag to scale induced surface charges:

   if ( scalec == 1) scalerf = .true.

   ! set saved grid options

   savh(nfocus) = 0
   savbcopt(nfocus) = 0

   do l = 1, nfocus
      savbcopt(l) = bcopt
   end do
   savh(nfocus) = space
   do l = nfocus - 1, 1, -1
      savh(l) = savh(l+1)*fscale
   end do

   ! if grid parameters are read in ... set up saved arrays here
   ! it is assumed that no focus is used for this special application

   if ( rxm /= 0 .and. rym /= 0 .and. rzm /= 0 ) then

      ! if grid parameter is read in then do not call setgrd

      pbgrid = .false.

      if ( nfocus > 1 ) then
         write(6,'(a,4i6)') ' PB Bomb in pb_read(): focus is invoked while attempting to read in grid parameters',&
                             rxm, rym, rzm, nfocus
         call mexit(6, 1)
      end if
      savxm(nfocus) = rxm
      savym(nfocus) = rym
      savzm(nfocus) = rzm
      savxmym(nfocus) = rxm*rym
      savxmymzm(nfocus) = rxm*rym*rzm

   ! else, set these to be minimum values ...

   else

      if ( .not. pbgrid ) then
         write(6,'(a)') ' PB Bomb in pb_read(): both grid dimension and grid origins are required if read in'
         call mexit(6, 1)
      end if

      do l = 1, nfocus
         savxm(l) = 1
         savym(l) = 1
         savzm(l) = 1
         savxmym(l) = 1
         savxmymzm(l) = 1
      end do

   end if

   ! if grid origin is read in ... set up saved arrays here
   ! it is assumed that no focus is used for this special application

   if ( rgox /= 0.0d0 .and. rgoy /= 0.0d0 .and. rgoz /= 0.0d0 ) then

      ! if grid parameter is read in then do not call setgrd

      pbgrid = .false.

      if ( nfocus > 1 ) then
         write(6,'(a,3f8.3,i6)') ' PB Bomb in pb_read(): focus is invoked while attempting to read in grid parameters',&
         rgox, rgoy, rgoz, nfocus
         call mexit(6, 1)
      end if
      savgox(nfocus) = rgox
      savgoy(nfocus) = rgoy
      savgoz(nfocus) = rgoz

   else
      if ( .not. pbgrid ) then
         write(6,'(a)') ' PB Bomb in pb_read(): both grid dimension and grid origins are required if read in'
         call mexit(6, 1)
      end if
   end if

   ! set up ligand focusing option

   ligand = .false.
   if ( len_trim(ligandmask) > 0 ) ligand = .true.
   if ( xmax-xmin > 0 .and. &
        ymax-ymin > 0 .and. &
        zmax-zmin > 0 ) ligand = .true.

   if (ligand) then
      if ( frcopt /= 0 ) then
         write(6,'(a)') ' PB Bomb in pb_read(): Ligandmask is not compatible with frcopt /= 0.'
         call mexit(6,1)
      else if ( eneopt == 2 ) then
         write(6,'(a)') ' PB Bomb in pb_read(): Ligandmask is not compatible with eneopt == 2.'
         call mexit(6,1)
      else if ( nfocus < 2 ) then
         write(6,'(a)') ' PB Bomb in pb_read(): You cannot use ligandmask without focusing.'
         call mexit(6,1)
      end if
      !buffer = max(space*nbuffer/2,buffer)
      !nbuffer = max(nbuffer,nint(buffer/space*2))
   end if

   ! set phimap output options when requested
   ! it also turns of nonpolar solvent model calculations

   if ( phiout == 1 ) then
      outphi = .true.
      if ( radiopt /= 2 ) write(6,'(a,a,i6)') &
      ' PB Info in pb_read(): radiopt is no longer coerced by inp mode.',&
      ' Please make sure atomic radii are consistent with those in the visualization program.',&
      radiopt
      inp = 0
      donpsa = .false.
      npopt = 0
   end if

end subroutine pb_read

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBSA residue masking function
subroutine myresmask ( maskstr, masklength, resmask, nres )
   ! Mengjuei Hsieh

   implicit none
   character(*) maskstr
   integer      resmask(*)
   integer      masklength, nres

   ! Local

   integer      strpos
   integer      i, range_start, range_end
   logical      assigned
   logical      res_mode, num_mode, rng_mode

   i=0
   res_mode=.false.
   num_mode=.false.
   rng_mode=.false.
   assigned=.false.

   ! warning elminator

   range_start = -1
   range_end = -1

   do strpos = 1, masklength
      if ( maskstr(strpos:strpos)==' ' ) cycle
      if ( maskstr(strpos:strpos)==':' .and. res_mode ) then
         print *,"Error parsing residue mask: dupe commas "; call mexit(6,1)
      else if ( maskstr(strpos:strpos) == ':' ) then
         res_mode = .true.
         cycle
      else if ( .not. res_mode ) then
         print *,"Error parsing residue mask: lack commas "; call mexit(6,1)
      else if ( maskstr(strpos:strpos) == '@' ) then
         print *,"mjhsieh: unsupported mask syntax "; call mexit(6,1)
      else if ( maskstr(strpos:strpos) >= 'A' .and. &
                maskstr(strpos:strpos) <= 'z' ) then
         print *,"mjhsieh: unsupported mask syntax "; call mexit(6,1)
      else if ( maskstr(strpos:strpos) >= '0' .and. &
                maskstr(strpos:strpos) <= '9' ) then
         num_mode = .true.
         assigned = .false.
         i = i*10 + ichar(maskstr(strpos:strpos))-ichar('0')
         if ( i > nres ) then
            print *,"mjhsieh: mask out of range"; call mexit(6,1)
         end if
      else if ( maskstr(strpos:strpos) == '-' ) then
         if ( rng_mode ) then
            print *,"mjhsieh: mask syntax error"
         else if ( num_mode ) then
            rng_mode = .true.
            num_mode = .false.
            range_start = i; i = 0
         else
            print *,"mjhsieh: mask syntax error"
         end if
      else if ( maskstr(strpos:strpos) == ',' ) then
         if ( rng_mode ) then
            range_end = i; i = 0
            resmask(range_start:range_end)=1
            rng_mode = .false.
         else if ( num_mode ) then
            resmask(i:i) = 1; i = 0
            assigned=.true.
         else
            print *,"mjhsieh; mask syntax error"
         end if
      else
         print *,"mjhsieh: unsupported mask syntax "; call mexit(6,1)
      end if
   end do
   if ( assigned ) return

   if ( .not. rng_mode ) range_start=i
   range_end=i
   if ( range_start > range_end ) then
      print *,"reverse number order in mask, are you sure?"; call mexit(6,1)
   end if
   resmask(range_start:range_end)=1
   !print *, range_start,"-",range_end
end subroutine myresmask

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBSA echoin function
subroutine myechoin (ilun,iout)
   ! Meng-Juei Hsieh

   implicit none
   integer ilun, iout

   logical boolinside, boolhastag
   character(len=80) mylnbuf
   boolinside=.false.
   boolhastag=.false.

   do while (.true.)
      read(ilun,'(80a)',end=1949) mylnbuf
      if (boolinside) then
         if (mylnbuf(1:14)=='!! BEGIN Amber') then
            write(iout,"('Error: Nested interface tags')")
            call mexit(iout,1)
         elseif (mylnbuf(1:14)=='!! END   Amber') then
            boolinside=.false.
         else
            write(iout,'(80a)') mylnbuf(1:79)
         endif
      elseif (mylnbuf(1:14)=='!! BEGIN Amber') then
         write(iout,*)
         write(iout,"(' The Interface script used to generate the input file:')")
         write(iout,*)
         boolinside=.true.
         boolhastag=.true.
      elseif (mylnbuf(1:14)=='!! END   Amber') then
         write(iout,"('Error: unclosed interface tag')")
         call mexit(iout,1)
      endif
   enddo

   1949 rewind(ilun)
   if (boolhastag) write(iout,"(79('-'))")
   write(iout,*)
   write(iout,"(' Here is the input file:')")
   write(iout,*)

   boolinside=.false.
   do while (.true.)
      read(ilun,'(80a)',end=9562) mylnbuf
      if (boolinside) then
         if (mylnbuf(1:14)=='!! END   Amber') then
            boolinside=.false.
         endif
      elseif (mylnbuf(1:14)=='!! BEGIN Amber') then
         boolinside=.true.
      else
         write(iout,'(80a)') mylnbuf(1:79)
      endif
   enddo

   9562 rewind(ilun)

   return
end subroutine myechoin

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBSA name list function
subroutine mynmlsrc(srchkey,ilun,ifound)
   ! Meng-Juei Hsieh

   implicit none
   character(len=*) srchkey
   integer ilun

   character(len=80) realkey
   character(len=80) mylnbuf
   integer ifound, ilen, ipos

   ilen = len(srchkey)
   realkey(2:ilen+1)=srchkey(1:ilen)
   ifound = 0
   ipos = 0

   do while (ifound == 0)
      read(ilun,'(80a)',end=824,err=9528) mylnbuf
      realkey(1:1)='&'
      ipos=index(mylnbuf, realkey(1:ilen+1), .false.)
      if (ipos>0) then
         ifound = 1
      else
         realkey(1:1)='$'
         ipos=index(mylnbuf, realkey(1:ilen+1), .false.)
         if (ipos>0) ifound = 1
      endif
   enddo

   824 continue
   !write(6,'(a)') realkey(1:ilen+1)

   if (ifound == 1 ) then
      backspace(ilun)
   else if (ifound == 0) then
      rewind(ilun)
   else
      write(6,'(a)') "mynmlsrc Error: exception"
      call mexit(6,1)
   endif
   return

   9528 rewind(ilun)
   return
end subroutine mynmlsrc
