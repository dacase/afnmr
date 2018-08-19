#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBMD cleanup routine
subroutine pb_free

   use poisson_boltzmann
   use solvent_accessibility

   implicit none
#  include "../include/md.h"
#ifdef SANDER
#  include "../sander/box.h"
#else
#  include "box.h"
#endif

   integer alloc_err(32)

   alloc_err(1:32) = 0

   if ( allocated(icrd   ) ) deallocate(   icrd, stat = alloc_err(1 ) )
   if ( allocated(grdcrg ) ) deallocate( grdcrg, stat = alloc_err(2 ) )
   if ( allocated(qgrdcrg) ) deallocate(qgrdcrg, stat = alloc_err(3 ) )

   if ( allocated(acrg) ) deallocate( acrg, stat = alloc_err(4 ) )
   if ( allocated(gcrg) ) deallocate( gcrg, stat = alloc_err(5 ) )

   if ( allocated(nshrt) ) deallocate(nshrt, stat = alloc_err(6 ) )
   if ( allocated(nex  ) ) deallocate(  nex, STAT = alloc_err(7 ) )
   if ( allocated(iex  ) ) deallocate(  iex, STAT = alloc_err(8 ) )

   if ( allocated(acrd) ) deallocate( acrd, stat = alloc_err(9 ) )
   if ( allocated(gcrd) ) deallocate( gcrd, stat = alloc_err(10) )

   if ( allocated(mdsig  ) ) deallocate(  mdsig, stat = alloc_err(11) )
   if ( allocated(rmin   ) ) deallocate(   rmin, stat = alloc_err(12) )
   if ( allocated(radi   ) ) deallocate(   radi, stat = alloc_err(13) )
   if ( allocated(radip  ) ) deallocate(  radip, stat = alloc_err(14) )
   if ( allocated(radip2 ) ) deallocate( radip2, stat = alloc_err(15) )
   if ( allocated(radip3 ) ) deallocate( radip3, stat = alloc_err(16) )
   if ( allocated(nzratm ) ) deallocate( nzratm, stat = alloc_err(17) )
   if ( allocated(nmax   ) ) deallocate(   nmax, stat = alloc_err(18) )
   if ( allocated(nexp   ) ) deallocate(   nexp, stat = alloc_err(19) )
   if ( allocated(sumnmax) ) deallocate(sumnmax, stat = alloc_err(20) )
   if ( allocated(sumnexp) ) deallocate(sumnexp, stat = alloc_err(21) )
   if ( allocated(avnmax ) ) deallocate( avnmax, stat = alloc_err(22) )
   if ( allocated(avnexp ) ) deallocate( avnexp, stat = alloc_err(23) )
   if ( allocated(scrd   ) ) deallocate(   scrd, stat = alloc_err(24) )

   if ( allocated(iar1pb ) ) deallocate( iar1pb, stat = alloc_err(25) )
   if ( allocated(iprshrt) ) deallocate(iprshrt, stat = alloc_err(26) )
   if ( allocated(cn1pb  ) ) deallocate(  cn1pb, stat = alloc_err(28) )
   if ( allocated(cn2pb  ) ) deallocate(  cn2pb, stat = alloc_err(29) )
   if ( allocated(cn3pb  ) ) deallocate(  cn3pb, stat = alloc_err(30) )

   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30)+&
        alloc_err(31)+alloc_err(32) /= 0 ) then
      write(6,'(a,i6)') 'PB Bomb in pb_free(): Deallocation aborted #1', alloc_err(1:32)
      call mexit(6, 1)
   end if

   alloc_err(1:32) = 0

   if ( allocated(phi    ) ) deallocate(    phi, stat = alloc_err(1 ) )
   if ( allocated(chgrd  ) ) deallocate(  chgrd, stat = alloc_err(2 ) )
   if ( allocated(epsx   ) ) deallocate(   epsx, stat = alloc_err(3 ) )
   if ( allocated(epsy   ) ) deallocate(   epsy, stat = alloc_err(4 ) )
   if ( allocated(epsz   ) ) deallocate(   epsz, stat = alloc_err(5 ) )
   if ( allocated(saltgrd) ) deallocate(saltgrd, stat = alloc_err(6 ) )
   if ( allocated(ioncrg ) ) deallocate( ioncrg, stat = alloc_err(7 ) )

   if ( allocated(bv) ) deallocate( bv, stat = alloc_err(8 ) )

   if ( allocated(insas ) ) deallocate( insas, stat = alloc_err(9 ) )
   if ( allocated(atmsas) ) deallocate(atmsas, stat = alloc_err(10) )
   if ( allocated(lvlset) ) deallocate(lvlset, stat = alloc_err(11) )
   if ( allocated(mlvlset)) deallocate(mlvlset,stat = alloc_err(12) )
   if ( allocated(zv    ) ) deallocate(    zv, stat = alloc_err(13) )

   if ( allocated(cphi) ) deallocate( cphi, stat = alloc_err(14) )

   if ( allocated(fedgex ) ) deallocate( fedgex, stat = alloc_err(15) )
   if ( allocated(fedgey ) ) deallocate( fedgey, stat = alloc_err(16) )
   if ( allocated(fedgez ) ) deallocate( fedgez, stat = alloc_err(17) )
   if ( allocated(iepsav ) ) deallocate( iepsav, stat = alloc_err(18) )
   if ( allocated(iepsavx) ) deallocate(iepsavx, stat = alloc_err(19) )
   if ( allocated(iepsavy) ) deallocate(iepsavy, stat = alloc_err(20) )
   if ( allocated(iepsavz) ) deallocate(iepsavz, stat = alloc_err(21) )

   if ( allocated(xs) ) deallocate( xs, stat = alloc_err(22) )

   if ( allocated(outflag    ) ) deallocate(    outflag, stat = alloc_err(23) )
   if ( allocated(outflagorig) ) deallocate(outflagorig, stat = alloc_err(24) )
   if ( allocated(mapout     ) ) deallocate(     mapout, stat = alloc_err(25) )

   if ( allocated(liveflag) ) deallocate( liveflag, stat = alloc_err(26) )
   if ( allocated(realflag) ) deallocate( realflag, stat = alloc_err(27) )

   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30)+&
        alloc_err(31)+alloc_err(32) /= 0 ) then
      write(6,'(a,i6)') 'PB Bomb in pb_free(): Deallocation aborted #2', alloc_err(1:32)
      call mexit(6, 1)
   end if

end subroutine pb_free
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBMD initialization routine
subroutine pb_init(ifcap, natom, nres, ntypes, nbonh, nbona, ipres, iac, ico, &
                   numex, natex, ibh, jbh, iba, jba, ibel, lbres, igraph, &
                   isymbl, cg, rin)

   ! Module variables

   use parms, only : cn1, cn2
   use poisson_boltzmann
   use solvent_accessibility
   use dispersion_cavity

   implicit none

   ! Common variables

#  include "../include/md.h"
#  include "pb_md.h"

   ! Passed variables

   integer ifcap,natom, nres, ntypes, nbonh, nbona
   integer ipres(*), iac(*), ico(*), numex(*), natex(*), ibh(*), jbh(*), iba(*), jba(*), ibel(*)
   character (len=4) :: lbres(*), igraph(*), isymbl(*)
   _REAL_ cg(natom), rin(natom)

   ! Local variables

   integer ires, iatm, jatm, maxmax, ic, i, j, jp, idum
   integer alloc_err(64)
   _REAL_ maxnba_l
   _REAL_ rinchk, crgchk
   _REAL_ ucrgh(natom), ucrga(natom)
   character (len=4) :: residue, resid(natom)
   integer focusresidue(nres)

   ! begin of code

   alloc_err(1:64) = 0

   ! allocate topology informations

   allocate(   icrd( 3,  natom), stat = alloc_err(1 ) )
   allocate( grdcrg( 3,8*natom), stat = alloc_err(2 ) )
   allocate(qgrdcrg(   8*natom), stat = alloc_err(3 ) )
   allocate(   acrg(     natom), stat = alloc_err(4 ) )
   allocate(   gcrg( 8,  natom), stat = alloc_err(5 ) )
   allocate(  nshrt( 0:  natom), stat = alloc_err(6 ) )
   allocate(    nex(     natom), STAT = alloc_err(7 ) )
   allocate(    iex(64,  natom), STAT = alloc_err(8 ) )

   if ( ifcap == 0 .or. (ifcap >= 3 .and. ifcap <= 5) ) then

   allocate(   acrd( 3,  natom), stat = alloc_err(9 ) )
   allocate(   gcrd( 3,  natom), stat = alloc_err(10) )

   ! allocate sas informations

   allocate(  mdsig(  natom  ), stat = alloc_err(11) )
   allocate(   rmin(  natom  ), stat = alloc_err(11) )
   allocate(   radi(  natom  ), stat = alloc_err(12) )
   allocate(  radip(  natom  ), stat = alloc_err(13) )
   allocate( radip2(  natom  ), stat = alloc_err(14) )
   allocate( radip3(  natom  ), stat = alloc_err(15) )
   allocate( nzratm(  natom  ), stat = alloc_err(16) )
   allocate(   nmax(  natom  ), stat = alloc_err(17) )
   allocate(   nexp(  natom  ), stat = alloc_err(18) )
   allocate(sumnmax(  natom  ), stat = alloc_err(19) )
   allocate(sumnexp(  natom  ), stat = alloc_err(20) )
   allocate( avnmax(  natom  ), stat = alloc_err(21) )
   allocate( avnexp(  natom  ), stat = alloc_err(22) )
   allocate(   scrd(3,maxsph ), stat = alloc_err(23) )

   else

   allocate(   acrd( 3, 0:natom), stat = alloc_err(9 ) )
   allocate(   gcrd( 3, 0:natom), stat = alloc_err(10) )

   allocate(   radi(    0: 0   ), stat = alloc_err(12) )
   allocate( radip3(    1: 1   ), stat = alloc_err(15) )
   allocate( nzratm(    1: 1   ), stat = alloc_err(16) )

   end if

   ! allocate pb nblists

   maxnba_l = dble(natom) * ( sqrt( max(cutnb,cutsa,cutfd) ) )**3 / 3.0d0
   if ( natom >= 65536 ) then
      write(6,'(a)') "PBSA Warning: natom**2 exceeds integer limit (2147483647)."
      maxmax = 2147483647
   else
      maxmax = ceiling(dble(natom)/2*dble(natom))
   end if
   if ( maxnba_l > maxmax ) then
      maxnba = maxmax
   else
      maxnba = int(maxnba_l)
   end if

   allocate( iar1pb (4,0:natom), stat = alloc_err(24) )
   allocate( iprshrt(  maxnba ), stat = alloc_err(25) )
   allocate( cn1pb  (  maxnba ), stat = alloc_err(27) )
   allocate( cn2pb  (  maxnba ), stat = alloc_err(28) )
   allocate( cn3pb  (  maxnba ), stat = alloc_err(29) )

   ! allocate ibelly and outflag

   if ( ifcap == 3 .or. ifcap == 4 ) then
      allocate( ibelly(natom), stat = alloc_err(30) )
   end if
   allocate( outflag(  natom  ), stat = alloc_err(31) )
   allocate( outflagorig(  natom  ), stat = alloc_err(32) )
   allocate( mapout(  natom  ), stat = alloc_err(33) )

   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30)+&
        alloc_err(31)+alloc_err(32)+alloc_err(33) /= 0 ) then
      write(6,'(a,i6)') 'PB Bomb in pb_init(): Allocation aborted', alloc_err(1:31)
      call mexit(6, 1)
   end if

   if ( pbverbose ) then
      write(6,'(a)')
      write(6,'(a)') '======== Implicit Solvent Initialization ========'
      write(6,'(a)')
      if ( pqropt == 0 ) write(6,'(5x,a,2i9)') 'Max Nonbonded Pairs:', maxnba, maxmax
      write(6,'(a)')
   end if

   ! the following setups are mostly for standard prmtop/inpcrd files

   if ( pqropt == 0 ) then

   ! getting some topology info into pb setup ...

   do ires = 1, nres
      write(residue,'(a4)') lbres(ires)
      do iatm = ipres(ires), ipres(ires+1) - 1
         resid(iatm) = residue
      enddo
   enddo

   ! atomic and total charge in electron for printing only

   totcrgp = ZERO
   totcrgn = ZERO
   do iatm = 1, natom

      acrg(iatm) = cg(iatm)*INV_AMBER_ELECTROSTATIC
      if ( acrg(iatm) > ZERO) then
         totcrgp = totcrgp + acrg(iatm)
      else
         totcrgn = totcrgn + acrg(iatm)
      end if

      ! get info about belly atoms for GBSP

      if(ifcap == 3 .or. ifcap == 4) then
        ibelly(iatm) = ibel(iatm)
      end if

   end do
   totcrg = totcrgp + totcrgn

   ! for pure implicit solvent, set up radii arrays

   if ( ifcap == 0 .or. (ifcap >= 3 .and. ifcap <= 5) ) then

      ! set up group charges for cavity radii analysis

      ucrgh(1:natom) = acrg(1:natom)
      do idum = 1, nbonh
         iatm = ibh(idum)/3 + 1
         jatm = jbh(idum)/3 + 1
         if (isymbl(iatm)(1:1) == 'H' ) then
            ucrgh(jatm) = ucrgh(jatm) + acrg(iatm)
         else
            ucrgh(iatm) = ucrgh(iatm) + acrg(jatm)
         endif
      enddo

      ucrga(1:natom) = ucrgh(1:natom)
      do idum = 1, nbona
         iatm = iba(idum)/3 + 1
         jatm = jba(idum)/3 + 1
         ucrga(iatm) = ucrga(iatm) + ucrgh(jatm)
         ucrga(jatm) = ucrga(jatm) + ucrgh(iatm)
      enddo

      do idum = 1, natom
         if (isymbl(idum)(1:1) == 'H' ) then
            ucrga(idum) = 0.0d0
            ucrgh(idum) = 0.0d0
         endif
      enddo

      ! van der Waals sigma radii for nonpolar solvation

      do iatm = 1, natom
         ic = ico(ntypes*(iac(iatm)-1) + iac(iatm))
         if (cn2(ic) /= ZERO) then
            mdsig(iatm) = (cn1(ic)/cn2(ic))**(SIXTH)/2 ! this is sigma
            rmin(iatm) = mdsig(iatm)*(2.0d0**(SIXTH)) ! this is Rmin
         else
            mdsig(iatm) = ZERO
            rmin(iatm) = ZERO
         endif
      end do

      ! cavity radii for polar solvation if not passing in from driver

      noradius = .false.
      if ( radiopt == 0 ) then
         rinchk = ZERO
         do iatm = 1, natom
            rinchk = rinchk+rin(iatm)
         end do
         if (rinchk == ZERO) then
            write(6,'(a)') ' PB Bomb in pb_init(): Requested radi to be read in, but found none'
            if ( ipb > 0 ) then
               write(6,'(a)') ' PB Bomb in pb_init(): No PB calculation done'
               write(6,'(a)') ' '
               noradius = .true.
            end if
            if ( inp == 1 ) then
               write(6,'(a)') ' PB Bomb in pb_init(): No NP calculation done'
               write(6,'(a)') ' '
               inp = 0
            end if
         end if
         radi = rin ! for pb
         mdsig = rin ! for np

         ! radiscale and protscale are default to ONE.
         ! if scale ligand and protein differently, please specify the ligand residue ID
         ! if scale all atoms in the same way, please use protscale
         if ( radiscale /= ONE .or. protscale /= ONE ) then
            if ( radires == '' )  then
               write(6,'(a)') 'PB Warning in pb_init(): radiscale will be used as no ligand is selected'
               write(6,'(a)') 'PB Warning in pb_init(): Scale all atomic radii with protscale'
               do iatm = 1, natom
                  radi(iatm) = radi(iatm)*protscale
               end do
            else
               do iatm = 1, natom
                  if ( resid(iatm) == radires(1:4) ) then
                  radi(iatm) = radi(iatm)*radiscale
                  else
                  radi(iatm) = radi(iatm)*protscale
                  end if
               end do
            end if
         end if
      else if ( radiopt == 1 ) then
         call pb_aaradi( natom, nbonh, ibh, jbh, radi, acrg, ucrgh, ucrga, resid, igraph, isymbl, rin )
         radi = radi*protscale
      else if ( radiopt == 2 ) then
         call phi_aaradi( natom, isymbl, radi )
      else
         write(6,'(a,i6)') 'PB Bomb in pb_init(): Unknown radi assigment option', radiopt
         if ( ipb > 0 ) then
            write(6,'(a)') 'PB Bomb in pb_init(): No PB calculation done'
            write(6,'(a)') ' '
            noradius = .true.
         end if
         if ( inp == 1 ) then
            write(6,'(a)') 'PB Bomb in pb_init(): No NP calculation done'
            write(6,'(a)') ' '
            inp = 0
         end if
      end if

      ! check if atomic charges are correctly input from prmtop file

      nocharge = .false.
      crgchk = ZERO
      do iatm = 1, natom
         crgchk = crgchk + abs(acrg(iatm))
      end do
      if ( crgchk < 1.0D-6 ) then
         write(6,'(a)') ' PB Info in pb_init(): Requested charges to be read in, but found none'
         write(6,'(a)') ' PB Info in pb_init(): No PB calculation done'
         write(6,'(a)') ' '
         nocharge = .true.
      end if

      ! let's do a summary of key parameters for pb before proceeding

      if ( pbverbose ) then
         write(6,'(a,i6)') ' no. of atoms processed in PB initialization:', natom
         write(6, '(a5,3a6,5a10)') 'NUM','RESI','NAME','TYPE',&
         'CHARGE', 'ATM CRG/H', 'GRP CRG', 'PB RADI', 'NP RADI'
         do iatm = 1, natom
            if ( use_rmin == 0 ) write(6, '(i5,3a6,5f10.6)') iatm,resid(iatm),igraph(iatm),isymbl(iatm),&
            real(acrg(iatm)),real(ucrgh(iatm)), real(ucrga(iatm)),&
            real(radi(iatm)), real(mdsig(iatm))
            if ( use_rmin == 1 ) write(6, '(i5,3a6,5f10.6)') iatm,resid(iatm),igraph(iatm),isymbl(iatm),&
            real(acrg(iatm)),real(ucrgh(iatm)), real(ucrga(iatm)),&
            real(radi(iatm)), real(rmin(iatm))
         end do
         write(6,*)
         write(6,'(a,3f14.4)') ' total system charges (+/-) for PB', totcrg, totcrgp, totcrgn
         write(6,'(a,f14.4,a,f14.4)') ' cavity_surften =', cavity_surften, ' cavity_offset =', cavity_offset
         write(6,*)
      end if

      ! initialization for sas surface

      if ( srsas ) then

         call sa_sphere(maxsph, scrd)
         if ( pbverbose ) write(6,'(2x,a,i6)') ' SAS Surface: surface dots generated: ', maxsph

         ! nmax and nexp accumulators

         sumnmax(1:natom) = ZERO
         sumnexp(1:natom) = ZERO

      ! initialization for for vdw surface

      else
         if ( pbverbose ) write(6,'(a)') ' VDW Surface: setting up working radii'
         radip3(1:natom) = radi(1:natom)
      end if

   else
      acrd(1:3,0:0) = ZERO; gcrd(1:3,0:0) = ZERO
   end if ! if ( ifcap == 0 .or. (ifcap >= 3 .and. ifcap <= 5) ) then

   ! assinging atom-based pointers to exclusion list

   nshrt(0) = 0
   do i = 1, natom
      nshrt(i) = nshrt(i-1) + numex(i)
   end do
   nex = 0
   do i = 1, natom-1
      do jp = nshrt(i-1) + 1, nshrt(i)
         j = natex(jp)
         if (j == 0)  cycle
         nex(i) = nex(i) + 1
         iex(nex(i),i) = j
         nex(j) = nex(j) + 1
         iex(nex(j),j) = i
      end do
   end do

   ! only the following setups are needed/can be done for the pqr file

   else

      ! apparently there is exclusion list

      nshrt = 0

      ! initialization for sas surface

      if ( srsas ) then

         call sa_sphere(maxsph, scrd)
         if ( pbverbose ) write(6,'(2x,a,i6)') ' SAS Surface: surface dots generated: ', maxsph

         ! nmax and nexp accumulators

         sumnmax(1:natom) = ZERO
         sumnexp(1:natom) = ZERO

      end if

   end if ! if ( pqropt == 0 ) then

   ! set up green 3-d array for pb_fdcoulomb

   call pb_green

   ! set up variables if using ligand

   if ( ligand ) then

      allocate( liveflag(natom), stat = alloc_err(1) )
      allocate( realflag(natom), stat = alloc_err(2) )
      if ( SUM(alloc_err(1:2)) /= 0 ) then
         write(6,'(a,i6)') ' PB Bomb in pb_init(): Allocation aborted', alloc_err(1:2)
         call mexit(6,1)
      end if

      ! Mengjuei -
      ! Basically liveflag HERE is a temporary switch storing the list of
      ! atoms labelled ligandmask. However if the geormetries of ligand
      ! box were specified, we don't need liveflag that early.
      ! Later liveflag will become the list of atom inside the block of
      ! rectangular bounding box.

      if ( xmax - xmin > ZERO .and. &
         ymax - ymin > ZERO .and. &
         zmax - zmin > ZERO       ) then
         ! If the geometries of bounding box are correctly specified,
         ! skip the ligandmask
         continue
      else if ( ligand .and. len_trim(ligandmask) > 0 ) then
         focusresidue = 0
         liveflag = 0
         call myresmask ( ligandmask, LEN_TRIM(ligandmask), focusresidue, nres )
         j = 0
         do i = 1, nres
            if ( focusresidue(i) == 0 ) cycle
            j = j + 1
            liveflag(ipres(i):ipres(i+1)-1)=1
         enddo
         write(6,'(a,a,a,i5)') ' Focusing Mask ', &
            ligandmask(1:LEN_TRIM(ligandmask)), ' matches residue', j
         xmax = ZERO; ymax = ZERO; zmax = ZERO
         xmin = ZERO; ymin = ZERO; zmin = ZERO
      else
         write(6,'(a)') 'wrong ligand box'
         call mexit(6,1)
      endif

   endif

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Setup table of precomputed FD Green's function.
subroutine pb_green

   implicit none

   ! Common variables

   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green

   ! Local variables

   integer i, j, k
   _REAL_ green_data

   do i = 0, 40
      do j = 0, i
         do k = 0, j
            green_data = green(i, j, k)
            green(i, k, j) = green_data
            green(j, i, k) = green_data
            green(k, i, j) = green_data
            green(j, k, i) = green_data
            green(k, j, i) = green_data
         end do
      end do
   end do

end subroutine pb_green

end subroutine pb_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Data for precomputed FD Green's function
block data blkgreen

   implicit none

   ! Common variables

   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green

   data green( 0, 0, 0) /   3.1759115351965534D+00  /
   data green( 1, 0, 0) /   1.0815164328033549D+00  /
   data green( 1, 1, 0) /   6.9355601002926537D-01  /
   data green( 1, 1, 1) /   5.4762175169703808D-01  /
   data green( 2, 0, 0) /   5.3896302150652231D-01  /
   data green( 2, 1, 0) /   4.5152984558740150D-01  /
   data green( 2, 1, 1) /   4.0168749336481041D-01  /
   data green( 2, 2, 0) /   3.5255169608092701D-01  /
   data green( 2, 2, 1) /   3.2872036419869677D-01  /
   data green( 2, 2, 2) /   2.8464545999514379D-01  /
   data green( 3, 0, 0) /   3.4614231388617117D-01  /
   data green( 3, 1, 0) /   3.2073335917807405D-01  /
   data green( 3, 1, 1) /   3.0200278891962756D-01  /
   data green( 3, 2, 0) /   2.7740487845668244D-01  /
   data green( 3, 2, 1) /   2.6587502119324352D-01  /
   data green( 3, 2, 2) /   2.4057055579159040D-01  /
   data green( 3, 3, 0) /   2.3517700631231578D-01  /
   data green( 3, 3, 1) /   2.2833414358944390D-01  /
   data green( 3, 3, 2) /   2.1176997261789529D-01  /
   data green( 3, 3, 3) /   1.9125147492997435D-01  /
   data green( 4, 0, 0) /   2.5495742509820746D-01  /
   data green( 4, 1, 0) /   2.4531753929893496D-01  /
   data green( 4, 1, 1) /   2.3711247941031846D-01  /
   data green( 4, 2, 0) /   2.2421716678228923D-01  /
   data green( 4, 2, 1) /   2.1821739620342098D-01  /
   data green( 4, 2, 2) /   2.0348788713212107D-01  /
   data green( 4, 3, 0) /   1.9979199689082094D-01  /
   data green( 4, 3, 1) /   1.9565392010998225D-01  /
   data green( 4, 3, 2) /   1.8494655280238623D-01  /
   data green( 4, 3, 3) /   1.7073297724205330D-01  /
   data green( 4, 4, 0) /   1.7649822175475413D-01  /
   data green( 4, 4, 1) /   1.7368106393811655D-01  /
   data green( 4, 4, 2) /   1.6610647283022761D-01  /
   data green( 4, 4, 3) /   1.5557916351723405D-01  /
   data green( 4, 4, 4) /   1.4383287095970104D-01  /
   data green( 5, 0, 0) /   2.0233207950733323D-01  /
   data green( 5, 1, 0) /   1.9777232591440136D-01  /
   data green( 5, 1, 1) /   1.9360221653757156D-01  /
   data green( 5, 2, 0) /   1.8635379364045454D-01  /
   data green( 5, 2, 1) /   1.8295790259257091D-01  /
   data green( 5, 2, 2) /   1.7402886898952086D-01  /
   data green( 5, 3, 0) /   1.7155174627560230D-01  /
   data green( 5, 3, 1) /   1.6895236723570453D-01  /
   data green( 5, 3, 2) /   1.6192808688203741D-01  /
   data green( 5, 3, 3) /   1.5209495588310476D-01  /
   data green( 5, 4, 0) /   1.5602160443532481D-01  /
   data green( 5, 4, 1) /   1.5408692441187632D-01  /
   data green( 5, 4, 2) /   1.4874275196062139D-01  /
   data green( 5, 4, 3) /   1.4103484141468445D-01  /
   data green( 5, 4, 4) /   1.3208657840216792D-01  /
   data green( 5, 5, 0) /   1.4126545371003910D-01  /
   data green( 5, 5, 1) /   1.3983695733689094D-01  /
   data green( 5, 5, 2) /   1.3582392088765388D-01  /
   data green( 5, 5, 3) /   1.2988672447183353D-01  /
   data green( 5, 5, 4) /   1.2279310690233598D-01  /
   data green( 5, 5, 5) /   1.1521211843231260D-01  /
   data green( 6, 0, 0) /   1.6794574828818579D-01  /
   data green( 6, 1, 0) /   1.6542610996454146D-01  /
   data green( 6, 1, 1) /   1.6304036280116682D-01  /
   data green( 6, 2, 0) /   1.5866571768529189D-01  /
   data green( 6, 2, 1) /   1.5659277294875293D-01  /
   data green( 6, 2, 2) /   1.5091334785578717D-01  /
   data green( 6, 3, 0) /   1.4923834821560472D-01  /
   data green( 6, 3, 1) /   1.4753562314215773D-01  /
   data green( 6, 3, 2) /   1.4280302442088691D-01  /
   data green( 6, 3, 3) /   1.3591090146313159D-01  /
   data green( 6, 4, 0) /   1.3864035604780006D-01  /
   data green( 6, 4, 1) /   1.3728680156459946D-01  /
   data green( 6, 4, 2) /   1.3347626533724810D-01  /
   data green( 6, 4, 3) /   1.2781887425314498D-01  /
   data green( 6, 4, 4) /   1.2103070281926526D-01  /
   data green( 6, 5, 0) /   1.2793779935790175D-01  /
   data green( 6, 5, 1) /   1.2687926029994986D-01  /
   data green( 6, 5, 2) /   1.2386716979797802D-01  /
   data green( 6, 5, 3) /   1.1931681810582123D-01  /
   data green( 6, 5, 4) /   1.1374332085276670D-01  /
   data green( 6, 5, 5) /   1.0763112996228928D-01  /
   data green( 6, 6, 0) /   1.1775701788045112D-01  /
   data green( 6, 6, 1) /   1.1693365376498031D-01  /
   data green( 6, 6, 2) /   1.1457015985192663D-01  /
   data green( 6, 6, 3) /   1.1094579267779280D-01  /
   data green( 6, 6, 4) /   1.0642359263985041D-01  /
   data green( 6, 6, 5) /   1.0136496444590526D-01  /
   data green( 6, 6, 6) /   9.6075969565313266D-02  /
   data green( 7, 0, 0) /   1.4363797036361628D-01  /
   data green( 7, 1, 0) /   1.4209214229703621D-01  /
   data green( 7, 1, 1) /   1.4060219444284022D-01  /
   data green( 7, 2, 0) /   1.3779050839364396D-01  /
   data green( 7, 2, 1) /   1.3644368361554321D-01  /
   data green( 7, 2, 2) /   1.3265962340592247D-01  /
   data green( 7, 3, 0) /   1.3150102300061872D-01  /
   data green( 7, 3, 1) /   1.3034042446739721D-01  /
   data green( 7, 3, 2) /   1.2705392184495914D-01  /
   data green( 7, 3, 3) /   1.2212665554762085D-01  /
   data green( 7, 4, 0) /   1.2407078114877013D-01  /
   data green( 7, 4, 1) /   1.2310238014856431D-01  /
   data green( 7, 4, 2) /   1.2033897002625792D-01  /
   data green( 7, 4, 3) /   1.1614371637871881D-01  /
   data green( 7, 4, 4) /   1.1097324830160045D-01  /
   data green( 7, 5, 0) /   1.1620544790922047D-01  /
   data green( 7, 5, 1) /   1.1541317997734810D-01  /
   data green( 7, 5, 2) /   1.1313659430526810D-01  /
   data green( 7, 5, 3) /   1.0963902658141100D-01  /
   data green( 7, 5, 4) /   1.0526457468703798D-01  /
   data green( 7, 5, 5) /   1.0035809074407886D-01  /
   data green( 7, 6, 0) /   1.0839960051847128D-01  /
   data green( 7, 6, 1) /   1.0775811212880239D-01  /
   data green( 7, 6, 2) /   1.0590358653641528D-01  /
   data green( 7, 6, 3) /   1.0302368368166874D-01  /
   data green( 7, 6, 4) /   9.9372078504935402D-02  /
   data green( 7, 6, 5) /   9.5213982272844708D-02  /
   data green( 7, 6, 6) /   9.0786974684721203D-02  /
   data green( 7, 7, 0) /   1.0095457742024190D-01  /
   data green( 7, 7, 1) /   1.0043709608501153D-01  /
   data green( 7, 7, 2) /   9.8933308963866914D-02  /
   data green( 7, 7, 3) /   9.6576010518451216D-02  /
   data green( 7, 7, 4) /   9.3550022573802100D-02  /
   data green( 7, 7, 5) /   9.0055966159822004D-02  /
   data green( 7, 7, 6) /   8.6281554868525373D-02  /
   data green( 7, 7, 7) /   8.2384815605608827D-02  /
   data green( 8, 0, 0) /   1.2551350470536657D-01  /
   data green( 8, 1, 0) /   1.2449387617473523D-01  /
   data green( 8, 1, 1) /   1.2350115203071543D-01  /
   data green( 8, 2, 0) /   1.2159680014783039D-01  /
   data green( 8, 2, 1) /   1.2067657803470230D-01  /
   data green( 8, 2, 2) /   1.1804918165874267D-01  /
   data green( 8, 3, 0) /   1.1722565131089910D-01  /
   data green( 8, 3, 1) /   1.1640591505254005D-01  /
   data green( 8, 3, 2) /   1.1405483320166934D-01  /
   data green( 8, 3, 3) /   1.1045375537523730D-01  /
   data green( 8, 4, 0) /   1.1187309963785286D-01  /
   data green( 8, 4, 1) /   1.1116412370701316D-01  /
   data green( 8, 4, 2) /   1.0912094214278861D-01  /
   data green( 8, 4, 3) /   1.0596552356227758D-01  /
   data green( 8, 4, 4) /   1.0199220485882379D-01  /
   data green( 8, 5, 0) /   1.0599814647548324D-01  /
   data green( 8, 5, 1) /   1.0539728507228295D-01  /
   data green( 8, 5, 2) /   1.0365763291219818D-01  /
   data green( 8, 5, 3) /   1.0094877232995098D-01  /
   data green( 8, 5, 4) /   9.7501683137435430D-02  /
   data green( 8, 5, 5) /   9.3560300582418598D-02  /
   data green( 8, 6, 0) /   9.9964335643309379D-02  /
   data green( 8, 6, 1) /   9.9461555890587694D-02  /
   data green( 8, 6, 2) /   9.7999660286958612D-02  /
   data green( 8, 6, 3) /   9.5705607271006451D-02  /
   data green( 8, 6, 4) /   9.2756615174408363D-02  /
   data green( 8, 6, 5) /   8.9345819097605389D-02  /
   data green( 8, 6, 6) /   8.5654804260273956D-02  /
   data green( 8, 7, 0) /   9.4027035657242714D-02  /
   data green( 8, 7, 1) /   9.3609232934177625D-02  /
   data green( 8, 7, 2) /   9.2389787053454075D-02  /
   data green( 8, 7, 3) /   9.0462682104850464D-02  /
   data green( 8, 7, 4) /   8.7962000877333968D-02  /
   data green( 8, 7, 5) /   8.5038127485457304D-02  /
   data green( 8, 7, 6) /   8.1837299038139272D-02  /
   data green( 8, 7, 7) /   7.8488076342692337D-02  /
   data green( 8, 8, 0) /   8.8347099219041922D-02  /
   data green( 8, 8, 1) /   8.8000815664067450D-02  /
   data green( 8, 8, 2) /   8.6986732384313065D-02  /
   data green( 8, 8, 3) /   8.5374106193392724D-02  /
   data green( 8, 8, 4) /   8.3263644590470240D-02  /
   data green( 8, 8, 5) /   8.0771236522732756D-02  /
   data green( 8, 8, 6) /   7.8013077927048174D-02  /
   data green( 8, 8, 7) /   7.5094908803896285D-02  /
   data green( 8, 8, 8) /   7.2105988413551828D-02  /
   data green( 9, 0, 0) /   1.1146755316964206D-01  /
   data green( 9, 1, 0) /   1.1075850583674679D-01  /
   data green( 9, 1, 1) /   1.1006380932257710D-01  /
   data green( 9, 2, 0) /   1.0871760893829935D-01  /
   data green( 9, 2, 1) /   1.0806273570284175D-01  /
   data green( 9, 2, 2) /   1.0617264407379035D-01  /
   data green( 9, 3, 0) /   1.0557115497401200D-01  /
   data green( 9, 3, 1) /   1.0497387959355930D-01  /
   data green( 9, 3, 2) /   1.0324528313574793D-01  /
   data green( 9, 3, 3) /   1.0055516317590933D-01  /
   data green( 9, 4, 0) /   1.0161577147793842D-01  /
   data green( 9, 4, 1) /   1.0108512018804984D-01  /
   data green( 9, 4, 2) /   9.9544569447315392D-02  /
   data green( 9, 4, 3) /   9.7133750288145901D-02  /
   data green( 9, 4, 4) /   9.4045567451916071D-02  /
   data green( 9, 5, 0) /   9.7151425517950751D-02  /
   data green( 9, 5, 1) /   9.6689071471067201D-02  /
   data green( 9, 5, 2) /   9.5342543335939861D-02  /
   data green( 9, 5, 3) /   9.3223160515376977D-02  /
   data green( 9, 5, 4) /   9.0487631191973089D-02  /
   data green( 9, 5, 5) /   8.7308708280350630D-02  /
   data green( 9, 6, 0) /   9.2438119427483756D-02  /
   data green( 9, 6, 1) /   9.2040709277995320D-02  /
   data green( 9, 6, 2) /   9.0879792058089631D-02  /
   data green( 9, 6, 3) /   8.9042230048201235D-02  /
   data green( 9, 6, 4) /   8.6652502158133263D-02  /
   data green( 9, 6, 5) /   8.3851084810229487D-02  /
   data green( 9, 6, 6) /   8.0775614605432935D-02  /
   data green( 9, 7, 0) /   8.7677735792507414D-02  /
   data green( 9, 7, 1) /   8.7339107254702550D-02  /
   data green( 9, 7, 2) /   8.6347105646557704D-02  /
   data green( 9, 7, 3) /   8.4768580715464142D-02  /
   data green( 9, 7, 4) /   8.2700913335015530D-02  /
   data green( 9, 7, 5) /   8.0256443217110243D-02  /
   data green( 9, 7, 6) /   7.7548153344838611D-02  /
   data green( 9, 7, 7) /   7.4679226766473927D-02  /
   data green( 9, 8, 0) /   8.3013446335815463D-02  /
   data green( 9, 8, 1) /   8.2726298256347072D-02  /
   data green( 9, 8, 2) /   8.1882949170755004D-02  /
   data green( 9, 8, 3) /   8.0534447987936020D-02  /
   data green( 9, 8, 4) /   7.8756261536013908D-02  /
   data green( 9, 8, 5) /   7.6637220823981778D-02  /
   data green( 9, 8, 6) /   7.4268862079690537D-02  /
   data green( 9, 8, 7) /   7.1737116898696524D-02  /
   data green( 9, 8, 8) /   6.9117068023207343D-02  /
   data green( 9, 9, 0) /   7.8538164723848558D-02  /
   data green( 9, 9, 1) /   7.8295129934218888D-02  /
   data green( 9, 9, 2) /   7.7579704572913216D-02  /
   data green( 9, 9, 3) /   7.6430747022821072D-02  /
   data green( 9, 9, 4) /   7.4906459007106949D-02  /
   data green( 9, 9, 5) /   7.3076604086992694D-02  /
   data green( 9, 9, 6) /   7.1014681444458105D-02  /
   data green( 9, 9, 7) /   6.8791462130493086D-02  /
   data green( 9, 9, 8) /   6.6470545288249888D-02  /
   data green( 9, 9, 9) /   6.4105946861472710D-02  /
   data green(10, 0, 0) /   1.0025779096549864D-01  /
   data green(10, 1, 0) /   9.9744378092650021D-02  /
   data green(10, 1, 1) /   9.9239220825569852D-02  /
   data green(10, 2, 0) /   9.8253721265523120D-02  /
   data green(10, 2, 1) /   9.7771894254122110D-02  /
   data green(10, 2, 2) /   9.6370645106820266D-02  /
   data green(10, 3, 0) /   9.5920138929816032D-02  /
   data green(10, 3, 1) /   9.5473068508164130D-02  /
   data green(10, 3, 2) /   9.4170609322243887D-02  /
   data green(10, 3, 3) /   9.2119156832430674D-02  /
   data green(10, 4, 0) /   9.2928708361714976D-02  /
   data green(10, 4, 1) /   9.2523305431405689D-02  /
   data green(10, 4, 2) /   9.1339777593220081D-02  /
   data green(10, 4, 3) /   8.9468517576079515D-02  /
   data green(10, 4, 4) /   8.7038436892435012D-02  /
   data green(10, 5, 0) /   8.9478372784664695D-02  /
   data green(10, 5, 1) /   8.9117345434184439D-02  /
   data green(10, 5, 2) /   8.8061033611591649D-02  /
   data green(10, 5, 3) /   8.6384035898050521D-02  /
   data green(10, 5, 4) /   8.4194165608626151D-02  /
   data green(10, 5, 5) /   8.1614517095279920D-02  /
   data green(10, 6, 0) /   8.5753801055144255D-02  /
   data green(10, 6, 1) /   8.5436609566040572D-02  /
   data green(10, 6, 2) /   8.4506503752885209D-02  /
   data green(10, 6, 3) /   8.3023737571137068D-02  /
   data green(10, 6, 4) /   8.1076538388971803D-02  /
   data green(10, 6, 5) /   7.8767421502744628D-02  /
   data green(10, 6, 6) /   7.6200407062187336D-02  /
   data green(10, 7, 0) /   8.1909598825097615D-02  /
   data green(10, 7, 1) /   8.1633561620629957D-02  /
   data green(10, 7, 2) /   8.0822417626880860D-02  /
   data green(10, 7, 3) /   7.9524105170223638D-02  /
   data green(10, 7, 4) /   7.7809691506037573D-02  /
   data green(10, 7, 5) /   7.5763159503138611D-02  /
   data green(10, 7, 6) /   7.3471474362184611D-02  /
   data green(10, 7, 7) /   7.1016743769080862D-02  /
   data green(10, 8, 0) /   7.8065081766800548D-02  /
   data green(10, 8, 1) /   7.7826341178523273D-02  /
   data green(10, 8, 2) /   7.7123406176462919D-02  /
   data green(10, 8, 3) /   7.5994043289169269D-02  /
   data green(10, 8, 4) /   7.4494883471572987D-02  /
   data green(10, 8, 5) /   7.2693917501350472D-02  /
   data green(10, 8, 6) /   7.0662922039119988D-02  /
   data green(10, 8, 7) /   6.8471173588417727D-02  /
   data green(10, 8, 8) /   6.6181095351799321D-02  /
   data green(10, 9, 0) /   7.4305917901511226D-02  /
   data green(10, 9, 1) /   7.4100156897928635D-02  /
   data green(10, 9, 2) /   7.3493226069464554D-02  /
   data green(10, 9, 3) /   7.2514711290517064D-02  /
   data green(10, 9, 4) /   7.1209439930399937D-02  /
   data green(10, 9, 5) /   6.9632021211213715D-02  /
   data green(10, 9, 6) /   6.7841149144940854D-02  /
   data green(10, 9, 7) /   6.5894656126428613D-02  /
   data green(10, 9, 8) /   6.3845863345559389D-02  /
   data green(10, 9, 9) /   6.1741348434695503D-02  /
   data green(10,10, 0) /   7.0689203118758787D-02  /
   data green(10,10, 1) /   7.0512114690102490D-02  /
   data green(10,10, 2) /   6.9988900943110416D-02  /
   data green(10,10, 3) /   6.9142684871975224D-02  /
   data green(10,10, 4) /   6.8008830793140365D-02  /
   data green(10,10, 5) /   6.6630994087524029D-02  /
   data green(10,10, 6) /   6.5056874411942120D-02  /
   data green(10,10, 7) /   6.3334377560641308D-02  /
   data green(10,10, 8) /   6.1508624198599454D-02  /
   data green(10,10, 9) /   5.9619959471520483D-02  /
   data green(10,10,10) /   5.7702901172333260D-02  /
   data green(11, 0, 0) /   9.1101680252749645D-02  /
   data green(11, 1, 0) /   9.0717808836991612D-02  /
   data green(11, 1, 1) /   9.0338970937297508D-02  /
   data green(11, 2, 0) /   8.9596413124129304D-02  /
   data green(11, 2, 1) /   8.9231974115813101D-02  /
   data green(11, 2, 2) /   8.8166219414399408D-02  /
   data green(11, 3, 0) /   8.7821111961317838D-02  /
   data green(11, 3, 1) /   8.7478583517837646D-02  /
   data green(11, 3, 2) /   8.6475724757080180D-02  /
   data green(11, 3, 3) /   8.4881524022027602D-02  /
   data green(11, 4, 0) /   8.5511356115059625D-02  /
   data green(11, 4, 1) /   8.5195812503100266D-02  /
   data green(11, 4, 2) /   8.4270630170684341D-02  /
   data green(11, 4, 3) /   8.2795947952195012D-02  /
   data green(11, 4, 4) /   8.0859687533282321D-02  /
   data green(11, 5, 0) /   8.2801610904808989D-02  /
   data green(11, 5, 1) /   8.2515679740336731D-02  /
   data green(11, 5, 2) /   8.1675995655269729D-02  /
   data green(11, 5, 3) /   8.0333600505491515D-02  /
   data green(11, 5, 4) /   7.8563834185046394D-02  /
   data green(11, 5, 5) /   7.6455220068587279D-02  /
   data green(11, 6, 0) /   7.9823496161538401D-02  /
   data green(11, 6, 1) /   7.9567736255404239D-02  /
   data green(11, 6, 2) /   7.8815432083571405D-02  /
   data green(11, 6, 3) /   7.7609012168489519D-02  /
   data green(11, 6, 4) /   7.6011711987152022D-02  /
   data green(11, 6, 5) /   7.4098822156660388D-02  /
   data green(11, 6, 6) /   7.1949036037832662D-02  /
   data green(11, 7, 0) /   7.6693851094873325D-02  /
   data green(11, 7, 1) /   7.6467295272534735D-02  /
   data green(11, 7, 2) /   7.5799823394525637D-02  /
   data green(11, 7, 3) /   7.4726160312652792D-02  /
   data green(11, 7, 4) /   7.3298549167302746D-02  /
   data green(11, 7, 5) /   7.1580008929404040D-02  /
   data green(11, 7, 6) /   6.9637460454742159D-02  /
   data green(11, 7, 7) /   6.7535939946806678D-02  /
   data green(11, 8, 0) /   7.3508845181332341D-02  /
   data green(11, 8, 1) /   7.3309542352970311D-02  /
   data green(11, 8, 2) /   7.2721459723984430D-02  /
   data green(11, 8, 3) /   7.1772705638302653D-02  /
   data green(11, 8, 4) /   7.0505947066466707D-02  /
   data green(11, 8, 5) /   6.8973297959075516D-02  /
   data green(11, 8, 6) /   6.7230955558135727D-02  /
   data green(11, 8, 7) /   6.5334507345380974D-02  /
   data green(11, 8, 8) /   6.3335430219634198D-02  /
   data green(11, 9, 0) /   7.0342744003802124D-02  /
   data green(11, 9, 1) /   7.0168211613751266D-02  /
   data green(11, 9, 2) /   6.9652476535854796D-02  /
   data green(11, 9, 3) /   6.8818126559272094D-02  /
   data green(11, 9, 4) /   6.7699733808848361D-02  /
   data green(11, 9, 5) /   6.6340022516074335D-02  /
   data green(11, 9, 6) /   6.4785739636482612D-02  /
   data green(11, 9, 7) /   6.3083910988519204D-02  /
   data green(11, 9, 8) /   6.1278910673583327D-02  /
   data green(11, 9, 9) /   5.9410498112540545D-02  /
   data green(11,10, 0) /   6.7249576764662714D-02  /
   data green(11,10, 1) /   6.7097135141444150D-02  /
   data green(11,10, 2) /   6.6646076978827865D-02  /
   data green(11,10, 3) /   6.5914477457283197D-02  /
   data green(11,10, 4) /   6.4930212969271553D-02  /
   data green(11,10, 5) /   6.3728108448816989D-02  /
   data green(11,10, 6) /   6.2346788266802802D-02  /
   data green(11,10, 7) /   6.0825727250224559D-02  /
   data green(11,10, 8) /   5.9202840734158095D-02  /
   data green(11,10, 9) /   5.7512767294399587D-02  /
   data green(11,10,10) /   5.5785842873145947D-02  /
   data green(11,11, 0) /   6.4266218318071694D-02  /
   data green(11,11, 1) /   6.4133213755402152D-02  /
   data green(11,11, 2) /   6.3739188115049750D-02  /
   data green(11,11, 3) /   6.3098579860750750D-02  /
   data green(11,11, 4) /   6.2233811287682454D-02  /
   data green(11,11, 5) /   6.1173173111136604D-02  /
   data green(11,11, 6) /   5.9948437638828365D-02  /
   data green(11,11, 7) /   5.8592557370455108D-02  /
   data green(11,11, 8) /   5.7137710728801874D-02  /
   data green(11,11, 9) /   5.5613833822526916D-02  /
   data green(11,11,10) /   5.4047665912661157D-02  /
   data green(11,11,11) /   5.2462257678254860D-02  /
   data green(12, 0, 0) /   8.3481055203032714D-02  /
   data green(12, 1, 0) /   8.3186439677825616D-02  /
   data green(12, 1, 1) /   8.2895038892605521D-02  /
   data green(12, 2, 0) /   8.2321888449316913D-02  /
   data green(12, 2, 1) /   8.2039763447092531D-02  /
   data green(12, 2, 2) /   8.1211273633789216D-02  /
   data green(12, 3, 0) /   8.0941596563226664D-02  /
   data green(12, 3, 1) /   8.0673809261550114D-02  /
   data green(12, 3, 2) /   7.9886782095288419D-02  /
   data green(12, 3, 3) /   7.8626641881184400D-02  /
   data green(12, 4, 0) /   7.9125080456315111D-02  /
   data green(12, 4, 1) /   7.8875320043277589D-02  /
   data green(12, 4, 2) /   7.8140522563240736D-02  /
   data green(12, 4, 3) /   7.6961727905604752D-02  /
   data green(12, 4, 4) /   7.5400124032775781D-02  /
   data green(12, 5, 0) /   7.6965080886917517D-02  /
   data green(12, 5, 1) /   7.6735577689252585D-02  /
   data green(12, 5, 2) /   7.6059597819942676D-02  /
   data green(12, 5, 3) /   7.4972777173897873D-02  /
   data green(12, 5, 4) /   7.3528619407139104D-02  /
   data green(12, 5, 5) /   7.1791490632830246D-02  /
   data green(12, 6, 0) /   7.4556241403595236D-02  /
   data green(12, 6, 1) /   7.4347904708403245D-02  /
   data green(12, 6, 2) /   7.3733521274853955D-02  /
   data green(12, 6, 3) /   7.2743430550932242D-02  /
   data green(12, 6, 4) /   7.1423515856441300D-02  /
   data green(12, 6, 5) /   6.9829534414241570D-02  /
   data green(12, 6, 6) /   6.8021243942003498D-02  /
   data green(12, 7, 0) /   7.1986575856202065D-02  /
   data green(12, 7, 1) /   7.1799256916804916D-02  /
   data green(12, 7, 2) /   7.1246175347529583D-02  /
   data green(12, 7, 3) /   7.0352766337072437D-02  /
   data green(12, 7, 4) /   6.9157775202103144D-02  /
   data green(12, 7, 5) /   6.7708764335504906D-02  /
   data green(12, 7, 6) /   6.6057347894089041D-02  /
   data green(12, 7, 7) /   6.4254960311512940D-02  /
   data green(12, 8, 0) /   6.9332309516577331D-02  /
   data green(12, 8, 1) /   6.9165101147695709D-02  /
   data green(12, 8, 2) /   6.8670804245790210D-02  /
   data green(12, 8, 3) /   6.7870496878270486D-02  /
   data green(12, 8, 4) /   6.6796512353697976D-02  /
   data green(12, 8, 5) /   6.5488936183021818D-02  /
   data green(12, 8, 6) /   6.3991805914012892D-02  /
   data green(12, 8, 7) /   6.2349633770772135D-02  /
   data green(12, 8, 8) /   6.0604649928077264D-02  /
   data green(12, 9, 0) /   6.6655700947803737D-02  /
   data green(12, 9, 1) /   6.6507214750507582D-02  /
   data green(12, 9, 2) /   6.6067758269828497D-02  /
   data green(12, 9, 3) /   6.5354654624826425D-02  /
   data green(12, 9, 4) /   6.4394653811605626D-02  /
   data green(12, 9, 5) /   6.3221234032008661D-02  /
   data green(12, 9, 6) /   6.1871611344422754D-02  /
   data green(12, 9, 7) /   6.0383924899015042D-02  /
   data green(12, 9, 8) /   5.8794920641088391D-02  /
   data green(12, 9, 9) /   5.7138284304581823D-02  /
   data green(12,10, 0) /   6.4005024864455293D-02  /
   data green(12,10, 1) /   6.3873617045918346D-02  /
   data green(12,10, 2) /   6.3484283680224740D-02  /
   data green(12,10, 3) /   6.2851183503601693D-02  /
   data green(12,10, 4) /   6.1996316019857840D-02  /
   data green(12,10, 5) /   6.0947459742092466D-02  /
   data green(12,10, 6) /   5.9735842214522342D-02  /
   data green(12,10, 7) /   5.8393888580770721D-02  /
   data green(12,10, 8) /   5.6953304259339789D-02  /
   data green(12,10, 9) /   5.5443628752505386D-02  /
   data green(12,10,10) /   5.3891289652420855D-02  /
   data green(12,11, 0) /   6.1415864434150272D-02  /
   data green(12,11, 1) /   6.1299802908201549D-02  /
   data green(12,11, 2) /   6.0955590558244795D-02  /
   data green(12,11, 3) /   6.0394762423602956D-02  /
   data green(12,11, 4) /   5.9635344407832178D-02  /
   data green(12,11, 5) /   5.8700286421337435D-02  /
   data green(12,11, 6) /   5.7615659408886306D-02  /
   data green(12,11, 7) /   5.6408870677325536D-02  /
   data green(12,11, 8) /   5.5107095855756572D-02  /
   data green(12,11, 9) /   5.3736045852449579D-02  /
   data green(12,11,10) /   5.2319109114446614D-02  /
   data green(12,11,11) /   5.0876849443848571D-02  /
   data green(12,12, 0) /   5.8913026545233144D-02  /
   data green(12,12, 1) /   5.8810604410957687D-02  /
   data green(12,12, 2) /   5.8506561014932144D-02  /
   data green(12,12, 3) /   5.8010280952088406D-02  /
   data green(12,12, 4) /   5.7336507013627164D-02  /
   data green(12,12, 5) /   5.6504152744411049D-02  /
   data green(12,12, 6) /   5.5534913350500176D-02  /
   data green(12,12, 7) /   5.4451860140137968D-02  /
   data green(12,12, 8) /   5.3278170125135089D-02  /
   data green(12,12, 9) /   5.2036089402015234D-02  /
   data green(12,12,10) /   5.0746173602959688D-02  /
   data green(12,12,11) /   4.9426803036115692D-02  /
   data green(12,12,12) /   4.8093939476936430D-02  /
   data green(13, 0, 0) /   7.7038892254143870D-02  /
   data green(13, 1, 0) /   7.6807807792401539D-02  /
   data green(13, 1, 1) /   7.6578856168499090D-02  /
   data green(13, 2, 0) /   7.6127354436534858D-02  /
   data green(13, 2, 1) /   7.5904596329480115D-02  /
   data green(13, 2, 2) /   7.5248331303574167D-02  /
   data green(13, 3, 0) /   7.5033879989309923D-02  /
   data green(13, 3, 1) /   7.4820809902577598D-02  /
   data green(13, 3, 2) /   7.4192720474885565D-02  /
   data green(13, 3, 3) /   7.3181307263292472D-02  /
   data green(13, 4, 0) /   7.3581809086131475D-02  /
   data green(13, 4, 1) /   7.3381117786207026D-02  /
   data green(13, 4, 2) /   7.2789077344646416D-02  /
   data green(13, 4, 3) /   7.1834353830334680D-02  /
   data green(13, 4, 4) /   7.0560362037884683D-02  /
   data green(13, 5, 0) /   7.1836397178280650D-02  /
   data green(13, 5, 1) /   7.1649882936637474D-02  /
   data green(13, 5, 2) /   7.1099192563140901D-02  /
   data green(13, 5, 3) /   7.0209686854276782D-02  /
   data green(13, 5, 4) /   6.9019974561842959D-02  /
   data green(13, 5, 5) /   6.7577416085632697D-02  /
   data green(13, 6, 0) /   6.9866486100106900D-02  /
   data green(13, 6, 1) /   6.9695094710508443D-02  /
   data green(13, 6, 2) /   6.9188587138744692D-02  /
   data green(13, 6, 3) /   6.8368990494838147D-02  /
   data green(13, 6, 4) /   6.7270023577079521D-02  /
   data green(13, 6, 5) /   6.5933369562009139D-02  /
   data green(13, 6, 6) /   6.4404662997527035D-02  /
   data green(13, 7, 0) /   6.7738539288556443D-02  /
   data green(13, 7, 1) /   6.7582489168464197D-02  /
   data green(13, 7, 2) /   6.7120879916130097D-02  /
   data green(13, 7, 3) /   6.6372559730946412D-02  /
   data green(13, 7, 4) /   6.5366543162599433D-02  /
   data green(13, 7, 5) /   6.4138983390169629D-02  /
   data green(13, 7, 6) /   6.2729852406757797D-02  /
   data green(13, 7, 7) /   6.1179858592548680D-02  /
   data green(13, 8, 0) /   6.5512532818734412D-02  /
   data green(13, 8, 1) /   6.5371479103523783D-02  /
   data green(13, 8, 2) /   6.4953834107432398D-02  /
   data green(13, 8, 3) /   6.4275538069933136D-02  /
   data green(13, 8, 4) /   6.3361264980720033D-02  /
   data green(13, 8, 5) /   6.2242002503830753D-02  /
   data green(13, 8, 6) /   6.0952350733635897D-02  /
   data green(13, 8, 7) /   5.9527954226633657D-02  /
   data green(13, 8, 8) /   5.8003360525108169D-02  /
   data green(13, 9, 0) /   6.3239697800972383D-02  /
   data green(13, 9, 1) /   6.3112899478047810D-02  /
   data green(13, 9, 2) /   6.2737115781767175D-02  /
   data green(13, 9, 3) /   6.2125708726380048D-02  /
   data green(13, 9, 4) /   6.1299472030394507D-02  /
   data green(13, 9, 5) /   6.0284720594834890D-02  /
   data green(13, 9, 6) /   5.9111121370494943D-02  /
   data green(13, 9, 7) /   5.7809584068516988D-02  /
   data green(13, 9, 8) /   5.6410449781933075D-02  /
   data green(13, 9, 9) /   5.4942108927762842D-02  /
   data green(13,10, 0) /   6.0961772948278335D-02  /
   data green(13,10, 1) /   6.0848240930676679D-02  /
   data green(13,10, 2) /   6.0511475724927176D-02  /
   data green(13,10, 3) /   5.9962606815814820D-02  /
   data green(13,10, 4) /   5.9219041684743220D-02  /
   data green(13,10, 5) /   5.8302971316011344D-02  /
   data green(13,10, 6) /   5.7239645944158896D-02  /
   data green(13,10, 7) /   5.6055662184196900D-02  /
   data green(13,10, 8) /   5.4777450991759673D-02  /
   data green(13,10, 9) /   5.3430081151840572D-02  /
   data green(13,10,10) /   5.2036419307475157D-02  /
   data green(13,11, 0) /   5.8711311060738298D-02  /
   data green(13,11, 1) /   5.8609927244535835D-02  /
   data green(13,11, 2) /   5.8308945207457677D-02  /
   data green(13,11, 3) /   5.7817595259099804D-02  /
   data green(13,11, 4) /   5.7150383280885314D-02  /
   data green(13,11, 5) /   5.6325929113665900D-02  /
   data green(13,11, 6) /   5.5365606150804020D-02  /
   data green(13,11, 7) /   5.4292162707946526D-02  /
   data green(13,11, 8) /   5.3128473491487534D-02  /
   data green(13,11, 9) /   5.1896518167446587D-02  /
   data green(13,11,10) /   5.0616630222339760D-02  /
   data green(13,11,11) /   4.9307014683711937D-02  /
   data green(13,12, 0) /   5.6512610790591403D-02  /
   data green(13,12, 1) /   5.6422216544588837D-02  /
   data green(13,12, 2) /   5.6153649805028460D-02  /
   data green(13,12, 3) /   5.5714546418382527D-02  /
   data green(13,12, 4) /   5.5116959784799566D-02  /
   data green(13,12, 5) /   5.4376461629832007D-02  /
   data green(13,12, 6) /   5.3511074200339570D-02  /
   data green(13,12, 7) /   5.2540168005270678D-02  /
   data green(13,12, 8) /   5.1483439748572064D-02  /
   data green(13,12, 9) /   5.0360050489548676D-02  /
   data green(13,12,10) /   4.9187965475366933D-02  /
   data green(13,12,11) /   4.7983503124550388D-02  /
   data green(13,12,12) /   4.6761075917757161D-02  /
   data green(13,13, 0) /   5.4382938153103688D-02  /
   data green(13,13, 1) /   5.4302395417780411D-02  /
   data green(13,13, 2) /   5.4062924450219284D-02  /
   data green(13,13, 3) /   5.3670835178183426D-02  /
   data green(13,13, 4) /   5.3136128915233095D-02  /
   data green(13,13, 5) /   5.2471803710961426D-02  /
   data green(13,13, 6) /   5.1693018901680793D-02  /
   data green(13,13, 7) /   5.0816217933487368D-02  /
   data green(13,13, 8) /   4.9858297217904107D-02  /
   data green(13,13, 9) /   4.8835885803172767D-02  /
   data green(13,13,10) /   4.7764773384341569D-02  /
   data green(13,13,11) /   4.6659498802682942D-02  /
   data green(13,13,12) /   4.5533091466540221D-02  /
   data green(13,13,13) /   4.4396945388845144D-02  /
   data green(14, 0, 0) /   7.1521067152224463D-02  /
   data green(14, 1, 0) /   7.1336448048906750D-02  /
   data green(14, 1, 1) /   7.1153289874625353D-02  /
   data green(14, 2, 0) /   7.0791357729220419D-02  /
   data green(14, 2, 1) /   7.0612462718602276D-02  /
   data green(14, 2, 2) /   7.0084080578924096D-02  /
   data green(14, 3, 0) /   6.9910900044811222D-02  /
   data green(14, 3, 1) /   6.9738735574032734D-02  /
   data green(14, 3, 2) /   6.9230014939934290D-02  /
   data green(14, 3, 3) /   6.8407053088129791D-02  /
   data green(14, 4, 0) /   6.8733261320469086D-02  /
   data green(14, 4, 1) /   6.8569807403971464D-02  /
   data green(14, 4, 2) /   6.8086556850069688D-02  /
   data green(14, 4, 3) /   6.7303961576302973D-02  /
   data green(14, 4, 4) /   6.6253391410176901D-02  /
   data green(14, 5, 0) /   6.7305241123252840D-02  /
   data green(14, 5, 1) /   6.7151917692435323D-02  /
   data green(14, 5, 2) /   6.6698323284597266D-02  /
   data green(14, 5, 3) /   6.5962832501606034D-02  /
   data green(14, 5, 4) /   6.4973739409044939D-02  /
   data green(14, 5, 5) /   6.3766317633261502D-02  /
   data green(14, 6, 0) /   6.5677549309192146D-02  /
   data green(14, 6, 1) /   6.5535218210693774D-02  /
   data green(14, 6, 2) /   6.5113843872996527D-02  /
   data green(14, 6, 3) /   6.4429655117049109D-02  /
   data green(14, 6, 4) /   6.3507747824746219D-02  /
   data green(14, 6, 5) /   6.2379596907404286D-02  /
   data green(14, 6, 6) /   6.1080290105624810D-02  /
   data green(14, 7, 0) /   6.3900662619366833D-02  /
   data green(14, 7, 1) /   6.3769685075261526D-02  /
   data green(14, 7, 2) /   6.3381634003663215D-02  /
   data green(14, 7, 3) /   6.2750640405105335D-02  /
   data green(14, 7, 4) /   6.1898652094577671D-02  /
   data green(14, 7, 5) /   6.0853368370315783D-02  /
   data green(14, 7, 6) /   5.9645910832576261D-02  /
   data green(14, 7, 7) /   5.8308577976996144D-02  /
   data green(14, 8, 0) /   6.2021692099252616D-02  /
   data green(14, 8, 1) /   6.1902017900768359D-02  /
   data green(14, 8, 2) /   6.1547187527449872D-02  /
   data green(14, 8, 3) /   6.0969363995849422D-02  /
   data green(14, 8, 4) /   6.0187521763864368D-02  /
   data green(14, 8, 5) /   5.9225759140602213D-02  /
   data green(14, 8, 6) /   5.8111367980085205D-02  /
   data green(14, 8, 7) /   5.6872937669220053D-02  /
   data green(14, 8, 8) /   5.5538705205438237D-02  /
   data green(14, 9, 0) /   6.0082381134922196D-02  /
   data green(14, 9, 1) /   5.9973648500838930D-02  /
   data green(14, 9, 2) /   5.9651018383987127D-02  /
   data green(14, 9, 3) /   5.9124865035543914D-02  /
   data green(14, 9, 4) /   5.8411442384083086D-02  /
   data green(14, 9, 5) /   5.7531522316268978D-02  /
   data green(14, 9, 6) /   5.6508815537399983D-02  /
   data green(14, 9, 7) /   5.5368391948828165D-02  /
   data green(14, 9, 8) /   5.4135273537362347D-02  /
   data green(14, 9, 9) /   5.2833307394447741D-02  /
   data green(14,10, 0) /   5.8118122102150595D-02  /
   data green(14,10, 1) /   5.8019753142352429D-02  /
   data green(14,10, 2) /   5.7727661933621999D-02  /
   data green(14,10, 3) /   5.7250635996137049D-02  /
   data green(14,10, 4) /   5.6602500645495532D-02  /
   data green(14,10, 5) /   5.5801030816572436D-02  /
   data green(14,10, 6) /   5.4866672428923818D-02  /
   data green(14,10, 7) /   5.3821240812028721D-02  /
   data green(14,10, 8) /   5.2686735081760222D-02  /
   data green(14,10, 9) /   5.1484360764093851D-02  /
   data green(14,10,10) /   5.0233803444069323D-02  /
   data green(14,11, 0) /   5.6157763702338125D-02  /
   data green(14,11, 1) /   5.6069046815551922D-02  /
   data green(14,11, 2) /   5.5805432652909714D-02  /
   data green(14,11, 3) /   5.5374327408455458D-02  /
   data green(14,11, 4) /   5.4787429435171273D-02  /
   data green(14,11, 5) /   5.4059865883125294D-02  /
   data green(14,11, 6) /   5.3209165529826900D-02  /
   data green(14,11, 7) /   5.2254195738594224D-02  /
   data green(14,11, 8) /   5.1214173477443621D-02  /
   data green(14,11, 9) /   5.0107827797013352D-02  /
   data green(14,11,10) /   4.8952754585591153D-02  /
   data green(14,11,11) /   4.7764971964642695D-02  /
   data green(14,12, 0) /   5.4223955895295460D-02  /
   data green(14,12, 1) /   5.4144111598639147D-02  /
   data green(14,12, 2) /   5.3906705194590412D-02  /
   data green(14,12, 3) /   5.3517957531095506D-02  /
   data green(14,12, 4) /   5.2987731450837086D-02  /
   data green(14,12, 5) /   5.2328850224814386D-02  /
   data green(14,12, 6) /   5.1556277163949858D-02  /
   data green(14,12, 7) /   5.0686253301140433D-02  /
   data green(14,12, 8) /   4.9735479162086260D-02  /
   data green(14,12, 9) /   4.8720404340718393D-02  /
   data green(14,12,10) /   4.7656662028461472D-02  /
   data green(14,12,11) /   4.6558660831667507D-02  /
   data green(14,12,12) /   4.5439326847425313D-02  /
   data green(14,13, 0) /   5.2333808250939189D-02  /
   data green(14,13, 1) /   5.2262038407090851D-02  /
   data green(14,13, 2) /   5.2048508247647382D-02  /
   data green(14,13, 3) /   5.1698432433441535D-02  /
   data green(14,13, 4) /   5.1220107516327225D-02  /
   data green(14,13, 5) /   5.0624375594595324D-02  /
   data green(14,13, 6) /   4.9923971682478391D-02  /
   data green(14,13, 7) /   4.9132827735398989D-02  /
   data green(14,13, 8) /   4.8265400036810166D-02  /
   data green(14,13, 9) /   4.7336071618846703D-02  /
   data green(14,13,10) /   4.6358662374729825D-02  /
   data green(14,13,11) /   4.5346060858057598D-02  /
   data green(14,13,12) /   4.4309976386099424D-02  /
   data green(14,13,13) /   4.3260799311150087D-02  /
   data green(14,14, 0) /   5.0499689218647954D-02  /
   data green(14,14, 1) /   5.0435211555844829D-02  /
   data green(14,14, 2) /   5.0243266365869094D-02  /
   data green(14,14, 3) /   4.9928221028707316D-02  /
   data green(14,14, 4) /   4.9497046724586005D-02  /
   data green(14,14, 5) /   4.8958897056806686D-02  /
   data green(14,14, 6) /   4.8324590135690376D-02  /
   data green(14,14, 7) /   4.7606048816327330D-02  /
   data green(14,14, 8) /   4.6815750472688170D-02  /
   data green(14,14, 9) /   4.5966227717712964D-02  /
   data green(14,14,10) /   4.5069648039435285D-02  /
   data green(14,14,11) /   4.4137486481467633D-02  /
   data green(14,14,12) /   4.3180293504941999D-02  /
   data green(14,14,13) /   4.2207551238163526D-02  /
   data green(14,14,14) /   4.1227605731002044D-02  /
   data green(15, 0, 0) /   6.6741718463575794D-02  /
   data green(15, 1, 0) /   6.6591875870343431D-02  /
   data green(15, 1, 1) /   6.6443061544234766D-02  /
   data green(15, 2, 0) /   6.6148518407864995D-02  /
   data green(15, 2, 1) /   6.6002716225330840D-02  /
   data green(15, 2, 2) /   6.5571196852897348D-02  /
   data green(15, 3, 0) /   6.5429430081802492D-02  /
   data green(15, 3, 1) /   6.5288418434299458D-02  /
   data green(15, 3, 2) /   6.4870943073563850D-02  /
   data green(15, 3, 3) /   6.4193058233011691D-02  /
   data green(15, 4, 0) /   6.4462002860675827D-02  /
   data green(15, 4, 1) /   6.4327255200614913D-02  /
   data green(15, 4, 2) /   6.3928156550965648D-02  /
   data green(15, 4, 3) /   6.3279581777500440D-02  /
   data green(15, 4, 4) /   6.2404584452480691D-02  /
   data green(15, 5, 0) /   6.3280403546704278D-02  /
   data green(15, 5, 1) /   6.3153033195458924D-02  /
   data green(15, 5, 2) /   6.2775596227334959D-02  /
   data green(15, 5, 3) /   6.2161628768365035D-02  /
   data green(15, 5, 4) /   6.1332172522636112D-02  /
   data green(15, 5, 5) /   6.0313817081037807D-02  /
   data green(15, 6, 0) /   6.1922469591038747D-02  /
   data green(15, 6, 1) /   6.1803218603768646D-02  /
   data green(15, 6, 2) /   6.1449645483230987D-02  /
   data green(15, 6, 3) /   6.0873875603002288D-02  /
   data green(15, 6, 4) /   6.0094819843321615D-02  /
   data green(15, 6, 5) /   5.9136487948468222D-02  /
   data green(15, 6, 6) /   5.8026062156260660D-02  /
   data green(15, 7, 0) /   6.0426824868676607D-02  /
   data green(15, 7, 1) /   6.0316088548612597D-02  /
   data green(15, 7, 2) /   5.9987567225035877D-02  /
   data green(15, 7, 3) /   5.9451977488546200D-02  /
   data green(15, 7, 4) /   5.8726091040834871D-02  /
   data green(15, 7, 5) /   5.7831307856564575D-02  /
   data green(15, 7, 6) /   5.6792008155677677D-02  /
   data green(15, 7, 7) /   5.5633912265835479D-02  /
   data green(15, 8, 0) /   5.8830540220955539D-02  /
   data green(15, 8, 1) /   5.8728415098283464D-02  /
   data green(15, 8, 2) /   5.8425256772998431D-02  /
   data green(15, 8, 3) /   5.7930431173199799D-02  /
   data green(15, 8, 4) /   5.7258647987353667D-02  /
   data green(15, 8, 5) /   5.6428771909247985D-02  /
   data green(15, 8, 6) /   5.5462433967076660D-02  /
   data green(15, 8, 7) /   5.4382628677338916D-02  /
   data green(15, 8, 8) /   5.3212448294356382D-02  /
   data green(15, 9, 0) /   5.7167477805479670D-02  /
   data green(15, 9, 1) /   5.7073820964955599D-02  /
   data green(15, 9, 2) /   5.6795631524700872D-02  /
   data green(15, 9, 3) /   5.6341020726826695D-02  /
   data green(15, 9, 4) /   5.5722772512931111D-02  /
   data green(15, 9, 5) /   5.4957365424121170D-02  /
   data green(15, 9, 6) /   5.4063817179798794D-02  /
   data green(15, 9, 7) /   5.3062500068440856D-02  /
   data green(15, 9, 8) /   5.1974051811766578D-02  /
   data green(15, 9, 9) /   5.0818466836011184D-02  /
   data green(15,10, 0) /   5.5467308542660090D-02  /
   data green(15,10, 1) /   5.5381798571274385D-02  /
   data green(15,10, 2) /   5.5127655701418199D-02  /
   data green(15,10, 3) /   5.4711854137890391D-02  /
   data green(15,10, 4) /   5.4145423556266052D-02  /
   data green(15,10, 5) /   5.3442652309609547D-02  /
   data green(15,10, 6) /   5.2620135933555967D-02  /
   data green(15,10, 7) /   5.1695787489868826D-02  /
   data green(15,10, 8) /   5.0687910907872946D-02  /
   data green(15,10, 9) /   4.9614409715431899D-02  /
   data green(15,10,10) /   4.8492170657570728D-02  /
   data green(15,11, 0) /   5.3755099524740356D-02  /
   data green(15,11, 1) /   5.3677292552536650D-02  /
   data green(15,11, 2) /   5.3445909357780715D-02  /
   data green(15,11, 3) /   5.3066913576319263D-02  /
   data green(15,11, 4) /   5.2549767942228660D-02  /
   data green(15,11, 5) /   5.1906790178700701D-02  /
   data green(15,11, 6) /   5.1152375813564001D-02  /
   data green(15,11, 7) /   5.0302178603179162D-02  /
   data green(15,11, 8) /   4.9372329593720225D-02  /
   data green(15,11, 9) /   4.8378755446786485D-02  /
   data green(15,11,10) /   4.7336632057020198D-02  /
   data green(15,11,11) /   4.6259986269626831D-02  /
   data green(15,12, 0) /   5.2051329430625591D-02  /
   data green(15,12, 1) /   5.1980706734717171D-02  /
   data green(15,12, 2) /   5.1770571332222168D-02  /
   data green(15,12, 3) /   5.1426002280866044D-02  /
   data green(15,12, 4) /   5.0955084212814535D-02  /
   data green(15,12, 5) /   5.0368389626546731D-02  /
   data green(15,12, 6) /   4.9678348045099401D-02  /
   data green(15,12, 7) /   4.8898572001542552D-02  /
   data green(15,12, 8) /   4.8043204067832697D-02  /
   data green(15,12, 9) /   4.7126334948353901D-02  /
   data green(15,12,10) /   4.6161524562695046D-02  /
   data green(15,12,11) /   4.5161440166867461D-02  /
   data green(15,12,12) /   4.4137610731260719D-02  /
   data green(15,13, 0) /   5.0372189424406265D-02  /
   data green(15,13, 1) /   5.0308195371694051D-02  /
   data green(15,13, 2) /   5.0117682634673089D-02  /
   data green(15,13, 3) /   4.9804965098688286D-02  /
   data green(15,13, 4) /   4.9376929979270241D-02  /
   data green(15,13, 5) /   4.8842623376183718D-02  /
   data green(15,13, 6) /   4.8212740563555061D-02  /
   data green(15,13, 7) /   4.7499074642150044D-02  /
   data green(15,13, 8) /   4.6713974013936822D-02  /
   data green(15,13, 9) /   4.5869849439935967D-02  /
   data green(15,13,10) /   4.4978758319236270D-02  /
   data green(15,13,11) /   4.4052080271698352D-02  /
   data green(15,13,12) /   4.3100286328481133D-02  /
   data green(15,13,13) /   4.2132795229529459D-02  /
   data green(15,14, 0) /   4.8730047849159908D-02  /
   data green(15,14, 1) /   4.8672118468185067D-02  /
   data green(15,14, 2) /   4.8499574557683757D-02  /
   data green(15,14, 3) /   4.8216074107452843D-02  /
   data green(15,14, 4) /   4.7827473614673807D-02  /
   data green(15,14, 5) /   4.7341497145686445D-02  /
   data green(15,14, 6) /   4.6767325788025720D-02  /
   data green(15,14, 7) /   4.6115148409393736D-02  /
   data green(15,14, 8) /   4.5395713114234210D-02  /
   data green(15,14, 9) /   4.4619912278230436D-02  /
   data green(15,14,10) /   4.3798424643985777D-02  /
   data green(15,14,11) /   4.2941427814156587D-02  /
   data green(15,14,12) /   4.2058385268910981D-02  /
   data green(15,14,13) /   4.1157904785368460D-02  /
   data green(15,14,14) /   4.0247660223840569D-02  /
   data green(15,15, 0) /   4.7133985983535000D-02  /
   data green(15,15, 1) /   4.7081569291026741D-02  /
   data green(15,15, 2) /   4.6925372109546448D-02  /
   data green(15,15, 3) /   4.6668494031307545D-02  /
   data green(15,15, 4) /   4.6315910613086750D-02  /
   data green(15,15, 5) /   4.5874209461567832D-02  /
   data green(15,15, 6) /   4.5351260178871068D-02  /
   data green(15,15, 7) /   4.4755849294566247D-02  /
   data green(15,15, 8) /   4.4097310784400240D-02  /
   data green(15,15, 9) /   4.3385178486215481D-02  /
   data green(15,15,10) /   4.2628880069532325D-02  /
   data green(15,15,11) /   4.1837484703420851D-02  /
   data green(15,15,12) /   4.1019509461561326D-02  /
   data green(15,15,13) /   4.0182783642530383D-02  /
   data green(15,15,14) /   3.9334365956389790D-02  /
   data green(15,15,15) /   3.8480506956412780D-02  /
   data green(16, 0, 0) /   6.2561740147856565D-02  /
   data green(16, 1, 0) /   6.2438447213243355D-02  /
   data green(16, 1, 1) /   6.2315895199434715D-02  /
   data green(16, 2, 0) /   6.2073014315161990D-02  /
   data green(16, 2, 1) /   6.1952639394085976D-02  /
   data green(16, 2, 2) /   6.1595781940670544D-02  /
   data green(16, 3, 0) /   6.1478322308863874D-02  /
   data green(16, 3, 1) /   6.1361430450451809D-02  /
   data green(16, 3, 2) /   6.1014813430274595D-02  /
   data green(16, 3, 3) /   6.0450246607811536D-02  /
   data green(16, 4, 0) /   6.0674411813849263D-02  /
   data green(16, 4, 1) /   6.0562112758318055D-02  /
   data green(16, 4, 2) /   6.0229006176710079D-02  /
   data green(16, 4, 3) /   5.9686101083876723D-02  /
   data green(16, 4, 4) /   5.8950606704434153D-02  /
   data green(16, 5, 0) /   5.9686641314340191D-02  /
   data green(16, 5, 1) /   5.9579807901895245D-02  /
   data green(16, 5, 2) /   5.9262790081391972D-02  /
   data green(16, 5, 3) /   5.8745713978110257D-02  /
   data green(16, 5, 4) /   5.8044445581566406D-02  /
   data green(16, 5, 5) /   5.7179263910756513D-02  /
   data green(16, 6, 0) /   5.8543602614121892D-02  /
   data green(16, 6, 1) /   5.8442856593576828D-02  /
   data green(16, 6, 2) /   5.8143771367247556D-02  /
   data green(16, 6, 3) /   5.7655526917500803D-02  /
   data green(16, 6, 4) /   5.6992544120241920D-02  /
   data green(16, 6, 5) /   5.6173323846220126D-02  /
   data green(16, 6, 6) /   5.5219090623647327D-02  /
   data green(16, 7, 0) /   5.7275099683473246D-02  /
   data green(16, 7, 1) /   5.7180820420649411D-02  /
   data green(16, 7, 2) /   5.6900801053163752D-02  /
   data green(16, 7, 3) /   5.6443259484098868D-02  /
   data green(16, 7, 4) /   5.5821140974645300D-02  /
   data green(16, 7, 5) /   5.5051119714842917D-02  /
   data green(16, 7, 6) /   5.4152421855752288D-02  /
   data green(16, 7, 7) /   5.3145621951983454D-02  /
   data green(16, 8, 0) /   5.5910416355757310D-02  /
   data green(16, 8, 1) /   5.5822766181410424D-02  /
   data green(16, 8, 2) /   5.5562308089320639D-02  /
   data green(16, 8, 3) /   5.5136320067624332D-02  /
   data green(16, 8, 4) /   5.4556299524043859D-02  /
   data green(16, 8, 5) /   5.3837117079769750D-02  /
   data green(16, 8, 6) /   5.2996009900311218D-02  /
   data green(16, 8, 7) /   5.2051539799104075D-02  /
   data green(16, 8, 8) /   5.1022623582489018D-02  /
   data green(16, 9, 0) /   5.4476995004428756D-02  /
   data green(16, 9, 1) /   5.4395954289156077D-02  /
   data green(16, 9, 2) /   5.4155016598019128D-02  /
   data green(16, 9, 3) /   5.3760569976693998D-02  /
   data green(16, 9, 4) /   5.3222734998936080D-02  /
   data green(16, 9, 5) /   5.2554656316870549D-02  /
   data green(16, 9, 6) /   5.1771652148197941D-02  /
   data green(16, 9, 7) /   5.0890323303043744D-02  /
   data green(16, 9, 8) /   4.9927711226555656D-02  /
   data green(16, 9, 9) /   4.8900570567222350D-02  /
   data green(16,10, 0) /   5.2999554681040915D-02  /
   data green(16,10, 1) /   5.2924960523723211D-02  /
   data green(16,10, 2) /   5.2703078683240832D-02  /
   data green(16,10, 3) /   5.2339475270375052D-02  /
   data green(16,10, 4) /   5.1842993789441016D-02  /
   data green(16,10, 5) /   5.1225167948440876D-02  /
   data green(16,10, 6) /   5.0499510379570824D-02  /
   data green(16,10, 7) /   4.9680758614135366D-02  /
   data green(16,10, 8) /   4.8784151754689892D-02  /
   data green(16,10, 9) /   4.7824793680256016D-02  /
   data green(16,10,10) /   4.6817136956451039D-02  /
   data green(16,11, 0) /   5.1499610367745051D-02  /
   data green(16,11, 1) /   5.1431194311155198D-02  /
   data green(16,11, 2) /   5.1227590331278050D-02  /
   data green(16,11, 3) /   5.0893620330694278D-02  /
   data green(16,11, 4) /   5.0436966694100284D-02  /
   data green(16,11, 5) /   4.9867689497129823D-02  /
   data green(16,11, 6) /   4.9197636591021769D-02  /
   data green(16,11, 7) /   4.8439810981784980D-02  /
   data green(16,11, 8) /   4.7607755059206336D-02  /
   data green(16,11, 9) /   4.6714998569179135D-02  /
   data green(16,11,10) /   4.5774600819850922D-02  /
   data green(16,11,11) /   4.4798801205342853D-02  /
   data green(16,12, 0) /   4.9995318269877270D-02  /
   data green(16,12, 1) /   4.9932740122585410D-02  /
   data green(16,12, 2) /   4.9746421790705403D-02  /
   data green(16,12, 3) /   4.9440521934056354D-02  /
   data green(16,12, 4) /   4.9021683997138472D-02  /
   data green(16,12, 5) /   4.8498641721667630D-02  /
   data green(16,12, 6) /   4.7881733101438104D-02  /
   data green(16,12, 7) /   4.7182373349853521D-02  /
   data green(16,12, 8) /   4.6412534687356510D-02  /
   data green(16,12, 9) /   4.5584271832154737D-02  /
   data green(16,12,10) /   4.4709319856231035D-02  /
   data green(16,12,11) /   4.3798778334256334D-02  /
   data green(16,12,12) /   4.2862884549441815D-02  /
   data green(16,13, 0) /   4.8501560272324572D-02  /
   data green(16,13, 1) /   4.8444436561091930D-02  /
   data green(16,13, 2) /   4.8274281200102770D-02  /
   data green(16,13, 3) /   4.7994669156425994D-02  /
   data green(16,13, 4) /   4.7611326056933703D-02  /
   data green(16,13, 5) /   4.7131807347448464D-02  /
   data green(16,13, 6) /   4.6565099847392981D-02  /
   data green(16,13, 7) /   4.5921185129073047D-02  /
   data green(16,13, 8) /   4.5210602782657740D-02  /
   data green(16,13, 9) /   4.4444045461011536D-02  /
   data green(16,13,10) /   4.3632008622372619D-02  /
   data green(16,13,11) /   4.2784508143390854D-02  /
   data green(16,13,12) /   4.1910870083387793D-02  /
   data green(16,13,13) /   4.1019589838327451D-02  /
   data green(16,14, 0) /   4.7030185532000078D-02  /
   data green(16,14, 1) /   4.6978112183701105D-02  /
   data green(16,14, 2) /   4.6822933660376070D-02  /
   data green(16,14, 3) /   4.6567716313656346D-02  /
   data green(16,14, 4) /   4.6217383117960656D-02  /
   data green(16,14, 5) /   4.5778453576860917D-02  /
   data green(16,14, 6) /   4.5258718294957333D-02  /
   data green(16,14, 7) /   4.4666878801058796D-02  /
   data green(16,14, 8) /   4.4012182726755958D-02  /
   data green(16,14, 9) /   4.3304080267297998D-02  /
   data green(16,14,10) /   4.2551921343323772D-02  /
   data green(16,14,11) /   4.1764705515455865D-02  /
   data green(16,14,12) /   4.0950889718956322D-02  /
   data green(16,14,13) /   4.0118253109235766D-02  /
   data green(16,14,14) /   3.9273814128524744D-02  /
   data green(16,15, 0) /   4.5590340810418399D-02  /
   data green(16,15, 1) /   4.5542910358354415D-02  /
   data green(16,15, 2) /   4.5401510109788351D-02  /
   data green(16,15, 3) /   4.5168766625153199D-02  /
   data green(16,15, 4) /   4.4848906478148695D-02  /
   data green(16,15, 5) /   4.4447545843038186D-02  /
   data green(16,15, 6) /   4.3971425370520470D-02  /
   data green(16,15, 7) /   4.3428113992669326D-02  /
   data green(16,15, 8) /   4.2825705348575514D-02  /
   data green(16,15, 9) /   4.2172527753449714D-02  /
   data green(16,15,10) /   4.1476883969793037D-02  /
   data green(16,15,11) /   4.0746831530559087D-02  /
   data green(16,15,12) /   3.9990008942797389D-02  /
   data green(16,15,13) /   3.9213508433247099D-02  /
   data green(16,15,14) /   3.8423792345857258D-02  /
   data green(16,15,15) /   3.7626647956435742D-02  /
   data green(16,16, 0) /   4.4188837684289509D-02  /
   data green(16,16, 1) /   4.4145651674248318D-02  /
   data green(16,16, 2) /   4.4016855661365192D-02  /
   data green(16,16, 3) /   4.3804698155456782D-02  /
   data green(16,16, 4) /   4.3512805185028113D-02  /
   data green(16,16, 5) /   4.3146010253974978D-02  /
   data green(16,16, 6) /   4.2710138669034052D-02  /
   data green(16,16, 7) /   4.2211764468318508D-02  /
   data green(16,16, 8) /   4.1657958528218166D-02  /
   data green(16,16, 9) /   4.1056044635992960D-02  /
   data green(16,16,10) /   4.0413377008987027D-02  /
   data green(16,16,11) /   3.9737148660461109D-02  /
   data green(16,16,12) /   3.9034235847315843D-02  /
   data green(16,16,13) /   3.8311080112766009D-02  /
   data green(16,16,14) /   3.7573606474798127D-02  /
   data green(16,16,15) /   3.6827174229235551D-02  /
   data green(16,16,16) /   3.6076555602481443D-02  /
   data green(17, 0, 0) /   5.8874933570590059D-02  /
   data green(17, 1, 0) /   5.8772262547228692D-02  /
   data green(17, 1, 1) /   5.8670136437714841D-02  /
   data green(17, 2, 0) /   5.8467519172827699D-02  /
   data green(17, 2, 1) /   5.8366998233465903D-02  /
   data green(17, 2, 2) /   5.8068589142404668D-02  /
   data green(17, 3, 0) /   5.7970216741466039D-02  /
   data green(17, 3, 1) /   5.7872276376868856D-02  /
   data green(17, 3, 2) /   5.7581472332439702D-02  /
   data green(17, 3, 3) /   5.7106592385554704D-02  /
   data green(17, 4, 0) /   5.7295278882579244D-02  /
   data green(17, 4, 1) /   5.7200765006386947D-02  /
   data green(17, 4, 2) /   5.6920063155433416D-02  /
   data green(17, 4, 3) /   5.6461451258693861D-02  /
   data green(17, 4, 4) /   5.5837962443237833D-02  /
   data green(17, 5, 0) /   5.6461814107575149D-02  /
   data green(17, 5, 1) /   5.6371413468285440D-02  /
   data green(17, 5, 2) /   5.6102844837053446D-02  /
   data green(17, 5, 3) /   5.5663791435960522D-02  /
   data green(17, 5, 4) /   5.5066372253219441D-02  /
   data green(17, 5, 5) /   5.4326227527928102D-02  /
   data green(17, 6, 0) /   5.5491691908725389D-02  /
   data green(17, 6, 1) /   5.5405918653778218D-02  /
   data green(17, 6, 2) /   5.5151008074621068D-02  /
   data green(17, 6, 3) /   5.4733996952303854D-02  /
   data green(17, 6, 4) /   5.4166007558197206D-02  /
   data green(17, 6, 5) /   5.3461436759363802D-02  /
   data green(17, 6, 6) /   5.2636990181678292D-02  /
   data green(17, 7, 0) /   5.4408113420984776D-02  /
   data green(17, 7, 1) /   5.4327310463659406D-02  /
   data green(17, 7, 2) /   5.4087079732630099D-02  /
   data green(17, 7, 3) /   5.3693790403112819D-02  /
   data green(17, 7, 4) /   5.3157531963809332D-02  /
   data green(17, 7, 5) /   5.2491406676105308D-02  /
   data green(17, 7, 6) /   5.1710680788051050D-02  /
   data green(17, 7, 7) /   5.0831896136352597D-02  /
   data green(17, 8, 0) /   5.3234330862865366D-02  /
   data green(17, 8, 1) /   5.3158682835295623D-02  /
   data green(17, 8, 2) /   5.2933687862707836D-02  /
   data green(17, 8, 3) /   5.2565052158388637D-02  /
   data green(17, 8, 4) /   5.2061836035933985D-02  /
   data green(17, 8, 5) /   5.1435845113302049D-02  /
   data green(17, 8, 6) /   5.0700894551966698D-02  /
   data green(17, 8, 7) /   4.9872031379458164D-02  /
   data green(17, 8, 8) /   4.8964791149258190D-02  /
   data green(17, 9, 0) /   5.1992612605982612D-02  /
   data green(17, 9, 1) /   5.1922166462399305D-02  /
   data green(17, 9, 2) /   5.1712557025002229D-02  /
   data green(17, 9, 3) /   5.1368852198382366D-02  /
   data green(17, 9, 4) /   5.0899117873636018D-02  /
   data green(17, 9, 5) /   5.0313900301757475D-02  /
   data green(17, 9, 6) /   4.9625595809592409D-02  /
   data green(17, 9, 7) /   4.8847777961828487D-02  /
   data green(17, 9, 8) /   4.7994546340122395D-02  /
   data green(17, 9, 9) /   4.7079946753699406D-02  /
   data green(17,10, 0) /   5.0703493123965042D-02  /
   data green(17,10, 1) /   5.0638182606471821D-02  /
   data green(17,10, 2) /   5.0443773674631383D-02  /
   data green(17,10, 3) /   5.0124734704289797D-02  /
   data green(17,10, 4) /   4.9688194268527740D-02  /
   data green(17,10, 5) /   4.9143505398023334D-02  /
   data green(17,10, 6) /   4.8501711042073022D-02  /
   data green(17,10, 7) /   4.7774967775853730D-02  /
   data green(17,10, 8) /   4.6975981040112878D-02  /
   data green(17,10, 9) /   4.6117494518561640D-02  /
   data green(17,10,10) /   4.5211862080921601D-02  /
   data green(17,11, 0) /   4.9385301108501443D-02  /
   data green(17,11, 1) /   4.9324971969062527D-02  /
   data green(17,11, 2) /   4.9145317514091826D-02  /
   data green(17,11, 3) /   4.8850254178036512D-02  /
   data green(17,11, 4) /   4.8446044607969543D-02  /
   data green(17,11, 5) /   4.7940933848847823D-02  /
   data green(17,11, 6) /   4.7344699772642808D-02  /
   data green(17,11, 7) /   4.6668163673313716D-02  /
   data green(17,11, 8) /   4.5922704768507114D-02  /
   data green(17,11, 9) /   4.5119814576820086D-02  /
   data green(17,11,10) /   4.4270716274881119D-02  /
   data green(17,11,11) /   4.3386062654215643D-02  /
   data green(17,12, 0) /   4.8053929303397538D-02  /
   data green(17,12, 1) /   4.7998363067965448D-02  /
   data green(17,12, 2) /   4.7832825823987679D-02  /
   data green(17,12, 3) /   4.7560734048507587D-02  /
   data green(17,12, 4) /   4.7187563363258370D-02  /
   data green(17,12, 5) /   4.6720546760304106D-02  /
   data green(17,12, 6) /   4.6168299053593577D-02  /
   data green(17,12, 7) /   4.5540404197925863D-02  /
   data green(17,12, 8) /   4.4847001032433879D-02  /
   data green(17,12, 9) /   4.4098397470796385D-02  /
   data green(17,12,10) /   4.3304734966056586D-02  /
   data green(17,12,11) /   4.2475716084264030D-02  /
   data green(17,12,12) /   4.1620399730101863D-02  /
   data green(17,13, 0) /   4.6722795285479859D-02  /
   data green(17,13, 1) /   4.6671730216143570D-02  /
   data green(17,13, 2) /   4.6519543397344092D-02  /
   data green(17,13, 3) /   4.6269204335118538D-02  /
   data green(17,13, 4) /   4.5925482743358476D-02  /
   data green(17,13, 5) /   4.5494699505651735D-02  /
   data green(17,13, 6) /   4.4984414647885773D-02  /
   data green(17,13, 7) /   4.4403081351325195D-02  /
   data green(17,13, 8) /   4.3759694677812482D-02  /
   data green(17,13, 9) /   4.3063459821650128D-02  /
   data green(17,13,10) /   4.2323498611042114D-02  /
   data green(17,13,11) /   4.1548606033174042D-02  /
   data green(17,13,12) /   4.0747061921729139D-02  /
   data green(17,13,13) /   3.9926497415187917D-02  /
   data green(17,14, 0) /   4.5402939892695351D-02  /
   data green(17,14, 1) /   4.5356088522198967D-02  /
   data green(17,14, 2) /   4.5216407597324017D-02  /
   data green(17,14, 3) /   4.4986471214569210D-02  /
   data green(17,14, 4) /   4.4670422667490516D-02  /
   data green(17,14, 5) /   4.4273769712074265D-02  /
   data green(17,14, 6) /   4.3803126385885294D-02  /
   data green(17,14, 7) /   4.3265924253503354D-02  /
   data green(17,14, 8) /   4.2670116046711302D-02  /
   data green(17,14, 9) /   4.2023892041016517D-02  /
   data green(17,14,10) /   4.1335425041037328D-02  /
   data green(17,14,11) /   4.0612654542348477D-02  /
   data green(17,14,12) /   3.9863115393950127D-02  /
   data green(17,14,13) /   3.9093811750990512D-02  /
   data green(17,14,14) /   3.8311133637121811D-02  /
   data green(17,15, 0) /   4.4103214945976923D-02  /
   data green(17,15, 1) /   4.4060278080943559D-02  /
   data green(17,15, 2) /   4.3932222243934609D-02  /
   data green(17,15, 3) /   4.3721274662561184D-02  /
   data green(17,15, 4) /   4.3431027484625156D-02  /
   data green(17,15, 5) /   4.3066269917156001D-02  /
   data green(17,15, 6) /   4.2632775244552702D-02  /
   data green(17,15, 7) /   4.2137060672976380D-02  /
   data green(17,15, 8) /   4.1586138305959711D-02  /
   data green(17,15, 9) /   4.0987273812823258D-02  /
   data green(17,15,10) /   4.0347766112906289D-02  /
   data green(17,15,11) /   3.9674757391426244D-02  /
   data green(17,15,12) /   3.8975078665144595D-02  /
   data green(17,15,13) /   3.8255132446295748D-02  /
   data green(17,15,14) /   3.7520811125747865D-02  /
   data green(17,15,15) /   3.6777447632016104D-02  /
   data green(17,16, 0) /   4.2830520568201830D-02  /
   data green(17,16, 1) /   4.2791197991563119D-02  /
   data green(17,16, 2) /   4.2673881959454631D-02  /
   data green(17,16, 3) /   4.2480497418020575D-02  /
   data green(17,16, 4) /   4.2214154872219793D-02  /
   data green(17,16, 5) /   4.1879012991855660D-02  /
   data green(17,16, 6) /   4.1480103275435026D-02  /
   data green(17,16, 7) /   4.1023130813660033D-02  /
   data green(17,16, 8) /   4.0514265683923176D-02  /
   data green(17,16, 9) /   3.9959938385926590D-02  /
   data green(17,16,10) /   3.9366650408941009D-02  /
   data green(17,16,11) /   3.8740808022672789D-02  /
   data green(17,16,12) /   3.8088584212536518D-02  /
   data green(17,16,13) /   3.7415810743993944D-02  /
   data green(17,16,14) /   3.6727899907536311D-02  /
   data green(17,16,15) /   3.6029793692631096D-02  /
   data green(17,16,16) /   3.5325936975727321D-02  /
   data green(17,17, 0) /   4.1590062056547604D-02  /
   data green(17,17, 1) /   4.1554060340986583D-02  /
   data green(17,17, 2) /   4.1446617655592627D-02  /
   data green(17,17, 3) /   4.1269396863973705D-02  /
   data green(17,17, 4) /   4.1025089779963879D-02  /
   data green(17,17, 5) /   4.0717304807443226D-02  /
   data green(17,17, 6) /   4.0350422769362924D-02  /
   data green(17,17, 7) /   3.9929431892376396D-02  /
   data green(17,17, 8) /   3.9459753450491464D-02  /
   data green(17,17, 9) /   3.8947068879520513D-02  /
   data green(17,17,10) /   3.8397157531699144D-02  /
   data green(17,17,11) /   3.7815752014836973D-02  /
   data green(17,17,12) /   3.7208415626905404D-02  /
   data green(17,17,13) /   3.6580444070910323D-02  /
   data green(17,17,14) /   3.5936791638475636D-02  /
   data green(17,17,15) /   3.5282020504753001D-02  /
   data green(17,17,16) /   3.4620270714305165D-02  /
   data green(17,17,17) /   3.3955247824246193D-02  /
   data green(18, 0, 0) /   5.5598811086768958D-02  /
   data green(18, 1, 0) /   5.5512402451281154D-02  /
   data green(18, 1, 1) /   5.5426401865465241D-02  /
   data green(18, 2, 0) /   5.5255624966177586D-02  /
   data green(18, 2, 1) /   5.5170828876893133D-02  /
   data green(18, 2, 2) /   5.4918811781946304D-02  /
   data green(18, 3, 0) /   5.4835627330787555D-02  /
   data green(18, 3, 1) /   5.4752775497002734D-02  /
   data green(18, 3, 2) /   5.4506499504101989D-02  /
   data green(18, 3, 3) /   5.4103460523249436D-02  /
   data green(18, 4, 0) /   5.4263700619811124D-02  /
   data green(18, 4, 1) /   5.4183445396836671D-02  /
   data green(18, 4, 2) /   5.3944839321316322D-02  /
   data green(18, 4, 3) /   5.3554197048099710D-02  /
   data green(18, 4, 4) /   5.3021520931166133D-02  /
   data green(18, 5, 0) /   5.3554445603235169D-02  /
   data green(18, 5, 1) /   5.3477330303023755D-02  /
   data green(18, 5, 2) /   5.3248002806628073D-02  /
   data green(18, 5, 3) /   5.2872369336382369D-02  /
   data green(18, 5, 4) /   5.2359798972426506D-02  /
   data green(18, 5, 5) /   5.1722483231645473D-02  /
   data green(18, 6, 0) /   5.2724784002113931D-02  /
   data green(18, 6, 1) /   5.2651231413801222D-02  /
   data green(18, 6, 2) /   5.2432436904713040D-02  /
   data green(18, 6, 3) /   5.2073857324430484D-02  /
   data green(18, 6, 4) /   5.1584163300244894D-02  /
   data green(18, 6, 5) /   5.0974664766053765D-02  /
   data green(18, 6, 6) /   5.0258615371592583D-02  /
   data green(18, 7, 0) /   5.1792937143525730D-02  /
   data green(18, 7, 1) /   5.1723247718618380D-02  /
   data green(18, 7, 2) /   5.1515880538515677D-02  /
   data green(18, 7, 3) /   5.1175822127446029D-02  /
   data green(18, 7, 4) /   5.0711010134861484D-02  /
   data green(18, 7, 5) /   5.0131825717262643D-02  /
   data green(18, 7, 6) /   4.9450475326450992D-02  /
   data green(18, 7, 7) /   4.8680330531113411D-02  /
   data green(18, 8, 0) /   5.0777477123876154D-02  /
   data green(18, 8, 1) /   5.0711835178731089D-02  /
   data green(18, 8, 2) /   5.0516447335609789D-02  /
   data green(18, 8, 3) /   5.0195826382570523D-02  /
   data green(18, 8, 4) /   4.9757169582423912D-02  /
   data green(18, 8, 5) /   4.9209916034279017D-02  /
   data green(18, 8, 6) /   4.8565204321085076D-02  /
   data green(18, 8, 7) /   4.7835288678238913D-02  /
   data green(18, 8, 8) /   4.7032967873898900D-02  /
   data green(18, 9, 0) /   4.9696523719837835D-02  /
   data green(18, 9, 1) /   4.9635009412487482D-02  /
   data green(18, 9, 2) /   4.9451845353873224D-02  /
   data green(18, 9, 3) /   4.9151081452283508D-02  /
   data green(18, 9, 4) /   4.8739189438278367D-02  /
   data green(18, 9, 5) /   4.8224681299120421D-02  /
   data green(18, 9, 6) /   4.7617638851730718D-02  /
   data green(18, 9, 7) /   4.6929203162900522D-02  /
   data green(18, 9, 8) /   4.6171069909279674D-02  /
   data green(18, 9, 9) /   4.5355028237605996D-02  /
   data green(18,10, 0) /   4.8567125135321562D-02  /
   data green(18,10, 1) /   4.8509729885049348D-02  /
   data green(18,10, 2) /   4.8338771514691792D-02  /
   data green(18,10, 3) /   4.8057858635785683D-02  /
   data green(18,10, 4) /   4.7672769237806646D-02  /
   data green(18,10, 5) /   4.7191124978493065D-02  /
   data green(18,10, 6) /   4.6621987116755123D-02  /
   data green(18,10, 7) /   4.5975414323658856D-02  /
   data green(18,10, 8) /   4.5262021082942451D-02  /
   data green(18,10, 9) /   4.4492568979559784D-02  /
   data green(18,10,10) /   4.3677613942192842D-02  /
   data green(18,11, 0) /   4.7404829917775927D-02  /
   data green(18,11, 1) /   4.7351473206189550D-02  /
   data green(18,11, 2) /   4.7192489107554733D-02  /
   data green(18,11, 3) /   4.6931073862665977D-02  /
   data green(18,11, 4) /   4.6572355295046393D-02  /
   data green(18,11, 5) /   4.6123117057017111D-02  /
   data green(18,11, 6) /   4.5591454427006801D-02  /
   data green(18,11, 7) /   4.4986394543167954D-02  /
   data green(18,11, 8) /   4.4317513229155656D-02  /
   data green(18,11, 9) /   4.3594575858995205D-02  /
   data green(18,11,10) /   4.2827222551421948D-02  /
   data green(18,11,11) /   4.2024710001660644D-02  /
   data green(18,12, 0) /   4.6223435020595692D-02  /
   data green(18,12, 1) /   4.6173980972615776D-02  /
   data green(18,12, 2) /   4.6026575125311706D-02  /
   data green(18,12, 3) /   4.5784034656588118D-02  /
   data green(18,12, 4) /   4.5450888022271818D-02  /
   data green(18,12, 5) /   4.5033143068805528D-02  /
   data green(18,12, 6) /   4.4537995841364642D-02  /
   data green(18,12, 7) /   4.3973506727035323D-02  /
   data green(18,12, 8) /   4.3348270392204721D-02  /
   data green(18,12, 9) /   4.2671102595662883D-02  /
   data green(18,12,10) /   4.1950761499124717D-02  /
   data green(18,12,11) /   4.1195714787779669D-02  /
   data green(18,12,12) /   4.0413957819183000D-02  /
   data green(18,13, 0) /   4.5034881812174660D-02  /
   data green(18,13, 1) /   4.4989154462781046D-02  /
   data green(18,13, 2) /   4.4852811211387851D-02  /
   data green(18,13, 3) /   4.4628325450505719D-02  /
   data green(18,13, 4) /   4.4319680531698012D-02  /
   data green(18,13, 5) /   4.3932175822839307D-02  /
   data green(18,13, 6) /   4.3472181743465850D-02  /
   data green(18,13, 7) /   4.2946865201750457D-02  /
   data green(18,13, 8) /   4.2363907032096623D-02  /
   data green(18,13, 9) /   4.1731230668221597D-02  /
   data green(18,13,10) /   4.1056757181961952D-02  /
   data green(18,13,11) /   4.0348196896269582D-02  /
   data green(18,13,12) /   3.9612882874573055D-02  /
   data green(18,13,13) /   3.8857647307360688D-02  /
   data green(18,14, 0) /   4.3849266548317202D-02  /
   data green(18,14, 1) /   4.3807063162386153D-02  /
   data green(18,14, 2) /   4.3681186545521035D-02  /
   data green(18,14, 3) /   4.3473801711264602D-02  /
   data green(18,14, 4) /   4.3188401732355247D-02  /
   data green(18,14, 5) /   4.2829646219401114D-02  /
   data green(18,14, 6) /   4.2403156162337981D-02  /
   data green(18,14, 7) /   4.1915282263063187D-02  /
   data green(18,14, 8) /   4.1372864275219752D-02  /
   data green(18,14, 9) /   4.0782997256578996D-02  /
   data green(18,14,10) /   4.0152817595586761D-02  /
   data green(18,14,11) /   3.9489317879047144D-02  /
   data green(18,14,12) /   3.8799195764531558D-02  /
   data green(18,14,13) /   3.8088738504151766D-02  /
   data green(18,14,14) /   3.7363741940729382D-02  /
   data green(18,15, 0) /   4.2674932242658614D-02  /
   data green(18,15, 1) /   4.2636034423633340D-02  /
   data green(18,15, 2) /   4.2519981053535896D-02  /
   data green(18,15, 3) /   4.2328662989064271D-02  /
   data green(18,15, 4) /   4.2065136310174708D-02  /
   data green(18,15, 5) /   4.1733488226789947D-02  /
   data green(18,15, 6) /   4.1338665845343077D-02  /
   data green(18,15, 7) /   4.0886281427513133D-02  /
   data green(18,15, 8) /   4.0382408270748495D-02  /
   data green(18,15, 9) /   3.9833380277680693D-02  /
   data green(18,15,10) /   3.9245606053416780D-02  /
   data green(18,15,11) /   3.8625405474926197D-02  /
   data green(18,15,12) /   3.7978873603861506D-02  /
   data green(18,15,13) /   3.7311773958650485D-02  /
   data green(18,15,14) /   3.6629460785659838D-02  /
   data green(18,15,15) /   3.5936828198902968D-02  /
   data green(18,16, 0) /   4.1518612739270655D-02  /
   data green(18,16, 1) /   4.1482795325543775D-02  /
   data green(18,16, 2) /   4.1375900786251585D-02  /
   data green(18,16, 3) /   4.1199577994457372D-02  /
   data green(18,16, 4) /   4.0956496373825278D-02  /
   data green(18,16, 5) /   4.0650234824904929D-02  /
   data green(18,16, 6) /   4.0285139164144688D-02  /
   data green(18,16, 7) /   3.9866158888930775D-02  /
   data green(18,16, 8) /   3.9398674619282946D-02  /
   data green(18,16, 9) /   3.8888326894358567D-02  /
   data green(18,16,10) /   3.8340855391454209D-02  /
   data green(18,16,11) /   3.7761955447834769D-02  /
   data green(18,16,12) /   3.7157156369186449D-02  /
   data green(18,16,13) /   3.6531723713918739D-02  /
   data green(18,16,14) /   3.5890585769571139D-02  /
   data green(18,16,15) /   3.5238282906518259D-02  /
   data green(18,16,16) /   3.4578937438009931D-02  /
   data green(18,17, 0) /   4.0385605260454399D-02  /
   data green(18,17, 1) /   4.0352643175326434D-02  /
   data green(18,17, 2) /   4.0254242404843099D-02  /
   data green(18,17, 3) /   4.0091839456122207D-02  /
   data green(18,17, 4) /   3.9867763631963329D-02  /
   data green(18,17, 5) /   3.9585145155810575D-02  /
   data green(18,17, 6) /   3.9247796682743905D-02  /
   data green(18,17, 7) /   3.8860076753541917D-02  /
   data green(18,17, 8) /   3.8426744281602768D-02  /
   data green(18,17, 9) /   3.7952812761539584D-02  /
   data green(18,17,10) /   3.7443411738977639D-02  /
   data green(18,17,11) /   3.6903661442535811D-02  /
   data green(18,17,12) /   3.6338564625306047D-02  /
   data green(18,17,13) /   3.5752917836046436D-02  /
   data green(18,17,14) /   3.5151242720058901D-02  /
   data green(18,17,15) /   3.4537736645237467D-02  /
   data green(18,17,16) /   3.3916241002688606D-02  /
   data green(18,17,17) /   3.3290224934187186D-02  /
   data green(18,18, 0) /   3.9279953469197305D-02  /
   data green(18,18, 1) /   3.9249626780149795D-02  /
   data green(18,18, 2) /   3.9159069173013601D-02  /
   data green(18,18, 3) /   3.9009531614803188D-02  /
   data green(18,18, 4) /   3.8803045597498816D-02  /
   data green(18,18, 5) /   3.8542347197556989D-02  /
   data green(18,18, 6) /   3.8230778691048399D-02  /
   data green(18,18, 7) /   3.7872174486157734D-02  /
   data green(18,18, 8) /   3.7470738636276348D-02  /
   data green(18,18, 9) /   3.7030920980509641D-02  /
   data green(18,18,10) /   3.6557298149146433D-02  /
   data green(18,18,11) /   3.6054464453977461D-02  /
   data green(18,18,12) /   3.5526936259992155D-02  /
   data green(18,18,13) /   3.4979071991671606D-02  /
   data green(18,18,14) /   3.4415008609302425D-02  /
   data green(18,18,15) /   3.3838614293498132D-02  /
   data green(18,18,16) /   3.3253456245018756D-02  /
   data green(18,18,17) /   3.2662781946532840D-02  /
   data green(18,18,18) /   3.2069511918433713D-02  /
   data green(19, 0, 0) /   5.2668323144899130D-02  /
   data green(19, 1, 0) /   5.2594912376581010D-02  /
   data green(19, 1, 1) /   5.2521812098727964D-02  /
   data green(19, 2, 0) /   5.2376543088382792D-02  /
   data green(19, 2, 1) /   5.2304360917300868D-02  /
   data green(19, 2, 2) /   5.2089624787282830D-02  /
   data green(19, 3, 0) /   5.2018670663264911D-02  /
   data green(19, 3, 1) /   5.1947975496527986D-02  /
   data green(19, 3, 2) /   5.1737637568657673D-02  /
   data green(19, 3, 3) /   5.1392777649538264D-02  /
   data green(19, 4, 0) /   5.1529961108591246D-02  /
   data green(19, 4, 1) /   5.1461261633479166D-02  /
   data green(19, 4, 2) /   5.1256828016798144D-02  /
   data green(19, 4, 3) /   5.0921540917789991D-02  /
   data green(19, 4, 4) /   5.0463171102706310D-02  /
   data green(19, 5, 0) /   5.0921714283863134D-02  /
   data green(19, 5, 1) /   5.0855443129355617D-02  /
   data green(19, 5, 2) /   5.0658196137279336D-02  /
   data green(19, 5, 3) /   5.0334568430749085D-02  /
   data green(19, 5, 4) /   4.9891884781900826D-02  /
   data green(19, 5, 5) /   4.9339744384983945D-02  /
   data green(19, 6, 0) /   5.0207166529594852D-02  /
   data green(19, 6, 1) /   5.0143670900560017D-02  /
   data green(19, 6, 2) /   4.9954641270281784D-02  /
   data green(19, 6, 3) /   4.9644355325492523D-02  /
   data green(19, 6, 4) /   4.9219641045499835D-02  /
   data green(19, 6, 5) /   4.8689464216213223D-02  /
   data green(19, 6, 6) /   4.8064421862867655D-02  /
   data green(19, 7, 0) /   4.9400752876942594D-02  /
   data green(19, 7, 1) /   4.9340291573476983D-02  /
   data green(19, 7, 2) /   4.9160249412076751D-02  /
   data green(19, 7, 3) /   4.8864567981185117D-02  /
   data green(19, 7, 4) /   4.8459548117981914D-02  /
   data green(19, 7, 5) /   4.7953481365825286D-02  /
   data green(19, 7, 6) /   4.7356195229601078D-02  /
   data green(19, 7, 7) /   4.6678559040948003D-02  /
   data green(19, 8, 0) /   4.8517400659565654D-02  /
   data green(19, 8, 1) /   4.8460146646499107D-02  /
   data green(19, 8, 2) /   4.8289608697260268D-02  /
   data green(19, 8, 3) /   4.8009385639271075D-02  /
   data green(19, 8, 4) /   4.7625239468620052D-02  /
   data green(19, 8, 5) /   4.7144770172480005D-02  /
   data green(19, 8, 6) /   4.6577012483844139D-02  /
   data green(19, 8, 7) /   4.5931994800977316D-02  /
   data green(19, 8, 8) /   4.5220298919097950D-02  /
   data green(19, 9, 0) /   4.7571908628871631D-02  /
   data green(19, 9, 1) /   4.7517955875033963D-02  /
   data green(19, 9, 2) /   4.7357205383164430D-02  /
   data green(19, 9, 3) /   4.7092916704810944D-02  /
   data green(19, 9, 4) /   4.6730317184399635D-02  /
   data green(19, 9, 5) /   4.6276318190183727D-02  /
   data green(19, 9, 6) /   4.5739161400930713D-02  /
   data green(19, 9, 7) /   4.5128029252666610D-02  /
   data green(19, 9, 8) /   4.4452652758207728D-02  /
   data green(19, 9, 9) /   4.3722944894257758D-02  /
   data green(19,10, 0) /   4.6578444280251832D-02  /
   data green(19,10, 1) /   4.6527817435133755D-02  /
   data green(19,10, 2) /   4.6376932431256380D-02  /
   data green(19,10, 3) /   4.6128721042976467D-02  /
   data green(19,10, 4) /   4.5787892810708497D-02  /
   data green(19,10, 5) /   4.5360689762235715D-02  /
   data green(19,10, 6) /   4.4854579077568311D-02  /
   data green(19,10, 7) /   4.4277912260333334D-02  /
   data green(19,10, 8) /   4.3639579015887914D-02  /
   data green(19,10, 9) /   4.2948680237060465D-02  /
   data green(19,10,10) /   4.2214238510271826D-02  /
   data green(19,11, 0) /   4.5550171829857544D-02  /
   data green(19,11, 1) /   4.5502837385078793D-02  /
   data green(19,11, 2) /   4.5361723422377362D-02  /
   data green(19,11, 3) /   4.5129451302984333D-02  /
   data green(19,11, 4) /   4.4810238982547330D-02  /
   data green(19,11, 5) /   4.4409690723902981D-02  /
   data green(19,11, 6) /   4.3934532231093164D-02  /
   data green(19,11, 7) /   4.3392314878837268D-02  /
   data green(19,11, 8) /   4.2791112729116461D-02  /
   data green(19,11, 9) /   4.2139233221350852D-02  /
   data green(19,11,10) /   4.1444957731677216D-02  /
   data green(19,11,11) /   4.0716322677344916D-02  /
   data green(19,12, 0) /   4.4499007144994322D-02  /
   data green(19,12, 1) /   4.4454884952851326D-02  /
   data green(19,12, 2) /   4.4323308979736027D-02  /
   data green(19,12, 3) /   4.4106611430265827D-02  /
   data green(19,12, 4) /   4.3808551218234383D-02  /
   data green(19,12, 5) /   4.3434134909036054D-02  /
   data green(19,12, 6) /   4.2989390028280734D-02  /
   data green(19,12, 7) /   4.2481110185798227D-02  /
   data green(19,12, 8) /   4.1916591736843888D-02  /
   data green(19,12, 9) /   4.1303379684634565D-02  /
   data green(19,12,10) /   4.0649036911865159D-02  /
   data green(19,12,11) /   3.9960946426176057D-02  /
   data green(19,12,12) /   3.9246151860290575D-02  /
   data green(19,13, 0) /   4.3435485093093112D-02  /
   data green(19,13, 1) /   4.3394459401978239D-02  /
   data green(19,13, 2) /   4.3272082286863509D-02  /
   data green(19,13, 3) /   4.3070420256977091D-02  /
   data green(19,13, 4) /   4.2792809418857436D-02  /
   data green(19,13, 5) /   4.2443703868013660D-02  /
   data green(19,13, 6) /   4.2028482784616870D-02  /
   data green(19,13, 7) /   4.1553232093516440D-02  /
   data green(19,13, 8) /   4.1024516977370598D-02  /
   data green(19,13, 9) /   4.0449160121378946D-02  /
   data green(19,13,10) /   3.9834037821526835D-02  /
   data green(19,13,11) /   3.9185902621081618D-02  /
   data green(19,13,12) /   3.8511237538364412D-02  /
   data green(19,13,13) /   3.7816143671526470D-02  /
   data green(19,14, 0) /   4.2368719017602287D-02  /
   data green(19,14, 1) /   4.2330648471865184D-02  /
   data green(19,14, 2) /   4.2217054537227655D-02  /
   data green(19,14, 3) /   4.2029762335572075D-02  /
   data green(19,14, 4) /   4.1771722954102461D-02  /
   data green(19,14, 5) /   4.1446885660009813D-02  /
   data green(19,14, 6) /   4.1060034516869293D-02  /
   data green(19,14, 7) /   4.0616602258054420D-02  /
   data green(19,14, 8) /   4.0122474782119895D-02  /
   data green(19,14, 9) /   3.9583798681748647D-02  /
   data green(19,14,10) /   3.9006802161478275D-02  /
   data green(19,14,11) /   3.8397637000620245D-02  /
   data green(19,14,12) /   3.7762246331605714D-02  /
   data green(19,14,13) /   3.7106260302647902D-02  /
   data green(19,14,14) /   3.6434919427631239D-02  /
   data green(19,15, 0) /   4.1306430375120294D-02  /
   data green(19,15, 1) /   4.1271156676731958D-02  /
   data green(19,15, 2) /   4.1165879332810518D-02  /
   data green(19,15, 3) /   4.0992206202391746D-02  /
   data green(19,15, 4) /   4.0752741054388419D-02  /
   data green(19,15, 5) /   4.0450976243759824D-02  /
   data green(19,15, 6) /   4.0091154846720004D-02  /
   data green(19,15, 7) /   3.9678112624016749D-02  /
   data green(19,15, 8) /   3.9217110718834727D-02  /
   data green(19,15, 9) /   3.8713669378157962D-02  /
   data green(19,15,10) /   3.8173411467946471D-02  /
   data green(19,15,11) /   3.7601922473970630D-02  /
   data green(19,15,12) /   3.7004631390729723D-02  /
   data green(19,15,13) /   3.6386714698015347D-02  /
   data green(19,15,14) /   3.5753023720357192D-02  /
   data green(19,15,15) /   3.5108034177045437D-02  /
   data green(19,16, 0) /   4.0255027713221564D-02  /
   data green(19,16, 1) /   4.0222382837217366D-02  /
   data green(19,16, 2) /   4.0124925979674736D-02  /
   data green(19,16, 3) /   3.9964070943460279D-02  /
   data green(19,16, 4) /   3.9742110609231426D-02  /
   data green(19,16, 5) /   3.9462127037003351D-02  /
   data green(19,16, 6) /   3.9127875467510412D-02  /
   data green(19,16, 7) /   3.8743650555441919D-02  /
   data green(19,16, 8) /   3.8314143696133666D-02  /
   data green(19,16, 9) /   3.7844299930267326D-02  /
   data green(19,16,10) /   3.7339181805196450D-02  /
   data green(19,16,11) /   3.6803845986233076D-02  /
   data green(19,16,12) /   3.6243236611661088D-02  /
   data green(19,16,13) /   3.5662097606063863D-02  /
   data green(19,16,14) /   3.5064904583734714D-02  /
   data green(19,16,15) /   3.4455815694756899D-02  /
   data green(19,16,16) /   3.3838639833918452D-02  /
   data green(19,17, 0) /   3.9219716947057906D-02  /
   data green(19,17, 1) /   3.9189528939980948D-02  /
   data green(19,17, 2) /   3.9099384182752162D-02  /
   data green(19,17, 3) /   3.8950524226692533D-02  /
   data green(19,17, 4) /   3.8744965428559118D-02  /
   data green(19,17, 5) /   3.8485423790251050D-02  /
   data green(19,17, 6) /   3.8175217562554813D-02  /
   data green(19,17, 7) /   3.7818154289439944D-02  /
   data green(19,17, 8) /   3.7418409468484223D-02  /
   data green(19,17, 9) /   3.6980403794268417D-02  /
   data green(19,17,10) /   3.6508685157490624D-02  /
   data green(19,17,11) /   3.6007820374281907D-02  /
   data green(19,17,12) /   3.5482300217170006D-02  /
   data green(19,17,13) /   3.4936459894413067D-02  /
   data green(19,17,14) /   3.4374415821720253D-02  /
   data green(19,17,15) /   3.3800018443907877D-02  /
   data green(19,17,16) /   3.3216820039373068D-02  /
   data green(19,17,17) /   3.2628055882433928D-02  /
   data green(19,18, 0) /   3.8204628366987695D-02  /
   data green(19,18, 1) /   3.8176725844017483D-02  /
   data green(19,18, 2) /   3.8093385916721199D-02  /
   data green(19,18, 3) /   3.7955698003031113D-02  /
   data green(19,18, 4) /   3.7765433754353026D-02  /
   data green(19,18, 5) /   3.7524984292586765D-02  /
   data green(19,18, 6) /   3.7237278548543924D-02  /
   data green(19,18, 7) /   3.6905688041268842D-02  /
   data green(19,18, 8) /   3.6533923893892595D-02  /
   data green(19,18, 9) /   3.6125931787277923D-02  /
   data green(19,18,10) /   3.5685789991218062D-02  /
   data green(19,18,11) /   3.5217614714827275D-02  /
   data green(19,18,12) /   3.4725475931845856D-02  /
   data green(19,18,13) /   3.4213325704321035D-02  /
   data green(19,18,14) /   3.3684939965263477D-02  /
   data green(19,18,15) /   3.3143873808096316D-02  /
   data green(19,18,16) /   3.2593429612352141D-02  /
   data green(19,18,17) /   3.2036636823685079D-02  /
   data green(19,18,18) /   3.1476241890334544D-02  /
   data green(19,19, 0) /   3.7212949332594600D-02  /
   data green(19,19, 1) /   3.7187164835181648D-02  /
   data green(19,19, 2) /   3.7110133613872449D-02  /
   data green(19,19, 3) /   3.6982811252481269D-02  /
   data green(19,19, 4) /   3.6806753731557097D-02  /
   data green(19,19, 5) /   3.6584065049335861D-02  /
   data green(19,19, 6) /   3.6317328808598927D-02  /
   data green(19,19, 7) /   3.6009528034620139D-02  /
   data green(19,19, 8) /   3.5663957897002642D-02  /
   data green(19,19, 9) /   3.5284135988917512D-02  /
   data green(19,19,10) /   3.4873714427077811D-02  /
   data green(19,19,11) /   3.4436397368034985D-02  /
   data green(19,19,12) /   3.3975866701409209D-02  /
   data green(19,19,13) /   3.3495717786836285D-02  /
   data green(19,19,14) /   3.2999406240858967D-02  /
   data green(19,19,15) /   3.2490206020413540D-02  /
   data green(19,19,16) /   3.1971178431909941D-02  /
   data green(19,19,17) /   3.1445151235591502D-02  /
   data green(19,19,18) /   3.0914706710725007D-02  /
   data green(19,19,19) /   3.0382177381536670D-02  /
   data green(20, 0, 0) /   5.0031478276301680D-02  /
   data green(20, 1, 0) /   4.9968581377467194D-02  /
   data green(20, 1, 1) /   4.9905924139138726D-02  /
   data green(20, 2, 0) /   4.9781328689671347D-02  /
   data green(20, 2, 1) /   4.9719381155990401D-02  /
   data green(20, 2, 2) /   4.9534939969833255D-02  /
   data green(20, 3, 0) /   4.9473941458771813D-02  /
   data green(20, 3, 1) /   4.9413146699462535D-02  /
   data green(20, 3, 2) /   4.9232119957696732D-02  /
   data green(20, 3, 3) /   4.8934848401084668D-02  /
   data green(20, 4, 0) /   4.9053157817650003D-02  /
   data green(20, 4, 1) /   4.8993916652765181D-02  /
   data green(20, 4, 2) /   4.8817492522266449D-02  /
   data green(20, 4, 3) /   4.8527703258848413D-02  /
   data green(20, 4, 4) /   4.8130654285689953D-02  /
   data green(20, 5, 0) /   4.8527826203046256D-02  /
   data green(20, 5, 1) /   4.8470485517928305D-02  /
   data green(20, 5, 2) /   4.8299693169863149D-02  /
   data green(20, 5, 3) /   4.8019064085649363D-02  /
   data green(20, 5, 4) /   4.7634384755039352D-02  /
   data green(20, 5, 5) /   4.7153285082029886D-02  /
   data green(20, 6, 0) /   4.7908406213529570D-02  /
   data green(20, 6, 1) /   4.7853251486849574D-02  /
   data green(20, 6, 2) /   4.7688938941569010D-02  /
   data green(20, 6, 3) /   4.7418855900808729D-02  /
   data green(20, 6, 4) /   4.7048430531165600D-02  /
   data green(20, 6, 5) /   4.6584831872048836D-02  /
   data green(20, 6, 6) /   4.6036596913984741D-02  /
   data green(20, 7, 0) /   4.7206429782015266D-02  /
   data green(20, 7, 1) /   4.7153681886164839D-02  /
   data green(20, 7, 2) /   4.6996506411740559D-02  /
   data green(20, 7, 3) /   4.6738047264842486D-02  /
   data green(20, 7, 4) /   4.6383348711899509D-02  /
   data green(20, 7, 5) /   4.5939084741412825D-02  /
   data green(20, 7, 6) /   4.5413221297670295D-02  /
   data green(20, 7, 7) /   4.4814643653417842D-02  /
   data green(20, 8, 0) /   4.6433972034705263D-02  /
   data green(20, 8, 1) /   4.6383787894926649D-02  /
   data green(20, 8, 2) /   4.6234217766940595D-02  /
   data green(20, 8, 3) /   4.5988154601179429D-02  /
   data green(20, 8, 4) /   4.5650246115163676D-02  /
   data green(20, 8, 5) /   4.5226653492127858D-02  /
   data green(20, 8, 6) /   4.4724748977990551D-02  /
   data green(20, 8, 7) /   4.4152780431068205D-02  /
   data green(20, 8, 8) /   4.3519530522318696D-02  /
   data green(20, 9, 0) /   4.5603171363506459D-02  /
   data green(20, 9, 1) /   4.5555647744047349D-02  /
   data green(20, 9, 2) /   4.5413973236751545D-02  /
   data green(20, 9, 3) /   4.5180789526770480D-02  /
   data green(20, 9, 4) /   4.4860346493796180D-02  /
   data green(20, 9, 5) /   4.4458289321935705D-02  /
   data green(20, 9, 6) /   4.3981390549590639D-02  /
   data green(20, 9, 7) /   4.3437251132649964D-02  /
   data green(20, 9, 8) /   4.2833994558056247D-02  /
   data green(20, 9, 9) /   4.2179975137404055D-02  /
   data green(20,10, 0) /   4.4725825217192705D-02  /
   data green(20,10, 1) /   4.4681004754132461D-02  /
   data green(20,10, 2) /   4.4547355789194297D-02  /
   data green(20,10, 3) /   4.4327274372312986D-02  /
   data green(20,10, 4) /   4.4024620654285140D-02  /
   data green(20,10, 5) /   4.3644532792557623D-02  /
   data green(20,10, 6) /   4.3193191694061731D-02  /
   data green(20,10, 7) /   4.2677557013381065D-02  /
   data green(20,10, 8) /   4.2105095027666996D-02  /
   data green(20,10, 9) /   4.1483516801034619D-02  /
   data green(20,10,10) /   4.0820541181962572D-02  /
   data green(20,11, 0) /   4.3813074865965529D-02  /
   data green(20,11, 1) /   4.3770953464063338D-02  /
   data green(20,11, 2) /   4.3645321327653946D-02  /
   data green(20,11, 3) /   4.3438339077072841D-02  /
   data green(20,11, 4) /   4.3153492544407347D-02  /
   data green(20,11, 5) /   4.2795431401488532D-02  /
   data green(20,11, 6) /   4.2369764250962734D-02  /
   data green(20,11, 7) /   4.1882827323514428D-02  /
   data green(20,11, 8) /   4.1341444292623111D-02  /
   data green(20,11, 9) /   4.0752693086621265D-02  /
   data green(20,11,10) /   4.0123692517808576D-02  /
   data green(20,11,11) /   3.9461417746702199D-02  /
   data green(20,12, 0) /   4.2875181020716933D-02  /
   data green(20,12, 1) /   4.2835715832704756D-02  /
   data green(20,12, 2) /   4.2717976660746289D-02  /
   data green(20,12, 3) /   4.2523902167074891D-02  /
   data green(20,12, 4) /   4.2256624546427743D-02  /
   data green(20,12, 5) /   4.1920330546979030D-02  /
   data green(20,12, 6) /   4.1520084217775456D-02  /
   data green(20,12, 7) /   4.1061625650275560D-02  /
   data green(20,12, 8) /   4.0551160451938670D-02  /
   data green(20,12, 9) /   3.9995153520705588D-02  /
   data green(20,12,10) /   3.9400138308051523D-02  /
   data green(20,12,11) /   3.8772549698694411D-02  /
   data green(20,12,12) /   3.8118585413479497D-02  /
   data green(20,13, 0) /   4.1921383779830929D-02  /
   data green(20,13, 1) /   4.1884501144415137D-02  /
   data green(20,13, 2) /   4.1774439333874085D-02  /
   data green(20,13, 3) /   4.1592930619797928D-02  /
   data green(20,13, 4) /   4.1342777684119043D-02  /
   data green(20,13, 5) /   4.1027734612722473D-02  /
   data green(20,13, 6) /   4.0652354457555252D-02  /
   data green(20,13, 7) /   4.0221815153507964D-02  /
   data green(20,13, 8) /   3.9741736098267713D-02  /
   data green(20,13, 9) /   3.9217996894771422D-02  /
   data green(20,13,10) /   3.8656567931394957D-02  /
   data green(20,13,11) /   3.8063360043532592D-02  /
   data green(20,13,12) /   3.7444097871108993D-02  /
   data green(20,13,13) /   3.6804219039773454D-02  /
   data green(20,14, 0) /   4.0959835145352744D-02  /
   data green(20,14, 1) /   4.0925438035264841D-02  /
   data green(20,14, 2) /   4.0822768250733520D-02  /
   data green(20,14, 3) /   4.0653368351468877D-02  /
   data green(20,14, 4) /   4.0419737523431656D-02  /
   data green(20,14, 5) /   4.0125230157912417D-02  /
   data green(20,14, 6) /   3.9773925389476633D-02  /
   data green(20,14, 7) /   3.9370477268740878D-02  /
   data green(20,14, 8) /   3.8919955781491190D-02  /
   data green(20,14, 9) /   3.8427688390777780D-02  /
   data green(20,14,10) /   3.7899110401440665D-02  /
   data green(20,14,11) /   3.7339630536538028D-02  /
   data green(20,14,12) /   3.6754515992740287D-02  /
   data green(20,14,13) /   3.6148799182956791D-02  /
   data green(20,14,14) /   3.5527206579047815D-02  /
   data green(20,15, 0) /   3.9997589923775272D-02  /
   data green(20,15, 1) /   3.9965564619744966D-02  /
   data green(20,15, 2) /   3.9869951547301101D-02  /
   data green(20,15, 3) /   3.9712120559054855D-02  /
   data green(20,15, 4) /   3.9494294006670344D-02  /
   data green(20,15, 5) /   3.9219460637647335D-02  /
   data green(20,15, 6) /   3.8891264382820742D-02  /
   data green(20,15, 7) /   3.8513875937536204D-02  /
   data green(20,15, 8) /   3.8091855561831506D-02  /
   data green(20,15, 9) /   3.7630015192469897D-02  /
   data green(20,15,10) /   3.7133286935458705D-02  /
   data green(20,15,11) /   3.6606603523368024D-02  /
   data green(20,15,12) /   3.6054794625263946D-02  /
   data green(20,15,13) /   3.5482501209642969D-02  /
   data green(20,15,14) /   3.4894108650056481D-02  /
   data green(20,15,15) /   3.4293698033141501D-02  /
   data green(20,16, 0) /   3.9040640543445741D-02  /
   data green(20,16, 1) /   3.9010862388151095D-02  /
   data green(20,16, 2) /   3.8921937795556423D-02  /
   data green(20,16, 3) /   3.8775080648313705D-02  /
   data green(20,16, 4) /   3.8572262818152125D-02  /
   data green(20,16, 5) /   3.8316141286362369D-02  /
   data green(20,16, 6) /   3.8009963639197633D-02  /
   data green(20,16, 7) /   3.7657458366619975D-02  /
   data green(20,16, 8) /   3.7262716884490850D-02  /
   data green(20,16, 9) /   3.6830074013488874D-02  /
   data green(20,16,10) /   3.6363992897787023D-02  /
   data green(20,16,11) /   3.5868959204453564D-02  /
   data green(20,16,12) /   3.5349388100583454D-02  /
   data green(20,16,13) /   3.4809546134640150D-02  /
   data green(20,16,14) /   3.4253488889938920D-02  /
   data green(20,16,15) /   3.3685014223416583D-02  /
   data green(20,16,16) /   3.3107630097240812D-02  /
   data green(20,17, 0) /   3.8093982461721847D-02  /
   data green(20,17, 1) /   3.8066320653514243D-02  /
   data green(20,17, 2) /   3.7983697628600349D-02  /
   data green(20,17, 3) /   3.7847187346230250D-02  /
   data green(20,17, 4) /   3.7658536558863372D-02  /
   data green(20,17, 5) /   3.7420103264991456D-02  /
   data green(20,17, 6) /   3.7134776596839572D-02  /
   data green(20,17, 7) /   3.6805883355347881D-02  /
   data green(20,17, 8) /   3.6437086855567992D-02  /
   data green(20,17, 9) /   3.6032283660550786D-02  /
   data green(20,17,10) /   3.5595503241001157D-02  /
   data green(20,17,11) /   3.5130814727434651D-02  /
   data green(20,17,12) /   3.4642243865511971D-02  /
   data green(20,17,13) /   3.4133702181156800D-02  /
   data green(20,17,14) /   3.3608929322943390D-02  /
   data green(20,17,15) /   3.3071448654263279D-02  /
   data green(20,17,16) /   3.2524535460937420D-02  /
   data green(20,17,17) /   3.1971196634300009D-02  /
   data green(20,18, 0) /   3.7161698765041280D-02  /
   data green(20,18, 1) /   3.7136020225083589D-02  /
   data green(20,18, 2) /   3.7059304683640278D-02  /
   data green(20,18, 3) /   3.6932501253135429D-02  /
   data green(20,18, 4) /   3.6757155472885272D-02  /
   data green(20,18, 5) /   3.6535357415479607D-02  /
   data green(20,18, 6) /   3.6269673895205766D-02  /
   data green(20,18, 7) /   3.5963068994958720D-02  /
   data green(20,18, 8) /   3.5618817533045459D-02  /
   data green(20,18, 9) /   3.5240416074861149D-02  /
   data green(20,18,10) /   3.4831495711488164D-02  /
   data green(20,18,11) /   3.4395740169605421D-02  /
   data green(20,18,12) /   3.3936811993355465D-02  /
   data green(20,18,13) /   3.3458288655895764D-02  /
   data green(20,18,14) /   3.2963609607281848D-02  /
   data green(20,18,15) /   3.2456034513142676D-02  /
   data green(20,18,16) /   3.1938612326029577D-02  /
   data green(20,18,17) /   3.1414160374865545D-02  /
   data green(20,18,18) /   3.0885252354753312D-02  /
   data green(20,19, 0) /   3.6247054795614408D-02  /
   data green(20,19, 1) /   3.6223227188293926D-02  /
   data green(20,19, 2) /   3.6152026881064656D-02  /
   data green(20,19, 3) /   3.6034292081697866D-02  /
   data green(20,19, 4) /   3.5871389289409673D-02  /
   data green(20,19, 5) /   3.5665169585342847D-02  /
   data green(20,19, 6) /   3.5417911335274777D-02  /
   data green(20,19, 7) /   3.5132252709790757D-02  /
   data green(20,19, 8) /   3.4811117785346492D-02  /
   data green(20,19, 9) /   3.4457640017434399D-02  /
   data green(20,19,10) /   3.4075086611539157D-02  /
   data green(20,19,11) /   3.3666786825034115D-02  /
   data green(20,19,12) /   3.3236066594946156D-02  /
   data green(20,19,13) /   3.2786191185053694D-02  /
   data green(20,19,14) /   3.2320316853688469D-02  /
   data green(20,19,15) /   3.1841451916759828D-02  /
   data green(20,19,16) /   3.1352427055375165D-02  /
   data green(20,19,17) /   3.0855874311771943D-02  /
   data green(20,19,18) /   3.0354213933276378D-02  /
   data green(20,19,19) /   2.9849648052348330D-02  /
   data green(20,20, 0) /   3.5352595837726264D-02  /
   data green(20,20, 1) /   3.5330489848635724D-02  /
   data green(20,20, 2) /   3.5264421170974856D-02  /
   data green(20,20, 3) /   3.5155129839317879D-02  /
   data green(20,20, 4) /   3.5003823615196196D-02  /
   data green(20,20, 5) /   3.4812141189764167D-02  /
   data green(20,20, 6) /   3.4582103789288152D-02  /
   data green(20,20, 7) /   3.4316057930865146D-02  /
   data green(20,20, 8) /   3.4016612388696729D-02  /
   data green(20,20, 9) /   3.3686572485710109D-02  /
   data green(20,20,10) /   3.3328874645925599D-02  /
   data green(20,20,11) /   3.2946523777306454D-02  /
   data green(20,20,12) /   3.2542535563384677D-02  /
   data green(20,20,13) /   3.2119885187902195D-02  /
   data green(20,20,14) /   3.1681463457453152D-02  /
   data green(20,20,15) /   3.1230040768611415D-02  /
   data green(20,20,16) /   3.0768238920105655D-02  /
   data green(20,20,17) /   3.0298510414872248D-02  /
   data green(20,20,18) /   2.9823124636472555D-02  /
   data green(20,20,19) /   2.9344160115006009D-02  /
   data green(20,20,20) /   2.8863502008443166D-02  /
   data green(21, 0, 0) /   4.7646221003042076D-02  /
   data green(21, 1, 0) /   4.7591920643971525D-02  /
   data green(21, 1, 1) /   4.7537807669189212D-02  /
   data green(21, 2, 0) /   4.7430143901425498D-02  /
   data green(21, 2, 1) /   4.7376586520535668D-02  /
   data green(21, 2, 2) /   4.7217012804342498D-02  /
   data green(21, 3, 0) /   4.7164198183119344D-02  /
   data green(21, 3, 1) /   4.7111545475023178D-02  /
   data green(21, 3, 2) /   4.6954654584875632D-02  /
   data green(21, 3, 3) /   4.6696666323879497D-02  /
   data green(21, 4, 0) /   4.6799384829960350D-02  /
   data green(21, 4, 1) /   4.6747955725804570D-02  /
   data green(21, 4, 2) /   4.6594694077627057D-02  /
   data green(21, 4, 3) /   4.6342619340609972D-02  /
   data green(21, 4, 4) /   4.5996578583657792D-02  /
   data green(21, 5, 0) /   4.6342707867378194D-02  /
   data green(21, 5, 1) /   4.6292782465689894D-02  /
   data green(21, 5, 2) /   4.6143981814486382D-02  /
   data green(21, 5, 3) /   4.5899178998587184D-02  /
   data green(21, 5, 4) /   4.5562989763800410D-02  /
   data green(21, 5, 5) /   4.5141532853018841D-02  /
   data green(21, 6, 0) /   4.5802511792821826D-02  /
   data green(21, 6, 1) /   4.5754325461345474D-02  /
   data green(21, 6, 2) /   4.5610685409870219D-02  /
   data green(21, 6, 3) /   4.5374299256133323D-02  /
   data green(21, 6, 4) /   4.5049520901697365D-02  /
   data green(21, 6, 5) /   4.4642129747486749D-02  /
   data green(21, 6, 6) /   4.4159053281602675D-02  /
   data green(21, 7, 0) /   4.5188083794584430D-02  /
   data green(21, 7, 1) /   4.5141824167979981D-02  /
   data green(21, 7, 2) /   4.5003903198849754D-02  /
   data green(21, 7, 3) /   4.4776849982241421D-02  /
   data green(21, 7, 4) /   4.4464735500830549D-02  /
   data green(21, 7, 5) /   4.4072971708905194D-02  /
   data green(21, 7, 6) /   4.3608058269614579D-02  /
   data green(21, 7, 7) /   4.3077299422081945D-02  /
   data green(21, 8, 0) /   4.4509254613290840D-02  /
   data green(21, 8, 1) /   4.4465061291202639D-02  /
   data green(21, 8, 2) /   4.4333275759784962D-02  /
   data green(21, 8, 3) /   4.4116241294088296D-02  /
   data green(21, 8, 4) /   4.3817733923358922D-02  /
   data green(21, 8, 5) /   4.3442781623784350D-02  /
   data green(21, 8, 6) /   4.2997435613642151D-02  /
   data green(21, 8, 7) /   4.2488513499054810D-02  /
   data green(21, 8, 8) /   4.1923334236565295D-02  /
   data green(21, 9, 0) /   4.3776026812174407D-02  /
   data green(21, 9, 1) /   4.3733993339933108D-02  /
   data green(21, 9, 2) /   4.3608623210392244D-02  /
   data green(21, 9, 3) /   4.3402071751771756D-02  /
   data green(21, 9, 4) /   4.3117816160222484D-02  /
   data green(21, 9, 5) /   4.2760494413358212D-02  /
   data green(21, 9, 6) /   4.2335700769975067D-02  /
   data green(21, 9, 7) /   4.1849754991136891D-02  /
   data green(21, 9, 8) /   4.1309462770090143D-02  /
   data green(21, 9, 9) /   4.0721883211984722D-02  /
   data green(21,10, 0) /   4.2998251285167426D-02  /
   data green(21,10, 1) /   4.2958428875163136D-02  /
   data green(21,10, 2) /   4.2839628613058560D-02  /
   data green(21,10, 3) /   4.2643820143578669D-02  /
   data green(21,10, 4) /   4.2374184911928264D-02  /
   data green(21,10, 5) /   4.2034973921338786D-02  /
   data green(21,10, 6) /   4.1631326480309976D-02  /
   data green(21,10, 7) /   4.1169064642059737D-02  /
   data green(21,10, 8) /   4.0654478485019040D-02  /
   data green(21,10, 9) /   4.0094116135492332D-02  /
   data green(21,10,10) /   3.9494589943817146D-02  /
   data green(21,11, 0) /   4.2185364199899218D-02  /
   data green(21,11, 1) /   4.2147766618844291D-02  /
   data green(21,11, 2) /   4.2035579552469497D-02  /
   data green(21,11, 3) /   4.1850592748003554D-02  /
   data green(21,11, 4) /   4.1595700604622372D-02  /
   data green(21,11, 5) /   4.1274777550121357D-02  /
   data green(21,11, 6) /   4.0892518637843024D-02  /
   data green(21,11, 7) /   4.0454257855006739D-02  /
   data green(21,11, 8) /   3.9965777136880716D-02  /
   data green(21,11, 9) /   3.9433118166204775D-02  /
   data green(21,11,10) /   3.8862407051836564D-02  /
   data green(21,11,11) /   3.8259699369862267D-02  /
   data green(21,12, 0) /   4.1346188668101178D-02  /
   data green(21,12, 1) /   4.1310797753435546D-02  /
   data green(21,12, 2) /   4.1205172323434056D-02  /
   data green(21,12, 3) /   4.1030930668138732D-02  /
   data green(21,12, 4) /   4.0790693117751693D-02  /
   data green(21,12, 5) /   4.0487973594423791D-02  /
   data green(21,12, 6) /   4.0127040372599440D-02  /
   data green(21,12, 7) /   3.9712756569118608D-02  /
   data green(21,12, 8) /   3.9250411412916048D-02  /
   data green(21,12, 9) /   3.8745552698216001D-02  /
   data green(21,12,10) /   3.8203829267840314D-02  /
   data green(21,12,11) /   3.7630850254224552D-02  /
   data green(21,12,12) /   3.7032065480979384D-02  /
   data green(21,13, 0) /   4.0488799130992403D-02  /
   data green(21,13, 1) /   4.0455570482837953D-02  /
   data green(21,13, 2) /   4.0356377040688014D-02  /
   data green(21,13, 3) /   4.0192675925273494D-02  /
   data green(21,13, 4) /   3.9966829383476926D-02  /
   data green(21,13, 5) /   3.9682010961755353D-02  /
   data green(21,13, 6) /   3.9342084587232075D-02  /
   data green(21,13, 7) /   3.8951465352691866D-02  /
   data green(21,13, 8) /   3.8514971330526386D-02  /
   data green(21,13, 9) /   3.8037675306103490D-02  /
   data green(21,13,10) /   3.7524764119046684D-02  /
   data green(21,13,11) /   3.6981411602377535D-02  /
   data green(21,13,12) /   3.6412669198763695D-02  /
   data green(21,13,13) /   3.5823376458982654D-02  /
   data green(21,14, 0) /   3.9620442080378182D-02  /
   data green(21,14, 1) /   3.9589310579477469D-02  /
   data green(21,14, 2) /   3.9496357699264503D-02  /
   data green(21,14, 3) /   3.9342890820223096D-02  /
   data green(21,14, 4) /   3.9131031986316865D-02  /
   data green(21,14, 5) /   3.8863637124186452D-02  /
   data green(21,14, 6) /   3.8544191552961113D-02  /
   data green(21,14, 7) /   3.8176689092378792D-02  /
   data green(21,14, 8) /   3.7765502587209329D-02  /
   data green(21,14, 9) /   3.7315253392744788D-02  /
   data green(21,14,10) /   3.6830686452996204D-02  /
   data green(21,14,11) /   3.6316556257526257D-02  /
   data green(21,14,12) /   3.5777527408968168D-02  /
   data green(21,14,13) /   3.5218091973888317D-02  /
   data green(21,14,14) /   3.4642504380629061D-02  /
   data green(21,15, 0) /   3.8747504239242822D-02  /
   data green(21,15, 1) /   3.8718389147245451D-02  /
   data green(21,15, 2) /   3.8631438725906249D-02  /
   data green(21,15, 3) /   3.8487822598183376D-02  /
   data green(21,15, 4) /   3.8289441447347569D-02  /
   data green(21,15, 5) /   3.8038857748358294D-02  /
   data green(21,15, 6) /   3.7739205846346581D-02  /
   data green(21,15, 7) /   3.7394087421187418D-02  /
   data green(21,15, 8) /   3.7007458856166131D-02  /
   data green(21,15, 9) /   3.6583516875104603D-02  /
   data green(21,15,10) /   3.6126588129740136D-02  /
   data green(21,15,11) /   3.5641027364523216D-02  /
   data green(21,15,12) /   3.5131127534519135D-02  /
   data green(21,15,13) /   3.4601043966925031D-02  /
   data green(21,15,14) /   3.4054733468210488D-02  /
   data green(21,15,15) /   3.3495908274857392D-02  /
   data green(21,16, 0) /   3.7875518385653582D-02  /
   data green(21,16, 1) /   3.7848327879427740D-02  /
   data green(21,16, 2) /   3.7767108581297573D-02  /
   data green(21,16, 3) /   3.7632904427428240D-02  /
   data green(21,16, 4) /   3.7447413799471502D-02  /
   data green(21,16, 5) /   3.7212930321182260D-02  /
   data green(21,16, 6) /   3.6932265735032742D-02  /
   data green(21,16, 7) /   3.6608659827705302D-02  /
   data green(21,16, 8) /   3.6245682813303086D-02  /
   data green(21,16, 9) /   3.5847135515367251D-02  /
   data green(21,16,10) /   3.5416952187123318D-02  /
   data green(21,16,11) /   3.4959109991315157D-02  /
   data green(21,16,12) /   3.4477548161969934D-02  /
   data green(21,16,13) /   3.3976098820454936D-02  /
   data green(21,16,14) /   3.3458430424842110D-02  /
   data green(21,16,15) /   3.2928003971157915D-02  /
   data green(21,16,16) /   3.2388041380818321D-02  /
   data green(21,17, 0) /   3.7009197207757596D-02  /
   data green(21,17, 1) /   3.6983832277547667D-02  /
   data green(21,17, 2) /   3.6908051109908684D-02  /
   data green(21,17, 3) /   3.6782783761775982D-02  /
   data green(21,17, 4) /   3.6609545022361925D-02  /
   data green(21,17, 5) /   3.6390383942152722D-02  /
   data green(21,17, 6) /   3.6127817863739892D-02  /
   data green(21,17, 7) /   3.5824755028661066D-02  /
   data green(21,17, 8) /   3.5484410231488718D-02  /
   data green(21,17, 9) /   3.5110217984117217D-02  /
   data green(21,17,10) /   3.4705747291255716D-02  /
   data green(21,17,11) /   3.4274621509753744D-02  /
   data green(21,17,12) /   3.3820445973371437D-02  /
   data green(21,17,13) /   3.3346745213536380D-02  /
   data green(21,17,14) /   3.2856910783299227D-02  /
   data green(21,17,15) /   3.2354159961231775D-02  /
   data green(21,17,16) /   3.1841505014417758D-02  /
   data green(21,17,17) /   3.1321732251760194D-02  /
   data green(21,18, 0) /   3.6152486515756560D-02  /
   data green(21,18, 1) /   3.6128844215994285D-02  /
   data green(21,18, 2) /   3.6058196197236443D-02  /
   data green(21,18, 3) /   3.5941369931327849D-02  /
   data green(21,18, 4) /   3.5779714566070549D-02  /
   data green(21,18, 5) /   3.5575057981865464D-02  /
   data green(21,18, 6) /   3.5329650480137927D-02  /
   data green(21,18, 7) /   3.5046098435093605D-02  /
   data green(21,18, 8) /   3.4727291593645813D-02  /
   data green(21,18, 9) /   3.4376327739370126D-02  /
   data green(21,18,10) /   3.3996438180704046D-02  /
   data green(21,18,11) /   3.3590917045492746D-02  /
   data green(21,18,12) /   3.3163056742327633D-02  /
   data green(21,18,13) /   3.2716091264205806D-02  /
   data green(21,18,14) /   3.2253148332757303D-02  /
   data green(21,18,15) /   3.1777210766425125D-02  /
   data green(21,18,16) /   3.1291086939504505D-02  /
   data green(21,18,17) /   3.0797389798653266D-02  /
   data green(21,18,18) /   3.0298523621901489D-02  /
   data green(21,19, 0) /   3.5308630461736351D-02  /
   data green(21,19, 1) /   3.5286606544183557D-02  /
   data green(21,19, 2) /   3.5220782547908489D-02  /
   data green(21,19, 3) /   3.5111893974778249D-02  /
   data green(21,19, 4) /   3.4961141249778789D-02  /
   data green(21,19, 5) /   3.4770153232792998D-02  /
   data green(21,19, 6) /   3.4540939223422142D-02  /
   data green(21,19, 7) /   3.4275832177679218D-02  /
   data green(21,19, 8) /   3.3977426166108902D-02  /
   data green(21,19, 9) /   3.3648511158231885D-02  /
   data green(21,19,10) /   3.3292008042274786D-02  /
   data green(21,19,11) /   3.2910906428772584D-02  /
   data green(21,19,12) /   3.2508207301439634D-02  /
   data green(21,19,13) /   3.2086872031053239D-02  /
   data green(21,19,14) /   3.1649778714723145D-02  /
   data green(21,19,15) /   3.1199686289327663D-02  /
   data green(21,19,16) /   3.0739206425673904D-02  /
   data green(21,19,17) /   3.0270782856650857D-02  /
   data green(21,19,18) /   2.9796677533587033D-02  /
   data green(21,19,19) /   2.9318962835988595D-02  /
   data green(21,20, 0) /   3.4480242868928632D-02  /
   data green(21,20, 1) /   3.4459733853262661D-02  /
   data green(21,20, 2) /   3.4398426787883088D-02  /
   data green(21,20, 3) /   3.4296975043170246D-02  /
   data green(21,20, 4) /   3.4156446041637853D-02  /
   data green(21,20, 5) /   3.3978290281707457D-02  /
   data green(21,20, 6) /   3.3764300472274994D-02  /
   data green(21,20, 7) /   3.3516562993812261D-02  /
   data green(21,20, 8) /   3.3237404172456027D-02  /
   data green(21,20, 9) /   3.2929333922384785D-02  /
   data green(21,20,10) /   3.2594989194729335D-02  /
   data green(21,20,11) /   3.2237079402230091D-02  /
   data green(21,20,12) /   3.1858335612603525D-02  /
   data green(21,20,13) /   3.1461464868233993D-02  /
   data green(21,20,14) /   3.1049110540414163D-02  /
   data green(21,20,15) /   3.0623819200295010D-02  /
   data green(21,20,16) /   3.0188014113199959D-02  /
   data green(21,20,17) /   2.9743975154555621D-02  /
   data green(21,20,18) /   2.9293824711202159D-02  /
   data green(21,20,19) /   2.8839518970211812D-02  /
   data green(21,20,20) /   2.8382843901880294D-02  /
   data green(21,21, 0) /   3.3669380173015184D-02  /
   data green(21,21, 1) /   3.3650284923269387D-02  /
   data green(21,21, 2) /   3.3593194446981299D-02  /
   data green(21,21, 3) /   3.3498688999806847D-02  /
   data green(21,21, 4) /   3.3367717442461621D-02  /
   data green(21,21, 5) /   3.3201570955338336D-02  /
   data green(21,21, 6) /   3.3001848263768452D-02  /
   data green(21,21, 7) /   3.2770414178193451D-02  /
   data green(21,21, 8) /   3.2509353487537258D-02  /
   data green(21,21, 9) /   3.2220922318836567D-02  /
   data green(21,21,10) /   3.1907499001962943D-02  /
   data green(21,21,11) /   3.1571536279450613D-02  /
   data green(21,21,12) /   3.1215516411443937D-02  /
   data green(21,21,13) /   3.0841910381529190D-02  /
   data green(21,21,14) /   3.0453142046237190D-02  /
   data green(21,21,15) /   3.0051557719707231D-02  /
   data green(21,21,16) /   2.9639401368498443D-02  /
   data green(21,21,17) /   2.9218795324887824D-02  /
   data green(21,21,18) /   2.8791726217648833D-02  /
   data green(21,21,19) /   2.8360035668281040D-02  /
   data green(21,21,20) /   2.7925415204162161D-02  /
   data green(21,21,21) /   2.7489404791148895D-02  /
   data green(22, 0, 0) /   4.5478165166064635D-02  /
   data green(22, 1, 0) /   4.5430962243515796D-02  /
   data green(22, 1, 1) /   4.5383907546982093D-02  /
   data green(22, 2, 0) /   4.5290242850719462D-02  /
   data green(22, 2, 1) /   4.5243628117243084D-02  /
   data green(22, 2, 2) /   4.5104654645399062D-02  /
   data green(22, 3, 0) /   4.5058627958511836D-02  /
   data green(22, 3, 1) /   4.5012731136341161D-02  /
   data green(22, 3, 2) /   4.4875888870684708D-02  /
   data green(22, 3, 3) /   4.4650601691221259D-02  /
   data green(22, 4, 0) /   4.4740333660005349D-02  /
   data green(22, 4, 1) /   4.4695410853761747D-02  /
   data green(22, 4, 2) /   4.4561460477719153D-02  /
   data green(22, 4, 3) /   4.4340894801059801D-02  /
   data green(22, 4, 4) /   4.4037599007435965D-02  /
   data green(22, 5, 0) /   4.4340959447060846D-02  /
   data green(22, 5, 1) /   4.4297238407196383D-02  /
   data green(22, 5, 2) /   4.4166856765280767D-02  /
   data green(22, 5, 3) /   4.3952119730843647D-02  /
   data green(22, 5, 4) /   4.3656742490801856D-02  /
   data green(22, 5, 5) /   4.3285673013508871D-02  /
   data green(22, 6, 0) /   4.3867221958747878D-02  /
   data green(22, 6, 1) /   4.3824897444861344D-02  /
   data green(22, 6, 2) /   4.3698663786837275D-02  /
   data green(22, 6, 3) /   4.3490704343594887D-02  /
   data green(22, 6, 4) /   4.3204540610767457D-02  /
   data green(22, 6, 5) /   4.2844867867647518D-02  /
   data green(22, 6, 6) /   4.2417346741428559D-02  /
   data green(22, 7, 0) /   4.3326658243418678D-02  /
   data green(22, 7, 1) /   4.3285889375732807D-02  /
   data green(22, 7, 2) /   4.3164277461481403D-02  /
   data green(22, 7, 3) /   4.2963873378703991D-02  /
   data green(22, 7, 4) /   4.2687987776880744D-02  /
   data green(22, 7, 5) /   4.2341040370302140D-02  /
   data green(22, 7, 6) /   4.1928368293785163D-02  /
   data green(22, 7, 7) /   4.1456009341734941D-02  /
   data green(22, 8, 0) /   4.2727322455875495D-02  /
   data green(22, 8, 1) /   4.2688231971300158D-02  /
   data green(22, 8, 2) /   4.2571607797236267D-02  /
   data green(22, 8, 3) /   4.2379361746193250D-02  /
   data green(22, 8, 4) /   4.2114582846064154D-02  /
   data green(22, 8, 5) /   4.1781400591313668D-02  /
   data green(22, 8, 6) /   4.1384810541433441D-02  /
   data green(22, 8, 7) /   4.0930476299834165D-02  /
   data green(22, 8, 8) /   4.0424522358783160D-02  /
   data green(22, 9, 0) /   4.2077496931215361D-02  /
   data green(22, 9, 1) /   4.2040172106618728D-02  /
   data green(22, 9, 2) /   4.1928796561053451D-02  /
   data green(22, 9, 3) /   4.1745140175578198D-02  /
   data green(22, 9, 4) /   4.1492065467121440D-02  /
   data green(22, 9, 5) /   4.1173404682892878D-02  /
   data green(22, 9, 6) /   4.0793802571812407D-02  /
   data green(22, 9, 7) /   4.0358537132991701D-02  /
   data green(22, 9, 8) /   3.9873331137778502D-02  /
   data green(22, 9, 9) /   3.9344166323339348D-02  /
   data green(22,10, 0) /   4.1385433731411854D-02  /
   data green(22,10, 1) /   4.1349928639842978D-02  /
   data green(22,10, 2) /   4.1243964107553530D-02  /
   data green(22,10, 3) /   4.1069168464396837D-02  /
   data green(22,10, 4) /   4.0828177987521974D-02  /
   data green(22,10, 5) /   4.0524527379757216D-02  /
   data green(22,10, 6) /   4.0162509216581502D-02  /
   data green(22,10, 7) /   3.9747013027504725D-02  /
   data green(22,10, 8) /   3.9283355197924238D-02  /
   data green(22,10, 9) /   3.8777110204893675D-02  /
   data green(22,10,10) /   3.8233952106282415D-02  /
   data green(22,11, 0) /   4.0659137142472580D-02  /
   data green(22,11, 1) /   4.0625475868035020D-02  /
   data green(22,11, 2) /   4.0524995683822444D-02  /
   data green(22,11, 3) /   4.0359186442139233D-02  /
   data green(22,11, 4) /   4.0130462755522005D-02  /
   data green(22,11, 5) /   3.9842067141011621D-02  /
   data green(22,11, 6) /   3.9497945318057961D-02  /
   data green(22,11, 7) /   3.9102602820623923D-02  /
   data green(22,11, 8) /   3.8660952609514555D-02  /
   data green(22,11, 9) /   3.8178162888181731D-02  /
   data green(22,11,10) /   3.7659513045486301D-02  /
   data green(22,11,11) /   3.7110263860349157D-02  /
   data green(22,12, 0) /   3.9906192150127348D-02  /
   data green(22,12, 1) /   3.9874372594690989D-02  /
   data green(22,12, 2) /   3.9779372265126088D-02  /
   data green(22,12, 3) /   3.9622547727294684D-02  /
   data green(22,12, 4) /   3.9406099909420539D-02  /
   data green(22,12, 5) /   3.9132989017335890D-02  /
   data green(22,12, 6) /   3.8806824629203630D-02  /
   data green(22,12, 7) /   3.8431738771221850D-02  /
   data green(22,12, 8) /   3.8012250290815931D-02  /
   data green(22,12, 9) /   3.7553128515525655D-02  /
   data green(22,12,10) /   3.7059263175666532D-02  /
   data green(22,12,11) /   3.6535546105593446D-02  /
   data green(22,12,12) /   3.5986768566420281D-02  /
   data green(22,13, 0) /   3.9133639291968281D-02  /
   data green(22,13, 1) /   3.9103637248018956D-02  /
   data green(22,13, 2) /   3.9014046479443981D-02  /
   data green(22,13, 3) /   3.8866097019316391D-02  /
   data green(22,13, 4) /   3.8661786625644956D-02  /
   data green(22,13, 5) /   3.8403806468490474D-02  /
   data green(22,13, 6) /   3.8095444825829390D-02  /
   data green(22,13, 7) /   3.7740475383387359D-02  /
   data green(22,13, 8) /   3.7343037225969808D-02  /
   data green(22,13, 9) /   3.6907513401315553D-02  /
   data green(22,13,10) /   3.6438414153567467D-02  /
   data green(22,13,11) /   3.5940269741171357D-02  /
   data green(22,13,12) /   3.5417536370165349D-02  /
   data green(22,13,13) /   3.4874517368818497D-02  /
   data green(22,14, 0) /   3.8347892807726137D-02  /
   data green(22,14, 1) /   3.8319666031873957D-02  /
   data green(22,14, 2) /   3.8235360778558679D-02  /
   data green(22,14, 3) /   3.8096088360831455D-02  /
   data green(22,14, 4) /   3.7903655619235449D-02  /
   data green(22,14, 5) /   3.7660500337814588D-02  /
   data green(22,14, 6) /   3.7369607278146179D-02  /
   data green(22,14, 7) /   3.7034410371482016D-02  /
   data green(22,14, 8) /   3.6658687069948626D-02  /
   data green(22,14, 9) /   3.6246450744277320D-02  /
   data green(22,14,10) /   3.5801846417478703D-02  /
   data green(22,14,11) /   3.5329054179754349D-02  /
   data green(22,14,12) /   3.4832203496371235D-02  /
   data green(22,14,13) /   3.4315300444868149D-02  /
   data green(22,14,14) /   3.3782168820528903D-02  /
   data green(22,15, 0) /   3.7554696751158884D-02  /
   data green(22,15, 1) /   3.7528188839673397D-02  /
   data green(22,15, 2) /   3.7449002782145448D-02  /
   data green(22,15, 3) /   3.7318139609140141D-02  /
   data green(22,15, 4) /   3.7137228545084888D-02  /
   data green(22,15, 5) /   3.6908471113439534D-02  /
   data green(22,15, 6) /   3.6634568237719087D-02  /
   data green(22,15, 7) /   3.6318634966991427D-02  /
   data green(22,15, 8) /   3.5964107878360831D-02  /
   data green(22,15, 9) /   3.5574650164139424D-02  /
   data green(22,15,10) /   3.5154058963234758D-02  /
   data green(22,15,11) /   3.4706178750670558D-02  /
   data green(22,15,12) /   3.4234823679464496D-02  /
   data green(22,15,13) /   3.3743710794834299D-02  /
   data green(22,15,14) /   3.3236405111952867D-02  /
   data green(22,15,15) /   3.2716276737265973D-02  /
   data green(22,16, 0) /   3.6759112564619757D-02  /
   data green(22,16, 1) /   3.6734256496671079D-02  /
   data green(22,16, 2) /   3.6659991549557978D-02  /
   data green(22,16, 3) /   3.6537217175527367D-02  /
   data green(22,16, 4) /   3.6367398760356814D-02  /
   data green(22,16, 5) /   3.6152519415715886D-02  /
   data green(22,16, 6) /   3.5895016912024715D-02  /
   data green(22,16, 7) /   3.5597709601427496D-02  /
   data green(22,16, 8) /   3.5263715564600111D-02  /
   data green(22,16, 9) /   3.4896369219066282D-02  /
   data green(22,16,10) /   3.4499139297274592D-02  /
   data green(22,16,11) /   3.4075551520067006D-02  /
   data green(22,16,12) /   3.3629118551575465D-02  /
   data green(22,16,13) /   3.3163279020815922D-02  /
   data green(22,16,14) /   3.2681346615991237D-02  /
   data green(22,16,15) /   3.2186469561781318D-02  /
   data green(22,16,16) /   3.1681600216517704D-02  /
   data green(22,17, 0) /   3.5965531328318277D-02  /
   data green(22,17, 1) /   3.5942252598683333D-02  /
   data green(22,17, 2) /   3.5872688212994081D-02  /
   data green(22,17, 3) /   3.5757644733398979D-02  /
   data green(22,17, 4) /   3.5598437505837313D-02  /
   data green(22,17, 5) /   3.5396849198775227D-02  /
   data green(22,17, 6) /   3.5155075399615252D-02  /
   data green(22,17, 7) /   3.4875660458590911D-02  /
   data green(22,17, 8) /   3.4561427113637098D-02  /
   data green(22,17, 9) /   3.4215403466670673D-02  /
   data green(22,17,10) /   3.3840750644834841D-02  /
   data green(22,17,11) /   3.3440694029652690D-02  /
   data green(22,17,12) /   3.3018460347128957D-02  /
   data green(22,17,13) /   3.2577222258730040D-02  /
   data green(22,17,14) /   3.2120051444484335D-02  /
   data green(22,17,15) /   3.1649880577827333D-02  /
   data green(22,17,16) /   3.1169474092254325D-02  /
   data green(22,17,17) /   3.0681407250119003D-02  /
   data green(22,18, 0) /   3.5177704228015423D-02  /
   data green(22,18, 1) /   3.5155923536157783D-02  /
   data green(22,18, 2) /   3.5090824694639020D-02  /
   data green(22,18, 3) /   3.4983129834970365D-02  /
   data green(22,18, 4) /   3.4834017738203944D-02  /
   data green(22,18, 5) /   3.4645088254558899D-02  /
   data green(22,18, 6) /   3.4418315481500684D-02  /
   data green(22,18, 7) /   3.4155992335478828D-02  /
   data green(22,18, 8) /   3.3860669456767953D-02  /
   data green(22,18, 9) /   3.3535091444660554D-02  /
   data green(22,18,10) /   3.3182133254342688D-02  /
   data green(22,18,11) /   3.2804739241793064D-02  /
   data green(22,18,12) /   3.2405866876100670D-02  /
   data green(22,18,13) /   3.1988436609664397D-02  /
   data green(22,18,14) /   3.1555288860608703D-02  /
   data green(22,18,15) /   3.1109148562586810D-02  /
   data green(22,18,16) /   3.0652597305827343D-02  /
   data green(22,18,17) /   3.0188052747236970D-02  /
   data green(22,18,18) /   2.9717754712175033D-02  /
   data green(22,19, 0) /   3.4398785501751267D-02  /
   data green(22,19, 1) /   3.4378420997905661D-02  /
   data green(22,19, 2) /   3.4317544902304944D-02  /
   data green(22,19, 3) /   3.4216802994786133D-02  /
   data green(22,19, 4) /   3.4077250393983391D-02  /
   data green(22,19, 5) /   3.3900321074641214D-02  /
   data green(22,19, 6) /   3.3687787642372909D-02  /
   data green(22,19, 7) /   3.3441713537847521D-02  /
   data green(22,19, 8) /   3.3164400109293897D-02  /
   data green(22,19, 9) /   3.2858331061818291D-02  /
   data green(22,19,10) /   3.2526116679671639D-02  /
   data green(22,19,11) /   3.2170439956164178D-02  /
   data green(22,19,12) /   3.1794006398934600D-02  /
   data green(22,19,13) /   3.1399498852663091D-02  /
   data green(22,19,14) /   3.0989538241097927D-02  /
   data green(22,19,15) /   3.0566650712088982D-02  /
   data green(22,19,16) /   3.0133241299985177D-02  /
   data green(22,19,17) /   2.9691573915663301D-02  /
   data green(22,19,18) /   2.9243757242502715D-02  /
   data green(22,19,19) /   2.8791735955985529D-02  /
   data green(22,20, 0) /   3.3631383034568682D-02  /
   data green(22,20, 1) /   3.3612352146675475D-02  /
   data green(22,20, 2) /   3.3555453665000921D-02  /
   data green(22,20, 3) /   3.3461264615597543D-02  /
   data green(22,20, 4) /   3.3330728617512802D-02  /
   data green(22,20, 5) /   3.3165129798436395D-02  /
   data green(22,20, 6) /   3.2966058281651403D-02  /
   data green(22,20, 7) /   3.2735369031404789D-02  /
   data green(22,20, 8) /   3.2475136076196141D-02  /
   data green(22,20, 9) /   3.2187604204344654D-02  /
   data green(22,20,10) /   3.1875140153597770D-02  /
   data green(22,20,11) /   3.1540185120518011D-02  /
   data green(22,20,12) /   3.1185210128888842D-02  /
   data green(22,20,13) /   3.0812675455901585D-02  /
   data green(22,20,14) /   3.0424994955542520D-02  /
   data green(22,20,15) /   3.0024505770509607D-02  /
   data green(22,20,16) /   2.9613443610071042D-02  /
   data green(22,20,17) /   2.9193923506520673D-02  /
   data green(22,20,18) /   2.8767925754737078D-02  /
   data green(22,20,19) /   2.8337286588912707D-02  /
   data green(22,20,20) /   2.7903693054090652D-02  /
   data green(22,21, 0) /   3.2877612726847506D-02  /
   data green(22,21, 1) /   3.2859833606547227D-02  /
   data green(22,21, 2) /   3.2806669591522668D-02  /
   data green(22,21, 3) /   3.2718636011528801D-02  /
   data green(22,21, 4) /   3.2596576308174405D-02  /
   data green(22,21, 5) /   3.2441639731192441D-02  /
   data green(22,21, 6) /   3.2255251752264430D-02  /
   data green(22,21, 7) /   3.2039078665115260D-02  /
   data green(22,21, 8) /   3.1794988041640750D-02  /
   data green(22,21, 9) /   3.1525006789374768D-02  /
   data green(22,21,10) /   3.1231278512015851D-02  /
   data green(22,21,11) /   3.0916021729418285D-02  /
   data green(22,21,12) /   3.0581490291238370D-02  /
   data green(22,21,13) /   3.0229937047512961D-02  /
   data green(22,21,14) /   2.9863581547679117D-02  /
   data green(22,21,15) /   2.9484582251458887D-02  /
   data green(22,21,16) /   2.9095013469997811D-02  /
   data green(22,21,17) /   2.8696847027034164D-02  /
   data green(22,21,18) /   2.8291938445159891D-02  /
   data green(22,21,19) /   2.7882017323725816D-02  /
   data green(22,21,20) /   2.7468681480891162D-02  /
   data green(22,21,21) /   2.7053394378135615D-02  /
   data green(22,22, 0) /   3.2139153671528850D-02  /
   data green(22,22, 1) /   3.2122546304849925D-02  /
   data green(22,22, 2) /   3.2072878917303524D-02  /
   data green(22,22, 3) /   3.1990611632156010D-02  /
   data green(22,22, 4) /   3.1876498145743268D-02  /
   data green(22,22, 5) /   3.1731566663610884D-02  /
   data green(22,22, 6) /   3.1557094550377703D-02  /
   data green(22,22, 7) /   3.1354577899530257D-02  /
   data green(22,22, 8) /   3.1125697402374006D-02  /
   data green(22,22, 9) /   3.0872281968549999D-02  /
   data green(22,22,10) /   3.0596271527686146D-02  /
   data green(22,22,11) /   3.0299680335394555D-02  /
   data green(22,22,12) /   2.9984561935551830D-02  /
   data green(22,22,13) /   2.9652976716303448D-02  /
   data green(22,22,14) /   2.9306962761459256D-02  /
   data green(22,22,15) /   2.8948510461753892D-02  /
   data green(22,22,16) /   2.8579541128011983D-02  /
   data green(22,22,17) /   2.8201889652307788D-02  /
   data green(22,22,18) /   2.7817291100876287D-02  /
   data green(22,22,19) /   2.7427370996825504D-02  /
   data green(22,22,20) /   2.7033638961220553D-02  /
   data green(22,22,21) /   2.6637485324984663D-02  /
   data green(22,22,22) /   2.6240180296811768D-02  /
   data green(23, 0, 0) /   4.3498921019282412D-02  /
   data green(23, 1, 0) /   4.3457629706374859D-02  /
   data green(23, 1, 1) /   4.3416456891185576D-02  /
   data green(23, 2, 0) /   4.3334466766377713D-02  /
   data green(23, 2, 1) /   4.3293646003481061D-02  /
   data green(23, 2, 2) /   4.3171881092196121D-02  /
   data green(23, 3, 0) /   4.3131530784544507D-02  /
   data green(23, 3, 1) /   4.3091285542822398D-02  /
   data green(23, 3, 2) /   4.2971230688552103D-02  /
   data green(23, 3, 3) /   4.2773376479958887D-02  /
   data green(23, 4, 0) /   4.2852208016975563D-02  /
   data green(23, 4, 1) /   4.2812745715503912D-02  /
   data green(23, 4, 2) /   4.2695017497900889D-02  /
   data green(23, 4, 3) /   4.2500968558528737D-02  /
   data green(23, 4, 4) /   4.2233740877234707D-02  /
   data green(23, 5, 0) /   4.2501016381840881D-02  /
   data green(23, 5, 1) /   4.2462523466523662D-02  /
   data green(23, 5, 2) /   4.2347676374601739D-02  /
   data green(23, 5, 3) /   4.2158340985737328D-02  /
   data green(23, 5, 4) /   4.1897532818454708D-02  /
   data green(23, 5, 5) /   4.1569284511135507D-02  /
   data green(23, 6, 0) /   4.2083407379463061D-02  /
   data green(23, 6, 1) /   4.2046045679308156D-02  /
   data green(23, 6, 2) /   4.1934561295935001D-02  /
   data green(23, 6, 3) /   4.1750729298283676D-02  /
   data green(23, 6, 4) /   4.1497420283982400D-02  /
   data green(23, 6, 5) /   4.1178476722391263D-02  /
   data green(23, 6, 6) /   4.0798554844103233D-02  /
   data green(23, 7, 0) /   4.1605542499838666D-02  /
   data green(23, 7, 1) /   4.1569446965355153D-02  /
   data green(23, 7, 2) /   4.1461727231528252D-02  /
   data green(23, 7, 3) /   4.1284058961832147D-02  /
   data green(23, 7, 4) /   4.1039153954616076D-02  /
   data green(23, 7, 5) /   4.0730645983280481D-02  /
   data green(23, 7, 6) /   4.0362944498197328D-02  /
   data green(23, 7, 7) /   3.9941067441088961D-02  /
   data green(23, 8, 0) /   4.1074061004727720D-02  /
   data green(23, 8, 1) /   4.1039338801135028D-02  /
   data green(23, 8, 2) /   4.0935703283604262D-02  /
   data green(23, 8, 3) /   4.0764724985488492D-02  /
   data green(23, 8, 4) /   4.0528947571516913D-02  /
   data green(23, 8, 5) /   4.0231783483404948D-02  /
   data green(23, 8, 6) /   3.9877379878213132D-02  /
   data green(23, 8, 7) /   3.9470464925006908D-02  /
   data green(23, 8, 8) /   3.9016185040908279D-02  /
   data green(23, 9, 0) /   4.0495854374592827D-02  /
   data green(23, 9, 1) /   4.0462585196367286D-02  /
   data green(23, 9, 2) /   4.0363271968941794D-02  /
   data green(23, 9, 3) /   4.0199377062932369D-02  /
   data green(23, 9, 4) /   3.9973270950448878D-02  /
   data green(23, 9, 5) /   3.9688137673994262D-02  /
   data green(23, 9, 6) /   3.9347853086999840D-02  /
   data green(23, 9, 7) /   3.8956844769883370D-02  /
   data green(23, 9, 8) /   3.8519943043542323D-02  /
   data green(23, 9, 9) /   3.8042232042707007D-02  /
   data green(23,10, 0) /   3.9877859749929752D-02  /
   data green(23,10, 1) /   3.9846097150275442D-02  /
   data green(23,10, 2) /   3.9751266683146748D-02  /
   data green(23,10, 3) /   3.9594721930009437D-02  /
   data green(23,10, 4) /   3.9378658946405944D-02  /
   data green(23,10, 5) /   3.9106031329196482D-02  /
   data green(23,10, 6) /   3.8780440522046621D-02  /
   data green(23,10, 7) /   3.8406009154847133D-02  /
   data green(23,10, 8) /   3.7987245722834795D-02  /
   data green(23,10, 9) /   3.7528908578141978D-02  /
   data green(23,10,10) /   3.7035876193117344D-02  /
   data green(23,11, 0) /   3.9226881037326883D-02  /
   data green(23,11, 1) /   3.9196654528536948D-02  /
   data green(23,11, 2) /   3.9106395867611235D-02  /
   data green(23,11, 3) /   3.8957351273795883D-02  /
   data green(23,11, 4) /   3.8751544448416185D-02  /
   data green(23,11, 5) /   3.8491700825275331D-02  /
   data green(23,11, 6) /   3.8181149463083951D-02  /
   data green(23,11, 7) /   3.7823709342437707D-02  /
   data green(23,11, 8) /   3.7423567322660722D-02  /
   data green(23,11, 9) /   3.6985154787465328D-02  /
   data green(23,11,10) /   3.6513029190601391D-02  /
   data green(23,11,11) /   3.6011765490073171D-02  /
   data green(23,12, 0) /   3.8549442608840134D-02  /
   data green(23,12, 1) /   3.8520760283402844D-02  /
   data green(23,12, 2) /   3.8435098782070340D-02  /
   data green(23,12, 3) /   3.8293600059627037D-02  /
   data green(23,12, 4) /   3.8098120212973971D-02  /
   data green(23,12, 5) /   3.7851162361465150D-02  /
   data green(23,12, 6) /   3.7555789470177223D-02  /
   data green(23,12, 7) /   3.7215522934181712D-02  /
   data green(23,12, 8) /   3.6834233209747533D-02  /
   data green(23,12, 9) /   3.6416028638958160D-02  /
   data green(23,12,10) /   3.5965147965985919D-02  /
   data green(23,12,11) /   3.5485861035728747D-02  /
   data green(23,12,12) /   3.4982380966024719D-02  /
   data green(23,13, 0) /   3.7851677166925869D-02  /
   data green(23,13, 1) /   3.7824528607298492D-02  /
   data green(23,13, 2) /   3.7743434524955648D-02  /
   data green(23,13, 3) /   3.7609436997409755D-02  /
   data green(23,13, 4) /   3.7424231353929964D-02  /
   data green(23,13, 5) /   3.7190107042562558D-02  /
   data green(23,13, 6) /   3.6909870608516519D-02  /
   data green(23,13, 7) /   3.6586755753129176D-02  /
   data green(23,13, 8) /   3.6224325879824985D-02  /
   data green(23,13, 9) /   3.5826374462449538D-02  /
   data green(23,13,10) /   3.5396828066725901D-02  /
   data green(23,13,11) /   3.4939656035569970D-02  /
   data green(23,13,12) /   3.4458789849447002D-02  /
   data green(23,13,13) /   3.3958054123861395D-02  /
   data green(23,14, 0) /   3.7139246659103571D-02  /
   data green(23,14, 1) /   3.7113605937788965D-02  /
   data green(23,14, 2) /   3.7037003317792715D-02  /
   data green(23,14, 3) /   3.6910386318514878D-02  /
   data green(23,14, 4) /   3.6735297859719789D-02  /
   data green(23,14, 5) /   3.6513824423389421D-02  /
   data green(23,14, 6) /   3.6248528343070796D-02  /
   data green(23,14, 7) /   3.5942368438039635D-02  /
   data green(23,14, 8) /   3.5598613612392527D-02  /
   data green(23,14, 9) /   3.5220754020036778D-02  /
   data green(23,14,10) /   3.4812414011042132D-02  /
   data green(23,14,11) /   3.4377270415307977D-02  /
   data green(23,14,12) /   3.3918978895006974D-02  /
   data green(23,14,13) /   3.3441110214767582D-02  /
   data green(23,14,14) /   3.2947097428902325D-02  /
   data green(23,15, 0) /   3.6417293216017707D-02  /
   data green(23,15, 1) /   3.6393121828945538D-02  /
   data green(23,15, 2) /   3.6320897190036165D-02  /
   data green(23,15, 3) /   3.6201478193068320D-02  /
   data green(23,15, 4) /   3.6036264720989857D-02  /
   data green(23,15, 5) /   3.5827152395944431D-02  /
   data green(23,15, 6) /   3.5576473309366016D-02  /
   data green(23,15, 7) /   3.5286926291771714D-02  /
   data green(23,15, 8) /   3.4961500648319210D-02  /
   data green(23,15, 9) /   3.4603397304792685D-02  /
   data green(23,15,10) /   3.4215951020105044D-02  /
   data green(23,15,11) /   3.3802556796979355D-02  /
   data green(23,15,12) /   3.3366602948816236D-02  /
   data green(23,15,13) /   3.2911412544979333D-02  /
   data green(23,15,14) /   3.2440194234886151D-02  /
   data green(23,15,15) /   3.1956002801270054D-02  /
   data green(23,16, 0) /   3.5690415929245628D-02  /
   data green(23,16, 1) /   3.5667665548064297D-02  /
   data green(23,16, 2) /   3.5599676048712350D-02  /
   data green(23,16, 3) /   3.5487223973282023D-02  /
   data green(23,16, 4) /   3.5331576120503944D-02  /
   data green(23,16, 5) /   3.5134450188516693D-02  /
   data green(23,16, 6) /   3.4897963082637824D-02  /
   data green(23,16, 7) /   3.4624569878652471D-02  /
   data green(23,16, 8) /   3.4316996761805815D-02  /
   data green(23,16, 9) /   3.3978171306345596D-02  /
   data green(23,16,10) /   3.3611153249321346D-02  /
   data green(23,16,11) /   3.3219068499913518D-02  /
   data green(23,16,12) /   3.2805048580006477D-02  /
   data green(23,16,13) /   3.2372177083309472D-02  /
   data green(23,16,14) /   3.1923444132070801D-02  /
   data green(23,16,15) /   3.1461709251927678D-02  /
   data green(23,16,16) /   3.0989672610216595D-02  /
   data green(23,17, 0) /   3.4962668772150186D-02  /
   data green(23,17, 1) /   3.4941283740410993D-02  /
   data green(23,17, 2) /   3.4877364591776476D-02  /
   data green(23,17, 3) /   3.4771611909288659D-02  /
   data green(23,17, 4) /   3.4625169581926976D-02  /
   data green(23,17, 5) /   3.4439590674771266D-02  /
   data green(23,17, 6) /   3.4216792483059907D-02  /
   data green(23,17, 7) /   3.3959003272725725D-02  /
   data green(23,17, 8) /   3.3668703503704149D-02  /
   data green(23,17, 9) /   3.3348564393707994D-02  /
   data green(23,17,10) /   3.3001386529812646D-02  /
   data green(23,17,11) /   3.2630040914338583D-02  /
   data green(23,17,12) /   3.2237414393343412D-02  /
   data green(23,17,13) /   3.1826360916750243D-02  /
   data green(23,17,14) /   3.1399659570449427D-02  /
   data green(23,17,15) /   3.0959979844625351D-02  /
   data green(23,17,16) /   3.0509854188816705D-02  /
   data green(23,17,17) /   3.0051657569971225D-02  /
   data green(23,18, 0) /   3.4237574949950923D-02  /
   data green(23,18, 1) /   3.4217494481708959D-02  /
   data green(23,18, 2) /   3.4157465484170486D-02  /
   data green(23,18, 3) /   3.4058118917466226D-02  /
   data green(23,18, 4) /   3.3920485873803169D-02  /
   data green(23,18, 5) /   3.3745968052366841D-02  /
   data green(23,18, 6) /   3.3536298776840258D-02  /
   data green(23,18, 7) /   3.3293496643072400D-02  /
   data green(23,18, 8) /   3.3019814143891546D-02  /
   data green(23,18, 9) /   3.2717683688993499D-02  /
   data green(23,18,10) /   3.2389663334391999D-02  /
   data green(23,18,11) /   3.2038384289005295D-02  /
   data green(23,18,12) /   3.1666501916755262D-02  /
   data green(23,18,13) /   3.1276651545678012D-02  /
   data green(23,18,14) /   3.0871409973061332D-02  /
   data green(23,18,15) /   3.0453263152743423D-02  /
   data green(23,18,16) /   3.0024580193396315D-02  /
   data green(23,18,17) /   2.9587593500983746D-02  /
   data green(23,18,18) /   2.9144384671669289D-02  /
   data green(23,19, 0) /   3.3518153290375781D-02  /
   data green(23,19, 1) /   3.3499313356360791D-02  /
   data green(23,19, 2) /   3.3442984513589555D-02  /
   data green(23,19, 3) /   3.3349734247082270D-02  /
   data green(23,19, 4) /   3.3220490688977394D-02  /
   data green(23,19, 5) /   3.3056517125702657D-02  /
   data green(23,19, 6) /   3.2859378255174453D-02  /
   data green(23,19, 7) /   3.2630899930855439D-02  /
   data green(23,19, 8) /   3.2373124357024601D-02  /
   data green(23,19, 9) /   3.2088262774707069D-02  /
   data green(23,19,10) /   3.1778647609832164D-02  /
   data green(23,19,11) /   3.1446685867295113D-02  /
   data green(23,19,12) /   3.1094815278351114D-02  /
   data green(23,19,13) /   3.0725464379326710D-02  /
   data green(23,19,14) /   3.0341017350961150D-02  /
   data green(23,19,15) /   2.9943784109026587D-02  /
   data green(23,19,16) /   2.9535975830586508D-02  /
   data green(23,19,17) /   2.9119685841083504D-02  /
   data green(23,19,18) /   2.8696875582868275D-02  /
   data green(23,19,19) /   2.8269365237093756D-02  /
   data green(23,20, 0) /   3.2806952816533709D-02  /
   data green(23,20, 1) /   3.2789287722767610D-02  /
   data green(23,20, 2) /   3.2736463946021768D-02  /
   data green(23,20, 3) /   3.2648991361586395D-02  /
   data green(23,20, 4) /   3.2527704547247203D-02  /
   data green(23,20, 5) /   3.2373740803913073D-02  /
   data green(23,20, 6) /   3.2188510993154902D-02  /
   data green(23,20, 7) /   3.1973664633806013D-02  /
   data green(23,20, 8) /   3.1731050898036760D-02  /
   data green(23,20, 9) /   3.1462677222696092D-02  /
   data green(23,20,10) /   3.1170667210307153D-02  /
   data green(23,20,11) /   3.0857219352808846D-02  /
   data green(23,20,12) /   3.0524567894136906D-02  /
   data green(23,20,13) /   3.0174946882568044D-02  /
   data green(23,20,14) /   2.9810558177652661D-02  /
   data green(23,20,15) /   2.9433543893601130D-02  /
   data green(23,20,16) /   2.9045963500213056D-02  /
   data green(23,20,17) /   2.8649775577062768D-02  /
   data green(23,20,18) /   2.8246824034124249D-02  /
   data green(23,20,19) /   2.7838828474725325D-02  /
   data green(23,20,20) /   2.7427378283055789D-02  /
   data green(23,21, 0) /   3.2106092268877812D-02  /
   data green(23,21, 1) /   3.2089535946118436D-02  /
   data green(23,21, 2) /   3.2040020901774123D-02  /
   data green(23,21, 3) /   3.1958004921915292D-02  /
   data green(23,21, 4) /   3.1844237900607406D-02  /
   data green(23,21, 5) /   3.1699742909330232D-02  /
   data green(23,21, 6) /   3.1525791021481307D-02  /
   data green(23,21, 7) /   3.1323871087657863D-02  /
   data green(23,21, 8) /   3.1095655829247001D-02  /
   data green(23,21, 9) /   3.0842965690860813D-02  /
   data green(23,21,10) /   3.0567731870055109D-02  /
   data green(23,21,11) /   3.0271959837892295D-02  /
   data green(23,21,12) /   2.9957694494614315D-02  /
   data green(23,21,13) /   2.9626987892425960D-02  /
   data green(23,21,14) /   2.9281870223863835D-02  /
   data green(23,21,15) /   2.8924324539105609D-02  /
   data green(23,21,16) /   2.8556265434912326D-02  /
   data green(23,21,17) /   2.8179521763331026D-02  /
   data green(23,21,18) /   2.7795823246937114D-02  /
   data green(23,21,19) /   2.7406790762284524D-02  /
   data green(23,21,20) /   2.7013929964012140D-02  /
   data green(23,21,21) /   2.6618627865913101D-02  /
   data green(23,22, 0) /   3.1417301982889063D-02  /
   data green(23,22, 1) /   3.1401789013586347D-02  /
   data green(23,22, 2) /   3.1355388191884992D-02  /
   data green(23,22, 3) /   3.1278510353415778D-02  /
   data green(23,22, 4) /   3.1171828981171944D-02  /
   data green(23,22, 5) /   3.1036263911579692D-02  /
   data green(23,22, 6) /   3.0872959617298062D-02  /
   data green(23,22, 7) /   3.0683259057099663D-02  /
   data green(23,22, 8) /   3.0468674231441117D-02  /
   data green(23,22, 9) /   3.0230854651245129D-02  /
   data green(23,22,10) /   2.9971554919070293D-02  /
   data green(23,22,11) /   2.9692602545146357D-02  /
   data green(23,22,12) /   2.9395866989568177D-02  /
   data green(23,22,13) /   2.9083230752891775D-02  /
   data green(23,22,14) /   2.8756563147670007D-02  /
   data green(23,22,15) /   2.8417697189067129D-02  /
   data green(23,22,16) /   2.8068409857007291D-02  /
   data green(23,22,17) /   2.7710405815444986D-02  /
   data green(23,22,18) /   2.7345304532902289D-02  /
   data green(23,22,19) /   2.6974630635702251D-02  /
   data green(23,22,20) /   2.6599807241865437D-02  /
   data green(23,22,21) /   2.6222151967802171D-02  /
   data green(23,22,22) /   2.5842875268638859D-02  /
   data green(23,23, 0) /   3.0741966130337195D-02  /
   data green(23,23, 1) /   3.0727432544024446D-02  /
   data green(23,23, 2) /   3.0683955639893304D-02  /
   data green(23,23, 3) /   3.0611904043348657D-02  /
   data green(23,23, 4) /   3.0511882470504653D-02  /
   data green(23,23, 5) /   3.0384717708280867D-02  /
   data green(23,23, 6) /   3.0231439888383758D-02  /
   data green(23,23, 7) /   3.0053259874960637D-02  /
   data green(23,23, 8) /   2.9851543713143146D-02  /
   data green(23,23, 9) /   2.9627785149506154D-02  /
   data green(23,23,10) /   2.9383577236634210D-02  /
   data green(23,23,11) /   2.9120583978884711D-02  /
   data green(23,23,12) /   2.8840512875477810D-02  /
   data green(23,23,13) /   2.8545089083056870D-02  /
   data green(23,23,14) /   2.8236031766422431D-02  /
   data green(23,23,15) /   2.7915033046131237D-02  /
   data green(23,23,16) /   2.7583739796179757D-02  /
   data green(23,23,17) /   2.7243738402722775D-02  /
   data green(23,23,18) /   2.6896542471694533D-02  /
   data green(23,23,19) /   2.6543583372661734D-02  /
   data green(23,23,20) /   2.6186203429321250D-02  /
   data green(23,23,21) /   2.5825651513008637D-02  /
   data green(23,23,22) /   2.5463080762339923D-02  /
   data green(23,23,23) /   2.5099548136777599D-02  /
   data green(24, 0, 0) /   4.1684842124130525D-02  /
   data green(24, 1, 0) /   4.1648514426701956D-02  /
   data green(24, 1, 1) /   4.1612282380419474D-02  /
   data green(24, 2, 0) /   4.1540105249665234D-02  /
   data green(24, 2, 1) /   4.1504157611061326D-02  /
   data green(24, 2, 2) /   4.1396878523711239D-02  /
   data green(24, 3, 0) /   4.1361310879756960D-02  /
   data green(24, 3, 1) /   4.1325828928511492D-02  /
   data green(24, 3, 2) /   4.1219934647749597D-02  /
   data green(24, 3, 3) /   4.1045258694370340D-02  /
   data green(24, 4, 0) /   4.1114875844454696D-02  /
   data green(24, 4, 1) /   4.1080028915039202D-02  /
   data green(24, 4, 2) /   4.0976023172499543D-02  /
   data green(24, 4, 3) /   4.0804440709280763D-02  /
   data green(24, 4, 4) /   4.0567843502005325D-02  /
   data green(24, 5, 0) /   4.0804476514498457D-02  /
   data green(24, 5, 1) /   4.0770418240690887D-02  /
   data green(24, 5, 2) /   4.0668758236232727D-02  /
   data green(24, 5, 3) /   4.0501019133711375D-02  /
   data green(24, 5, 4) /   4.0269667761836402D-02  /
   data green(24, 5, 5) /   3.9978014971612144D-02  /
   data green(24, 6, 0) /   4.0434572077734598D-02  /
   data green(24, 6, 1) /   4.0401437523710546D-02  /
   data green(24, 6, 2) /   4.0302525405050783D-02  /
   data green(24, 6, 3) /   4.0139289918620161D-02  /
   data green(24, 6, 4) /   3.9914088299381242D-02  /
   data green(24, 6, 5) /   3.9630086844198428D-02  /
   data green(24, 6, 6) /   3.9291139882013652D-02  /
   data green(24, 7, 0) /   4.0010234440712208D-02  /
   data green(24, 7, 1) /   3.9978138204587947D-02  /
   data green(24, 7, 2) /   3.9882315420961510D-02  /
   data green(24, 7, 3) /   3.9724144922372265D-02  /
   data green(24, 7, 4) /   3.9505863150203734D-02  /
   data green(24, 7, 5) /   3.9230476870770997D-02  /
   data green(24, 7, 6) /   3.8901650548712939D-02  /
   data green(24, 7, 7) /   3.8523576458390323D-02  /
   data green(24, 8, 0) /   3.9536969095789187D-02  /
   data green(24, 8, 1) /   3.9506004385455575D-02  /
   data green(24, 8, 2) /   3.9413548917295683D-02  /
   data green(24, 8, 3) /   3.9260901286851890D-02  /
   data green(24, 8, 4) /   3.9050169209078951D-02  /
   data green(24, 8, 5) /   3.8784189202111262D-02  /
   data green(24, 8, 6) /   3.8466422734236172D-02  /
   data green(24, 8, 7) /   3.8100836120113550D-02  /
   data green(24, 8, 8) /   3.7691771949567944D-02  /
   data green(24, 9, 0) /   3.9020538168949528D-02  /
   data green(24, 9, 1) /   3.8990776776639882D-02  /
   data green(24, 9, 2) /   3.8901903026546579D-02  /
   data green(24, 9, 3) /   3.8755132367127382D-02  /
   data green(24, 9, 4) /   3.8552438980722310D-02  /
   data green(24, 9, 5) /   3.8296482511022419D-02  /
   data green(24, 9, 6) /   3.7990513106049215D-02  /
   data green(24, 9, 7) /   3.7638261275912320D-02  /
   data green(24, 9, 8) /   3.7243819547141961D-02  /
   data green(24, 9, 9) /   3.6811522689534103D-02  /
   data green(24,10, 0) /   3.8466795055696001D-02  /
   data green(24,10, 1) /   3.8438288103828974D-02  /
   data green(24,10, 2) /   3.8353149074489043D-02  /
   data green(24,10, 3) /   3.8212509149378744D-02  /
   data green(24,10, 4) /   3.8018207032842581D-02  /
   data green(24,10, 5) /   3.7772722627699526D-02  /
   data green(24,10, 6) /   3.7479090881570795D-02  /
   data green(24,10, 7) /   3.7140801544375497D-02  /
   data green(24,10, 8) /   3.6761691039892415D-02  /
   data green(24,10, 9) /   3.6345832517833625D-02  /
   data green(24,10,10) /   3.5897429514934916D-02  /
   data green(24,11, 0) /   3.7881537665644865D-02  /
   data green(24,11, 1) /   3.7854316964570152D-02  /
   data green(24,11, 2) /   3.7773008254295026D-02  /
   data green(24,11, 3) /   3.7638658894972081D-02  /
   data green(24,11, 4) /   3.7452972676523921D-02  /
   data green(24,11, 5) /   3.7218250208478534D-02  /
   data green(24,11, 6) /   3.6937311300508889D-02  /
   data green(24,11, 7) /   3.6613404359228670D-02  /
   data green(24,11, 8) /   3.6250108263964478D-02  /
   data green(24,11, 9) /   3.5851232106247966D-02  /
   data green(24,11,10) /   3.5420717661480243D-02  /
   data green(24,11,11) /   3.4962548627429667D-02  /
   data green(24,12, 0) /   3.7270384731854983D-02  /
   data green(24,12, 1) /   3.7244464578980160D-02  /
   data green(24,12, 2) /   3.7167029691699098D-02  /
   data green(24,12, 3) /   3.7039045364217531D-02  /
   data green(24,12, 4) /   3.6862083144984882D-02  /
   data green(24,12, 5) /   3.6638267600465964D-02  /
   data green(24,12, 6) /   3.6370206824612297D-02  /
   data green(24,12, 7) /   3.6060911058376770D-02  /
   data green(24,12, 8) /   3.5713704192043651D-02  /
   data green(24,12, 9) /   3.5332132892574940D-02  /
   data green(24,12,10) /   3.4919877688234728D-02  /
   data green(24,12,11) /   3.4480669651125187D-02  /
   data green(24,12,12) /   3.4018215459376534D-02  /
   data green(24,13, 0) /   3.6638677227046217D-02  /
   data green(24,13, 1) /   3.6614056482698540D-02  /
   data green(24,13, 2) /   3.6540492965718568D-02  /
   data green(24,13, 3) /   3.6418872708114523D-02  /
   data green(24,13, 4) /   3.6250639385268693D-02  /
   data green(24,13, 5) /   3.6037747039583776D-02  /
   data green(24,13, 6) /   3.5782598216329918D-02  /
   data green(24,13, 7) /   3.5487971274824760D-02  /
   data green(24,13, 8) /   3.5156941015261296D-02  /
   data green(24,13, 9) /   3.4792796767835754D-02  /
   data green(24,13,10) /   3.4398961771740316D-02  /
   data green(24,13,11) /   3.3978917105038794D-02  /
   data green(24,13,12) /   3.3536132706053524D-02  /
   data green(24,13,13) /   3.3074007245920675D-02  /
   data green(24,14, 0) /   3.5991404888373789D-02  /
   data green(24,14, 1) /   3.5968069181719423D-02  /
   data green(24,14, 2) /   3.5898335156902002D-02  /
   data green(24,14, 3) /   3.5783013182267301D-02  /
   data green(24,14, 4) /   3.5623424722259156D-02  /
   data green(24,14, 5) /   3.5421360561224420D-02  /
   data green(24,14, 6) /   3.5179026000967062D-02  /
   data green(24,14, 7) /   3.4898976256391483D-02  /
   data green(24,14, 8) /   3.4584045618185892D-02  /
   data green(24,14, 9) /   3.4237273985266402D-02  /
   data green(24,14,10) /   3.3861834126598343D-02  /
   data green(24,14,11) /   3.3460962573495075D-02  /
   data green(24,14,12) /   3.3037896445331769D-02  /
   data green(24,14,13) /   3.2595817850987245D-02  /
   data green(24,14,14) /   3.2137806853577623D-02  /
   data green(24,15, 0) /   3.5333156298707091D-02  /
   data green(24,15, 1) /   3.5311080242092605D-02  /
   data green(24,15, 2) /   3.5245100969552537D-02  /
   data green(24,15, 3) /   3.5135957346446769D-02  /
   data green(24,15, 4) /   3.4984855211617649D-02  /
   data green(24,15, 5) /   3.4793430619965024D-02  /
   data green(24,15, 6) /   3.4563701505052176D-02  /
   data green(24,15, 7) /   3.4298010509261470D-02  /
   data green(24,15, 8) /   3.3998962040791558D-02  /
   data green(24,15, 9) /   3.3669356669809983D-02  /
   data green(24,15,10) /   3.3312125795259925D-02  /
   data green(24,15,11) /   3.2930269147062721D-02  /
   data green(24,15,12) /   3.2526797196460747D-02  /
   data green(24,15,13) /   3.2104679993262226D-02  /
   data green(24,15,14) /   3.1666803390141537D-02  /
   data green(24,15,15) /   3.1215933096726659D-02  /
   data green(24,16, 0) /   3.4668089926557538D-02  /
   data green(24,16, 1) /   3.4647239244400045D-02  /
   data green(24,16, 2) /   3.4584913439557025D-02  /
   data green(24,16, 3) /   3.4481784392591447D-02  /
   data green(24,16, 4) /   3.4338949497951278D-02  /
   data green(24,16, 5) /   3.4157899441526800D-02  /
   data green(24,16, 6) /   3.3940475724207032D-02  /
   data green(24,16, 7) /   3.3688820261546203D-02  /
   data green(24,16, 8) /   3.3405319669213306D-02  /
   data green(24,16, 9) /   3.3092546909379406D-02  /
   data green(24,16,10) /   3.2753202842476646D-02  /
   data green(24,16,11) /   3.2390059938768298D-02  /
   data green(24,16,12) /   3.2005910003080645D-02  /
   data green(24,16,13) /   3.1603517305233968D-02  /
   data green(24,16,14) /   3.1185578035860792D-02  /
   data green(24,16,15) /   3.0754686561601897D-02  /
   data green(24,16,16) /   3.0313308563292996D-02  /
   data green(24,17, 0) /   3.3999922944564290D-02  /
   data green(24,17, 1) /   3.3980256450082631D-02  /
   data green(24,17, 2) /   3.3921462155082203D-02  /
   data green(24,17, 3) /   3.3824149657881164D-02  /
   data green(24,17, 4) /   3.3689315407357444D-02  /
   data green(24,17, 5) /   3.3518314543981895D-02  /
   data green(24,17, 6) /   3.3312823691769140D-02  /
   data green(24,17, 7) /   3.3074796669274396D-02  /
   data green(24,17, 8) /   3.2806415336456639D-02  /
   data green(24,17, 9) /   3.2510037866721408D-02  /
   data green(24,17,10) /   3.2188146642281069D-02  /
   data green(24,17,11) /   3.1843297744303979D-02  /
   data green(24,17,12) /   3.1478073685080876D-02  /
   data green(24,17,13) /   3.1095040648991077D-02  /
   data green(24,17,14) /   3.0696711111704526D-02  /
   data green(24,17,15) /   3.0285512325987508D-02  /
   data green(24,17,16) /   2.9863760822436344D-02  /
   data green(24,17,17) /   2.9433642790107364D-02  /
   data green(24,18, 0) /   3.3331934445746149D-02  /
   data green(24,18, 1) /   3.3313405823202760D-02  /
   data green(24,18, 2) /   3.3258005705842628D-02  /
   data green(24,18, 3) /   3.3166286155482463D-02  /
   data green(24,18, 4) /   3.3039150263877506D-02  /
   data green(24,18, 5) /   3.2877827608524771D-02  /
   data green(24,18, 6) /   3.2683841745867262D-02  /
   data green(24,18, 7) /   3.2458971398642622D-02  /
   data green(24,18, 8) /   3.2205207213786578D-02  /
   data green(24,18, 9) /   3.1924706042601850D-02  /
   data green(24,18,10) /   3.1619744634365662D-02  /
   data green(24,18,11) /   3.1292674459457714D-02  /
   data green(24,18,12) /   3.0945879118052983D-02  /
   data green(24,18,13) /   3.0581735478510102D-02  /
   data green(24,18,14) /   3.0202579357927226D-02  /
   data green(24,18,15) /   2.9810676233764186D-02  /
   data green(24,18,16) /   2.9408197181420150D-02  /
   data green(24,18,17) /   2.8997199982545116D-02  /
   data green(24,18,18) /   2.8579615150136611D-02  /
   data green(24,19, 0) /   3.2666979761297218D-02  /
   data green(24,19, 1) /   3.2649539131817111D-02  /
   data green(24,19, 2) /   3.2597385145596972D-02  /
   data green(24,19, 3) /   3.2511017006087890D-02  /
   data green(24,19, 4) /   3.2391251946045613D-02  /
   data green(24,19, 5) /   3.2239203879142884D-02  /
   data green(24,19, 6) /   3.2056255062120544D-02  /
   data green(24,19, 7) /   3.1844022158207613D-02  /
   data green(24,19, 8) /   3.1604318285362830D-02  /
   data green(24,19, 9) /   3.1339112707877717D-02  /
   data green(24,19,10) /   3.1050489792619949D-02  /
   data green(24,19,11) /   3.0740608717608954D-02  /
   data green(24,19,12) /   3.0411665213658134D-02  /
   data green(24,19,13) /   3.0065856365738854D-02  /
   data green(24,19,14) /   2.9705349225601616D-02  /
   data green(24,19,15) /   2.9332253714178243D-02  /
   data green(24,19,16) /   2.8948600039814368D-02  /
   data green(24,19,17) /   2.8556320639336320D-02  /
   data green(24,19,18) /   2.8157236470736056D-02  /
   data green(24,19,19) /   2.7753047351389788D-02  /
   data green(24,20, 0) /   3.2007512859844739D-02  /
   data green(24,20, 1) /   3.1991108124895480D-02  /
   data green(24,20, 2) /   3.1942045511411922D-02  /
   data green(24,20, 3) /   3.1860775891654293D-02  /
   data green(24,20, 4) /   3.1748037910886102D-02  /
   data green(24,20, 5) /   3.1604839449606910D-02  /
   data green(24,20, 6) /   3.1432432962903045D-02  /
   data green(24,20, 7) /   3.1232285861726317D-02  /
   data green(24,20, 8) /   3.1006047269250658D-02  /
   data green(24,20, 9) /   3.0755512557920188D-02  /
   data green(24,20,10) /   3.0482587052852877D-02  /
   data green(24,20,11) /   3.0189250186703525D-02  /
   data green(24,20,12) /   2.9877521227590253D-02  /
   data green(24,20,13) /   2.9549427495964436D-02  /
   data green(24,20,14) /   2.9206975759379351D-02  /
   data green(24,20,15) /   2.8852127265099239D-02  /
   data green(24,20,16) /   2.8486776655044500D-02  /
   data green(24,20,17) /   2.8112734817104083D-02  /
   data green(24,20,18) /   2.7731715568414986D-02  /
   data green(24,20,19) /   2.7345325942880808D-02  /
   data green(24,20,20) /   2.6955059766769221D-02  /
   data green(24,21, 0) /   3.1355614194759650D-02  /
   data green(24,21, 1) /   3.1340192163157379D-02  /
   data green(24,21, 2) /   3.1294062813181533D-02  /
   data green(24,21, 3) /   3.1217633002579143D-02  /
   data green(24,21, 4) /   3.1111569735805419D-02  /
   data green(24,21, 5) /   3.0976784087207446D-02  /
   data green(24,21, 6) /   3.0814409769182309D-02  /
   data green(24,21, 7) /   3.0625777319197935D-02  /
   data green(24,21, 8) /   3.0412385025844663D-02  /
   data green(24,21, 9) /   3.0175867782546672D-02  /
   data green(24,21,10) /   2.9917965050184136D-02  /
   data green(24,21,11) /   2.9640489035310892D-02  /
   data green(24,21,12) /   2.9345294062424129D-02  /
   data green(24,21,13) /   2.9034247953104751D-02  /
   data green(24,21,14) /   2.8709206038649602D-02  /
   data green(24,21,15) /   2.8371988241730321D-02  /
   data green(24,21,16) /   2.8024359479819204D-02  /
   data green(24,21,17) /   2.7668013478594766D-02  /
   data green(24,21,18) /   2.7304559943820687D-02  /
   data green(24,21,19) /   2.6935514928604447D-02  /
   data green(24,21,20) /   2.6562294150062771D-02  /
   data green(24,21,21) /   2.6186208953714332D-02  /
   data green(24,22, 0) /   3.0713021799417789D-02  /
   data green(24,22, 1) /   3.0698529111751242D-02  /
   data green(24,22, 2) /   3.0655174325336806D-02  /
   data green(24,22, 3) /   3.0583324350017715D-02  /
   data green(24,22, 4) /   3.0483581105180874D-02  /
   data green(24,22, 5) /   3.0356767589786116D-02  /
   data green(24,22, 6) /   3.0203909274866200D-02  /
   data green(24,22, 7) /   3.0026211631709990D-02  /
   data green(24,22, 8) /   2.9825034735537742D-02  /
   data green(24,22, 9) /   2.9601865948042424D-02  /
   data green(24,22,10) /   2.9358291683654780D-02  /
   data green(24,22,11) /   2.9095969210068011D-02  /
   data green(24,22,12) /   2.8816599333726888D-02  /
   data green(24,22,13) /   2.8521900688326236D-02  /
   data green(24,22,14) /   2.8213586192315523D-02  /
   data green(24,22,15) /   2.7893342082734772D-02  /
   data green(24,22,16) /   2.7562809778427497D-02  /
   data green(24,22,17) /   2.7223570684398676D-02  /
   data green(24,22,18) /   2.6877133926758570D-02  /
   data green(24,22,19) /   2.6524926907674010D-02  /
   data green(24,22,20) /   2.6168288493134226D-02  /
   data green(24,22,21) /   2.5808464592402311D-02  /
   data green(24,22,22) /   2.5446605854737160D-02  /
   data green(24,23, 0) /   3.0081163864098075D-02  /
   data green(24,23, 1) /   3.0067547733371679D-02  /
   data green(24,23, 2) /   3.0026810434108351D-02  /
   data green(24,23, 3) /   2.9959282721431158D-02  /
   data green(24,23, 4) /   2.9865507554527250D-02  /
   data green(24,23, 5) /   2.9746228033818680D-02  /
   data green(24,23, 6) /   2.9602371256232447D-02  /
   data green(24,23, 7) /   2.9435028767018766D-02  /
   data green(24,23, 8) /   2.9245434395754940D-02  /
   data green(24,23, 9) /   2.9034940322384597D-02  /
   data green(24,23,10) /   2.8804992226636888D-02  /
   data green(24,23,11) /   2.8557104335451777D-02  /
   data green(24,23,12) /   2.8292835105894462D-02  /
   data green(24,23,13) /   2.8013764175328752D-02  /
   data green(24,23,14) /   2.7721471087003180D-02  /
   data green(24,23,15) /   2.7417516168025444D-02  /
   data green(24,23,16) /   2.7103423807104975D-02  /
   data green(24,23,17) /   2.6780668258786142D-02  /
   data green(24,23,18) /   2.6450661994488998D-02  /
   data green(24,23,19) /   2.6114746531775097D-02  /
   data green(24,23,20) /   2.5774185603263115D-02  /
   data green(24,23,21) /   2.5430160475393180D-02  /
   data green(24,23,22) /   2.5083767193487767D-02  /
   data green(24,23,23) /   2.4736015511215261D-02  /
   data green(24,24, 0) /   2.9461191432220299D-02  /
   data green(24,24, 1) /   2.9448400221982588D-02  /
   data green(24,24, 2) /   2.9410126688015823D-02  /
   data green(24,24, 3) /   2.9346668939890887D-02  /
   data green(24,24, 4) /   2.9258516651404071D-02  /
   data green(24,24, 5) /   2.9146340617625114D-02  /
   data green(24,24, 6) /   2.9010978752128481D-02  /
   data green(24,24, 7) /   2.8853419090136414D-02  /
   data green(24,24, 8) /   2.8674780457209920D-02  /
   data green(24,24, 9) /   2.8476291515847466D-02  /
   data green(24,24,10) /   2.8259268913705161D-02  /
   data green(24,24,11) /   2.8025095230329915D-02  /
   data green(24,24,12) /   2.7775197360154043D-02  /
   data green(24,24,13) /   2.7511025885668963D-02  /
   data green(24,24,14) /   2.7234035894599714D-02  /
   data green(24,24,15) /   2.6945669586887898D-02  /
   data green(24,24,16) /   2.6647340908893376D-02  /
   data green(24,24,17) /   2.6340422349680641D-02  /
   data green(24,24,18) /   2.6026233942242436D-02  /
   data green(24,24,19) /   2.5706034434083126D-02  /
   data green(24,24,20) /   2.5381014528293756D-02  /
   data green(24,24,21) /   2.5052292048361756D-02  /
   data green(24,24,22) /   2.4720908846732524D-02  /
   data green(24,24,23) /   2.4387829257131639D-02  /
   data green(24,24,24) /   2.4053939881981679D-02  /
   data green(25, 0, 0) /   4.0016074018692789D-02  /
   data green(25, 1, 0) /   3.9983944719202150D-02  /
   data green(25, 1, 1) /   3.9951893315804662D-02  /
   data green(25, 2, 0) /   3.9888024203032044D-02  /
   data green(25, 2, 1) /   3.9856204580579438D-02  /
   data green(25, 2, 2) /   3.9761205532449453D-02  /
   data green(25, 3, 0) /   3.9729695542854375D-02  /
   data green(25, 3, 1) /   3.9698255974639433D-02  /
   data green(25, 3, 2) /   3.9604387878852752D-02  /
   data green(25, 3, 3) /   3.9449424972202517D-02  /
   data green(25, 4, 0) /   3.9511201825418708D-02  /
   data green(25, 4, 1) /   3.9480281588574535D-02  /
   data green(25, 4, 2) /   3.9387959028794023D-02  /
   data green(25, 4, 3) /   3.9235531194569170D-02  /
   data green(25, 4, 4) /   3.9025103192562798D-02  /
   data green(25, 5, 0) /   3.9235558301578738D-02  /
   data green(25, 5, 1) /   3.9205284788140622D-02  /
   data green(25, 5, 2) /   3.9114887090842008D-02  /
   data green(25, 5, 3) /   3.8965617190560874D-02  /
   data green(25, 5, 4) /   3.8759507845853525D-02  /
   data green(25, 5, 5) /   3.8499296106467615D-02  /
   data green(25, 6, 0) /   3.8906439084312607D-02  /
   data green(25, 6, 1) /   3.8876925534890951D-02  /
   data green(25, 6, 2) /   3.8788790034844603D-02  /
   data green(25, 6, 3) /   3.8643232452921672D-02  /
   data green(25, 6, 4) /   3.8442201837446226D-02  /
   data green(25, 6, 5) /   3.8188324319021208D-02  /
   data green(25, 6, 6) /   3.7884809662155808D-02  /
   data green(25, 7, 0) /   3.8528046561734836D-02  /
   data green(25, 7, 1) /   3.8499390491332691D-02  /
   data green(25, 7, 2) /   3.8413807844934075D-02  /
   data green(25, 7, 3) /   3.8272440795764112D-02  /
   data green(25, 7, 4) /   3.8077145645002808D-02  /
   data green(25, 7, 5) /   3.7830425496119112D-02  /
   data green(25, 7, 6) /   3.7535342848669148D-02  /
   data green(25, 7, 7) /   3.7195417971599962D-02  /
   data green(25, 8, 0) /   3.8104972189434455D-02  /
   data green(25, 8, 1) /   3.8077254517285629D-02  /
   data green(25, 8, 2) /   3.7994466100354299D-02  /
   data green(25, 8, 3) /   3.7857687319748576D-02  /
   data green(25, 8, 4) /   3.7668675063067446D-02  /
   data green(25, 8, 5) /   3.7429800404154036D-02  /
   data green(25, 8, 6) /   3.7143967550216789D-02  /
   data green(25, 8, 7) /   3.6814519377567620D-02  /
   data green(25, 8, 8) /   3.6445135321988252D-02  /
   data green(25, 9, 0) /   3.7642056934339273D-02  /
   data green(25, 9, 1) /   3.7615341778691337D-02  /
   data green(25, 9, 2) /   3.7535539054785710D-02  /
   data green(25, 9, 3) /   3.7403664696332273D-02  /
   data green(25, 9, 4) /   3.7221371813813582D-02  /
   data green(25, 9, 5) /   3.6990893475557923D-02  /
   data green(25, 9, 6) /   3.6714968146553825D-02  /
   data green(25, 9, 7) /   3.6396752567910175D-02  /
   data green(25, 9, 8) /   3.6039727284402551D-02  /
   data green(25, 9, 9) /   3.5647599964546473D-02  /
   data green(25,10, 0) /   3.7144258541993790D-02  /
   data green(25,10, 1) /   3.7118593601303242D-02  /
   data green(25,10, 2) /   3.7041919229738172D-02  /
   data green(25,10, 3) /   3.6915185596831893D-02  /
   data green(25,10, 4) /   3.6739939816325000D-02  /
   data green(25,10, 5) /   3.6518273803086303D-02  /
   data green(25,10, 6) /   3.6252756188744940D-02  /
   data green(25,10, 7) /   3.5946352554801587D-02  /
   data green(25,10, 8) /   3.5602338643204022D-02  /
   data green(25,10, 9) /   3.5224211178250336D-02  /
   data green(25,10,10) /   3.4815600537864401D-02  /
   data green(25,11, 0) /   3.6616531239851001D-02  /
   data green(25,11, 1) /   3.6591948656134868D-02  /
   data green(25,11, 2) /   3.6518499032428521D-02  /
   data green(25,11, 3) /   3.6397066651621311D-02  /
   data green(25,11, 4) /   3.6229092329449278D-02  /
   data green(25,11, 5) /   3.6016526220397513D-02  /
   data green(25,11, 6) /   3.5761766066078983D-02  /
   data green(25,11, 7) /   3.5467584645708629D-02  /
   data green(25,11, 8) /   3.5137050563713378D-02  /
   data green(25,11, 9) /   3.4773446514169026D-02  /
   data green(25,11,10) /   3.4380188841432721D-02  /
   data green(25,11,11) /   3.3960751649293935D-02  /
   data green(25,12, 0) /   3.6063721731638437D-02  /
   data green(25,12, 1) /   3.6040239319655200D-02  /
   data green(25,12, 2) /   3.5970068204912969D-02  /
   data green(25,12, 3) /   3.5854027685907480D-02  /
   data green(25,12, 4) /   3.5693453630459109D-02  /
   data green(25,12, 5) /   3.5490156023671206D-02  /
   data green(25,12, 6) /   3.5246363301814930D-02  /
   data green(25,12, 7) /   3.4964656765369542D-02  /
   data green(25,12, 8) /   3.4647898712336790D-02  /
   data green(25,12, 9) /   3.4299157962129380D-02  /
   data green(25,12,10) /   3.3921636186501715D-02  /
   data green(25,12,11) /   3.3518597990942631D-02  /
   data green(25,12,12) /   3.3093307075877001D-02  /
   data green(25,13, 0) /   3.5490483609725593D-02  /
   data green(25,13, 1) /   3.5468106335428351D-02  /
   data green(25,13, 2) /   3.5401229229941475D-02  /
   data green(25,13, 3) /   3.5290608353805188D-02  /
   data green(25,13, 4) /   3.5137477342739865D-02  /
   data green(25,13, 5) /   3.4943509431650928D-02  /
   data green(25,13, 6) /   3.4710767549475087D-02  /
   data green(25,13, 7) /   3.4441645349459898D-02  /
   data green(25,13, 8) /   3.4138802358852686D-02  /
   data green(25,13, 9) /   3.3805096479721922D-02  /
   data green(25,13,10) /   3.3443516876008425D-02  /
   data green(25,13,11) /   3.3057119892248674D-02  /
   data green(25,13,12) /   3.2648970131206344D-02  /
   data green(25,13,13) /   3.2222088237581034D-02  /
   data green(25,14, 0) /   3.4901210781946863D-02  /
   data green(25,14, 1) /   3.4879932382460621D-02  /
   data green(25,14, 2) /   3.4816331324361433D-02  /
   data green(25,14, 3) /   3.4711102841366345D-02  /
   data green(25,14, 4) /   3.4565382133456994D-02  /
   data green(25,14, 5) /   3.4380710561182033D-02  /
   data green(25,14, 6) /   3.4158991123733655D-02  /
   data green(25,14, 7) /   3.3902435697070112D-02  /
   data green(25,14, 8) /   3.3613506799011983D-02  /
   data green(25,14, 9) /   3.3294856709131605D-02  /
   data green(25,14,10) /   3.2949266622786140D-02  /
   data green(25,14,11) /   3.2579588201630816D-02  /
   data green(25,14,12) /   3.2188689449986888D-02  /
   data green(25,14,13) /   3.1779406353063526D-02  /
   data green(25,14,14) /   3.1354501210305821D-02  /
   data green(25,15, 0) /   3.4299989277108166D-02  /
   data green(25,15, 1) /   3.4279793929230938D-02  /
   data green(25,15, 2) /   3.4219422442280614D-02  /
   data green(25,15, 3) /   3.4119512129583404D-02  /
   data green(25,15, 4) /   3.3981104362093693D-02  /
   data green(25,15, 5) /   3.3805614604424654D-02  /
   data green(25,15, 6) /   3.3594792866546347D-02  /
   data green(25,15, 7) /   3.3350676700015673D-02  /
   data green(25,15, 8) /   3.3075539129959468D-02  /
   data green(25,15, 9) /   3.2771833983369913D-02  /
   data green(25,15,10) /   3.2442140965506759D-02  /
   data green(25,15,11) /   3.2089112581412879D-02  /
   data green(25,15,12) /   3.1715424641210854D-02  /
   data green(25,15,13) /   3.1323731671770487D-02  /
   data green(25,15,14) /   3.0916628126535634D-02  /
   data green(25,15,15) /   3.0496615875603084D-02  /
   data green(25,16, 0) /   3.3690565898028185D-02  /
   data green(25,16, 1) /   3.3671429860046072D-02  /
   data green(25,16, 2) /   3.3614217827003572D-02  /
   data green(25,16, 3) /   3.3519512440430346D-02  /
   data green(25,16, 4) /   3.3388266414110399D-02  /
   data green(25,16, 5) /   3.3221776074538772D-02  /
   data green(25,16, 6) /   3.3021646362710044D-02  /
   data green(25,16, 7) /   3.2789749118668582D-02  /
   data green(25,16, 8) /   3.2528176705300130D-02  /
   data green(25,16, 9) /   3.2239193101709467D-02  /
   data green(25,16,10) /   3.1925184519849804D-02  /
   data green(25,16,11) /   3.1588611395772147D-02  /
   data green(25,16,12) /   3.1231963312933417D-02  /
   data green(25,16,13) /   3.0857718066899475D-02  /
   data green(25,16,14) /   3.0468305714412006D-02  /
   data green(25,16,15) /   3.0066078095815681D-02  /
   data green(25,16,16) /   2.9653284001464979D-02  /
   data green(25,17, 0) /   3.3076331622766464D-02  /
   data green(25,17, 1) /   3.3058224792835447D-02  /
   data green(25,17, 2) /   3.3004083085353235D-02  /
   data green(25,17, 3) /   3.2914437927484788D-02  /
   data green(25,17, 4) /   3.2790158898525761D-02  /
   data green(25,17, 5) /   3.2632430439941969D-02  /
   data green(25,17, 6) /   3.2442720984224377D-02  /
   data green(25,17, 7) /   3.2222746054505987D-02  /
   data green(25,17, 8) /   3.1974427096039940D-02  /
   data green(25,17, 9) /   3.1699847875901452D-02  /
   data green(25,17,10) /   3.1401210236005984D-02  /
   data green(25,17,11) /   3.1080790825897254D-02  /
   data green(25,17,12) /   3.0740900202713219D-02  /
   data green(25,17,13) /   3.0383845396666647D-02  /
   data green(25,17,14) /   3.0011896731010915D-02  /
   data green(25,17,15) /   2.9627259381792736D-02  /
   data green(25,17,16) /   2.9232049884993224D-02  /
   data green(25,17,17) /   2.8828277560710082D-02  /
   data green(25,18, 0) /   3.2460317372258921D-02  /
   data green(25,18, 1) /   3.2443204724019056D-02  /
   data green(25,18, 2) /   3.2392029471520936D-02  /
   data green(25,18, 3) /   3.2307275381739241D-02  /
   data green(25,18, 4) /   3.2189734592051482D-02  /
   data green(25,18, 5) /   3.2040487165912249D-02  /
   data green(25,18, 6) /   3.1860873937306210D-02  /
   data green(25,18, 7) /   3.1652463961647496D-02  /
   data green(25,18, 8) /   3.1417018075763893D-02  /
   data green(25,18, 9) /   3.1156450143866036D-02  /
   data green(25,18,10) /   3.0872787534841316D-02  /
   data green(25,18,11) /   3.0568132253409389D-02  /
   data green(25,18,12) /   3.0244623954855700D-02  /
   data green(25,18,13) /   2.9904405834672534D-02  /
   data green(25,18,14) /   2.9549594124921589D-02  /
   data green(25,18,15) /   2.9182251670328518D-02  /
   data green(25,18,16) /   2.8804365816564464D-02  /
   data green(25,18,17) /   2.8417830633286446D-02  /
   data green(25,18,18) /   2.8024433322588006D-02  /
   data green(25,19, 0) /   3.1845199708182323D-02  /
   data green(25,19, 1) /   3.1829042579549424D-02  /
   data green(25,19, 2) /   3.1780719004832758D-02  /
   data green(25,19, 3) /   3.1700668650665544D-02  /
   data green(25,19, 4) /   3.1589611927301878D-02  /
   data green(25,19, 5) /   3.1448532082856806D-02  /
   data green(25,19, 6) /   3.1278651371427955D-02  /
   data green(25,19, 7) /   3.1081402410537910D-02  /
   data green(25,19, 8) /   3.0858396006029701D-02  /
   data green(25,19, 9) /   3.0611386794054424D-02  /
   data green(25,19,10) /   3.0342238033182301D-02  /
   data green(25,19,11) /   3.0052886785919233D-02  /
   data green(25,19,12) /   2.9745310574606575D-02  /
   data green(25,19,13) /   2.9421496401372097D-02  /
   data green(25,19,14) /   2.9083412805424789D-02  /
   data green(25,19,15) /   2.8732985411763425D-02  /
   data green(25,19,16) /   2.8372076218320396D-02  /
   data green(25,19,17) /   2.8002466684734766D-02  /
   data green(25,19,18) /   2.7625844532270377D-02  /
   data green(25,19,19) /   2.7243794044011244D-02  /
   data green(25,20, 0) /   3.1233314136686834D-02  /
   data green(25,20, 1) /   3.1218071360374140D-02  /
   data green(25,20, 2) /   3.1172477147121441D-02  /
   data green(25,20, 3) /   3.1096930557374183D-02  /
   data green(25,20, 4) /   3.0992085894957201D-02  /
   data green(25,20, 5) /   3.0858837053588748D-02  /
   data green(25,20, 6) /   3.0698296641627203D-02  /
   data green(25,20, 7) /   3.0511770826992614D-02  /
   data green(25,20, 8) /   3.0300730986613185D-02  /
   data green(25,20, 9) /   3.0066783312297056D-02  /
   data green(25,20,10) /   2.9811637519382397D-02  /
   data green(25,20,11) /   2.9537075734049373D-02  /
   data green(25,20,12) /   2.9244922512654325D-02  /
   data green(25,20,13) /   2.8937016787405328D-02  /
   data green(25,20,14) /   2.8615186353308474D-02  /
   data green(25,20,15) /   2.8281225326661850D-02  /
   data green(25,20,16) /   2.7936874828217110D-02  /
   data green(25,20,17) /   2.7583806984171129D-02  /
   data green(25,20,18) /   2.7223612201824014D-02  /
   data green(25,20,19) /   2.6857789567381093D-02  /
   data green(25,20,20) /   2.6487740131672322D-02  /
   data green(25,21, 0) /   3.0626673914102735D-02  /
   data green(25,21, 1) /   3.0612302788237900D-02  /
   data green(25,21, 2) /   3.0569310974829757D-02  /
   data green(25,21, 3) /   3.0498060302900691D-02  /
   data green(25,21, 4) /   3.0399144408371482D-02  /
   data green(25,21, 5) /   3.0273375069533683D-02  /
   data green(25,21, 6) /   3.0121763949437823D-02  /
   data green(25,21, 7) /   2.9945500539066430D-02  /
   data green(25,21, 8) /   2.9745927219287996D-02  /
   data green(25,21, 9) /   2.9524512422427700D-02  /
   data green(25,21,10) /   2.9282822876684451D-02  /
   data green(25,21,11) /   2.9022495864593174D-02  /
   data green(25,21,12) /   2.8745212330197677D-02  /
   data green(25,21,13) /   2.8452671540838092D-02  /
   data green(25,21,14) /   2.8146567861503820D-02  /
   data green(25,21,15) /   2.7828570044973421D-02  /
   data green(25,21,16) /   2.7500303290205790D-02  /
   data green(25,21,17) /   2.7163334183094901D-02  /
   data green(25,21,18) /   2.6819158513614116D-02  /
   data green(25,21,19) /   2.6469191864903888D-02  /
   data green(25,21,20) /   2.6114762794142182D-02  /
   data green(25,21,21) /   2.5757108371442631D-02  /
   data green(25,22, 0) /   3.0026992531257315D-02  /
   data green(25,22, 1) /   3.0013449635637367D-02  /
   data green(25,22, 2) /   2.9972931051076963D-02  /
   data green(25,22, 3) /   2.9905764592162480D-02  /
   data green(25,22, 4) /   2.9812488419776646D-02  /
   data green(25,22, 5) /   2.9693839126063722D-02  /
   data green(25,22, 6) /   2.9550735784988198D-02  /
   data green(25,22, 7) /   2.9384260636539653D-02  /
   data green(25,22, 8) /   2.9195637180433203D-02  /
   data green(25,22, 9) /   2.8986206512885561D-02  /
   data green(25,22,10) /   2.8757402747926895D-02  /
   data green(25,22,11) /   2.8510728327117294D-02  /
   data green(25,22,12) /   2.8247729946080317D-02  /
   data green(25,22,13) /   2.7969975722589763D-02  /
   data green(25,22,14) /   2.7679034109509327D-02  /
   data green(25,22,15) /   2.7376454926842691D-02  /
   data green(25,22,16) /   2.7063752759499897D-02  /
   data green(25,22,17) /   2.6742392848380106D-02  /
   data green(25,22,18) /   2.6413779497266799D-02  /
   data green(25,22,19) /   2.6079246930069293D-02  /
   data green(25,22,20) /   2.5740052463537628D-02  /
   data green(25,22,21) /   2.5397371809632741D-02  /
   data green(25,22,22) /   2.5052296288003893D-02  /
   data green(25,23, 0) /   2.9435708355869807D-02  /
   data green(25,23, 1) /   2.9422950224265342D-02  /
   data green(25,23, 2) /   2.9384775496601331D-02  /
   data green(25,23, 3) /   2.9321481006694065D-02  /
   data green(25,23, 4) /   2.9233554344824070D-02  /
   data green(25,23, 5) /   2.9121663476460214D-02  /
   data green(25,23, 6) /   2.8986642821178726D-02  /
   data green(25,23, 7) /   2.8829476353318166D-02  /
   data green(25,23, 8) /   2.8651278379235402D-02  /
   data green(25,23, 9) /   2.8453272698519612D-02  /
   data green(25,23,10) /   2.8236770867990684D-02  /
   data green(25,23,11) /   2.8003150260896684D-02  /
   data green(25,23,12) /   2.7753832555227426D-02  /
   data green(25,23,13) /   2.7490263202022733D-02  /
   data green(25,23,14) /   2.7213892325327178D-02  /
   data green(25,23,15) /   2.6926157398290535D-02  /
   data green(25,23,16) /   2.6628467932317591D-02  /
   data green(25,23,17) /   2.6322192314320841D-02  /
   data green(25,23,18) /   2.6008646835677197D-02  /
   data green(25,23,19) /   2.5689086878479483D-02  /
   data green(25,23,20) /   2.5364700161661188D-02  /
   data green(25,23,21) /   2.5036601901835463D-02  /
   data green(25,23,22) /   2.4705831710508549D-02  /
   data green(25,23,23) /   2.4373352029275090D-02  /
   data green(25,24, 0) /   2.8854010210580247D-02  /
   data green(25,24, 1) /   2.8841993872457963D-02  /
   data green(25,24, 2) /   2.8806035049002317D-02  /
   data green(25,24, 3) /   2.8746402428531530D-02  /
   data green(25,24, 4) /   2.8663537620926959D-02  /
   data green(25,24, 5) /   2.8558046117290362D-02  /
   data green(25,24, 6) /   2.8430685146272245D-02  /
   data green(25,24, 7) /   2.8282348898721239D-02  /
   data green(25,24, 8) /   2.8114051672882879D-02  /
   data green(25,24, 9) /   2.7926909539700222D-02  /
   data green(25,24,10) /   2.7722121141389944D-02  /
   data green(25,24,11) /   2.7500948218608354D-02  /
   data green(25,24,12) /   2.7264696416568222D-02  /
   data green(25,24,13) /   2.7014696854301280D-02  /
   data green(25,24,14) /   2.6752288860517528D-02  /
   data green(25,24,15) /   2.6478804190891724D-02  /
   data green(25,24,16) /   2.6195552951290835D-02  /
   data green(25,24,17) /   2.5903811364687827D-02  /
   data green(25,24,18) /   2.5604811440356430D-02  /
   data green(25,24,19) /   2.5299732535206159D-02  /
   data green(25,24,20) /   2.4989694740395724D-02  /
   data green(25,24,21) /   2.4675753982178979D-02  /
   data green(25,24,22) /   2.4358898693963089D-02  /
   data green(25,24,23) /   2.4040047895822592D-02  /
   data green(25,24,24) /   2.3720050506831687D-02  /
   data green(25,25, 0) /   2.8282862939895504D-02  /
   data green(25,25, 1) /   2.8271546346389360D-02  /
   data green(25,25, 2) /   2.8237678169478343D-02  /
   data green(25,25, 3) /   2.8181501582085332D-02  /
   data green(25,25, 4) /   2.8103416477371326D-02  /
   data green(25,25, 5) /   2.8003971598112187D-02  /
   data green(25,25, 6) /   2.7883853946896922D-02  /
   data green(25,25, 7) /   2.7743875873120806D-02  /
   data green(25,25, 8) /   2.7584960302167175D-02  /
   data green(25,25, 9) /   2.7408124614533437D-02  /
   data green(25,25,10) /   2.7214463697365326D-02  /
   data green(25,25,11) /   2.7005132679462807D-02  /
   data green(25,25,12) /   2.6781329826608913D-02  /
   data green(25,25,13) /   2.6544280021617343D-02  /
   data green(25,25,14) /   2.6295219188049831D-02  /
   data green(25,25,15) /   2.6035379943524315D-02  /
   data green(25,25,16) /   2.5765978693008521D-02  /
   data green(25,25,17) /   2.5488204298881146D-02  /
   data green(25,25,18) /   2.5203208396365750D-02  /
   data green(25,25,19) /   2.4912097362737484D-02  /
   data green(25,25,20) /   2.4615925898026778D-02  /
   data green(25,25,21) /   2.4315692134482505D-02  /
   data green(25,25,22) /   2.4012334161755203D-02  /
   data green(25,25,23) /   2.3706727834007058D-02  /
   data green(25,25,24) /   2.3399685712941570D-02  /
   data green(25,25,25) /   2.3091956995837153D-02  /
   data green(26, 0, 0) /   3.8475823111217562D-02  /
   data green(26, 1, 0) /   3.8447269035176722D-02  /
   data green(26, 1, 1) /   3.8418778914845284D-02  /
   data green(26, 2, 0) /   3.8361990545311545D-02  /
   data green(26, 2, 1) /   3.8333690846489674D-02  /
   data green(26, 2, 2) /   3.8249169752120983D-02  /
   data green(26, 3, 0) /   3.8221124399639546D-02  /
   data green(26, 3, 1) /   3.8193137328464077D-02  /
   data green(26, 3, 2) /   3.8109547117281387D-02  /
   data green(26, 3, 3) /   3.7971452992000873D-02  /
   data green(26, 4, 0) /   3.8026518086475454D-02  /
   data green(26, 4, 1) /   3.7998958999415119D-02  /
   data green(26, 4, 2) /   3.7916643247426128D-02  /
   data green(26, 4, 3) /   3.7780642074014055D-02  /
   data green(26, 4, 4) /   3.7592697572526056D-02  /
   data green(26, 5, 0) /   3.7780662808961335D-02  /
   data green(26, 5, 1) /   3.7753637972266656D-02  /
   data green(26, 5, 2) /   3.7672913266479237D-02  /
   data green(26, 5, 3) /   3.7539525425467443D-02  /
   data green(26, 5, 4) /   3.7355160986247180D-02  /
   data green(26, 5, 5) /   3.7122097337444163D-02  /
   data green(26, 6, 0) /   3.7486606495045378D-02  /
   data green(26, 6, 1) /   3.7460211287004622D-02  /
   data green(26, 6, 2) /   3.7381361880428066D-02  /
   data green(26, 6, 3) /   3.7251054940293980D-02  /
   data green(26, 6, 4) /   3.7070912462496927D-02  /
   data green(26, 6, 5) /   3.6843125967739962D-02  /
   data green(26, 6, 6) /   3.6570383755540503D-02  /
   data green(26, 7, 0) /   3.7147852673284229D-02  /
   data green(26, 7, 1) /   3.7122170284562779D-02  /
   data green(26, 7, 2) /   3.7045444226347188D-02  /
   data green(26, 7, 3) /   3.6918626589605125D-02  /
   data green(26, 7, 4) /   3.6743267527416155D-02  /
   data green(26, 7, 5) /   3.6521462889096419D-02  /
   data green(26, 7, 6) /   3.6255785863210259D-02  /
   data green(26, 7, 7) /   3.5949206918735876D-02  /
   data green(26, 8, 0) /   3.6768251510172163D-02  /
   data green(26, 8, 1) /   3.6743352158445299D-02  /
   data green(26, 8, 2) /   3.6668958948076105D-02  /
   data green(26, 8, 3) /   3.6545975976121463D-02  /
   data green(26, 8, 4) /   3.6375875986606639D-02  /
   data green(26, 8, 5) /   3.6160651637851669D-02  /
   data green(26, 8, 6) /   3.5902751790119848D-02  /
   data green(26, 8, 7) /   3.5605006733576972D-02  /
   data green(26, 8, 8) /   3.5270546658421158D-02  /
   data green(26, 9, 0) /   3.6351889148275203D-02  /
   data green(26, 9, 1) /   3.6327829787794183D-02  /
   data green(26, 9, 2) /   3.6255939497051567D-02  /
   data green(26, 9, 3) /   3.6137072025686362D-02  /
   data green(26, 9, 4) /   3.5972618850876451D-02  /
   data green(26, 9, 5) /   3.5764464174717431D-02  /
   data green(26, 9, 6) /   3.5514925990843896D-02  /
   data green(26, 9, 7) /   3.5226686768223048D-02  /
   data green(26, 9, 8) /   3.4902717661624456D-02  /
   data green(26, 9, 9) /   3.4546200172438861D-02  /
   data green(26,10, 0) /   3.5902980819469806D-02  /
   data green(26,10, 1) /   3.5879805297432288D-02  /
   data green(26,10, 2) /   3.5810549018590543D-02  /
   data green(26,10, 3) /   3.5696014037595841D-02  /
   data green(26,10, 4) /   3.5537508321926328D-02  /
   data green(26,10, 5) /   3.5336804489792900D-02  /
   data green(26,10, 6) /   3.5096085680378107D-02  /
   data green(26,10, 7) /   3.4817881738866116D-02  /
   data green(26,10, 8) /   3.4504999238163872D-02  /
   data green(26,10, 9) /   3.4160448891884444D-02  /
   data green(26,10,10) /   3.3787373672885319D-02  /
   data green(26,11, 0) /   3.5425772187559240D-02  /
   data green(26,11, 1) /   3.5403511779001098D-02  /
   data green(26,11, 2) /   3.5336983197868765D-02  /
   data green(26,11, 3) /   3.5226936370138559D-02  /
   data green(26,11, 4) /   3.5074594981368827D-02  /
   data green(26,11, 5) /   3.4881618891620725D-02  /
   data green(26,11, 6) /   3.4650054739298962D-02  /
   data green(26,11, 7) /   3.4382277565059648D-02  /
   data green(26,11, 8) /   3.4080926602897237D-02  /
   data green(26,11, 9) /   3.3748838433240370D-02  /
   data green(26,11,10) /   3.3388980499286972D-02  /
   data green(26,11,11) /   3.3004387603583221D-02  /
   data green(26,12, 0) /   3.4924452169088606D-02  /
   data green(26,12, 1) /   3.4903126410836441D-02  /
   data green(26,12, 2) /   3.4839384269845945D-02  /
   data green(26,12, 3) /   3.4733923910428809D-02  /
   data green(26,12, 4) /   3.4587885256001823D-02  /
   data green(26,12, 5) /   3.4402815957238819D-02  /
   data green(26,12, 6) /   3.4180626581682409D-02  /
   data green(26,12, 7) /   3.3923537524520174D-02  /
   data green(26,12, 8) /   3.3634020431912104D-02  /
   data green(26,12, 9) /   3.3314736987471720D-02  /
   data green(26,12,10) /   3.2968477760262349D-02  /
   data green(26,12,11) /   3.2598103490609182D-02  /
   data green(26,12,12) /   3.2206490751587519D-02  /
   data green(26,13, 0) /   3.4403079246865274D-02  /
   data green(26,13, 1) /   3.4382696988088685D-02  /
   data green(26,13, 2) /   3.4321768195422271D-02  /
   data green(26,13, 3) /   3.4220940314761333D-02  /
   data green(26,13, 4) /   3.4081271121798241D-02  /
   data green(26,13, 5) /   3.3904198073253573D-02  /
   data green(26,13, 6) /   3.3691497873861130D-02  /
   data green(26,13, 7) /   3.3445238451167210D-02  /
   data green(26,13, 8) /   3.3167725797324132D-02  /
   data green(26,13, 9) /   3.2861448204373700D-02  /
   data green(26,13,10) /   3.2529020303051707D-02  /
   data green(26,13,11) /   3.2173129048665108D-02  /
   data green(26,13,12) /   3.1796483425490810D-02  /
   data green(26,13,13) /   3.1401769211025790D-02  /
   data green(26,14, 0) /   3.3865522151552384D-02  /
   data green(26,14, 1) /   3.3846082742076594D-02  /
   data green(26,14, 2) /   3.3787965893217392D-02  /
   data green(26,14, 3) /   3.3691769924723843D-02  /
   data green(26,14, 4) /   3.3558472971100754D-02  /
   data green(26,14, 5) /   3.3389405512601535D-02  /
   data green(26,14, 6) /   3.3186214067161240D-02  /
   data green(26,14, 7) /   3.2950817953807977D-02  /
   data green(26,14, 8) /   3.2685361280872065D-02  /
   data green(26,14, 9) /   3.2392162384633258D-02  /
   data green(26,14,10) /   3.2073662857840808D-02  /
   data green(26,14,11) /   3.1732378089855143D-02  /
   data green(26,14,12) /   3.1370850927478135D-02  /
   data green(26,14,13) /   3.0991609697749559D-02  /
   data green(26,14,14) /   3.0597131449059026D-02  /
   data green(26,15, 0) /   3.3315414825505020D-02  /
   data green(26,15, 1) /   3.3296909371397454D-02  /
   data green(26,15, 2) /   3.3241578473951885D-02  /
   data green(26,15, 3) /   3.3149973344882518D-02  /
   data green(26,15, 4) /   3.3022995679369074D-02  /
   data green(26,15, 5) /   3.2861873142222012D-02  /
   data green(26,15, 6) /   3.2668126903341879D-02  /
   data green(26,15, 7) /   3.2443532878587991D-02  /
   data green(26,15, 8) /   3.2190078551267605D-02  /
   data green(26,15, 9) /   3.1909917324102165D-02  /
   data green(26,15,10) /   3.1605322290361809D-02  /
   data green(26,15,11) /   3.1278641137293985D-02  /
   data green(26,15,12) /   3.0932253634700646D-02  /
   data green(26,15,13) /   3.0568532849651209D-02  /
   data green(26,15,14) /   3.0189810896980886D-02  /
   data green(26,15,15) /   2.9798349712189153D-02  /
   data green(26,16, 0) /   3.2756124841644731D-02  /
   data green(26,16, 1) /   3.2738537468778181D-02  /
   data green(26,16, 2) /   3.2685945694354068D-02  /
   data green(26,16, 3) /   3.2598855951808461D-02  /
   data green(26,16, 4) /   3.2478097211122457D-02  /
   data green(26,16, 5) /   3.2324799184518693D-02  /
   data green(26,16, 6) /   3.2140363408075044D-02  /
   data green(26,16, 7) /   3.1926428627933386D-02  /
   data green(26,16, 8) /   3.1684832116210093D-02  /
   data green(26,16, 9) /   3.1417568616456011D-02  /
   data green(26,16,10) /   3.1126748577627799D-02  /
   data green(26,16,11) /   3.0814557195771176D-02  /
   data green(26,16,12) /   3.0483215567924146D-02  /
   data green(26,16,13) /   3.0134945000380262D-02  /
   data green(26,16,14) /   2.9771935230349489D-02  /
   data green(26,16,15) /   2.9396317040019323D-02  /
   data green(26,16,16) /   2.9010139483879059D-02  /
   data green(26,17, 0) /   3.2190733936076521D-02  /
   data green(26,17, 1) /   3.2174043014745218D-02  /
   data green(26,17, 2) /   3.2124126338192442D-02  /
   data green(26,17, 3) /   3.2041448100978941D-02  /
   data green(26,17, 4) /   3.1926768610208406D-02  /
   data green(26,17, 5) /   3.1781124972468745D-02  /
   data green(26,17, 6) /   3.1605805419112919D-02  /
   data green(26,17, 7) /   3.1402318497181017D-02  /
   data green(26,17, 8) /   3.1172358528311524D-02  /
   data green(26,17, 9) /   3.0917768811065772D-02  /
   data green(26,17,10) /   3.0640504017265008D-02  /
   data green(26,17,11) /   3.0342593123178708D-02  /
   data green(26,17,12) /   3.0026104040845387D-02  /
   data green(26,17,13) /   2.9693110895712582D-02  /
   data green(26,17,14) /   2.9345664656567984D-02  /
   data green(26,17,15) /   2.8985767582620510D-02  /
   data green(26,17,16) /   2.8615351726990812D-02  /
   data green(26,17,17) /   2.8236261537593812D-02  /
   data green(26,18, 0) /   3.1622029008820551D-02  /
   data green(26,18, 1) /   3.1606208304746738D-02  /
   data green(26,18, 2) /   3.1558888927338768D-02  /
   data green(26,18, 3) /   3.1480495493230294D-02  /
   data green(26,18, 4) /   3.1371723914952236D-02  /
   data green(26,18, 5) /   3.1233524334792179D-02  /
   data green(26,18, 6) /   3.1067078394757947D-02  /
   data green(26,18, 7) /   3.0873771893128260D-02  /
   data green(26,18, 8) /   3.0655164033213594D-02  /
   data green(26,18, 9) /   3.0412954540033323D-02  /
   data green(26,18,10) /   3.0148949908218467D-02  /
   data green(26,18,11) /   2.9865029959485076D-02  /
   data green(26,18,12) /   2.9563115745679453D-02  /
   data green(26,18,13) /   2.9245139651708989D-02  /
   data green(26,18,14) /   2.8913018350165541D-02  /
   data green(26,18,15) /   2.8568629053164738D-02  /
   data green(26,18,16) /   2.8213789311037980D-02  /
   data green(26,18,17) /   2.7850240432576171D-02  /
   data green(26,18,18) /   2.7479634454277771D-02  /
   data green(26,19, 0) /   3.1052501819752162D-02  /
   data green(26,19, 1) /   3.1037521548071014D-02  /
   data green(26,19, 2) /   3.0992711034542153D-02  /
   data green(26,19, 3) /   3.0918458026657301D-02  /
   data green(26,19, 4) /   3.0815398397234527D-02  /
   data green(26,19, 5) /   3.0684401099767042D-02  /
   data green(26,19, 6) /   3.0526548094119033D-02  /
   data green(26,19, 7) /   3.0343110138922082D-02  /
   data green(26,19, 8) /   3.0135519483845956D-02  /
   data green(26,19, 9) /   2.9905340561073621D-02  /
   data green(26,19,10) /   2.9654239772276457D-02  /
   data green(26,19,11) /   2.9383955402658795D-02  /
   data green(26,19,12) /   2.9096268579179883D-02  /
   data green(26,19,13) /   2.8792976040384510D-02  /
   data green(26,19,14) /   2.8475865315581550D-02  /
   data green(26,19,15) /   2.8146692735666664D-02  /
   data green(26,19,16) /   2.7807164528828213D-02  /
   data green(26,19,17) /   2.7458921101023905D-02  /
   data green(26,19,18) /   2.7103524469728045D-02  /
   data green(26,19,19) /   2.6742448713374840D-02  /
   data green(26,20, 0) /   3.0484355617242862D-02  /
   data green(26,20, 1) /   3.0470183385753725D-02  /
   data green(26,20, 2) /   3.0427785473905792D-02  /
   data green(26,20, 3) /   3.0357515456945833D-02  /
   data green(26,20, 4) /   3.0259953512220744D-02  /
   data green(26,20, 5) /   3.0135893182950610D-02  /
   data green(26,20, 6) /   2.9986323685412979D-02  /
   data green(26,20, 7) /   2.9812408522384551D-02  /
   data green(26,20, 8) /   2.9615461285821133D-02  /
   data green(26,20, 9) /   2.9396919593384391D-02  /
   data green(26,20,10) /   2.9158318107228194D-02  /
   data green(26,20,11) /   2.8901261535043580D-02  /
   data green(26,20,12) /   2.8627398422076771D-02  /
   data green(26,20,13) /   2.8338396420294558D-02  /
   data green(26,20,14) /   2.8035919579475738D-02  /
   data green(26,20,15) /   2.7721608056609365D-02  /
   data green(26,20,16) /   2.7397060494898880D-02  /
   data green(26,20,17) /   2.7063819190051890D-02  /
   data green(26,20,18) /   2.6723358045092339D-02  /
   data green(26,20,19) /   2.6377073218994247D-02  /
   data green(26,20,20) /   2.6026276300218178D-02  /
   data green(26,21, 0) /   2.9919517045436821D-02  /
   data green(26,21, 1) /   2.9906118681325932D-02  /
   data green(26,21, 2) /   2.9866031746459975D-02  /
   data green(26,21, 3) /   2.9799578282087090D-02  /
   data green(26,21, 4) /   2.9707287027255200D-02  /
   data green(26,21, 5) /   2.9589881792532854D-02  /
   data green(26,21, 6) /   2.9448265892229149D-02  /
   data green(26,21, 7) /   2.9283503282942467D-02  /
   data green(26,21, 8) /   2.9096797161342707D-02  /
   data green(26,21, 9) /   2.8889466830864366D-02  /
   data green(26,21,10) /   2.8662923655592344D-02  /
   data green(26,21,11) /   2.8418646884199320D-02  /
   data green(26,21,12) /   2.8158160054595972D-02  /
   data green(26,21,13) /   2.7883008590227142D-02  /
   data green(26,21,14) /   2.7594739081743910D-02  /
   data green(26,21,15) /   2.7294880622896069D-02  /
   data green(26,21,16) /   2.6984928445630153D-02  /
   data green(26,21,17) /   2.6666329983603510D-02  /
   data green(26,21,18) /   2.6340473390774447D-02  /
   data green(26,21,19) /   2.6008678455612253D-02  /
   data green(26,21,20) /   2.5672189783233771D-02  /
   data green(26,21,21) /   2.5332172067391592D-02  /
   data green(26,22, 0) /   2.9359651846878811D-02  /
   data green(26,22, 1) /   2.9346992107235466D-02  /
   data green(26,22, 2) /   2.9309111281893944D-02  /
   data green(26,22, 3) /   2.9246302422508714D-02  /
   data green(26,22, 4) /   2.9159046942057280D-02  /
   data green(26,22, 5) /   2.9048004415837481D-02  /
   data green(26,22, 6) /   2.8913998901842996D-02  /
   data green(26,22, 7) /   2.8758002329721873D-02  /
   data green(26,22, 8) /   2.8581115599112861D-02  /
   data green(26,22, 9) /   2.8384548079963566D-02  /
   data green(26,22,10) /   2.8169596219228581D-02  /
   data green(26,22,11) /   2.7937621933138687D-02  /
   data green(26,22,12) /   2.7690031407622846D-02  /
   data green(26,22,13) /   2.7428254848761876D-02  /
   data green(26,22,14) /   2.7153727628476954D-02  /
   data green(26,22,15) /   2.6867873166048183D-02  /
   data green(26,22,16) /   2.6572087780825773D-02  /
   data green(26,22,17) /   2.6267727651699407D-02  /
   data green(26,22,18) /   2.5956097929101491D-02  /
   data green(26,22,19) /   2.5638443968553869D-02  /
   data green(26,22,20) /   2.5315944592586133D-02  /
   data green(26,22,21) /   2.4989707240574480D-02  /
   data green(26,22,22) /   2.4660764833003619D-02  /
   data green(26,23, 0) /   2.8806183080752411D-02  /
   data green(26,23, 1) /   2.8794226251653888D-02  /
   data green(26,23, 2) /   2.8758445214460963D-02  /
   data green(26,23, 3) /   2.8699106456613912D-02  /
   data green(26,23, 4) /   2.8616647990559178D-02  /
   data green(26,23, 5) /   2.8511670415585728D-02  /
   data green(26,23, 6) /   2.8384924909801063D-02  /
   data green(26,23, 7) /   2.8237298617215141D-02  /
   data green(26,23, 8) /   2.8069797974503576D-02  /
   data green(26,23, 9) /   2.7883530568921230D-02  /
   data green(26,23,10) /   2.7679686132573988D-02  /
   data green(26,23,11) /   2.7459517260984535D-02  /
   data green(26,23,12) /   2.7224320399902038D-02  /
   data green(26,23,13) /   2.6975417579361994D-02  /
   data green(26,23,14) /   2.6714139294619754D-02  /
   data green(26,23,15) /   2.6441808846338573D-02  /
   data green(26,23,16) /   2.6159728363398405D-02  /
   data green(26,23,17) /   2.5869166646076174D-02  /
   data green(26,23,18) /   2.5571348889150615D-02  /
   data green(26,23,19) /   2.5267448276488060D-02  /
   data green(26,23,20) /   2.4958579382455673D-02  /
   data green(26,23,21) /   2.4645793271638090D-02  /
   data green(26,23,22) /   2.4330074156485880D-02  /
   data green(26,23,23) /   2.4012337451773016D-02  /
   data green(26,24, 0) /   2.8260310790579887D-02  /
   data green(26,24, 1) /   2.8249021182527909D-02  /
   data green(26,24, 2) /   2.8215233638928857D-02  /
   data green(26,24, 3) /   2.8159190372589602D-02  /
   data green(26,24, 4) /   2.8081289706140291D-02  /
   data green(26,24, 5) /   2.7982078244345381D-02  /
   data green(26,24, 6) /   2.7862240341417752D-02  /
   data green(26,24, 7) /   2.7722585256596939D-02  /
   data green(26,24, 8) /   2.7564032460263337D-02  /
   data green(26,24, 9) /   2.7387595595027986D-02  /
   data green(26,24,10) /   2.7194365610969936D-02  /
   data green(26,24,11) /   2.6985493583002469D-02  /
   data green(26,24,12) /   2.6762173684509265D-02  /
   data green(26,24,13) /   2.6525626739412840D-02  /
   data green(26,24,14) /   2.6277084709935315D-02  /
   data green(26,24,15) /   2.6017776404839228D-02  /
   data green(26,24,16) /   2.5748914617945887D-02  /
   data green(26,24,17) /   2.5471684833597099D-02  /
   data green(26,24,18) /   2.5187235567959158D-02  /
   data green(26,24,19) /   2.4896670355184659D-02  /
   data green(26,24,20) /   2.4601041337007398D-02  /
   data green(26,24,21) /   2.4301344374035257D-02  /
   data green(26,24,22) /   2.3998515566780631D-02  /
   data green(26,24,23) /   2.3693429053726938D-02  /
   data green(26,24,24) /   2.3386895941480083D-02  /
   data green(26,25, 0) /   2.7723032262716881D-02  /
   data green(26,25, 1) /   2.7712374612023170D-02  /
   data green(26,25, 2) /   2.7680475495195402D-02  /
   data green(26,25, 3) /   2.7627554994299638D-02  /
   data green(26,25, 4) /   2.7553975221088257D-02  /
   data green(26,25, 5) /   2.7460233464912047D-02  /
   data green(26,25, 6) /   2.7346952958801986D-02  /
   data green(26,25, 7) /   2.7214871596109056D-02  /
   data green(26,25, 8) /   2.7064828989791522D-02  /
   data green(26,25, 9) /   2.6897752304133792D-02  /
   data green(26,25,10) /   2.6714641303707883D-02  /
   data green(26,25,11) /   2.6516553057792977D-02  /
   data green(26,25,12) /   2.6304586712718413D-02  /
   data green(26,25,13) /   2.6079868703221348D-02  /
   data green(26,25,14) /   2.5843538721061144D-02  /
   data green(26,25,15) /   2.5596736699152002D-02  /
   data green(26,25,16) /   2.5340591006531986D-02  /
   data green(26,25,17) /   2.5076207987268468D-02  /
   data green(26,25,18) /   2.4804662917931489D-02  /
   data green(26,25,19) /   2.4526992405809984D-02  /
   data green(26,25,20) /   2.4244188205074617D-02  /
   data green(26,25,21) /   2.3957192391377516D-02  /
   data green(26,25,22) /   2.3666893807057714D-02  /
   data green(26,25,23) /   2.3374125668850221D-02  /
   data green(26,25,24) /   2.3079664217070847D-02  /
   data green(26,25,25) /   2.2784228278732752D-02  /
   data green(26,26, 0) /   2.7195162207886057D-02  /
   data green(26,26, 1) /   2.7185101994485178D-02  /
   data green(26,26, 2) /   2.7154988417377755D-02  /
   data green(26,26, 3) /   2.7105021421689982D-02  /
   data green(26,26, 4) /   2.7035530147461420D-02  /
   data green(26,26, 5) /   2.6946966932907809D-02  /
   data green(26,26, 6) /   2.6839899221399139D-02  /
   data green(26,26, 7) /   2.6714999652982140D-02  /
   data green(26,26, 8) /   2.6573034672834814D-02  /
   data green(26,26, 9) /   2.6414852022567721D-02  /
   data green(26,26,10) /   2.6241367495112369D-02  /
   data green(26,26,11) /   2.6053551330748333D-02  /
   data green(26,26,12) /   2.5852414612479988D-02  /
   data green(26,26,13) /   2.5638995986221061D-02  /
   data green(26,26,14) /   2.5414348988378358D-02  /
   data green(26,26,15) /   2.5179530213953104D-02  /
   data green(26,26,16) /   2.4935588505623971D-02  /
   data green(26,26,17) /   2.4683555291533984D-02  /
   data green(26,26,18) /   2.4424436149263689D-02  /
   data green(26,26,19) /   2.4159203627738643D-02  /
   data green(26,26,20) /   2.3888791318964459D-02  /
   data green(26,26,21) /   2.3614089138298901D-02  /
   data green(26,26,22) /   2.3335939745725402D-02  /
   data green(26,26,23) /   2.3055136021131534D-02  /
   data green(26,26,24) /   2.2772419493440537D-02  /
   data green(26,26,25) /   2.2488479615899608D-02  /
   data green(26,26,26) /   2.2203953777072540D-02  /
   data green(27, 0, 0) /   3.7049788507905627D-02  /
   data green(27, 1, 0) /   3.7024298005638452D-02  /
   data green(27, 1, 1) /   3.6998860409934130D-02  /
   data green(27, 2, 0) /   3.6948143941041560D-02  /
   data green(27, 2, 1) /   3.6922863957616694D-02  /
   data green(27, 2, 2) /   3.6847337052734130D-02  /
   data green(27, 3, 0) /   3.6822267566267726D-02  /
   data green(27, 3, 1) /   3.6797246633319319D-02  /
   data green(27, 3, 2) /   3.6722491504823436D-02  /
   data green(27, 3, 3) /   3.6598914597211762D-02  /
   data green(27, 4, 0) /   3.6648201486002910D-02  /
   data green(27, 4, 1) /   3.6623535773283832D-02  /
   data green(27, 4, 2) /   3.6549838998572830D-02  /
   data green(27, 4, 3) /   3.6428002012094646D-02  /
   data green(27, 4, 4) /   3.6259476122070942D-02  /
   data green(27, 5, 0) /   3.6428018026135017D-02  /
   data green(27, 5, 1) /   3.6403796683598874D-02  /
   data green(27, 5, 2) /   3.6331423982445030D-02  /
   data green(27, 5, 3) /   3.6211764095209271D-02  /
   data green(27, 5, 4) /   3.6046225273694998D-02  /
   data green(27, 5, 5) /   3.5836714010222960D-02  /
   data green(27, 6, 0) /   3.6164261829704773D-02  /
   data green(27, 6, 1) /   3.6140565554833859D-02  /
   data green(27, 6, 2) /   3.6069757527598757D-02  /
   data green(27, 6, 3) /   3.5952670830844552D-02  /
   data green(27, 6, 4) /   3.5790663515838032D-02  /
   data green(27, 6, 5) /   3.5585575042840566D-02  /
   data green(27, 6, 6) /   3.5339669209186686D-02  /
   data green(27, 7, 0) /   3.5859870903627390D-02  /
   data green(27, 7, 1) /   3.5836770870962618D-02  /
   data green(27, 7, 2) /   3.5767739810476984D-02  /
   data green(27, 7, 3) /   3.5653576071687794D-02  /
   data green(27, 7, 4) /   3.5495581591689009D-02  /
   data green(27, 7, 5) /   3.5295520842241281D-02  /
   data green(27, 7, 6) /   3.5055566977099752D-02  /
   data green(27, 7, 7) /   3.4778238347240809D-02  /
   data green(27, 8, 0) /   3.5518090733148378D-02  /
   data green(27, 8, 1) /   3.5495647902780945D-02  /
   data green(27, 8, 2) /   3.5428575730136747D-02  /
   data green(27, 8, 3) /   3.5317634987005866D-02  /
   data green(27, 8, 4) /   3.5164066864306599D-02  /
   data green(27, 8, 5) /   3.4969554582415567D-02  /
   data green(27, 8, 6) /   3.4736172965019400D-02  /
   data green(27, 8, 7) /   3.4466328888394210D-02  /
   data green(27, 8, 8) /   3.4162695838135824D-02  /
   data green(27, 9, 0) /   3.5142386050081513D-02  /
   data green(27, 9, 1) /   3.5120650846869418D-02  /
   data green(27, 9, 2) /   3.5055688147376432D-02  /
   data green(27, 9, 3) /   3.4948219096140577D-02  /
   data green(27, 9, 4) /   3.4799420782508320D-02  /
   data green(27, 9, 5) /   3.4610890603381676D-02  /
   data green(27, 9, 6) /   3.4384599385071171D-02  /
   data green(27, 9, 7) /   3.4122835916516679D-02  /
   data green(27, 9, 8) /   3.3828145848097166D-02  /
   data green(27, 9, 9) /   3.3503267963068802D-02  /
   data green(27,10, 0) /   3.4736354444125983D-02  /
   data green(27,10, 1) /   3.4715366778434877D-02  /
   data green(27,10, 2) /   3.4652632851856581D-02  /
   data green(27,10, 3) /   3.4548832892401266D-02  /
   data green(27,10, 4) /   3.4405077755598833D-02  /
   data green(27,10, 5) /   3.4222876067028528D-02  /
   data green(27,10, 6) /   3.4004090934721752D-02  /
   data green(27,10, 7) /   3.3750888626570358D-02  /
   data green(27,10, 8) /   3.3465681890506908D-02  /
   data green(27,10, 9) /   3.3151070656327934D-02  /
   data green(27,10,10) /   3.2809782717104587D-02  /
   data green(27,11, 0) /   3.4303645338943846D-02  /
   data green(27,11, 1) /   3.4283434924174926D-02  /
   data green(27,11, 2) /   3.4223018717207933D-02  /
   data green(27,11, 3) /   3.4123035441947711D-02  /
   data green(27,11, 4) /   3.3984528719076181D-02  /
   data green(27,11, 5) /   3.3808916961627238D-02  /
   data green(27,11, 6) /   3.3597953650973941D-02  /
   data green(27,11, 7) /   3.3353680139066699D-02  /
   data green(27,11, 8) /   3.3078373385293999D-02  /
   data green(27,11, 9) /   3.2774491103732803D-02  /
   data green(27,11,10) /   3.2444616684317792D-02  /
   data green(27,11,11) /   3.2091405992413072D-02  /
   data green(27,12, 0) /   3.3847887026795621D-02  /
   data green(27,12, 1) /   3.3828473939339090D-02  /
   data green(27,12, 2) /   3.3770435699606353D-02  /
   data green(27,12, 3) /   3.3674369565917676D-02  /
   data green(27,12, 4) /   3.3541251934717006D-02  /
   data green(27,12, 5) /   3.3372410917203207D-02  /
   data green(27,12, 6) /   3.3169490093360385D-02  /
   data green(27,12, 7) /   3.2934405351930089D-02  /
   data green(27,12, 8) /   3.2669296966922556D-02  /
   data green(27,12, 9) /   3.2376479132912368D-02  /
   data green(27,12,10) /   3.2058389094652778D-02  /
   data green(27,12,11) /   3.1717537788614206D-02  /
   data green(27,12,12) /   3.1356463601448052D-02  /
   data green(27,13, 0) /   3.3372623574647632D-02  /
   data green(27,13, 1) /   3.3354018997903093D-02  /
   data green(27,13, 2) /   3.3298392476678762D-02  /
   data green(27,13, 3) /   3.3206300382389670D-02  /
   data green(27,13, 4) /   3.3078652772932135D-02  /
   data green(27,13, 5) /   3.2916688542370752D-02  /
   data green(27,13, 6) /   3.2721942520427177D-02  /
   data green(27,13, 7) /   3.2496206208029954D-02  /
   data green(27,13, 8) /   3.2241484056767002D-02  /
   data green(27,13, 9) /   3.1959947274039439D-02  /
   data green(27,13,10) /   3.1653887071159764D-02  /
   data green(27,13,11) /   3.1325669090735100D-02  /
   data green(27,13,12) /   3.0977690482981997D-02  /
   data green(27,13,13) /   3.0612340782092951D-02  /
   data green(27,14, 0) /   3.2881262570843978D-02  /
   data green(27,14, 1) /   3.2863469665743023D-02  /
   data green(27,14, 2) /   3.2810264698768271D-02  /
   data green(27,14, 3) /   3.2722164183014678D-02  /
   data green(27,14, 4) /   3.2600013454654714D-02  /
   data green(27,14, 5) /   3.2444964260689567D-02  /
   data green(27,14, 6) /   3.2258445035621253D-02  /
   data green(27,14, 7) /   3.2042125347989177D-02  /
   data green(27,14, 8) /   3.1797876199187419D-02  /
   data green(27,14, 9) /   3.1527727931479184D-02  /
   data green(27,14,10) /   3.1233827456356739D-02  /
   data green(27,14,11) /   3.0918396366222047D-02  /
   data green(27,14,12) /   3.0583691267085635D-02  /
   data green(27,14,13) /   3.0231967396219705D-02  /
   data green(27,14,14) /   2.9865446294587282D-02  /
   data green(27,15, 0) /   3.2377033939929932D-02  /
   data green(27,15, 1) /   3.2360048788842029D-02  /
   data green(27,15, 2) /   3.2309254097579239D-02  /
   data green(27,15, 3) /   3.2225127909858417D-02  /
   data green(27,15, 4) /   3.2108453044793060D-02  /
   data green(27,15, 5) /   3.1960296969076117D-02  /
   data green(27,15, 6) /   3.1781985057458655D-02  /
   data green(27,15, 7) /   3.1575068535161420D-02  /
   data green(27,15, 8) /   3.1341288577873758D-02  /
   data green(27,15, 9) /   3.1082538118524339D-02  /
   data green(27,15,10) /   3.0800822879799229D-02  /
   data green(27,15,11) /   3.0498223031662317D-02  /
   data green(27,15,12) /   3.0176856684645627D-02  /
   data green(27,15,13) /   2.9838846196325375D-02  /
   data green(27,15,14) /   2.9486288014100817D-02  /
   data green(27,15,15) /   2.9121226523531381D-02  /
   data green(27,16, 0) /   3.1862959452702261D-02  /
   data green(27,16, 1) /   3.1846772030481550D-02  /
   data green(27,16, 2) /   3.1798358106389826D-02  /
   data green(27,16, 3) /   3.1718158919082397D-02  /
   data green(27,16, 4) /   3.1606897426719668D-02  /
   data green(27,16, 5) /   3.1465560298685080D-02  /
   data green(27,16, 6) /   3.1295373950833316D-02  /
   data green(27,16, 7) /   3.1097775748877555D-02  /
   data green(27,16, 8) /   3.0874381667991876D-02  /
   data green(27,16, 9) /   3.0626951768020671D-02  /
   data green(27,16,10) /   3.0357354826062895D-02  /
   data green(27,16,11) /   3.0067533372830170D-02  /
   data green(27,16,12) /   2.9759470222913940D-02  /
   data green(27,16,13) /   2.9435157391744574D-02  /
   data green(27,16,14) /   2.9096568073736413D-02  /
   data green(27,16,15) /   2.8745632135262003D-02  /
   data green(27,16,16) /   2.8384215367789071D-02  /
   data green(27,17, 0) /   3.1341832113736985D-02  /
   data green(27,17, 1) /   3.1326427247841911D-02  /
   data green(27,17, 2) /   3.1280349206384395D-02  /
   data green(27,17, 3) /   3.1204004284949204D-02  /
   data green(27,17, 4) /   3.1098058563202285D-02  /
   data green(27,17, 5) /   3.0963421846238303D-02  /
   data green(27,17, 6) /   3.0801226257970396D-02  /
   data green(27,17, 7) /   3.0612800460093872D-02  /
   data green(27,17, 8) /   3.0399640616158692D-02  /
   data green(27,17, 9) /   3.0163379288427241D-02  /
   data green(27,17,10) /   2.9905753447493348D-02  /
   data green(27,17,11) /   2.9628572699808377D-02  /
   data green(27,17,12) /   2.9333688709864119D-02  /
   data green(27,17,13) /   2.9022966628106276D-02  /
   data green(27,17,14) /   2.8698259149548884D-02  /
   data green(27,17,15) /   2.8361383637187355D-02  /
   data green(27,17,16) /   2.8014102561820291D-02  /
   data green(27,17,17) /   2.7658107345718750D-02  /
   data green(27,18, 0) /   3.0816204315342117D-02  /
   data green(27,18, 1) /   3.0801562605485756D-02  /
   data green(27,18, 2) /   3.0757762921800002D-02  /
   data green(27,18, 3) /   3.0685178607715229D-02  /
   data green(27,18, 4) /   3.0584422062196425D-02  /
   data green(27,18, 5) /   3.0456330460894781D-02  /
   data green(27,18, 6) /   3.0301946690088920D-02  /
   data green(27,18, 7) /   3.0122496333047411D-02  /
   data green(27,18, 8) /   2.9919361678198631D-02  /
   data green(27,18, 9) /   2.9694053782762422D-02  /
   data green(27,18,10) /   2.9448183625409543D-02  /
   data green(27,18,11) /   2.9183433323765627D-02  /
   data green(27,18,12) /   2.8901528288001653D-02  /
   data green(27,18,13) /   2.8604211043639204D-02  /
   data green(27,18,14) /   2.8293217299048434D-02  /
   data green(27,18,15) /   2.7970254669169130D-02  /
   data green(27,18,16) /   2.7636984308103487D-02  /
   data green(27,18,17) /   2.7295005558237159D-02  /
   data green(27,18,18) /   2.6945843598470135D-02  /
   data green(27,19, 0) /   3.0288383488125230D-02  /
   data green(27,19, 1) /   3.0274482164081830D-02  /
   data green(27,19, 2) /   3.0232893226447267D-02  /
   data green(27,19, 3) /   3.0163959127325418D-02  /
   data green(27,19, 4) /   3.0068241902507932D-02  /
   data green(27,19, 5) /   2.9946510506649093D-02  /
   data green(27,19, 6) /   2.9799723874426136D-02  /
   data green(27,19, 7) /   2.9629010429516763D-02  /
   data green(27,19, 8) /   2.9435644878015491D-02  /
   data green(27,19, 9) /   2.9221023182847153D-02  /
   data green(27,19,10) /   2.8986636621297406D-02  /
   data green(27,19,11) /   2.8734045784048479D-02  /
   data green(27,19,12) /   2.8464855289673126D-02  /
   data green(27,19,13) /   2.8180689874169947D-02  /
   data green(27,19,14) /   2.7883172382372080D-02  /
   data green(27,19,15) /   2.7573904048052624D-02  /
   data green(27,19,16) /   2.7254447312021306D-02  /
   data green(27,19,17) /   2.6926311300224364D-02  /
   data green(27,19,18) /   2.6590939972328937D-02  /
   data green(27,19,19) /   2.6249702858793161D-02  /
   data green(27,20, 0) /   2.9760433930073889D-02  /
   data green(27,20, 1) /   2.9747247633602533D-02  /
   data green(27,20, 2) /   2.9707794072611534D-02  /
   data green(27,20, 3) /   2.9642386889429856D-02  /
   data green(27,20, 4) /   2.9551541113981065D-02  /
   data green(27,20, 5) /   2.9435961954181272D-02  /
   data green(27,20, 6) /   2.9296529779167312D-02  /
   data green(27,20, 7) /   2.9134281914216045D-02  /
   data green(27,20, 8) /   2.8950391967355972D-02  /
   data green(27,20, 9) /   2.8746147463021978D-02  /
   data green(27,20,10) /   2.8522926567689941D-02  /
   data green(27,20,11) /   2.8282174660049019D-02  /
   data green(27,20,12) /   2.8025381430692248D-02  /
   data green(27,20,13) /   2.7754059102197734D-02  /
   data green(27,20,14) /   2.7469722249316503D-02  /
   data green(27,20,15) /   2.7173869580057037D-02  /
   data green(27,20,16) /   2.6867967920056576D-02  /
   data green(27,20,17) /   2.6553438531521551D-02  /
   data green(27,20,18) /   2.6231645799181334D-02  /
   data green(27,20,19) /   2.5903888232286718D-02  /
   data green(27,20,20) /   2.5571391665180701D-02  /
   data green(27,21, 0) /   2.9234183531744676D-02  /
   data green(27,21, 1) /   2.9221685014831775D-02  /
   data green(27,21, 2) /   2.9184285784717261D-02  /
   data green(27,21, 3) /   2.9122272736452126D-02  /
   data green(27,21, 4) /   2.9036117226261756D-02  /
   data green(27,21, 5) /   2.8926465167391016D-02  /
   data green(27,21, 6) /   2.8794123741205720D-02  /
   data green(27,21, 7) /   2.8640045252910015D-02  /
   data green(27,21, 8) /   2.8465308750027521D-02  /
   data green(27,21, 9) /   2.8271100072475505D-02  /
   data green(27,21,10) /   2.8058691015348986D-02  /
   data green(27,21,11) /   2.7829418262232144D-02  /
   data green(27,21,12) /   2.7584662693252125D-02  /
   data green(27,21,13) /   2.7325829595128390D-02  /
   data green(27,21,14) /   2.7054330207883698D-02  /
   data green(27,21,15) /   2.6771564942371286D-02  /
   data green(27,21,16) /   2.6478908501350790D-02  /
   data green(27,21,17) /   2.6177697040370183D-02  /
   data green(27,21,18) /   2.5869217417622912D-02  /
   data green(27,21,19) /   2.5554698507213217D-02  /
   data green(27,21,20) /   2.5235304489452291D-02  /
   data green(27,21,21) /   2.4912129985290413D-02  /
   data green(27,22, 0) /   2.8711234209355421D-02  /
   data green(27,22, 1) /   2.8699394946022840D-02  /
   data green(27,22, 2) /   2.8663965149621506D-02  /
   data green(27,22, 3) /   2.8605206980237496D-02  /
   data green(27,22, 4) /   2.8523551376406470D-02  /
   data green(27,22, 5) /   2.8419589316942256D-02  /
   data green(27,22, 6) /   2.8294060078480156D-02  /
   data green(27,22, 7) /   2.8147836940678007D-02  /
   data green(27,22, 8) /   2.7981910868712104D-02  /
   data green(27,22, 9) /   2.7797372748768843D-02  /
   data green(27,22,10) /   2.7595394766175854D-02  /
   data green(27,22,11) /   2.7377211499679522D-02  /
   data green(27,22,12) /   2.7144101263258082D-02  /
   data green(27,22,13) /   2.6897368164292494D-02  /
   data green(27,22,14) /   2.6638325270178675D-02  /
   data green(27,22,15) /   2.6368279190909037D-02  /
   data green(27,22,16) /   2.6088516298678522D-02  /
   data green(27,22,17) /   2.5800290722209317D-02  /
   data green(27,22,18) /   2.5504814177163829D-02  /
   data green(27,22,19) /   2.5203247627465859D-02  /
   data green(27,22,20) /   2.4896694717161400D-02  /
   data green(27,22,21) /   2.4586196869194641D-02  /
   data green(27,22,22) /   2.4272729915897073D-02  /
   data green(27,23, 0) /   2.8192974987878067D-02  /
   data green(27,23, 1) /   2.8181765700681247D-02  /
   data green(27,23, 2) /   2.8148218161073959D-02  /
   data green(27,23, 3) /   2.8092571732870848D-02  /
   data green(27,23, 4) /   2.8015220078133653D-02  /
   data green(27,23, 5) /   2.7916703456511000D-02  /
   data green(27,23, 6) /   2.7797698361565978D-02  /
   data green(27,23, 7) /   2.7659004879349282D-02  /
   data green(27,23, 8) /   2.7501532222273473D-02  /
   data green(27,23, 9) /   2.7326282932938654D-02  /
   data green(27,23,10) /   2.7134336267348892D-02  /
   data green(27,23,11) /   2.6926831256393273D-02  /
   data green(27,23,12) /   2.6704949911706126D-02  /
   data green(27,23,13) /   2.6469900991452668D-02  /
   data green(27,23,14) /   2.6222904678278477D-02  /
   data green(27,23,15) /   2.5965178450835315D-02  /
   data green(27,23,16) /   2.5697924356886566D-02  /
   data green(27,23,17) /   2.5422317824290559D-02  /
   data green(27,23,18) /   2.5139498079601634D-02  /
   data green(27,23,19) /   2.4850560185104009D-02  /
   data green(27,23,20) /   2.4556548655353175D-02  /
   data green(27,23,21) /   2.4258452574441688D-02  /
   data green(27,23,22) /   2.3957202105211361D-02  /
   data green(27,23,23) /   2.3653666260937324D-02  /
   data green(27,24, 0) /   2.7680596824373974D-02  /
   data green(27,24, 1) /   2.7669987929523597D-02  /
   data green(27,24, 2) /   2.7638234519796866D-02  /
   data green(27,24, 3) /   2.7585555011023305D-02  /
   data green(27,24, 4) /   2.7512308787332348D-02  /
   data green(27,24, 5) /   2.7418989420726080D-02  /
   data green(27,24, 6) /   2.7306215532688829D-02  /
   data green(27,24, 7) /   2.7174719625855197D-02  /
   data green(27,24, 8) /   2.7025335272777083D-02  /
   data green(27,24, 9) /   2.6858983086179393D-02  /
   data green(27,24,10) /   2.6676655910117263D-02  /
   data green(27,24,11) /   2.6479403665149757D-02  /
   data green(27,24,12) /   2.6268318255451541D-02  /
   data green(27,24,13) /   2.6044518905147838D-02  /
   data green(27,24,14) /   2.5809138239161290D-02  /
   data green(27,24,15) /   2.5563309364771834D-02  /
   data green(27,24,16) /   2.5308154148017747D-02  /
   data green(27,24,17) /   2.5044772817645056D-02  /
   data green(27,24,18) /   2.4774234971534612D-02  /
   data green(27,24,19) /   2.4497572008637163D-02  /
   data green(27,24,20) /   2.4215770964898457D-02  /
   data green(27,24,21) /   2.3929769695228937D-02  /
   data green(27,24,22) /   2.3640453315414797D-02  /
   data green(27,24,23) /   2.3348651797655166D-02  /
   data green(27,24,24) /   2.3055138600453138D-02  /
   data green(27,25, 0) /   2.7175108413893354D-02  /
   data green(27,25, 1) /   2.7165070390824404D-02  /
   data green(27,25, 2) /   2.7135023139064681D-02  /
   data green(27,25, 3) /   2.7085165873149176D-02  /
   data green(27,25, 4) /   2.7015826536344739D-02  /
   data green(27,25, 5) /   2.6927455834216556D-02  /
   data green(27,25, 6) /   2.6820619182077007D-02  /
   data green(27,25, 7) /   2.6695986845360824D-02  /
   data green(27,25, 8) /   2.6554322603240955D-02  /
   data green(27,25, 9) /   2.6396471299174204D-02  /
   data green(27,25,10) /   2.6223345656872902D-02  /
   data green(27,25,11) /   2.6035912737117972D-02  /
   data green(27,25,12) /   2.5835180391697913D-02  /
   data green(27,25,13) /   2.5622184038297337D-02  /
   data green(27,25,14) /   2.5397974037630069D-02  /
   data green(27,25,15) /   2.5163603905002170D-02  /
   data green(27,25,16) /   2.4920119536193146D-02  /
   data green(27,25,17) /   2.4668549575134962D-02  /
   data green(27,25,18) /   2.4409897000921891D-02  /
   data green(27,25,19) /   2.4145131966193053D-02  /
   data green(27,25,20) /   2.3875185879261542D-02  /
   data green(27,25,21) /   2.3600946689316112D-02  /
   data green(27,25,22) /   2.3323255307857206D-02  /
   data green(27,25,23) /   2.3042903080107282D-02  /
   data green(27,25,24) /   2.2760630206979864D-02  /
   data green(27,25,25) /   2.2477125010618534D-02  /
   data green(27,26, 0) /   2.6677352366456082D-02  /
   data green(27,26, 1) /   2.6667856058800444D-02  /
   data green(27,26, 2) /   2.6639428048850285D-02  /
   data green(27,26, 3) /   2.6592249988350684D-02  /
   data green(27,26, 4) /   2.6526621043997091D-02  /
   data green(27,26, 5) /   2.6442952649381124D-02  /
   data green(27,26, 6) /   2.6341761412450433D-02  /
   data green(27,26, 7) /   2.6223660415720405D-02  /
   data green(27,26, 8) /   2.6089349190937965D-02  /
   data green(27,26, 9) /   2.5939602679595768D-02  /
   data green(27,26,10) /   2.5775259504971153D-02  /
   data green(27,26,11) /   2.5597209880655825D-02  /
   data green(27,26,12) /   2.5406383466236821D-02  /
   data green(27,26,13) /   2.5203737455012625D-02  /
   data green(27,26,14) /   2.4990245143986862D-02  /
   data green(27,26,15) /   2.4766885195706093D-02  /
   data green(27,26,16) /   2.4534631757596334D-02  /
   data green(27,26,17) /   2.4294445559889653D-02  /
   data green(27,26,18) /   2.4047266070223242D-02  /
   data green(27,26,19) /   2.3794004743291814D-02  /
   data green(27,26,20) /   2.3535539368799964D-02  /
   data green(27,26,21) /   2.3272709491174236D-02  /
   data green(27,26,22) /   2.3006312850403287D-02  /
   data green(27,26,23) /   2.2737102774961399D-02  /
   data green(27,26,24) /   2.2465786444735177D-02  /
   data green(27,26,25) /   2.2193023933709505D-02  /
   data green(27,26,26) /   2.1919427938245462D-02  /
   data green(27,27, 0) /   2.6188021277610565D-02  /
   data green(27,27, 1) /   2.6179038133066251D-02  /
   data green(27,27, 2) /   2.6152144223210386D-02  /
   data green(27,27, 3) /   2.6107505164414305D-02  /
   data green(27,27, 4) /   2.6045393833035592D-02  /
   data green(27,27, 5) /   2.5966185750166924D-02  /
   data green(27,27, 6) /   2.5870352836206902D-02  /
   data green(27,27, 7) /   2.5758455736879040D-02  /
   data green(27,27, 8) /   2.5631134960829201D-02  /
   data green(27,27, 9) /   2.5489101095263599D-02  /
   data green(27,27,10) /   2.5333124379622741D-02  /
   data green(27,27,11) /   2.5164023918285282D-02  /
   data green(27,27,12) /   2.4982656802798476D-02  /
   data green(27,27,13) /   2.4789907393791881D-02  /
   data green(27,27,14) /   2.4586676984627038D-02  /
   data green(27,27,15) /   2.4373874035253699D-02  /
   data green(27,27,16) /   2.4152405127986824D-02  /
   data green(27,27,17) /   2.3923166759148327D-02  /
   data green(27,27,18) /   2.3687038043610050D-02  /
   data green(27,27,19) /   2.3444874374768276D-02  /
   data green(27,27,20) /   2.3197502051505673D-02  /
   data green(27,27,21) /   2.2945713857012898D-02  /
   data green(27,27,22) /   2.2690265552330866D-02  /
   data green(27,27,23) /   2.2431873230220004D-02  /
   data green(27,27,24) /   2.2171211462289252D-02  /
   data green(27,27,25) /   2.1908912163868689D-02  /
   data green(27,27,26) /   2.1645564096402503D-02  /
   data green(27,27,27) /   2.1381712925620213D-02  /
   data green(28, 0, 0) /   3.5725715913662262D-02  /
   data green(28, 1, 0) /   3.5702865729838593D-02  /
   data green(28, 1, 1) /   3.5680059618249202D-02  /
   data green(28, 2, 0) /   3.5634579613798165D-02  /
   data green(28, 2, 1) /   3.5611904862181319D-02  /
   data green(28, 2, 2) /   3.5544141639403608D-02  /
   data green(28, 3, 0) /   3.5521642304283604D-02  /
   data green(28, 3, 1) /   3.5499183669460116D-02  /
   data green(28, 3, 2) /   3.5432064629821153D-02  /
   data green(28, 3, 3) /   3.5321047557433559D-02  /
   data green(28, 4, 0) /   3.5365333690571585D-02  /
   data green(28, 4, 1) /   3.5343171838793937D-02  /
   data green(28, 4, 2) /   3.5276937471363762D-02  /
   data green(28, 4, 3) /   3.5167376185488985D-02  /
   data green(28, 4, 4) /   3.5015704588320269D-02  /
   data green(28, 5, 0) /   3.5167388664943310D-02  /
   data green(28, 5, 1) /   3.5145598792628810D-02  /
   data green(28, 5, 2) /   3.5080473323211243D-02  /
   data green(28, 5, 3) /   3.4972737046708890D-02  /
   data green(28, 5, 4) /   3.4823572912581535D-02  /
   data green(28, 5, 5) /   3.4634586090822386D-02  /
   data green(28, 6, 0) /   3.4929944443753115D-02  /
   data green(28, 6, 1) /   3.4908595130133557D-02  /
   data green(28, 6, 2) /   3.4844783106563935D-02  /
   data green(28, 6, 3) /   3.4739208834439414D-02  /
   data green(28, 6, 4) /   3.4593015893462128D-02  /
   data green(28, 6, 5) /   3.4407756711814577D-02  /
   data green(28, 6, 6) /   3.4185347459698825D-02  /
   data green(28, 7, 0) /   3.4655478443701690D-02  /
   data green(28, 7, 1) /   3.4634630769493684D-02  /
   data green(28, 7, 2) /   3.4572314436128793D-02  /
   data green(28, 7, 3) /   3.4469202620505207D-02  /
   data green(28, 7, 4) /   3.4326394728644176D-02  /
   data green(28, 7, 5) /   3.4145383970306378D-02  /
   data green(28, 7, 6) /   3.3928014635700110D-02  /
   data green(28, 7, 7) /   3.3676431433720962D-02  /
   data green(28, 8, 0) /   3.4346740129447197D-02  /
   data green(28, 8, 1) /   3.4326447077123165D-02  /
   data green(28, 8, 2) /   3.4265784585104070D-02  /
   data green(28, 8, 3) /   3.4165396183641922D-02  /
   data green(28, 8, 4) /   3.4026333255614112D-02  /
   data green(28, 8, 5) /   3.3850024581692756D-02  /
   data green(28, 8, 6) /   3.3638236167015785D-02  /
   data green(28, 8, 7) /   3.3393023529875505D-02  /
   data green(28, 8, 8) /   3.3116678897411093D-02  /
   data green(28, 9, 0) /   3.4006680281200664D-02  /
   data green(28, 9, 1) /   3.3986986414748566D-02  /
   data green(28, 9, 2) /   3.3928110862203639D-02  /
   data green(28, 9, 3) /   3.3830665741865135D-02  /
   data green(28, 9, 4) /   3.3695651524745733D-02  /
   data green(28, 9, 5) /   3.3524428628548902D-02  /
   data green(28, 9, 6) /   3.3318679899943603D-02  /
   data green(28, 9, 7) /   3.3080365982744116D-02  /
   data green(28, 9, 8) /   3.2811675818730270D-02  /
   data green(28, 9, 9) /   3.2514974597123741D-02  /
   data green(28,10, 0) /   3.3638380899390953D-02  /
   data green(28,10, 1) /   3.3619322306149980D-02  /
   data green(28,10, 2) /   3.3562341557128381D-02  /
   data green(28,10, 3) /   3.3468018171267953D-02  /
   data green(28,10, 4) /   3.3337299750652347D-02  /
   data green(28,10, 5) /   3.3171475657048741D-02  /
   data green(28,10, 6) /   3.2972142198308393D-02  /
   data green(28,10, 7) /   3.2741161139743895D-02  /
   data green(28,10, 8) /   3.2480613588588114D-02  /
   data green(28,10, 9) /   3.2192751371670020D-02  /
   data green(28,10,10) /   3.1879947948450685D-02  /
   data green(28,11, 0) /   3.3244988526832324D-02  /
   data green(28,11, 1) /   3.3226592992122723D-02  /
   data green(28,11, 2) /   3.3171590187793122D-02  /
   data green(28,11, 3) /   3.3080526386944656D-02  /
   data green(28,11, 4) /   3.2954295239197406D-02  /
   data green(28,11, 5) /   3.2794113523860732D-02  /
   data green(28,11, 6) /   3.2601489037768668D-02  /
   data green(28,11, 7) /   3.2378182254572035D-02  /
   data green(28,11, 8) /   3.2126163608637751D-02  /
   data green(28,11, 9) /   3.1847568330304304D-02  /
   data green(28,11,10) /   3.1544650698716474D-02  /
   data green(28,11,11) /   3.1219739405031127D-02  /
   data green(28,12, 0) /   3.2829653199415415D-02  /
   data green(28,12, 1) /   3.2811940576718109D-02  /
   data green(28,12, 2) /   3.2758975228648679D-02  /
   data green(28,12, 3) /   3.2671270026416517D-02  /
   data green(28,12, 4) /   3.2549664377171002D-02  /
   data green(28,12, 5) /   3.2395302013905045D-02  /
   data green(28,12, 6) /   3.2209601537945450D-02  /
   data green(28,12, 7) /   3.1994221179680767D-02  /
   data green(28,12, 8) /   3.1751019442719700D-02  /
   data green(28,12, 9) /   3.1482013370654933D-02  /
   data green(28,12,10) /   3.1189336130650096D-02  /
   data green(28,12,11) /   3.0875195461827036D-02  /
   data green(28,12,12) /   3.0541834313908386D-02  /
   data green(28,13, 0) /   3.2395474607574706D-02  /
   data green(28,13, 1) /   3.2378457342921232D-02  /
   data green(28,13, 2) /   3.2327566885982935D-02  /
   data green(28,13, 3) /   3.2243282981033329D-02  /
   data green(28,13, 4) /   3.2126391201662371D-02  /
   data green(28,13, 5) /   3.1977962709718775D-02  /
   data green(28,13, 6) /   3.1799327369319570D-02  /
   data green(28,13, 7) /   3.1592041519898889D-02  /
   data green(28,13, 8) /   3.1357851895098592D-02  /
   data green(28,13, 9) /   3.1098657247544480D-02  /
   data green(28,13,10) /   3.0816469208122768D-02  /
   data green(28,13,11) /   3.0513373786767311D-02  /
   data green(28,13,12) /   3.0191494731039405D-02  /
   data green(28,13,13) /   2.9852959723128531D-02  /
   data green(28,14, 0) /   3.1945456427447859D-02  /
   data green(28,14, 1) /   3.1929140196024192D-02  /
   data green(28,14, 2) /   3.1880341876376485D-02  /
   data green(28,14, 3) /   3.1799508727693208D-02  /
   data green(28,14, 4) /   3.1687373495397873D-02  /
   data green(28,14, 5) /   3.1544936049812941D-02  /
   data green(28,14, 6) /   3.1373438960001596D-02  /
   data green(28,14, 7) /   3.1174338156126944D-02  /
   data green(28,14, 8) /   3.0949270000143286D-02  /
   data green(28,14, 9) /   3.0700016156133811D-02  /
   data green(28,14,10) /   3.0428467631639394D-02  /
   data green(28,14,11) /   3.0136589261637328D-02  /
   data green(28,14,12) /   2.9826385744966225D-02  /
   data green(28,14,13) /   2.9499870139477427D-02  /
   data green(28,14,14) /   2.9159035497823661D-02  /
   data green(28,15, 0) /   3.1482469212844130D-02  /
   data green(28,15, 1) /   3.1466853627920975D-02  /
   data green(28,15, 2) /   3.1420146607664935D-02  /
   data green(28,15, 3) /   3.1342763869798604D-02  /
   data green(28,15, 4) /   3.1235386829080291D-02  /
   data green(28,15, 5) /   3.1098946010608242D-02  /
   data green(28,15, 6) /   3.0934598950717918D-02  /
   data green(28,15, 7) /   3.0743703600181384D-02  /
   data green(28,15, 8) /   3.0527788395109889D-02  /
   data green(28,15, 9) /   3.0288520229871083D-02  /
   data green(28,15,10) /   3.0027671555827147D-02  /
   data green(28,15,11) /   2.9747087749182825D-02  /
   data green(28,15,12) /   2.9448655755185855D-02  /
   data green(28,15,13) /   2.9134274841590216D-02  /
   data green(28,15,14) /   2.8805830099443562D-02  /
   data green(28,15,15) /   2.8465169130273506D-02  /
   data green(28,16, 0) /   3.1009221759938716D-02  /
   data green(28,16, 1) /   3.0994301118335021D-02  /
   data green(28,16, 2) /   3.0949668690457347D-02  /
   data green(28,16, 3) /   3.0875709834768714D-02  /
   data green(28,16, 4) /   3.0773056523432655D-02  /
   data green(28,16, 5) /   3.0642572414724491D-02  /
   data green(28,16, 6) /   3.0485332933933195D-02  /
   data green(28,16, 7) /   3.0302601251251463D-02  /
   data green(28,16, 8) /   3.0095801180810493D-02  /
   data green(28,16, 9) /   2.9866488090661645D-02  /
   data green(28,16,10) /   2.9616318910606211D-02  /
   data green(28,16,11) /   2.9347022260762211D-02  /
   data green(28,16,12) /   2.9060369610474937D-02  /
   data green(28,16,13) /   2.8758148229005177D-02  /
   data green(28,16,14) /   2.8442136521412708D-02  /
   data green(28,16,15) /   2.8114082169308378D-02  /
   data green(28,16,16) /   2.7775683328690803D-02  /
   data green(28,17, 0) /   3.0528240482617118D-02  /
   data green(28,17, 1) /   3.0514004516217445D-02  /
   data green(28,17, 2) /   3.0471416339132928D-02  /
   data green(28,17, 3) /   3.0400832312331923D-02  /
   data green(28,17, 4) /   3.0302837148901755D-02  /
   data green(28,17, 5) /   3.0178230524208486D-02  /
   data green(28,17, 6) /   3.0028009181455002D-02  /
   data green(28,17, 7) /   2.9853345307328066D-02  /
   data green(28,17, 8) /   2.9655562073928995D-02  /
   data green(28,17, 9) /   2.9436107305062340D-02  /
   data green(28,17,10) /   2.9196526227986991D-02  /
   data green(28,17,11) /   2.8938434221718247D-02  /
   data green(28,17,12) /   2.8663490379509077D-02  /
   data green(28,17,13) /   2.8373372578128414D-02  /
   data green(28,17,14) /   2.8069754602646792D-02  /
   data green(28,17,15) /   2.7754285724703277D-02  /
   data green(28,17,16) /   2.7428572985132207D-02  /
   data green(28,17,17) /   2.7094166296603855D-02  /
   data green(28,18, 0) /   3.0041856070398259D-02  /
   data green(28,18, 1) /   3.0028290679101986D-02  /
   data green(28,18, 2) /   2.9987704957428481D-02  /
   data green(28,18, 3) /   2.9920427756789984D-02  /
   data green(28,18, 4) /   2.9826998923906042D-02  /
   data green(28,18, 5) /   2.9708157325403693D-02  /
   data green(28,18, 6) /   2.9564824819436893D-02  /
   data green(28,18, 7) /   2.9398086847257962D-02  /
   data green(28,18, 8) /   2.9209170425994167D-02  /
   data green(28,18, 9) /   2.8999420381658619D-02  /
   data green(28,18,10) /   2.8770274668919911D-02  /
   data green(28,18,11) /   2.8523239585840597D-02  /
   data green(28,18,12) /   2.8259865615388315D-02  /
   data green(28,18,13) /   2.7981724520799812D-02  /
   data green(28,18,14) /   2.7690388199395764D-02  /
   data green(28,18,15) /   2.7387409669458043D-02  /
   data green(28,18,16) /   2.7074306436334929D-02  /
   data green(28,18,17) /   2.6752546364329992D-02  /
   data green(28,18,18) /   2.6423536075410807D-02  /
   data green(28,19, 0) /   2.9552196535419458D-02  /
   data green(28,19, 1) /   2.9539284482759061D-02  /
   data green(28,19, 2) /   2.9500650038322765D-02  /
   data green(28,19, 3) /   2.9436596111194938D-02  /
   data green(28,19, 4) /   2.9347620207661063D-02  /
   data green(28,19, 5) /   2.9234403748117324D-02  /
   data green(28,19, 6) /   2.9097797747015673D-02  /
   data green(28,19, 7) /   2.8938805438473308D-02  /
   data green(28,19, 8) /   2.8758562526328475D-02  /
   data green(28,19, 9) /   2.8558315790911958D-02  /
   data green(28,19,10) /   2.8339400795512840D-02  /
   data green(28,19,11) /   2.8103219406846791D-02  /
   data green(28,19,12) /   2.7851217781946516D-02  /
   data green(28,19,13) /   2.7584865386752963D-02  /
   data green(28,19,14) /   2.7305635508063534D-02  /
   data green(28,19,15) /   2.7014987609029515D-02  /
   data green(28,19,16) /   2.6714351766862647D-02  /
   data green(28,19,17) /   2.6405115326213365D-02  /
   data green(28,19,18) /   2.6088611807576580D-02  /
   data green(28,19,19) /   2.5766112030152739D-02  /
   data green(28,20, 0) /   2.9061185676125409D-02  /
   data green(28,20, 1) /   2.9048907234262407D-02  /
   data green(28,20, 2) /   2.9012165427566457D-02  /
   data green(28,20, 3) /   2.8951238829263124D-02  /
   data green(28,20, 4) /   2.8866585199284648D-02  /
   data green(28,20, 5) /   2.8758831974948478D-02  /
   data green(28,20, 6) /   2.8628763505561660D-02  /
   data green(28,20, 7) /   2.8477305533961682D-02  /
   data green(28,20, 8) /   2.8305507513033606D-02  /
   data green(28,20, 9) /   2.8114523394378928D-02  /
   data green(28,20,10) /   2.7905591539194041D-02  /
   data green(28,20,11) /   2.7680014380587769D-02  /
   data green(28,20,12) /   2.7439138416904676D-02  /
   data green(28,20,13) /   2.7184335043584699D-02  /
   data green(28,20,14) /   2.6916982643912631D-02  /
   data green(28,20,15) /   2.6638450263935801D-02  /
   data green(28,20,16) /   2.6350083100489908D-02  /
   data green(28,20,17) /   2.6053189939244981D-02  /
   data green(28,20,18) /   2.5749032596235427D-02  /
   data green(28,20,19) /   2.5438817344357607D-02  /
   data green(28,20,20) /   2.5123688247387980D-02  /
   data green(28,21, 0) /   2.8570545975938378D-02  /
   data green(28,21, 1) /   2.8558879511577345D-02  /
   data green(28,21, 2) /   2.8523965988326620D-02  /
   data green(28,21, 3) /   2.8466061255979325D-02  /
   data green(28,21, 4) /   2.8385585936084607D-02  /
   data green(28,21, 5) /   2.8283116973222288D-02  /
   data green(28,21, 6) /   2.8159376277056662D-02  /
   data green(28,21, 7) /   2.8015216888390242D-02  /
   data green(28,21, 8) /   2.7851607177368762D-02  /
   data green(28,21, 9) /   2.7669613626821327D-02  /
   data green(28,21,10) /   2.7470382767928014D-02  /
   data green(28,21,11) /   2.7255122820863797D-02  /
   data green(28,21,12) /   2.7025085553605888D-02  /
   data green(28,21,13) /   2.6781548812917100D-02  /
   data green(28,21,14) /   2.6525800108563438D-02  /
   data green(28,21,15) /   2.6259121551131064D-02  /
   data green(28,21,16) /   2.5982776360997922D-02  /
   data green(28,21,17) /   2.5697997085913051D-02  /
   data green(28,21,18) /   2.5405975591034389D-02  /
   data green(28,21,19) /   2.5107854820839305D-02  /
   data green(28,21,20) /   2.4804722278634238D-02  /
   data green(28,21,21) /   2.4497605127057013D-02  /
   data green(28,22, 0) /   2.8081804997585186D-02  /
   data green(28,22, 1) /   2.8070727494411588D-02  /
   data green(28,22, 2) /   2.8037573743783451D-02  /
   data green(28,22, 3) /   2.7982578463565229D-02  /
   data green(28,22, 4) /   2.7906127714806328D-02  /
   data green(28,22, 5) /   2.7808751407027342D-02  /
   data green(28,22, 6) /   2.7691113208645944D-02  /
   data green(28,22, 7) /   2.7553998234894527D-02  /
   data green(28,22, 8) /   2.7398298951411909D-02  /
   data green(28,22, 9) /   2.7224999772347366D-02  /
   data green(28,22,10) /   2.7035160846680208D-02  /
   data green(28,22,11) /   2.6829901516879102D-02  /
   data green(28,22,12) /   2.6610383902995360D-02  /
   data green(28,22,13) /   2.6377797016975221D-02  /
   data green(28,22,14) /   2.6133341751231411D-02  /
   data green(28,22,15) /   2.5878217017342178D-02  /
   data green(28,22,16) /   2.5613607239889653D-02  /
   data green(28,22,17) /   2.5340671341053427D-02  /
   data green(28,22,18) /   2.5060533286981716D-02  /
   data green(28,22,19) /   2.4774274209598851D-02  /
   data green(28,22,20) /   2.4482926068916232D-02  /
   data green(28,22,21) /   2.4187466781802756D-02  /
   data green(28,22,22) /   2.3888816713566772D-02  /
   data green(28,23, 0) /   2.7596304411424077D-02  /
   data green(28,23, 1) /   2.7585791927935113D-02  /
   data green(28,23, 2) /   2.7554326649012379D-02  /
   data green(28,23, 3) /   2.7502123710142785D-02  /
   data green(28,23, 4) /   2.7429537125122038D-02  /
   data green(28,23, 5) /   2.7337053146112268D-02  /
   data green(28,23, 6) /   2.7225281312565488D-02  /
   data green(28,23, 7) /   2.7094943508507864D-02  /
   data green(28,23, 8) /   2.6946861405360185D-02  /
   data green(28,23, 9) /   2.6781942704140055D-02  /
   data green(28,23,10) /   2.6601166605894254D-02  /
   data green(28,23,11) /   2.6405568933490785D-02  /
   data green(28,23,12) /   2.6196227303779015D-02  /
   data green(28,23,13) /   2.5974246709929021D-02  /
   data green(28,23,14) /   2.5740745823423141D-02  /
   data green(28,23,15) /   2.5496844267827414D-02  /
   data green(28,23,16) /   2.5243651056098763D-02  /
   data green(28,23,17) /   2.4982254323324536D-02  /
   data green(28,23,18) /   2.4713712430366224D-02  /
   data green(28,23,19) /   2.4439046463078144D-02  /
   data green(28,23,20) /   2.4159234108057773D-02  /
   data green(28,23,21) /   2.3875204850023875D-02  /
   data green(28,23,22) /   2.3587836408091346D-02  /
   data green(28,23,23) /   2.3297952308117806D-02  /
   data green(28,24, 0) /   2.7115210894845327D-02  /
   data green(28,24, 1) /   2.7105238958937149D-02  /
   data green(28,24, 2) /   2.7075389239166654D-02  /
   data green(28,24, 3) /   2.7025858780400952D-02  /
   data green(28,24, 4) /   2.6956971971626005D-02  /
   data green(28,24, 5) /   2.6869174669262381D-02  /
   data green(28,24, 6) /   2.6763026264490931D-02  /
   data green(28,24, 7) /   2.6639189968358200D-02  /
   data green(28,24, 8) /   2.6498421638850005D-02  /
   data green(28,24, 9) /   2.6341557507041223D-02  /
   data green(28,24,10) /   2.6169501174182722D-02  /
   data green(28,24,11) /   2.5983210248816006D-02  /
   data green(28,24,12) /   2.5783682974498293D-02  /
   data green(28,24,13) /   2.5571945167111351D-02  /
   data green(28,24,14) /   2.5349037739204109D-02  /
   data green(28,24,15) /   2.5116005040775199D-02  /
   data green(28,24,16) /   2.4873884194664068D-02  /
   data green(28,24,17) /   2.4623695553295303D-02  /
   data green(28,24,18) /   2.4366434354442686D-02  /
   data green(28,24,19) /   2.4103063608908067D-02  /
   data green(28,24,20) /   2.3834508213902553D-02  /
   data green(28,24,21) /   2.3561650253267251D-02  /
   data green(28,24,22) /   2.3285325419755466D-02  /
   data green(28,24,23) /   2.3006320475291418D-02  /
   data green(28,24,24) /   2.2725371651968633D-02  /
   data green(28,25, 0) /   2.6639528248164370D-02  /
   data green(28,25, 1) /   2.6630072191640996D-02  /
   data green(28,25, 2) /   2.6601764506571954D-02  /
   data green(28,25, 3) /   2.6554785569812011D-02  /
   data green(28,25, 4) /   2.6489432458285000D-02  /
   data green(28,25, 5) /   2.6406113751858300D-02  /
   data green(28,25, 6) /   2.6305342508943369D-02  /
   data green(28,25, 7) /   2.6187727649162312D-02  /
   data green(28,25, 8) /   2.6053964021404120D-02  /
   data green(28,25, 9) /   2.5904821465022293D-02  /
   data green(28,25,10) /   2.5741133186148919D-02  /
   data green(28,25,11) /   2.5563783770538355D-02  /
   data green(28,25,12) /   2.5373697140365345D-02  /
   data green(28,25,13) /   2.5171824737074130D-02  /
   data green(28,25,14) /   2.4959134178271603D-02  /
   data green(28,25,15) /   2.4736598596559849D-02  /
   data green(28,25,16) /   2.4505186824875690D-02  /
   data green(28,25,17) /   2.4265854548891443D-02  /
   data green(28,25,18) /   2.4019536504514022D-02  /
   data green(28,25,19) /   2.3767139759235823D-02  /
   data green(28,25,20) /   2.3509538081287022D-02  /
   data green(28,25,21) /   2.3247567370997233D-02  /
   data green(28,25,22) /   2.2982022104844069D-02  /
   data green(28,25,23) /   2.2713652724339777D-02  /
   data green(28,25,24) /   2.2443163888894252D-02  /
   data green(28,25,25) /   2.2171213503599612D-02  /
   data green(28,26, 0) /   2.6170110181745643D-02  /
   data green(28,26, 1) /   2.6161145419120381D-02  /
   data green(28,26, 2) /   2.6134306466297712D-02  /
   data green(28,26, 3) /   2.6089758378003154D-02  /
   data green(28,26, 4) /   2.6027773109409003D-02  /
   data green(28,26, 5) /   2.5948724922547968D-02  /
   data green(28,26, 6) /   2.5853084169918054D-02  /
   data green(28,26, 7) /   2.5741409655712023D-02  /
   data green(28,26, 8) /   2.5614339813406630D-02  /
   data green(28,26, 9) /   2.5472582964659967D-02  /
   data green(28,26,10) /   2.5316906937967364D-02  /
   data green(28,26,11) /   2.5148128326575335D-02  /
   data green(28,26,12) /   2.4967101654776138D-02  /
   data green(28,26,13) /   2.4774708701541710D-02  /
   data green(28,26,14) /   2.4571848202566932D-02  /
   data green(28,26,15) /   2.4359426118444466D-02  /
   data green(28,26,16) /   2.4138346620178316D-02  /
   data green(28,26,17) /   2.3909503905701096D-02  /
   data green(28,26,18) /   2.3673774924362300D-02  /
   data green(28,26,19) /   2.3432013052027672D-02  /
   data green(28,26,20) /   2.3185042728601994D-02  /
   data green(28,26,21) /   2.2933655043214186D-02  /
   data green(28,26,22) /   2.2678604230370512D-02  /
   data green(28,26,23) /   2.2420605023171062D-02  /
   data green(28,26,24) /   2.2160330797030434D-02  /
   data green(28,26,25) /   2.1898412428889572D-02  /
   data green(28,26,26) /   2.1635437792176221D-02  /
   data green(28,27, 0) /   2.5707673333309321D-02  /
   data green(28,27, 1) /   2.5699175589987855D-02  /
   data green(28,27, 2) /   2.5673732972040626D-02  /
   data green(28,27, 3) /   2.5631496476769249D-02  /
   data green(28,27, 4) /   2.5572714997819053D-02  /
   data green(28,27, 5) /   2.5497731266498309D-02  /
   data green(28,27, 6) /   2.5406976352647298D-02  /
   data green(28,27, 7) /   2.5300962896398670D-02  /
   data green(28,27, 8) /   2.5180277275478266D-02  /
   data green(28,27, 9) /   2.5045570935969024D-02  /
   data green(28,27,10) /   2.4897551127122573D-02  /
   data green(28,27,11) /   2.4736971282989410D-02  /
   data green(28,27,12) /   2.4564621286120022D-02  /
   data green(28,27,13) /   2.4381317832650219D-02  /
   data green(28,27,14) /   2.4187895095371486D-02  /
   data green(28,27,15) /   2.3985195853748031D-02  /
   data green(28,27,16) /   2.3774063229163145D-02  /
   data green(28,27,17) /   2.3555333131756812D-02  /
   data green(28,27,18) /   2.3329827493648581D-02  /
   data green(28,27,19) /   2.3098348333455129D-02  /
   data green(28,27,20) /   2.2861672669826479D-02  /
   data green(28,27,21) /   2.2620548277946202D-02  /
   data green(28,27,22) /   2.2375690262972849D-02  /
   data green(28,27,23) /   2.2127778408388546D-02  /
   data green(28,27,24) /   2.1877455245088218D-02  /
   data green(28,27,25) /   2.1625324778550663D-02  /
   data green(28,27,26) /   2.1371951806217600D-02  /
   data green(28,27,27) /   2.1117861754837895D-02  /
   data green(28,28, 0) /   2.5252810169591959D-02  /
   data green(28,28, 1) /   2.5244755664857170D-02  /
   data green(28,28, 2) /   2.5220638437838030D-02  /
   data green(28,28, 3) /   2.5180596610000942D-02  /
   data green(28,28, 4) /   2.5124857940072128D-02  /
   data green(28,28, 5) /   2.5053736238612463D-02  /
   data green(28,28, 6) /   2.4967626504461882D-02  /
   data green(28,28, 7) /   2.4866998929481000D-02  /
   data green(28,28, 8) /   2.4752391946949094D-02  /
   data green(28,28, 9) /   2.4624404519585516D-02  /
   data green(28,28,10) /   2.4483687874915200D-02  /
   data green(28,28,11) /   2.4330936898653094D-02  /
   data green(28,28,12) /   2.4166881391506342D-02  /
   data green(28,28,13) /   2.3992277382288388D-02  /
   data green(28,28,14) /   2.3807898671820944D-02  /
   data green(28,28,15) /   2.3614528759263741D-02  /
   data green(28,28,16) /   2.3412953276812823D-02  /
   data green(28,28,17) /   2.3203953031643573D-02  /
   data green(28,28,18) /   2.2988297726905183D-02  /
   data green(28,28,19) /   2.2766740407644731D-02  /
   data green(28,28,20) /   2.2540012653661217D-02  /
   data green(28,28,21) /   2.2308820520113393D-02  /
   data green(28,28,22) /   2.2073841208637899D-02  /
   data green(28,28,23) /   2.1835720436956422D-02  /
   data green(28,28,24) /   2.1595070463461104D-02  /
   data green(28,28,25) /   2.1352468714908139D-02  /
   data green(28,28,26) /   2.1108456959858527D-02  /
   data green(28,28,27) /   2.0863540967543993D-02  /
   data green(28,28,28) /   2.0618190591030910D-02  /
   data green(29, 0, 0) /   3.4493044054713561D-02  /
   data green(29, 1, 0) /   3.4472481609434176D-02  /
   data green(29, 1, 1) /   3.4451956115521244D-02  /
   data green(29, 2, 0) /   3.4411015983262509D-02  /
   data green(29, 2, 1) /   3.4390600674560061D-02  /
   data green(29, 2, 2) /   3.4329573799682531D-02  /
   data green(29, 3, 0) /   3.4309305616143898D-02  /
   data green(29, 3, 1) /   3.4289071748361280D-02  /
   data green(29, 3, 2) /   3.4228585936442374D-02  /
   data green(29, 3, 3) /   3.4128489116769176D-02  /
   data green(29, 4, 0) /   3.4168426010611756D-02  /
   data green(29, 4, 1) /   3.4148441635455495D-02  /
   data green(29, 4, 2) /   3.4088699852294409D-02  /
   data green(29, 4, 3) /   3.3989828437012842D-02  /
   data green(29, 4, 4) /   3.3852853211709555D-02  /
   data green(29, 5, 0) /   3.3989838243942527D-02  /
   data green(29, 5, 1) /   3.3970167115091744D-02  /
   data green(29, 5, 2) /   3.3911359539557026D-02  /
   data green(29, 5, 3) /   3.3814026929322720D-02  /
   data green(29, 5, 4) /   3.3679168582480462D-02  /
   data green(29, 5, 5) /   3.3508143285919077D-02  /
   data green(29, 6, 0) /   3.3775347463901692D-02  /
   data green(29, 6, 1) /   3.3756048113527920D-02  /
   data green(29, 6, 2) /   3.3698349387871807D-02  /
   data green(29, 6, 3) /   3.3602843508551854D-02  /
   data green(29, 6, 4) /   3.3470498657454960D-02  /
   data green(29, 6, 5) /   3.3302631813757116D-02  /
   data green(29, 6, 6) /   3.3100872853976872D-02  /
   data green(29, 7, 0) /   3.3527053646395008D-02  /
   data green(29, 7, 1) /   3.3508178658912285D-02  /
   data green(29, 7, 2) /   3.3451745724628834D-02  /
   data green(29, 7, 3) /   3.3358325468488999D-02  /
   data green(29, 7, 4) /   3.3228851040288249D-02  /
   data green(29, 7, 5) /   3.3064592321745336D-02  /
   data green(29, 7, 6) /   3.2867121806358819D-02  /
   data green(29, 7, 7) /   3.2638273923933668D-02  /
   data green(29, 8, 0) /   3.3247297164386104D-02  /
   data green(29, 8, 1) /   3.3228892661164439D-02  /
   data green(29, 8, 2) /   3.3173863221390161D-02  /
   data green(29, 8, 3) /   3.3082755911757071D-02  /
   data green(29, 8, 4) /   3.2956465650653459D-02  /
   data green(29, 8, 5) /   3.2796210886255811D-02  /
   data green(29, 8, 6) /   3.2603501389863211D-02  /
   data green(29, 8, 7) /   3.2380099809966749D-02  /
   data green(29, 8, 8) /   3.2127978849119103D-02  /
   data green(29, 9, 0) /   3.2938601778786998D-02  /
   data green(29, 9, 1) /   3.2920707114944506D-02  /
   data green(29, 9, 2) /   3.2867198726999218D-02  /
   data green(29, 9, 3) /   3.2778598613190969D-02  /
   data green(29, 9, 4) /   3.2655760989285576D-02  /
   data green(29, 9, 5) /   3.2499849504480822D-02  /
   data green(29, 9, 6) /   3.2312307037973310D-02  /
   data green(29, 9, 7) /   3.2094819591654625D-02  /
   data green(29, 9, 8) /   3.1849275998417394D-02  /
   data green(29, 9, 9) /   3.1577725238873040D-02  /
   data green(29,10, 0) /   3.2603617531886613D-02  /
   data green(29,10, 1) /   3.2586265195074356D-02  /
   data green(29,10, 2) /   3.2534374963498966D-02  /
   data green(29,10, 3) /   3.2448442698615833D-02  /
   data green(29,10, 4) /   3.2329280156055386D-02  /
   data green(29,10, 5) /   3.2177993773893521D-02  /
   data green(29,10, 6) /   3.1995956520623724D-02  /
   data green(29,10, 7) /   3.1784774187680463D-02  /
   data green(29,10, 8) /   3.1546247702239745D-02  /
   data green(29,10, 9) /   3.1282333109225353D-02  /
   data green(29,10,10) /   3.0995100832826469D-02  /
   data green(29,11, 0) /   3.2245065738998303D-02  /
   data green(29,11, 1) /   3.2228281431067742D-02  /
   data green(29,11, 2) /   3.2178086244706257D-02  /
   data green(29,11, 3) /   3.2094949255045295D-02  /
   data green(29,11, 4) /   3.1979638677479459D-02  /
   data green(29,11, 5) /   3.1833202233617221D-02  /
   data green(29,11, 6) /   3.1656941060951273D-02  /
   data green(29,11, 7) /   3.1452378422534427D-02  /
   data green(29,11, 8) /   3.1221224650348336D-02  /
   data green(29,11, 9) /   3.0965339828413780D-02  /
   data green(29,11,10) /   3.0686695693544794D-02  /
   data green(29,11,11) /   3.0387338116686643D-02  /
   data green(29,12, 0) /   3.1865687881853580D-02  /
   data green(29,12, 1) /   3.1849490757861414D-02  /
   data green(29,12, 2) /   3.1801047995375004D-02  /
   data green(29,12, 3) /   3.1720801618783642D-02  /
   data green(29,12, 4) /   3.1609475847127692D-02  /
   data green(29,12, 5) /   3.1468059017531060D-02  /
   data green(29,12, 6) /   3.1297779533638242D-02  /
   data green(29,12, 7) /   3.1100076971018433D-02  /
   data green(29,12, 8) /   3.0876569635323595D-02  /
   data green(29,12, 9) /   3.0629019939798597D-02  /
   data green(29,12,10) /   3.0359298949926578D-02  /
   data green(29,12,11) /   3.0069351345991067D-02  /
   data green(29,12,12) /   2.9761161896269289D-02  /
   data green(29,13, 0) /   3.1468199758094825D-02  /
   data green(29,13, 1) /   3.1452602793324347D-02  /
   data green(29,13, 2) /   3.1405951410239073D-02  /
   data green(29,13, 3) /   3.1328660662055273D-02  /
   data green(29,13, 4) /   3.1221410873721097D-02  /
   data green(29,13, 5) /   3.1085131081241946D-02  /
   data green(29,13, 6) /   3.0920976967925430D-02  /
   data green(29,13, 7) /   3.0730304311137435D-02  /
   data green(29,13, 8) /   3.0514639103518208D-02  /
   data green(29,13, 9) /   3.0275645581217391D-02  /
   data green(29,13,10) /   3.0015093380975539D-02  /
   data green(29,13,11) /   2.9734824967242245D-02  /
   data green(29,13,12) /   2.9436724334483900D-02  /
   data green(29,13,13) /   2.9122687815644494D-02  /
   data green(29,14, 0) /   3.1055251781375862D-02  /
   data green(29,14, 1) /   3.1040262235735593D-02  /
   data green(29,14, 2) /   3.0995424142125305D-02  /
   data green(29,14, 3) /   3.0921125960538278D-02  /
   data green(29,14, 4) /   3.0818004709483551D-02  /
   data green(29,14, 5) /   3.0686930862461530D-02  /
   data green(29,14, 6) /   3.0528988198410979D-02  /
   data green(29,14, 7) /   3.0345449508547218D-02  /
   data green(29,14, 8) /   3.0137749199203027D-02  /
   data green(29,14, 9) /   2.9907453896125393D-02  /
   data green(29,14,10) /   2.9656232151758482D-02  /
   data green(29,14,11) /   2.9385824291046126D-02  /
   data green(29,14,12) /   2.9098013315371589D-02  /
   data green(29,14,13) /   2.8794597633136167D-02  /
   data green(29,14,14) /   2.8477366214512634D-02  /
   data green(29,15, 0) /   3.0629395893906241D-02  /
   data green(29,15, 1) /   3.0615015843815600D-02  /
   data green(29,15, 2) /   3.0571997483856882D-02  /
   data green(29,15, 3) /   3.0500703309726079D-02  /
   data green(29,15, 4) /   3.0401728030451363D-02  /
   data green(29,15, 5) /   3.0275884850237723D-02  /
   data green(29,15, 6) /   3.0124187142124361D-02  /
   data green(29,15, 7) /   2.9947826312720627D-02  /
   data green(29,15, 8) /   2.9748146781779273D-02  /
   data green(29,15, 9) /   2.9526619062969689D-02  /
   data green(29,15,10) /   2.9284811933864015D-02  /
   data green(29,15,11) /   2.9024364630022129D-02  /
   data green(29,15,12) /   2.8746959900255276D-02  /
   data green(29,15,13) /   2.8454298630103909D-02  /
   data green(29,15,14) /   2.8148076591460409D-02  /
   data green(29,15,15) /   2.7829963720605717D-02  /
   data green(29,16, 0) /   3.0193059174798711D-02  /
   data green(29,16, 1) /   3.0179286084994057D-02  /
   data green(29,16, 2) /   3.0138080136452569D-02  /
   data green(29,16, 3) /   3.0069778693509401D-02  /
   data green(29,16, 4) /   2.9974935486400910D-02  /
   data green(29,16, 5) /   2.9854308197479220D-02  /
   data green(29,16, 6) /   2.9708841854616974D-02  /
   data green(29,16, 7) /   2.9539648736378009D-02  /
   data green(29,16, 8) /   2.9347985605919061D-02  /
   data green(29,16, 9) /   2.9135229149599064D-02  /
   data green(29,16,10) /   2.8902850502336375D-02  /
   data green(29,16,11) /   2.8652389699760823D-02  /
   data green(29,16,12) /   2.8385430815473272D-02  /
   data green(29,16,13) /   2.8103578430680169D-02  /
   data green(29,16,14) /   2.7808435954335818D-02  /
   data green(29,16,15) /   2.7501586175507849D-02  /
   data green(29,16,16) /   2.7184574295474582D-02  /
   data green(29,17, 0) /   2.9748523919193852D-02  /
   data green(29,17, 1) /   2.9735351230275688D-02  /
   data green(29,17, 2) /   2.9695938351978014D-02  /
   data green(29,17, 3) /   2.9630598509448911D-02  /
   data green(29,17, 4) /   2.9539846046329046D-02  /
   data green(29,17, 5) /   2.9424385228527583D-02  /
   data green(29,17, 6) /   2.9285095245852960D-02  /
   data green(29,17, 7) /   2.9123012029981077D-02  /
   data green(29,17, 8) /   2.8939307608220086D-02  /
   data green(29,17, 9) /   2.8735267767710473D-02  /
   data green(29,17,10) /   2.8512268814121815D-02  /
   data green(29,17,11) /   2.8271754176402191D-02  /
   data green(29,17,12) /   2.8015211541480289D-02  /
   data green(29,17,13) /   2.7744151108703350D-02  /
   data green(29,17,14) /   2.7460085442691629D-02  /
   data green(29,17,15) /   2.7164511284486853D-02  /
   data green(29,17,16) /   2.6858893562640058D-02  /
   data green(29,17,17) /   2.6544651734979878D-02  /
   data green(29,18, 0) /   2.9297913730806899D-02  /
   data green(29,18, 1) /   2.9285331442322798D-02  /
   data green(29,18, 2) /   2.9247682009423136D-02  /
   data green(29,18, 3) /   2.9185255628163303D-02  /
   data green(29,18, 4) /   2.9098529042483251D-02  /
   data green(29,18, 5) /   2.8988155475858611D-02  /
   data green(29,18, 6) /   2.8854951125400062D-02  /
   data green(29,18, 7) /   2.8699878759267911D-02  /
   data green(29,18, 8) /   2.8524029048592390D-02  /
   data green(29,18, 9) /   2.8328600316300846D-02  /
   data green(29,18,10) /   2.8114877397110941D-02  /
   data green(29,18,11) /   2.7884210278404584D-02  /
   data green(29,18,12) /   2.7637993136232237D-02  /
   data green(29,18,13) /   2.7377644301494142D-02  /
   data green(29,18,14) /   2.7104587596357959D-02  /
   data green(29,18,15) /   2.6820235378115603D-02  /
   data green(29,18,16) /   2.6525973524123177D-02  /
   data green(29,18,17) /   2.6223148493179839D-02  /
   data green(29,18,18) /   2.5913056510181585D-02  /
   data green(29,19, 0) /   2.8843185012349653D-02  /
   data green(29,19, 1) /   2.8831180245365864D-02  /
   data green(29,19, 2) /   2.8795256024540368D-02  /
   data green(29,19, 3) /   2.8735680707807108D-02  /
   data green(29,19, 4) /   2.8652895360955579D-02  /
   data green(29,19, 5) /   2.8547504727025844D-02  /
   data green(29,19, 6) /   2.8420265096078669D-02  /
   data green(29,19, 7) /   2.8272069546759312D-02  /
   data green(29,19, 8) /   2.8103931111542264D-02  /
   data green(29,19, 9) /   2.7916964464745727D-02  /
   data green(29,19,10) /   2.7712366745906856D-02  /
   data green(29,19,11) /   2.7491398113144563D-02  /
   data green(29,19,12) /   2.7255362576113256D-02  /
   data green(29,19,13) /   2.7005589591953184D-02  /
   data green(29,19,14) /   2.6743416826918286D-02  /
   data green(29,19,15) /   2.6470174397804341D-02  /
   data green(29,19,16) /   2.6187170817086811D-02  /
   data green(29,19,17) /   2.5895680779041464D-02  /
   data green(29,19,18) /   2.5596934845118102D-02  /
   data green(29,19,19) /   2.5292111018254931D-02  /
   data green(29,20, 0) /   2.8386123146795852D-02  /
   data green(29,20, 1) /   2.8374680673943595D-02  /
   data green(29,20, 2) /   2.8340436402612391D-02  /
   data green(29,20, 3) /   2.8283638092123468D-02  /
   data green(29,20, 4) /   2.8204693133769450D-02  /
   data green(29,20, 5) /   2.8104160469323659D-02  /
   data green(29,20, 6) /   2.7982739721220128D-02  /
   data green(29,20, 7) /   2.7841257944095185D-02  /
   data green(29,20, 8) /   2.7680654478807798D-02  /
   data green(29,20, 9) /   2.7501964433290536D-02  /
   data green(29,20,10) /   2.7306301329066630D-02  /
   data green(29,20,11) /   2.7094839439668266D-02  /
   data green(29,20,12) /   2.6868796311010832D-02  /
   data green(29,20,13) /   2.6629415898823009D-02  /
   data green(29,20,14) /   2.6377952690011801D-02  /
   data green(29,20,15) /   2.6115657098994671D-02  /
   data green(29,20,16) /   2.5843762351841437D-02  /
   data green(29,20,17) /   2.5563472995096610D-02  /
   data green(29,20,18) /   2.5275955096017710D-02  /
   data green(29,20,19) /   2.4982328139243379D-02  /
   data green(29,20,20) /   2.4683658573163486D-02  /
   data green(29,21, 0) /   2.7928342627020323D-02  /
   data green(29,21, 1) /   2.7917445361693279D-02  /
   data green(29,21, 2) /   2.7884830206335844D-02  /
   data green(29,21, 3) /   2.7830725582184215D-02  /
   data green(29,21, 4) /   2.7755507246953318D-02  /
   data green(29,21, 5) /   2.7659691076825508D-02  /
   data green(29,21, 6) /   2.7543923345314039D-02  /
   data green(29,21, 7) /   2.7408968854149743D-02  /
   data green(29,21, 8) /   2.7255697334527873D-02  /
   data green(29,21, 9) /   2.7085068576429373D-02  /
   data green(29,21,10) /   2.6898116758659770D-02  /
   data green(29,21,11) /   2.6695934443949805D-02  /
   data green(29,21,12) /   2.6479656674702316D-02  /
   data green(29,21,13) /   2.6250445559645029D-02  /
   data green(29,21,14) /   2.6009475684304669D-02  /
   data green(29,21,15) /   2.5757920613575718D-02  /
   data green(29,21,16) /   2.5496940687213029D-02  /
   data green(29,21,17) /   2.5227672242777374D-02  /
   data green(29,21,18) /   2.4951218338613874D-02  /
   data green(29,21,19) /   2.4668640994197504D-02  /
   data green(29,21,20) /   2.4380954918152595D-02  /
   data green(29,21,21) /   2.4089122656177717D-02  /
   data green(29,22, 0) /   2.7471290399970021D-02  /
   data green(29,22, 1) /   2.7460919839565539D-02  /
   data green(29,22, 2) /   2.7429878717763447D-02  /
   data green(29,22, 3) /   2.7378377376441963D-02  /
   data green(29,22, 4) /   2.7306761980632335D-02  /
   data green(29,22, 5) /   2.7215508082435051D-02  /
   data green(29,22, 6) /   2.7105211941851462D-02  /
   data green(29,22, 7) /   2.6976579911733224D-02  /
   data green(29,22, 8) /   2.6830416249788540D-02  /
   data green(29,22, 9) /   2.6667609756261857D-02  /
   data green(29,22,10) /   2.6489119650856565D-02  /
   data green(29,22,11) /   2.6295961097564841D-02  /
   data green(29,22,12) /   2.6089190763474844D-02  /
   data green(29,22,13) /   2.5869892760485898D-02  /
   data green(29,22,14) /   2.5639165270905661D-02  /
   data green(29,22,15) /   2.5398108103064504D-02  /
   data green(29,22,16) /   2.5147811365166972D-02  /
   data green(29,22,17) /   2.4889345388002217D-02  /
   data green(29,22,18) /   2.4623751972673521D-02  /
   data green(29,22,19) /   2.4352036990311764D-02  /
   data green(29,22,20) /   2.4075164318242359D-02  /
   data green(29,22,21) /   2.3794051062057973D-02  /
   data green(29,22,22) /   2.3509563985715375D-02  /
   data green(29,23, 0) /   2.7016251732365670D-02  /
   data green(29,23, 1) /   2.7006388353144187D-02  /
   data green(29,23, 2) /   2.6976863111972220D-02  /
   data green(29,23, 3) /   2.6927869509885224D-02  /
   data green(29,23, 4) /   2.6859726129911079D-02  /
   data green(29,23, 5) /   2.6772870906185374D-02  /
   data green(29,23, 6) /   2.6667853386069876D-02  /
   data green(29,23, 7) /   2.6545325250519423D-02  /
   data green(29,23, 8) /   2.6406029406977696D-02  /
   data green(29,23, 9) /   2.6250788001258574D-02  /
   data green(29,23,10) /   2.6080489709522828D-02  /
   data green(29,23,11) /   2.5896076669183001D-02  /
   data green(29,23,12) /   2.5698531390054472D-02  /
   data green(29,23,13) /   2.5488863956832749D-02  /
   data green(29,23,14) /   2.5268099794068348D-02  /
   data green(29,23,15) /   2.5037268218489912D-02  /
   data green(29,23,16) /   2.4797391954000290D-02  /
   data green(29,23,17) /   2.4549477734842885D-02  /
   data green(29,23,18) /   2.4294508074768531D-02  /
   data green(29,23,19) /   2.4033434236433966D-02  /
   data green(29,23,20) /   2.3767170397072664D-02  /
   data green(29,23,21) /   2.3496588974482410D-02  /
   data green(29,23,22) /   2.3222517051872841D-02  /
   data green(29,23,23) /   2.2945733821003952D-02  /
   data green(29,24, 0) /   2.6564357967235171D-02  /
   data green(29,24, 1) /   2.6554981570511321D-02  /
   data green(29,24, 2) /   2.6526912020280548D-02  /
   data green(29,24, 3) /   2.6480327180634929D-02  /
   data green(29,24, 4) /   2.6415520009353332D-02  /
   data green(29,24, 5) /   2.6332893460760648D-02  /
   data green(29,24, 6) /   2.6232953595127236D-02  /
   data green(29,24, 7) /   2.6116301123282920D-02  /
   data green(29,24, 8) /   2.5983621658159203D-02  /
   data green(29,24, 9) /   2.5835674973872881D-02  /
   data green(29,24,10) /   2.5673283587078564D-02  /
   data green(29,24,11) /   2.5497320975036108D-02  /
   data green(29,24,12) /   2.5308699731466447D-02  /
   data green(29,24,13) /   2.5108359936814671D-02  /
   data green(29,24,14) /   2.4897257986482048D-02  /
   data green(29,24,15) /   2.4676356081623871D-02  /
   data green(29,24,16) /   2.4446612544921739D-02  /
   data green(29,24,17) /   2.4208973080804011D-02  /
   data green(29,24,18) /   2.3964363058037914D-02  /
   data green(29,24,19) /   2.3713680854151986D-02  /
   data green(29,24,20) /   2.3457792266996805D-02  /
   data green(29,24,21) /   2.3197525969695411D-02  /
   data green(29,24,22) /   2.2933669961623850D-02  /
   data green(29,24,23) /   2.2666968949911680D-02  /
   data green(29,24,24) /   2.2398122582987331D-02  /
   data green(29,25, 0) /   2.6116595615219812D-02  /
   data green(29,25, 1) /   2.6107685626227697D-02  /
   data green(29,25, 2) /   2.6081010433449615D-02  /
   data green(29,25, 3) /   2.6036733422461777D-02  /
   data green(29,25, 4) /   2.5975123810659834D-02  /
   data green(29,25, 5) /   2.5896552117894541D-02  /
   data green(29,25, 6) /   2.5801484036153552D-02  /
   data green(29,25, 7) /   2.5690472895195222D-02  /
   data green(29,25, 8) /   2.5564150958742485D-02  /
   data green(29,25, 9) /   2.5423219811705287D-02  /
   data green(29,25,10) /   2.5268440112309824D-02  /
   data green(29,25,11) /   2.5100620984206526D-02  /
   data green(29,25,12) /   2.4920609313607199D-02  /
   data green(29,25,13) /   2.4729279196857430D-02  /
   data green(29,25,14) /   2.4527521756594551D-02  /
   data green(29,25,15) /   2.4316235511989859D-02  /
   data green(29,25,16) /   2.4096317452767267D-02  /
   data green(29,25,17) /   2.3868654929827558D-02  /
   data green(29,25,18) /   2.3634118439229957D-02  /
   data green(29,25,19) /   2.3393555342485113D-02  /
   data green(29,25,20) /   2.3147784535723002D-02  /
   data green(29,25,21) /   2.2897592054054779D-02  /
   data green(29,25,22) /   2.2643727575744184D-02  /
   data green(29,25,23) /   2.2386901773730591D-02  /
   data green(29,25,24) /   2.2127784449447213D-02  /
   data green(29,25,25) /   2.1867003375411399D-02  /
   data green(29,26, 0) /   2.5673816304303231D-02  /
   data green(29,26, 1) /   2.5665352026249630D-02  /
   data green(29,26, 2) /   2.5640009473199887D-02  /
   data green(29,26, 3) /   2.5597938657380259D-02  /
   data green(29,26, 4) /   2.5539386855801753D-02  /
   data green(29,26, 5) /   2.5464694588222935D-02  /
   data green(29,26, 6) /   2.5374290167207256D-02  /
   data green(29,26, 7) /   2.5268682989666020D-02  /
   data green(29,26, 8) /   2.5148455772247336D-02  /
   data green(29,26, 9) /   2.5014255955998672D-02  /
   data green(29,26,10) /   2.4866786518326141D-02  /
   data green(29,26,11) /   2.4706796432524911D-02  /
   data green(29,26,12) /   2.4535071007817530D-02  /
   data green(29,26,13) /   2.4352422327170189D-02  /
   data green(29,26,14) /   2.4159679977785443D-02  /
   data green(29,26,15) /   2.3957682241907474D-02  /
   data green(29,26,16) /   2.3747267885289247D-02  /
   data green(29,26,17) /   2.3529268649128000D-02  /
   data green(29,26,18) /   2.3304502520059157D-02  /
   data green(29,26,19) /   2.3073767823218985D-02  /
   data green(29,26,20) /   2.2837838156456548D-02  /
   data green(29,26,21) /   2.2597458160194893D-02  /
   data green(29,26,22) /   2.2353340097617583D-02  /
   data green(29,26,23) /   2.2106161203935690D-02  /
   data green(29,26,24) /   2.1856561751404317D-02  /
   data green(29,26,25) /   2.1605143768270926D-02  /
   data green(29,26,26) /   2.1352470344597488D-02  /
   data green(29,27, 0) /   2.5236747190932027D-02  /
   data green(29,27, 1) /   2.5228708017533397D-02  /
   data green(29,27, 2) /   2.5204636638140446D-02  /
   data green(29,27, 3) /   2.5164670738337322D-02  /
   data green(29,27, 4) /   2.5109037361130026D-02  /
   data green(29,27, 5) /   2.5038049337196090D-02  /
   data green(29,27, 6) /   2.4952100442399868D-02  /
   data green(29,27, 7) /   2.4851659428194381D-02  /
   data green(29,27, 8) /   2.4737263099317040D-02  /
   data green(29,27, 9) /   2.4609508633704275D-02  /
   data green(29,27,10) /   2.4469045351271698D-02  /
   data green(29,27,11) /   2.4316566141180072D-02  /
   data green(29,27,12) /   2.4152798751999535D-02  /
   data green(29,27,13) /   2.3978497136787794D-02  /
   data green(29,27,14) /   2.3794433026815717D-02  /
   data green(29,27,15) /   2.3601387884991547D-02  /
   data green(29,27,16) /   2.3400145364496017D-02  /
   data green(29,27,17) /   2.3191484371236127D-02  /
   data green(29,27,18) /   2.2976172801801985D-02  /
   data green(29,27,19) /   2.2754962002814939D-02  /
   data green(29,27,20) /   2.2528581973788674D-02  /
   data green(29,27,21) /   2.2297737314537356D-02  /
   data green(29,27,22) /   2.2063103900163045D-02  /
   data green(29,27,23) /   2.1825326251922680D-02  /
   data green(29,27,24) /   2.1585015560809297D-02  /
   data green(29,27,25) /   2.1342748312331756D-02  /
   data green(29,27,26) /   2.1099065455479807D-02  /
   data green(29,27,27) /   2.0854472055883956D-02  /
   data green(29,28, 0) /   2.4806001510609361D-02  /
   data green(29,28, 1) /   2.4798367100868628D-02  /
   data green(29,28, 2) /   2.4775506204044388D-02  /
   data green(29,28, 3) /   2.4737545164278507D-02  /
   data green(29,28, 4) /   2.4684692398090621D-02  /
   data green(29,28, 5) /   2.4617235227072043D-02  /
   data green(29,28, 6) /   2.4535535576691650D-02  /
   data green(29,28, 7) /   2.4440024666338784D-02  /
   data green(29,28, 8) /   2.4331196840835699D-02  /
   data green(29,28, 9) /   2.4209602711855377D-02  /
   data green(29,28,10) /   2.4075841788503697D-02  /
   data green(29,28,11) /   2.3930554779759073D-02  /
   data green(29,28,12) /   2.3774415747928247D-02  /
   data green(29,28,13) /   2.3608124282551320D-02  /
   data green(29,28,14) /   2.3432397849315233D-02  /
   data green(29,28,15) /   2.3247964449726297D-02  /
   data green(29,28,16) /   2.3055555705821650D-02  /
   data green(29,28,17) /   2.2855900461314874D-02  /
   data green(29,28,18) /   2.2649718967422777D-02  /
   data green(29,28,19) /   2.2437717699195895D-02  /
   data green(29,28,20) /   2.2220584827278084D-02  /
   data green(29,28,21) /   2.1998986351244400D-02  /
   data green(29,28,22) /   2.1773562884405902D-02  /
   data green(29,28,23) /   2.1544927066431203D-02  /
   data green(29,28,24) /   2.1313661569362787D-02  /
   data green(29,28,25) /   2.1080317654513923D-02  /
   data green(29,28,26) /   2.0845414232131883D-02  /
   data green(29,28,27) /   2.0609437372349354D-02  /
   data green(29,28,28) /   2.0372840214517803D-02  /
   data green(29,29, 0) /   2.4382089014518094D-02  /
   data green(29,29, 1) /   2.4374839433704845D-02  /
   data green(29,29, 2) /   2.4353129525869582D-02  /
   data green(29,29, 3) /   2.4317075216093196D-02  /
   data green(29,29, 4) /   2.4266867800675763D-02  /
   data green(29,29, 5) /   2.4202771137352905D-02  /
   data green(29,29, 6) /   2.4125117825635894D-02  /
   data green(29,29, 7) /   2.4034304484761911D-02  /
   data green(29,29, 8) /   2.3930786258616617D-02  /
   data green(29,29, 9) /   2.3815070693107600D-02  /
   data green(29,29,10) /   2.3687711141393994D-02  /
   data green(29,29,11) /   2.3549299856071126D-02  /
   data green(29,29,12) /   2.3400460925174005D-02  /
   data green(29,29,13) /   2.3241843201306376D-02  /
   data green(29,29,14) /   2.3074113361158721D-02  /
   data green(29,29,15) /   2.2897949217130273D-02  /
   data green(29,29,16) /   2.2714033384769822D-02  /
   data green(29,29,17) /   2.2523047390338075D-02  /
   data green(29,29,18) /   2.2325666282944853D-02  /
   data green(29,29,19) /   2.2122553796280954D-02  /
   data green(29,29,20) /   2.1914358086652080D-02  /
   data green(29,29,21) /   2.1701708057370421D-02  /
   data green(29,29,22) /   2.1485210264942715D-02  /
   data green(29,29,23) /   2.1265446390128062D-02  /
   data green(29,29,24) /   2.1042971246903356D-02  /
   data green(29,29,25) /   2.0818311294631277D-02  /
   data green(29,29,26) /   2.0591963613147482D-02  /
   data green(29,29,27) /   2.0364395296876125D-02  /
   data green(29,29,28) /   2.0136043222209968D-02  /
   data green(29,29,29) /   1.9907314141995373D-02  /
   data green(30, 0, 0) /   3.3342621976882383D-02  /
   data green(30, 1, 0) /   3.3324051657747863D-02  /
   data green(30, 1, 1) /   3.3305512506889727D-02  /
   data green(30, 2, 0) /   3.3268527711078681D-02  /
   data green(30, 2, 1) /   3.3250081538351475D-02  /
   data green(30, 2, 2) /   3.3194927936686658D-02  /
   data green(30, 3, 0) /   3.3176605901982847D-02  /
   data green(30, 3, 1) /   3.3158312958105694D-02  /
   data green(30, 3, 2) /   3.3103616471725633D-02  /
   data green(30, 3, 3) /   3.3013058396271146D-02  /
   data green(30, 4, 0) /   3.3049195242101578D-02  /
   data green(30, 4, 1) /   3.3031113247579798D-02  /
   data green(30, 4, 2) /   3.2977046093934856D-02  /
   data green(30, 4, 3) /   3.2887525326492056D-02  /
   data green(30, 4, 4) /   3.2763420642950479D-02  /
   data green(30, 5, 0) /   3.2887533094014848D-02  /
   data green(30, 5, 1) /   3.2869716365438770D-02  /
   data green(30, 5, 2) /   3.2816440629550157D-02  /
   data green(30, 5, 3) /   3.2728224461625301D-02  /
   data green(30, 5, 4) /   3.2605916497894877D-02  /
   data green(30, 5, 5) /   3.2450672832216872D-02  /
   data green(30, 6, 0) /   3.2693152222263641D-02  /
   data green(30, 6, 1) /   3.2675650925256323D-02  /
   data green(30, 6, 2) /   3.2623316334401159D-02  /
   data green(30, 6, 3) /   3.2536651773733244D-02  /
   data green(30, 6, 4) /   3.2416481106189936D-02  /
   data green(30, 6, 5) /   3.2263927051631711D-02  /
   data green(30, 6, 6) /   3.2080382423930597D-02  /
   data green(30, 7, 0) /   3.2467841488555985D-02  /
   data green(30, 7, 1) /   3.2450701038263870D-02  /
   data green(30, 7, 2) /   3.2399443174980891D-02  /
   data green(30, 7, 3) /   3.2314554005202785D-02  /
   data green(30, 7, 4) /   3.2196829414742507D-02  /
   data green(30, 7, 5) /   3.2047354413505660D-02  /
   data green(30, 7, 6) /   3.1867475712933691D-02  /
   data green(30, 7, 7) /   3.1658768877229800D-02  /
   data green(30, 8, 0) /   3.2213602109358516D-02  /
   data green(30, 8, 1) /   3.2196862730230341D-02  /
   data green(30, 8, 2) /   3.2146801718687251D-02  /
   data green(30, 8, 3) /   3.2063886333176951D-02  /
   data green(30, 8, 4) /   3.1948881820719859D-02  /
   data green(30, 8, 5) /   3.1802831869099282D-02  /
   data green(30, 8, 6) /   3.1627032631608723D-02  /
   data green(30, 8, 7) /   3.1423001575354430D-02  /
   data green(30, 8, 8) /   3.1192442580535114D-02  /
   data green(30, 9, 0) /   3.1932601465359534D-02  /
   data green(30, 9, 1) /   3.1916297912893396D-02  /
   data green(30, 9, 2) /   3.1867537586767115D-02  /
   data green(30, 9, 3) /   3.1786767610622900D-02  /
   data green(30, 9, 4) /   3.1674720486587059D-02  /
   data green(30, 9, 5) /   3.1532395710927724D-02  /
   data green(30, 9, 6) /   3.1361035321273788D-02  /
   data green(30, 9, 7) /   3.1162094533145633D-02  /
   data green(30, 9, 8) /   3.0937208789887545D-02  /
   data green(30, 9, 9) /   3.0688158620828965D-02  /
   data green(30,10, 0) /   3.1627126383994661D-02  /
   data green(30,10, 1) /   3.1611287822898258D-02  /
   data green(30,10, 2) /   3.1563915358469682D-02  /
   data green(30,10, 3) /   3.1485435032636407D-02  /
   data green(30,10, 4) /   3.1376545046405496D-02  /
   data green(30,10, 5) /   3.1238198571535220D-02  /
   data green(30,10, 6) /   3.1071580864935406D-02  /
   data green(30,10, 7) /   3.0878081749286291D-02  /
   data green(30,10, 8) /   3.0659264679178806D-02  /
   data green(30,10, 9) /   3.0416833681329008D-02  /
   data green(30,10,10) /   3.0152599442967883D-02  /
   data green(30,11, 0) /   3.1299537631281849D-02  /
   data green(30,11, 1) /   3.1284187657643361D-02  /
   data green(30,11, 2) /   3.1238273635457500D-02  /
   data green(30,11, 3) /   3.1162199903741892D-02  /
   data green(30,11, 4) /   3.1056629333833720D-02  /
   data green(30,11, 5) /   3.0922467347987242D-02  /
   data green(30,11, 6) /   3.0760840617525292D-02  /
   data green(30,11, 7) /   3.0573071410635987D-02  /
   data green(30,11, 8) /   3.0360648704940811D-02  /
   data green(30,11, 9) /   3.0125197247261183D-02  /
   data green(30,11,10) /   2.9868445734698726D-02  /
   data green(30,11,11) /   2.9592195216017053D-02  /
   data green(30,12, 0) /   3.0952227078890077D-02  /
   data green(30,12, 1) /   3.0937383868829622D-02  /
   data green(30,12, 2) /   3.0892982712011029D-02  /
   data green(30,12, 3) /   3.0819405926682069D-02  /
   data green(30,12, 4) /   3.0717280518079766D-02  /
   data green(30,12, 5) /   3.0587463395656274D-02  /
   data green(30,12, 6) /   3.0431021646457751D-02  /
   data green(30,12, 7) /   3.0249208743796070D-02  /
   data green(30,12, 8) /   3.0043437704538229D-02  /
   data green(30,12, 9) /   2.9815252273255186D-02  /
   data green(30,12,10) /   2.9566297208599226D-02  /
   data green(30,12,11) /   2.9298288683994547D-02  /
   data green(30,12,12) /   2.9012985702757493D-02  /
   data green(30,13, 0) /   3.0587578691116134D-02  /
   data green(30,13, 1) /   3.0573255255093867D-02  /
   data green(30,13, 2) /   3.0530405982571479D-02  /
   data green(30,13, 3) /   3.0459391128016215D-02  /
   data green(30,13, 4) /   3.0360801740755713D-02  /
   data green(30,13, 5) /   3.0235446056093837D-02  /
   data green(30,13, 6) /   3.0084331313804353D-02  /
   data green(30,13, 7) /   2.9908641795916333D-02  /
   data green(30,13, 8) /   2.9709713999129084D-02  /
   data green(30,13, 9) /   2.9489009919342224D-02  /
   data green(30,13,10) /   2.9248089427585746D-02  /
   data green(30,13,11) /   2.8988582664189563D-02  /
   data green(30,13,12) /   2.8712163281336259D-02  /
   data green(30,13,13) /   2.8420523235498381D-02  /
   data green(30,14, 0) /   3.0207934137334982D-02  /
   data green(30,14, 1) /   3.0194138657748205D-02  /
   data green(30,14, 2) /   3.0152865886005506D-02  /
   data green(30,14, 3) /   3.0084454212146255D-02  /
   data green(30,14, 4) /   2.9989459034331215D-02  /
   data green(30,14, 5) /   2.9868640285581956D-02  /
   data green(30,14, 6) /   2.9722945749405760D-02  /
   data green(30,14, 7) /   2.9553490873684210D-02  /
   data green(30,14, 8) /   2.9361535905104735D-02  /
   data green(30,14, 9) /   2.9148461225469883D-02  /
   data green(30,14,10) /   2.8915741776900414D-02  /
   data green(30,14,11) /   2.8664921420244873D-02  /
   data green(30,14,12) /   2.8397587988341893D-02  /
   data green(30,14,13) /   2.8115349683706832D-02  /
   data green(30,14,14) /   2.7819813340058869D-02  /
   data green(30,15, 0) /   2.9815563506787431D-02  /
   data green(30,15, 1) /   2.9802299736479822D-02  /
   data green(30,15, 2) /   2.9762614863356754D-02  /
   data green(30,15, 3) /   2.9696825820201875D-02  /
   data green(30,15, 4) /   2.9605452997779610D-02  /
   data green(30,15, 5) /   2.9489208858301572D-02  /
   data green(30,15, 6) /   2.9348982686041931D-02  /
   data green(30,15, 7) /   2.9185822107313387D-02  /
   data green(30,15, 8) /   2.9000912114753308D-02  /
   data green(30,15, 9) /   2.8795552386579298D-02  /
   data green(30,15,10) /   2.8571133700270163D-02  /
   data green(30,15,11) /   2.8329114206023746D-02  /
   data green(30,15,12) /   2.8070996255374885D-02  /
   data green(30,15,13) /   2.7798304383501119D-02  /
   data green(30,15,14) /   2.7512564929760855D-02  /
   data green(30,15,15) /   2.7215287659424294D-02  /
   data green(30,16, 0) /   2.9412641305765323D-02  /
   data green(30,16, 1) /   2.9399909006286702D-02  /
   data green(30,16, 2) /   2.9361811513919703D-02  /
   data green(30,16, 3) /   2.9298644884259192D-02  /
   data green(30,16, 4) /   2.9210895427203763D-02  /
   data green(30,16, 5) /   2.9099229350367649D-02  /
   data green(30,16, 6) /   2.8964478871934079D-02  /
   data green(30,16, 7) /   2.8807625363778844D-02  /
   data green(30,16, 8) /   2.8629780178727297D-02  /
   data green(30,16, 9) /   2.8432163867997108D-02  /
   data green(30,16,10) /   2.8216084506066302D-02  /
   data green(30,16,11) /   2.7982915813568831D-02  /
   data green(30,16,12) /   2.7734075710188103D-02  /
   data green(30,16,13) /   2.7471005846459506D-02  /
   data green(30,16,14) /   2.7195152564262103D-02  /
   data green(30,16,15) /   2.6907949628835717D-02  /
   data green(30,16,16) /   2.6610802967860763D-02  /
   data green(30,17, 0) /   2.9001227666388983D-02  /
   data green(30,17, 1) /   2.8989023066947894D-02  /
   data green(30,17, 2) /   2.8952501887134767D-02  /
   data green(30,17, 3) /   2.8891940024381783D-02  /
   data green(30,17, 4) /   2.8807790862211861D-02  /
   data green(30,17, 5) /   2.8700675881437078D-02  /
   data green(30,17, 6) /   2.8571372055136993D-02  /
   data green(30,17, 7) /   2.8420796522839339D-02  /
   data green(30,17, 8) /   2.8249989123188532D-02  /
   data green(30,17, 9) /   2.8060093412958592D-02  /
   data green(30,17,10) /   2.7852336813183897D-02  /
   data green(30,17,11) /   2.7628010502927260D-02  /
   data green(30,17,12) /   2.7388449632561700D-02  /
   data green(30,17,13) /   2.7135014357745447D-02  /
   data green(30,17,14) /   2.6869072109618985D-02  /
   data green(30,17,15) /   2.6591981423262672D-02  /
   data green(30,17,16) /   2.6305077551643640D-02  /
   data green(30,17,17) /   2.6009660001635567D-02  /
   data green(30,18, 0) /   2.8583254498254020D-02  /
   data green(30,18, 1) /   2.8571570758963116D-02  /
   data green(30,18, 2) /   2.8536605652105917D-02  /
   data green(30,18, 3) /   2.8478615743027385D-02  /
   data green(30,18, 4) /   2.8398022819686848D-02  /
   data green(30,18, 5) /   2.8295405406311192D-02  /
   data green(30,18, 6) /   2.8171487355905304D-02  /
   data green(30,18, 7) /   2.8027123957616560D-02  /
   data green(30,18, 8) /   2.7863286070229108D-02  /
   data green(30,18, 9) /   2.7681042837986884D-02  /
   data green(30,18,10) /   2.7481543559011688D-02  /
   data green(30,18,11) /   2.7265999261696976D-02  /
   data green(30,18,12) /   2.7035664504512821D-02  /
   data green(30,18,13) /   2.6791819854918183D-02  /
   data green(30,18,14) /   2.6535755429532348D-02  /
   data green(30,18,15) /   2.6268755796463213D-02  /
   data green(30,18,16) /   2.5992086457381710D-02  /
   data green(30,18,17) /   2.5706982046422924D-02  /
   data green(30,18,18) /   2.5414636309082812D-02  /
   data green(30,19, 0) /   2.8160516170344031D-02  /
   data green(30,19, 1) /   2.8149343836279647D-02  /
   data green(30,19, 2) /   2.8115906743710940D-02  /
   data green(30,19, 3) /   2.8060443029865025D-02  /
   data green(30,19, 4) /   2.7983344346986704D-02  /
   data green(30,19, 5) /   2.7885148211821201D-02  /
   data green(30,19, 6) /   2.7766527709050914D-02  /
   data green(30,19, 7) /   2.7628278931098527D-02  /
   data green(30,19, 8) /   2.7471306604019766D-02  /
   data green(30,19, 9) /   2.7296608390521927D-02  /
   data green(30,19,10) /   2.7105258375860350D-02  /
   data green(30,19,11) /   2.6898390231927572D-02  /
   data green(30,19,12) /   2.6677180522392166D-02  /
   data green(30,19,13) /   2.6442832561617467D-02  /
   data green(30,19,14) /   2.6196561177318764D-02  /
   data green(30,19,15) /   2.5939578656681152D-02  /
   data green(30,19,16) /   2.5673082082847801D-02  /
   data green(30,19,17) /   2.5398242197554103D-02  /
   data green(30,19,18) /   2.5116193859636367D-02  /
   data green(30,19,19) /   2.4828028110653835D-02  /
   data green(30,20, 0) /   2.7734664217392440D-02  /
   data green(30,20, 1) /   2.7723991652931757D-02  /
   data green(30,20, 2) /   2.7692047991164488D-02  /
   data green(30,20, 3) /   2.7639053897104481D-02  /
   data green(30,20, 4) /   2.7565372433975915D-02  /
   data green(30,20, 5) /   2.7471502182152548D-02  /
   data green(30,20, 6) /   2.7358067966947543D-02  /
   data green(30,20, 7) /   2.7225809529672382D-02  /
   data green(30,20, 8) /   2.7075568536357304D-02  /
   data green(30,20, 9) /   2.6908274356314786D-02  /
   data green(30,20,10) /   2.6724929057680355D-02  /
   data green(30,20,11) /   2.6526592060249973D-02  /
   data green(30,20,12) /   2.6314364859853379D-02  /
   data green(30,20,13) /   2.6089376196732468D-02  /
   data green(30,20,14) /   2.5852767987117480D-02  /
   data green(30,20,15) /   2.5605682276798894D-02  /
   data green(30,20,16) /   2.5349249412167513D-02  /
   data green(30,20,17) /   2.5084577561656617D-02  /
   data green(30,20,18) /   2.4812743661798885D-02  /
   data green(30,20,19) /   2.4534785809468950D-02  /
   data green(30,20,20) /   2.4251697076800961D-02  /
   data green(30,21, 0) /   2.7307205516031142D-02  /
   data green(30,21, 1) /   2.7297019311717026D-02  /
   data green(30,21, 2) /   2.7266529185435082D-02  /
   data green(30,21, 3) /   2.7215939315271380D-02  /
   data green(30,21, 4) /   2.7145585772223722D-02  /
   data green(30,21, 5) /   2.7055930343704553D-02  /
   data green(30,21, 6) /   2.6947552200780749D-02  /
   data green(30,21, 7) /   2.6821137700837876D-02  /
   data green(30,21, 8) /   2.6677468670622966D-02  /
   data green(30,21, 9) /   2.6517409549014888D-02  /
   data green(30,21,10) /   2.6341893783728230D-02  /
   data green(30,21,11) /   2.6151909872239879D-02  /
   data green(30,21,12) /   2.5948487416527478D-02  /
   data green(30,21,13) /   2.5732683526637181D-02  /
   data green(30,21,14) /   2.5505569863126260D-02  /
   data green(30,21,15) /   2.5268220556746383D-02  /
   data green(30,21,16) /   2.5021701188918746D-02  /
   data green(30,21,17) /   2.4767058961825412D-02  /
   data green(30,21,18) /   2.4505314134982741D-02  /
   data green(30,21,19) /   2.4237452758024098D-02  /
   data green(30,21,20) /   2.3964420688500174D-02  /
   data green(30,21,21) /   2.3687118849588181D-02  /
   data green(30,22, 0) /   2.6879503363717826D-02  /
   data green(30,22, 1) /   2.6869788710410632D-02  /
   data green(30,22, 2) /   2.6840708028481573D-02  /
   data green(30,22, 3) /   2.6792450004621257D-02  /
   data green(30,22, 4) /   2.6725325333246282D-02  /
   data green(30,22, 5) /   2.6639761182088276D-02  /
   data green(30,22, 6) /   2.6536293716910613D-02  /
   data green(30,22, 7) /   2.6415558939195633D-02  /
   data green(30,22, 8) /   2.6278282137818677D-02  /
   data green(30,22, 9) /   2.6125266286890582D-02  /
   data green(30,22,10) /   2.5957379736449841D-02  /
   data green(30,22,11) /   2.5775543541045734D-02  /
   data green(30,22,12) /   2.5580718755046077D-02  /
   data green(30,22,13) /   2.5373893995081928D-02  /
   data green(30,22,14) /   2.5156073532279142D-02  /
   data green(30,22,15) /   2.4928266132906503D-02  /
   data green(30,22,16) /   2.4691474818832080D-02  /
   data green(30,22,17) /   2.4446687671499052D-02  /
   data green(30,22,18) /   2.4194869757363013D-02  /
   data green(30,22,19) /   2.3936956210724350D-02  /
   data green(30,22,20) /   2.3673846472942874D-02  /
   data green(30,22,21) /   2.3406399655927283D-02  /
   data green(30,22,22) /   2.3135430972863791D-02  /
   data green(30,23, 0) /   2.6452780909276309D-02  /
   data green(30,23, 1) /   2.6443521936515130D-02  /
   data green(30,23, 2) /   2.6415803421747520D-02  /
   data green(30,23, 3) /   2.6369799550208333D-02  /
   data green(30,23, 4) /   2.6305797248288122D-02  /
   data green(30,23, 5) /   2.6224191231823304D-02  /
   data green(30,23, 6) /   2.6125477310170207D-02  /
   data green(30,23, 7) /   2.6010244166544922D-02  /
   data green(30,23, 8) /   2.5879163876780201D-02  /
   data green(30,23, 9) /   2.5732981456776226D-02  /
   data green(30,23,10) /   2.5572503742865991D-02  /
   data green(30,23,11) /   2.5398587909428994D-02  /
   data green(30,23,12) /   2.5212129915590697D-02  /
   data green(30,23,13) /   2.5014053149644066D-02  /
   data green(30,23,14) /   2.4805297508276477D-02  /
   data green(30,23,15) /   2.4586809110354997D-02  /
   data green(30,23,16) /   2.4359530804481395D-02  /
   data green(30,23,17) /   2.4124393588157749D-02  /
   data green(30,23,18) /   2.3882309016256632D-02  /
   data green(30,23,19) /   2.3634162639220687D-02  /
   data green(30,23,20) /   2.3380808478222624D-02  /
   data green(30,23,21) /   2.3123064516171601D-02  /
   data green(30,23,22) /   2.2861709160319963D-02  /
   data green(30,23,23) /   2.2597478614336912D-02  /
   data green(30,24, 0) /   2.6028126419957508D-02  /
   data green(30,24, 1) /   2.6019306497243176D-02  /
   data green(30,24, 2) /   2.5992900585948577D-02  /
   data green(30,24, 3) /   2.5949069341427697D-02  /
   data green(30,24, 4) /   2.5888077502527466D-02  /
   data green(30,24, 5) /   2.5810289466740993D-02  /
   data green(30,24, 6) /   2.5716163300005423D-02  /
   data green(30,24, 7) /   2.5606243372338210D-02  /
   data green(30,24, 8) /   2.5481151847229194D-02  /
   data green(30,24, 9) /   2.5341579277994387D-02  /
   data green(30,24,10) /   2.5188274577547042D-02  /
   data green(30,24,11) /   2.5022034629466151D-02  /
   data green(30,24,12) /   2.4843693798787939D-02  /
   data green(30,24,13) /   2.4654113582137937D-02  /
   data green(30,24,14) /   2.4454172610586713D-02  /
   data green(30,24,15) /   2.4244757187084430D-02  /
   data green(30,24,16) /   2.4026752505670924D-02  /
   data green(30,24,17) /   2.3801034663898529D-02  /
   data green(30,24,18) /   2.3568463544830335D-02  /
   data green(30,24,19) /   2.3329876612049952D-02  /
   data green(30,24,20) /   2.3086083631435158D-02  /
   data green(30,24,21) /   2.2837862307747332D-02  /
   data green(30,24,22) /   2.2585954802763462D-02  /
   data green(30,24,23) /   2.2331065084832884D-02  /
   data green(30,24,24) /   2.2073857047237511D-02  /
   data green(30,25, 0) /   2.5606499919160666D-02  /
   data green(30,25, 1) /   2.5598101920294766D-02  /
   data green(30,25, 2) /   2.5572957551955728D-02  /
   data green(30,25, 3) /   2.5531214882833948D-02  /
   data green(30,25, 4) /   2.5473118000162508D-02  /
   data green(30,25, 5) /   2.5399003059711916D-02  /
   data green(30,25, 6) /   2.5309292932553559D-02  /
   data green(30,25, 7) /   2.5204490614164042D-02  /
   data green(30,25, 8) /   2.5085171593743689D-02  /
   data green(30,25, 9) /   2.4951975404285554D-02  /
   data green(30,25,10) /   2.4805596586393423D-02  /
   data green(30,25,11) /   2.4646775301222695D-02  /
   data green(30,25,12) /   2.4476287820929928D-02  /
   data green(30,25,13) /   2.4294937109883893D-02  /
   data green(30,25,14) /   2.4103543688180920D-02  /
   data green(30,25,15) /   2.3902936942486186D-02  /
   data green(30,25,16) /   2.3693947019699357D-02  /
   data green(30,25,17) /   2.3477397408144640D-02  /
   data green(30,25,18) /   2.3254098280455996D-02  /
   data green(30,25,19) /   2.3024840643350902D-02  /
   data green(30,25,20) /   2.2790391313057691D-02  /
   data green(30,25,21) /   2.2551488711974009D-02  /
   data green(30,25,22) /   2.2308839462594193D-02  /
   data green(30,25,23) /   2.2063115739005015D-02  /
   data green(30,25,24) /   2.1814953324255336D-02  /
   data green(30,25,25) /   2.1564950313432491D-02  /
   data green(30,26, 0) /   2.5188740785422631D-02  /
   data green(30,26, 1) /   2.5180747317113152D-02  /
   data green(30,26, 2) /   2.5156812617681636D-02  /
   data green(30,26, 3) /   2.5117073076477692D-02  /
   data green(30,26, 4) /   2.5061753608008396D-02  /
   data green(30,26, 5) /   2.4991164128689999D-02  /
   data green(30,26, 6) /   2.4905694776883126D-02  /
   data green(30,26, 7) /   2.4805810019439836D-02  /
   data green(30,26, 8) /   2.4692041816353136D-02  /
   data green(30,26, 9) /   2.4564982035349119D-02  /
   data green(30,26,10) /   2.4425274319884348D-02  /
   data green(30,26,11) /   2.4273605617043826D-02  /
   data green(30,26,12) /   2.4110697566827096D-02  /
   data green(30,26,13) /   2.3937297942231223D-02  /
   data green(30,26,14) /   2.3754172311657752D-02  /
   data green(30,26,15) /   2.3562096072944311D-02  /
   data green(30,26,16) /   2.3361846983258397D-02  /
   data green(30,26,17) /   2.3154198282654807D-02  /
   data green(30,26,18) /   2.2939912482613699D-02  /
   data green(30,26,19) /   2.2719735865470443D-02  /
   data green(30,26,20) /   2.2494393717211673D-02  /
   data green(30,26,21) /   2.2264586295288855D-02  /
   data green(30,26,22) /   2.2030985515297218D-02  /
   data green(30,26,23) /   2.1794232325767913D-02  /
   data green(30,26,24) /   2.1554934728932270D-02  /
   data green(30,26,25) /   2.1313666396991032D-02  /
   data green(30,26,26) /   2.1070965827907209D-02  /
   data green(30,27, 0) /   2.4775575962303415D-02  /
   data green(30,27, 1) /   2.4767969559021698D-02  /
   data green(30,27, 2) /   2.4745192423687040D-02  /
   data green(30,27, 3) /   2.4707370132325240D-02  /
   data green(30,27, 4) /   2.4654709839535351D-02  /
   data green(30,27, 5) /   2.4587497137853378D-02  /
   data green(30,27, 6) /   2.4506091792462496D-02  /
   data green(30,27, 7) /   2.4410922475045865D-02  /
   data green(30,27, 8) /   2.4302480645442198D-02  /
   data green(30,27, 9) /   2.4181313747813755D-02  /
   data green(30,27,10) /   2.4048017898793423D-02  /
   data green(30,27,11) /   2.3903230248535828D-02  /
   data green(30,27,12) /   2.3747621192163518D-02  /
   data green(30,27,13) /   2.3581886599539706D-02  /
   data green(30,27,14) /   2.3406740216642827D-02  /
   data green(30,27,15) /   2.3222906373255730D-02  /
   data green(30,27,16) /   2.3031113110474323D-02  /
   data green(30,27,17) /   2.2832085818919108D-02  /
   data green(30,27,18) /   2.2626541455630349D-02  /
   data green(30,27,19) /   2.2415183385428937D-02  /
   data green(30,27,20) /   2.2198696871818618D-02  /
   data green(30,27,21) /   2.1977745223886888D-02  /
   data green(30,27,22) /   2.1752966589521874D-02  /
   data green(30,27,23) /   2.1524971371808265D-02  /
   data green(30,27,24) /   2.1294340234745957D-02  /
   data green(30,27,25) /   2.1061622656365856D-02  /
   data green(30,27,26) /   2.0827335981716155D-02  /
   data green(30,27,27) /   2.0591964924807526D-02  /
   data green(30,28, 0) /   2.4367628486876812D-02  /
   data green(30,28, 1) /   2.4360391774462593D-02  /
   data green(30,28, 2) /   2.4338720357271072D-02  /
   data green(30,28, 3) /   2.4302729819104561D-02  /
   data green(30,28, 4) /   2.4252610895315233D-02  /
   data green(30,28, 5) /   2.4188626674488456D-02  /
   data green(30,28, 6) /   2.4111108794241380D-02  /
   data green(30,28, 7) /   2.4020452738068101D-02  /
   data green(30,28, 8) /   2.3917112361937255D-02  /
   data green(30,28, 9) /   2.3801593795395450D-02  /
   data green(30,28,10) /   2.3674448871826886D-02  /
   data green(30,28,11) /   2.3536268246218111D-02  /
   data green(30,28,12) /   2.3387674356579187D-02  /
   data green(30,28,13) /   2.3229314377681939D-02  /
   data green(30,28,14) /   2.3061853303818328D-02  /
   data green(30,28,15) /   2.2885967281835277D-02  /
   data green(30,28,16) /   2.2702337297809988D-02  /
   data green(30,28,17) /   2.2511643301427006D-02  /
   data green(30,28,18) /   2.2314558832373852D-02  /
   data green(30,28,19) /   2.2111746193733829D-02  /
   data green(30,28,20) /   2.1903852199126182D-02  /
   data green(30,28,21) /   2.1691504503761228D-02  /
   data green(30,28,22) /   2.1475308515016136D-02  /
   data green(30,28,23) /   2.1255844865811304D-02  /
   data green(30,28,24) /   2.1033667424057833D-02  /
   data green(30,28,25) /   2.0809301803717671D-02  /
   data green(30,28,26) /   2.0583244337442200D-02  /
   data green(30,28,27) /   2.0355961467142324D-02  /
   data green(30,28,28) /   2.0127889506957267D-02  /
   data green(30,29, 0) /   2.3965426099240024D-02  /
   data green(30,29, 1) /   2.3958541930052060D-02  /
   data green(30,29, 2) /   2.3937925048665324D-02  /
   data green(30,29, 3) /   2.3903681820728360D-02  /
   data green(30,29, 4) /   2.3855987827213630D-02  /
   data green(30,29, 5) /   2.3795085371830844D-02  /
   data green(30,29, 6) /   2.3721280089158636D-02  /
   data green(30,29, 7) /   2.3634936745820671D-02  /
   data green(30,29, 8) /   2.3536474346079379D-02  /
   data green(30,29, 9) /   2.3426360667462127D-02  /
   data green(30,29,10) /   2.3305106361088887D-02  /
   data green(30,29,11) /   2.3173258755170310D-02  /
   data green(30,29,12) /   2.3031395498904969D-02  /
   data green(30,29,13) /   2.2880118178201406D-02  /
   data green(30,29,14) /   2.2720046024942604D-02  /
   data green(30,29,15) /   2.2551809828700257D-02  /
   data green(30,29,16) /   2.2376046144753643D-02  /
   data green(30,29,17) /   2.2193391875841992D-02  /
   data green(30,29,18) /   2.2004479288102213D-02  /
   data green(30,29,19) /   2.1809931504848466D-02  /
   data green(30,29,20) /   2.1610358505852464D-02  /
   data green(30,29,21) /   2.1406353645069463D-02  /
   data green(30,29,22) /   2.1198490686672993D-02  /
   data green(30,29,23) /   2.0987321348029941D-02  /
   data green(30,29,24) /   2.0773373328967611D-02  /
   data green(30,29,25) /   2.0557148799354472D-02  /
   data green(30,29,26) /   2.0339123311556843D-02  /
   data green(30,29,27) /   2.0119745100600300D-02  /
   data green(30,29,28) /   1.9899434732676330D-02  /
   data green(30,29,29) /   1.9678585061780763D-02  /
   data green(30,30, 0) /   2.3569409745247336D-02  /
   data green(30,30, 1) /   2.3562861307786519D-02  /
   data green(30,30, 2) /   2.3543248772067642D-02  /
   data green(30,30, 3) /   2.3510670011728248D-02  /
   data green(30,30, 4) /   2.3465286640798277D-02  /
   data green(30,30, 5) /   2.3407321793759529D-02  /
   data green(30,30, 6) /   2.3337057101848985D-02  /
   data green(30,30, 7) /   2.3254828945301540D-02  /
   data green(30,30, 8) /   2.3161024077878102D-02  /
   data green(30,30, 9) /   2.3056074732639067D-02  /
   data green(30,30,10) /   2.2940453326171593D-02  /
   data green(30,30,11) /   2.2814666882269904D-02  /
   data green(30,30,12) /   2.2679251295562572D-02  /
   data green(30,30,13) /   2.2534765551142373D-02  /
   data green(30,30,14) /   2.2381786008415224D-02  /
   data green(30,30,15) /   2.2220900846795809D-02  /
   data green(30,30,16) /   2.2052704758247291D-02  /
   data green(30,30,17) /   2.1877793957720808D-02  /
   data green(30,30,18) /   2.1696761567989099D-02  /
   data green(30,30,19) /   2.1510193420820476D-02  /
   data green(30,30,20) /   2.1318664302441271D-02  /
   data green(30,30,21) /   2.1122734658222807D-02  /
   data green(30,30,22) /   2.0922947759823333D-02  /
   data green(30,30,23) /   2.0719827327831147D-02  /
   data green(30,30,24) /   2.0513875594405078D-02  /
   data green(30,30,25) /   2.0305571783518540D-02  /
   data green(30,30,26) /   2.0095370981137464D-02  /
   data green(30,30,27) /   1.9883703363898952D-02  /
   data green(30,30,28) /   1.9670973752464932D-02  /
   data green(30,30,29) /   1.9457561454536421D-02  /
   data green(30,30,30) /   1.9243820362350835D-02  /
   data green(31, 0, 0) /   3.2266481175589140D-02  /
   data green(31, 1, 0) /   3.2249653635312424D-02  /
   data green(31, 1, 1) /   3.2232852533618420D-02  /
   data green(31, 2, 0) /   3.2199329646775902D-02  /
   data green(31, 2, 1) /   3.2182607442788003D-02  /
   data green(31, 2, 2) /   3.2132597800283214D-02  /
   data green(31, 3, 0) /   3.2115980926361483D-02  /
   data green(31, 3, 1) /   3.2099388840633070D-02  /
   data green(31, 3, 2) /   3.2049767508912932D-02  /
   data green(31, 3, 3) /   3.1967577664422313D-02  /
   data green(31, 4, 0) /   3.2000379950840395D-02  /
   data green(31, 4, 1) /   3.1983967190442324D-02  /
   data green(31, 4, 2) /   3.1934881035967075D-02  /
   data green(31, 4, 3) /   3.1853573927157708D-02  /
   data green(31, 4, 4) /   3.1740786997219363D-02  /
   data green(31, 5, 0) /   3.1853580124903642D-02  /
   data green(31, 5, 1) /   3.1837393181139713D-02  /
   data green(31, 5, 2) /   3.1788980982343680D-02  /
   data green(31, 5, 3) /   3.1708785612758694D-02  /
   data green(31, 5, 4) /   3.1597531361906100D-02  /
   data green(31, 5, 5) /   3.1456206608329013D-02  /
   data green(31, 6, 0) /   3.1676889436596557D-02  /
   data green(31, 6, 1) /   3.1660971477642484D-02  /
   data green(31, 6, 2) /   3.1613362115014504D-02  /
   data green(31, 6, 3) /   3.1534491226428374D-02  /
   data green(31, 6, 4) /   3.1425063241682261D-02  /
   data green(31, 6, 5) /   3.1286039720190027D-02  /
   data green(31, 6, 6) /   3.1118616160475870D-02  /
   data green(31, 7, 0) /   3.1471838876791003D-02  /
   data green(31, 7, 1) /   3.1456229251647365D-02  /
   data green(31, 7, 2) /   3.1409540228701328D-02  /
   data green(31, 7, 3) /   3.1332187866094124D-02  /
   data green(31, 7, 4) /   3.1224854102548560D-02  /
   data green(31, 7, 5) /   3.1088470110881353D-02  /
   data green(31, 7, 6) /   3.0924194124968545D-02  /
   data green(31, 7, 7) /   3.0733384762868861D-02  /
   data green(31, 8, 0) /   3.1240147077388827D-02  /
   data green(31, 8, 1) /   3.1224880941014506D-02  /
   data green(31, 8, 2) /   3.1179217265578010D-02  /
   data green(31, 8, 3) /   3.1103556932071806D-02  /
   data green(31, 8, 4) /   3.0998557170059820D-02  /
   data green(31, 8, 5) /   3.0865115751577955D-02  /
   data green(31, 8, 6) /   3.0704349921127898D-02  /
   data green(31, 8, 7) /   3.0517571019640505D-02  /
   data green(31, 8, 8) /   3.0306255903607644D-02  /
   data green(31, 9, 0) /   3.0983682694230318D-02  /
   data green(31, 9, 1) /   3.0968790757160487D-02  /
   data green(31, 9, 2) /   3.0924244192930297D-02  /
   data green(31, 9, 3) /   3.0850427611378763D-02  /
   data green(31, 9, 4) /   3.0747971741560676D-02  /
   data green(31, 9, 5) /   3.0617738512590167D-02  /
   data green(31, 9, 6) /   3.0460801149051871D-02  /
   data green(31, 9, 7) /   3.0278420171417072D-02  /
   data green(31, 9, 8) /   3.0072016327219300D-02  /
   data green(31, 9, 9) /   2.9843141543667601D-02  /
   data green(31,10, 0) /   3.0704426029643364D-02  /
   data green(31,10, 1) /   3.0689934429314100D-02  /
   data green(31,10, 2) /   3.0646583109559859D-02  /
   data green(31,10, 3) /   3.0574739577962654D-02  /
   data green(31,10, 4) /   3.0475006697785235D-02  /
   data green(31,10, 5) /   3.0348208685061911D-02  /
   data green(31,10, 6) /   3.0195372409368128D-02  /
   data green(31,10, 7) /   3.0017704820141388D-02  /
   data green(31,10, 8) /   2.9816567447389430D-02  /
   data green(31,10, 9) /   2.9593448988511899D-02  /
   data green(31,10,10) /   2.9349936992925374D-02  /
   data green(31,11, 0) /   3.0404431270521375D-02  /
   data green(31,11, 1) /   3.0390361556325157D-02  /
   data green(31,11, 2) /   3.0348269936172678D-02  /
   data green(31,11, 3) /   3.0278506238796261D-02  /
   data green(31,11, 4) /   3.0181644509308415D-02  /
   data green(31,11, 5) /   3.0058469935755651D-02  /
   data green(31,11, 6) /   2.9909961374184110D-02  /
   data green(31,11, 7) /   2.9737270225733008D-02  /
   data green(31,11, 8) /   2.9541696537682283D-02  /
   data green(31,11, 9) /   2.9324663260929688D-02  /
   data green(31,11,10) /   2.9087689599802118D-02  /
   data green(31,11,11) /   2.8832364342029029D-02  /
   data green(31,12, 0) /   3.0085790531429620D-02  /
   data green(31,12, 1) /   3.0072159751477989D-02  /
   data green(31,12, 2) /   3.0031378863150351D-02  /
   data green(31,12, 3) /   2.9963779679459839D-02  /
   data green(31,12, 4) /   2.9869906864423038D-02  /
   data green(31,12, 5) /   2.9750505787787943D-02  /
   data green(31,12, 6) /   2.9606506274326250D-02  /
   data green(31,12, 7) /   2.9439002934209646D-02  /
   data green(31,12, 8) /   2.9249232870784583D-02  /
   data green(31,12, 9) /   2.9038551619991505D-02  /
   data green(31,12,10) /   2.8808408182134541D-02  /
   data green(31,12,11) /   2.8560319966412839D-02  /
   data green(31,12,12) /   2.8295848389613999D-02  /
   data green(31,13, 0) /   2.9750600662189163D-02  /
   data green(31,13, 1) /   2.9737421536973310D-02  /
   data green(31,13, 2) /   2.9697989504063153D-02  /
   data green(31,13, 3) /   2.9632618243886411D-02  /
   data green(31,13, 4) /   2.9541822834292100D-02  /
   data green(31,13, 5) /   2.9426308519522704D-02  /
   data green(31,13, 6) /   2.9286955667026944D-02  /
   data green(31,13, 7) /   2.9124801533946830D-02  /
   data green(31,13, 8) /   2.8941019566354684D-02  /
   data green(31,13, 9) /   2.8736897009395981D-02  /
   data green(31,13,10) /   2.8513811615507492D-02  /
   data green(31,13,11) /   2.8273208204733766D-02  /
   data green(31,13,12) /   2.8016575762746306D-02  /
   data green(31,13,13) /   2.7745425667259606D-02  /
   data green(31,14, 0) /   2.9400933529233994D-02  /
   data green(31,14, 1) /   2.9388214695839431D-02  /
   data green(31,14, 2) /   2.9350157458084856D-02  /
   data green(31,14, 3) /   2.9287057443784379D-02  /
   data green(31,14, 4) /   2.9199400260240149D-02  /
   data green(31,14, 5) /   2.9087851152897779D-02  /
   data green(31,14, 6) /   2.8953241138911175D-02  /
   data green(31,14, 7) /   2.8796550175817667D-02  /
   data green(31,14, 8) /   2.8618888018388862D-02  /
   data green(31,14, 9) /   2.8421473468767302D-02  /
   data green(31,14,10) /   2.8205612736073280D-02  /
   data green(31,14,11) /   2.7972677594967514D-02  /
   data green(31,14,12) /   2.7724083974016860D-02  /
   data green(31,14,13) /   2.7461271521704508D-02  /
   data green(31,14,14) /   2.7185684598905170D-02  /
   data green(31,15, 0) /   2.9038810230758411D-02  /
   data green(31,15, 1) /   2.9026556540884103D-02  /
   data green(31,15, 2) /   2.8989888739676777D-02  /
   data green(31,15, 3) /   2.8929084653943383D-02  /
   data green(31,15, 4) /   2.8844600816187862D-02  /
   data green(31,15, 5) /   2.8737062979800554D-02  /
   data green(31,15, 6) /   2.8607253387172411D-02  /
   data green(31,15, 7) /   2.8456095292901374D-02  /
   data green(31,15, 8) /   2.8284635329015820D-02  /
   data green(31,15, 9) /   2.8094024348015525D-02  /
   data green(31,15,10) /   2.7885497392187226D-02  /
   data green(31,15,11) /   2.7660353416661492D-02  /
   data green(31,15,12) /   2.7419935343939120D-02  /
   data green(31,15,13) /   2.7165610955600702D-02  /
   data green(31,15,14) /   2.6898755039858287D-02  /
   data green(31,15,15) /   2.6620733118746895D-02  /
   data green(31,16, 0) /   2.8666179474043433D-02  /
   data green(31,16, 1) /   2.8654392329613514D-02  /
   data green(31,16, 2) /   2.8619118306028184D-02  /
   data green(31,16, 3) /   2.8560617826338543D-02  /
   data green(31,16, 4) /   2.8479318982203218D-02  /
   data green(31,16, 5) /   2.8375808865850144D-02  /
   data green(31,16, 6) /   2.8250821921662034D-02  /
   data green(31,16, 7) /   2.8105225765480989D-02  /
   data green(31,16, 8) /   2.7940004996726903D-02  /
   data green(31,16, 9) /   2.7756243574052073D-02  /
   data green(31,16,10) /   2.7555106339041440D-02  /
   data green(31,16,11) /   2.7337820256446758D-02  /
   data green(31,16,12) /   2.7105655897690348D-02  /
   data green(31,16,13) /   2.6859909632380084D-02  /
   data green(31,16,14) /   2.6601886916561616D-02  /
   data green(31,16,15) /   2.6332886982696611D-02  /
   data green(31,16,16) /   2.6054189150731225D-02  /
   data green(31,17, 0) /   2.8284900141224905D-02  /
   data green(31,17, 1) /   2.8273577852637982D-02  /
   data green(31,17, 2) /   2.8239692713475300D-02  /
   data green(31,17, 3) /   2.8183488260208580D-02  /
   data green(31,17, 4) /   2.8105364974232654D-02  /
   data green(31,17, 5) /   2.8005872386067256D-02  /
   data green(31,17, 6) /   2.7885698452853214D-02  /
   data green(31,17, 7) /   2.7745656607333975D-02  /
   data green(31,17, 8) /   2.7586670946156727D-02  /
   data green(31,17, 9) /   2.7409760067684648D-02  /
   data green(31,17,10) /   2.7216020084017686D-02  /
   data green(31,17,11) /   2.7006607320149943D-02  /
   data green(31,17,12) /   2.6782721178516194D-02  /
   data green(31,17,13) /   2.6545587594210884D-02  /
   data green(31,17,14) /   2.6296443440219638D-02  /
   data green(31,17,15) /   2.6036522168527615D-02  /
   data green(31,17,16) /   2.5767040897081036D-02  /
   data green(31,17,17) /   2.5489189078700295D-02  /
   data green(31,18, 0) /   2.7896727904057950D-02  /
   data green(31,18, 1) /   2.7885866057868375D-02  /
   data green(31,18, 2) /   2.7853356770376099D-02  /
   data green(31,18, 3) /   2.7799427303961413D-02  /
   data green(31,18, 4) /   2.7724451517100675D-02  /
   data green(31,18, 5) /   2.7628942693158142D-02  /
   data green(31,18, 6) /   2.7513543881915965D-02  /
   data green(31,18, 7) /   2.7379016106359219D-02  /
   data green(31,18, 8) /   2.7226224849970487D-02  /
   data green(31,18, 9) /   2.7056125278899102D-02  /
   data green(31,18,10) /   2.6869746668231092D-02  /
   data green(31,18,11) /   2.6668176493397817D-02  /
   data green(31,18,12) /   2.6452544619275715D-02  /
   data green(31,18,13) /   2.6224007974606840D-02  /
   data green(31,18,14) /   2.5983736042516980D-02  /
   data green(31,18,15) /   2.5732897433805749D-02  /
   data green(31,18,16) /   2.5472647742789463D-02  /
   data green(31,18,17) /   2.5204118819703490D-02  /
   data green(31,18,18) /   2.4928409532196702D-02  /
   data green(31,19, 0) /   2.7503305621508725D-02  /
   data green(31,19, 1) /   2.7492897446362092D-02  /
   data green(31,19, 2) /   2.7461743928310189D-02  /
   data green(31,19, 3) /   2.7410056740553565D-02  /
   data green(31,19, 4) /   2.7338184225615656D-02  /
   data green(31,19, 5) /   2.7246604899400017D-02  /
   data green(31,19, 6) /   2.7135918692454311D-02  /
   data green(31,19, 7) /   2.7006836239472160D-02  /
   data green(31,19, 8) /   2.6860166584369426D-02  /
   data green(31,19, 9) /   2.6696803704204002D-02  /
   data green(31,19,10) /   2.6517712270113694D-02  /
   data green(31,19,11) /   2.6323913058221418D-02  /
   data green(31,19,12) /   2.6116468400328365D-02  /
   data green(31,19,13) /   2.5896468026389948D-02  /
   data green(31,19,14) /   2.5665015602045853D-02  /
   data green(31,19,15) /   2.5423216208853860D-02  /
   data green(31,19,16) /   2.5172164956215447D-02  /
   data green(31,19,17) /   2.4912936855719440D-02  /
   data green(31,19,18) /   2.4646578033610356D-02  /
   data green(31,19,19) /   2.4374098307457442D-02  /
   data green(31,20, 0) /   2.7106157165320145D-02  /
   data green(31,20, 1) /   2.7096193887093278D-02  /
   data green(31,20, 2) /   2.7066370065192202D-02  /
   data green(31,20, 3) /   2.7016882520226607D-02  /
   data green(31,20, 4) /   2.6948055271618556D-02  /
   data green(31,20, 5) /   2.6860333667142269D-02  /
   data green(31,20, 6) /   2.6754276458808474D-02  /
   data green(31,20, 7) /   2.6630546098697829D-02  /
   data green(31,20, 8) /   2.6489897578706227D-02  /
   data green(31,20, 9) /   2.6333166171023614D-02  /
   data green(31,20,10) /   2.6161254440862119D-02  /
   data green(31,20,11) /   2.5975118900130368D-02  /
   data green(31,20,12) /   2.5775756652207469D-02  /
   data green(31,20,13) /   2.5564192346346229D-02  /
   data green(31,20,14) /   2.5341465718716670D-02  /
   data green(31,20,15) /   2.5108619949086169D-02  /
   data green(31,20,16) /   2.4866691010941498D-02  /
   data green(31,20,17) /   2.4616698141497195D-02  /
   data green(31,20,18) /   2.4359635509030936D-02  /
   data green(31,20,19) /   2.4096465110292432D-02  /
   data green(31,20,20) /   2.3828110891704016D-02  /
   data green(31,21, 0) /   2.6706684264622166D-02  /
   data green(31,21, 1) /   2.6697155443800139D-02  /
   data green(31,21, 2) /   2.6668630259640138D-02  /
   data green(31,21, 3) /   2.6621291450059481D-02  /
   data green(31,21, 4) /   2.6555439960190861D-02  /
   data green(31,21, 5) /   2.6471489648156446D-02  /
   data green(31,21, 6) /   2.6369960130969754D-02  /
   data green(31,21, 7) /   2.6251468010605693D-02  /
   data green(31,21, 8) /   2.6116716765181118D-02  /
   data green(31,21, 9) /   2.5966485620103265D-02  /
   data green(31,21,10) /   2.5801617728324632D-02  /
   data green(31,21,11) /   2.5623007987937929D-02  /
   data green(31,21,12) /   2.5431590810686026D-02  /
   data green(31,21,13) /   2.5228328128709857D-02  /
   data green(31,21,14) /   2.5014197891672630D-02  /
   data green(31,21,15) /   2.4790183265152171D-02  /
   data green(31,21,16) /   2.4557262696728019D-02  /
   data green(31,21,17) /   2.4316400971117855D-02  /
   data green(31,21,18) /   2.4068541332271088D-02  /
   data green(31,21,19) /   2.3814598710270859D-02  /
   data green(31,21,20) /   2.3555454055492350D-02  /
   data green(31,21,21) /   2.3291949752496342D-02  /
   data green(31,22, 0) /   2.6306165936208149D-02  /
   data green(31,22, 1) /   2.6297059782466756D-02  /
   data green(31,22, 2) /   2.6269798130911470D-02  /
   data green(31,22, 3) /   2.6224550424077966D-02  /
   data green(31,22, 4) /   2.6161595811623976D-02  /
   data green(31,22, 5) /   2.6081318384409899D-02  /
   data green(31,22, 6) /   2.5984200727377289D-02  /
   data green(31,22, 7) /   2.5870816001328484D-02  /
   data green(31,22, 8) /   2.5741818803634076D-02  /
   data green(31,22, 9) /   2.5597935085022065D-02  /
   data green(31,22,10) /   2.5439951413311863D-02  /
   data green(31,22,11) /   2.5268703875544727D-02  /
   data green(31,22,12) /   2.5085066898555736D-02  /
   data green(31,22,13) /   2.4889942246399118D-02  /
   data green(31,22,14) /   2.4684248423377995D-02  /
   data green(31,22,15) /   2.4468910676161852D-02  /
   data green(31,22,16) /   2.4244851750019846D-02  /
   data green(31,22,17) /   2.4012983514813766D-02  /
   data green(31,22,18) /   2.3774199538041756D-02  /
   data green(31,22,19) /   2.3529368646483643D-02  /
   data green(31,22,20) /   2.3279329486040397D-02  /
   data green(31,22,21) /   2.3024886061939198D-02  /
   data green(31,22,22) /   2.2766804218972852D-02  /
   data green(31,23, 0) /   2.5905760066586590D-02  /
   data green(31,23, 1) /   2.5897063727268844D-02  /
   data green(31,23, 2) /   2.5871027317359213D-02  /
   data green(31,23, 3) /   2.5827807775280185D-02  /
   data green(31,23, 4) /   2.5767663742012229D-02  /
   data green(31,23, 5) /   2.5690951277466810D-02  /
   data green(31,23, 6) /   2.5598118059667076D-02  /
   data green(31,23, 7) /   2.5489696250265873D-02  /
   data green(31,23, 8) /   2.5366294245334525D-02  /
   data green(31,23, 9) /   2.5228587554867562D-02  /
   data green(31,23,10) /   2.5077309067471028D-02  /
   data green(31,23,11) /   2.4913238958422335D-02  /
   data green(31,23,12) /   2.4737194490582603D-02  /
   data green(31,23,13) /   2.4550019939944594D-02  /
   data green(31,23,14) /   2.4352576852725683D-02  /
   data green(31,23,15) /   2.4145734810891190D-02  /
   data green(31,23,16) /   2.3930362849872350D-02  /
   data green(31,23,17) /   2.3707321637967963D-02  /
   data green(31,23,18) /   2.3477456493199478D-02  /
   data green(31,23,19) /   2.3241591281636567D-02  /
   data green(31,23,20) /   2.3000523212492692D-02  /
   data green(31,23,21) /   2.2755018520329944D-02  /
   data green(31,23,22) /   2.2505809003911175D-02  /
   data green(31,23,23) /   2.2253589374711787D-02  /
   data green(31,24, 0) /   2.5506506729586466D-02  /
   data green(31,24, 1) /   2.5498206550231762D-02  /
   data green(31,24, 2) /   2.5473354683036686D-02  /
   data green(31,24, 3) /   2.5432096346412806D-02  /
   data green(31,24, 4) /   2.5374670949192175D-02  /
   data green(31,24, 5) /   2.5301408245617284D-02  /
   data green(31,24, 6) /   2.5212723123102292D-02  /
   data green(31,24, 7) /   2.5109109182802663D-02  /
   data green(31,24, 8) /   2.4991131304359430D-02  /
   data green(31,24, 9) /   2.4859417408255431D-02  /
   data green(31,24,10) /   2.4714649641483736D-02  /
   data green(31,24,11) /   2.4557555214774113D-02  /
   data green(31,24,12) /   2.4388897113136371D-02  /
   data green(31,24,13) /   2.4209464887110340D-02  /
   data green(31,24,14) /   2.4020065711358465D-02  /
   data green(31,24,15) /   2.3821515871783779D-02  /
   data green(31,24,16) /   2.3614632813940106D-02  /
   data green(31,24,17) /   2.3400227855783520D-02  /
   data green(31,24,18) /   2.3179099638283032D-02  /
   data green(31,24,19) /   2.2952028359310629D-02  /
   data green(31,24,20) /   2.2719770810536541D-02  /
   data green(31,24,21) /   2.2483056214444325D-02  /
   data green(31,24,22) /   2.2242582839462524D-02  /
   data green(31,24,23) /   2.1999015355742724D-02  /
   data green(31,24,24) /   2.1752982882261282D-02  /
   data green(31,25, 0) /   2.5109332853774431D-02  /
   data green(31,25, 1) /   2.5101414610068179D-02  /
   data green(31,25, 2) /   2.5077704871525778D-02  /
   data green(31,25, 3) /   2.5038337904518260D-02  /
   data green(31,25, 4) /   2.4983535137233427D-02  /
   data green(31,25, 5) /   2.4913601712229828D-02  /
   data green(31,25, 6) /   2.4828921808403286D-02  /
   data green(31,25, 7) /   2.4729952871713680D-02  /
   data green(31,25, 8) /   2.4617218921687684D-02  /
   data green(31,25, 9) /   2.4491303120527438D-02  /
   data green(31,25,10) /   2.4352839803111072D-02  /
   data green(31,25,11) /   2.4202506169296232D-02  /
   data green(31,25,12) /   2.4041013835250709D-02  /
   data green(31,25,13) /   2.3869100428965902D-02  /
   data green(31,25,14) /   2.3687521397876325D-02  /
   data green(31,25,15) /   2.3497042175018263D-02  /
   data green(31,25,16) /   2.3298430825868743D-02  /
   data green(31,25,17) /   2.3092451272331581D-02  /
   data green(31,25,18) /   2.2879857164566468D-02  /
   data green(31,25,19) /   2.2661386446586211D-02  /
   data green(31,25,20) /   2.2437756638651329D-02  /
   data green(31,25,21) /   2.2209660839101138D-02  /
   data green(31,25,22) /   2.1977764430781227D-02  /
   data green(31,25,23) /   2.1742702462849147D-02  /
   data green(31,25,24) /   2.1505077667477456D-02  /
   data green(31,25,25) /   2.1265459062690773D-02  /
   data green(31,26, 0) /   2.4715057892542042D-02  /
   data green(31,26, 1) /   2.4707506994008626D-02  /
   data green(31,26, 2) /   2.4684895863656253D-02  /
   data green(31,26, 3) /   2.4647348560636613D-02  /
   data green(31,26, 4) /   2.4595069747382986D-02  /
   data green(31,26, 5) /   2.4528341601460248D-02  /
   data green(31,26, 6) /   2.4447519620945525D-02  /
   data green(31,26, 7) /   2.4353027444526797D-02  /
   data green(31,26, 8) /   2.4245350831896582D-02  /
   data green(31,26, 9) /   2.4125030967759122D-02  /
   data green(31,26,10) /   2.3992657263400181D-02  /
   data green(31,26,11) /   2.3848859833267986D-02  /
   data green(31,26,12) /   2.3694301820776517D-02  /
   data green(31,26,13) /   2.3529671738308673D-02  /
   data green(31,26,14) /   2.3355675972161767D-02  /
   data green(31,26,15) /   2.3173031585100274D-02  /
   data green(31,26,16) /   2.2982459528488302D-02  /
   data green(31,26,17) /   2.2784678353864919D-02  /
   data green(31,26,18) /   2.2580398491411323D-02  /
   data green(31,26,19) /   2.2370317140998477D-02  /
   data green(31,26,20) /   2.2155113801177843D-02  /
   data green(31,26,21) /   2.1935446443168454D-02  /
   data green(31,26,22) /   2.1711948320992863D-02  /
   data green(31,26,23) /   2.1485225395628929D-02  /
   data green(31,26,24) /   2.1255854340429134D-02  /
   data green(31,26,25) /   2.1024381087037379D-02  /
   data green(31,26,26) /   2.0791319865431365D-02  /
   data green(31,27, 0) /   2.4324400192545632D-02  /
   data green(31,27, 1) /   2.4317201859030614D-02  /
   data green(31,27, 2) /   2.4295645237682111D-02  /
   data green(31,27, 3) /   2.4259844896809459D-02  /
   data green(31,27, 4) /   2.4209989902579925D-02  /
   data green(31,27, 5) /   2.4146341054747795D-02  /
   data green(31,27, 6) /   2.4069227128351346D-02  /
   data green(31,27, 7) /   2.3979040226668073D-02  /
   data green(31,27, 8) /   2.3876230372186114D-02  /
   data green(31,27, 9) /   2.3761299478198033D-02  /
   data green(31,27,10) /   2.3634794853427902D-02  /
   data green(31,27,11) /   2.3497302395816006D-02  /
   data green(31,27,12) /   2.3349439629499716D-02  /
   data green(31,27,13) /   2.3191848731730899D-02  /
   data green(31,27,14) /   2.3025189684769774D-02  /
   data green(31,27,15) /   2.2850133672646027D-02  /
   data green(31,27,16) /   2.2667356825106726D-02  /
   data green(31,27,17) /   2.2477534392091988D-02  /
   data green(31,27,18) /   2.2281335412644460D-02  /
   data green(31,27,19) /   2.2079417923105419D-02  /
   data green(31,27,20) /   2.1872424731469293D-02  /
   data green(31,27,21) /   2.1660979768393370D-02  /
   data green(31,27,22) /   2.1445685010959680D-02  /
   data green(31,27,23) /   2.1227117963079845D-02  /
   data green(31,27,24) /   2.1005829666502194D-02  /
   data green(31,27,25) /   2.0782343208692534D-02  /
   data green(31,27,26) /   2.0557152688294323D-02  /
   data green(31,27,27) /   2.0330722595244228D-02  /
   data green(31,28, 0) /   2.3937983800182808D-02  /
   data green(31,28, 1) /   2.3931123212685392D-02  /
   data green(31,28, 2) /   2.3910576873662506D-02  /
   data green(31,28, 3) /   2.3876450544708969D-02  /
   data green(31,28, 4) /   2.3828918813458728D-02  /
   data green(31,28, 5) /   2.3768222620617896D-02  /
   data green(31,28, 6) /   2.3694665894578942D-02  /
   data green(31,28, 7) /   2.3608611385024635D-02  /
   data green(31,28, 8) /   2.3510475805802653D-02  /
   data green(31,28, 9) /   2.3400724411477197D-02  /
   data green(31,28,10) /   2.3279865140961704D-02  /
   data green(31,28,11) /   2.3148442465437388D-02  /
   data green(31,28,12) /   2.3007031076578328D-02  /
   data green(31,28,13) /   2.2856229545401691D-02  /
   data green(31,28,14) /   2.2696654072492075D-02  /
   data green(31,28,15) /   2.2528932437701085D-02  /
   data green(31,28,16) /   2.2353698242548008D-02  /
   data green(31,28,17) /   2.2171585522302219D-02  /
   data green(31,28,18) /   2.1983223787926890D-02  /
   data green(31,28,19) /   2.1789233541429635D-02  /
   data green(31,28,20) /   2.1590222292312878D-02  /
   data green(31,28,21) /   2.1386781088224262D-02  /
   data green(31,28,22) /   2.1179481559923486D-02  /
   data green(31,28,23) /   2.0968873469524394D-02  /
   data green(31,28,24) /   2.0755482741741667D-02  /
   data green(31,28,25) /   2.0539809950571697D-02  /
   data green(31,28,26) /   2.0322329228388294D-02  /
   data green(31,28,27) /   2.0103487560697287D-02  /
   data green(31,28,28) /   1.9883704427588477D-02  /
   data green(31,29, 0) /   2.3556345488693734D-02  /
   data green(31,29, 1) /   2.3549807916452968D-02  /
   data green(31,29, 2) /   2.3530227886003226D-02  /
   data green(31,29, 3) /   2.3497703001565222D-02  /
   data green(31,29, 4) /   2.3452394433933318D-02  /
   data green(31,29, 5) /   2.3394524709011877D-02  /
   data green(31,29, 6) /   2.3324374695574055D-02  /
   data green(31,29, 7) /   2.3242279871554389D-02  /
   data green(31,29, 8) /   2.3148625964761568D-02  /
   data green(31,29, 9) /   2.3043844076462362D-02  /
   data green(31,29,10) /   2.2928405404508481D-02  /
   data green(31,29,11) /   2.2802815686468785D-02  /
   data green(31,29,12) /   2.2667609482742275D-02  /
   data green(31,29,13) /   2.2523344415230195D-02  /
   data green(31,29,14) /   2.2370595469361654D-02  /
   data green(31,29,15) /   2.2209949456743854D-02  /
   data green(31,29,16) /   2.2041999723152478D-02  /
   data green(31,29,17) /   2.1867341172710206D-02  /
   data green(31,29,18) /   2.1686565664615038D-02  /
   data green(31,29,19) /   2.1500257824300869D-02  /
   data green(31,29,20) /   2.1308991296977282D-02  /
   data green(31,29,21) /   2.1113325458536839D-02  /
   data green(31,29,22) /   2.0913802587156319D-02  /
   data green(31,29,23) /   2.0710945488768541D-02  /
   data green(31,29,24) /   2.0505255561054932D-02  /
   data green(31,29,25) /   2.0297211273734884D-02  /
   data green(31,29,26) /   2.0087267037659087D-02  /
   data green(31,29,27) /   1.9875852431451187D-02  /
   data green(31,29,28) /   1.9663371752044710D-02  /
   data green(31,29,29) /   1.9450203854263690D-02  /
   data green(31,30, 0) /   2.3179941828715477D-02  /
   data green(31,30, 1) /   2.3173712734650000D-02  /
   data green(31,30, 2) /   2.3155055607780223D-02  /
   data green(31,30, 3) /   2.3124060508023417D-02  /
   data green(31,30, 4) /   2.3080876192437230D-02  /
   data green(31,30, 5) /   2.3025708138124133D-02  /
   data green(31,30, 6) /   2.2958815846857737D-02  /
   data green(31,30, 7) /   2.2880509500220417D-02  /
   data green(31,30, 8) /   2.2791146048584572D-02  /
   data green(31,30, 9) /   2.2691124828430196D-02  /
   data green(31,30,10) /   2.2580882809971392D-02  /
   data green(31,30,11) /   2.2460889580772303D-02  /
   data green(31,30,12) /   2.2331642171076634D-02  /
   data green(31,30,13) /   2.2193659823236775D-02  /
   data green(31,30,14) /   2.2047478801333967D-02  /
   data green(31,30,15) /   2.1893647328355893D-02  /
   data green(31,30,16) /   2.1732720727729908D-02  /
   data green(31,30,17) /   2.1565256834202195D-02  /
   data green(31,30,18) /   2.1391811726594417D-02  /
   data green(31,30,19) /   2.1212935822397753D-02  /
   data green(31,30,20) /   2.1029170361949685D-02  /
   data green(31,30,21) /   2.0841044298466650D-02  /
   data green(31,30,22) /   2.0649071599770041D-02  /
   data green(31,30,23) /   2.0453748958349286D-02  /
   data green(31,30,24) /   2.0255553898572771D-02  /
   data green(31,30,25) /   2.0054943263429863D-02  /
   data green(31,30,26) /   1.9852352058146779D-02  /
   data green(31,30,27) /   1.9648192624295322D-02  /
   data green(31,30,28) /   1.9442854115500736D-02  /
   data green(31,30,29) /   1.9236702244420621D-02  /
   data green(31,30,30) /   1.9030079270165236D-02  /
   data green(31,31, 0) /   2.2809156161699377D-02  /
   data green(31,31, 1) /   2.2803221288383826D-02  /
   data green(31,31, 2) /   2.2785444486593500D-02  /
   data green(31,31, 3) /   2.2755908848124361D-02  /
   data green(31,31, 4) /   2.2714751659677655D-02  /
   data green(31,31, 5) /   2.2662162635523914D-02  /
   data green(31,31, 6) /   2.2598381506044506D-02  /
   data green(31,31, 7) /   2.2523695021801470D-02  /
   data green(31,31, 8) /   2.2438433445537746D-02  /
   data green(31,31, 9) /   2.2342966614408155D-02  /
   data green(31,31,10) /   2.2237699661524769D-02  /
   data green(31,31,11) /   2.2123068489473294D-02  /
   data green(31,31,12) /   2.1999535088893962D-02  /
   data green(31,31,13) /   2.1867582792737931D-02  /
   data green(31,31,14) /   2.1727711551750909D-02  /
   data green(31,31,15) /   2.1580433309527268D-02  /
   data green(31,31,16) /   2.1426267546606228D-02  /
   data green(31,31,17) /   2.1265737053051667D-02  /
   data green(31,31,18) /   2.1099363978267050D-02  /
   data green(31,31,19) /   2.0927666195907666D-02  /
   data green(31,31,20) /   2.0751154011074354D-02  /
   data green(31,31,21) /   2.0570327226845857D-02  /
   data green(31,31,22) /   2.0385672577900128D-02  /
   data green(31,31,23) /   2.0197661530680484D-02  /
   data green(31,31,24) /   2.0006748442401589D-02  /
   data green(31,31,25) /   1.9813369065220236D-02  /
   data green(31,31,26) /   1.9617939377116398D-02  /
   data green(31,31,27) /   1.9420854717397882D-02  /
   data green(31,31,28) /   1.9222489202177496D-02  /
   data green(31,31,29) /   1.9023195393571037D-02  /
   data green(31,31,30) /   1.8823304195610065D-02  /
   data green(31,31,31) /   1.8623124949829411D-02  /
   data green(32, 0, 0) /   3.1257650535402713D-02  /
   data green(32, 1, 0) /   3.1242354264524749D-02  /
   data green(32, 1, 1) /   3.1227080538619983D-02  /
   data green(32, 2, 0) /   3.1196600722326710D-02  /
   data green(32, 2, 1) /   3.1181394297065838D-02  /
   data green(32, 2, 2) /   3.1135908961610694D-02  /
   data green(32, 3, 0) /   3.1120792377303606D-02  /
   data green(32, 3, 1) /   3.1105697017187971D-02  /
   data green(32, 3, 2) /   3.1060543240446294D-02  /
   data green(32, 3, 3) /   3.0985724718121391D-02  /
   data green(32, 4, 0) /   3.1015589030790960D-02  /
   data green(32, 4, 1) /   3.1000646886493885D-02  /
   data green(32, 4, 2) /   3.0955950513010842D-02  /
   data green(32, 4, 3) /   3.0881886926086729D-02  /
   data green(32, 4, 4) /   3.0779090762238101D-02  /
   data green(32, 5, 0) /   3.0881891905690708D-02  /
   data green(32, 5, 1) /   3.0867142946067357D-02  /
   data green(32, 5, 2) /   3.0823023319631996D-02  /
   data green(32, 5, 3) /   3.0749911717091044D-02  /
   data green(32, 5, 4) /   3.0648429213552389D-02  /
   data green(32, 5, 5) /   3.0519424653564945D-02  /
   data green(32, 6, 0) /   3.0720822440336065D-02  /
   data green(32, 6, 1) /   3.0706303956200328D-02  /
   data green(32, 6, 2) /   3.0662872440570012D-02  /
   data green(32, 6, 3) /   3.0590896749287343D-02  /
   data green(32, 6, 4) /   3.0490981932830445D-02  /
   data green(32, 6, 5) /   3.0363955148139953D-02  /
   data green(32, 6, 6) /   3.0210846848607458D-02  /
   data green(32, 7, 0) /   3.0533696754909849D-02  /
   data green(32, 7, 1) /   3.0519442947470960D-02  /
   data green(32, 7, 2) /   3.0476801698892966D-02  /
   data green(32, 7, 3) /   3.0406130701611814D-02  /
   data green(32, 7, 4) /   3.0308016811831252D-02  /
   data green(32, 7, 5) /   3.0183262552497342D-02  /
   data green(32, 7, 6) /   3.0032868081523570D-02  /
   data green(32, 7, 7) /   2.9858009410765097D-02  /
   data green(32, 8, 0) /   3.0321996901924095D-02  /
   data green(32, 8, 1) /   3.0308038564081985D-02  /
   data green(32, 8, 2) /   3.0266279580062807D-02  /
   data green(32, 8, 3) /   3.0197065346143182D-02  /
   data green(32, 8, 4) /   3.0100962671880040D-02  /
   data green(32, 8, 5) /   2.9978746925709138D-02  /
   data green(32, 8, 6) /   2.9831384849919800D-02  /
   data green(32, 8, 7) /   2.9660013783467099D-02  /
   data green(32, 8, 8) /   2.9465918147391162D-02  /
   data green(32, 9, 0) /   3.0087340078669096D-02  /
   data green(32, 9, 1) /   3.0073704372580369D-02  /
   data green(32, 9, 2) /   3.0032908827137498D-02  /
   data green(32, 9, 3) /   2.9965285613124202D-02  /
   data green(32, 9, 4) /   2.9871379970963001D-02  /
   data green(32, 9, 5) /   2.9751938037360825D-02  /
   data green(32, 9, 6) /   2.9607890558534092D-02  /
   data green(32, 9, 7) /   2.9440333179303762D-02  /
   data green(32, 9, 8) /   2.9250504107346476D-02  /
   data green(32, 9, 9) /   2.9039760009714120D-02  /
   data green(32,10, 0) /   2.9831446970485561D-02  /
   data green(32,10, 1) /   2.9818157300297368D-02  /
   data green(32,10, 2) /   2.9778395162509720D-02  /
   data green(32,10, 3) /   2.9712478777619324D-02  /
   data green(32,10, 4) /   2.9620930626412200D-02  /
   data green(32,10, 5) /   2.9504465983336906D-02  /
   data green(32,10, 6) /   2.9363977562834050D-02  /
   data green(32,10, 7) /   2.9200516917654358D-02  /
   data green(32,10, 8) /   2.9015273331602785D-02  /
   data green(32,10, 9) /   2.8809551004830222D-02  /
   data green(32,10,10) /   2.8584745337956290D-02  /
   data green(32,11, 0) /   2.9556110318122969D-02  /
   data green(32,11, 1) /   2.9543186292821433D-02  /
   data green(32,11, 2) /   2.9504516213746876D-02  /
   data green(32,11, 3) /   2.9440403826132021D-02  /
   data green(32,11, 4) /   2.9351347985256546D-02  /
   data green(32,11, 5) /   2.9238031910204215D-02  /
   data green(32,11, 6) /   2.9101308782396244D-02  /
   data green(32,11, 7) /   2.8942184277544591D-02  /
   data green(32,11, 8) /   2.8761796716316197D-02  /
   data green(32,11, 9) /   2.8561395572329059D-02  /
   data green(32,11,10) /   2.8342319086095344D-02  /
   data green(32,11,11) /   2.8105971703727126D-02  /
   data green(32,12, 0) /   2.9263164674021066D-02  /
   data green(32,12, 1) /   2.9250622152159827D-02  /
   data green(32,12, 2) /   2.9213091595717354D-02  /
   data green(32,12, 3) /   2.9150861939820891D-02  /
   data green(32,12, 4) /   2.9064407857610206D-02  /
   data green(32,12, 5) /   2.8954379737043769D-02  /
   data green(32,12, 6) /   2.8821590236291023D-02  /
   data green(32,12, 7) /   2.8666997956671065D-02  /
   data green(32,12, 8) /   2.8491688861931166D-02  /
   data green(32,12, 9) /   2.8296856123449017D-02  /
   data green(32,12,10) /   2.8083779082494060D-02  /
   data green(32,12,11) /   2.7853801995971109D-02  /
   data green(32,12,12) /   2.7608313176608213D-02  /
   data green(32,13, 0) /   2.8954458147408477D-02  /
   data green(32,13, 1) /   2.8942309353176082D-02  /
   data green(32,13, 2) /   2.8905954939712569D-02  /
   data green(32,13, 3) /   2.8845668873702754D-02  /
   data green(32,13, 4) /   2.8761901376924540D-02  /
   data green(32,13, 5) /   2.8655269619037572D-02  /
   data green(32,13, 6) /   2.8526545221650313D-02  /
   data green(32,13, 7) /   2.8376639064355783D-02  /
   data green(32,13, 8) /   2.8206583966482696D-02  /
   data green(32,13, 9) /   2.8017515866412612D-02  /
   data green(32,13,10) /   2.7810654133121601D-02  /
   data green(32,13,11) /   2.7587281624578821D-02  /
   data green(32,13,12) /   2.7348725059517277D-02  /
   data green(32,13,13) /   2.7096336199157540D-02  /
   data green(32,14, 0) /   2.8631826753442582D-02  /
   data green(32,14, 1) /   2.8620080452112080D-02  /
   data green(32,14, 2) /   2.8584928479139855D-02  /
   data green(32,14, 3) /   2.8526629834405227D-02  /
   data green(32,14, 4) /   2.8445610279947534D-02  /
   data green(32,14, 5) /   2.8342453733330018D-02  /
   data green(32,14, 6) /   2.8217890701146435D-02  /
   data green(32,14, 7) /   2.8072784197073458D-02  /
   data green(32,14, 8) /   2.7908113665272929D-02  /
   data green(32,14, 9) /   2.7724957475260305D-02  /
   data green(32,14,10) /   2.7524474568109723D-02  /
   data green(32,14,11) /   2.7307885818074789D-02  /
   data green(32,14,12) /   2.7076455632401721D-02  /
   data green(32,14,13) /   2.6831474250737827D-02  /
   data green(32,14,14) /   2.6574241130246531D-02  /
   data green(32,15, 0) /   2.8297071792717291D-02  /
   data green(32,15, 1) /   2.8285733512936603D-02  /
   data green(32,15, 2) /   2.8251800615763288D-02  /
   data green(32,15, 3) /   2.8195517277470817D-02  /
   data green(32,15, 4) /   2.8117285023160268D-02  /
   data green(32,15, 5) /   2.8017654798393504D-02  /
   data green(32,15, 6) /   2.7897316303717400D-02  /
   data green(32,15, 7) /   2.7757084992607928D-02  /
   data green(32,15, 8) /   2.7597887203308930D-02  /
   data green(32,15, 9) /   2.7420743937491450D-02  /
   data green(32,15,10) /   2.7226753813061329D-02  /
   data green(32,15,11) /   2.7017075706404495D-02  /
   data green(32,15,12) /   2.6792911564290502D-02  /
   data green(32,15,13) /   2.6555489812221052D-02  /
   data green(32,15,14) /   2.6306049719574526D-02  /
   data green(32,15,15) /   2.6045827007947196D-02  /
   data green(32,16, 0) /   2.7951940507284962D-02  /
   data green(32,16, 1) /   2.7941012797800646D-02  /
   data green(32,16, 2) /   2.7908306713145134D-02  /
   data green(32,16, 3) /   2.7854051871388719D-02  /
   data green(32,16, 4) /   2.7778625983406240D-02  /
   data green(32,16, 5) /   2.7682547575000119D-02  /
   data green(32,16, 6) /   2.7566466186681338D-02  /
   data green(32,16, 7) /   2.7431150410482767D-02  /
   data green(32,16, 8) /   2.7277474186928493D-02  /
   data green(32,16, 9) /   2.7106401824846756D-02  /
   data green(32,16,10) /   2.6918972221478554D-02  /
   data green(32,16,11) /   2.6716282751568444D-02  /
   data green(32,16,12) /   2.6499473264671877D-02  /
   data green(32,16,13) /   2.6269710583757488D-02  /
   data green(32,16,14) /   2.6028173839952921D-02  /
   data green(32,16,15) /   2.5776040912776611D-02  /
   data green(32,16,16) /   2.5514476176971303D-02  /
   data green(32,17, 0) /   2.7598110097583093D-02  /
   data green(32,17, 1) /   2.7587592806697904D-02  /
   data green(32,17, 2) /   2.7556113204466284D-02  /
   data green(32,17, 3) /   2.7503886718861691D-02  /
   data green(32,17, 4) /   2.7431267837604290D-02  /
   data green(32,17, 5) /   2.7338743448872247D-02  /
   data green(32,17, 6) /   2.7226923865002948D-02  /
   data green(32,17, 7) /   2.7096531850314361D-02  /
   data green(32,17, 8) /   2.6948390032035651D-02  /
   data green(32,17, 9) /   2.6783407110023735D-02  /
   data green(32,17,10) /   2.6602563295815000D-02  /
   data green(32,17,11) /   2.6406895405593878D-02  /
   data green(32,17,12) /   2.6197482007208481D-02  /
   data green(32,17,13) /   2.5975428981797135D-02  /
   data green(32,17,14) /   2.5741855809881680D-02  /
   data green(32,17,15) /   2.5497882834099910D-02  /
   data green(32,17,16) /   2.5244619690093933D-02  /
   data green(32,17,17) /   2.4983155036997093D-02  /
   data green(32,18, 0) /   2.7237175047623314D-02  /
   data green(32,18, 1) /   2.7227065614812945D-02  /
   data green(32,18, 2) /   2.7196804966535395D-02  /
   data green(32,18, 3) /   2.7146594792502095D-02  /
   data green(32,18, 4) /   2.7076767085949239D-02  /
   data green(32,18, 5) /   2.6987778068153705D-02  /
   data green(32,18, 6) /   2.6880199990765463D-02  /
   data green(32,18, 7) /   2.6754711101846192D-02  /
   data green(32,18, 8) /   2.6612084113809324D-02  /
   data green(32,18, 9) /   2.6453173545317392D-02  /
   data green(32,18,10) /   2.6278902323946477D-02  /
   data green(32,18,11) /   2.6090248032811739D-02  /
   data green(32,18,12) /   2.5888229164292120D-02  /
   data green(32,18,13) /   2.5673891710329310D-02  /
   data green(32,18,14) /   2.5448296374891360D-02  /
   data green(32,18,15) /   2.5212506643683380D-02  /
   data green(32,18,16) /   2.4967577892549310D-02  /
   data green(32,18,17) /   2.4714547662392044D-02  /
   data green(32,18,18) /   2.4454427177469661D-02  /
   data green(32,19, 0) /   2.6870637596606010D-02  /
   data green(32,19, 1) /   2.6860931347112382D-02  /
   data green(32,19, 2) /   2.6831875803666166D-02  /
   data green(32,19, 3) /   2.6783659435342496D-02  /
   data green(32,19, 4) /   2.6716592578034406D-02  /
   data green(32,19, 5) /   2.6631101906208526D-02  /
   data green(32,19, 6) /   2.6527722966078200D-02  /
   data green(32,19, 7) /   2.6407091023853621D-02  /
   data green(32,19, 8) /   2.6269930529843895D-02  /
   data green(32,19, 9) /   2.6117043530296188D-02  /
   data green(32,19,10) /   2.5949297373303117D-02  /
   data green(32,19,11) /   2.5767612053430686D-02  /
   data green(32,19,12) /   2.5572947523483481D-02  /
   data green(32,19,13) /   2.5366291273394913D-02  /
   data green(32,19,14) /   2.5148646438478840D-02  /
   data green(32,19,15) /   2.4921020655288725D-02  /
   data green(32,19,16) /   2.4684415836140651D-02  /
   data green(32,19,17) /   2.4439818985736026D-02  /
   data green(32,19,18) /   2.4188194137621230D-02  /
   data green(32,19,19) /   2.3930475446285165D-02  /
   data green(32,20, 0) /   2.6499901114210894D-02  /
   data green(32,20, 1) /   2.6490591548953311D-02  /
   data green(32,20, 2) /   2.6462721804718502D-02  /
   data green(32,20, 3) /   2.6416467696831258D-02  /
   data green(32,20, 4) /   2.6352118822560111D-02  /
   data green(32,20, 5) /   2.6270073542717629D-02  /
   data green(32,20, 6) /   2.6170832196639099D-02  /
   data green(32,20, 7) /   2.6054988774921992D-02  /
   data green(32,20, 8) /   2.5923221316608059D-02  /
   data green(32,20, 9) /   2.5776281325951313D-02  /
   data green(32,20,10) /   2.5614982517899981D-02  /
   data green(32,20,11) /   2.5440189201303195D-02  /
   data green(32,20,12) /   2.5252804595900406D-02  /
   data green(32,20,13) /   2.5053759355320969D-02  /
   data green(32,20,14) /   2.4844000536031627D-02  /
   data green(32,20,15) /   2.4624481214053938D-02  /
   data green(32,20,16) /   2.4396150909954655D-02  /
   data green(32,20,17) /   2.4159946940516819D-02  /
   data green(32,20,18) /   2.3916786774715686D-02  /
   data green(32,20,19) /   2.3667561433822423D-02  /
   data green(32,20,20) /   2.3413129941853524D-02  /
   data green(32,21, 0) /   2.6126266082573249D-02  /
   data green(32,21, 1) /   2.6117345157261461D-02  /
   data green(32,21, 2) /   2.6090637282442378D-02  /
   data green(32,21, 3) /   2.6046306220949896D-02  /
   data green(32,21, 4) /   2.5984621807462936D-02  /
   data green(32,21, 5) /   2.5905955402521323D-02  /
   data green(32,21, 6) /   2.5810773740089769D-02  /
   data green(32,21, 7) /   2.5699631366619128D-02  /
   data green(32,21, 8) /   2.5573161907414359D-02  /
   data green(32,21, 9) /   2.5432068422053332D-02  /
   data green(32,21,10) /   2.5277113124004373D-02  /
   data green(32,21,11) /   2.5109106740701874D-02  /
   data green(32,21,12) /   2.4928897780177640D-02  /
   data green(32,21,13) /   2.4737361950517959D-02  /
   data green(32,21,14) /   2.4535391950952785D-02  /
   data green(32,21,15) /   2.4323887820517923D-02  /
   data green(32,21,16) /   2.4103747994218036D-02  /
   data green(32,21,17) /   2.3875861179571609D-02  /
   data green(32,21,18) /   2.3641099130182360D-02  /
   data green(32,21,19) /   2.3400310359061500D-02  /
   data green(32,21,20) /   2.3154314803942359D-02  /
   data green(32,21,21) /   2.2903899430526706D-02  /
   data green(32,22, 0) /   2.5750928357388722D-02  /
   data green(32,22, 1) /   2.5742386746201187D-02  /
   data green(32,22, 2) /   2.5716812973443109D-02  /
   data green(32,22, 3) /   2.5674359371971416D-02  /
   data green(32,22, 4) /   2.5615277025806515D-02  /
   data green(32,22, 5) /   2.5539911659746613D-02  /
   data green(32,22, 6) /   2.5448698070977868D-02  /
   data green(32,22, 7) /   2.5342153276892328D-02  /
   data green(32,22, 8) /   2.5220868587119496D-02  /
   data green(32,22, 9) /   2.5085500831324877D-02  /
   data green(32,22,10) /   2.4936762987058922D-02  /
   data green(32,22,11) /   2.4775414453994751D-02  /
   data green(32,22,12) /   2.4602251213075889D-02  /
   data green(32,22,13) /   2.4418096092724609D-02  /
   data green(32,22,14) /   2.4223789341029464D-02  /
   data green(32,22,15) /   2.4020179674623384D-02  /
   data green(32,22,16) /   2.3808115943710958D-02  /
   data green(32,22,17) /   2.3588439520236117D-02  /
   data green(32,22,18) /   2.3361977484119504D-02  /
   data green(32,22,19) /   2.3129536652187983D-02  /
   data green(32,22,20) /   2.2891898466891630D-02  /
   data green(32,22,21) /   2.2649814737868385D-02  /
   data green(32,22,22) /   2.2404004209272591D-02  /
   data green(32,23, 0) /   2.5374979369910976D-02  /
   data green(32,23, 1) /   2.5366806710453587D-02  /
   data green(32,23, 2) /   2.5342336165910615D-02  /
   data green(32,23, 3) /   2.5301709271610458D-02  /
   data green(32,23, 4) /   2.5245159390222170D-02  /
   data green(32,23, 5) /   2.5173008001271011D-02  /
   data green(32,23, 6) /   2.5085659669619936D-02  /
   data green(32,23, 7) /   2.4983595845917490D-02  /
   data green(32,23, 8) /   2.4867367682100020D-02  /
   data green(32,23, 9) /   2.4737588066345989D-02  /
   data green(32,23,10) /   2.4594923093874593D-02  /
   data green(32,23,11) /   2.4440083192732477D-02  /
   data green(32,23,12) /   2.4273814117845856D-02  /
   data green(32,23,13) /   2.4096888013205704D-02  /
   data green(32,23,14) /   2.3910094722505364D-02  /
   data green(32,23,15) /   2.3714233504448432D-02  /
   data green(32,23,16) /   2.3510105281933476D-02  /
   data green(32,23,17) /   2.3298505525980881D-02  /
   data green(32,23,18) /   2.3080217847010958D-02  /
   data green(32,23,19) /   2.2856008339112256D-02  /
   data green(32,23,20) /   2.2626620698190006D-02  /
   data green(32,23,21) /   2.2392772113020661D-02  /
   data green(32,23,22) /   2.2155149909669946D-02  /
   data green(32,23,23) /   2.1914408914625958D-02  /
   data green(32,24, 0) /   2.4999407936736707D-02  /
   data green(32,24, 1) /   2.4991593054187231D-02  /
   data green(32,24, 2) /   2.4968192426741972D-02  /
   data green(32,24, 3) /   2.4929337425021819D-02  /
   data green(32,24, 4) /   2.4875244721349749D-02  /
   data green(32,24, 5) /   2.4806212944971484D-02  /
   data green(32,24, 6) /   2.4722618142117934D-02  /
   data green(32,24, 7) /   2.4624908175036559D-02  /
   data green(32,24, 8) /   2.4513596220846982D-02  /
   data green(32,24, 9) /   2.4389253550300020D-02  /
   data green(32,24,10) /   2.4252501777743730D-02  /
   data green(32,24,11) /   2.4104004776839858D-02  /
   data green(32,24,12) /   2.3944460452312504D-02  /
   data green(32,24,13) /   2.3774592547118750D-02  /
   data green(32,24,14) /   2.3595142648067908D-02  /
   data green(32,24,15) /   2.3406862532410179D-02  /
   data green(32,24,16) /   2.3210506974661239D-02  /
   data green(32,24,17) /   2.3006827108279847D-02  /
   data green(32,24,18) /   2.2796564412007666D-02  /
   data green(32,24,19) /   2.2580445366771527D-02  /
   data green(32,24,20) /   2.2359176806885032D-02  /
   data green(32,24,21) /   2.2133441969488424D-02  /
   data green(32,24,22) /   2.1903897229132236D-02  /
   data green(32,24,23) /   2.1671169490338628D-02  /
   data green(32,24,24) /   2.1435854199889873D-02  /
   data green(32,25, 0) /   2.4625103361221107D-02  /
   data green(32,25, 1) /   2.4617634470573638D-02  /
   data green(32,25, 2) /   2.4595268615919541D-02  /
   data green(32,25, 3) /   2.4558127628467041D-02  /
   data green(32,25, 4) /   2.4506412509914723D-02  /
   data green(32,25, 5) /   2.4440400420952788D-02  /
   data green(32,25, 6) /   2.4360440589874750D-02  /
   data green(32,25, 7) /   2.4266949258697539D-02  /
   data green(32,25, 8) /   2.4160403807885353D-02  /
   data green(32,25, 9) /   2.4041336218065711D-02  /
   data green(32,25,10) /   2.3910326037565419D-02  /
   data green(32,25,11) /   2.3767993028150688D-02  /
   data green(32,25,12) /   2.3614989658399279D-02  /
   data green(32,25,13) /   2.3451993605365418D-02  /
   data green(32,25,14) /   2.3279700411572596D-02  /
   data green(32,25,15) /   2.3098816426994247D-02  /
   data green(32,25,16) /   2.2910052145734758D-02  /
   data green(32,25,17) /   2.2714116025761202D-02  /
   data green(32,25,18) /   2.2511708858330683D-02  /
   data green(32,25,19) /   2.2303518732639412D-02  /
   data green(32,25,20) /   2.2090216621448503D-02  /
   data green(32,25,21) /   2.1872452595587524D-02  /
   data green(32,25,22) /   2.1650852659687414D-02  /
   data green(32,25,23) /   2.1426016188459521D-02  /
   data green(32,25,24) /   2.1198513932379065D-02  /
   data green(32,25,25) /   2.0968886553682408D-02  /
   data green(32,26, 0) /   2.4252859535492187D-02  /
   data green(32,26, 1) /   2.4245724421641524D-02  /
   data green(32,26, 2) /   2.4224356900402734D-02  /
   data green(32,26, 3) /   2.4188869874974955D-02  /
   data green(32,26, 4) /   2.4139449674379344D-02  /
   data green(32,26, 5) /   2.4076353344765371D-02  /
   data green(32,26, 6) /   2.3999904966048227D-02  /
   data green(32,26, 7) /   2.3910491096497091D-02  /
   data green(32,26, 8) /   2.3808555468866687D-02  /
   data green(32,26, 9) /   2.3694593077183346D-02  /
   data green(32,26,10) /   2.3569143802950651D-02  /
   data green(32,26,11) /   2.3432785733275212D-02  /
   data green(32,26,12) /   2.3286128321504843D-02  /
   data green(32,26,13) /   2.3129805533985724D-02  /
   data green(32,26,14) /   2.2964469115257728D-02  /
   data green(32,26,15) /   2.2790782089342917D-02  /
   data green(32,26,16) /   2.2609412597730739D-02  /
   data green(32,26,17) /   2.2421028156211539D-02  /
   data green(32,26,18) /   2.2226290393779845D-02  /
   data green(32,26,19) /   2.2025850318239604D-02  /
   data green(32,26,20) /   2.1820344135567826D-02  /
   data green(32,26,21) /   2.1610389634056634D-02  /
   data green(32,26,22) /   2.1396583130121653D-02  /
   data green(32,26,23) /   2.1179496960654628D-02  /
   data green(32,26,24) /   2.0959677496996516D-02  /
   data green(32,26,25) /   2.0737643647989430D-02  /
   data green(32,26,26) /   2.0513885814017554D-02  /
   data green(32,27, 0) /   2.3883379782184237D-02  /
   data green(32,27, 1) /   2.3876565958240111D-02  /
   data green(32,27, 2) /   2.3856159509246726D-02  /
   data green(32,27, 3) /   2.3822265002923788D-02  /
   data green(32,27, 4) /   2.3775055063545228D-02  /
   data green(32,27, 5) /   2.3714767937623921D-02  /
   data green(32,27, 6) /   2.3641704180705221D-02  /
   data green(32,27, 7) /   2.3556222554873645D-02  /
   data green(32,27, 8) /   2.3458735245109066D-02  /
   data green(32,27, 9) /   2.3349702516524011D-02  /
   data green(32,27,10) /   2.3229626943398002D-02  /
   data green(32,27,11) /   2.3099047344727265D-02  /
   data green(32,27,12) /   2.2958532559932934D-02  /
   data green(32,27,13) /   2.2808675192865717D-02  /
   data green(32,27,14) /   2.2650085442944983D-02  /
   data green(32,27,15) /   2.2483385129942534D-02  /
   data green(32,27,16) /   2.2309202004391699D-02  /
   data green(32,27,17) /   2.2128164419714436D-02  /
   data green(32,27,18) /   2.1940896425700736D-02  /
   data green(32,27,19) /   2.1748013326661667D-02  /
   data green(32,27,20) /   2.1550117732007636D-02  /
   data green(32,27,21) /   2.1347796112651633D-02  /
   data green(32,27,22) /   2.1141615863846604D-02  /
   data green(32,27,23) /   2.0932122864055592D-02  /
   data green(32,27,24) /   2.0719839510323999D-02  /
   data green(32,27,25) /   2.0505263203383722D-02  /
   data green(32,27,26) /   2.0288865250293324D-02  /
   data green(32,27,27) /   2.0071090148674576D-02  /
   data green(32,28, 0) /   2.3517282207609902D-02  /
   data green(32,28, 1) /   2.3510777052320762D-02  /
   data green(32,28, 2) /   2.3491294003624263D-02  /
   data green(32,28, 3) /   2.3458929863653298D-02  /
   data green(32,28, 4) /   2.3413844483597034D-02  /
   data green(32,28, 5) /   2.3356258577421424D-02  /
   data green(32,28, 6) /   2.3286450743664224D-02  /
   data green(32,28, 7) /   2.3204753773475625D-02  /
   data green(32,28, 8) /   2.3111550339429195D-02  /
   data green(32,28, 9) /   2.3007268172042905D-02  /
   data green(32,28,10) /   2.2892374839092341D-02  /
   data green(32,28,11) /   2.2767372246581350D-02  /
   data green(32,28,12) /   2.2632790979809712D-02  /
   data green(32,28,13) /   2.2489184598696737D-02  /
   data green(32,28,14) /   2.2337123993899890D-02  /
   data green(32,28,15) /   2.2177191899941277D-02  /
   data green(32,28,16) /   2.2009977649215477D-02  /
   data green(32,28,17) /   2.1836072237109157D-02  /
   data green(32,28,18) /   2.1656063754196112D-02  /
   data green(32,28,19) /   2.1470533227197892D-02  /
   data green(32,28,20) /   2.1280050896650561D-02  /
   data green(32,28,21) /   2.1085172946417746D-02  /
   data green(32,28,22) /   2.0886438688660108D-02  /
   data green(32,28,23) /   2.0684368197821486D-02  /
   data green(32,28,24) /   2.0479460378738931D-02  /
   data green(32,28,25) /   2.0272191447155110D-02  /
   data green(32,28,26) /   2.0063013795665163D-02  /
   data green(32,28,27) /   1.9852355214369167D-02  /
   data green(32,28,28) /   1.9640618433089625D-02  /
   data green(32,29, 0) /   2.3155105371118166D-02  /
   data green(32,29, 1) /   2.3148896246633334D-02  /
   data green(32,29, 2) /   2.3130298867893029D-02  /
   data green(32,29, 3) /   2.3099402815994000D-02  /
   data green(32,29, 4) /   2.3056356059913210D-02  /
   data green(32,29, 5) /   2.3001362993991154D-02  /
   data green(32,29, 6) /   2.2934681762282718D-02  /
   data green(32,29, 7) /   2.2856620937924926D-02  /
   data green(32,29, 8) /   2.2767535640086056D-02  /
   data green(32,29, 9) /   2.2667823182134460D-02  /
   data green(32,29,10) /   2.2557918352097684D-02  /
   data green(32,29,11) /   2.2438288430181922D-02  /
   data green(32,29,12) /   2.2309428048194680D-02  /
   data green(32,29,13) /   2.2171853992437346D-02  /
   data green(32,29,14) /   2.2026100045427183D-02  /
   data green(32,29,15) /   2.1872711953191730D-02  /
   data green(32,29,16) /   2.1712242594429287D-02  /
   data green(32,29,17) /   2.1545247416147304D-02  /
   data green(32,29,18) /   2.1372280188055599D-02  /
   data green(32,29,19) /   2.1193889115536955D-02  /
   data green(32,29,20) /   2.1010613338910938D-02  /
   data green(32,29,21) /   2.0822979835327000D-02  /
   data green(32,29,22) /   2.0631500729266013D-02  /
   data green(32,29,23) /   2.0436671008496367D-02  /
   data green(32,29,24) /   2.0238966634544116D-02  /
   data green(32,29,25) /   2.0038843030339232D-02  /
   data green(32,29,26) /   1.9836733922676424D-02  /
   data green(32,29,27) /   1.9633050513410438D-02  /
   data green(32,29,28) /   1.9428180950787861D-02  /
   data green(32,29,29) /   1.9222490070870679D-02  /
   data green(32,30, 0) /   2.2797314107352366D-02  /
   data green(32,30, 1) /   2.2791388458780943D-02  /
   data green(32,30, 2) /   2.2773639259343472D-02  /
   data green(32,30, 3) /   2.2744149386505232D-02  /
   data green(32,30, 4) /   2.2703055774066588D-02  /
   data green(32,30, 5) /   2.2650547651154495D-02  /
   data green(32,30, 6) /   2.2586864139334306D-02  /
   data green(32,30, 7) /   2.2512291267222730D-02  /
   data green(32,30, 8) /   2.2427158474679326D-02  /
   data green(32,30, 9) /   2.2331834688515604D-02  /
   data green(32,30,10) /   2.2226724058420980D-02  /
   data green(32,30,11) /   2.2112261445373777D-02  /
   data green(32,30,12) /   2.1988907755251873D-02  /
   data green(32,30,13) /   2.1857145207899564D-02  /
   data green(32,30,14) /   2.1717472626883318D-02  /
   data green(32,30,15) /   2.1570400828004526D-02  /
   data green(32,30,16) /   2.1416448175815350D-02  /
   data green(32,30,17) /   2.1256136367406159D-02  /
   data green(32,30,18) /   2.1089986492095363D-02  /
   data green(32,30,19) /   2.0918515404813369D-02  /
   data green(32,30,20) /   2.0742232440340791D-02  /
   data green(32,30,21) /   2.0561636485474672D-02  /
   data green(32,30,22) /   2.0377213416924476D-02  /
   data green(32,30,23) /   2.0189433904472684D-02  /
   data green(32,30,24) /   1.9998751571795828D-02  /
   data green(32,30,25) /   1.9805601501385963D-02  /
   data green(32,30,26) /   1.9610399065242569D-02  /
   data green(32,30,27) /   1.9413539059376404D-02  /
   data green(32,30,28) /   1.9215395117601331D-02  /
   data green(32,30,29) /   1.9016319378486558D-02  /
   data green(32,30,30) /   1.8816642378579145D-02  /
   data green(32,31, 0) /   2.2444305367998820D-02  /
   data green(32,31, 1) /   2.2438650806355029D-02  /
   data green(32,31, 2) /   2.2421712783746155D-02  /
   data green(32,31, 3) /   2.2393567963214071D-02  /
   data green(32,31, 4) /   2.2354343044771550D-02  /
   data green(32,31, 5) /   2.2304213185586515D-02  /
   data green(32,31, 6) /   2.2243399842613070D-02  /
   data green(32,31, 7) /   2.2172168089392857D-02  /
   data green(32,31, 8) /   2.2090823469923825D-02  /
   data green(32,31, 9) /   2.1999708461263010D-02  /
   data green(32,31,10) /   2.1899198622662180D-02  /
   data green(32,31,11) /   2.1789698512438224D-02  /
   data green(32,31,12) /   2.1671637454499636D-02  /
   data green(32,31,13) /   2.1545465234654556D-02  /
   data green(32,31,14) /   2.1411647802786125D-02  /
   data green(32,31,15) /   2.1270663051047318D-02  /
   data green(32,31,16) /   2.1122996730799273D-02  /
   data green(32,31,17) /   2.0969138562516153D-02  /
   data green(32,31,18) /   2.0809578583727031D-02  /
   data green(32,31,19) /   2.0644803770654526D-02  /
   data green(32,31,20) /   2.0475294959896612D-02  /
   data green(32,31,21) /   2.0301524087583660D-02  /
   data green(32,31,22) /   2.0123951755167176D-02  /
   data green(32,31,23) /   1.9943025123541281D-02  /
   data green(32,31,24) /   1.9759176130681625D-02  /
   data green(32,31,25) /   1.9572820022471866D-02  /
   data green(32,31,26) /   1.9384354181893337D-02  /
   data green(32,31,27) /   1.9194157238251356D-02  /
   data green(32,31,28) /   1.9002588435547297D-02  /
   data green(32,31,29) /   1.8809987237398695D-02  /
   data green(32,31,30) /   1.8616673144964756D-02  /
   data green(32,31,31) /   1.8422945704048736D-02  /
   data green(32,32, 0) /   2.2096413976758506D-02  /
   data green(32,32, 1) /   2.2091018346935720D-02  /
   data green(32,32, 2) /   2.2074855190646267D-02  /
   data green(32,32, 3) /   2.2047995417088708D-02  /
   data green(32,32, 4) /   2.2010556247440773D-02  /
   data green(32,32, 5) /   2.1962699797759763D-02  /
   data green(32,32, 6) /   2.1904631142321770D-02  /
   data green(32,32, 7) /   2.1836595902444220D-02  /
   data green(32,32, 8) /   2.1758877415664516D-02  /
   data green(32,32, 9) /   2.1671793547936326D-02  /
   data green(32,32,10) /   2.1575693217052990D-02  /
   data green(32,32,11) /   2.1470952698722987D-02  /
   data green(32,32,12) /   2.1357971787633451D-02  /
   data green(32,32,13) /   2.1237169884572436D-02  /
   data green(32,32,14) /   2.1108982077457369D-02  /
   data green(32,32,15) /   2.0973855279222135D-02  /
   data green(32,32,16) /   2.0832244479277031D-02  /
   data green(32,32,17) /   2.0684609158030418D-02  /
   data green(32,32,18) /   2.0531409906103210D-02  /
   data green(32,32,19) /   2.0373105281713849D-02  /
   data green(32,32,20) /   2.0210148931573677D-02  /
   data green(32,32,21) /   2.0042986992772089D-02  /
   data green(32,32,22) /   1.9872055785767655D-02  /
   data green(32,32,23) /   1.9697779801899604D-02  /
   data green(32,32,24) /   1.9520569982909491D-02  /
   data green(32,32,25) /   1.9340822284884110D-02  /
   data green(32,32,26) /   1.9158916514822103D-02  /
   data green(32,32,27) /   1.8975215424679767D-02  /
   data green(32,32,28) /   1.8790064045226446D-02  /
   data green(32,32,29) /   1.8603789240276076D-02  /
   data green(32,32,30) /   1.8416699460783389D-02  /
   data green(32,32,31) /   1.8229084677816226D-02  /
   data green(32,32,32) /   1.8041216473449240D-02  /
   data green(33, 0, 0) /   3.0310004978728163D-02  /
   data green(33, 1, 0) /   3.0296059616866558D-02  /
   data green(33, 1, 1) /   3.0282133574920291D-02  /
   data green(33, 2, 0) /   3.0254339451224275D-02  /
   data green(33, 2, 1) /   3.0240471099861707D-02  /
   data green(33, 2, 2) /   3.0198980894356627D-02  /
   data green(33, 3, 0) /   3.0185189549966605D-02  /
   data green(33, 3, 1) /   3.0171416461184986D-02  /
   data green(33, 3, 2) /   3.0130210723833963D-02  /
   data green(33, 3, 3) /   3.0061910311239878D-02  /
   data green(33, 4, 0) /   3.0089176177923373D-02  /
   data green(33, 4, 1) /   3.0075534621463815D-02  /
   data green(33, 4, 2) /   3.0034721669438973D-02  /
   data green(33, 4, 3) /   2.9967069918901319D-02  /
   data green(33, 4, 4) /   2.9873125296930859D-02  /
   data green(33, 5, 0) /   2.9967073945978757D-02  /
   data green(33, 5, 1) /   2.9953598427247433D-02  /
   data green(33, 5, 2) /   2.9913281318709000D-02  /
   data green(33, 5, 3) /   2.9846448481229037D-02  /
   data green(33, 5, 4) /   2.9753634853683628D-02  /
   data green(33, 5, 5) /   2.9635572589675908D-02  /
   data green(33, 6, 0) /   2.9819848632418672D-02  /
   data green(33, 6, 1) /   2.9806571485115028D-02  /
   data green(33, 6, 2) /   2.9766846804392986D-02  /
   data green(33, 6, 3) /   2.9700992477192344D-02  /
   data green(33, 6, 4) /   2.9609530432489418D-02  /
   data green(33, 6, 5) /   2.9493175181149422D-02  /
   data green(33, 6, 6) /   2.9352818471841794D-02  /
   data green(33, 7, 0) /   2.9648636415465959D-02  /
   data green(33, 7, 1) /   2.9635587459093247D-02  /
   data green(33, 7, 2) /   2.9596544294940866D-02  /
   data green(33, 7, 3) /   2.9531815737421908D-02  /
   data green(33, 7, 4) /   2.9441908909619216D-02  /
   data green(33, 7, 5) /   2.9327518236898747D-02  /
   data green(33, 7, 6) /   2.9189510702383147D-02  /
   data green(33, 7, 7) /   2.9028907971740320D-02  /
   data green(33, 8, 0) /   2.9454720372412833D-02  /
   data green(33, 8, 1) /   2.9441926641439159D-02  /
   data green(33, 8, 2) /   2.9403645778543170D-02  /
   data green(33, 8, 3) /   2.9340176578108411D-02  /
   data green(33, 8, 4) /   2.9252009806573824D-02  /
   data green(33, 8, 5) /   2.9139817691018872D-02  /
   data green(33, 8, 6) /   2.9004439829156941D-02  /
   data green(33, 8, 7) /   2.8846866093782193D-02  /
   data green(33, 8, 8) /   2.8668217199112205D-02  /
   data green(33, 9, 0) /   2.9239505160213827D-02  /
   data green(33, 9, 1) /   2.9226990708135715D-02  /
   data green(33, 9, 2) /   2.9189544041617634D-02  /
   data green(33, 9, 3) /   2.9127453145503381D-02  /
   data green(33, 9, 4) /   2.9041191135439976D-02  /
   data green(33, 9, 5) /   2.8931406273031625D-02  /
   data green(33, 9, 6) /   2.8798908572734164D-02  /
   data green(33, 9, 7) /   2.8644653537403444D-02  /
   data green(33, 9, 8) /   2.8469723648847735D-02  /
   data green(33, 9, 9) /   2.8275308290263680D-02  /
   data green(33,10, 0) /   2.9004490795883207D-02  /
   data green(33,10, 1) /   2.8992276574072991D-02  /
   data green(33,10, 2) /   2.8955726746697302D-02  /
   data green(33,10, 3) /   2.8895117859575173D-02  /
   data green(33,10, 4) /   2.8810904343512085D-02  /
   data green(33,10, 5) /   2.8703709078148211D-02  /
   data green(33,10, 6) /   2.8574310725714498D-02  /
   data green(33,10, 7) /   2.8423628334499525D-02  /
   data green(33,10, 8) /   2.8252703796079914D-02  /
   data green(33,10, 9) /   2.8062682788867182D-02  /
   data green(33,10,10) /   2.7854794852961207D-02  /
   data green(33,11, 0) /   2.8751246408066908D-02  /
   data green(33,11, 1) /   2.8739350216276421D-02  /
   data green(33,11, 2) /   2.8703750469128052D-02  /
   data green(33,11, 3) /   2.8644711801552203D-02  /
   data green(33,11, 4) /   2.8562669181872204D-02  /
   data green(33,11, 5) /   2.8458219037436130D-02  /
   data green(33,11, 6) /   2.8332107333319447D-02  /
   data green(33,11, 7) /   2.8185215066496686D-02  /
   data green(33,11, 8) /   2.8018541716807187D-02  /
   data green(33,11, 9) /   2.7833187242353934D-02  /
   data green(33,11,10) /   2.7630333220263356D-02  /
   data green(33,11,11) /   2.7411223716200781D-02  /
   data green(33,12, 0) /   2.8481384742845641D-02  /
   data green(33,12, 1) /   2.8469821245745042D-02  /
   data green(33,12, 2) /   2.8435215465713540D-02  /
   data green(33,12, 3) /   2.8377819806303196D-02  /
   data green(33,12, 4) /   2.8298049242192373D-02  /
   data green(33,12, 5) /   2.8196473011331654D-02  /
   data green(33,12, 6) /   2.8073803445658456D-02  /
   data green(33,12, 7) /   2.7930882365694226D-02  /
   data green(33,12, 8) /   2.7768665537883385D-02  /
   data green(33,12, 9) /   2.7588205737535679D-02  /
   data green(33,12,10) /   2.7390634974192654D-02  /
   data green(33,12,11) /   2.7177146422005458D-02  /
   data green(33,12,12) /   2.6948976559058545D-02  /
   data green(33,13, 0) /   2.8196538088445851D-02  /
   data green(33,13, 1) /   2.8185318890690262D-02  /
   data green(33,13, 2) /   2.8151741832476151D-02  /
   data green(33,13, 3) /   2.8096046907466955D-02  /
   data green(33,13, 4) /   2.8018628796956969D-02  /
   data green(33,13, 5) /   2.7920029125754125D-02  /
   data green(33,13, 6) /   2.7800926042044018D-02  /
   data green(33,13, 7) /   2.7662121510310323D-02  /
   data green(33,13, 8) /   2.7504526774568980D-02  /
   data green(33,13, 9) /   2.7329146490766013D-02  /
   data green(33,13,10) /   2.7137062041626872D-02  /
   data green(33,13,11) /   2.6929414536054391D-02  /
   data green(33,13,12) /   2.6707387961610964D-02  /
   data green(33,13,13) /   2.6472192907175514D-02  /
   data green(33,14, 0) /   2.7898336147071514D-02  /
   data green(33,14, 1) /   2.7887469918137939D-02  /
   data green(33,14, 2) /   2.7854947574761011D-02  /
   data green(33,14, 3) /   2.7800996652385876D-02  /
   data green(33,14, 4) /   2.7725991451624920D-02  /
   data green(33,14, 5) /   2.7630445848557218D-02  /
   data green(33,14, 6) /   2.7515003612196155D-02  /
   data green(33,14, 7) /   2.7380426583239935D-02  /
   data green(33,14, 8) /   2.7227581131123304D-02  /
   data green(33,14, 9) /   2.7057423345507724D-02  /
   data green(33,14,10) /   2.6870983433067153D-02  /
   data green(33,14,11) /   2.6669349781986446D-02  /
   data green(33,14,12) /   2.6453653127772998D-02  /
   data green(33,14,13) /   2.6225051208695553D-02  /
   data green(33,14,14) /   2.5984714241949283D-02  /
   data green(33,15, 0) /   2.7588386238944648D-02  /
   data green(33,15, 1) /   2.7577878878342191D-02  /
   data green(33,15, 2) /   2.7546428972210556D-02  /
   data green(33,15, 3) /   2.7494251666164059D-02  /
   data green(33,15, 4) /   2.7421700983555643D-02  /
   data green(33,15, 5) /   2.7329263175352649D-02  /
   data green(33,15, 6) /   2.7217547756302734D-02  /
   data green(33,15, 7) /   2.7087276548163730D-02  /
   data green(33,15, 8) /   2.6939271108536864D-02  /
   data green(33,15, 9) /   2.6774438960455824D-02  /
   data green(33,15,10) /   2.6593759052696512D-02  /
   data green(33,15,11) /   2.6398266874770451D-02  /
   data green(33,15,12) /   2.6189039626104665D-02  /
   data green(33,15,13) /   2.5967181799365194D-02  /
   data green(33,15,14) /   2.5733811487221148D-02  /
   data green(33,15,15) /   2.5490047664233913D-02  /
   data green(33,16, 0) /   2.7268256083764601D-02  /
   data green(33,16, 1) /   2.7258110917125683D-02  /
   data green(33,16, 2) /   2.7227743483423713D-02  /
   data green(33,16, 3) /   2.7177356709109843D-02  /
   data green(33,16, 4) /   2.7107284611080806D-02  /
   data green(33,16, 5) /   2.7017986166797250D-02  /
   data green(33,16, 6) /   2.6910037044222703D-02  /
   data green(33,16, 7) /   2.6784119480883394D-02  /
   data green(33,16, 8) /   2.6641010654169944D-02  /
   data green(33,16, 9) /   2.6481569919106211D-02  /
   data green(33,16,10) /   2.6306725304538352D-02  /
   data green(33,16,11) /   2.6117459654815097D-02  /
   data green(33,16,12) /   2.5914796783515943D-02  /
   data green(33,16,13) /   2.5699787971521889D-02  /
   data green(33,16,14) /   2.5473499097165557D-02  /
   data green(33,16,15) /   2.5236998634991736D-02  /
   data green(33,16,16) /   2.4991346705355497D-02  /
   data green(33,17, 0) /   2.6939459275969493D-02  /
   data green(33,17, 1) /   2.6929677272886442D-02  /
   data green(33,17, 2) /   2.6900395308082237D-02  /
   data green(33,17, 3) /   2.6851804347000196D-02  /
   data green(33,17, 4) /   2.6784218814303571D-02  /
   data green(33,17, 5) /   2.6698070961405210D-02  /
   data green(33,17, 6) /   2.6593903260531012D-02  /
   data green(33,17, 7) /   2.6472359085184589D-02  /
   data green(33,17, 8) /   2.6334171984981242D-02  /
   data green(33,17, 9) /   2.6180153894442920D-02  /
   data green(33,17,10) /   2.6011182629829651D-02  /
   data green(33,17,11) /   2.5828189026009644D-02  /
   data green(33,17,12) /   2.5632144048379643D-02  /
   data green(33,17,13) /   2.5424046185394920D-02  /
   data green(33,17,14) /   2.5204909388329092D-02  /
   data green(33,17,15) /   2.4975751779636294D-02  /
   data green(33,17,16) /   2.4737585302864962D-02  /
   data green(33,17,17) /   2.4491406438310267D-02  /
   data green(33,18, 0) /   2.6603443457866961D-02  /
   data green(33,18, 1) /   2.6594023463040194D-02  /
   data green(33,18, 2) /   2.6565823613388651D-02  /
   data green(33,18, 3) /   2.6519023244362221D-02  /
   data green(33,18, 4) /   2.6453917722300224D-02  /
   data green(33,18, 5) /   2.6370913283968640D-02  /
   data green(33,18, 6) /   2.6270520061595786D-02  /
   data green(33,18, 7) /   2.6153343525975085D-02  /
   data green(33,18, 8) /   2.6020074623842298D-02  /
   data green(33,18, 9) /   2.5871478914929497D-02  /
   data green(33,18,10) /   2.5708385028200463D-02  /
   data green(33,18,11) /   2.5531672756209333D-02  /
   data green(33,18,12) /   2.5342261092644022D-02  /
   data green(33,18,13) /   2.5141096492993391D-02  /
   data green(33,18,14) /   2.4929141604457958D-02  /
   data green(33,18,15) /   2.4707364671465233D-02  /
   data green(33,18,16) /   2.4476729780196358D-02  /
   data green(33,18,17) /   2.4238188061896575D-02  /
   data green(33,18,18) /   2.3992669932594725D-02  /
   data green(33,19, 0) /   2.6261581102068252D-02  /
   data green(33,19, 1) /   2.6252520072273772D-02  /
   data green(33,19, 2) /   2.6225393339978063D-02  /
   data green(33,19, 3) /   2.6180369000467427D-02  /
   data green(33,19, 4) /   2.6117723992530258D-02  /
   data green(33,19, 5) /   2.6037839382867092D-02  /
   data green(33,19, 6) /   2.5941193986548223D-02  /
   data green(33,19, 7) /   2.5828356530959259D-02  /
   data green(33,19, 8) /   2.5699976610126702D-02  /
   data green(33,19, 9) /   2.5556774703157380D-02  /
   data green(33,19,10) /   2.5399531544131656D-02  /
   data green(33,19,11) /   2.5229077131461230D-02  /
   data green(33,19,12) /   2.5046279653554413D-02  /
   data green(33,19,13) /   2.4852034586366919D-02  /
   data green(33,19,14) /   2.4647254189220476D-02  /
   data green(33,19,15) /   2.4432857590521592D-02  /
   data green(33,19,16) /   2.4209761617099684D-02  /
   data green(33,19,17) /   2.3978872482025917D-02  /
   data green(33,19,18) /   2.3741078407910499D-02  /
   data green(33,19,19) /   2.3497243227366255D-02  /
   data green(33,20, 0) /   2.5915162742859280D-02  /
   data green(33,20, 1) /   2.5906455983323324D-02  /
   data green(33,20, 2) /   2.5880388431225670D-02  /
   data green(33,20, 3) /   2.5837117377189984D-02  /
   data green(33,20, 4) /   2.5776902038695811D-02  /
   data green(33,20, 5) /   2.5700099261234344D-02  /
   data green(33,20, 6) /   2.5607157697218486D-02  /
   data green(33,20, 7) /   2.5498610647114176D-02  /
   data green(33,20, 8) /   2.5375067782810495D-02  /
   data green(33,20, 9) /   2.5237205997826626D-02  /
   data green(33,20,10) /   2.5085759641975690D-02  /
   data green(33,20,11) /   2.4921510399755785D-02  /
   data green(33,20,12) /   2.4745277062909670D-02  /
   data green(33,20,13) /   2.4557905429734638D-02  /
   data green(33,20,14) /   2.4360258538666558D-02  /
   data green(33,20,15) /   2.4153207413444476D-02  /
   data green(33,20,16) /   2.3937622463856893D-02  /
   data green(33,20,17) /   2.3714365651625699D-02  /
   data green(33,20,18) /   2.3484283497120315D-02  /
   data green(33,20,19) /   2.3248200970726128D-02  /
   data green(33,20,20) /   2.3006916283887510D-02  /
   data green(33,21, 0) /   2.5565392444694773D-02  /
   data green(33,21, 1) /   2.5557033839598495D-02  /
   data green(33,21, 2) /   2.5532007278641130D-02  /
   data green(33,21, 3) /   2.5490459716931806D-02  /
   data green(33,21, 4) /   2.5432633412748887D-02  /
   data green(33,21, 5) /   2.5358862016954604D-02  /
   data green(33,21, 6) /   2.5269565272811347D-02  /
   data green(33,21, 7) /   2.5165242489790681D-02  /
   data green(33,21, 8) /   2.5046464986904941D-02  /
   data green(33,21, 9) /   2.4913867723521811D-02  /
   data green(33,21,10) /   2.4768140347987515D-02  /
   data green(33,21,11) /   2.4610017896793363D-02  /
   data green(33,21,12) /   2.4440271370183647D-02  /
   data green(33,21,13) /   2.4259698395221803D-02  /
   data green(33,21,14) /   2.4069114165947154D-02  /
   data green(33,21,15) /   2.3869342824107260D-02  /
   data green(33,21,16) /   2.3661209414825024D-02  /
   data green(33,21,17) /   2.3445532521158420D-02  /
   data green(33,21,18) /   2.3223117651354756D-02  /
   data green(33,21,19) /   2.2994751423962938D-02  /
   data green(33,21,20) /   2.2761196569828314D-02  /
   data green(33,21,21) /   2.2523187747042459D-02  /
   data green(33,22, 0) /   2.5213385263237555D-02  /
   data green(33,22, 1) /   2.5205367496193454D-02  /
   data green(33,22, 2) /   2.5181360143221563D-02  /
   data green(33,22, 3) /   2.5141500315940566D-02  /
   data green(33,22, 4) /   2.5086014113811891D-02  /
   data green(33,22, 5) /   2.5015213073493001D-02  /
   data green(33,22, 6) /   2.4929489352141128D-02  /
   data green(33,22, 7) /   2.4829309789391515D-02  /
   data green(33,22, 8) /   2.4715209021351350D-02  /
   data green(33,22, 9) /   2.4587781840349498D-02  /
   data green(33,22,10) /   2.4447675005843026D-02  /
   data green(33,22,11) /   2.4295578714854642D-02  /
   data green(33,22,12) /   2.4132217935156759D-02  /
   data green(33,22,13) /   2.3958343792119450D-02  /
   data green(33,22,14) /   2.3774725181992615D-02  /
   data green(33,22,15) /   2.3582140761871633D-02  /
   data green(33,22,16) /   2.3381371441234956D-02  /
   data green(33,22,17) /   2.3173193473219941D-02  /
   data green(33,22,18) /   2.2958372217057890D-02  /
   data green(33,22,19) /   2.2737656617459343D-02  /
   data green(33,22,20) /   2.2511774423120627D-02  /
   data green(33,22,21) /   2.2281428145559491D-02  /
   data green(33,22,22) /   2.2047291741586018D-02  /
   data green(33,23, 0) /   2.4860166437846679D-02  /
   data green(33,23, 1) /   2.4852481199242719D-02  /
   data green(33,23, 2) /   2.4829468295855233D-02  /
   data green(33,23, 3) /   2.4791255501256580D-02  /
   data green(33,23, 4) /   2.4738053579283055D-02  /
   data green(33,23, 5) /   2.4670153065599047D-02  /
   data green(33,23, 6) /   2.4587919897768290D-02  /
   data green(33,23, 7) /   2.4491790021590239D-02  /
   data green(33,23, 8) /   2.4382263127035616D-02  /
   data green(33,23, 9) /   2.4259895685608875D-02  /
   data green(33,23,10) /   2.4125293471895338D-02  /
   data green(33,23,11) /   2.3979103755417494D-02  /
   data green(33,23,12) /   2.3822007345165954D-02  /
   data green(33,23,13) /   2.3654710659094965D-02  /
   data green(33,23,14) /   2.3477937975555075D-02  /
   data green(33,23,15) /   2.3292424004326986D-02  /
   data green(33,23,16) /   2.3098906892926877D-02  /
   data green(33,23,17) /   2.2898121760456906D-02  /
   data green(33,23,18) /   2.2690794827645890D-02  /
   data green(33,23,19) /   2.2477638188876561D-02  /
   data green(33,23,20) /   2.2259345250737680D-02  /
   data green(33,23,21) /   2.2036586842577241D-02  /
   data green(33,23,22) /   2.1810007988057031D-02  /
   data green(33,23,23) /   2.1580225313026738D-02  /
   data green(33,24, 0) /   2.4506672051327217D-02  /
   data green(33,24, 1) /   2.4499310230385721D-02  /
   data green(33,24, 2) /   2.4477264616375852D-02  /
   data green(33,24, 3) /   2.4440654155548799D-02  /
   data green(33,24, 4) /   2.4389675108776099D-02  /
   data green(33,24, 5) /   2.4324598138520191D-02  /
   data green(33,24, 6) /   2.4245764350102558D-02  /
   data green(33,24, 7) /   2.4153580399836626D-02  /
   data green(33,24, 8) /   2.4048512805400470D-02  /
   data green(33,24, 9) /   2.3931081610542328D-02  /
   data green(33,24,10) /   2.3801853566398742D-02  /
   data green(33,24,11) /   2.3661434995325626D-02  /
   data green(33,24,12) /   2.3510464500534855D-02  /
   data green(33,24,13) /   2.3349605676650589D-02  /
   data green(33,24,14) /   2.3179539963442030D-02  /
   data green(33,24,15) /   2.3000959768505466D-02  /
   data green(33,24,16) /   2.2814561965669083D-02  /
   data green(33,24,17) /   2.2621041855484546D-02  /
   data green(33,24,18) /   2.2421087653369937D-02  /
   data green(33,24,19) /   2.2215375550674120D-02  /
   data green(33,24,20) /   2.2004565374875162D-02  /
   data green(33,24,21) /   2.1789296857860768D-02  /
   data green(33,24,22) /   2.1570186506146435D-02  /
   data green(33,24,23) /   2.1347825054181438D-02  /
   data green(33,24,24) /   2.1122775471642503D-02  /
   data green(33,25, 0) /   2.4153750900176012D-02  /
   data green(33,25, 1) /   2.4146702760404225D-02  /
   data green(33,25, 2) /   2.4125595397806121D-02  /
   data green(33,25, 3) /   2.4090539440452880D-02  /
   data green(33,25, 4) /   2.4041717477105869D-02  /
   data green(33,25, 5) /   2.3979381423960541D-02  /
   data green(33,25, 6) /   2.3903848943028608D-02  /
   data green(33,25, 7) /   2.3815499011177743D-02  /
   data green(33,25, 8) /   2.3714766759147535D-02  /
   data green(33,25, 9) /   2.3602137714932657D-02  /
   data green(33,25,10) /   2.3478141595370518D-02  /
   data green(33,25,11) /   2.3343345793528057D-02  /
   data green(33,25,12) /   2.3198348707811494D-02  /
   data green(33,25,13) /   2.3043773052150253D-02  /
   data green(33,25,14) /   2.2880259275873914D-02  /
   data green(33,25,15) /   2.2708459207886753D-02  /
   data green(33,25,16) /   2.2529030023392294D-02  /
   data green(33,25,17) /   2.2342628613678743D-02  /
   data green(33,25,18) /   2.2149906421229455D-02  /
   data green(33,25,19) /   2.1951504784459952D-02  /
   data green(33,25,20) /   2.1748050819359904D-02  /
   data green(33,25,21) /   2.1540153849743046D-02  /
   data green(33,25,22) /   2.1328402384042291D-02  /
   data green(33,25,23) /   2.1113361624848250D-02  /
   data green(33,25,24) /   2.0895571487768581D-02  /
   data green(33,25,25) /   2.0675545098666634D-02  /
   data green(33,26, 0) /   2.3802167333722681D-02  /
   data green(33,26, 1) /   2.3795422671131841D-02  /
   data green(33,26, 2) /   2.3775223116977294D-02  /
   data green(33,26, 3) /   2.3741671483040197D-02  /
   data green(33,26, 4) /   2.3694937505692744D-02  /
   data green(33,26, 5) /   2.3635255468127662D-02  /
   data green(33,26, 6) /   2.3562920963501380D-02  /
   data green(33,26, 7) /   2.3478286885969597D-02  /
   data green(33,26, 8) /   2.3381758754628730D-02  /
   data green(33,26, 9) /   2.3273789488933841D-02  /
   data green(33,26,10) /   2.3154873762881757D-02  /
   data green(33,26,11) /   2.3025542069049845D-02  /
   data green(33,26,12) /   2.2886354622659373D-02  /
   data green(33,26,13) /   2.2737895230611999D-02  /
   data green(33,26,14) /   2.2580765241538455D-02  /
   data green(33,26,15) /   2.2415577681031907D-02  /
   data green(33,26,16) /   2.2242951662215210D-02  /
   data green(33,26,17) /   2.2063507146418044D-02  /
   data green(33,26,18) /   2.1877860112785131D-02  /
   data green(33,26,19) /   2.1686618179790371D-02  /
   data green(33,26,20) /   2.1490376706476692D-02  /
   data green(33,26,21) /   2.1289715387242641D-02  /
   data green(33,26,22) /   2.1085195341491777D-02  /
   data green(33,26,23) /   2.0877356688665527D-02  /
   data green(33,26,24) /   2.0666716590202805D-02  /
   data green(33,26,25) /   2.0453767732818925D-02  /
   data green(33,26,26) /   2.0238977222108406D-02  /
   data green(33,27, 0) /   2.3452604840977432D-02  /
   data green(33,27, 1) /   2.3446153125016772D-02  /
   data green(33,27, 2) /   2.3426829952607347D-02  /
   data green(33,27, 3) /   2.3394730809313102D-02  /
   data green(33,27, 4) /   2.3350013380167293D-02  /
   data green(33,27, 5) /   2.3292895404558374D-02  /
   data green(33,27, 6) /   2.3223651753669982D-02  /
   data green(33,27, 7) /   2.3142610806786781D-02  /
   data green(33,27, 8) /   2.3050150218774670D-02  /
   data green(33,27, 9) /   2.2946692183212738D-02  /
   data green(33,27,10) /   2.2832698303665879D-02  /
   data green(33,27,11) /   2.2708664189360031D-02  /
   data green(33,27,12) /   2.2575113891190369D-02  /
   data green(33,27,13) /   2.2432594289902894D-02  /
   data green(33,27,14) /   2.2281669540934249D-02  /
   data green(33,27,15) /   2.2122915670388255D-02  /
   data green(33,27,16) /   2.1956915404640211D-02  /
   data green(33,27,17) /   2.1784253302781467D-02  /
   data green(33,27,18) /   2.1605511247207890D-02  /
   data green(33,27,19) /   2.1421264333718659D-02  /
   data green(33,27,20) /   2.1232077189044848D-02  /
   data green(33,27,21) /   2.1038500731187747D-02  /
   data green(33,27,22) /   2.0841069376630907D-02  /
   data green(33,27,23) /   2.0640298688606947D-02  /
   data green(33,27,24) /   2.0436683452267006D-02  /
   data green(33,27,25) /   2.0230696155847907D-02  /
   data green(33,27,26) /   2.0022785851724567D-02  /
   data green(33,27,27) /   1.9813377367478301D-02  /
   data green(33,28, 0) /   2.3105670187532654D-02  /
   data green(33,28, 1) /   2.3099500685131533D-02  /
   data green(33,28, 2) /   2.3081021854969293D-02  /
   data green(33,28, 3) /   2.3050322331071633D-02  /
   data green(33,28, 4) /   2.3007548523590351D-02  /
   data green(33,28, 5) /   2.2952902685034349D-02  /
   data green(33,28, 6) /   2.2886640273521344D-02  /
   data green(33,28, 7) /   2.2809066679937123D-02  /
   data green(33,28, 8) /   2.2720533400058820D-02  /
   data green(33,28, 9) /   2.2621433743600147D-02  /
   data green(33,28,10) /   2.2512198179472351D-02  /
   data green(33,28,11) /   2.2393289420239441D-02  /
   data green(33,28,12) /   2.2265197348874225D-02  /
   data green(33,28,13) /   2.2128433887766044D-02  /
   data green(33,28,14) /   2.1983527903897061D-02  /
   data green(33,28,15) /   2.1831020235696948D-02  /
   data green(33,28,16) /   2.1671458916873443D-02  /
   data green(33,28,17) /   2.1505394661079404D-02  /
   data green(33,28,18) /   2.1333376659186342D-02  /
   data green(33,28,19) /   2.1155948728712438D-02  /
   data green(33,28,20) /   2.0973645843056252D-02  /
   data green(33,28,21) /   2.0786991056992877D-02  /
   data green(33,28,22) /   2.0596492834685301D-02  /
   data green(33,28,23) /   2.0402642777453486D-02  /
   data green(33,28,24) /   2.0205913740847140D-02  /
   data green(33,28,25) /   2.0006758324231912D-02  /
   data green(33,28,26) /   1.9805607711108609D-02  /
   data green(33,28,27) /   1.9602870834677891D-02  /
   data green(33,28,28) /   1.9398933840635200D-02  /
   data green(33,29, 0) /   2.2761897929786393D-02  /
   data green(33,29, 1) /   2.2755999813234132D-02  /
   data green(33,29, 2) /   2.2738332995759839D-02  /
   data green(33,29, 3) /   2.2708979716433857D-02  /
   data green(33,29, 4) /   2.2668075857897249D-02  /
   data green(33,29, 5) /   2.2615809204163152D-02  /
   data green(33,29, 6) /   2.2552417063207542D-02  /
   data green(33,29, 7) /   2.2478183312928016D-02  /
   data green(33,29, 8) /   2.2393434941586823D-02  /
   data green(33,29, 9) /   2.2298538163602159D-02  /
   data green(33,29,10) /   2.2193894198247931D-02  /
   data green(33,29,11) /   2.2079934802375197D-02  /
   data green(33,29,12) /   2.1957117648744941D-02  /
   data green(33,29,13) /   2.1825921639175672D-02  /
   data green(33,29,14) /   2.1686842236789166D-02  /
   data green(33,29,15) /   2.1540386894604264D-02  /
   data green(33,29,16) /   2.1387070649053312D-02  /
   data green(33,29,17) /   2.1227411937173403D-02  /
   data green(33,29,18) /   2.1061928685742779D-02  /
   data green(33,29,19) /   2.0891134709943041D-02  /
   data green(33,29,20) /   2.0715536448632992D-02  /
   data green(33,29,21) /   2.0535630053355812D-02  /
   data green(33,29,22) /   2.0351898839031792D-02  /
   data green(33,29,23) /   2.0164811096105351D-02  /
   data green(33,29,24) /   1.9974818256839395D-02  /
   data green(33,29,25) /   1.9782353402538869D-02  /
   data green(33,29,26) /   1.9587830093742035D-02  /
   data green(33,29,27) /   1.9391641501801535D-02  /
   data green(33,29,28) /   1.9194159817710346D-02  /
   data green(33,29,29) /   1.8995735912411518D-02  /
   data green(33,30, 0) /   2.2421755158719806D-02  /
   data green(33,30, 1) /   2.2416117598351444D-02  /
   data green(33,30, 2) /   2.2399230451355202D-02  /
   data green(33,30, 3) /   2.2371169998389843D-02  /
   data green(33,30, 4) /   2.2332062309617769D-02  /
   data green(33,30, 5) /   2.2282081675824180D-02  /
   data green(33,30, 6) /   2.2221448465875118D-02  /
   data green(33,30, 7) /   2.2150426461784525D-02  /
   data green(33,30, 8) /   2.2069319733743171D-02  /
   data green(33,30, 9) /   2.1978469126165634D-02  /
   data green(33,30,10) /   2.1878248431905181D-02  /
   data green(33,30,11) /   2.1769060335177365D-02  /
   data green(33,30,12) /   2.1651332204466918D-02  /
   data green(33,30,13) /   2.1525511814933486D-02  /
   data green(33,30,14) /   2.1392063075848497D-02  /
   data green(33,30,15) /   2.1251461832733510D-02  /
   data green(33,30,16) /   2.1104191806522891D-02  /
   data green(33,30,17) /   2.0950740723660610D-02  /
   data green(33,30,18) /   2.0791596681975603D-02  /
   data green(33,30,19) /   2.0627244787854812D-02  /
   data green(33,30,20) /   2.0458164090999407D-02  /
   data green(33,30,21) /   2.0284824834205420D-02  /
   data green(33,30,22) /   2.0107686027396283D-02  /
   data green(33,30,23) /   1.9927193347728859D-02  /
   data green(33,30,24) /   1.9743777361117779D-02  /
   data green(33,30,25) /   1.9557852055036344D-02  /
   data green(33,30,26) /   1.9369813667976500D-02  /
   data green(33,30,27) /   1.9180039797457322D-02  /
   data green(33,30,28) /   1.8988888765909134D-02  /
   data green(33,30,29) /   1.8796699222048845D-02  /
   data green(33,30,30) /   1.8603789954406976D-02  /
   data green(33,31, 0) /   2.2085646349472604D-02  /
   data green(33,31, 1) /   2.2080258592284675D-02  /
   data green(33,31, 2) /   2.2064118996324562D-02  /
   data green(33,31, 3) /   2.2037298299048380D-02  /
   data green(33,31, 4) /   2.1999913438643683D-02  /
   data green(33,31, 5) /   2.1952126141696251D-02  /
   data green(33,31, 6) /   2.1894140992998447D-02  /
   data green(33,31, 7) /   2.1826203032351821D-02  /
   data green(33,31, 8) /   2.1748594933005450D-02  /
   data green(33,31, 9) /   2.1661633824131932D-02  /
   data green(33,31,10) /   2.1565667825273107D-02  /
   data green(33,31,11) /   2.1461072363897448D-02  /
   data green(33,31,12) /   2.1348246348125734D-02  /
   data green(33,31,13) /   2.1227608265431661D-02  /
   data green(33,31,14) /   2.1099592274923251D-02  /
   data green(33,31,15) /   2.0964644355944598D-02  /
   data green(33,31,16) /   2.0823218569533553D-02  /
   data green(33,31,17) /   2.0675773482082296D-02  /
   data green(33,31,18) /   2.0522768792725875D-02  /
   data green(33,31,19) /   2.0364662197868585D-02  /
   data green(33,31,20) /   2.0201906518152617D-02  /
   data green(33,31,21) /   2.0034947105345566D-02  /
   data green(33,31,22) /   1.9864219539285796D-02  /
   data green(33,31,23) /   1.9690147618346119D-02  /
   data green(33,31,24) /   1.9513141640969678D-02  /
   data green(33,31,25) /   1.9333596970765882D-02  /
   data green(33,31,26) /   1.9151892873455721D-02  /
   data green(33,31,27) /   1.8968391610613426D-02  /
   data green(33,31,28) /   1.8783437772628395D-02  /
   data green(33,31,29) /   1.8597357831546469D-02  /
   data green(33,31,30) /   1.8410459893368446D-02  /
   data green(33,31,31) /   1.8223033628901042D-02  /
   data green(33,32, 0) /   2.1753918215340954D-02  /
   data green(33,32, 1) /   2.1748769650749740D-02  /
   data green(33,32, 2) /   2.1733345906180417D-02  /
   data green(33,32, 3) /   2.1707712569008542D-02  /
   data green(33,32, 4) /   2.1671978090126505D-02  /
   data green(33,32, 5) /   2.1626292512811510D-02  /
   data green(33,32, 6) /   2.1570845734250218D-02  /
   data green(33,32, 7) /   2.1505865338946630D-02  /
   data green(33,32, 8) /   2.1431614051879433D-02  /
   data green(33,32, 9) /   2.1348386866187216D-02  /
   data green(33,32,10) /   2.1256507905167130D-02  /
   data green(33,32,11) /   2.1156327081387510D-02  /
   data green(33,32,12) /   2.1048216616753022D-02  /
   data green(33,32,13) /   2.0932567486517319D-02  /
   data green(33,32,14) /   2.0809785847688660D-02  /
   data green(33,32,15) /   2.0680289508251896D-02  /
   data green(33,32,16) /   2.0544504488405516D-02  /
   data green(33,32,17) /   2.0402861718884974D-02  /
   data green(33,32,18) /   2.0255793914710427D-02  /
   data green(33,32,19) /   2.0103732655648562D-02  /
   data green(33,32,20) /   1.9947105697581444D-02  /
   data green(33,32,21) /   1.9786334532061925D-02  /
   data green(33,32,22) /   1.9621832204799947D-02  /
   data green(33,32,23) /   1.9454001397818933D-02  /
   data green(33,32,24) /   1.9283232774654973D-02  /
   data green(33,32,25) /   1.9109903583314665D-02  /
   data green(33,32,26) /   1.8934376507791024D-02  /
   data green(33,32,27) /   1.8756998755763660D-02  /
   data green(33,32,28) /   1.8578101367654136D-02  /
   data green(33,32,29) /   1.8397998730424577D-02  /
   data green(33,32,30) /   1.8216988278339240D-02  /
   data green(33,32,31) /   1.8035350362283601D-02  /
   data green(33,32,32) /   1.7853348269082241D-02  /
   data green(33,33, 0) /   2.1426864485145178D-02  /
   data green(33,33, 1) /   2.1421944699133742D-02  /
   data green(33,33, 2) /   2.1407205688451794D-02  /
   data green(33,33, 3) /   2.1382708261126426D-02  /
   data green(33,33, 4) /   2.1348552989053988D-02  /
   data green(33,33, 5) /   2.1304879064074689D-02  /
   data green(33,33, 6) /   2.1251862732372541D-02  /
   data green(33,33, 7) /   2.1189715341500860D-02  /
   data green(33,33, 8) /   2.1118681041954933D-02  /
   data green(33,33, 9) /   2.1039034191364478D-02  /
   data green(33,33,10) /   2.0951076513901649D-02  /
   data green(33,33,11) /   2.0855134070317072D-02  /
   data green(33,33,12) /   2.0751554095126763D-02  /
   data green(33,33,13) /   2.0640701756954066D-02  /
   data green(33,33,14) /   2.0522956896019812D-02  /
   data green(33,33,15) /   2.0398710789463197D-02  /
   data green(33,33,16) /   2.0268362990789400D-02  /
   data green(33,33,17) /   2.0132318284528748D-02  /
   data green(33,33,18) /   1.9990983791407343D-02  /
   data green(33,33,19) /   1.9844766253221457D-02  /
   data green(33,33,20) /   1.9694069520405651D-02  /
   data green(33,33,21) /   1.9539292259194534D-02  /
   data green(33,33,22) /   1.9380825889470564D-02  /
   data green(33,33,23) /   1.9219052759007337D-02  /
   data green(33,33,24) /   1.9054344554957983D-02  /
   data green(33,33,25) /   1.8887060949174680D-02  /
   data green(33,33,26) /   1.8717548470311360D-02  /
   data green(33,33,27) /   1.8546139592669444D-02  /
   data green(33,33,28) /   1.8373152029380516D-02  /
   data green(33,33,29) /   1.8198888215748266D-02  /
   data green(33,33,30) /   1.8023634967349739D-02  /
   data green(33,33,31) /   1.7847663296763926D-02  /
   data green(33,33,32) /   1.7671228372497679D-02  /
   data green(33,33,33) /   1.7494569603745493D-02  /
   data green(34, 0, 0) /   2.9418140869500054D-02  /
   data green(34, 1, 0) /   2.9405391856881482D-02  /
   data green(34, 1, 1) /   2.9392659477445165D-02  /
   data green(34, 2, 0) /   2.9367244618462383D-02  /
   data green(34, 2, 1) /   2.9354561920418180D-02  /
   data green(34, 2, 2) /   2.9316612757137726D-02  /
   data green(34, 3, 0) /   2.9303996370978296D-02  /
   data green(34, 3, 1) /   2.9291395754795814D-02  /
   data green(34, 3, 2) /   2.9253691766336954D-02  /
   data green(34, 3, 3) /   2.9191175863847247D-02  /
   data green(34, 4, 0) /   2.9216135297876193D-02  /
   data green(34, 4, 1) /   2.9203648106494125D-02  /
   data green(34, 4, 2) /   2.9166282920714881D-02  /
   data green(34, 4, 3) /   2.9104326828482387D-02  /
   data green(34, 4, 4) /   2.9018251474177113D-02  /
   data green(34, 5, 0) /   2.9104330105344835D-02  /
   data green(34, 5, 1) /   2.9091986246150628D-02  /
   data green(34, 5, 2) /   2.9055049210313513D-02  /
   data green(34, 5, 3) /   2.8993800601796904D-02  /
   data green(34, 5, 4) /   2.8908703108224090D-02  /
   data green(34, 5, 5) /   2.8800390814824401D-02  /
   data green(34, 6, 0) /   2.8969416022501068D-02  /
   data green(34, 6, 1) /   2.8957243631337514D-02  /
   data green(34, 6, 2) /   2.8920818809830540D-02  /
   data green(34, 6, 3) /   2.8860416658333315D-02  /
   data green(34, 6, 4) /   2.8776489240461486D-02  /
   data green(34, 6, 5) /   2.8669656207850724D-02  /
   data green(34, 6, 6) /   2.8540692215378149D-02  /
   data green(34, 7, 0) /   2.8812377814867885D-02  /
   data green(34, 7, 1) /   2.8800402970127356D-02  /
   data green(34, 7, 2) /   2.8764568291300852D-02  /
   data green(34, 7, 3) /   2.8705141463058706D-02  /
   data green(34, 7, 4) /   2.8622562432500189D-02  /
   data green(34, 7, 5) /   2.8517434384724406D-02  /
   data green(34, 7, 6) /   2.8390511623137497D-02  /
   data green(34, 7, 7) /   2.8242684827346060D-02  /
   data green(34, 8, 0) /   2.8634330473994767D-02  /
   data green(34, 8, 1) /   2.8622576966367979D-02  /
   data green(34, 8, 2) /   2.8587403535090160D-02  /
   data green(34, 8, 3) /   2.8529069654464976D-02  /
   data green(34, 8, 4) /   2.8448001853376301D-02  /
   data green(34, 8, 5) /   2.8344785074742924D-02  /
   data green(34, 8, 6) /   2.8220151065103404D-02  /
   data green(34, 8, 7) /   2.8074964241813147D-02  /
   data green(34, 8, 8) /   2.7910205562022190D-02  /
   data green(34, 9, 0) /   2.8436498298046356D-02  /
   data green(34, 9, 1) /   2.8424987458890331D-02  /
   data green(34, 9, 2) /   2.8390539043688615D-02  /
   data green(34, 9, 3) /   2.8333403645154891D-02  /
   data green(34, 9, 4) /   2.8253993273055934D-02  /
   data green(34, 9, 5) /   2.8152873123487649D-02  /
   data green(34, 9, 6) /   2.8030750512564390D-02  /
   data green(34, 9, 7) /   2.7888461395253215D-02  /
   data green(34, 9, 8) /   2.7726954962880634D-02  /
   data green(34, 9, 9) /   2.7547276856438048D-02  /
   data green(34,10, 0) /   2.8220193088386906D-02  /
   data green(34,10, 1) /   2.8208943677147890D-02  /
   data green(34,10, 2) /   2.8175276373280329D-02  /
   data green(34,10, 3) /   2.8119432342566641D-02  /
   data green(34,10, 4) /   2.8041808179624757D-02  /
   data green(34,10, 5) /   2.7942948105857984D-02  /
   data green(34,10, 6) /   2.7823533472751610D-02  /
   data green(34,10, 7) /   2.7684369963648228D-02  /
   data green(34,10, 8) /   2.7526372955855097D-02  /
   data green(34,10, 9) /   2.7350551546714195D-02  /
   data green(34,10,10) /   2.7157991761549802D-02  /
   data green(34,11, 0) /   2.7986792158996743D-02  /
   data green(34,11, 1) /   2.7975820307824105D-02  /
   data green(34,11, 2) /   2.7942982370781941D-02  /
   data green(34,11, 3) /   2.7888509666302568D-02  /
   data green(34,11, 4) /   2.7812782681283895D-02  /
   data green(34,11, 5) /   2.7716323709741031D-02  /
   data green(34,11, 6) /   2.7599786942214578D-02  /
   data green(34,11, 7) /   2.7463946371115128D-02  /
   data green(34,11, 8) /   2.7309681941712956D-02  /
   data green(34,11, 9) /   2.7137964418320986D-02  /
   data green(34,11,10) /   2.6949839449776188D-02  /
   data green(34,11,11) /   2.6746411308939934D-02  /
   data green(34,12, 0) /   2.7737716795049910D-02  /
   data green(34,12, 1) /   2.7727036006784625D-02  /
   data green(34,12, 2) /   2.7695067844911332D-02  /
   data green(34,12, 3) /   2.7642033481073171D-02  /
   data green(34,12, 4) /   2.7568296799080041D-02  /
   data green(34,12, 5) /   2.7474357479905037D-02  /
   data green(34,12, 6) /   2.7360841685270396D-02  /
   data green(34,12, 7) /   2.7228490677145369D-02  /
   data green(34,12, 8) /   2.7078147770762919D-02  /
   data green(34,12, 9) /   2.6910744056569087D-02  /
   data green(34,12,10) /   2.6727283341230562D-02  /
   data green(34,12,11) /   2.6528826750555286D-02  /
   data green(34,12,12) /   2.6316477410510242D-02  /
   data green(34,13, 0) /   2.7474411711968894D-02  /
   data green(34,13, 1) /   2.7464032906160486D-02  /
   data green(34,13, 2) /   2.7432967216512515D-02  /
   data green(34,13, 3) /   2.7381425482976670D-02  /
   data green(34,13, 4) /   2.7309754677778825D-02  /
   data green(34,13, 5) /   2.7218431436597280D-02  /
   data green(34,13, 6) /   2.7108053336694672D-02  /
   data green(34,13, 7) /   2.6979328231958973D-02  /
   data green(34,13, 8) /   2.6833062010848160D-02  /
   data green(34,13, 9) /   2.6670145178944196D-02  /
   data green(34,13,10) /   2.6491538682559414D-02  /
   data green(34,13,11) /   2.6298259384517662D-02  /
   data green(34,13,12) /   2.6091365580087199D-02  /
   data green(34,13,13) /   2.5871942903282423D-02  /
   data green(34,14, 0) /   2.7198325965320062D-02  /
   data green(34,14, 1) /   2.7188257565850533D-02  /
   data green(34,14, 2) /   2.7158119594215654D-02  /
   data green(34,14, 3) /   2.7108112479893042D-02  /
   data green(34,14, 4) /   2.7038566148346200D-02  /
   data green(34,14, 5) /   2.6949933993085436D-02  /
   data green(34,14, 6) /   2.6842784741886637D-02  /
   data green(34,14, 7) /   2.6717792500572581D-02  /
   data green(34,14, 8) /   2.6575725309613329D-02  /
   data green(34,14, 9) /   2.6417432582373779D-02  /
   data green(34,14,10) /   2.6243831808475677D-02  /
   data green(34,14,11) /   2.6055894902178769D-02  /
   data green(34,14,12) /   2.5854634555838595D-02  /
   data green(34,14,13) /   2.5641090925172486D-02  /
   data green(34,14,14) /   2.5416318929615770D-02  /
   data green(34,15, 0) /   2.6910895653430023D-02  /
   data green(34,15, 1) /   2.6901143710697590D-02  /
   data green(34,15, 2) /   2.6871951614809174D-02  /
   data green(34,15, 3) /   2.6823509402251516D-02  /
   data green(34,15, 4) /   2.6756129973951180D-02  /
   data green(34,15, 5) /   2.6670243498509459D-02  /
   data green(34,15, 6) /   2.6566389854163736D-02  /
   data green(34,15, 7) /   2.6445209367411480D-02  /
   data green(34,15, 8) /   2.6307432153999454D-02  /
   data green(34,15, 9) /   2.6153866399396052D-02  /
   data green(34,15,10) /   2.5985385930285862D-02  /
   data green(34,15,11) /   2.5802917426615491D-02  /
   data green(34,15,12) /   2.5607427606912862D-02  /
   data green(34,15,13) /   2.5399910690426897D-02  /
   data green(34,15,14) /   2.5181376401038260D-02  /
   data green(34,15,15) /   2.4952838733030457D-02  /
   data green(34,16, 0) /   2.6613528646137093D-02  /
   data green(34,16, 1) /   2.6604096986536516D-02  /
   data green(34,16, 2) /   2.6575862280868778D-02  /
   data green(34,16, 3) /   2.6529004275601585D-02  /
   data green(34,16, 4) /   2.6463819009312217D-02  /
   data green(34,16, 5) /   2.6380713633721912D-02  /
   data green(34,16, 6) /   2.6280199414140525D-02  /
   data green(34,16, 7) /   2.6162883143076597D-02  /
   data green(34,16, 8) /   2.6029457244583380D-02  /
   data green(34,16, 9) /   2.5880688876183531D-02  /
   data green(34,16,10) /   2.5717408349304027D-02  /
   data green(34,16,11) /   2.5540497188487679D-02  /
   data green(34,16,12) /   2.5350876135602350D-02  /
   data green(34,16,13) /   2.5149493379932174D-02  /
   data green(34,16,14) /   2.4937313260976507D-02  /
   data green(34,16,15) /   2.4715305650782456D-02  /
   data green(34,16,16) /   2.4484436179448253D-02  /
   data green(34,17, 0) /   2.6307591470829436D-02  /
   data green(34,17, 1) /   2.6298481866403083D-02  /
   data green(34,17, 2) /   2.6271209927328124D-02  /
   data green(34,17, 3) /   2.6225945287281525D-02  /
   data green(34,17, 4) /   2.6162967406430756D-02  /
   data green(34,17, 5) /   2.6082660793958466D-02  /
   data green(34,17, 6) /   2.5985508545774819D-02  /
   data green(34,17, 7) /   2.5872084408422352D-02  /
   data green(34,17, 8) /   2.5743043620212021D-02  /
   data green(34,17, 9) /   2.5599112807787146D-02  /
   data green(34,17,10) /   2.5441079229971447D-02  /
   data green(34,17,11) /   2.5269779661230337D-02  /
   data green(34,17,12) /   2.5086089195504771D-02  /
   data green(34,17,13) /   2.4890910229348334D-02  /
   data green(34,17,14) /   2.4685161853438133D-02  /
   data green(34,17,15) /   2.4469769846066771D-02  /
   data green(34,17,16) /   2.4245657423597380D-02  /
   data green(34,17,17) /   2.4013736863341399D-02  /
   data green(34,18, 0) /   2.5994398395460357D-02  /
   data green(34,18, 1) /   2.5985610747012333D-02  /
   data green(34,18, 2) /   2.5959301358333676D-02  /
   data green(34,18, 3) /   2.5915629990514649D-02  /
   data green(34,18, 4) /   2.5854859912687392D-02  /
   data green(34,18, 5) /   2.5777353507489772D-02  /
   data green(34,18, 6) /   2.5683566321786342D-02  /
   data green(34,18, 7) /   2.5574039752422385D-02  /
   data green(34,18, 8) /   2.5449392593231996D-02  /
   data green(34,18, 9) /   2.5310311694616540D-02  /
   data green(34,18,10) /   2.5157542000156134D-02  /
   data green(34,18,11) /   2.4991876226128893D-02  /
   data green(34,18,12) /   2.4814144440435153D-02  /
   data green(34,18,13) /   2.4625203778767292D-02  /
   data green(34,18,14) /   2.4425928509848188D-02  /
   data green(34,18,15) /   2.4217200630295779D-02  /
   data green(34,18,16) /   2.3999901135302325D-02  /
   data green(34,18,17) /   2.3774902075860079D-02  /
   data green(34,18,18) /   2.3543059478484535D-02  /
   data green(34,19, 0) /   2.5675202670529796D-02  /
   data green(34,19, 1) /   2.5666735198120322D-02  /
   data green(34,19, 2) /   2.5641383118846638D-02  /
   data green(34,19, 3) /   2.5599296613401604D-02  /
   data green(34,19, 4) /   2.5540723232816531D-02  /
   data green(34,19, 5) /   2.5466003866712460D-02  /
   data green(34,19, 6) /   2.5375567280570561D-02  /
   data green(34,19, 7) /   2.5269923392137662D-02  /
   data green(34,19, 8) /   2.5149655490146892D-02  /
   data green(34,19, 9) /   2.5015411621633597D-02  /
   data green(34,19,10) /   2.4867895386692013D-02  /
   data green(34,19,11) /   2.4707856381685434D-02  /
   data green(34,19,12) /   2.4536080524460997D-02  /
   data green(34,19,13) /   2.4353380479303716D-02  /
   data green(34,19,14) /   2.4160586376830974D-02  /
   data green(34,19,15) /   2.3958536996610956D-02  /
   data green(34,19,16) /   2.3748071549856718D-02  /
   data green(34,19,17) /   2.3530022167887028D-02  /
   data green(34,19,18) /   2.3305207170734508D-02  /
   data green(34,19,19) /   2.3074425160639022D-02  /
   data green(34,20, 0) /   2.5351189829534958D-02  /
   data green(34,20, 1) /   2.5343039265029482D-02  /
   data green(34,20, 2) /   2.5318634803502941D-02  /
   data green(34,20, 3) /   2.5278117378987854D-02  /
   data green(34,20, 4) /   2.5221719365911242D-02  /
   data green(34,20, 5) /   2.5149760888952436D-02  /
   data green(34,20, 6) /   2.5062644818963600D-02  /
   data green(34,20, 7) /   2.4960850606984213D-02  /
   data green(34,20, 8) /   2.4844927138282489D-02  /
   data green(34,20, 9) /   2.4715484809542956D-02  /
   data green(34,20,10) /   2.4573187044252549D-02  /
   data green(34,20,11) /   2.4418741464091620D-02  /
   data green(34,20,12) /   2.4252890928329045D-02  /
   data green(34,20,13) /   2.4076404639921868D-02  /
   data green(34,20,14) /   2.3890069497621000D-02  /
   data green(34,20,15) /   2.3694681849460629D-02  /
   data green(34,20,16) /   2.3491039776191814D-02  /
   data green(34,20,17) /   2.3279936005075833D-02  /
   data green(34,20,18) /   2.3062151526389126D-02  /
   data green(34,20,19) /   2.2838449958197345D-02  /
   data green(34,20,20) /   2.2609572680362518D-02  /
   data green(34,21, 0) /   2.5023472900301463D-02  /
   data green(34,21, 1) /   2.5015634677476833D-02  /
   data green(34,21, 2) /   2.4992164258426854D-02  /
   data green(34,21, 3) /   2.4953193696120330D-02  /
   data green(34,21, 4) /   2.4898940782636287D-02  /
   data green(34,21, 5) /   2.4829705678918765D-02  /
   data green(34,21, 6) /   2.4745866340673457D-02  /
   data green(34,21, 7) /   2.4647872875902966D-02  /
   data green(34,21, 8) /   2.4536240996540920D-02  /
   data green(34,21, 9) /   2.4411544746008811D-02  /
   data green(34,21,10) /   2.4274408695786853D-02  /
   data green(34,21,11) /   2.4125499807276747D-02  /
   data green(34,21,12) /   2.3965519150842593D-02  /
   data green(34,21,13) /   2.3795193662827889D-02  /
   data green(34,21,14) /   2.3615268104741807D-02  /
   data green(34,21,15) /   2.3426497368037241D-02  /
   data green(34,21,16) /   2.3229639244374524D-02  /
   data green(34,21,17) /   2.3025447756353509D-02  /
   data green(34,21,18) /   2.2814667118646631D-02  /
   data green(34,21,19) /   2.2598026375347517D-02  /
   data green(34,21,20) /   2.2376234737013979D-02  /
   data green(34,21,21) /   2.2149977620952409D-02  /
   data green(34,22, 0) /   2.4693089347108149D-02  /
   data green(34,22, 1) /   2.4685557785659087D-02  /
   data green(34,22, 2) /   2.4663004499255840D-02  /
   data green(34,22, 3) /   2.4625553048450155D-02  /
   data green(34,22, 4) /   2.4573407275599281D-02  /
   data green(34,22, 5) /   2.4506848232704770D-02  /
   data green(34,22, 6) /   2.4426230008404784D-02  /
   data green(34,22, 7) /   2.4331974574583186D-02  /
   data green(34,22, 8) /   2.4224565797307061D-02  /
   data green(34,22, 9) /   2.4104542774447027D-02  /
   data green(34,22,10) /   2.3972492672912196D-02  /
   data green(34,22,11) /   2.3829043241922483D-02  /
   data green(34,22,12) /   2.3674855175540910D-02  /
   data green(34,22,13) /   2.3510614488525906D-02  /
   data green(34,22,14) /   2.3337025055432892D-02  /
   data green(34,22,15) /   2.3154801444944521D-02  /
   data green(34,22,16) /   2.2964662160855305D-02  /
   data green(34,22,17) /   2.2767323379175364D-02  /
   data green(34,22,18) /   2.2563493248547773D-02  /
   data green(34,22,19) /   2.2353866799550027D-02  /
   data green(34,22,20) /   2.2139121488247301D-02  /
   data green(34,22,21) /   2.1919913381162143D-02  /
   data green(34,22,22) /   2.1696873973010504D-02  /
   data green(34,23, 0) /   2.4360999544118800D-02  /
   data green(34,23, 1) /   2.4353768024721698D-02  /
   data green(34,23, 2) /   2.4332112149124143D-02  /
   data green(34,23, 3) /   2.4296147389301363D-02  /
   data green(34,23, 4) /   2.4246064296032533D-02  /
   data green(34,23, 5) /   2.4182125703258686D-02  /
   data green(34,23, 6) /   2.4104662927556770D-02  /
   data green(34,23, 7) /   2.4014071069591802D-02  /
   data green(34,23, 8) /   2.3910803546162668D-02  /
   data green(34,23, 9) /   2.3795365997484389D-02  /
   data green(34,23,10) /   2.3668309724229300D-02  /
   data green(34,23,11) /   2.3530224812530840D-02  /
   data green(34,23,12) /   2.3381733102945758D-02  /
   data green(34,23,13) /   2.3223481151872994D-02  /
   data green(34,23,14) /   2.3056133321968458D-02  /
   data green(34,23,15) /   2.2880365122654386D-02  /
   data green(34,23,16) /   2.2696856903939813D-02  /
   data green(34,23,17) /   2.2506287987483151D-02  /
   data green(34,23,18) /   2.2309331299103081D-02  /
   data green(34,23,19) /   2.2106648547630045D-02  /
   data green(34,23,20) /   2.1898885976786537D-02  /
   data green(34,23,21) /   2.1686670700227759D-02  /
   data green(34,23,22) /   2.1470607615335788D-02  /
   data green(34,23,23) /   2.1251276879057521D-02  /
   data green(34,24, 0) /   2.4028086572432394D-02  /
   data green(34,24, 1) /   2.4021147700777120D-02  /
   data green(34,24, 2) /   2.4000367191917227D-02  /
   data green(34,24, 3) /   2.3965852841409504D-02  /
   data green(34,24, 4) /   2.3917782580848920D-02  /
   data green(34,24, 5) /   2.3856401937711384D-02  /
   data green(34,24, 6) /   2.3782020579343706D-02  /
   data green(34,24, 7) /   2.3695008035712204D-02  /
   data green(34,24, 8) /   2.3595788714993710D-02  /
   data green(34,24, 9) /   2.3484836340613095D-02  /
   data green(34,24,10) /   2.3362667947514921D-02  /
   data green(34,24,11) /   2.3229837579234722D-02  /
   data green(34,24,12) /   2.3086929825942969D-02  /
   data green(34,24,13) /   2.2934553337562599D-02  /
   data green(34,24,14) /   2.2773334435999153D-02  /
   data green(34,24,15) /   2.2603910937297720D-02  /
   data green(34,24,16) /   2.2426926279044038D-02  /
   data green(34,24,17) /   2.2243024031452754D-02  /
   data green(34,24,18) /   2.2052842853177934D-02  /
   data green(34,24,19) /   2.1857011935691555D-02  /
   data green(34,24,20) /   2.1656146963733438D-02  /
   data green(34,24,21) /   2.1450846604334280D-02  /
   data green(34,24,22) /   2.1241689523604820D-02  /
   data green(34,24,23) /   2.1029231919086060D-02  /
   data green(34,24,24) /   2.0814005546065054D-02  /
   data green(34,25, 0) /   2.3695157133976626D-02  /
   data green(34,25, 1) /   2.3688502892352019D-02  /
   data green(34,25, 2) /   2.3668573836706935D-02  /
   data green(34,25, 3) /   2.3635470500749194D-02  /
   data green(34,25, 4) /   2.3589358873838197D-02  /
   data green(34,25, 5) /   2.3530468096028146D-02  /
   data green(34,25, 6) /   2.3459087319554671D-02  /
   data green(34,25, 7) /   2.3375561820386543D-02  /
   data green(34,25, 8) /   2.3280288460860173D-02  /
   data green(34,25, 9) /   2.3173710617536026D-02  /
   data green(34,25,10) /   2.3056312696916464D-02  /
   data green(34,25,11) /   2.2928614365460141D-02  /
   data green(34,25,12) /   2.2791164619597151D-02  /
   data green(34,25,13) /   2.2644535816588009D-02  /
   data green(34,25,14) /   2.2489317778653407D-02  /
   data green(34,25,15) /   2.2326112071522682D-02  /
   data green(34,25,16) /   2.2155526545169179D-02  /
   data green(34,25,17) /   2.1978170209786940D-02  /
   data green(34,25,18) /   2.1794648504752261D-02  /
   data green(34,25,19) /   2.1605559003066405D-02  /
   data green(34,25,20) /   2.1411487579156056D-02  /
   data green(34,25,21) /   2.1213005054365134D-02  /
   data green(34,25,22) /   2.1010664322336801D-02  /
   data green(34,25,23) /   2.0804997945972077D-02  /
   data green(34,25,24) /   2.0596516208872227D-02  /
   data green(34,25,25) /   2.0385705597142395D-02  /
   data green(34,26, 0) /   2.3362943383426656D-02  /
   data green(34,26, 1) /   2.3356565269028597D-02  /
   data green(34,26, 2) /   2.3337462296875430D-02  /
   data green(34,26, 3) /   2.3305728150830235D-02  /
   data green(34,26, 4) /   2.3261517551336108D-02  /
   data green(34,26, 5) /   2.3205044166287568D-02  /
   data green(34,26, 6) /   2.3136577764164182D-02  /
   data green(34,26, 7) /   2.3056440683225835D-02  /
   data green(34,26, 8) /   2.2965003706079987D-02  /
   data green(34,26, 9) /   2.2862681440763807D-02  /
   data green(34,26,10) /   2.2749927317319757D-02  /
   data green(34,26,11) /   2.2627228312594626D-02  /
   data green(34,26,12) /   2.2495099515787625D-02  /
   data green(34,26,13) /   2.2354078643435256D-02  /
   data green(34,26,14) /   2.2204720605520856D-02  /
   data green(34,26,15) /   2.2047592214819816D-02  /
   data green(34,26,16) /   2.1883267120078095D-02  /
   data green(34,26,17) /   2.1712321030836146D-02  /
   data green(34,26,18) /   2.1535327288285149D-02  /
   data green(34,26,19) /   2.1352852823062200D-02  /
   data green(34,26,20) /   2.1165454527854486D-02  /
   data green(34,26,21) /   2.0973676060499921D-02  /
   data green(34,26,22) /   2.0778045082247603D-02  /
   data green(34,26,23) /   2.0579070926188776D-02  /
   data green(34,26,24) /   2.0377242682700303D-02  /
   data green(34,26,25) /   2.0173027682098346D-02  /
   data green(34,26,26) /   1.9966870349545876D-02  /
   data green(34,27, 0) /   2.3032105492391385D-02  /
   data green(34,27, 1) /   2.3025994642012441D-02  /
   data green(34,27, 2) /   2.3007691300120828D-02  /
   data green(34,27, 3) /   2.2977282706068299D-02  /
   data green(34,27, 4) /   2.2934912974303914D-02  /
   data green(34,27, 5) /   2.2880781202727012D-02  /
   data green(34,27, 6) /   2.2815138892946816D-02  /
   data green(34,27, 7) /   2.2738286747495579D-02  /
   data green(34,27, 8) /   2.2650570922851911D-02  /
   data green(34,27, 9) /   2.2552378827777778D-02  /
   data green(34,27,10) /   2.2444134563670372D-02  /
   data green(34,27,11) /   2.2326294107287304D-02  /
   data green(34,27,12) /   2.2199340336412715D-02  /
   data green(34,27,13) /   2.2063777996048958D-02  /
   data green(34,27,14) /   2.1920128696933817D-02  /
   data green(34,27,15) /   2.1768926030083673D-02  /
   data green(34,27,16) /   2.1610710871191156D-02  /
   data green(34,27,17) /   2.1446026937628836D-02  /
   data green(34,27,18) /   2.1275416649074962D-02  /
   data green(34,27,19) /   2.1099417330894713D-02  /
   data green(34,27,20) /   2.0918557787822083D-02  /
   data green(34,27,21) /   2.0733355264563574D-02  /
   data green(34,27,22) /   2.0544312799967004D-02  /
   data green(34,27,23) /   2.0351916972569179D-02  /
   data green(34,27,24) /   2.0156636027773178D-02  /
   data green(34,27,25) /   1.9958918370661254D-02  /
   data green(34,27,26) /   1.9759191403510786D-02  /
   data green(34,27,27) /   1.9557860683390330D-02  /
   data green(34,28, 0) /   2.2703234776559059D-02  /
   data green(34,28, 1) /   2.2697382077715496D-02  /
   data green(34,28, 2) /   2.2679851161621162D-02  /
   data green(34,28, 3) /   2.2650723218469900D-02  /
   data green(34,28, 4) /   2.2610132403774526D-02  /
   data green(34,28, 5) /   2.2558264126951368D-02  /
   data green(34,28, 6) /   2.2495352715614914D-02  /
   data green(34,28, 7) /   2.2421678512852105D-02  /
   data green(34,28, 8) /   2.2337564477024945D-02  /
   data green(34,28, 9) /   2.2243372363211864D-02  /
   data green(34,28,10) /   2.2139498571988339D-02  /
   data green(34,28,11) /   2.2026369754773522D-02  /
   data green(34,28,12) /   2.1904438265494774D-02  /
   data green(34,28,13) /   2.1774177546049683D-02  /
   data green(34,28,14) /   2.1636077528296074D-02  /
   data green(34,28,15) /   2.1490640128477328D-02  /
   data green(34,28,16) /   2.1338374901555251D-02  /
   data green(34,28,17) /   2.1179794913352538D-02  /
   data green(34,28,18) /   2.1015412878179405D-02  /
   data green(34,28,19) /   2.0845737599172385D-02  /
   data green(34,28,20) /   2.0671270738303747D-02  /
   data green(34,28,21) /   2.0492503933254413D-02  /
   data green(34,28,22) /   2.0309916269342618D-02  /
   data green(34,28,23) /   2.0123972106654633D-02  /
   data green(34,28,24) /   1.9935119255552160D-02  /
   data green(34,28,25) /   1.9743787487893804D-02  /
   data green(34,28,26) /   1.9550387366610092D-02  /
   data green(34,28,27) /   1.9355309372674462D-02  /
   data green(34,28,28) /   1.9158923305945093D-02  /
   data green(34,29, 0) /   2.2376857234879418D-02  /
   data green(34,29, 1) /   2.2371253423742195D-02  /
   data green(34,29, 2) /   2.2354467270673507D-02  /
   data green(34,29, 3) /   2.2326574299490535D-02  /
   data green(34,29, 4) /   2.2287699333665100D-02  /
   data green(34,29, 5) /   2.2238014949024398D-02  /
   data green(34,29, 6) /   2.2177739360474884D-02  /
   data green(34,29, 7) /   2.2107133793127124D-02  /
   data green(34,29, 8) /   2.2026499399102697D-02  /
   data green(34,29, 9) /   2.1936173789877920D-02  /
   data green(34,29,10) /   2.1836527260035016D-02  /
   data green(34,29,11) /   2.1727958781659524D-02  /
   data green(34,29,12) /   2.1610891849382911D-02  /
   data green(34,29,13) /   2.1485770254383023D-02  /
   data green(34,29,14) /   2.1353053861782240D-02  /
   data green(34,29,15) /   2.1213214460160921D-02  /
   data green(34,29,16) /   2.1066731744716585D-02  /
   data green(34,29,17) /   2.0914089487356957D-02  /
   data green(34,29,18) /   2.0755771938122673D-02  /
   data green(34,29,19) /   2.0592260493178270D-02  /
   data green(34,29,20) /   2.0424030655532444D-02  /
   data green(34,29,21) /   2.0251549305944715D-02  /
   data green(34,29,22) /   2.0075272293381949D-02  /
   data green(34,29,23) /   1.9895642347082153D-02  /
   data green(34,29,24) /   1.9713087305883118D-02  /
   data green(34,29,25) /   1.9528018655044287D-02  /
   data green(34,29,26) /   1.9340830356350237D-02  /
   data green(34,29,27) /   1.9151897953811173D-02  /
   data green(34,29,28) /   1.8961577934716832D-02  /
   data green(34,29,29) /   1.8770207324080013D-02  /
   data green(34,30, 0) /   2.2053437369004558D-02  /
   data green(34,30, 1) /   2.2048073115733896D-02  /
   data green(34,30, 2) /   2.2032003859962006D-02  /
   data green(34,30, 3) /   2.2005299827378574D-02  /
   data green(34,30, 4) /   2.1968077112885027D-02  /
   data green(34,30, 5) /   2.1920496282438265D-02  /
   data green(34,30, 6) /   2.1862760462101651D-02  /
   data green(34,30, 7) /   2.1795112958586278D-02  /
   data green(34,30, 8) /   2.1717834465237191D-02  /
   data green(34,30, 9) /   2.1631239915095708D-02  /
   data green(34,30,10) /   2.1535675048146023D-02  /
   data green(34,30,11) /   2.1431512763045609D-02  /
   data green(34,30,12) /   2.1319149324568108D-02  /
   data green(34,30,13) /   2.1199000496778517D-02  /
   data green(34,30,14) /   2.1071497668828189D-02  /
   data green(34,30,15) /   2.0937084035476252D-02  /
   data green(34,30,16) /   2.0796210888340976D-02  /
   data green(34,30,17) /   2.0649334066803311D-02  /
   data green(34,30,18) /   2.0496910609774141D-02  /
   data green(34,30,19) /   2.0339395641528826D-02  /
   data green(34,30,20) /   2.0177239516809754D-02  /
   data green(34,30,21) /   2.0010885242660783D-02  /
   data green(34,30,22) /   1.9840766187201268D-02  /
   data green(34,30,23) /   1.9667304078934921D-02  /
   data green(34,30,24) /   1.9490907294336519D-02  /
   data green(34,30,25) /   1.9311969426433047D-02  /
   data green(34,30,26) /   1.9130868122925009D-02  /
   data green(34,30,27) /   1.8947964179066853D-02  /
   data green(34,30,28) /   1.8763600868008563D-02  /
   data green(34,30,29) /   1.8578103489532382D-02  /
   data green(34,30,30) /   1.8391779117028093D-02  /
   data green(34,31, 0) /   2.1733382170206646D-02  /
   data green(34,31, 1) /   2.1728248152454645D-02  /
   data green(34,31, 2) /   2.1712867945332486D-02  /
   data green(34,31, 3) /   2.1687306828709597D-02  /
   data green(34,31, 4) /   2.1651672746601631D-02  /
   data green(34,31, 5) /   2.1606115044313170D-02  /
   data green(34,31, 6) /   2.1550822741204137D-02  /
   data green(34,31, 7) /   2.1486022377983003D-02  /
   data green(34,31, 8) /   2.1411975486002575D-02  /
   data green(34,31, 9) /   2.1328975732897115D-02  /
   data green(34,31,10) /   2.1237345803874746D-02  /
   data green(34,31,11) /   2.1137434080982721D-02  /
   data green(34,31,12) /   2.1029611183705662D-02  /
   data green(34,31,13) /   2.0914266433435573D-02  /
   data green(34,31,14) /   2.0791804301840008D-02  /
   data green(34,31,15) /   2.0662640899178034D-02  /
   data green(34,31,16) /   2.0527200553446725D-02  /
   data green(34,31,17) /   2.0385912525172543D-02  /
   data green(34,31,18) /   2.0239207895991346D-02  /
   data green(34,31,19) /   2.0087516662175110D-02  /
   data green(34,31,20) /   1.9931265057224107D-02  /
   data green(34,31,21) /   1.9770873120783971D-02  /
   data green(34,31,22) /   1.9606752524659633D-02  /
   data green(34,31,23) /   1.9439304660732161D-02  /
   data green(34,31,24) /   1.9268918990251693D-02  /
   data green(34,31,25) /   1.9095971649347010D-02  /
   data green(34,31,26) /   1.8920824301694127D-02  /
   data green(34,31,27) /   1.8743823226124078D-02  /
   data green(34,31,28) /   1.8565298624499873D-02  /
   data green(34,31,29) /   1.8385564133409808D-02  /
   data green(34,31,30) /   1.8204916522052177D-02  /
   data green(34,31,31) /   1.8023635558053419D-02  /
   data green(34,32, 0) /   2.1417045179169942D-02  /
   data green(34,32, 1) /   2.1412132144622908D-02  /
   data green(34,32, 2) /   2.1397413341901567D-02  /
   data green(34,32, 3) /   2.1372949440480779D-02  /
   data green(34,32, 4) /   2.1338840783800481D-02  /
   data green(34,32, 5) /   2.1295226248961619D-02  /
   data green(34,32, 6) /   2.1242281686050402D-02  /
   data green(34,32, 7) /   2.1180217971253249D-02  /
   data green(34,32, 8) /   2.1109278715517769D-02  /
   data green(34,32, 9) /   2.1029737676643998D-02  /
   data green(34,32,10) /   2.0941895927200277D-02  /
   data green(34,32,11) /   2.0846078833467330D-02  /
   data green(34,32,12) /   2.0742632901727384D-02  /
   data green(34,32,13) /   2.0631922547704015D-02  /
   data green(34,32,14) /   2.0514326842962286D-02  /
   data green(34,32,15) /   2.0390236288787288D-02  /
   data green(34,32,16) /   2.0260049663696211D-02  /
   data green(34,32,17) /   2.0124170985552407D-02  /
   data green(34,32,18) /   1.9983006623492554D-02  /
   data green(34,32,19) /   1.9836962588795578D-02  /
   data green(34,32,20) /   1.9686442027646221D-02  /
   data green(34,32,21) /   1.9531842932677939D-02  /
   data green(34,32,22) /   1.9373556084394814D-02  /
   data green(34,32,23) /   1.9211963228205581D-02  /
   data green(34,32,24) /   1.9047435487959025D-02  /
   data green(34,32,25) /   1.8880332012617285D-02  /
   data green(34,32,26) /   1.8710998849078595D-02  /
   data green(34,32,27) /   1.8539768031174155D-02  /
   data green(34,32,28) /   1.8366956872501179D-02  /
   data green(34,32,29) /   1.8192867448983286D-02  /
   data green(34,32,30) /   1.8017786255825698D-02  /
   data green(34,32,31) /   1.7841984022798917D-02  /
   data green(34,32,32) /   1.7665715671481598D-02  /
   data green(34,33, 0) /   2.1104730540960834D-02  /
   data green(34,33, 1) /   2.1100029359852984D-02  /
   data green(34,33, 2) /   2.1085944679044861D-02  /
   data green(34,33, 3) /   2.1062532875617844D-02  /
   data green(34,33, 4) /   2.1029887214434884D-02  /
   data green(34,33, 5) /   2.0988136818699283D-02  /
   data green(34,33, 6) /   2.0937445260079610D-02  /
   data green(34,33, 7) /   2.0878008798392186D-02  /
   data green(34,33, 8) /   2.0810054307552717D-02  /
   data green(34,33, 9) /   2.0733836929977906D-02  /
   data green(34,33,10) /   2.0649637505697035D-02  /
   data green(34,33,11) /   2.0557759825049503D-02  /
   data green(34,33,12) /   2.0458527754991702D-02  /
   data green(34,33,13) /   2.0352282288771565D-02  /
   data green(34,33,14) /   2.0239378567162142D-02  /
   data green(34,33,15) /   2.0120182916733061D-02  /
   data green(34,33,16) /   1.9995069946966697D-02  /
   data green(34,33,17) /   1.9864419743602868D-02  /
   data green(34,33,18) /   1.9728615190636505D-02  /
   data green(34,33,19) /   1.9588039448109280D-02  /
   data green(34,33,20) /   1.9443073607427523D-02  /
   data green(34,33,21) /   1.9294094540583555D-02  /
   data green(34,33,22) /   1.9141472954510810D-02  /
   data green(34,33,23) /   1.8985571656988789D-02  /
   data green(34,33,24) /   1.8826744036127975D-02  /
   data green(34,33,25) /   1.8665332751574705D-02  /
   data green(34,33,26) /   1.8501668632220986D-02  /
   data green(34,33,27) /   1.8336069772398756D-02  /
   data green(34,33,28) /   1.8168840816278529D-02  /
   data green(34,33,29) /   1.8000272418455105D-02  /
   data green(34,33,30) /   1.7830640867453854D-02  /
   data green(34,33,31) /   1.7660207858084472D-02  /
   data green(34,33,32) /   1.7489220398156070D-02  /
   data green(34,33,33) /   1.7317910834993313D-02  /
   data green(34,34, 0) /   2.0796696992829343D-02  /
   data green(34,34, 1) /   2.0792198701381048D-02  /
   data green(34,34, 2) /   2.0778721352006854D-02  /
   data green(34,34, 3) /   2.0756317329756250D-02  /
   data green(34,34, 4) /   2.0725073314207622D-02  /
   data green(34,34, 5) /   2.0685109350202135D-02  /
   data green(34,34, 6) /   2.0636577574346830D-02  /
   data green(34,34, 7) /   2.0579660623612137D-02  /
   data green(34,34, 8) /   2.0514569758294736D-02  /
   data green(34,34, 9) /   2.0441542736492593D-02  /
   data green(34,34,10) /   2.0360841480922855D-02  /
   data green(34,34,11) /   2.0272749581337957D-02  /
   data green(34,34,12) /   2.0177569676950047D-02  /
   data green(34,34,13) /   2.0075620763200908D-02  /
   data green(34,34,14) /   1.9967235466002049D-02  /
   data green(34,34,15) /   1.9852757324345393D-02  /
   data green(34,34,16) /   1.9732538119103686D-02  /
   data green(34,34,17) /   1.9606935282076840D-02  /
   data green(34,34,18) /   1.9476309415077071D-02  /
   data green(34,34,19) /   1.9341021944264698D-02  /
   data green(34,34,20) /   1.9201432930219260D-02  /
   data green(34,34,21) /   1.9057899049514146D-02  /
   data green(34,34,22) /   1.8910771758994329D-02  /
   data green(34,34,23) /   1.8760395649648574D-02  /
   data green(34,34,24) /   1.8607106993007556D-02  /
   data green(34,34,25) /   1.8451232479453648D-02  /
   data green(34,34,26) /   1.8293088144734192D-02  /
   data green(34,34,27) /   1.8132978478350617D-02  /
   data green(34,34,28) /   1.7971195705349080D-02  /
   data green(34,34,29) /   1.7808019231351793D-02  /
   data green(34,34,30) /   1.7643715239415395D-02  /
   data green(34,34,31) /   1.7478536426447890D-02  /
   data green(34,34,32) /   1.7312721866417930D-02  /
   data green(34,34,33) /   1.7146496987405089D-02  /
   data green(34,34,34) /   1.6980073649621577D-02  /
   data green(35, 0, 0) /   2.8577272810746190D-02  /
   data green(35, 1, 0) /   2.8565587081569536D-02  /
   data green(35, 1, 1) /   2.8553915735151342D-02  /
   data green(35, 2, 0) /   2.8530616190853869D-02  /
   data green(35, 2, 1) /   2.8518987814806232D-02  /
   data green(35, 2, 2) /   2.8484188274959368D-02  /
   data green(35, 3, 0) /   2.8472617249973037D-02  /
   data green(35, 3, 1) /   2.8461059903362340D-02  /
   data green(35, 3, 2) /   2.8426472577692055D-02  /
   data green(35, 3, 3) /   2.8369107682204974D-02  /
   data green(35, 4, 0) /   2.8392012920022344D-02  /
   data green(35, 4, 1) /   2.8380553797963299D-02  /
   data green(35, 4, 2) /   2.8346259943223338D-02  /
   data green(35, 4, 3) /   2.8289380191456812D-02  /
   data green(35, 4, 4) /   2.8210323674718887D-02  /
   data green(35, 5, 0) /   2.8289382873411709D-02  /
   data green(35, 5, 1) /   2.8278047996166338D-02  /
   data green(35, 5, 2) /   2.8244125364679139D-02  /
   data green(35, 5, 3) /   2.8187859324199031D-02  /
   data green(35, 5, 4) /   2.8109651664400978D-02  /
   data green(35, 5, 5) /   2.8010053667120820D-02  /
   data green(35, 6, 0) /   2.8165452319699958D-02  /
   data green(35, 6, 1) /   2.8154266254300422D-02  /
   data green(35, 6, 2) /   2.8120788263305042D-02  /
   data green(35, 6, 3) /   2.8065257357659877D-02  /
   data green(35, 6, 4) /   2.7988066603371119D-02  /
   data green(35, 6, 5) /   2.7889755410566531D-02  /
   data green(35, 6, 6) /   2.7770999158450632D-02  /
   data green(35, 7, 0) /   2.8021078036990863D-02  /
   data green(35, 7, 1) /   2.8010063657796690D-02  /
   data green(35, 7, 2) /   2.7977098674757415D-02  /
   data green(35, 7, 3) /   2.7922416004331117D-02  /
   data green(35, 7, 4) /   2.7846398743760915D-02  /
   data green(35, 7, 5) /   2.7749572733216231D-02  /
   data green(35, 7, 6) /   2.7632596543889814D-02  /
   data green(35, 7, 7) /   2.7496249262434679D-02  /
   data green(35, 8, 0) /   2.7857232425905539D-02  /
   data green(35, 8, 1) /   2.7846410718666035D-02  /
   data green(35, 8, 2) /   2.7814021476175413D-02  /
   data green(35, 8, 3) /   2.7760290852001317D-02  /
   data green(35, 8, 4) /   2.7685590878919943D-02  /
   data green(35, 8, 5) /   2.7590432330746893D-02  /
   data green(35, 8, 6) /   2.7475455109205491D-02  /
   data green(35, 8, 7) /   2.7341416507371819D-02  /
   data green(35, 8, 8) /   2.7189177763633322D-02  /
   data green(35, 9, 0) /   2.7674986147901950D-02  /
   data green(35, 9, 1) /   2.7664376059955410D-02  /
   data green(35, 9, 2) /   2.7632619208098268D-02  /
   data green(35, 9, 3) /   2.7579934411649739D-02  /
   data green(35, 9, 4) /   2.7506681701251961D-02  /
   data green(35, 9, 5) /   2.7413355501672954D-02  /
   data green(35, 9, 6) /   2.7300575446056307D-02  /
   data green(35, 9, 7) /   2.7169075153209328D-02  /
   data green(35, 9, 8) /   2.7019689358867591D-02  /
   data green(35, 9, 9) /   2.6853339829174927D-02  /
   data green(35,10, 0) /   2.7475489923099304D-02  /
   data green(35,10, 1) /   2.7465108260432577D-02  /
   data green(35,10, 2) /   2.7434034058799499D-02  /
   data green(35,10, 3) /   2.7382478331462082D-02  /
   data green(35,10, 4) /   2.7310788331471950D-02  /
   data green(35,10, 5) /   2.7219441071394512D-02  /
   data green(35,10, 6) /   2.7109034586509911D-02  /
   data green(35,10, 7) /   2.6980277252414805D-02  /
   data green(35,10, 8) /   2.6833975524094531D-02  /
   data green(35,10, 9) /   2.6671020499253992D-02  /
   data green(35,10,10) /   2.6492373723356832D-02  /
   data green(35,11, 0) /   2.7259956046828522D-02  /
   data green(35,11, 1) /   2.7249817416956906D-02  /
   data green(35,11, 2) /   2.7219469563245339D-02  /
   data green(35,11, 3) /   2.7169115320557494D-02  /
   data green(35,11, 4) /   2.7099088551082718D-02  /
   data green(35,11, 5) /   2.7009848011748518D-02  /
   data green(35,11, 6) /   2.6901969081089848D-02  /
   data green(35,11, 7) /   2.6776133635472866D-02  /
   data green(35,11, 8) /   2.6633118417416422D-02  /
   data green(35,11, 9) /   2.6473782272799579D-02  /
   data green(35,11,10) /   2.6299052648352528D-02  /
   data green(35,11,11) /   2.6109911736775768D-02  /
   data green(35,12, 0) /   2.7029640142918827D-02  /
   data green(35,12, 1) /   2.7019756941016874D-02  /
   data green(35,12, 2) /   2.6990172528602155D-02  /
   data green(35,12, 3) /   2.6941081286865211D-02  /
   data green(35,12, 4) /   2.6872803232246816D-02  /
   data green(35,12, 5) /   2.6785778237409800D-02  /
   data green(35,12, 6) /   2.6680558230004172D-02  /
   data green(35,12, 7) /   2.6557797638070528D-02  /
   data green(35,12, 8) /   2.6418242400418498D-02  /
   data green(35,12, 9) /   2.6262717892620181D-02  /
   data green(35,12,10) /   2.6092116133730539D-02  /
   data green(35,12,11) /   2.5907382636127813D-02  /
   data green(35,12,12) /   2.5709503242717788D-02  /
   data green(35,13, 0) /   2.6785823610676635D-02  /
   data green(35,13, 1) /   2.6776206045156029D-02  /
   data green(35,13, 2) /   2.6747415638334736D-02  /
   data green(35,13, 3) /   2.6699638135135442D-02  /
   data green(35,13, 4) /   2.6633179402715815D-02  /
   data green(35,13, 5) /   2.6548460006365578D-02  /
   data green(35,13, 6) /   2.6446007882410756D-02  /
   data green(35,13, 7) /   2.6326449356182672D-02  /
   data green(35,13, 8) /   2.6190498799240486D-02  /
   data green(35,13, 9) /   2.6038947250548665D-02  /
   data green(35,13,10) /   2.5872650340561486D-02  /
   data green(35,13,11) /   2.5692515855670870D-02  /
   data green(35,13,12) /   2.5499491264763287D-02  /
   data green(35,13,13) /   2.5294551501999687D-02  /
   data green(35,14, 0) /   2.6529797147748894D-02  /
   data green(35,14, 1) /   2.6520453300571321D-02  /
   data green(35,14, 2) /   2.6492481113467647D-02  /
   data green(35,14, 3) /   2.6446057599182362D-02  /
   data green(35,14, 4) /   2.6381474313743828D-02  /
   data green(35,14, 5) /   2.6299132284615809D-02  /
   data green(35,14, 6) /   2.6199535154607170D-02  /
   data green(35,14, 7) /   2.6083280769325094D-02  /
   data green(35,14, 8) /   2.5951051478762654D-02  /
   data green(35,14, 9) /   2.5803603452305601D-02  /
   data green(35,14,10) /   2.5641755320389013D-02  /
   data green(35,14,11) /   2.5466376455638783D-02  /
   data green(35,14,12) /   2.5278375192907234D-02  /
   data green(35,14,13) /   2.5078687263175498D-02  /
   data green(35,14,14) /   2.4868264683323822D-02  /
   data green(35,15, 0) /   2.6262845648783104D-02  /
   data green(35,15, 1) /   2.6253781565217044D-02  /
   data green(35,15, 2) /   2.6226645728610928D-02  /
   data green(35,15, 3) /   2.6181606403090009D-02  /
   data green(35,15, 4) /   2.6118940801731837D-02  /
   data green(35,15, 5) /   2.6039030360781871D-02  /
   data green(35,15, 6) /   2.5942354346731488D-02  /
   data green(35,15, 7) /   2.5829482004492810D-02  /
   data green(35,15, 8) /   2.5701063494455566D-02  /
   data green(35,15, 9) /   2.5557819893077899D-02  /
   data green(35,15,10) /   2.5400532545227381D-02  /
   data green(35,15,11) /   2.5230032057057247D-02  /
   data green(35,15,12) /   2.5047187206889164D-02  /
   data green(35,15,13) /   2.4852894030140484D-02  /
   data green(35,15,14) /   2.4648065304958753D-02  /
   data green(35,15,15) /   2.4433620630307351D-02  /
   data green(35,16, 0) /   2.5986234695725374D-02  /
   data green(35,16, 1) /   2.5977454497986885D-02  /
   data green(35,16, 2) /   2.5951167397513544D-02  /
   data green(35,16, 3) /   2.5907532964785651D-02  /
   data green(35,16, 4) /   2.5846814155087050D-02  /
   data green(35,16, 5) /   2.5769372919613464D-02  /
   data green(35,16, 6) /   2.5675664263883343D-02  /
   data green(35,16, 7) /   2.5566228943018500D-02  /
   data green(35,16, 8) /   2.5441685019858669D-02  /
   data green(35,16, 9) /   2.5302718536924344D-02  /
   data green(35,16,10) /   2.5150073566357251D-02  /
   data green(35,16,11) /   2.4984541903358734D-02  /
   data green(35,16,12) /   2.4806952659260771D-02  /
   data green(35,16,13) /   2.4618161991716989D-02  /
   data green(35,16,14) /   2.4419043183502420D-02  /
   data green(35,16,15) /   2.4210477250180981D-02  /
   data green(35,16,16) /   2.3993344222574309D-02  /
   data green(35,17, 0) /   2.5701198774603522D-02  /
   data green(35,17, 1) /   2.5692704793825540D-02  /
   data green(35,17, 2) /   2.5667273462999497D-02  /
   data green(35,17, 3) /   2.5625055776813770D-02  /
   data green(35,17, 4) /   2.5566300621041427D-02  /
   data green(35,17, 5) /   2.5491350708928422D-02  /
   data green(35,17, 6) /   2.5400637075810167D-02  /
   data green(35,17, 7) /   2.5294672303863650D-02  /
   data green(35,17, 8) /   2.5174042682265984D-02  /
   data green(35,17, 9) /   2.5039399531296445D-02  /
   data green(35,17,10) /   2.4891449931521342D-02  /
   data green(35,17,11) /   2.4730947101279616D-02  /
   data green(35,17,12) /   2.4558680658032742D-02  /
   data green(35,17,13) /   2.4375466983052690D-02  /
   data green(35,17,14) /   2.4182139886060043D-02  /
   data green(35,17,15) /   2.3979541738650492D-02  /
   data green(35,17,16) /   2.3768515214560620D-02  /
   data green(35,17,17) /   2.3549895742823165D-02  /
   data green(35,18, 0) /   2.5408931279511198D-02  /
   data green(35,18, 1) /   2.5400724200716311D-02  /
   data green(35,18, 2) /   2.5376150752911570D-02  /
   data green(35,18, 3) /   2.5335353527021449D-02  /
   data green(35,18, 4) /   2.5278567616572403D-02  /
   data green(35,18, 5) /   2.5206116865825214D-02  /
   data green(35,18, 6) /   2.5118408782864653D-02  /
   data green(35,18, 7) /   2.5015928272980919D-02  /
   data green(35,18, 8) /   2.4899230378151747D-02  /
   data green(35,18, 9) /   2.4768932229960822D-02  /
   data green(35,18,10) /   2.4625704435327345D-02  /
   data green(35,18,11) /   2.4470262117056947D-02  /
   data green(35,18,12) /   2.4303355825104934D-02  /
   data green(35,18,13) /   2.4125762520674915D-02  /
   data green(35,18,14) /   2.3938276815299062D-02  /
   data green(35,18,15) /   2.3741702622481164D-02  /
   data green(35,18,16) /   2.3536845352007556D-02  /
   data green(35,18,17) /   2.3324504748248553D-02  /
   data green(35,18,18) /   2.3105468445123335D-02  /
   data green(35,19, 0) /   2.5110576299874513D-02  /
   data green(35,19, 1) /   2.5102655315029858D-02  /
   data green(35,19, 2) /   2.5078937399743206D-02  /
   data green(35,19, 3) /   2.5039556958776438D-02  /
   data green(35,19, 4) /   2.4984735645656203D-02  /
   data green(35,19, 5) /   2.4914778907578310D-02  /
   data green(35,19, 6) /   2.4830071297275047D-02  /
   data green(35,19, 7) /   2.4731070691742640D-02  /
   data green(35,19, 8) /   2.4618301585468814D-02  /
   data green(35,19, 9) /   2.4492347645645819D-02  /
   data green(35,19,10) /   2.4353843728292725D-02  /
   data green(35,19,11) /   2.4203467557277793D-02  /
   data green(35,19,12) /   2.4041931263458299D-02  /
   data green(35,19,13) /   2.3869972969474160D-02  /
   data green(35,19,14) /   2.3688348588381414D-02  /
   data green(35,19,15) /   2.3497823982699947D-02  /
   data green(35,19,16) /   2.3299167606048447D-02  /
   data green(35,19,17) /   2.3093143723769052D-02  /
   data green(35,19,18) /   2.2880506283096746D-02  /
   data green(35,19,19) /   2.2661993478604123D-02  /
   data green(35,20, 0) /   2.4807222133460188D-02  /
   data green(35,20, 1) /   2.4799585098218523D-02  /
   data green(35,20, 2) /   2.4776716368501144D-02  /
   data green(35,20, 3) /   2.4738742417800944D-02  /
   data green(35,20, 4) /   2.4685871873378534D-02  /
   data green(35,20, 5) /   2.4618392341974195D-02  /
   data green(35,20, 6) /   2.4536666099382484D-02  /
   data green(35,20, 7) /   2.4441124769504326D-02  /
   data green(35,20, 8) /   2.4332263143669458D-02  /
   data green(35,20, 9) /   2.4210632309253587D-02  /
   data green(35,20,10) /   2.4076832267426049D-02  /
   data green(35,20,11) /   2.3931504223250136D-02  /
   data green(35,20,12) /   2.3775322727747444D-02  /
   data green(35,20,13) /   2.3608987841714828D-02  /
   data green(35,20,14) /   2.3433217476104145D-02  /
   data green(35,20,15) /   2.3248740044858340D-02  /
   data green(35,20,16) /   2.3056287544526311D-02  /
   data green(35,20,17) /   2.2856589152007875D-02  /
   data green(35,20,18) /   2.2650365408560098D-02  /
   data green(35,20,19) /   2.2438323035719725D-02  /
   data green(35,20,20) /   2.2221150407864917D-02  /
   data green(35,21, 0) /   2.4499896425517205D-02  /
   data green(35,21, 1) /   2.4492540015845610D-02  /
   data green(35,21, 2) /   2.4470510595564078D-02  /
   data green(35,21, 3) /   2.4433926991289014D-02  /
   data green(35,21, 4) /   2.4382985266519247D-02  /
   data green(35,21, 5) /   2.4317955811590977D-02  /
   data green(35,21, 6) /   2.4239179389039130D-02  /
   data green(35,21, 7) /   2.4147062246845367D-02  /
   data green(35,21, 8) /   2.4042070434839156D-02  /
   data green(35,21, 9) /   2.3924723476213366D-02  /
   data green(35,21,10) /   2.3795587556283339D-02  /
   data green(35,21,11) /   2.3655268394223505D-02  /
   data green(35,21,12) /   2.3504403960897317D-02  /
   data green(35,21,13) /   2.3343657197713284D-02  /
   data green(35,21,14) /   2.3173708878584696D-02  /
   data green(35,21,15) /   2.2995250740594685D-02  /
   data green(35,21,16) /   2.2808978989984165D-02  /
   data green(35,21,17) /   2.2615588269690221D-02  /
   data green(35,21,18) /   2.2415766153887073D-02  /
   data green(35,21,19) /   2.2210188214714151D-02  /
   data green(35,21,20) /   2.1999513687345804D-02  /
   data green(35,21,21) /   2.1784381742319722D-02  /
   data green(35,22, 0) /   2.4189562803672827D-02  /
   data green(35,22, 1) /   2.4182482669198511D-02  /
   data green(35,22, 2) /   2.4161279610653140D-02  /
   data green(35,22, 3) /   2.4126065114483404D-02  /
   data green(35,22, 4) /   2.4077023179960039D-02  /
   data green(35,22, 5) /   2.4014407656554122D-02  /
   data green(35,22, 6) /   2.3938538622769331D-02  /
   data green(35,22, 7) /   2.3849797906900939D-02  /
   data green(35,22, 8) /   2.3748623870757105D-02  /
   data green(35,22, 9) /   2.3635505592620211D-02  /
   data green(35,22,10) /   2.3510976595244402D-02  /
   data green(35,22,11) /   2.3375608268419487D-02  /
   data green(35,22,12) /   2.3230003133851967D-02  /
   data green(35,22,13) /   2.3074788093361290D-02  /
   data green(35,22,14) /   2.2910607790424080D-02  /
   data green(35,22,15) /   2.2738118200815706D-02  /
   data green(35,22,16) /   2.2557980551462579D-02  /
   data green(35,22,17) /   2.2370855648592470D-02  /
   data green(35,22,18) /   2.2177398677753646D-02  /
   data green(35,22,19) /   2.1978254520068140D-02  /
   data green(35,22,20) /   2.1774053611850423D-02  /
   data green(35,22,21) /   2.1565408358975399D-02  /
   data green(35,22,22) /   2.1352910103481112D-02  /
   data green(35,23, 0) /   2.3877118857882174D-02  /
   data green(35,23, 1) /   2.3870309769408200D-02  /
   data green(35,23, 2) /   2.3849917493693498D-02  /
   data green(35,23, 3) /   2.3816046499535313D-02  /
   data green(35,23, 4) /   2.3768869247903897D-02  /
   data green(35,23, 5) /   2.3708623759947586D-02  /
   data green(35,23, 6) /   2.3635610306973309D-02  /
   data green(35,23, 7) /   2.3550187311945827D-02  /
   data green(35,23, 8) /   2.3452766570563357D-02  /
   data green(35,23, 9) /   2.3343807913845306D-02  /
   data green(35,23,10) /   2.3223813443038178D-02  /
   data green(35,23,11) /   2.3093321471435231D-02  /
   data green(35,23,12) /   2.2952900306620828D-02  /
   data green(35,23,13) /   2.2803142001140191D-02  /
   data green(35,23,14) /   2.2644656190296156D-02  /
   data green(35,23,15) /   2.2478064123448750D-02  /
   data green(35,23,16) /   2.2303992980675119D-02  /
   data green(35,23,17) /   2.2123070550770986D-02  /
   data green(35,23,18) /   2.1935920330133644D-02  /
   data green(35,23,19) /   2.1743157085772444D-02  /
   data green(35,23,20) /   2.1545382910143008D-02  /
   data green(35,23,21) /   2.1343183781170547D-02  /
   data green(35,23,22) /   2.1137126628057046D-02  /
   data green(35,23,23) /   2.0927756892474680D-02  /
   data green(35,24, 0) /   2.3563395303617420D-02  /
   data green(35,24, 1) /   2.3556851292853651D-02  /
   data green(35,24, 2) /   2.3537252007109735D-02  /
   data green(35,24, 3) /   2.3504695230091469D-02  /
   data green(35,24, 4) /   2.3459342427325721D-02  /
   data green(35,24, 5) /   2.3401416528268579D-02  /
   data green(35,24, 6) /   2.3331198905424580D-02  /
   data green(35,24, 7) /   2.3249025630120809D-02  /
   data green(35,24, 8) /   2.3155283101213588D-02  /
   data green(35,24, 9) /   2.3050403155607160D-02  /
   data green(35,24,10) /   2.2934857777697222D-02  /
   data green(35,24,11) /   2.2809153528633806D-02  /
   data green(35,24,12) /   2.2673825815782665D-02  /
   data green(35,24,13) /   2.2529433118321871D-02  /
   data green(35,24,14) /   2.2376551277070705D-02  /
   data green(35,24,15) /   2.2215767946060550D-02  /
   data green(35,24,16) /   2.2047677290735673D-02  /
   data green(35,24,17) /   2.1872875003739918D-02  /
   data green(35,24,18) /   2.1691953694697973D-02  /
   data green(35,24,19) /   2.1505498695867337D-02  /
   data green(35,24,20) /   2.1314084311556996D-02  /
   data green(35,24,21) /   2.1118270526213744D-02  /
   data green(35,24,22) /   2.0918600174389564D-02  /
   data green(35,24,23) /   2.0715596565635413D-02  /
   data green(35,24,24) /   2.0509761548831249D-02  /
   data green(35,25, 0) /   2.3249156163120541D-02  /
   data green(35,25, 1) /   2.3242870653218568D-02  /
   data green(35,25, 2) /   2.3224044740541593D-02  /
   data green(35,25, 3) /   2.3192769861257394D-02  /
   data green(35,25, 4) /   2.3149197036960834D-02  /
   data green(35,25, 5) /   2.3093534854816546D-02  /
   data green(35,25, 6) /   2.3026046714376830D-02  /
   data green(35,25, 7) /   2.2947047411788622D-02  /
   data green(35,25, 8) /   2.2856899147017225D-02  /
   data green(35,25, 9) /   2.2756007051129980D-02  /
   data green(35,25,10) /   2.2644814338297441D-02  /
   data green(35,25,11) /   2.2523797190889872D-02  /
   data green(35,25,12) /   2.2393459485992612D-02  /
   data green(35,25,13) /   2.2254327468129418D-02  /
   data green(35,25,14) /   2.2106944466415798D-02  /
   data green(35,25,15) /   2.1951865745309194D-02  /
   data green(35,25,16) /   2.1789653567190981D-02  /
   data green(35,25,17) /   2.1620872532832531D-02  /
   data green(35,25,18) /   2.1446085252967678D-02  /
   data green(35,25,19) /   2.1265848391276428D-02  /
   data green(35,25,20) /   2.1080709106556930D-02  /
   data green(35,25,21) /   2.0891201910120696D-02  /
   data green(35,25,22) /   2.0697845943788856D-02  /
   data green(35,25,23) /   2.0501142674500334D-02  /
   data green(35,25,24) /   2.0301573993584909D-02  /
   data green(35,25,25) /   2.0099600702246569D-02  /
   data green(35,26, 0) /   2.2935099802411981D-02  /
   data green(35,26, 1) /   2.2929065728373169D-02  /
   data green(35,26, 2) /   2.2910992107588513D-02  /
   data green(35,26, 3) /   2.2880964366912098D-02  /
   data green(35,26, 4) /   2.2839123637064074D-02  /
   data green(35,26, 5) /   2.2785664915342327D-02  /
   data green(35,26, 6) /   2.2720834559468800D-02  /
   data green(35,26, 7) /   2.2644927175259094D-02  /
   data green(35,26, 8) /   2.2558281974149465D-02  /
   data green(35,26, 9) /   2.2461278686935374D-02  /
   data green(35,26,10) /   2.2354333127091484D-02  /
   data green(35,26,11) /   2.2237892500663044D-02  /
   data green(35,26,12) /   2.2112430560026670D-02  /
   data green(35,26,13) /   2.1978442696054014D-02  /
   data green(35,26,14) /   2.1836441057744387D-02  /
   data green(35,26,15) /   2.1686949780681673D-02  /
   data green(35,26,16) /   2.1530500396237068D-02  /
   data green(35,26,17) /   2.1367627482819784D-02  /
   data green(35,26,18) /   2.1198864609200230D-02  /
   data green(35,26,19) /   2.1024740608482036D-02  /
   data green(35,26,20) /   2.0845776210109988D-02  /
   data green(35,26,21) /   2.0662481046726022D-02  /
   data green(35,26,22) /   2.0475351043001314D-02  /
   data green(35,26,23) /   2.0284866184977938D-02  /
   data green(35,26,24) /   2.0091488661066408D-02  /
   data green(35,26,25) /   1.9895661359721328D-02  /
   data green(35,26,26) /   1.9697806703948578D-02  /
   data green(35,27, 0) /   2.2621860669360231D-02  /
   data green(35,27, 1) /   2.2616070587801535D-02  /
   data green(35,27, 2) /   2.2598727041540277D-02  /
   data green(35,27, 3) /   2.2569909783371728D-02  /
   data green(35,27, 4) /   2.2529750601750199D-02  /
   data green(35,27, 5) /   2.2478431651313996D-02  /
   data green(35,27, 6) /   2.2416183174009192D-02  /
   data green(35,27, 7) /   2.2343280666309986D-02  /
   data green(35,27, 8) /   2.2260041559958470D-02  /
   data green(35,27, 9) /   2.2166821492956009D-02  /
   data green(35,27,10) /   2.2064010253983132D-02  /
   data green(35,27,11) /   2.1952027486912488D-02  /
   data green(35,27,12) /   2.1831318242667267D-02  /
   data green(35,27,13) /   2.1702348463559361D-02  /
   data green(35,27,14) /   2.1565600480719104D-02  /
   data green(35,27,15) /   2.1421568598691577D-02  /
   data green(35,27,16) /   2.1270754833160849D-02  /
   data green(35,27,17) /   2.1113664858536715D-02  /
   data green(35,27,18) /   2.0950804212253756D-02  /
   data green(35,27,19) /   2.0782674792517934D-02  /
   data green(35,27,20) /   2.0609771676271108D-02  /
   data green(35,27,21) /   2.0432580274650249D-02  /
   data green(35,27,22) /   2.0251573834448131D-02  /
   data green(35,27,23) /   2.0067211286224490D-02  /
   data green(35,27,24) /   1.9879935432889153D-02  /
   data green(35,27,25) /   1.9690171466843463D-02  /
   data green(35,27,26) /   1.9498325799132629D-02  /
   data green(35,27,27) /   1.9304785180493115D-02  /
   data green(35,28, 0) /   2.2310011589119907D-02  /
   data green(35,28, 1) /   2.2304457777226564D-02  /
   data green(35,28, 2) /   2.2287821247777898D-02  /
   data green(35,28, 3) /   2.2260176408793280D-02  /
   data green(35,28, 4) /   2.2221646245666426D-02  /
   data green(35,28, 5) /   2.2172400805532944D-02  /
   data green(35,28, 6) /   2.2112655126942930D-02  /
   data green(35,28, 7) /   2.2042666663912917D-02  /
   data green(35,28, 8) /   2.1962732264072266D-02  /
   data green(35,28, 9) /   2.1873184769002092D-02  /
   data green(35,28,10) /   2.1774389310766917D-02  /
   data green(35,28,11) /   2.1666739381971795D-02  /
   data green(35,28,12) /   2.1550652757475516D-02  /
   data green(35,28,13) /   2.1426567344309189D-02  /
   data green(35,28,14) /   2.1294937032636227D-02  /
   data green(35,28,15) /   2.1156227615071129D-02  /
   data green(35,28,16) /   2.1010912834720390D-02  /
   data green(35,28,17) /   2.0859470614315340D-02  /
   data green(35,28,18) /   2.0702379510167502D-02  /
   data green(35,28,19) /   2.0540115425765702D-02  /
   data green(35,28,20) /   2.0373148610984902D-02  /
   data green(35,28,21) /   2.0201940964378926D-02  /
   data green(35,28,22) /   2.0026943648112396D-02  /
   data green(35,28,23) /   1.9848595017928190D-02  /
   data green(35,28,24) /   1.9667318864261023D-02  /
   data green(35,28,25) /   1.9483522955263108D-02  /
   data green(35,28,26) /   1.9297597868122605D-02  /
   data green(35,28,27) /   1.9109916091612170D-02  /
   data green(35,28,28) /   1.8920831380252769D-02  /
   data green(35,29, 0) /   2.2000066486442097D-02  /
   data green(35,29, 1) /   2.1994741030216682D-02  /
   data green(35,29, 2) /   2.1978787883465259D-02  /
   data green(35,29, 3) /   2.1952276430322236D-02  /
   data green(35,29, 4) /   2.1915321378918832D-02  /
   data green(35,29, 5) /   2.1868081386453596D-02  /
   data green(35,29, 6) /   2.1810757179773665D-02  /
   data green(35,29, 7) /   2.1743589214818712D-02  /
   data green(35,29, 8) /   2.1666854927762151D-02  /
   data green(35,29, 9) /   2.1580865638220040D-02  /
   data green(35,29,10) /   2.1485963170290320D-02  /
   data green(35,29,11) /   2.1382516260344833D-02  /
   data green(35,29,12) /   2.1270916821447081D-02  /
   data green(35,29,13) /   2.1151576133129069D-02  /
   data green(35,29,14) /   2.1024921022236080D-02  /
   data green(35,29,15) /   2.0891390095908815D-02  /
   data green(35,29,16) /   2.0751430081832046D-02  /
   data green(35,29,17) /   2.0605492323973219D-02  /
   data green(35,29,18) /   2.0454029474504500D-02  /
   data green(35,29,19) /   2.0297492414770218D-02  /
   data green(35,29,20) /   2.0136327430325156D-02  /
   data green(35,29,21) /   1.9970973657482880D-02  /
   data green(35,29,22) /   1.9801860811689123D-02  /
   data green(35,29,23) /   1.9629407201532900D-02  /
   data green(35,29,24) /   1.9454018026444135D-02  /
   data green(35,29,25) /   1.9276083951166579D-02  /
   data green(35,29,26) /   1.9095979945968806D-02  /
   data green(35,29,27) /   1.8914064378257083D-02  /
   data green(35,29,28) /   1.8730678338745775D-02  /
   data green(35,29,29) /   1.8546145183570083D-02  /
   data green(35,30, 0) /   2.1692483418753634D-02  /
   data green(35,30, 1) /   2.1687378290888539D-02  /
   data green(35,30, 2) /   2.1672084549298351D-02  /
   data green(35,30, 3) /   2.1646666864834428D-02  /
   data green(35,30, 4) /   2.1611232177608828D-02  /
   data green(35,30, 5) /   2.1565928450481105D-02  /
   data green(35,30, 6) /   2.1510942964031242D-02  /
   data green(35,30, 7) /   2.1446500191284162D-02  /
   data green(35,30, 8) /   2.1372859298892672D-02  /
   data green(35,30, 9) /   2.1290311328250350D-02  /
   data green(35,30,10) /   2.1199176114919845D-02  /
   data green(35,30,11) /   2.1099799007739974D-02  /
   data green(35,30,12) /   2.0992547450029014D-02  /
   data green(35,30,13) /   2.0877807484522684D-02  /
   data green(35,30,14) /   2.0755980241243605D-02  /
   data green(35,30,15) /   2.0627478463615880D-02  /
   data green(35,30,16) /   2.0492723123080044D-02  /
   data green(35,30,17) /   2.0352140166514585D-02  /
   data green(35,30,18) /   2.0206157434223070D-02  /
   data green(35,30,19) /   2.0055201779380888D-02  /
   data green(35,30,20) /   1.9899696412912941D-02  /
   data green(35,30,21) /   1.9740058491019570D-02  /
   data green(35,30,22) /   1.9576696956173978D-02  /
   data green(35,30,23) /   1.9410010636528554D-02  /
   data green(35,30,24) /   1.9240386603398517D-02  /
   data green(35,30,25) /   1.9068198781909112D-02  /
   data green(35,30,26) /   1.8893806806029242D-02  /
   data green(35,30,27) /   1.8717555106074958D-02  /
   data green(35,30,28) /   1.8539772214326258D-02  /
   data green(35,30,29) /   1.8360770272618968D-02  /
   data green(35,30,30) /   1.8180844724592419D-02  /
   data green(35,31, 0) /   2.1387667818683440D-02  /
   data green(35,31, 1) /   2.1382774946547230D-02  /
   data green(35,31, 2) /   2.1368116492642541D-02  /
   data green(35,31, 3) /   2.1343752713415733D-02  /
   data green(35,31, 4) /   2.1309783271257807D-02  /
   data green(35,31, 5) /   2.1266346104977112D-02  /
   data green(35,31, 6) /   2.1213615883778115D-02  /
   data green(35,31, 7) /   2.1151802078499888D-02  /
   data green(35,31, 8) /   2.1081146691374897D-02  /
   data green(35,31, 9) /   2.1001921691633695D-02  /
   data green(35,31,10) /   2.0914426208749174D-02  /
   data green(35,31,11) /   2.0818983537905546D-02  /
   data green(35,31,12) /   2.0715938013394386D-02  /
   data green(35,31,13) /   2.0605651805153589D-02  /
   data green(35,31,14) /   2.0488501691712693D-02  /
   data green(35,31,15) /   2.0364875859573261D-02  /
   data green(35,31,16) /   2.0235170774759066D-02  /
   data green(35,31,17) /   2.0099788167159134D-02  /
   data green(35,31,18) /   1.9959132162607782D-02  /
   data green(35,31,19) /   1.9813606591642212D-02  /
   data green(35,31,20) /   1.9663612497776900D-02  /
   data green(35,31,21) /   1.9509545862135769D-02  /
   data green(35,31,22) /   1.9351795555559798D-02  /
   data green(35,31,23) /   1.9190741523995025D-02  /
   data green(35,31,24) /   1.9026753208165773D-02  /
   data green(35,31,25) /   1.8860188194319998D-02  /
   data green(35,31,26) /   1.8691391089234278D-02  /
   data green(35,31,27) /   1.8520692609696066D-02  /
   data green(35,31,28) /   1.8348408874327191D-02  /
   data green(35,31,29) /   1.8174840883844672D-02  /
   data green(35,31,30) /   1.8000274174627552D-02  /
   data green(35,31,31) /   1.7824978629717287D-02  /
   data green(35,32, 0) /   2.1085975859265364D-02  /
   data green(35,32, 1) /   2.1081287183608566D-02  /
   data green(35,32, 2) /   2.1067239935747894D-02  /
   data green(35,32, 3) /   2.1043890243846636D-02  /
   data green(35,32, 4) /   2.1011330962197414D-02  /
   data green(35,32, 5) /   2.0969690648094844D-02  /
   data green(35,32, 6) /   2.0919132160553498D-02  /
   data green(35,32, 7) /   2.0859850910629488D-02  /
   data green(35,32, 8) /   2.0792072799774568D-02  /
   data green(35,32, 9) /   2.0716051888083679D-02  /
   data green(35,32,10) /   2.0632067838351427D-02  /
   data green(35,32,11) /   2.0540423184456549D-02  /
   data green(35,32,12) /   2.0441440473742561D-02  /
   data green(35,32,13) /   2.0335459332809905D-02  /
   data green(35,32,14) /   2.0222833504591538D-02  /
   data green(35,32,15) /   2.0103927901902234D-02  /
   data green(35,32,16) /   1.9979115719018576D-02  /
   data green(35,32,17) /   1.9848775638465276D-02  /
   data green(35,32,18) /   1.9713289165268991D-02  /
   data green(35,32,19) /   1.9573038115701710D-02  /
   data green(35,32,20) /   1.9428402282170687D-02  /
   data green(35,32,21) /   1.9279757290597090D-02  /
   data green(35,32,22) /   1.9127472661514929D-02  /
   data green(35,32,23) /   1.8971910081339714D-02  /
   data green(35,32,24) /   1.8813421885896627D-02  /
   data green(35,32,25) /   1.8652349754429701D-02  /
   data green(35,32,26) /   1.8489023608973944D-02  /
   data green(35,32,27) /   1.8323760711178635D-02  /
   data green(35,32,28) /   1.8156864946417078D-02  /
   data green(35,32,29) /   1.7988626283283324D-02  /
   data green(35,32,30) /   1.7819320395326660D-02  /
   data green(35,32,31) /   1.7649208431064725D-02  /
   data green(35,32,32) /   1.7478536917897330D-02  /
   data green(35,33, 0) /   2.0787717868914507D-02  /
   data green(35,33, 1) /   2.0783225393974487D-02  /
   data green(35,33, 2) /   2.0769765456438084D-02  /
   data green(35,33, 3) /   2.0747390328863816D-02  /
   data green(35,33, 4) /   2.0716186505230074D-02  /
   data green(35,33, 5) /   2.0676273774442734D-02  /
   data green(35,33, 6) /   2.0627803950616384D-02  /
   data green(35,33, 7) /   2.0570959286354509D-02  /
   data green(35,33, 8) /   2.0505950601178751D-02  /
   data green(35,33, 9) /   2.0433015162116624D-02  /
   data green(35,33,10) /   2.0352414357130012D-02  /
   data green(35,33,11) /   2.0264431204485858D-02  /
   data green(35,33,12) /   2.0169367742324900D-02  /
   data green(35,33,13) /   2.0067542342616564D-02  /
   data green(35,33,14) /   1.9959286992484074D-02  /
   data green(35,33,15) /   1.9844944583673623D-02  /
   data green(35,33,16) /   1.9724866247874947D-02  /
   data green(35,33,17) /   1.9599408771855991D-02  /
   data green(35,33,18) /   1.9468932122129835D-02  /
   data green(35,33,19) /   1.9333797104309897D-02  /
   data green(35,33,20) /   1.9194363177601172D-02  /
   data green(35,33,21) /   1.9050986440176342D-02  /
   data green(35,33,22) /   1.8904017796632758D-02  /
   data green(35,33,23) /   1.8753801314432491D-02  /
   data green(35,33,24) /   1.8600672772279763D-02  /
   data green(35,33,25) /   1.8444958399853616D-02  /
   data green(35,33,26) /   1.8286973805228309D-02  /
   data green(35,33,27) /   1.8127023083698807D-02  /
   data green(35,33,28) /   1.7965398099586496D-02  /
   data green(35,33,29) /   1.7802377930914867D-02  /
   data green(35,33,30) /   1.7638228465592716D-02  /
   data green(35,33,31) /   1.7473202136886166D-02  /
   data green(35,33,32) /   1.7307537785461423D-02  /
   data green(35,33,33) /   1.7141460635092014D-02  /
   data green(35,34, 0) /   2.0493161736146136D-02  /
   data green(35,34, 1) /   2.0488857571872039D-02  /
   data green(35,34, 2) /   2.0475961361407034D-02  /
   data green(35,34, 3) /   2.0454521780543654D-02  /
   data green(35,34, 4) /   2.0424619388208769D-02  /
   data green(35,34, 5) /   2.0386365787629885D-02  /
   data green(35,34, 6) /   2.0339902476053735D-02  /
   data green(35,34, 7) /   2.0285399406123440D-02  /
   data green(35,34, 8) /   2.0223053287279125D-02  /
   data green(35,34, 9) /   2.0153085659891067D-02  /
   data green(35,34,10) /   2.0075740778156246D-02  /
   data green(35,34,11) /   1.9991283340027897D-02  /
   data green(35,34,12) /   1.9899996103588998D-02  /
   data green(35,34,13) /   1.9802177429355076D-02  /
   data green(35,34,14) /   1.9698138787070826D-02  /
   data green(35,34,15) /   1.9588202263750236D-02  /
   data green(35,34,16) /   1.9472698107133247D-02  /
   data green(35,34,17) /   1.9351962335537289D-02  /
   data green(35,34,18) /   1.9226334441423926D-02  /
   data green(35,34,19) /   1.9096155212036635D-02  /
   data green(35,34,20) /   1.8961764686340841D-02  /
   data green(35,34,21) /   1.8823500263352076D-02  /
   data green(35,34,22) /   1.8681694972890826D-02  /
   data green(35,34,23) /   1.8536675915955966D-02  /
   data green(35,34,24) /   1.8388762878343574D-02  /
   data green(35,34,25) /   1.8238267117915354D-02  /
   data green(35,34,26) /   1.8085490323079434D-02  /
   data green(35,34,27) /   1.7930723737611447D-02  /
   data green(35,34,28) /   1.7774247444917509D-02  /
   data green(35,34,29) /   1.7616329803218013D-02  /
   data green(35,34,30) /   1.7457227021892476D-02  /
   data green(35,34,31) /   1.7297182868342541D-02  /
   data green(35,34,32) /   1.7136428494171203D-02  /
   data green(35,34,33) /   1.6975182369202185D-02  /
   data green(35,34,34) /   1.6813650311838079D-02  /
   data green(35,35, 0) /   2.0202536255696596D-02  /
   data green(35,35, 1) /   2.0198412652828776D-02  /
   data green(35,35, 2) /   2.0186057003839228D-02  /
   data green(35,35, 3) /   2.0165514632600596D-02  /
   data green(35,35, 4) /   2.0136860565446880D-02  /
   data green(35,35, 5) /   2.0100198771772500D-02  /
   data green(35,35, 6) /   2.0055661122063990D-02  /
   data green(35,35, 7) /   2.0003406083719238D-02  /
   data green(35,35, 8) /   1.9943617179675625D-02  /
   data green(35,35, 9) /   1.9876501238748245D-02  /
   data green(35,35,10) /   1.9802286469580273D-02  /
   data green(35,35,11) /   1.9721220392172351D-02  /
   data green(35,35,12) /   1.9633567662070840D-02  /
   data green(35,35,13) /   1.9539607822478121D-02  /
   data green(35,35,14) /   1.9439633018858669D-02  /
   data green(35,35,15) /   1.9333945709135207D-02  /
   data green(35,35,16) /   1.9222856400409643D-02  /
   data green(35,35,17) /   1.9106681440425034D-02  /
   data green(35,35,18) /   1.8985740888839014D-02  /
   data green(35,35,19) /   1.8860356489939379D-02  /
   data green(35,35,20) /   1.8730849764826665D-02  /
   data green(35,35,21) /   1.8597540237437590D-02  /
   data green(35,35,22) /   1.8460743805192855D-02  /
   data green(35,35,23) /   1.8320771261614622D-02  /
   data green(35,35,24) /   1.8177926975047332D-02  /
   data green(35,35,25) /   1.8032507724686817D-02  /
   data green(35,35,26) /   1.7884801692516362D-02  /
   data green(35,35,27) /   1.7735087607488851D-02  /
   data green(35,35,28) /   1.7583634036389184D-02  /
   data green(35,35,29) /   1.7430698814260345D-02  /
   data green(35,35,30) /   1.7276528606063546D-02  /
   data green(35,35,31) /   1.7121358590350393D-02  /
   data green(35,35,32) /   1.6965412255122925D-02  /
   data green(35,35,33) /   1.6808901295717531D-02  /
   data green(35,35,34) /   1.6652025604437490D-02  /
   data green(35,35,35) /   1.6494973341742009D-02  /
   data green(36, 0, 0) /   2.7783147668698868D-02  /
   data green(36, 1, 0) /   2.7772410160632893D-02  /
   data green(36, 1, 1) /   2.7761685140711419D-02  /
   data green(36, 2, 0) /   2.7740272565505816D-02  /
   data green(36, 2, 1) /   2.7729584864092199D-02  /
   data green(36, 2, 2) /   2.7697596107621875D-02  /
   data green(36, 3, 0) /   2.7686958211259029D-02  /
   data green(36, 3, 1) /   2.7676332224943524D-02  /
   data green(36, 3, 2) /   2.7644527896065377D-02  /
   data green(36, 3, 3) /   2.7591764691084836D-02  /
   data green(36, 4, 0) /   2.7612834502946503D-02  /
   data green(36, 4, 1) /   2.7602293918511222D-02  /
   data green(36, 4, 2) /   2.7570744806833841D-02  /
   data green(36, 4, 3) /   2.7518403695912304D-02  /
   data green(36, 4, 4) /   2.7445626862420528D-02  /
   data green(36, 5, 0) /   2.7518405903070354D-02  /
   data green(36, 5, 1) /   2.7507973440492797D-02  /
   data green(36, 5, 2) /   2.7476747450867509D-02  /
   data green(36, 5, 3) /   2.7424940765200384D-02  /
   data green(36, 5, 4) /   2.7352903608771902D-02  /
   data green(36, 5, 5) /   2.7261117037965465D-02  /
   data green(36, 6, 0) /   2.7404304476695226D-02  /
   data green(36, 6, 1) /   2.7394001657496975D-02  /
   data green(36, 6, 2) /   2.7363163118602783D-02  /
   data green(36, 6, 3) /   2.7311997292419576D-02  /
   data green(36, 6, 4) /   2.7240847203376866D-02  /
   data green(36, 6, 5) /   2.7150184093389609D-02  /
   data green(36, 6, 6) /   2.7040598826412937D-02  /
   data green(36, 7, 0) /   2.7271278345878384D-02  /
   data green(36, 7, 1) /   2.7261125291938002D-02  /
   data green(36, 7, 2) /   2.7230734355635285D-02  /
   data green(36, 7, 3) /   2.7180308934748470D-02  /
   data green(36, 7, 4) /   2.7110183810226923D-02  /
   data green(36, 7, 5) /   2.7020818985608767D-02  /
   data green(36, 7, 6) /   2.6912791376894261D-02  /
   data green(36, 7, 7) /   2.6786784644738781D-02  /
   data green(36, 8, 0) /   2.7120178459213544D-02  /
   data green(36, 8, 1) /   2.7110193725795172D-02  /
   data green(36, 8, 2) /   2.7080305868439236D-02  /
   data green(36, 8, 3) /   2.7030712686466635D-02  /
   data green(36, 8, 4) /   2.6961739792382283D-02  /
   data green(36, 8, 5) /   2.6873834686723830D-02  /
   data green(36, 8, 6) /   2.6767558762064596D-02  /
   data green(36, 8, 7) /   2.6643577513934976D-02  /
   data green(36, 8, 8) /   2.6502649287298924D-02  /
   data green(36, 9, 0) /   2.6951944120449631D-02  /
   data green(36, 9, 1) /   2.6942144565743280D-02  /
   data green(36, 9, 2) /   2.6912810198320918D-02  /
   data green(36, 9, 3) /   2.6864132731929818D-02  /
   data green(36, 9, 4) /   2.6796427810741166D-02  /
   data green(36, 9, 5) /   2.6710129337100319D-02  /
   data green(36, 9, 6) /   2.6605781813175730D-02  /
   data green(36, 9, 7) /   2.6484030959292221D-02  /
   data green(36, 9, 8) /   2.6345612920212771D-02  /
   data green(36, 9, 9) /   2.6191342402368145D-02  /
   data green(36,10, 0) /   2.6767587734613282D-02  /
   data green(36,10, 1) /   2.6757988426636469D-02  /
   data green(36,10, 2) /   2.6729252616278384D-02  /
   data green(36,10, 3) /   2.6681565523727142D-02  /
   data green(36,10, 4) /   2.6615232154015657D-02  /
   data green(36,10, 5) /   2.6530671891105816D-02  /
   data green(36,10, 6) /   2.6428411195352313D-02  /
   data green(36,10, 7) /   2.6309074651553948D-02  /
   data green(36,10, 8) /   2.6173374660759273D-02  /
   data green(36,10, 9) /   2.6022100099383862D-02  /
   data green(36,10,10) /   2.5856104283378027D-02  /
   data green(36,11, 0) /   2.6568179222042472D-02  /
   data green(36,11, 1) /   2.6558793382394048D-02  /
   data green(36,11, 2) /   2.6530695683773982D-02  /
   data green(36,11, 3) /   2.6484064524387067D-02  /
   data green(36,11, 4) /   2.6419193729187503D-02  /
   data green(36,11, 5) /   2.6336487419773246D-02  /
   data green(36,11, 6) /   2.6236453080589001D-02  /
   data green(36,11, 7) /   2.6119693052730542D-02  /
   data green(36,11, 8) /   2.5986894730000040D-02  /
   data green(36,11, 9) /   2.5838819760833368D-02  /
   data green(36,11,10) /   2.5676292573676238D-02  /
   data green(36,11,11) /   2.5500188542754033D-02  /
   data green(36,12, 0) /   2.6354830522924151D-02  /
   data green(36,12, 1) /   2.6345669505682592D-02  /
   data green(36,12, 2) /   2.6318243897239371D-02  /
   data green(36,12, 3) /   2.6272725023576176D-02  /
   data green(36,12, 4) /   2.6209395116327264D-02  /
   data green(36,12, 5) /   2.6128642464188705D-02  /
   data green(36,12, 6) /   2.6030954855773641D-02  /
   data green(36,12, 7) /   2.5916911529199543D-02  /
   data green(36,12, 8) /   2.5787173884400372D-02  /
   data green(36,12, 9) /   2.5642475241654611D-02  /
   data green(36,12,10) /   2.5483609943490665D-02  /
   data green(36,12,11) /   2.5311422097316599D-02  /
   data green(36,12,12) /   2.5126794244014365D-02  /
   data green(36,13, 0) /   2.6128680571111122D-02  /
   data green(36,13, 1) /   2.6119753874176015D-02  /
   data green(36,13, 2) /   2.6093028791134672D-02  /
   data green(36,13, 3) /   2.6048669400737912D-02  /
   data green(36,13, 4) /   2.5986946051024411D-02  /
   data green(36,13, 5) /   2.5908230794443921D-02  /
   data green(36,13, 6) /   2.5812991210610172D-02  /
   data green(36,13, 7) /   2.5701782816090230D-02  /
   data green(36,13, 8) /   2.5575240298682189D-02  /
   data green(36,13, 9) /   2.5434067839620100D-02  /
   data green(36,13,10) /   2.5279028800470421D-02  /
   data green(36,13,11) /   2.5110935052416131D-02  /
   data green(36,13,12) /   2.4930636215196925D-02  /
   data green(36,13,13) /   2.4739009052838080D-02  /
   data green(36,14, 0) /   2.5890881060570956D-02  /
   data green(36,14, 1) /   2.5882196365987778D-02  /
   data green(36,14, 2) /   2.5856194819890786D-02  /
   data green(36,14, 3) /   2.5813033149764141D-02  /
   data green(36,14, 4) /   2.5752969645870986D-02  /
   data green(36,14, 5) /   2.5676359879110915D-02  /
   data green(36,14, 6) /   2.5583650902673196D-02  /
   data green(36,14, 7) /   2.5475374121332651D-02  /
   data green(36,14, 8) /   2.5352137047635840D-02  /
   data green(36,14, 9) /   2.5214614188681546D-02  /
   data green(36,14,10) /   2.5063537320125145D-02  /
   data green(36,14,11) /   2.4899685405629476D-02  /
   data green(36,14,12) /   2.4723874411138023D-02  /
   data green(36,14,13) /   2.4536947245509217D-02  /
   data green(36,14,14) /   2.4339764034058714D-02  /
   data green(36,15, 0) /   2.5642583265360275D-02  /
   data green(36,15, 1) /   2.5634146504652365D-02  /
   data green(36,15, 2) /   2.5608886277568077D-02  /
   data green(36,15, 3) /   2.5566951921977700D-02  /
   data green(36,15, 4) /   2.5508589603737097D-02  /
   data green(36,15, 5) /   2.5434138313489129D-02  /
   data green(36,15, 6) /   2.5344024442459912D-02  /
   data green(36,15, 7) /   2.5238755106014636D-02  /
   data green(36,15, 8) /   2.5118910416542004D-02  /
   data green(36,15, 9) /   2.4985134930158324D-02  /
   data green(36,15,10) /   2.4838128504196973D-02  /
   data green(36,15,11) /   2.4678636804613935D-02  /
   data green(36,15,12) /   2.4507441695056330D-02  /
   data green(36,15,13) /   2.4325351723675566D-02  /
   data green(36,15,14) /   2.4133192901440156D-02  /
   data green(36,15,15) /   2.3931799938534143D-02  /
   data green(36,16, 0) /   2.5384926108854830D-02  /
   data green(36,16, 1) /   2.5376741549103222D-02  /
   data green(36,16, 2) /   2.5352235449829404D-02  /
   data green(36,16, 3) /   2.5311549780607980D-02  /
   data green(36,16, 4) /   2.5254918614037607D-02  /
   data green(36,16, 5) /   2.5182664395278137D-02  /
   data green(36,16, 6) /   2.5095192883985783D-02  /
   data green(36,16, 7) /   2.4992986922935979D-02  /
   data green(36,16, 8) /   2.4876599217904225D-02  /
   data green(36,16, 9) /   2.4746644334772236D-02  /
   data green(36,16,10) /   2.4603790131807711D-02  /
   data green(36,16,11) /   2.4448748847709780D-02  /
   data green(36,16,12) /   2.4282268059964612D-02  /
   data green(36,16,13) /   2.4105121714413320D-02  /
   data green(36,16,14) /   2.3918101407121240D-02  /
   data green(36,16,15) /   2.3722008075268738D-02  /
   data green(36,16,16) /   2.3517644226514419D-02  /
   data green(36,17, 0) /   2.5119025613903983D-02  /
   data green(36,17, 1) /   2.5111095960243898D-02  /
   data green(36,17, 2) /   2.5087352129604331D-02  /
   data green(36,17, 3) /   2.5047928797753052D-02  /
   data green(36,17, 4) /   2.4993048062416173D-02  /
   data green(36,17, 5) /   2.4923015977321775D-02  /
   data green(36,17, 6) /   2.4838217849546033D-02  /
   data green(36,17, 7) /   2.4739112440683986D-02  /
   data green(36,17, 8) /   2.4626225240213254D-02  /
   data green(36,17, 9) /   2.4500140999318940D-02  /
   data green(36,17,10) /   2.4361495724895960D-02  /
   data green(36,17,11) /   2.4210968336477554D-02  /
   data green(36,17,12) /   2.4049272183993622D-02  /
   data green(36,17,13) /   2.3877146612483174D-02  /
   data green(36,17,14) /   2.3695348742417396D-02  /
   data green(36,17,15) /   2.3504645612553315D-02  /
   data green(36,17,16) /   2.3305806807710803D-02  /
   data green(36,17,17) /   2.3099597667979185D-02  /
   data green(36,18, 0) /   2.4845965805696153D-02  /
   data green(36,18, 1) /   2.4838292316007371D-02  /
   data green(36,18, 2) /   2.4815314568655329D-02  /
   data green(36,18, 3) /   2.4777160066539801D-02  /
   data green(36,18, 4) /   2.4724039127202703D-02  /
   data green(36,18, 5) /   2.4656241671517792D-02  /
   data green(36,18, 6) /   2.4574132863510238D-02  /
   data green(36,18, 7) /   2.4478147728840376D-02  /
   data green(36,18, 8) /   2.4368784905001887D-02  /
   data green(36,18, 9) /   2.4246599694726999D-02  /
   data green(36,18,10) /   2.4112196604976023D-02  /
   data green(36,18,11) /   2.3966221557223149D-02  /
   data green(36,18,12) /   2.3809353950971533D-02  /
   data green(36,18,13) /   2.3642298752351267D-02  /
   data green(36,18,14) /   2.3465778764348640D-02  /
   data green(36,18,15) /   2.3280527215934141D-02  /
   data green(36,18,16) /   2.3087280785404188D-02  /
   data green(36,18,17) /   2.2886773149908161D-02  /
   data green(36,18,18) /   2.2679729129564914D-02  /
   data green(36,19, 0) /   2.4566791085686168D-02  /
   data green(36,19, 1) /   2.4559373693506286D-02  /
   data green(36,19, 2) /   2.4537161884393553D-02  /
   data green(36,19, 3) /   2.4500276149035220D-02  /
   data green(36,19, 4) /   2.4448915284814961D-02  /
   data green(36,19, 5) /   2.4383353428026659D-02  /
   data green(36,19, 6) /   2.4303936021511614D-02  /
   data green(36,19, 7) /   2.4211074833089024D-02  /
   data green(36,19, 8) /   2.4105242163456202D-02  /
   data green(36,19, 9) /   2.3986964399265388D-02  /
   data green(36,19,10) /   2.3856815077387346D-02  /
   data green(36,19,11) /   2.3715407629923232D-02  /
   data green(36,19,12) /   2.3563387976684436D-02  /
   data green(36,19,13) /   2.3401427123311818D-02  /
   data green(36,19,14) /   2.3230213909880160D-02  /
   data green(36,19,15) /   2.3050448037819266D-02  /
   data green(36,19,16) /   2.2862833483431089D-02  /
   data green(36,19,17) /   2.2668072385325599D-02  /
   data green(36,19,18) /   2.2466859471789394D-02  /
   data green(36,19,19) /   2.2259877073352759D-02  /
   data green(36,20, 0) /   2.4282500049397408D-02  /
   data green(36,20, 1) /   2.4275337491444785D-02  /
   data green(36,20, 2) /   2.4253887896177179D-02  /
   data green(36,20, 3) /   2.4218264935872615D-02  /
   data green(36,20, 4) /   2.4168656202409329D-02  /
   data green(36,20, 5) /   2.4105320470962394D-02  /
   data green(36,20, 6) /   2.4028583979538529D-02  /
   data green(36,20, 7) /   2.3938835828401756D-02  /
   data green(36,20, 8) /   2.3836522624668378D-02  /
   data green(36,20, 9) /   2.3722142513023872D-02  /
   data green(36,20,10) /   2.3596238743223873D-02  /
   data green(36,20,11) /   2.3459392928734394D-02  /
   data green(36,20,12) /   2.3312218148835086D-02  /
   data green(36,20,13) /   2.3155352039328073D-02  /
   data green(36,20,14) /   2.2989450005464540D-02  /
   data green(36,20,15) /   2.2815178675764262D-02  /
   data green(36,20,16) /   2.2633209698067130D-02  /
   data green(36,20,17) /   2.2444213960425705D-02  /
   data green(36,20,18) /   2.2248856300259964D-02  /
   data green(36,20,19) /   2.2047790746377680D-02  /
   data green(36,20,20) /   2.1841656320695943D-02  /
   data green(36,21, 0) /   2.3994040683977518D-02  /
   data green(36,21, 1) /   2.3987130629098508D-02  /
   data green(36,21, 2) /   2.3966436328668538D-02  /
   data green(36,21, 3) /   2.3932064857246110D-02  /
   data green(36,21, 4) /   2.3884192960260526D-02  /
   data green(36,21, 5) /   2.3823064536540366D-02  /
   data green(36,21, 6) /   2.3748987212973194D-02  /
   data green(36,21, 7) /   2.3662328104885683D-02  /
   data green(36,21, 8) /   2.3563508875008680D-02  /
   data green(36,21, 9) /   2.3453000218274931D-02  /
   data green(36,21,10) /   2.3331315908805838D-02  /
   data green(36,21,11) /   2.3199006549214023D-02  /
   data green(36,21,12) /   2.3056653161005112D-02  /
   data green(36,21,13) /   2.2904860748893651D-02  /
   data green(36,21,14) /   2.2744251961930126D-02  /
   data green(36,21,15) /   2.2575460961287930D-02  /
   data green(36,21,16) /   2.2399127589256700D-02  /
   data green(36,21,17) /   2.2215891917316224D-02  /
   data green(36,21,18) /   2.2026389233957618D-02  /
   data green(36,21,19) /   2.1831245515916618D-02  /
   data green(36,21,20) /   2.1631073410311549D-02  /
   data green(36,21,21) /   2.1426468740323508D-02  /
   data green(36,22, 0) /   2.3702306853132342D-02  /
   data green(36,22, 1) /   2.3695646029952212D-02  /
   data green(36,22, 2) /   2.3675697291723492D-02  /
   data green(36,22, 3) /   2.3642561357012763D-02  /
   data green(36,22, 4) /   2.3596404518700244D-02  /
   data green(36,22, 5) /   2.3537456332352019D-02  /
   data green(36,22, 6) /   2.3466006468743590D-02  /
   data green(36,22, 7) /   2.3382400814504739D-02  /
   data green(36,22, 8) /   2.3287036922311856D-02  /
   data green(36,22, 9) /   2.3180358925214153D-02  /
   data green(36,22,10) /   2.3062852038193003D-02  /
   data green(36,22,11) /   2.2935036773839255D-02  /
   data green(36,22,12) /   2.2797462998271862D-02  /
   data green(36,22,13) /   2.2650703948512235D-02  /
   data green(36,22,14) /   2.2495350324053785D-02  /
   data green(36,22,15) /   2.2332004554019558D-02  /
   data green(36,22,16) /   2.2161275327852663D-02  /
   data green(36,22,17) /   2.1983772462701983D-02  /
   data green(36,22,18) /   2.1800102165292743D-02  /
   data green(36,22,19) /   2.1610862730768166D-02  /
   data green(36,22,20) /   2.1416640706322813D-02  /
   data green(36,22,21) /   2.1218007533868426D-02  /
   data green(36,22,22) /   2.1015516673811279D-02  /
   data green(36,23, 0) /   2.3408135957067644D-02  /
   data green(36,23, 1) /   2.3401720278099641D-02  /
   data green(36,23, 2) /   2.3382504926330417D-02  /
   data green(36,23, 3) /   2.3350584521738122D-02  /
   data green(36,23, 4) /   2.3306115324622175D-02  /
   data green(36,23, 5) /   2.3249313116726861D-02  /
   data green(36,23, 6) /   2.3180450314195737D-02  /
   data green(36,23, 7) /   2.3099852387524773D-02  /
   data green(36,23, 8) /   2.3007893679455681D-02  /
   data green(36,23, 9) /   2.2904992723758472D-02  /
   data green(36,23,10) /   2.2791607175777606D-02  /
   data green(36,23,11) /   2.2668228469368407D-02  /
   data green(36,23,12) /   2.2535376314569117D-02  /
   data green(36,23,13) /   2.2393593146367890D-02  /
   data green(36,23,14) /   2.2243438627724717D-02  /
   data green(36,23,15) /   2.2085484300190559D-02  /
   data green(36,23,16) /   2.1920308463692917D-02  /
   data green(36,23,17) /   2.1748491354001590D-02  /
   data green(36,23,18) /   2.1570610672703723D-02  /
   data green(36,23,19) /   2.1387237510792433D-02  /
   data green(36,23,20) /   2.1198932693721063D-02  /
   data green(36,23,21) /   2.1006243563406270D-02  /
   data green(36,23,22) /   2.0809701201490546D-02  /
   data green(36,23,23) /   2.0609818088405642D-02  /
   data green(36,24, 0) /   2.3112307642562016D-02  /
   data green(36,24, 1) /   2.3106132322990934D-02  /
   data green(36,24, 2) /   2.3087636093561036D-02  /
   data green(36,24, 3) /   2.3056907743911226D-02  /
   data green(36,24, 4) /   2.3014093939880631D-02  /
   data green(36,24, 5) /   2.2959397284385639D-02  /
   data green(36,24, 6) /   2.2893073673464191D-02  /
   data green(36,24, 7) /   2.2815429014639978D-02  /
   data green(36,24, 8) /   2.2726815388979284D-02  /
   data green(36,24, 9) /   2.2627626749143610D-02  /
   data green(36,24,10) /   2.2518294253091715D-02  /
   data green(36,24,11) /   2.2399281336763083D-02  /
   data green(36,24,12) /   2.2271078629183868D-02  /
   data green(36,24,13) /   2.2134198810245621D-02  /
   data green(36,24,14) /   2.1989171505330689D-02  /
   data green(36,24,15) /   2.1836538302501225D-02  /
   data green(36,24,16) /   2.1676847967703392D-02  /
   data green(36,24,17) /   2.1510651921949564D-02  /
   data green(36,24,18) /   2.1338500032301290D-02  /
   data green(36,24,19) /   2.1160936756208608D-02  /
   data green(36,24,20) /   2.0978497666827502D-02  /
   data green(36,24,21) /   2.0791706375710348D-02  /
   data green(36,24,22) /   2.0601071859037475D-02  /
   data green(36,24,23) /   2.0407086184530555D-02  /
   data green(36,24,24) /   2.0210222628481819D-02  /
   data green(36,25, 0) /   2.2815543432280077D-02  /
   data green(36,25, 1) /   2.2809603102070398D-02  /
   data green(36,25, 2) /   2.2791809977368294D-02  /
   data green(36,25, 3) /   2.2762247292289130D-02  /
   data green(36,25, 4) /   2.2721052567462981D-02  /
   data green(36,25, 5) /   2.2668415837922503D-02  /
   data green(36,25, 6) /   2.2604577235207687D-02  /
   data green(36,25, 7) /   2.2529823983571250D-02  /
   data green(36,25, 8) /   2.2444486882961492D-02  /
   data green(36,25, 9) /   2.2348936361386660D-02  /
   data green(36,25,10) /   2.2243578186059560D-02  /
   data green(36,25,11) /   2.2128848926292188D-02  /
   data green(36,25,12) /   2.2005211261529864D-02  /
   data green(36,25,13) /   2.1873149225404131D-02  /
   data green(36,25,14) /   2.1733163471587668D-02  /
   data green(36,25,15) /   2.1585766639983378D-02  /
   data green(36,25,16) /   2.1431478892862232D-02  /
   data green(36,25,17) /   2.1270823680489853D-02  /
   data green(36,25,18) /   2.1104323785046560D-02  /
   data green(36,25,19) /   2.0932497680718122D-02  /
   data green(36,25,20) /   2.0755856237121392D-02  /
   data green(36,25,21) /   2.0574899783073461D-02  /
   data green(36,25,22) /   2.0390115538384393D-02  /
   data green(36,25,23) /   2.0201975413042771D-02  /
   data green(36,25,24) /   2.0010934165992585D-02  /
   data green(36,25,25) /   1.9817427909724541D-02  /
   data green(36,26, 0) /   2.2518507141818064D-02  /
   data green(36,26, 1) /   2.2512795950189818D-02  /
   data green(36,26, 2) /   2.2495688471288577D-02  /
   data green(36,26, 3) /   2.2467262661360596D-02  /
   data green(36,26, 4) /   2.2427647350082902D-02  /
   data green(36,26, 5) /   2.2377020623102938D-02  /
   data green(36,26, 6) /   2.2315607613661151D-02  /
   data green(36,26, 7) /   2.2243677756611829D-02  /
   data green(36,26, 8) /   2.2161541569646594D-02  /
   data green(36,26, 9) /   2.2069547035521456D-02  /
   data green(36,26,10) /   2.1968075665350081D-02  /
   data green(36,26,11) /   2.1857538326463082D-02  /
   data green(36,26,12) /   2.1738370918995401D-02  /
   data green(36,26,13) /   2.1611029983428967D-02  /
   data green(36,26,14) /   2.1475988317074850D-02  /
   data green(36,26,15) /   2.1333730671287986D-02  /
   data green(36,26,16) /   2.1184749593490940D-02  /
   data green(36,26,17) /   2.1029541469275957D-02  /
   data green(36,26,18) /   2.0868602810392940D-02  /
   data green(36,26,19) /   2.0702426824725357D-02  /
   data green(36,26,20) /   2.0531500294769332D-02  /
   data green(36,26,21) /   2.0356300781973907D-02  /
   data green(36,26,22) /   2.0177294165819298D-02  /
   data green(36,26,23) /   1.9994932518886285D-02  /
   data green(36,26,24) /   1.9809652312524840D-02  /
   data green(36,26,25) /   1.9621872942124547D-02  /
   data green(36,26,26) /   1.9431995556437643D-02  /
   data green(36,27, 0) /   2.2221805956635008D-02  /
   data green(36,27, 1) /   2.2216317668296497D-02  /
   data green(36,27, 2) /   2.2199877222581158D-02  /
   data green(36,27, 3) /   2.2172557575166157D-02  /
   data green(36,27, 4) /   2.2134479318781056D-02  /
   data green(36,27, 5) /   2.2085809208522284D-02  /
   data green(36,27, 6) /   2.2026758147072598D-02  /
   data green(36,27, 7) /   2.1957578677224653D-02  /
   data green(36,27, 8) /   2.1878562039411152D-02  /
   data green(36,27, 9) /   2.1790034860079136D-02  /
   data green(36,27,10) /   2.1692355542501558D-02  /
   data green(36,27,11) /   2.1585910434902392D-02  /
   data green(36,27,12) /   2.1471109851616852D-02  /
   data green(36,27,13) /   2.1348384021557577D-02  /
   data green(36,27,14) /   2.1218179034749221D-02  /
   data green(36,27,15) /   2.1080952852433020D-02  /
   data green(36,27,16) /   2.0937171439588141D-02  /
   data green(36,27,17) /   2.0787305071041672D-02  /
   data green(36,27,18) /   2.0631824854025219D-02  /
   data green(36,27,19) /   2.0471199501440288D-02  /
   data green(36,27,20) /   2.0305892381541506D-02  /
   data green(36,27,21) /   2.0136358861513725D-02  /
   data green(36,27,22) /   1.9963043954733299D-02  /
   data green(36,27,23) /   1.9786380274534384D-02  /
   data green(36,27,24) /   1.9606786291166280D-02  /
   data green(36,27,25) /   1.9424664883393297D-02  /
   data green(36,27,26) /   1.9240402171877224D-02  /
   data green(36,27,27) /   1.9054366618078709D-02  /
   data green(36,28, 0) /   2.1925992047904838D-02  /
   data green(36,28, 1) /   2.1920720130727860D-02  /
   data green(36,28, 2) /   2.1904927214020813D-02  /
   data green(36,28, 3) /   2.1878681527151483D-02  /
   data green(36,28, 4) /   2.1842095875228749D-02  /
   data green(36,28, 5) /   2.1795326295869305D-02  /
   data green(36,28, 6) /   2.1738570222813954D-02  /
   data green(36,28, 7) /   2.1672064198481465D-02  /
   data green(36,28, 8) /   2.1596081186772996D-02  /
   data green(36,28, 9) /   2.1510927544785369D-02  /
   data green(36,28,10) /   2.1416939717365781D-02  /
   data green(36,28,11) /   2.1314480721557411D-02  /
   data green(36,28,12) /   2.1203936488962988D-02  /
   data green(36,28,13) /   2.1085712133005297D-02  /
   data green(36,28,14) /   2.0960228205185755D-02  /
   data green(36,28,15) /   2.0827916999992362D-02  /
   data green(36,28,16) /   2.0689218962387731D-02  /
   data green(36,28,17) /   2.0544579245141654D-02  /
   data green(36,28,18) /   2.0394444455986298D-02  /
   data green(36,28,19) /   2.0239259626981202D-02  /
   data green(36,28,20) /   2.0079465430864769D-02  /
   data green(36,28,21) /   1.9915495661788692D-02  /
   data green(36,28,22) /   1.9747774990887340D-02  /
   data green(36,28,23) /   1.9576717000783610D-02  /
   data green(36,28,24) /   1.9402722497489373D-02  /
   data green(36,28,25) /   1.9226178093291180D-02  /
   data green(36,28,26) /   1.9047455050148796D-02  /
   data green(36,28,27) /   1.8866908369872926D-02  /
   data green(36,28,28) /   1.8684876114855625D-02  /
   data green(36,29, 0) /   2.1631564615466170D-02  /
   data green(36,29, 1) /   2.1626502319535414D-02  /
   data green(36,29, 2) /   2.1611336772502838D-02  /
   data green(36,29, 3) /   2.1586131746431040D-02  /
   data green(36,29, 4) /   2.1550992699796781D-02  /
   data green(36,29, 5) /   2.1506065554990628D-02  /
   data green(36,29, 6) /   2.1451535025920589D-02  /
   data green(36,29, 7) /   2.1387622533052233D-02  /
   data green(36,29, 8) /   2.1314583751466450D-02  /
   data green(36,29, 9) /   2.1232705844137372D-02  /
   data green(36,29,10) /   2.1142304437455264D-02  /
   data green(36,29,11) /   2.1043720398960249D-02  /
   data green(36,29,12) /   2.0937316478321085D-02  /
   data green(36,29,13) /   2.0823473871876365D-02  /
   data green(36,29,14) /   2.0702588768716498D-02  /
   data green(36,29,15) /   2.0575068932536834D-02  /
   data green(36,29,16) /   2.0441330368593236D-02  /
   data green(36,29,17) /   2.0301794119315889D-02  /
   data green(36,29,18) /   2.0156883225770264D-02  /
   data green(36,29,19) /   2.0007019885466781D-02  /
   data green(36,29,20) /   1.9852622830267484D-02  /
   data green(36,29,21) /   1.9694104941539812D-02  /
   data green(36,29,22) /   1.9531871113450602D-02  /
   data green(36,29,23) /   1.9366316369525204D-02  /
   data green(36,29,24) /   1.9197824232422698D-02  /
   data green(36,29,25) /   1.9026765342369988D-02  /
   data green(36,29,26) /   1.8853496315887071D-02  /
   data green(36,29,27) /   1.8678358833329575D-02  /
   data green(36,29,28) /   1.8501678941351608D-02  /
   data green(36,29,29) /   1.8323766554610967D-02  /
   data green(36,30, 0) /   2.1338972256614644D-02  /
   data green(36,30, 1) /   2.1334112684781376D-02  /
   data green(36,30, 2) /   2.1319553903997354D-02  /
   data green(36,30, 3) /   2.1295355490982847D-02  /
   data green(36,30, 4) /   2.1261615987275700D-02  /
   data green(36,30, 5) /   2.1218471787377572D-02  /
   data green(36,30, 6) /   2.1166095616768694D-02  /
   data green(36,30, 7) /   2.1104694632876161D-02  /
   data green(36,30, 8) /   2.1034508189447243D-02  /
   data green(36,30, 9) /   2.0955805310740153D-02  /
   data green(36,30,10) /   2.0868881926343195D-02  /
   data green(36,30,11) /   2.0774057920194967D-02  /
   data green(36,30,12) /   2.0671674048501818D-02  /
   data green(36,30,13) /   2.0562088780802294D-02  /
   data green(36,30,14) /   2.0445675116546078D-02  /
   data green(36,30,15) /   2.0322817426413266D-02  /
   data green(36,30,16) /   2.0193908363417674D-02  /
   data green(36,30,17) /   2.0059345883848716D-02  /
   data green(36,30,18) /   1.9919530412556480D-02  /
   data green(36,30,19) /   1.9774862181208008D-02  /
   data green(36,30,20) /   1.9625738762165328D-02  /
   data green(36,30,21) /   1.9472552814751017D-02  /
   data green(36,30,22) /   1.9315690055045579D-02  /
   data green(36,30,23) /   1.9155527455135948D-02  /
   data green(36,30,24) /   1.8992431673006970D-02  /
   data green(36,30,25) /   1.8826757710107257D-02  /
   data green(36,30,26) /   1.8658847790063271D-02  /
   data green(36,30,27) /   1.8489030449074215D-02  /
   data green(36,30,28) /   1.8317619826182095D-02  /
   data green(36,30,29) /   1.8144915139847961D-02  /
   data green(36,30,30) /   1.7971200336033329D-02  /
   data green(36,31, 0) /   2.1048615570780482D-02  /
   data green(36,31, 1) /   2.1043951741005599D-02  /
   data green(36,31, 2) /   2.1029978865513534D-02  /
   data green(36,31, 3) /   2.1006752579203342D-02  /
   data green(36,31, 4) /   2.0974364922746139D-02  /
   data green(36,31, 5) /   2.0932943331937681D-02  /
   data green(36,31, 6) /   2.0882649253402788D-02  /
   data green(36,31, 7) /   2.0823676415949625D-02  /
   data green(36,31, 8) /   2.0756248793445963D-02  /
   data green(36,31, 9) /   2.0680618300446946D-02  /
   data green(36,31,10) /   2.0597062265809714D-02  /
   data green(36,31,11) /   2.0505880732110432D-02  /
   data green(36,31,12) /   2.0407393629829945D-02  /
   data green(36,31,13) /   2.0301937875046264D-02  /
   data green(36,31,14) /   2.0189864437874083D-02  /
   data green(36,31,15) /   2.0071535426271610D-02  /
   data green(36,31,16) /   1.9947321226276644D-02  /
   data green(36,31,17) /   1.9817597735435542D-02  /
   data green(36,31,18) /   1.9682743721361869D-02  /
   data green(36,31,19) /   1.9543138332210844D-02  /
   data green(36,31,20) /   1.9399158780575729D-02  /
   data green(36,31,21) /   1.9251178217077263D-02  /
   data green(36,31,22) /   1.9099563804879399D-02  /
   data green(36,31,23) /   1.8944675001644110D-02  /
   data green(36,31,24) /   1.8786862051132773D-02  /
   data green(36,31,25) /   1.8626464682834108D-02  /
   data green(36,31,26) /   1.8463811014692296D-02  /
   data green(36,31,27) /   1.8299216651237241D-02  /
   data green(36,31,28) /   1.8132983967179155D-02  /
   data green(36,31,29) /   1.7965401564801173D-02  /
   data green(36,31,30) /   1.7796743892232108D-02  /
   data green(36,31,31) /   1.7627271008865700D-02  /
   data green(36,32, 0) /   2.0760849921607181D-02  /
   data green(36,32, 1) /   2.0756374821493489D-02  /
   data green(36,32, 2) /   2.0742966896049951D-02  /
   data green(36,32, 3) /   2.0720678082374135D-02  /
   data green(36,32, 4) /   2.0689594320954643D-02  /
   data green(36,32, 5) /   2.0649834637436658D-02  /
   data green(36,32, 6) /   2.0601549884151783D-02  /
   data green(36,32, 7) /   2.0544921167341270D-02  /
   data green(36,32, 8) /   2.0480157991862757D-02  /
   data green(36,32, 9) /   2.0407496159981726D-02  /
   data green(36,32,10) /   2.0327195464488856D-02  /
   data green(36,32,11) /   2.0239537218786557D-02  /
   data green(36,32,12) /   2.0144821667742171D-02  /
   data green(36,32,13) /   2.0043365323051159D-02  /
   data green(36,32,14) /   1.9935498265678010D-02  /
   data green(36,32,15) /   1.9821561455769084D-02  /
   data green(36,32,16) /   1.9701904087413707D-02  /
   data green(36,32,17) /   1.9576881021936606D-02  /
   data green(36,32,18) /   1.9446850329216822D-02  /
   data green(36,32,19) /   1.9312170962022871D-02  /
   data green(36,32,20) /   1.9173200583701040D-02  /
   data green(36,32,21) /   1.9030293564906837D-02  /
   data green(36,32,22) /   1.8883799160565401D-02  /
   data green(36,32,23) /   1.8734059873993604D-02  /
   data green(36,32,24) /   1.8581410011205751D-02  /
   data green(36,32,25) /   1.8426174424916723D-02  /
   data green(36,32,26) /   1.8268667444694082D-02  /
   data green(36,32,27) /   1.8109191987111684D-02  /
   data green(36,32,28) /   1.7948038837625597D-02  /
   data green(36,32,29) /   1.7785486094213373D-02  /
   data green(36,32,30) /   1.7621798761565888D-02  /
   data green(36,32,31) /   1.7457228483761935D-02  /
   data green(36,32,32) /   1.7292013402850050D-02  /
   data green(36,33, 0) /   2.0475988289165675D-02  /
   data green(36,33, 1) /   2.0471694923160706D-02  /
   data green(36,33, 2) /   2.0458831039590408D-02  /
   data green(36,33, 3) /   2.0437445111506601D-02  /
   data green(36,33, 4) /   2.0407617363232790D-02  /
   data green(36,33, 5) /   2.0369458936385941D-02  /
   data green(36,33, 6) /   2.0323110746214200D-02  /
   data green(36,33, 7) /   2.0268742051186776D-02  /
   data green(36,33, 8) /   2.0206548763994948D-02  /
   data green(36,33, 9) /   2.0136751536438287D-02  /
   data green(36,33,10) /   2.0059593653972866D-02  /
   data green(36,33,11) /   1.9975338777926194D-02  /
   data green(36,33,12) /   1.9884268574523797D-02  /
   data green(36,33,13) /   1.9786680269953831D-02  /
   data green(36,33,14) /   1.9682884169789697D-02  /
   data green(36,33,15) /   1.9573201179297117D-02  /
   data green(36,33,16) /   1.9457960358601475D-02  /
   data green(36,33,17) /   1.9337496543525696D-02  /
   data green(36,33,18) /   1.9212148059283646D-02  /
   data green(36,33,19) /   1.9082254550280679D-02  /
   data green(36,33,20) /   1.8948154945181737D-02  /
   data green(36,33,21) /   1.8810185572291378D-02  /
   data green(36,33,22) /   1.8668678436271118D-02  /
   data green(36,33,23) /   1.8523959663397975D-02  /
   data green(36,33,24) /   1.8376348119024277D-02  /
   data green(36,33,25) /   1.8226154197693856D-02  /
   data green(36,33,26) /   1.8073678783543004D-02  /
   data green(36,33,27) /   1.7919212376189146D-02  /
   data green(36,33,28) /   1.7763034375292159D-02  /
   data green(36,33,29) /   1.7605412515353538D-02  /
   data green(36,33,30) /   1.7446602441082273D-02  /
   data green(36,33,31) /   1.7286847412771089D-02  /
   data green(36,33,32) /   1.7126378130565771D-02  /
   data green(36,33,33) /   1.6965412666231510D-02  /
   data green(36,34, 0) /   2.0194304155692240D-02  /
   data green(36,34, 1) /   2.0190185585494742D-02  /
   data green(36,34, 2) /   2.0177845003742301D-02  /
   data green(36,34, 3) /   2.0157327642425435D-02  /
   data green(36,34, 4) /   2.0128708376194477D-02  /
   data green(36,34, 5) /   2.0092090965099396D-02  /
   data green(36,34, 6) /   2.0047607015541879D-02  /
   data green(36,34, 7) /   1.9995414679721876D-02  /
   data green(36,34, 8) /   1.9935697118511099D-02  /
   data green(36,34, 9) /   1.9868660756553524D-02  /
   data green(36,34,10) /   1.9794533361385334D-02  /
   data green(36,34,11) /   1.9713561980426002D-02  /
   data green(36,34,12) /   1.9626010770805138D-02  /
   data green(36,34,13) /   1.9532158757175028D-02  /
   data green(36,34,14) /   1.9432297551974887D-02  /
   data green(36,34,15) /   1.9326729071143071D-02  /
   data green(36,34,16) /   1.9215763276123690D-02  /
   data green(36,34,17) /   1.9099715970308649D-02  /
   data green(36,34,18) /   1.8978906674923725D-02  /
   data green(36,34,19) /   1.8853656605941082D-02  /
   data green(36,34,20) /   1.8724286770009198D-02  /
   data green(36,34,21) /   1.8591116193752654D-02  /
   data green(36,34,22) /   1.8454460297216958D-02  /
   data green(36,34,23) /   1.8314629418805631D-02  /
   data green(36,34,24) /   1.8171927495855507D-02  /
   data green(36,34,25) /   1.8026650902075038D-02  /
   data green(36,34,26) /   1.7879087440470887D-02  /
   data green(36,34,27) /   1.7729515488133465D-02  /
   data green(36,34,28) /   1.7578203287350835D-02  /
   data green(36,34,29) /   1.7425408375971108D-02  /
   data green(36,34,30) /   1.7271377148722584D-02  /
   data green(36,34,31) /   1.7116344540307128D-02  /
   data green(36,34,32) /   1.6960533820480208D-02  /
   data green(36,34,33) /   1.6804156490989170D-02  /
   data green(36,34,34) /   1.6647412274127506D-02  /
   data green(36,35, 0) /   1.9916034378114861D-02  /
   data green(36,35, 1) /   1.9912083756846392D-02  /
   data green(36,35, 2) /   1.9900246007395959D-02  /
   data green(36,35, 3) /   1.9880563332615075D-02  /
   data green(36,35, 4) /   1.9853105605945286D-02  /
   data green(36,35, 5) /   1.9817969683932173D-02  /
   data green(36,35, 6) /   1.9775278462392350D-02  /
   data green(36,35, 7) /   1.9725179694164454D-02  /
   data green(36,35, 8) /   1.9667844590514010D-02  /
   data green(36,35, 9) /   1.9603466231725694D-02  /
   data green(36,35,10) /   1.9532257815124233D-02  /
   data green(36,35,11) /   1.9454450770663625D-02  /
   data green(36,35,12) /   1.9370292775298302D-02  /
   data green(36,35,13) /   1.9280045697614521D-02  /
   data green(36,35,14) /   1.9183983503698491D-02  /
   data green(36,35,15) /   1.9082390154021224D-02  /
   data green(36,35,16) /   1.8975557519315556D-02  /
   data green(36,35,17) /   1.8863783341113449D-02  /
   data green(36,35,18) /   1.8747369259910870D-02  /
   data green(36,35,19) /   1.8626618930948646D-02  /
   data green(36,35,20) /   1.8501836244450651D-02  /
   data green(36,35,21) /   1.8373323663950936D-02  /
   data green(36,35,22) /   1.8241380693161641D-02  /
   data green(36,35,23) /   1.8106302478767815D-02  /
   data green(36,35,24) /   1.7968378553647694D-02  /
   data green(36,35,25) /   1.7827891722363230D-02  /
   data green(36,35,26) /   1.7685117088381826D-02  /
   data green(36,35,27) /   1.7540321220402302D-02  /
   data green(36,35,28) /   1.7393761453375432D-02  /
   data green(36,35,29) /   1.7245685318336621D-02  /
   data green(36,35,30) /   1.7096330093992775D-02  /
   data green(36,35,31) /   1.6945922472115391D-02  /
   data green(36,35,32) /   1.6794678328163583D-02  /
   data green(36,35,33) /   1.6642802588170191D-02  /
   data green(36,35,34) /   1.6490489182744612D-02  /
   data green(36,35,35) /   1.6337921079046507D-02  /
   data green(36,36, 0) /   1.9641382009603400D-02  /
   data green(36,36, 1) /   1.9637592610315203D-02  /
   data green(36,36, 2) /   1.9626237579676520D-02  /
   data green(36,36, 3) /   1.9607356292149915D-02  /
   data green(36,36, 4) /   1.9581013950201519D-02  /
   data green(36,36, 5) /   1.9547300960211299D-02  /
   data green(36,36, 6) /   1.9506332075212816D-02  /
   data green(36,36, 7) /   1.9458245320313046D-02  /
   data green(36,36, 8) /   1.9403200720324535D-02  /
   data green(36,36, 9) /   1.9341378852243254D-02  /
   data green(36,36,10) /   1.9272979247651659D-02  /
   data green(36,36,11) /   1.9198218671872814D-02  /
   data green(36,36,12) /   1.9117329307729459D-02  /
   data green(36,36,13) /   1.9030556872083177D-02  /
   data green(36,36,14) /   1.8938158692976607D-02  /
   data green(36,36,15) /   1.8840401774234943D-02  /
   data green(36,36,16) /   1.8737560872873787D-02  /
   data green(36,36,17) /   1.8629916612697094D-02  /
   data green(36,36,18) /   1.8517753655145634D-02  /
   data green(36,36,19) /   1.8401358945870173D-02  /
   data green(36,36,20) /   1.8281020052751163D-02  /
   data green(36,36,21) /   1.8157023608259503D-02  /
   data green(36,36,22) /   1.8029653866232970D-02  /
   data green(36,36,23) /   1.7899191380404675D-02  /
   data green(36,36,24) /   1.7765911809422653D-02  /
   data green(36,36,25) /   1.7630084850694044D-02  /
   data green(36,36,26) /   1.7491973303207769D-02  /
   data green(36,36,27) /   1.7351832257562389D-02  /
   data green(36,36,28) /   1.7209908409762300D-02  /
   data green(36,36,29) /   1.7066439493951945D-02  /
   data green(36,36,30) /   1.6921653828126389D-02  /
   data green(36,36,31) /   1.6775769965979955D-02  /
   data green(36,36,32) /   1.6628996447412212D-02  /
   data green(36,36,33) /   1.6481531639785246D-02  /
   data green(36,36,34) /   1.6333563661793300D-02  /
   data green(36,36,35) /   1.6185270381741015D-02  /
   data green(36,36,36) /   1.6036819482107730D-02  /
   data green(37, 0, 0) /   2.7031972558915429D-02  /
   data green(37, 1, 0) /   2.7022083366600225D-02  /
   data green(37, 1, 1) /   2.7012205059666911D-02  /
   data green(37, 2, 0) /   2.6992481102104716D-02  /
   data green(37, 2, 1) /   2.6982635330964316D-02  /
   data green(37, 2, 2) /   2.6953162850456697D-02  /
   data green(37, 3, 0) /   2.6943360499241724D-02  /
   data green(37, 3, 1) /   2.6933568556370895D-02  /
   data green(37, 3, 2) /   2.6904256968216057D-02  /
   data green(37, 3, 3) /   2.6855617280348597D-02  /
   data green(37, 4, 0) /   2.6875042146304880D-02  /
   data green(37, 4, 1) /   2.6865324737887303D-02  /
   data green(37, 4, 2) /   2.6836235936423234D-02  /
   data green(37, 4, 3) /   2.6787964858477387D-02  /
   data green(37, 4, 4) /   2.6720822890435789D-02  /
   data green(37, 5, 0) /   2.6787966684383113D-02  /
   data green(37, 5, 1) /   2.6778343716844323D-02  /
   data green(37, 5, 2) /   2.6749537209396160D-02  /
   data green(37, 5, 3) /   2.6701733219031971D-02  /
   data green(37, 5, 4) /   2.6635238119267185D-02  /
   data green(37, 5, 5) /   2.6550473156348935D-02  /
   data green(37, 6, 0) /   2.6682686976528758D-02  /
   data green(37, 6, 1) /   2.6673177362952501D-02  /
   data green(37, 6, 2) /   2.6644709691892320D-02  /
   data green(37, 6, 3) /   2.6597466374929036D-02  /
   data green(37, 6, 4) /   2.6531747812081962D-02  /
   data green(37, 6, 5) /   2.6447967096407076D-02  /
   data green(37, 6, 6) /   2.6346642859459313D-02  /
   data green(37, 7, 0) /   2.6559858518494527D-02  /
   data green(37, 7, 1) /   2.6550480009025566D-02  /
   data green(37, 7, 2) /   2.6522404245325745D-02  /
   data green(37, 7, 3) /   2.6475809459411200D-02  /
   data green(37, 7, 4) /   2.6410989201484183D-02  /
   data green(37, 7, 5) /   2.6328347213201726D-02  /
   data green(37, 7, 6) /   2.6228390498650649D-02  /
   data green(37, 7, 7) /   2.6111720824339432D-02  /
   data green(37, 8, 0) /   2.6420228411457348D-02  /
   data green(37, 8, 1) /   2.6410997450770909D-02  /
   data green(37, 8, 2) /   2.6383362768241929D-02  /
   data green(37, 8, 3) /   2.6337497939298612D-02  /
   data green(37, 8, 4) /   2.6273688881215260D-02  /
   data green(37, 8, 5) /   2.6192328912440069D-02  /
   data green(37, 8, 6) /   2.6093912072453276D-02  /
   data green(37, 8, 7) /   2.5979024922843579D-02  /
   data green(37, 8, 8) /   2.5848337091864564D-02  /
   data green(37, 9, 0) /   2.6264623249482363D-02  /
   data green(37, 9, 1) /   2.6255554863302031D-02  /
   data green(37, 9, 2) /   2.6228406199436428D-02  /
   data green(37, 9, 3) /   2.6183345760673306D-02  /
   data green(37, 9, 4) /   2.6120651147766915D-02  /
   data green(37, 9, 5) /   2.6040704319182369D-02  /
   data green(37, 9, 6) /   2.5943985179188576D-02  /
   data green(37, 9, 7) /   2.5831063703666571D-02  /
   data green(37, 9, 8) /   2.5702590852690475D-02  /
   data green(37, 9, 9) /   2.5559288545840709D-02  /
   data green(37,10, 0) /   2.6093936288815327D-02  /
   data green(37,10, 1) /   2.6085044000357136D-02  /
   data green(37,10, 2) /   2.6058421806412348D-02  /
   data green(37,10, 3) /   2.6014232784289867D-02  /
   data green(37,10, 4) /   2.5952745637860448D-02  /
   data green(37,10, 5) /   2.5874330168998840D-02  /
   data green(37,10, 6) /   2.5779451149179414D-02  /
   data green(37,10, 7) /   2.5668660788774476D-02  /
   data green(37,10, 8) /   2.5542590039310416D-02  /
   data green(37,10, 9) /   2.5401938989710399D-02  /
   data green(37,10,10) /   2.5247466630791023D-02  /
   data green(37,11, 0) /   2.5909114263100719D-02  /
   data green(37,11, 1) /   2.5900410039271828D-02  /
   data green(37,11, 2) /   2.5874350119099679D-02  /
   data green(37,11, 3) /   2.5831091865499964D-02  /
   data green(37,11, 4) /   2.5770894609538981D-02  /
   data green(37,11, 5) /   2.5694115341819951D-02  /
   data green(37,11, 6) /   2.5601202878814375D-02  /
   data green(37,11, 7) /   2.5492690689567818D-02  /
   data green(37,11, 8) /   2.5369188603860265D-02  /
   data green(37,11, 9) /   2.5231373647485778D-02  /
   data green(37,11,10) /   2.5079980263248723D-02  /
   data green(37,11,11) /   2.4915790177762718D-02  /
   data green(37,12, 0) /   2.5711144190107216D-02  /
   data green(37,12, 1) /   2.5702638416345156D-02  /
   data green(37,12, 2) /   2.5677171850666595D-02  /
   data green(37,12, 3) /   2.5634895915900197D-02  /
   data green(37,12, 4) /   2.5576060197739926D-02  /
   data green(37,12, 5) /   2.5501008361404221D-02  /
   data green(37,12, 6) /   2.5410172620050117D-02  /
   data green(37,12, 7) /   2.5304066928131981D-02  /
   data green(37,12, 8) /   2.5183279106447361D-02  /
   data green(37,12, 9) /   2.5048462128962813D-02  /
   data green(37,12,10) /   2.4900324814095560D-02  /
   data green(37,12,11) /   2.4739622165096557D-02  /
   data green(37,12,12) /   2.4567145596341260D-02  /
   data green(37,13, 0) /   2.5501040484142926D-02  /
   data green(37,13, 1) /   2.5492741965983876D-02  /
   data green(37,13, 2) /   2.5467895116429176D-02  /
   data green(37,13, 3) /   2.5426645253792603D-02  /
   data green(37,13, 4) /   2.5369231946050470D-02  /
   data green(37,13, 5) /   2.5295985155363777D-02  /
   data green(37,13, 6) /   2.5207320012269316D-02  /
   data green(37,13, 7) /   2.5103730380534071D-02  /
   data green(37,13, 8) /   2.4985781405106124D-02  /
   data green(37,13, 9) /   2.4854101257683104D-02  /
   data green(37,13,10) /   2.4709372306608999D-02  /
   data green(37,13,11) /   2.4552321940212441D-02  /
   data green(37,13,12) /   2.4383713266011652D-02  /
   data green(37,13,13) /   2.4204335893616549D-02  /
   data green(37,14, 0) /   2.5279832647229931D-02  /
   data green(37,14, 1) /   2.5271748636065149D-02  /
   data green(37,14, 2) /   2.5247543221422454D-02  /
   data green(37,14, 3) /   2.5207355510925018D-02  /
   data green(37,14, 4) /   2.5151414877845446D-02  /
   data green(37,14, 5) /   2.5080037333572412D-02  /
   data green(37,14, 6) /   2.4993620607918321D-02  /
   data green(37,14, 7) /   2.4892638086256837D-02  /
   data green(37,14, 8) /   2.4777631781813962D-02  /
   data green(37,14, 9) /   2.4649204542244278D-02  /
   data green(37,14,10) /   2.4508011701383343D-02  /
   data green(37,14,11) /   2.4354752389844886D-02  /
   data green(37,14,12) /   2.4190160712528913D-02  /
   data green(37,14,13) /   2.4014996988169404D-02  /
   data green(37,14,14) /   2.3830039227129514D-02  /
   data green(37,15, 0) /   2.5048553764647988D-02  /
   data green(37,15, 1) /   2.5040690004677769D-02  /
   data green(37,15, 2) /   2.5017143240447227D-02  /
   data green(37,15, 3) /   2.4978046317098882D-02  /
   data green(37,15, 4) /   2.4923618325315366D-02  /
   data green(37,15, 5) /   2.4854161199566852D-02  /
   data green(37,15, 6) /   2.4770055101865317D-02  /
   data green(37,15, 7) /   2.4671752728324475D-02  /
   data green(37,15, 8) /   2.4559772703083353D-02  /
   data green(37,15, 9) /   2.4434692243679227D-02  /
   data green(37,15,10) /   2.4297139293249281D-02  /
   data green(37,15,11) /   2.4147784318033626D-02  /
   data green(37,15,12) /   2.3987331964056637D-02  /
   data green(37,15,13) /   2.3816512755493943D-02  /
   data green(37,15,14) /   2.3636075000292554D-02  /
   data green(37,15,15) /   2.3446777047479609D-02  /
   data green(37,16, 0) /   2.4808229979932799D-02  /
   data green(37,16, 1) /   2.4800590773051957D-02  /
   data green(37,16, 2) /   2.4777715564579233D-02  /
   data green(37,16, 3) /   2.4739730935264360D-02  /
   data green(37,16, 4) /   2.4686845687099085D-02  /
   data green(37,16, 5) /   2.4619347663221047D-02  /
   data green(37,16, 6) /   2.4537599429811204D-02  /
   data green(37,16, 7) /   2.4442032946008680D-02  /
   data green(37,16, 8) /   2.4333143373103144D-02  /
   data green(37,16, 9) /   2.4211482192519881D-02  /
   data green(37,16,10) /   2.4077649812914142D-02  /
   data green(37,16,11) /   2.3932287850036069D-02  /
   data green(37,16,12) /   2.3776071259353819D-02  /
   data green(37,16,13) /   2.3609700491518195D-02  /
   data green(37,16,14) /   2.3433893825685372D-02  /
   data green(37,16,15) /   2.3249380016708301D-02  /
   data green(37,16,16) /   2.3056891370553122D-02  /
   data green(37,17, 0) /   2.4559871073781523D-02  /
   data green(37,17, 1) /   2.4552459359018952D-02  /
   data green(37,17, 2) /   2.4530264538144841D-02  /
   data green(37,17, 3) /   2.4493406970536345D-02  /
   data green(37,17, 4) /   2.4442085237140466D-02  /
   data green(37,17, 5) /   2.4376573176244108D-02  /
   data green(37,17, 6) /   2.4297215855964255D-02  /
   data green(37,17, 7) /   2.4204424598704640D-02  /
   data green(37,17, 8) /   2.4098671196104475D-02  /
   data green(37,17, 9) /   2.3980481470008735D-02  /
   data green(37,17,10) /   2.3850428345274118D-02  /
   data green(37,17,11) /   2.3709124603763156D-02  /
   data green(37,17,12) /   2.3557215486032058D-02  /
   data green(37,17,13) /   2.3395371298670687D-02  /
   data green(37,17,14) /   2.3224280171938002D-02  /
   data green(37,17,15) /   2.3044641095338332D-02  /
   data green(37,17,16) /   2.2857157339253024D-02  /
   data green(37,17,17) /   2.2662530349814056D-02  /
   data green(37,18, 0) /   2.4304462223060749D-02  /
   data green(37,18, 1) /   2.4297279667226244D-02  /
   data green(37,18, 2) /   2.4275770262475155D-02  /
   data green(37,18, 3) /   2.4240048229570890D-02  /
   data green(37,18, 4) /   2.4190302061355147D-02  /
   data green(37,18, 5) /   2.4126791767220112D-02  /
   data green(37,18, 6) /   2.4049845126780951D-02  /
   data green(37,18, 7) /   2.3959853057776120D-02  /
   data green(37,18, 8) /   2.3857264224622768D-02  /
   data green(37,18, 9) /   2.3742579029838829D-02  /
   data green(37,18,10) /   2.3616343140295301D-02  /
   data green(37,18,11) /   2.3479140703933518D-02  /
   data green(37,18,12) /   2.3331587410471782D-02  /
   data green(37,18,13) /   2.3174323542317441D-02  /
   data green(37,18,14) /   2.3008007150209741D-02  /
   data green(37,18,15) /   2.2833307472998207D-02  /
   data green(37,18,16) /   2.2650898703433345D-02  /
   data green(37,18,17) /   2.2461454182926562D-02  /
   data green(37,18,18) /   2.2265641088870990D-02  /
   data green(37,19, 0) /   2.4042956972136320D-02  /
   data green(37,19, 1) /   2.4036004068475968D-02  /
   data green(37,19, 2) /   2.4015181599244064D-02  /
   data green(37,19, 3) /   2.3980597763813873D-02  /
   data green(37,19, 4) /   2.3932431156559689D-02  /
   data green(37,19, 5) /   2.3870928211774792D-02  /
   data green(37,19, 6) /   2.3796399727630136D-02  /
   data green(37,19, 7) /   2.3709216564581533D-02  /
   data green(37,19, 8) /   2.3609804633243657D-02  /
   data green(37,19, 9) /   2.3498639301351997D-02  /
   data green(37,19,10) /   2.3376239358642841D-02  /
   data green(37,19,11) /   2.3243160682232242D-02  /
   data green(37,19,12) /   2.3099989743606598D-02  /
   data green(37,19,13) /   2.2947337092152853D-02  /
   data green(37,19,14) /   2.2785830939955190D-02  /
   data green(37,19,15) /   2.2616110959205914D-02  /
   data green(37,19,16) /   2.2438822387921907D-02  /
   data green(37,19,17) /   2.2254610522630128D-02  /
   data green(37,19,18) /   2.2064115659136307D-02  /
   data green(37,19,19) /   2.1867968525178257D-02  /
   data green(37,20, 0) /   2.3776271410370958D-02  /
   data green(37,20, 1) /   2.3769547582270810D-02  /
   data green(37,20, 2) /   2.3749410368182410D-02  /
   data green(37,20, 3) /   2.3715962092566909D-02  /
   data green(37,20, 4) /   2.3669371689166826D-02  /
   data green(37,20, 5) /   2.3609872337285186D-02  /
   data green(37,20, 6) /   2.3537758243999677D-02  /
   data green(37,20, 7) /   2.3453380658724609D-02  /
   data green(37,20, 8) /   2.3357143224450257D-02  /
   data green(37,20, 9) /   2.3249496783457120D-02  /
   data green(37,20,10) /   2.3130933763965692D-02  /
   data green(37,20,11) /   2.3001982277959927D-02  /
   data green(37,20,12) /   2.2863200059511042D-02  /
   data green(37,20,13) /   2.2715168367748482D-02  /
   data green(37,20,14) /   2.2558485969780450D-02  /
   data green(37,20,15) /   2.2393763307088307D-02  /
   data green(37,20,16) /   2.2221616934998759D-02  /
   data green(37,20,17) /   2.2042664309577385D-02  /
   data green(37,20,18) /   2.1857518980449280D-02  /
   data green(37,20,19) /   2.1666786232321052D-02  /
   data green(37,20,20) /   2.1471059202932258D-02  /
   data green(37,21, 0) /   2.3505279517621134D-02  /
   data green(37,21, 1) /   2.3498783224702384D-02  /
   data green(37,21, 2) /   2.3479326702201798D-02  /
   data green(37,21, 3) /   2.3447006570373335D-02  /
   data green(37,21, 4) /   2.3401982380147780D-02  /
   data green(37,21, 5) /   2.3344474431103043D-02  /
   data green(37,21, 6) /   2.3274760799091920D-02  /
   data green(37,21, 7) /   2.3193173651580280D-02  /
   data green(37,21, 8) /   2.3100094945072074D-02  /
   data green(37,21, 9) /   2.2995951611383635D-02  /
   data green(37,21,10) /   2.2881210347645953D-02  /
   data green(37,21,11) /   2.2756372128675994D-02  /
   data green(37,21,12) /   2.2621966559918730D-02  /
   data green(37,21,13) /   2.2478546184873095D-02  /
   data green(37,21,14) /   2.2326680853296102D-02  /
   data green(37,21,15) /   2.2166952246162205D-02  /
   data green(37,21,16) /   2.1999948641032052D-02  /
   data green(37,21,17) /   2.1826259987865113D-02  /
   data green(37,21,18) /   2.1646473351073087D-02  /
   data green(37,21,19) /   2.1461168759370534D-02  /
   data green(37,21,20) /   2.1270915491264579D-02  /
   data green(37,21,21) /   2.1076268811261376D-02  /
   data green(37,22, 0) /   2.3230809614171517D-02  /
   data green(37,22, 1) /   2.3224538458460692D-02  /
   data green(37,22, 2) /   2.3205755497723812D-02  /
   data green(37,22, 3) /   2.3174551838185169D-02  /
   data green(37,22, 4) /   2.3131077957993965D-02  /
   data green(37,22, 5) /   2.3075541696846912D-02  /
   data green(37,22, 6) /   2.3008205515666520D-02  /
   data green(37,22, 7) /   2.2929383096661596D-02  /
   data green(37,22, 8) /   2.2839435368930731D-02  /
   data green(37,22, 9) /   2.2738766056126382D-02  /
   data green(37,22,10) /   2.2627816850276677D-02  /
   data green(37,22,11) /   2.2507062319568800D-02  /
   data green(37,22,12) /   2.2377004657853371D-02  /
   data green(37,22,13) /   2.2238168380124927D-02  /
   data green(37,22,14) /   2.2091095061711953D-02  /
   data green(37,22,15) /   2.1936338209916750D-02  /
   data green(37,22,16) /   2.1774458345982201D-02  /
   data green(37,22,17) /   2.1606018363156163D-02  /
   data green(37,22,18) /   2.1431579213871264D-02  /
   data green(37,22,19) /   2.1251695966216209D-02  /
   data green(37,22,20) /   2.1066914257417264D-02  /
   data green(37,22,21) /   2.0877767160371222D-02  /
   data green(37,22,22) /   2.0684772468668632D-02  /
   data green(37,23, 0) /   2.2953641832629997D-02  /
   data green(37,23, 1) /   2.2947592662848428D-02  /
   data green(37,23, 2) /   2.2929473879166652D-02  /
   data green(37,23, 3) /   2.2899371279016836D-02  /
   data green(37,23, 4) /   2.2857426602783192D-02  /
   data green(37,23, 5) /   2.2803835684858040D-02  /
   data green(37,23, 6) /   2.2738845931741667D-02  /
   data green(37,23, 7) /   2.2662753190406656D-02  /
   data green(37,23, 8) /   2.2575898083596335D-02  /
   data green(37,23, 9) /   2.2478661899114499D-02  /
   data green(37,23,10) /   2.2371462127215767D-02  /
   data green(37,23,11) /   2.2254747743826069D-02  /
   data green(37,23,12) /   2.2128994337601841D-02  /
   data green(37,23,13) /   2.1994699176015493D-02  /
   data green(37,23,14) /   2.1852376300109177D-02  /
   data green(37,23,15) /   2.1702551729756183D-02  /
   data green(37,23,16) /   2.1545758851734154D-02  /
   data green(37,23,17) /   2.1382534052190334D-02  /
   data green(37,23,18) /   2.1213412643700583D-02  /
   data green(37,23,19) /   2.1038925125580588D-02  /
   data green(37,23,20) /   2.0859593804834320D-02  /
   data green(37,23,21) /   2.0675929794476688D-02  /
   data green(37,23,22) /   2.0488430396225523D-02  /
   data green(37,23,23) /   2.0297576865916971D-02  /
   data green(37,24, 0) /   2.2674506516425073D-02  /
   data green(37,24, 1) /   2.2668675528798816D-02  /
   data green(37,24, 2) /   2.2651209583655567D-02  /
   data green(37,24, 3) /   2.2622189385906984D-02  /
   data green(37,24, 4) /   2.2581748291575953D-02  /
   data green(37,24, 5) /   2.2530070610051033D-02  /
   data green(37,24, 6) /   2.2467389286931507D-02  /
   data green(37,24, 7) /   2.2393983024179512D-02  /
   data green(37,24, 8) /   2.2310172906461312D-02  /
   data green(37,24, 9) /   2.2216318612038408D-02  /
   data green(37,24,10) /   2.2112814293109160D-02  /
   data green(37,24,11) /   2.2000084214008508D-02  /
   data green(37,24,12) /   2.1878578236212802D-02  /
   data green(37,24,13) /   2.1748767236865238D-02  /
   data green(37,24,14) /   2.1611138542854177D-02  /
   data green(37,24,15) /   2.1466191455738783D-02  /
   data green(37,24,16) /   2.1314432934478687D-02  /
   data green(37,24,17) /   2.1156373493461372D-02  /
   data green(37,24,18) /   2.0992523363201277D-02  /
   data green(37,24,19) /   2.0823388950744905D-02  /
   data green(37,24,20) /   2.0649469626646547D-02  /
   data green(37,24,21) /   2.0471254855703596D-02  /
   data green(37,24,22) /   2.0289221679719409D-02  /
   data green(37,24,23) /   2.0103832552580164D-02  /
   data green(37,24,24) /   1.9915533521013388D-02  /
   data green(37,25, 0) /   2.2394083442038987D-02  /
   data green(37,25, 1) /   2.2388466276374700D-02  /
   data green(37,25, 2) /   2.2371640164458990D-02  /
   data green(37,25, 3) /   2.2343680942374222D-02  /
   data green(37,25, 4) /   2.2304713947641875D-02  /
   data green(37,25, 5) /   2.2254912462559211D-02  /
   data green(37,25, 6) /   2.2194495588250197D-02  /
   data green(37,25, 7) /   2.2123725600217845D-02  /
   data green(37,25, 8) /   2.2042904847167911D-02  /
   data green(37,25, 9) /   2.1952372263503828D-02  /
   data green(37,25,10) /   2.1852499571939247D-02  /
   data green(37,25,11) /   2.1743687256047636D-02  /
   data green(37,25,12) /   2.1626360383310984D-02  /
   data green(37,25,13) /   2.1500964357503205D-02  /
   data green(37,25,14) /   2.1367960675317182D-02  /
   data green(37,25,15) /   2.1227822756351881D-02  /
   data green(37,25,16) /   2.1081031908314837D-02  /
   data green(37,25,17) /   2.0928073480972299D-02  /
   data green(37,25,18) /   2.0769433253409436D-02  /
   data green(37,25,19) /   2.0605594089930346D-02  /
   data green(37,25,20) /   2.0437032890782972D-02  /
   data green(37,25,21) /   2.0264217855130034D-02  /
   data green(37,25,22) /   2.0087606065544494D-02  /
   data green(37,25,23) /   1.9907641395962416D-02  /
   data green(37,25,24) /   1.9724752738596665D-02  /
   data green(37,25,25) /   1.9539352539866401D-02  /
   data green(37,26, 0) /   2.2113001759201643D-02  /
   data green(37,26, 1) /   2.2107593589292155D-02  /
   data green(37,26, 2) /   2.2091392908643041D-02  /
   data green(37,26, 3) /   2.2064470912424691D-02  /
   data green(37,26, 4) /   2.2026945292725784D-02  /
   data green(37,26, 5) /   2.1978978813086410D-02  /
   data green(37,26, 6) /   2.1920777360503043D-02  /
   data green(37,26, 7) /   2.1852587520308178D-02  /
   data green(37,26, 8) /   2.1774693729224141D-02  /
   data green(37,26, 9) /   2.1687415069730872D-02  /
   data green(37,26,10) /   2.1591101774463252D-02  /
   data green(37,26,11) /   2.1486131512575397D-02  /
   data green(37,26,12) /   2.1372905530906892D-02  /
   data green(37,26,13) /   2.1251844721487816D-02  /
   data green(37,26,14) /   2.1123385683650926D-02  /
   data green(37,26,15) /   2.0987976844064032D-02  /
   data green(37,26,16) /   2.0846074691694243D-02  /
   data green(37,26,17) /   2.0698140177420530D-02  /
   data green(37,26,18) /   2.0544635320084308D-02  /
   data green(37,26,19) /   2.0386020052549435D-02  /
   data green(37,26,20) /   2.0222749333143799D-02  /
   data green(37,26,21) /   2.0055270539941585D-02  /
   data green(37,26,22) /   1.9884021157936544D-02  /
   data green(37,26,23) /   1.9709426762418474D-02  /
   data green(37,26,24) /   1.9531899295912883D-02  /
   data green(37,26,25) /   1.9351835630945585D-02  /
   data green(37,26,26) /   1.9169616406673674D-02  /
   data green(37,27, 0) /   2.1831840544133910D-02  /
   data green(37,27, 1) /   2.1826636161843615D-02  /
   data green(37,27, 2) /   2.1811045365174614D-02  /
   data green(37,27, 3) /   2.1785134937750923D-02  /
   data green(37,27, 4) /   2.1749015301935980D-02  /
   data green(37,27, 5) /   2.1702839214993823D-02  /
   data green(37,27, 6) /   2.1646799986204269D-02  /
   data green(37,27, 7) /   2.1581129255460837D-02  /
   data green(37,27, 8) /   2.1506094382785019D-02  /
   data green(37,27, 9) /   2.1421995505299253D-02  /
   data green(37,27,10) /   2.1329162323328776D-02  /
   data green(37,27,11) /   2.1227950680362939D-02  /
   data green(37,27,12) /   2.1118739002615423D-02  /
   data green(37,27,13) /   2.1001924662985810D-02  /
   data green(37,27,14) /   2.0877920331524989D-02  /
   data green(37,27,15) /   2.0747150370288767D-02  /
   data green(37,27,16) /   2.0610047325014598D-02  /
   data green(37,27,17) /   2.0467048559682333D-02  /
   data green(37,27,18) /   2.0318593073036285D-02  /
   data green(37,27,19) /   2.0165118528850442D-02  /
   data green(37,27,20) /   2.0007058524389758D-02  /
   data green(37,27,21) /   1.9844840114394687D-02  /
   data green(37,27,22) /   1.9678881601196914D-02  /
   data green(37,27,23) /   1.9509590595412317D-02  /
   data green(37,27,24) /   1.9337362346166603D-02  /
   data green(37,27,25) /   1.9162578335057062D-02  /
   data green(37,27,26) /   1.8985605124072262D-02  /
   data green(37,27,27) /   1.8806793444478789D-02  /
   data green(37,28, 0) /   2.1551129864752112D-02  /
   data green(37,28, 1) /   2.1546123757383049D-02  /
   data green(37,28, 2) /   2.1531126383383656D-02  /
   data green(37,28, 3) /   2.1506200343268875D-02  /
   data green(37,28, 4) /   2.1471449164107384D-02  /
   data green(37,28, 5) /   2.1427016108127247D-02  /
   data green(37,28, 6) /   2.1373082542596813D-02  /
   data green(37,28, 7) /   2.1309865907111990D-02  /
   data green(37,28, 8) /   2.1237617322421257D-02  /
   data green(37,28, 9) /   2.1156618891354867D-02  /
   data green(37,28,10) /   2.1067180747128159D-02  /
   data green(37,28,11) /   2.0969637907181229D-02  /
   data green(37,28,12) /   2.0864346991801744D-02  /
   data green(37,28,13) /   2.0751682866139827D-02  /
   data green(37,28,14) /   2.0632035262014927D-02  /
   data green(37,28,15) /   2.0505805432339666D-02  /
   data green(37,28,16) /   2.0373402886290571D-02  /
   data green(37,28,17) /   2.0235242247802963D-02  /
   data green(37,28,18) /   2.0091740273831923D-02  /
   data green(37,28,19) /   1.9943313062363400D-02  /
   data green(37,28,20) /   1.9790373473624816D-02  /
   data green(37,28,21) /   1.9633328781547540D-02  /
   data green(37,28,22) /   1.9472578566455429D-02  /
   data green(37,28,23) /   1.9308512854337160D-02  /
   data green(37,28,24) /   1.9141510503011487D-02  /
   data green(37,28,25) /   1.8971937831082520D-02  /
   data green(37,28,26) /   1.8800147481841778D-02  /
   data green(37,28,27) /   1.8626477511212668D-02  /
   data green(37,28,28) /   1.8451250686431879D-02  /
   data green(37,29, 0) /   2.1271352262764553D-02  /
   data green(37,29, 1) /   2.1266538683517539D-02  /
   data green(37,29, 2) /   2.1252117567567105D-02  /
   data green(37,29, 3) /   2.1228147557830052D-02  /
   data green(37,29, 4) /   2.1194725655935671D-02  /
   data green(37,29, 5) /   2.1151986134525857D-02  /
   data green(37,29, 6) /   2.1100099048124321D-02  /
   data green(37,29, 7) /   2.1039268374750009D-02  /
   data green(37,29, 8) /   2.0969729827626654D-02  /
   data green(37,29, 9) /   2.0891748382156950D-02  /
   data green(37,29,10) /   2.0805615567634633D-02  /
   data green(37,29,11) /   2.0711646575887949D-02  /
   data green(37,29,12) /   2.0610177240177999D-02  /
   data green(37,29,13) /   2.0501560937283926D-02  /
   data green(37,29,14) /   2.0386165463917819D-02  /
   data green(37,29,15) /   2.0264369935596813D-02  /
   data green(37,29,16) /   2.0136561752069199D-02  /
   data green(37,29,17) /   2.0003133668568186D-02  /
   data green(37,29,18) /   1.9864481006791602D-02  /
   data green(37,29,19) /   1.9720999033803463D-02  /
   data green(37,29,20) /   1.9573080531243053D-02  /
   data green(37,29,21) /   1.9421113571498158D-02  /
   data green(37,29,22) /   1.9265479512016556D-02  /
   data green(37,29,23) /   1.9106551213825419D-02  /
   data green(37,29,24) /   1.8944691485700431D-02  /
   data green(37,29,25) /   1.8780251751345130D-02  /
   data green(37,29,26) /   1.8613570933441938D-02  /
   data green(37,29,27) /   1.8444974545534554D-02  /
   data green(37,29,28) /   1.8274773980385545D-02  /
   data green(37,29,29) /   1.8103265981696617D-02  /
   data green(37,30, 0) /   2.0992944565124780D-02  /
   data green(37,30, 1) /   2.0988317596646711D-02  /
   data green(37,30, 2) /   2.0974455060905152D-02  /
   data green(37,30, 3) /   2.0951411864155153D-02  /
   data green(37,30, 4) /   2.0919278845142018D-02  /
   data green(37,30, 5) /   2.0878181782811544D-02  /
   data green(37,30, 6) /   2.0828280037003825D-02  /
   data green(37,30, 7) /   2.0769764850754980D-02  /
   data green(37,30, 8) /   2.0702857349262017D-02  /
   data green(37,30, 9) /   2.0627806275815789D-02  /
   data green(37,30,10) /   2.0544885508939143D-02  /
   data green(37,30,11) /   2.0454391407514136D-02  /
   data green(37,30,12) /   2.0356640031833583D-02  /
   data green(37,30,13) /   2.0251964288320524D-02  /
   data green(37,30,14) /   2.0140711044226677D-02  /
   data green(37,30,15) /   2.0023238256091560D-02  /
   data green(37,30,16) /   1.9899912152294103D-02  /
   data green(37,30,17) /   1.9771104505852143D-02  /
   data green(37,30,18) /   1.9637190028926960D-02  /
   data green(37,30,19) /   1.9498543915467726D-02  /
   data green(37,30,20) /   1.9355539553276797D-02  /
   data green(37,30,21) /   1.9208546421658496D-02  /
   data green(37,30,22) /   1.9057928185882509D-02  /
   data green(37,30,23) /   1.8904040995065222D-02  /
   data green(37,30,24) /   1.8747231985844597D-02  /
   data green(37,30,25) /   1.8587837990460054D-02  /
   data green(37,30,26) /   1.8426184444589516D-02  /
   data green(37,30,27) /   1.8262584487558148D-02  /
   data green(37,30,28) /   1.8097338245313333D-02  /
   data green(37,30,29) /   1.7930732284841218D-02  /
   data green(37,30,30) /   1.7763039227447409D-02  /
   data green(37,31, 0) /   2.0716299945766434D-02  /
   data green(37,31, 1) /   2.0711853556917469D-02  /
   data green(37,31, 2) /   2.0698531580182380D-02  /
   data green(37,31, 3) /   2.0676385400187621D-02  /
   data green(37,31, 4) /   2.0645500045847615D-02  /
   data green(37,31, 5) /   2.0605993285685804D-02  /
   data green(37,31, 6) /   2.0558014387830774D-02  /
   data green(37,31, 7) /   2.0501742570131678D-02  /
   data green(37,31, 8) /   2.0437385171594241D-02  /
   data green(37,31, 9) /   2.0365175581070399D-02  /
   data green(37,31,10) /   2.0285370962719665D-02  /
   data green(37,31,11) /   2.0198249820135818D-02  /
   data green(37,31,12) /   2.0104109442184527D-02  /
   data green(37,31,13) /   2.0003263273566513D-02  /
   data green(37,31,14) /   1.9896038251989809D-02  /
   data green(37,31,15) /   1.9782772151723263D-02  /
   data green(37,31,16) /   1.9663810970362249D-02  /
   data green(37,31,17) /   1.9539506392030281D-02  /
   data green(37,31,18) /   1.9410213356143704D-02  /
   data green(37,31,19) /   1.9276287756454364D-02  /
   data green(37,31,20) /   1.9138084290522949D-02  /
   data green(37,31,21) /   1.8995954475214856D-02  /
   data green(37,31,22) /   1.8850244839384179D-02  /
   data green(37,31,23) /   1.8701295300727866D-02  /
   data green(37,31,24) /   1.8549437729939916D-02  /
   data green(37,31,25) /   1.8394994701835587D-02  /
   data green(37,31,26) /   1.8238278430090786D-02  /
   data green(37,31,27) /   1.8079589879669976D-02  /
   data green(37,31,28) /   1.7919218048901654D-02  /
   data green(37,31,29) /   1.7757439411489745D-02  /
   data green(37,31,30) /   1.7594517507498995D-02  /
   data green(37,31,31) /   1.7430702671488791D-02  /
   data green(37,32, 0) /   2.0441770167444549D-02  /
   data green(37,32, 1) /   2.0437498263528909D-02  /
   data green(37,32, 2) /   2.0424698631580213D-02  /
   data green(37,32, 3) /   2.0403419342683641D-02  /
   data green(37,32, 4) /   2.0373739957740732D-02  /
   data green(37,32, 5) /   2.0335770703094991D-02  /
   data green(37,32, 6) /   2.0289651339962203D-02  /
   data green(37,32, 7) /   2.0235549750267149D-02  /
   data green(37,32, 8) /   2.0173660266638017D-02  /
   data green(37,32, 9) /   2.0104201778569780D-02  /
   data green(37,32,10) /   2.0027415650030817D-02  /
   data green(37,32,11) /   1.9943563485995131D-02  /
   data green(37,32,12) /   1.9852924786518990D-02  /
   data green(37,32,13) /   1.9755794527076771D-02  /
   data green(37,32,14) /   1.9652480702992452D-02  /
   data green(37,32,15) /   1.9543301874051767D-02  /
   data green(37,32,16) /   1.9428584742879833D-02  /
   data green(37,32,17) /   1.9308661797562552D-02  /
   data green(37,32,18) /   1.9183869045426870D-02  /
   data green(37,32,19) /   1.9054543861026129D-02  /
   data green(37,32,20) /   1.8921022967348312D-02  /
   data green(37,32,21) /   1.8783640565208818D-02  /
   data green(37,32,22) /   1.8642726621826491D-02  /
   data green(37,32,23) /   1.8498605325808676D-02  /
   data green(37,32,24) /   1.8351593712270436D-02  /
   data green(37,32,25) /   1.8202000458642793D-02  /
   data green(37,32,26) /   1.8050124848926839D-02  /
   data green(37,32,27) /   1.7896255901745404D-02  /
   data green(37,32,28) /   1.7740671655540077D-02  /
   data green(37,32,29) /   1.7583638602650623D-02  /
   data green(37,32,30) /   1.7425411262778970D-02  /
   data green(37,32,31) /   1.7266231885454168D-02  /
   data green(37,32,32) /   1.7106330270547521D-02  /
   data green(37,33, 0) /   2.0169667942458669D-02  /
   data green(37,33, 1) /   2.0165564409245416D-02  /
   data green(37,33, 2) /   2.0153268846644763D-02  /
   data green(37,33, 3) /   2.0132826212552985D-02  /
   data green(37,33, 4) /   2.0104310929124906D-02  /
   data green(37,33, 5) /   2.0067826131889854D-02  /
   data green(37,33, 6) /   2.0023502639402436D-02  /
   data green(37,33, 7) /   1.9971497663493838D-02  /
   data green(37,33, 8) /   1.9911993284791984D-02  /
   data green(37,33, 9) /   1.9845194722009993D-02  /
   data green(37,33,10) /   1.9771328426468476D-02  /
   data green(37,33,11) /   1.9690640035362091D-02  /
   data green(37,33,12) /   1.9603392218390459D-02  /
   data green(37,33,13) /   1.9509862452566713D-02  /
   data green(37,33,14) /   1.9410340759350237D-02  /
   data green(37,33,15) /   1.9305127436805736D-02  /
   data green(37,33,16) /   1.9194530817373594D-02  /
   data green(37,33,17) /   1.9078865079167798D-02  /
   data green(37,33,18) /   1.8958448135625128D-02  /
   data green(37,33,19) /   1.8833599624944786D-02  /
   data green(37,33,20) /   1.8704639017206929D-02  /
   data green(37,33,21) /   1.8571883853459537D-02  /
   data green(37,33,22) /   1.8435648127522217D-02  /
   data green(37,33,23) /   1.8296240817860659D-02  /
   data green(37,33,24) /   1.8153964573712753D-02  /
   data green(37,33,25) /   1.8009114556750472D-02  /
   data green(37,33,26) /   1.7861977436981678D-02  /
   data green(37,33,27) /   1.7712830539355795D-02  /
   data green(37,33,28) /   1.7561941135647269D-02  /
   data green(37,33,29) /   1.7409565874647467D-02  /
   data green(37,33,30) /   1.7255950342487794D-02  /
   data green(37,33,31) /   1.7101328744023238D-02  /
   data green(37,33,32) /   1.6945923695600294D-02  /
   data green(37,33,33) /   1.6789946119187143D-02  /
   data green(37,34, 0) /   1.9900269359737245D-02  /
   data green(37,34, 1) /   1.9896328101654735D-02  /
   data green(37,34, 2) /   1.9884518386140179D-02  /
   data green(37,34, 3) /   1.9864882249950495D-02  /
   data green(37,34, 4) /   1.9837489292255175D-02  /
   data green(37,34, 5) /   1.9802435990912005D-02  /
   data green(37,34, 6) /   1.9759844763769650D-02  /
   data green(37,34, 7) /   1.9709862792803571D-02  /
   data green(37,34, 8) /   1.9652660633003069D-02  /
   data green(37,34, 9) /   1.9588430631369714D-02  /
   data green(37,34,10) /   1.9517385184079090D-02  /
   data green(37,34,11) /   1.9439754861747763D-02  /
   data green(37,34,12) /   1.9355786433818706D-02  /
   data green(37,34,13) /   1.9265740823346634D-02  /
   data green(37,34,14) /   1.9169891022972153D-02  /
   data green(37,34,15) /   1.9068520001691281D-02  /
   data green(37,34,16) /   1.8961918630240163D-02  /
   data green(37,34,17) /   1.8850383650628071D-02  /
   data green(37,34,18) /   1.8734215712674147D-02  /
   data green(37,34,19) /   1.8613717497447640D-02  /
   data green(37,34,20) /   1.8489191944388160D-02  /
   data green(37,34,21) /   1.8360940595695360D-02  /
   data green(37,34,22) /   1.8229262068419842D-02  /
   data green(37,34,23) /   1.8094450661639554D-02  /
   data green(37,34,24) /   1.7956795103236858D-02  /
   data green(37,34,25) /   1.7816577438151383D-02  /
   data green(37,34,26) /   1.7674072057612532D-02  /
   data green(37,34,27) /   1.7529544866776176D-02  /
   data green(37,34,28) /   1.7383252586415335D-02  /
   data green(37,34,29) /   1.7235442182845041D-02  /
   data green(37,34,30) /   1.7086350419089708D-02  /
   data green(37,34,31) /   1.6936203519410876D-02  /
   data green(37,34,32) /   1.6785216938684407D-02  /
   data green(37,34,33) /   1.6633595227723386D-02  /
   data green(37,34,34) /   1.6481531985459374D-02  /
   data green(37,35, 0) /   1.9633816334004135D-02  /
   data green(37,35, 1) /   1.9630031306928752D-02  /
   data green(37,35, 2) /   1.9618689367656199D-02  /
   data green(37,35, 3) /   1.9599829815173241D-02  /
   data green(37,35, 4) /   1.9573517727281588D-02  /
   data green(37,35, 5) /   1.9539843338172234D-02  /
   data green(37,35, 6) /   1.9498921183438739D-02  /
   data green(37,35, 7) /   1.9450889028326226D-02  /
   data green(37,35, 8) /   1.9395906598682656D-02  /
   data green(37,35, 9) /   1.9334154137170854D-02  /
   data green(37,35,10) /   1.9265830809738829D-02  /
   data green(37,35,11) /   1.9191152989088005D-02  /
   data green(37,35,12) /   1.9110352442906209D-02  /
   data green(37,35,13) /   1.9023674454953995D-02  /
   data green(37,35,14) /   1.8931375906745061D-02  /
   data green(37,35,15) /   1.8833723346600051D-02  /
   data green(37,35,16) /   1.8730991071351514D-02  /
   data green(37,35,17) /   1.8623459244023477D-02  /
   data green(37,35,18) /   1.8511412068494739D-02  /
   data green(37,35,19) /   1.8395136039579624D-02  /
   data green(37,35,20) /   1.8274918284217324D-02  /
   data green(37,35,21) /   1.8151045006643558D-02  /
   data green(37,35,22) /   1.8023800047608257D-02  /
   data green(37,35,23) /   1.7893463564972577D-02  /
   data green(37,35,24) /   1.7760310840429606D-02  /
   data green(37,35,25) /   1.7624611214693865D-02  /
   data green(37,35,26) /   1.7486627151330376D-02  /
   data green(37,35,27) /   1.7346613427471791D-02  /
   data green(37,35,28) /   1.7204816448011324D-02  /
   data green(37,35,29) /   1.7061473678468123D-02  /
   data green(37,35,30) /   1.6916813190592100D-02  /
   data green(37,35,31) /   1.6771053313898526D-02  /
   data green(37,35,32) /   1.6624402385680470D-02  /
   data green(37,35,33) /   1.6477058591620998D-02  /
   data green(37,35,34) /   1.6329209888892655D-02  /
   data green(37,35,35) /   1.6181034003565797D-02  /
   data green(37,36, 0) /   1.9370519040380119D-02  /
   data green(37,36, 1) /   1.9366884279459262D-02  /
   data green(37,36, 2) /   1.9355992280401045D-02  /
   data green(37,36, 3) /   1.9337879778895603D-02  /
   data green(37,36, 4) /   1.9312607618478632D-02  /
   data green(37,36, 5) /   1.9280260183994532D-02  /
   data green(37,36, 6) /   1.9240944622983924D-02  /
   data green(37,36, 7) /   1.9194789869005982D-02  /
   data green(37,36, 8) /   1.9141945484181428D-02  /
   data green(37,36, 9) /   1.9082580341015946D-02  /
   data green(37,36,10) /   1.9016881165772690D-02  /
   data green(37,36,11) /   1.8945050967264239D-02  /
   data green(37,36,12) /   1.8867307375912069D-02  /
   data green(37,36,13) /   1.8783880918281953D-02  /
   data green(37,36,14) /   1.8695013252072221D-02  /
   data green(37,36,15) /   1.8600955385758393D-02  /
   data green(37,36,16) /   1.8501965905839780D-02  /
   data green(37,36,17) /   1.8398309232968115D-02  /
   data green(37,36,18) /   1.8290253926242398D-02  /
   data green(37,36,19) /   1.8178071052713452D-02  /
   data green(37,36,20) /   1.8062032636737967D-02  /
   data green(37,36,21) /   1.7942410201335482D-02  /
   data green(37,36,22) /   1.7819473411205183D-02  /
   data green(37,36,23) /   1.7693488824618405D-02  /
   data green(37,36,24) /   1.7564718759070914D-02  /
   data green(37,36,25) /   1.7433420273403687D-02  /
   data green(37,36,26) /   1.7299844267113281D-02  /
   data green(37,36,27) /   1.7164234695799818D-02  /
   data green(37,36,28) /   1.7026827900154310D-02  /
   data green(37,36,29) /   1.6887852044574849D-02  /
   data green(37,36,30) /   1.6747526660420436D-02  /
   data green(37,36,31) /   1.6606062288055164D-02  /
   data green(37,36,32) /   1.6463660211190431D-02  /
   data green(37,36,33) /   1.6320512276582804D-02  /
   data green(37,36,34) /   1.6176800791872125D-02  /
   data green(37,36,35) /   1.6032698494225998D-02  /
   data green(37,36,36) /   1.5888368582474430D-02  /
   data green(37,37, 0) /   1.9110558304707377D-02  /
   data green(37,37, 1) /   1.9107067947662364D-02  /
   data green(37,37, 2) /   1.9096608357490923D-02  /
   data green(37,37, 3) /   1.9079213872084968D-02  /
   data green(37,37, 4) /   1.9054941373152764D-02  /
   data green(37,37, 5) /   1.9023869770591315D-02  /
   data green(37,37, 6) /   1.8986099293488570D-02  /
   data green(37,37, 7) /   1.8941750600180385D-02  /
   data green(37,37, 8) /   1.8890963722708969D-02  /
   data green(37,37, 9) /   1.8833896863517715D-02  /
   data green(37,37,10) /   1.8770725064214743D-02  /
   data green(37,37,11) /   1.8701638767708355D-02  /
   data green(37,37,12) /   1.8626842295942074D-02  /
   data green(37,37,13) /   1.8546552265842106D-02  /
   data green(37,37,14) /   1.8460995965953689D-02  /
   data green(37,37,15) /   1.8370409715626108D-02  /
   data green(37,37,16) /   1.8275037227557776D-02  /
   data green(37,37,17) /   1.8175127993096257D-02  /
   data green(37,37,18) /   1.8070935707970838D-02  /
   data green(37,37,19) /   1.7962716754188784D-02  /
   data green(37,37,20) /   1.7850728751723496D-02  /
   data green(37,37,21) /   1.7735229191431636D-02  /
   data green(37,37,22) /   1.7616474158422674D-02  /
   data green(37,37,23) /   1.7494717152924805D-02  /
   data green(37,37,24) /   1.7370208013596719D-02  /
   data green(37,37,25) /   1.7243191946264766D-02  /
   data green(37,37,26) /   1.7113908659254152D-02  /
   data green(37,37,27) /   1.6982591604851324D-02  /
   data green(37,37,28) /   1.6849467325000114D-02  /
   data green(37,37,29) /   1.6714754898101994D-02  /
   data green(37,37,30) /   1.6578665482762527D-02  /
   data green(37,37,31) /   1.6441401953497352D-02  /
   data green(37,37,32) /   1.6303158622769969D-02  /
   data green(37,37,33) /   1.6164121043270637D-02  /
   data green(37,37,34) /   1.6024465884042157D-02  /
   data green(37,37,35) /   1.5884360873898527D-02  /
   data green(37,37,36) /   1.5743964805549095D-02  /
   data green(37,37,37) /   1.5603427593913725D-02  /
   data green(38, 0, 0) /   2.6320354218392742D-02  /
   data green(38, 1, 0) /   2.6311226258614571D-02  /
   data green(38, 1, 1) /   2.6302107822160914D-02  /
   data green(38, 2, 0) /   2.6283899519351780D-02  /
   data green(38, 2, 1) /   2.6274809553094487D-02  /
   data green(38, 2, 2) /   2.6247596396757531D-02  /
   data green(38, 3, 0) /   2.6238544423039904D-02  /
   data green(38, 3, 1) /   2.6229501576972441D-02  /
   data green(38, 3, 2) /   2.6202429289631474D-02  /
   data green(38, 3, 3) /   2.6157495337619883D-02  /
   data green(38, 4, 0) /   2.6175441715483277D-02  /
   data green(38, 4, 1) /   2.6166464152869266D-02  /
   data green(38, 4, 2) /   2.6139587037728677D-02  /
   data green(38, 4, 3) /   2.6094976128712331D-02  /
   data green(38, 4, 4) /   2.6032904524705042D-02  /
   data green(38, 5, 0) /   2.6094977646706017D-02  /
   data green(38, 5, 1) /   2.6086082865954048D-02  /
   data green(38, 5, 2) /   2.6059453241317553D-02  /
   data green(38, 5, 3) /   2.6015251986921638D-02  /
   data green(38, 5, 4) /   2.5953748028932525D-02  /
   data green(38, 5, 5) /   2.5875311468779565D-02  /
   data green(38, 6, 0) /   2.5997637453694673D-02  /
   data green(38, 6, 1) /   2.5988842125927003D-02  /
   data green(38, 6, 2) /   2.5962509840147606D-02  /
   data green(38, 6, 3) /   2.5918800774737214D-02  /
   data green(38, 6, 4) /   2.5857978877027392D-02  /
   data green(38, 6, 5) /   2.5780407443960884D-02  /
   data green(38, 6, 6) /   2.5686543140227391D-02  /
   data green(38, 7, 0) /   2.5883997359051551D-02  /
   data green(38, 7, 1) /   2.5875317184671696D-02  /
   data green(38, 7, 2) /   2.5849329187748186D-02  /
   data green(38, 7, 3) /   2.5806190060681133D-02  /
   data green(38, 7, 4) /   2.5746158032767974D-02  /
   data green(38, 7, 5) /   2.5669588584619551D-02  /
   data green(38, 7, 6) /   2.5576928645555861D-02  /
   data green(38, 7, 7) /   2.5468709458309275D-02  /
   data green(38, 8, 0) /   2.5754715340011794D-02  /
   data green(38, 8, 1) /   2.5746164926803334D-02  /
   data green(38, 8, 2) /   2.5720564906180597D-02  /
   data green(38, 8, 3) /   2.5678068079783324D-02  /
   data green(38, 8, 4) /   2.5618926293919507D-02  /
   data green(38, 8, 5) /   2.5543486301863896D-02  /
   data green(38, 8, 6) /   2.5452184159532060D-02  /
   data green(38, 8, 7) /   2.5345538330802605D-02  /
   data green(38, 8, 8) /   2.5224141712820343D-02  /
   data green(38, 9, 0) /   2.5610520949567794D-02  /
   data green(38, 9, 1) /   2.5602113714022053D-02  /
   data green(38, 9, 2) /   2.5576941799667980D-02  /
   data green(38, 9, 3) /   2.5535153761318111D-02  /
   data green(38, 9, 4) /   2.5476994476928944D-02  /
   data green(38, 9, 5) /   2.5402801169599486D-02  /
   data green(38, 9, 6) /   2.5312998017474095D-02  /
   data green(38, 9, 7) /   2.5208089519210112D-02  /
   data green(38, 9, 8) /   2.5088652815247758D-02  /
   data green(38, 9, 9) /   2.4955329187874391D-02  /
   data green(38,10, 0) /   2.5452204484981290D-02  /
   data green(38,10, 1) /   2.5443952577704796D-02  /
   data green(38,10, 2) /   2.5419245119012536D-02  /
   data green(38,10, 3) /   2.5378226111565948D-02  /
   data green(38,10, 4) /   2.5321132962552392D-02  /
   data green(38,10, 5) /   2.5248292674844956D-02  /
   data green(38,10, 6) /   2.5160116683947825D-02  /
   data green(38,10, 7) /   2.5057094499368614D-02  /
   data green(38,10, 8) /   2.4939786340067591D-02  /
   data green(38,10, 9) /   2.4808814975450649D-02  /
   data green(38,10,10) /   2.4664856995449857D-02  /
   data green(38,11, 0) /   2.5280605799095576D-02  /
   data green(38,11, 1) /   2.5272520054334204D-02  /
   data green(38,11, 2) /   2.5248309468973301D-02  /
   data green(38,11, 3) /   2.5208113239784002D-02  /
   data green(38,11, 4) /   2.5152160885126108D-02  /
   data green(38,11, 5) /   2.5080768612390003D-02  /
   data green(38,11, 6) /   2.4994334391679932D-02  /
   data green(38,11, 7) /   2.4893331885095259D-02  /
   data green(38,11, 8) /   2.4778303410350148D-02  /
   data green(38,11, 9) /   2.4649852138299060D-02  /
   data green(38,11,10) /   2.4508633735681015D-02  /
   data green(38,11,11) /   2.4355347667131686D-02  /
   data green(38,12, 0) /   2.5096603037785097D-02  /
   data green(38,12, 1) /   2.5088692946358820D-02  /
   data green(38,12, 2) /   2.5065007638986052D-02  /
   data green(38,12, 3) /   2.5025681304125864D-02  /
   data green(38,12, 4) /   2.4970935237218411D-02  /
   data green(38,12, 5) /   2.4901074389262909D-02  /
   data green(38,12, 6) /   2.4816482683907085D-02  /
   data green(38,12, 7) /   2.4717617242992873D-02  /
   data green(38,12, 8) /   2.4605001688222543D-02  /
   data green(38,12, 9) /   2.4479218706410431D-02  /
   data green(38,12,10) /   2.4340902077165643D-02  /
   data green(38,12,11) /   2.4190728364850679D-02  /
   data green(38,12,12) /   2.4029408471816734D-02  /
   data green(38,13, 0) /   2.4901101564441575D-02  /
   data green(38,13, 1) /   2.4893375268744772D-02  /
   data green(38,13, 2) /   2.4870239615574757D-02  /
   data green(38,13, 3) /   2.4831823632712788D-02  /
   data green(38,13, 4) /   2.4778340140536651D-02  /
   data green(38,13, 5) /   2.4710082484442258D-02  /
   data green(38,13, 6) /   2.4627420099139360D-02  /
   data green(38,13, 7) /   2.4530793035349868D-02  /
   data green(38,13, 8) /   2.4420705605476045D-02  /
   data green(38,13, 9) /   2.4297719323556229D-02  /
   data green(38,13,10) /   2.4162445325809138D-02  /
   data green(38,13,11) /   2.4015536461296390D-02  /
   data green(38,13,12) /   2.3857679238173776D-02  /
   data green(38,13,13) /   2.3689585800499093D-02  /
   data green(38,14, 0) /   2.4695023301887394D-02  /
   data green(38,14, 1) /   2.4687487611088962D-02  /
   data green(38,14, 2) /   2.4664922004777275D-02  /
   data green(38,14, 3) /   2.4627450245626566D-02  /
   data green(38,14, 4) /   2.4575276505338375D-02  /
   data green(38,14, 5) /   2.4508682281629206D-02  /
   data green(38,14, 6) /   2.4428022210872837D-02  /
   data green(38,14, 7) /   2.4333718897617488D-02  /
   data green(38,14, 8) /   2.4226256906557331D-02  /
   data green(38,14, 9) /   2.4106176080224414D-02  /
   data green(38,14,10) /   2.3974064356227495D-02  /
   data green(38,14,11) /   2.3830550261281395D-02  /
   data green(38,14,12) /   2.3676295255952856D-02  /
   data green(38,14,13) /   2.3511986094738285D-02  /
   data green(38,14,14) /   2.3338327351794377D-02  /
   data green(38,15, 0) /   2.4479296686009459D-02  /
   data green(38,15, 1) /   2.4471957109201924D-02  /
   data green(38,15, 2) /   2.4449978057336924D-02  /
   data green(38,15, 3) /   2.4413477968663506D-02  /
   data green(38,15, 4) /   2.4362652266544830D-02  /
   data green(38,15, 5) /   2.4297770459937786D-02  /
   data green(38,15, 6) /   2.4219172203111105D-02  /
   data green(38,15, 7) /   2.4127262426717963D-02  /
   data green(38,15, 8) /   2.4022505675037369D-02  /
   data green(38,15, 9) /   2.3905419800820261D-02  /
   data green(38,15,10) /   2.3776569179288349D-02  /
   data green(38,15,11) /   2.3636557606400949D-02  /
   data green(38,15,12) /   2.3486021043873189D-02  /
   data green(38,15,13) /   2.3325620365251338D-02  /
   data green(38,15,14) /   2.3156034244526728D-02  /
   data green(38,15,15) /   2.2977952312341725D-02  /
   data green(38,16, 0) /   2.4254847386208610D-02  /
   data green(38,16, 1) /   2.4247708180999722D-02  /
   data green(38,16, 2) /   2.4226328450737668D-02  /
   data green(38,16, 3) /   2.4190821291664673D-02  /
   data green(38,16, 4) /   2.4141373347615663D-02  /
   data green(38,16, 5) /   2.4078242091326838D-02  /
   data green(38,16, 6) /   2.4001752127822139D-02  /
   data green(38,16, 7) /   2.3912290623172686D-02  /
   data green(38,16, 8) /   2.3810301982998126D-02  /
   data green(38,16, 9) /   2.3696281920641802D-02  /
   data green(38,16,10) /   2.3570771064597731D-02  /
   data green(38,16,11) /   2.3434348258441876D-02  /
   data green(38,16,12) /   2.3287623704515325D-02  /
   data green(38,16,13) /   2.3131232095492005D-02  /
   data green(38,16,14) /   2.2965825866533889D-02  /
   data green(38,16,15) /   2.2792068685924635D-02  /
   data green(38,16,16) /   2.2610629284881674D-02  /
   data green(38,17, 0) /   2.4022589907753682D-02  /
   data green(38,17, 1) /   2.4015654141665132D-02  /
   data green(38,17, 2) /   2.3994882942654963D-02  /
   data green(38,17, 3) /   2.3960384085344417D-02  /
   data green(38,17, 4) /   2.3912335465191843D-02  /
   data green(38,17, 5) /   2.3850982556596942D-02  /
   data green(38,17, 6) /   2.3776634954698582D-02  /
   data green(38,17, 7) /   2.3689662095690314D-02  /
   data green(38,17, 8) /   2.3590488269974338D-02  /
   data green(38,17, 9) /   2.3479587056996189D-02  /
   data green(38,17,10) /   2.3357475319767406D-02  /
   data green(38,17,11) /   2.3224706900825610D-02  /
   data green(38,17,12) /   2.3081866159939273D-02  /
   data green(38,17,13) /   2.2929561487735159D-02  /
   data green(38,17,14) /   2.2768418919306459D-02  /
   data green(38,17,15) /   2.2599075958579139D-02  /
   data green(38,17,16) /   2.2422175708668467D-02  /
   data green(38,17,17) /   2.2238361386545887D-02  /
   data green(38,18, 0) /   2.3783420152297952D-02  /
   data green(38,18, 1) /   2.3776689774319328D-02  /
   data green(38,18, 2) /   2.3756532972009574D-02  /
   data green(38,18, 3) /   2.3723052252704963D-02  /
   data green(38,18, 4) /   2.3676416850436956D-02  /
   data green(38,18, 5) /   2.3616860355647831D-02  /
   data green(38,18, 6) /   2.3544677488584829D-02  /
   data green(38,18, 7) /   2.3460220103126438D-02  /
   data green(38,18, 8) /   2.3363892525771559D-02  /
   data green(38,18, 9) /   2.3256146348027237D-02  /
   data green(38,18,10) /   2.3137474799106425D-02  /
   data green(38,18,11) /   2.3008406829615410D-02  /
   data green(38,18,12) /   2.2869501035969481D-02  /
   data green(38,18,13) /   2.2721339550048257D-02  /
   data green(38,18,14) /   2.2564522009700931D-02  /
   data green(38,18,15) /   2.2399659713867783D-02  /
   data green(38,18,16) /   2.2227370052096088D-02  /
   data green(38,18,17) /   2.2048271282902657D-02  /
   data green(38,18,18) /   2.1862977719535277D-02  /
   data green(38,19, 0) /   2.3538208976748107D-02  /
   data green(38,19, 1) /   2.3531684896472015D-02  /
   data green(38,19, 2) /   2.3512145248123478D-02  /
   data green(38,19, 3) /   2.3479687355906446D-02  /
   data green(38,19, 4) /   2.3434471928432454D-02  /
   data green(38,19, 5) /   2.3376720853927020D-02  /
   data green(38,19, 6) /   2.3306714197132223D-02  /
   data green(38,19, 7) /   2.3224786477025566D-02  /
   data green(38,19, 8) /   2.3131322320999178D-02  /
   data green(38,19, 9) /   2.3026751603664089D-02  /
   data green(38,19,10) /   2.2911544186624420D-02  /
   data green(38,19,11) /   2.2786204379327279D-02  /
   data green(38,19,12) /   2.2651265240587207D-02  /
   data green(38,19,13) /   2.2507282835977551D-02  /
   data green(38,19,14) /   2.2354830558502035D-02  /
   data green(38,19,15) /   2.2194493609452581D-02  /
   data green(38,19,16) /   2.2026863723832182D-02  /
   data green(38,19,17) /   2.1852534210892997D-02  /
   data green(38,19,18) /   2.1672095365899761D-02  /
   data green(38,19,19) /   2.1486130294802046D-02  /
   data green(38,20, 0) /   2.3287796758529104D-02  /
   data green(38,20, 1) /   2.3281478930448449D-02  /
   data green(38,20, 2) /   2.3262556336633721D-02  /
   data green(38,20, 3) /   2.3231121227992448D-02  /
   data green(38,20, 4) /   2.3187325966032082D-02  /
   data green(38,20, 5) /   2.3131380976704377D-02  /
   data green(38,20, 6) /   2.3063551961727612D-02  /
   data green(38,20, 7) /   2.2984156439334255D-02  /
   data green(38,20, 8) /   2.2893559701535712D-02  /
   data green(38,20, 9) /   2.2792170286567204D-02  /
   data green(38,20,10) /   2.2680435072864429D-02  /
   data green(38,20,11) /   2.2558834104640137D-02  /
   data green(38,20,12) /   2.2427875258997365D-02  /
   data green(38,20,13) /   2.2288088860845340D-02  /
   data green(38,20,14) /   2.2140022345130024D-02  /
   data green(38,20,15) /   2.1984235056618208D-02  /
   data green(38,20,16) /   2.1821293266305733D-02  /
   data green(38,20,17) /   2.1651765471095263D-02  /
   data green(38,20,18) /   2.1476218030327864D-02  /
   data green(38,20,19) /   2.1295211179618273D-02  /
   data green(38,20,20) /   2.1109295449726336D-02  /
   data green(38,21, 0) /   2.3032988947802042D-02  /
   data green(38,21, 1) /   2.3026876458561302D-02  /
   data green(38,21, 2) /   2.3008568223560198D-02  /
   data green(38,21, 3) /   2.2978151551892267D-02  /
   data green(38,21, 4) /   2.2935770671988983D-02  /
   data green(38,21, 5) /   2.2881624836705982D-02  /
   data green(38,21, 6) /   2.2815965739228786D-02  /
   data green(38,21, 7) /   2.2739094305045792D-02  /
   data green(38,21, 8) /   2.2651356939078736D-02  /
   data green(38,21, 9) /   2.2553141317725272D-02  /
   data green(38,21,10) /   2.2444871822767824D-02  /
   data green(38,21,11) /   2.2327004717748476D-02  /
   data green(38,21,12) /   2.2200023167593765D-02  /
   data green(38,21,13) /   2.2064432199256669D-02  /
   data green(38,21,14) /   2.1920753695318734D-02  /
   data green(38,21,15) /   2.1769521504352077D-02  /
   data green(38,21,16) /   2.1611276741927293D-02  /
   data green(38,21,17) /   2.1446563345035748D-02  /
   data green(38,21,18) /   2.1275923930924690D-02  /
   data green(38,21,19) /   2.1099895999431673D-02  /
   data green(38,21,20) /   2.0919008506294457D-02  /
   data green(38,21,21) /   2.0733778823973103D-02  /
   data green(38,22, 0) /   2.2774552564724229D-02  /
   data green(38,22, 1) /   2.2768643721365753D-02  /
   data green(38,22, 2) /   2.2750944816605046D-02  /
   data green(38,22, 3) /   2.2721538366990202D-02  /
   data green(38,22, 4) /   2.2680560711300455D-02  /
   data green(38,22, 5) /   2.2628200259107924D-02  /
   data green(38,22, 6) /   2.2564695100913350D-02  /
   data green(38,22, 7) /   2.2490330038880683D-02  /
   data green(38,22, 8) /   2.2405433109816079D-02  /
   data green(38,22, 9) /   2.2310371681838574D-02  /
   data green(38,22,10) /   2.2205548212910181D-02  /
   data green(38,22,11) /   2.2091395762941418D-02  /
   data green(38,22,12) /   2.1968373351634110D-02  /
   data green(38,22,13) /   2.1836961251783329D-02  /
   data green(38,22,14) /   2.1697656302771012D-02  /
   data green(38,22,15) /   2.1550967321868311D-02  /
   data green(38,22,16) /   2.1397410682201433D-02  /
   data green(38,22,17) /   2.1237506116326015D-02  /
   data green(38,22,18) /   2.1071772793788759D-02  /
   data green(38,22,19) /   2.0900725710289392D-02  /
   data green(38,22,20) /   2.0724872415494460D-02  /
   data green(38,22,21) /   2.0544710096534913D-02  /
   data green(38,22,22) /   2.0360723025006997D-02  /
   data green(38,23, 0) /   2.2513213582418867D-02  /
   data green(38,23, 1) /   2.2507505999934761D-02  /
   data green(38,23, 2) /   2.2490409325424839D-02  /
   data green(38,23, 3) /   2.2462001446320910D-02  /
   data green(38,23, 4) /   2.2422411078632170D-02  /
   data green(38,23, 5) /   2.2371816150998553D-02  /
   data green(38,23, 6) /   2.2310441598391458D-02  /
   data green(38,23, 7) /   2.2238556618736059D-02  /
   data green(38,23, 8) /   2.2156471457209094D-02  /
   data green(38,23, 9) /   2.2064533791951616D-02  /
   data green(38,23,10) /   2.1963124801190562D-02  /
   data green(38,23,11) /   2.1852654995193010D-02  /
   data green(38,23,12) /   2.1733559897134144D-02  /
   data green(38,23,13) /   2.1606295655023929D-02  /
   data green(38,23,14) /   2.1471334662592510D-02  /
   data green(38,23,15) /   2.1329161260847652D-02  /
   data green(38,23,16) /   2.1180267584304545D-02  /
   data green(38,23,17) /   2.1025149607088147D-02  /
   data green(38,23,18) /   2.0864303434656260D-02  /
   data green(38,23,19) /   2.0698221877195049D-02  /
   data green(38,23,20) /   2.0527391331163727D-02  /
   data green(38,23,21) /   2.0352288986319165D-02  /
   data green(38,23,22) /   2.0173380367080867D-02  /
   data green(38,23,23) /   1.9991117209484793D-02  /
   data green(38,24, 0) /   2.2249655123721762D-02  /
   data green(38,24, 1) /   2.2244145810498158D-02  /
   data green(38,24, 2) /   2.2227642450040920D-02  /
   data green(38,24, 3) /   2.2200218474908062D-02  /
   data green(38,24, 4) /   2.2161995263191969D-02  /
   data green(38,24, 5) /   2.2113140649995821D-02  /
   data green(38,24, 6) /   2.2053866893902409D-02  /
   data green(38,24, 7) /   2.1984428146419769D-02  /
   data green(38,24, 8) /   2.1905117482806345D-02  /
   data green(38,24, 9) /   2.1816263560898046D-02  /
   data green(38,24,10) /   2.1718226980361282D-02  /
   data green(38,24,11) /   2.1611396418092267D-02  /
   data green(38,24,12) /   2.1496184616306362D-02  /
   data green(38,24,13) /   2.1373024298360067D-02  /
   data green(38,24,14) /   2.1242364083763957D-02  /
   data green(38,24,15) /   2.1104664468490530D-02  /
   data green(38,24,16) /   2.0960393929919541D-02  /
   data green(38,24,17) /   2.0810025207976045D-02  /
   data green(38,24,18) /   2.0654031805590056D-02  /
   data green(38,24,19) /   2.0492884742902100D-02  /
   data green(38,24,20) /   2.0327049590985958D-02  /
   data green(38,24,21) /   2.0156983802538546D-02  /
   data green(38,24,22) /   1.9983134349225166D-02  /
   data green(38,24,23) /   1.9805935668338222D-02  /
   data green(38,24,24) /   1.9625807915244889D-02  /
   data green(38,25, 0) /   2.1984516391577721D-02  /
   data green(38,25, 1) /   2.1979201831588868D-02  /
   data green(38,25, 2) /   2.1963281298338053D-02  /
   data green(38,25, 3) /   2.1936823951523625D-02  /
   data green(38,25, 4) /   2.1899944129153078D-02  /
   data green(38,25, 5) /   2.1852799978403191D-02  /
   data green(38,25, 6) /   2.1795591584081876D-02  /
   data green(38,25, 7) /   2.1728558637830006D-02  /
   data green(38,25, 8) /   2.1651977700638814D-02  /
   data green(38,25, 9) /   2.1566159118759833D-02  /
   data green(38,25,10) /   2.1471443658451942D-02  /
   data green(38,25,11) /   2.1368198928159540D-02  /
   data green(38,25,12) /   2.1256815657665448D-02  /
   data green(38,25,13) /   2.1137703902633889D-02  /
   data green(38,25,14) /   2.1011289239955218D-02  /
   data green(38,25,15) /   2.0878009014693093D-02  /
   data green(38,25,16) /   2.0738308693529623D-02  /
   data green(38,25,17) /   2.0592638372737771D-02  /
   data green(38,25,18) /   2.0441449481221816D-02  /
   data green(38,25,19) /   2.0285191711377187D-02  /
   data green(38,25,20) /   2.0124310202725704D-02  /
   data green(38,25,21) /   1.9959242995734063D-02  /
   data green(38,25,22) /   1.9790418766134153D-02  /
   data green(38,25,23) /   1.9618254843591851D-02  /
   data green(38,25,24) /   1.9443155512832242D-02  /
   data green(38,25,25) /   1.9265510590389362D-02  /
   data green(38,26, 0) /   2.1718392248634571D-02  /
   data green(38,26, 1) /   2.1713268479500072D-02  /
   data green(38,26, 2) /   2.1697918949219182D-02  /
   data green(38,26, 3) /   2.1672408731693488D-02  /
   data green(38,26, 4) /   2.1636845431182869D-02  /
   data green(38,26, 5) /   2.1591377924633632D-02  /
   data green(38,26, 6) /   2.1536194641507991D-02  /
   data green(38,26, 7) /   2.1471521419831371D-02  /
   data green(38,26, 8) /   2.1397618985706232D-02  /
   data green(38,26, 9) /   2.1314780110373303D-02  /
   data green(38,26,10) /   2.1223326503855094D-02  /
   data green(38,26,11) /   2.1123605507208605D-02  /
   data green(38,26,12) /   2.1015986646456240D-02  /
   data green(38,26,13) /   2.0900858110451045D-02  /
   data green(38,26,14) /   2.0778623212436641D-02  /
   data green(38,26,15) /   2.0649696891110379D-02  /
   data green(38,26,16) /   2.0514502301860516D-02  /
   data green(38,26,17) /   2.0373467542813998D-02  /
   data green(38,26,18) /   2.0227022553697217D-02  /
   data green(38,26,19) /   2.0075596218562357D-02  /
   data green(38,26,20) /   1.9919613696429691D-02  /
   data green(38,26,21) /   1.9759493997070478D-02  /
   data green(38,26,22) /   1.9595647812698475D-02  /
   data green(38,26,23) /   1.9428475610400366D-02  /
   data green(38,26,24) /   1.9258365984825131D-02  /
   data green(38,26,25) /   1.9085694266038910D-02  /
   data green(38,26,26) /   1.8910821373568660D-02  /
   data green(38,27, 0) /   2.1451833360527455D-02  /
   data green(38,27, 1) /   2.1446896046781395D-02  /
   data green(38,27, 2) /   2.1432104576845301D-02  /
   data green(38,27, 3) /   2.1407520128535196D-02  /
   data green(38,27, 4) /   2.1373243883256884D-02  /
   data green(38,27, 5) /   2.1329415872086720D-02  /
   data green(38,27, 6) /   2.1276213396598461D-02  /
   data green(38,27, 7) /   2.1213849059130940D-02  /
   data green(38,27, 8) /   2.1142568444893452D-02  /
   data green(38,27, 9) /   2.1062647504516776D-02  /
   data green(38,27,10) /   2.0974389690217497D-02  /
   data green(38,27,11) /   2.0878122901574411D-02  /
   data green(38,27,12) /   2.0774196298018278D-02  /
   data green(38,27,13) /   2.0662977034589209D-02  /
   data green(38,27,14) /   2.0544846975460211D-02  /
   data green(38,27,15) /   2.0420199436356260D-02  /
   data green(38,27,16) /   2.0289436002543473D-02  /
   data green(38,27,17) /   2.0152963463777960D-02  /
   data green(38,27,18) /   2.0011190901743411D-02  /
   data green(38,27,19) /   1.9864526959323493D-02  /
   data green(38,27,20) /   1.9713377314783273D-02  /
   data green(38,27,21) /   1.9558142377778573D-02  /
   data green(38,27,22) /   1.9399215218249190D-02  /
   data green(38,27,23) /   1.9236979733820374D-02  /
   data green(38,27,24) /   1.9071809056439581D-02  /
   data green(38,27,25) /   1.8904064194682092D-02  /
   data green(38,27,26) /   1.8734092904505048D-02  /
   data green(38,27,27) /   1.8562228778224117D-02  /
   data green(38,28, 0) /   2.1185346818943265D-02  /
   data green(38,28, 1) /   2.1180591320073505D-02  /
   data green(38,28, 2) /   2.1166344052887413D-02  /
   data green(38,28, 3) /   2.1142662489389723D-02  /
   data green(38,28, 4) /   2.1109641700147723D-02  /
   data green(38,28, 5) /   2.1067413296670300D-02  /
   data green(38,28, 6) /   2.1016143983199120D-02  /
   data green(38,28, 7) /   2.0956033748961531D-02  /
   data green(38,28, 8) /   2.0887313738875986D-02  /
   data green(38,28, 9) /   2.0810243846338148D-02  /
   data green(38,28,10) /   2.0725110075903611D-02  /
   data green(38,28,11) /   2.0632221726349143D-02  /
   data green(38,28,12) /   2.0531908445732942D-02  /
   data green(38,28,13) /   2.0424517209747216D-02  /
   data green(38,28,14) /   2.0310409272981413D-02  /
   data green(38,28,15) /   2.0189957139854560D-02  /
   data green(38,28,16) /   2.0063541598129307D-02  /
   data green(38,28,17) /   1.9931548853303088D-02  /
   data green(38,28,18) /   1.9794367797010968D-02  /
   data green(38,28,19) /   1.9652387437088550D-02  /
   data green(38,28,20) /   1.9505994511340354D-02  /
   data green(38,28,21) /   1.9355571301523437D-02  /
   data green(38,28,22) /   1.9201493658747029D-02  /
   data green(38,28,23) /   1.9044129246534643D-02  /
   data green(38,28,24) /   1.8883836003292815D-02  /
   data green(38,28,25) /   1.8720960821948445D-02  /
   data green(38,28,26) /   1.8555838441092466D-02  /
   data green(38,28,27) /   1.8388790539116077D-02  /
   data green(38,28,28) /   1.8220125020539241D-02  /
   data green(38,29, 0) /   2.0919397164209198D-02  /
   data green(38,29, 1) /   2.0914818597208347D-02  /
   data green(38,29, 2) /   2.0901100947263350D-02  /
   data green(38,29, 3) /   2.0878298169622455D-02  /
   data green(38,29, 4) /   2.0846499534211971D-02  /
   data green(38,29, 5) /   2.0805828657165713D-02  /
   data green(38,29, 6) /   2.0756442173948810D-02  /
   data green(38,29, 7) /   2.0698528081829825D-02  /
   data green(38,29, 8) /   2.0632303785703274D-02  /
   data green(38,29, 9) /   2.0558013886372347D-02  /
   data green(38,29,10) /   2.0475927754240279D-02  /
   data green(38,29,11) /   2.0386336933859463D-02  /
   data green(38,29,12) /   2.0289552425939657D-02  /
   data green(38,29,13) /   2.0185901893271057D-02  /
   data green(38,29,14) /   2.0075726835668037D-02  /
   data green(38,29,15) /   1.9959379776625801D-02  /
   data green(38,29,16) /   1.9837221501072224D-02  /
   data green(38,29,17) /   1.9709618379577278D-02  /
   data green(38,29,18) /   1.9576939809848766D-02  /
   data green(38,29,19) /   1.9439555801488182D-02  /
   data green(38,29,20) /   1.9297834724987622D-02  /
   data green(38,29,21) /   1.9152141240983467D-02  /
   data green(38,29,22) /   1.9002834420987241D-02  /
   data green(38,29,23) /   1.8850266066307925D-02  /
   data green(38,29,24) /   1.8694779227753242D-02  /
   data green(38,29,25) /   1.8536706925015838D-02  /
   data green(38,29,26) /   1.8376371061453539D-02  /
   data green(38,29,27) /   1.8214081527279387D-02  /
   data green(38,29,28) /   1.8050135481985238D-02  /
   data green(38,29,29) /   1.7884816805115174D-02  /
   data green(38,30, 0) /   2.0654407732309577D-02  /
   data green(38,30, 1) /   2.0650001028633923D-02  /
   data green(38,30, 2) /   2.0636797852882156D-02  /
   data green(38,30, 3) /   2.0614848829883129D-02  /
   data green(38,30, 4) /   2.0584237734826392D-02  /
   data green(38,30, 5) /   2.0545080607134135D-02  /
   data green(38,30, 6) /   2.0497524535732614D-02  /
   data green(38,30, 7) /   2.0441746140506183D-02  /
   data green(38,30, 8) /   2.0377949780333211D-02  /
   data green(38,30, 9) /   2.0306365522726036D-02  /
   data green(38,30,10) /   2.0227246913607368D-02  /
   data green(38,30,11) /   2.0140868588093306D-02  /
   data green(38,30,12) /   2.0047523764302508D-02  /
   data green(38,30,13) /   1.9947521662210126D-02  /
   data green(38,30,14) /   1.9841184888494186D-02  /
   data green(38,30,15) /   1.9728846826295221D-02  /
   data green(38,30,16) /   1.9610849065971780D-02  /
   data green(38,30,17) /   1.9487538909444643D-02  /
   data green(38,30,18) /   1.9359266976750063D-02  /
   data green(38,30,19) /   1.9226384939136785D-02  /
   data green(38,30,20) /   1.9089243398603187D-02  /
   data green(38,30,21) /   1.8948189929327605D-02  /
   data green(38,30,22) /   1.8803567292124967D-02  /
   data green(38,30,23) /   1.8655711828974963D-02  /
   data green(38,30,24) /   1.8504952040895004D-02  /
   data green(38,30,25) /   1.8351607349038169D-02  /
   data green(38,30,26) /   1.8195987035922883D-02  /
   data green(38,30,27) /   1.8038389361167282D-02  /
   data green(38,30,28) /   1.7879100844011332D-02  /
   data green(38,30,29) /   1.7718395703252194D-02  /
   data green(38,30,30) /   1.7556535443970664D-02  /
   data green(38,31, 0) /   2.0390762257413798D-02  /
   data green(38,31, 1) /   2.0386522214374710D-02  /
   data green(38,31, 2) /   2.0373817965990269D-02  /
   data green(38,31, 3) /   2.0352696989053540D-02  /
   data green(38,31, 4) /   2.0323237863583369D-02  /
   data green(38,31, 5) /   2.0285549462592180D-02  /
   data green(38,31, 6) /   2.0239769840798321D-02  /
   data green(38,31, 7) /   2.0186064844393219D-02  /
   data green(38,31, 8) /   2.0124626469017342D-02  /
   data green(38,31, 9) /   2.0055670997275950D-02  /
   data green(38,31,10) /   1.9979436950332097D-02  /
   data green(38,31,11) /   1.9896182890291009D-02  /
   data green(38,31,12) /   1.9806185111222244D-02  /
   data green(38,31,13) /   1.9709735256781186D-02  /
   data green(38,31,14) /   1.9607137901555865D-02  /
   data green(38,31,15) /   1.9498708131572562D-02  /
   data green(38,31,16) /   1.9384769156969364D-02  /
   data green(38,31,17) /   1.9265649986825450D-02  /
   data green(38,31,18) /   1.9141683192661813D-02  /
   data green(38,31,19) /   1.9013202783354801D-02  /
   data green(38,31,20) /   1.8880542210267602D-02  /
   data green(38,31,21) /   1.8744032517437399D-02  /
   data green(38,31,22) /   1.8604000647773906D-02  /
   data green(38,31,23) /   1.8460767912525086D-02  /
   data green(38,31,24) /   1.8314648627828202D-02  /
   data green(38,31,25) /   1.8165948919045861D-02  /
   data green(38,31,26) /   1.8014965690830487D-02  /
   data green(38,31,27) /   1.7861985758486642D-02  /
   data green(38,31,28) /   1.7707285134217612D-02  /
   data green(38,31,29) /   1.7551128460244772D-02  /
   data green(38,31,30) /   1.7393768579556993D-02  /
   data green(38,31,31) /   1.7235446234160682D-02  /
   data green(38,32, 0) /   2.0128806667777163D-02  /
   data green(38,32, 1) /   2.0124727994492261D-02  /
   data green(38,32, 2) /   2.0112506860391646D-02  /
   data green(38,32, 3) /   2.0092187771666166D-02  /
   data green(38,32, 4) /   2.0063844404738591D-02  /
   data green(38,32, 5) /   2.0027578865854663D-02  /
   data green(38,32, 6) /   1.9983520675026085D-02  /
   data green(38,32, 7) /   1.9931825494035854D-02  /
   data green(38,32, 8) /   1.9872673622742209D-02  /
   data green(38,32, 9) /   1.9806268291687704D-02  /
   data green(38,32,10) /   1.9732833781942998D-02  /
   data green(38,32,11) /   1.9652613405136504D-02  /
   data green(38,32,12) /   1.9565867377724866D-02  /
   data green(38,32,13) /   1.9472870623764727D-02  /
   data green(38,32,14) /   1.9373910539808081D-02  /
   data green(38,32,15) /   1.9269284754140255D-02  /
   data green(38,32,16) /   1.9159298910515151D-02  /
   data green(38,32,17) /   1.9044264503933952D-02  /
   data green(38,32,18) /   1.8924496792986828D-02  /
   data green(38,32,19) /   1.8800312809959533D-02  /
   data green(38,32,20) /   1.8672029486423981D-02  /
   data green(38,32,21) /   1.8539961908496855D-02  /
   data green(38,32,22) /   1.8404421712469626D-02  /
   data green(38,32,23) /   1.8265715628173005D-02  /
   data green(38,32,24) /   1.8124144174312691D-02  /
   data green(38,32,25) /   1.7980000507156656D-02  /
   data green(38,32,26) /   1.7833569421406278D-02  /
   data green(38,32,27) /   1.7685126499868021D-02  /
   data green(38,32,28) /   1.7534937406669830D-02  /
   data green(38,32,29) /   1.7383257317234106D-02  /
   data green(38,32,30) /   1.7230330477016314D-02  /
   data green(38,32,31) /   1.7076389880124523D-02  /
   data green(38,32,32) /   1.6921657058326163D-02  /
   data green(38,33, 0) /   1.9868851019913706D-02  /
   data green(38,33, 1) /   1.9864928378024671D-02  /
   data green(38,33, 2) /   1.9853174400759295D-02  /
   data green(38,33, 3) /   1.9833630795407450D-02  /
   data green(38,33, 4) /   1.9806366617077895D-02  /
   data green(38,33, 5) /   1.9771477592418808D-02  /
   data green(38,33, 6) /   1.9729085191084818D-02  /
   data green(38,33, 7) /   1.9679335462511098D-02  /
   data green(38,33, 8) /   1.9622397659611962D-02  /
   data green(38,33, 9) /   1.9558462674421691D-02  /
   data green(38,33,10) /   1.9487741313355934D-02  /
   data green(38,33,11) /   1.9410462441644450D-02  /
   data green(38,33,12) /   1.9326871027552418D-02  /
   data green(38,33,13) /   1.9237226117282319D-02  /
   data green(38,33,14) /   1.9141798770974657D-02  /
   data green(38,33,15) /   1.9040869989070348D-02  /
   data green(38,33,16) /   1.8934728656546530D-02  /
   data green(38,33,17) /   1.8823669530291735D-02  /
   data green(38,33,18) /   1.8707991292253503D-02  /
   data green(38,33,19) /   1.8587994688082237D-02  /
   data green(38,33,20) /   1.8463980767919004D-02  /
   data green(38,33,21) /   1.8336249242832513D-02  /
   data green(38,33,22) /   1.8205096967295647D-02  /
   data green(38,33,23) /   1.8070816555082754D-02  /
   data green(38,33,24) /   1.7933695133133771D-02  /
   data green(38,33,25) /   1.7794013235320337D-02  /
   data green(38,33,26) /   1.7652043835701379D-02  /
   data green(38,33,27) /   1.7508051518795063D-02  /
   data green(38,33,28) /   1.7362291782632794D-02  /
   data green(38,33,29) /   1.7215010468900512D-02  /
   data green(38,33,30) /   1.7066443313305112D-02  /
   data green(38,33,31) /   1.6916815608415157D-02  /
   data green(38,33,32) /   1.6766341970593682D-02  /
   data green(38,33,33) /   1.6615226202243955D-02  /
   data green(38,34, 0) /   1.9611171522958964D-02  /
   data green(38,34, 1) /   1.9607399562382093D-02  /
   data green(38,34, 2) /   1.9596096747192593D-02  /
   data green(38,34, 3) /   1.9577302151155912D-02  /
   data green(38,34, 4) /   1.9551080480067526D-02  /
   data green(38,34, 5) /   1.9517521454285702D-02  /
   data green(38,34, 6) /   1.9476738960519224D-02  /
   data green(38,34, 7) /   1.9428869988506719D-02  /
   data green(38,34, 8) /   1.9374073371859379D-02  /
   data green(38,34, 9) /   1.9312528355401717D-02  /
   data green(38,34,10) /   1.9244433013764446D-02  /
   data green(38,34,11) /   1.9170002547712638D-02  /
   data green(38,34,12) /   1.9089467485716022D-02  /
   data green(38,34,13) /   1.9003071818593155D-02  /
   data green(38,34,14) /   1.8911071094724740D-02  /
   data green(38,34,15) /   1.8813730502386498D-02  /
   data green(38,34,16) /   1.8711322964272813D-02  /
   data green(38,34,17) /   1.8604127267354158D-02  /
   data green(38,34,18) /   1.8492426248925563D-02  /
   data green(38,34,19) /   1.8376505057158019D-02  /
   data green(38,34,20) /   1.8256649501752473D-02  /
   data green(38,34,21) /   1.8133144507508404D-02  /
   data green(38,34,22) /   1.8006272680836665D-02  /
   data green(38,34,23) /   1.7876312996541740D-02  /
   data green(38,34,24) /   1.7743539609632302D-02  /
   data green(38,34,25) /   1.7608220794539518D-02  /
   data green(38,34,26) /   1.7470618011964680D-02  /
   data green(38,34,27) /   1.7330985101668143D-02  /
   data green(38,34,28) /   1.7189567597861347D-02  /
   data green(38,34,29) /   1.7046602162478518D-02  /
   data green(38,34,30) /   1.6902316130479829D-02  /
   data green(38,34,31) /   1.6756927160462214D-02  /
   data green(38,34,32) /   1.6610642983211205D-02  /
   data green(38,34,33) /   1.6463661240399174D-02  /
   data green(38,34,34) /   1.6316169405396642D-02  /
   data green(38,35, 0) /   1.9356012611935071D-02  /
   data green(38,35, 1) /   1.9352386001951764D-02  /
   data green(38,35, 2) /   1.9341518409898031D-02  /
   data green(38,35, 3) /   1.9323446434640498D-02  /
   data green(38,35, 4) /   1.9298230693664976D-02  /
   data green(38,35, 5) /   1.9265955259474311D-02  /
   data green(38,35, 6) /   1.9226726884988041D-02  /
   data green(38,35, 7) /   1.9180674031861907D-02  /
   data green(38,35, 8) /   1.9127945718900311D-02  /
   data green(38,35, 9) /   1.9068710210492300D-02  /
   data green(38,35,10) /   1.9003153567198065D-02  /
   data green(38,35,11) /   1.8931478082207357D-02  /
   data green(38,35,12) /   1.8853900628366182D-02  /
   data green(38,35,13) /   1.8770650940829622D-02  /
   data green(38,35,14) /   1.8681969860173470D-02  /
   data green(38,35,15) /   1.8588107560032773D-02  /
   data green(38,35,16) /   1.8489321782090042D-02  /
   data green(38,35,17) /   1.8385876099584964D-02  /
   data green(38,35,18) /   1.8278038228537907D-02  /
   data green(38,35,19) /   1.8166078403655955D-02  /
   data green(38,35,20) /   1.8050267833503938D-02  /
   data green(38,35,21) /   1.7930877247053952D-02  /
   data green(38,35,22) /   1.7808175541246737D-02  /
   data green(38,35,23) /   1.7682428536771769D-02  /
   data green(38,35,24) /   1.7553897846955734D-02  /
   data green(38,35,25) /   1.7422839862484829D-02  /
   data green(38,35,26) /   1.7289504852708935D-02  /
   data green(38,35,27) /   1.7154136182510740D-02  /
   data green(38,35,28) /   1.7016969642182963D-02  /
   data green(38,35,29) /   1.6878232886448748D-02  /
   data green(38,35,30) /   1.6738144977683007D-02  /
   data green(38,35,31) /   1.6596916027537125D-02  /
   data green(38,35,32) /   1.6454746930524851D-02  /
   data green(38,35,33) /   1.6311829182676492D-02  /
   data green(38,35,34) /   1.6168344778093004D-02  /
   data green(38,35,35) /   1.6024466176110960D-02  /
   data green(38,36, 0) /   1.9103589035047264D-02  /
   data green(38,36, 1) /   1.9100102491068058D-02  /
   data green(38,36, 2) /   1.9089654319227725D-02  /
   data green(38,36, 3) /   1.9072278795085797D-02  /
   data green(38,36, 4) /   1.9048032697345736D-02  /
   data green(38,36, 5) /   1.9016994793529765D-02  /
   data green(38,36, 6) /   1.8979265132762886D-02  /
   data green(38,36, 7) /   1.8934964158050890D-02  /
   data green(38,36, 8) /   1.8884231653350438D-02  /
   data green(38,36, 9) /   1.8827225543209699D-02  /
   data green(38,36,10) /   1.8764120564750748D-02  /
   data green(38,36,11) /   1.8695106833231463D-02  /
   data green(38,36,12) /   1.8620388323348502D-02  /
   data green(38,36,13) /   1.8540181288828138D-02  /
   data green(38,36,14) /   1.8454712642717631D-02  /
   data green(38,36,15) /   1.8364218320177240D-02  /
   data green(38,36,16) /   1.8268941644529076D-02  /
   data green(38,36,17) /   1.8169131715909635D-02  /
   data green(38,36,18) /   1.8065041840161583D-02  /
   data green(38,36,19) /   1.7956928013661774D-02  /
   data green(38,36,20) /   1.7845047477686931D-02  /
   data green(38,36,21) /   1.7729657353735002D-02  /
   data green(38,36,22) /   1.7611013369013308D-02  /
   data green(38,36,23) /   1.7489368679132228D-02  /
   data green(38,36,24) /   1.7364972792954386D-02  /
   data green(38,36,25) /   1.7238070602585236D-02  /
   data green(38,36,26) /   1.7108901519683851D-02  /
   data green(38,36,27) /   1.6977698717645758D-02  /
   data green(38,36,28) /   1.6844688477777418D-02  /
   data green(38,36,29) /   1.6710089636352284D-02  /
   data green(38,36,30) /   1.6574113128411533D-02  /
   data green(38,36,31) /   1.6436961623344289D-02  /
   data green(38,36,32) /   1.6298829246641951D-02  /
   data green(38,36,33) /   1.6159901381757370D-02  /
   data green(38,36,34) /   1.6020354545695822D-02  /
   data green(38,36,35) /   1.5880356331804069D-02  /
   data green(38,36,36) /   1.5740065413188665D-02  /
   data green(38,37, 0) /   1.8854087926079652D-02  /
   data green(38,37, 1) /   1.8850736232428677D-02  /
   data green(38,37, 2) /   1.8840691882198026D-02  /
   data green(38,37, 3) /   1.8823986972037448D-02  /
   data green(38,37, 4) /   1.8800674679641513D-02  /
   data green(38,37, 5) /   1.8770828794458738D-02  /
   data green(38,37, 6) /   1.8734543072095919D-02  /
   data green(38,37, 7) /   1.8691930423436372D-02  /
   data green(38,37, 8) /   1.8643121952096452D-02  /
   data green(38,37, 9) /   1.8588265856075314D-02  /
   data green(38,37,10) /   1.8527526211258475D-02  /
   data green(38,37,11) /   1.8461081655782387D-02  /
   data green(38,37,12) /   1.8389123995138913D-02  /
   data green(38,37,13) /   1.8311856748296454D-02  /
   data green(38,37,14) /   1.8229493655054701D-02  /
   data green(38,37,15) /   1.8142257164364167D-02  /
   data green(38,37,16) /   1.8050376922472341D-02  /
   data green(38,37,17) /   1.7954088278556348D-02  /
   data green(38,37,18) /   1.7853630824027574D-02  /
   data green(38,37,19) /   1.7749246980005712D-02  /
   data green(38,37,20) /   1.7641180645622295D-02  /
   data green(38,37,21) /   1.7529675917886345D-02  /
   data green(38,37,22) /   1.7414975891884572D-02  /
   data green(38,37,23) /   1.7297321548146311D-02  /
   data green(38,37,24) /   1.7176950732124451D-02  /
   data green(38,37,25) /   1.7054097228965177D-02  /
   data green(38,37,26) /   1.6928989935091111D-02  /
   data green(38,37,27) /   1.6801852126626989D-02  /
   data green(38,37,28) /   1.6672900823369374D-02  /
   data green(38,37,29) /   1.6542346245849784D-02  /
   data green(38,37,30) /   1.6410391362067477D-02  /
   data green(38,37,31) /   1.6277231519670633D-02  /
   data green(38,37,32) /   1.6143054158735459D-02  /
   data green(38,37,33) /   1.6008038599823051D-02  /
   data green(38,37,34) /   1.5872355901669761D-02  /
   data green(38,37,35) /   1.5736168782673952D-02  /
   data green(38,37,36) /   1.5599631600266700D-02  /
   data green(38,37,37) /   1.5462890382278355D-02  /
   data green(38,38, 0) /   1.8607670838356079D-02  /
   data green(38,38, 1) /   1.8604448867426539D-02  /
   data green(38,38, 2) /   1.8594793001968792D-02  /
   data green(38,38, 3) /   1.8578733296867902D-02  /
   data green(38,38, 4) /   1.8556319554721912D-02  /
   data green(38,38, 5) /   1.8527620897671619D-02  /
   data green(38,38, 6) /   1.8492725178104363D-02  /
   data green(38,38, 7) /   1.8451738238033370D-02  /
   data green(38,38, 8) /   1.8404783029285457D-02  /
   data green(38,38, 9) /   1.8351998608633738D-02  /
   data green(38,38,10) /   1.8293539023647356D-02  /
   data green(38,38,11) /   1.8229572106264594D-02  /
   data green(38,38,12) /   1.8160278191914670D-02  /
   data green(38,38,13) /   1.8085848782416095D-02  /
   data green(38,38,14) /   1.8006485170878519D-02  /
   data green(38,38,15) /   1.7922397046456470D-02  /
   data green(38,38,16) /   1.7833801096081316D-02  /
   data green(38,38,17) /   1.7740919619277904D-02  /
   data green(38,38,18) /   1.7643979170902508D-02  /
   data green(38,38,19) /   1.7543209245172789D-02  /
   data green(38,38,20) /   1.7438841012751137D-02  /
   data green(38,38,21) /   1.7331106120944701D-02  /
   data green(38,38,22) /   1.7220235565345958D-02  /
   data green(38,38,23) /   1.7106458639505162D-02  /
   data green(38,38,24) /   1.6990001967537929D-02  /
   data green(38,38,25) /   1.6871088622964422D-02  /
   data green(38,38,26) /   1.6749937335576500D-02  /
   data green(38,38,27) /   1.6626761786757951D-02  /
   data green(38,38,28) /   1.6501769992455652D-02  /
   data green(38,38,29) /   1.6375163771924474D-02  /
   data green(38,38,30) /   1.6247138299450029D-02  /
   data green(38,38,31) /   1.6117881735490096D-02  /
   data green(38,38,32) /   1.5987574933062843D-02  /
   data green(38,38,33) /   1.5856391214739966D-02  /
   data green(38,38,34) /   1.5724496215265549D-02  /
   data green(38,38,35) /   1.5592047784604315D-02  /
   data green(38,38,36) /   1.5459195946114277D-02  /
   data green(38,38,37) /   1.5326082904523405D-02  /
   data green(38,38,38) /   1.5192843098456147D-02  /
   data green(39, 0, 0) /   2.5645247716982612D-02  /
   data green(39, 1, 0) /   2.5636804803020830D-02  /
   data green(39, 1, 1) /   2.5628370249880468D-02  /
   data green(39, 2, 0) /   2.5611526226162473D-02  /
   data green(39, 2, 1) /   2.5603116672359906D-02  /
   data green(39, 2, 2) /   2.5577937844636509D-02  /
   data green(39, 3, 0) /   2.5569561650217652D-02  /
   data green(39, 3, 1) /   2.5561193486828532D-02  /
   data green(39, 3, 2) /   2.5536138420494282D-02  /
   data green(39, 3, 3) /   2.5494543908682967D-02  /
   data green(39, 4, 0) /   2.5511157771110238D-02  /
   data green(39, 4, 1) /   2.5502846983189852D-02  /
   data green(39, 4, 2) /   2.5477963477418143D-02  /
   data green(39, 4, 3) /   2.5436653026821333D-02  /
   data green(39, 4, 4) /   2.5379155942504648D-02  /
   data green(39, 5, 0) /   2.5436654294766887D-02  /
   data green(39, 5, 1) /   2.5428416312060134D-02  /
   data green(39, 5, 2) /   2.5403750507757221D-02  /
   data green(39, 5, 3) /   2.5362800528798253D-02  /
   data green(39, 5, 4) /   2.5305803196894298D-02  /
   data green(39, 5, 5) /   2.5233084710541580D-02  /
   data green(39, 6, 0) /   2.5346478488027662D-02  /
   data green(39, 6, 1) /   2.5338328048141387D-02  /
   data green(39, 6, 2) /   2.5313924019263394D-02  /
   data green(39, 6, 3) /   2.5273407508716421D-02  /
   data green(39, 6, 4) /   2.5217011169683771D-02  /
   data green(39, 6, 5) /   2.5145055496704288D-02  /
   data green(39, 6, 6) /   2.5057943802871553D-02  /
   data green(39, 7, 0) /   2.5241138472764831D-02  /
   data green(39, 7, 1) /   2.5233089499474514D-02  /
   data green(39, 7, 2) /   2.5208988889482215D-02  /
   data green(39, 7, 3) /   2.5168974829638927D-02  /
   data green(39, 7, 4) /   2.5113275178876016D-02  /
   data green(39, 7, 5) /   2.5042203870366920D-02  /
   data green(39, 7, 6) /   2.4956156031996270D-02  /
   data green(39, 7, 7) /   2.4855601972799259D-02  /
   data green(39, 8, 0) /   2.5121215466387285D-02  /
   data green(39, 8, 1) /   2.5113280965162875D-02  /
   data green(39, 8, 2) /   2.5089522674838789D-02  /
   data green(39, 8, 3) /   2.5050075517301904D-02  /
   data green(39, 8, 4) /   2.4995161990957641D-02  /
   data green(39, 8, 5) /   2.4925088691072583D-02  /
   data green(39, 8, 6) /   2.4840241589042591D-02  /
   data green(39, 8, 7) /   2.4741080212100322D-02  /
   data green(39, 8, 8) /   2.4628130892956689D-02  /
   data green(39, 9, 0) /   2.4987355194887196D-02  /
   data green(39, 9, 1) /   2.4979547167086353D-02  /
   data green(39, 9, 2) /   2.4956167098038193D-02  /
   data green(39, 9, 3) /   2.4917346339289090D-02  /
   data green(39, 9, 4) /   2.4863301526417230D-02  /
   data green(39, 9, 5) /   2.4794331227302700D-02  /
   data green(39, 9, 6) /   2.4710811393366512D-02  /
   data green(39, 9, 7) /   2.4613189748701025D-02  /
   data green(39, 9, 8) /   2.4501979278823544D-02  /
   data green(39, 9, 9) /   2.4377751000008812D-02  /
   data green(39,10, 0) /   2.4840258716999311D-02  /
   data green(39,10, 1) /   2.4832588093521599D-02  /
   data green(39,10, 2) /   2.4809618949750787D-02  /
   data green(39,10, 3) /   2.4771478802438684D-02  /
   data green(39,10, 4) /   2.4718377988988054D-02  /
   data green(39,10, 5) /   2.4650606451581097D-02  /
   data green(39,10, 6) /   2.4568529371139909D-02  /
   data green(39,10, 7) /   2.4472581779116369D-02  /
   data green(39,10, 8) /   2.4363262300677892D-02  /
   data green(39,10, 9) /   2.4241126201302553D-02  /
   data green(39,10,10) /   2.4106777919644762D-02  /
   data green(39,11, 0) /   2.4680672900038034D-02  /
   data green(39,11, 1) /   2.4673149494600789D-02  /
   data green(39,11, 2) /   2.4650620642623217D-02  /
   data green(39,11, 3) /   2.4613209803412896D-02  /
   data green(39,11, 4) /   2.4561120649272874D-02  /
   data green(39,11, 5) /   2.4494633991606122D-02  /
   data green(39,11, 6) /   2.4414103605925032D-02  /
   data green(39,11, 7) /   2.4319951076612164D-02  /
   data green(39,11, 8) /   2.4212659806556104D-02  /
   data green(39,11, 9) /   2.4092768354416382D-02  /
   data green(39,11,10) /   2.3960863272791053D-02  /
   data green(39,11,11) /   2.3817571623963956D-02  /
   data green(39,12, 0) /   2.4509380780348516D-02  /
   data green(39,12, 1) /   2.4502013261957688D-02  /
   data green(39,12, 2) /   2.4479950648216942D-02  /
   data green(39,12, 3) /   2.4443312160153677D-02  /
   data green(39,12, 4) /   2.4392294506518936D-02  /
   data green(39,12, 5) /   2.4327168956215433D-02  /
   data green(39,12, 6) /   2.4248277360317306D-02  /
   data green(39,12, 7) /   2.4156027237250570D-02  /
   data green(39,12, 8) /   2.4050886057658389D-02  /
   data green(39,12, 9) /   2.3933374882256178D-02  /
   data green(39,12,10) /   2.3804061516147012D-02  /
   data green(39,12,11) /   2.3663553346597144D-02  /
   data green(39,12,12) /   2.3512490028510094D-02  /
   data green(39,13, 0) /   2.4327192025344439D-02  /
   data green(39,13, 1) /   2.4319987909020603D-02  /
   data green(39,13, 2) /   2.4298414031798363D-02  /
   data green(39,13, 3) /   2.4262585236620289D-02  /
   data green(39,13, 4) /   2.4212691037457507D-02  /
   data green(39,13, 5) /   2.4148992840721639D-02  /
   data green(39,13, 6) /   2.4071820167994754D-02  /
   data green(39,13, 7) /   2.3981565986339298D-02  /
   data green(39,13, 8) /   2.3878681274064131D-02  /
   data green(39,13, 9) /   2.3763668965734213D-02  /
   data green(39,13,10) /   2.3637077429999993D-02  /
   data green(39,13,11) /   2.3499493637450852D-02  /
   data green(39,13,12) /   2.3351536173465949D-02  /
   data green(39,13,13) /   2.3193848243553805D-02  /
   data green(39,14, 0) /   2.4134933691465441D-02  /
   data green(39,14, 1) /   2.4127899345857250D-02  /
   data green(39,14, 2) /   2.4106833277613961D-02  /
   data green(39,14, 3) /   2.4071845851342505D-02  /
   data green(39,14, 4) /   2.4023119219847504D-02  /
   data green(39,14, 5) /   2.3960904695611563D-02  /
   data green(39,14, 6) /   2.3885519175821487D-02  /
   data green(39,14, 7) /   2.3797340719950138D-02  /
   data green(39,14, 8) /   2.3696803399174592D-02  /
   data green(39,14, 9) /   2.3584391551940871D-02  /
   data green(39,14,10) /   2.3460633589378256D-02  /
   data green(39,14,11) /   2.3326095497965760D-02  /
   data green(39,14,12) /   2.3181374185121482D-02  /
   data green(39,14,13) /   2.3027090806762580D-02  /
   data green(39,14,14) /   2.2863884205106655D-02  /
   data green(39,15, 0) /   2.3933441444908776D-02  /
   data green(39,15, 1) /   2.3926582115098675D-02  /
   data green(39,15, 2) /   2.3906039570193906D-02  /
   data green(39,15, 3) /   2.3871919633709059D-02  /
   data green(39,15, 4) /   2.3824396992398161D-02  /
   data green(39,15, 5) /   2.3763712717447850D-02  /
   data green(39,15, 6) /   2.3690170891450556D-02  /
   data green(39,15, 7) /   2.3604134433044641D-02  /
   data green(39,15, 8) /   2.3506020230047095D-02  /
   data green(39,15, 9) /   2.3396293706050380D-02  /
   data green(39,15,10) /   2.3275462954434313D-02  /
   data green(39,15,11) /   2.3144072577487254D-02  /
   data green(39,15,12) /   2.3002697367061967D-02  /
   data green(39,15,13) /   2.2851935957383832D-02  /
   data green(39,15,14) /   2.2692404570946439D-02  /
   data green(39,15,15) /   2.2524730965668052D-02  /
   data green(39,16, 0) /   2.3723551381556247D-02  /
   data green(39,16, 1) /   2.3716871225132927D-02  /
   data green(39,16, 2) /   2.3696864667190447D-02  /
   data green(39,16, 3) /   2.3663632962362426D-02  /
   data green(39,16, 4) /   2.3617343283866631D-02  /
   data green(39,16, 5) /   2.3558226392767424D-02  /
   data green(39,16, 6) /   2.3486573464812317D-02  /
   data green(39,16, 7) /   2.3402732159798878D-02  /
   data green(39,16, 8) /   2.3307102036059336D-02  /
   data green(39,16, 9) /   2.3200129425918613D-02  /
   data green(39,16,10) /   2.3082301896532820D-02  /
   data green(39,16,11) /   2.2954142424275556D-02  /
   data green(39,16,12) /   2.2816203409991704D-02  /
   data green(39,16,13) /   2.2669060657398168D-02  /
   data green(39,16,14) /   2.2513307428268136D-02  /
   data green(39,16,15) /   2.2349548676503052D-02  /
   data green(39,16,16) /   2.2178395549550624D-02  /
   data green(39,17, 0) /   2.3506092550903742D-02  /
   data green(39,17, 1) /   2.3499594685244152D-02  /
   data green(39,17, 2) /   2.3480133468028117D-02  /
   data green(39,17, 3) /   2.3447805589313690D-02  /
   data green(39,17, 4) /   2.3402770714016554D-02  /
   data green(39,17, 5) /   2.3345249296472453D-02  /
   data green(39,17, 6) /   2.3275519603532902D-02  /
   data green(39,17, 7) /   2.3193914024465068D-02  /
   data green(39,17, 8) /   2.3100814762285357D-02  /
   data green(39,17, 9) /   2.2996649013557590D-02  /
   data green(39,17,10) /   2.2881883751804343D-02  /
   data green(39,17,11) /   2.2757020233426566D-02  /
   data green(39,17,12) /   2.2622588344558048D-02  /
   data green(39,17,13) /   2.2479140902954255D-02  /
   data green(39,17,14) /   2.2327248021351623D-02  /
   data green(39,17,15) /   2.2167491628369149D-02  /
   data green(39,17,16) /   2.2000460230654943D-02  /
   data green(39,17,17) /   2.1826743986318987D-02  /
   data green(39,18, 0) /   2.3281880257586447D-02  /
   data green(39,18, 1) /   2.3275566816244962D-02  /
   data green(39,18, 2) /   2.3256657351779508D-02  /
   data green(39,18, 3) /   2.3225244022961541D-02  /
   data green(39,18, 4) /   2.3181479039289521D-02  /
   data green(39,18, 5) /   2.3125572617121143D-02  /
   data green(39,18, 6) /   2.3057790194122903D-02  /
   data green(39,18, 7) /   2.2978448973910305D-02  /
   data green(39,18, 8) /   2.2887913887879377D-02  /
   data green(39,18, 9) /   2.2786593072786323D-02  /
   data green(39,18,10) /   2.2674932970308710D-02  /
   data green(39,18,11) /   2.2553413158530090D-02  /
   data green(39,18,12) /   2.2422541025154903D-02  /
   data green(39,18,13) /   2.2282846388588913D-02  /
   data green(39,18,14) /   2.2134876166271293D-02  /
   data green(39,18,15) /   2.1979189180379697D-02  /
   data green(39,18,16) /   2.1816351179872084D-02  /
   data green(39,18,17) /   2.1646930145419123D-02  /
   data green(39,18,18) /   2.1471491930735772D-02  /
   data green(39,19, 0) /   2.3051710184581180D-02  /
   data green(39,19, 1) /   2.3045582380716759D-02  /
   data green(39,19, 2) /   2.3027228328474935D-02  /
   data green(39,19, 3) /   2.2996735714371368D-02  /
   data green(39,19, 4) /   2.2954249387732484D-02  /
   data green(39,19, 5) /   2.2899969453870447D-02  /
   data green(39,19, 6) /   2.2834148673898245D-02  /
   data green(39,19, 7) /   2.2757089236979684D-02  /
   data green(39,19, 8) /   2.2669138984754492D-02  /
   data green(39,19, 9) /   2.2570687178414438D-02  /
   data green(39,19,10) /   2.2462159906141414D-02  /
   data green(39,19,11) /   2.2344015232264276D-02  /
   data green(39,19,12) /   2.2216738189644938D-02  /
   data green(39,19,13) /   2.2080835713729624D-02  /
   data green(39,19,14) /   2.1936831610795963D-02  /
   data green(39,19,15) /   2.1785261644689361D-02  /
   data green(39,19,16) /   2.1626668816323778D-02  /
   data green(39,19,17) /   2.1461598898997923D-02  /
   data green(39,19,18) /   2.1290596280704087D-02  /
   data green(39,19,19) /   2.1114200152597878D-02  /
   data green(39,20, 0) /   2.2816353355356603D-02  /
   data green(39,20, 1) /   2.2810411550223630D-02  /
   data green(39,20, 2) /   2.2792614021495381D-02  /
   data green(39,20, 3) /   2.2763044064923178D-02  /
   data green(39,20, 4) /   2.2721839301907423D-02  /
   data green(39,20, 5) /   2.2669189904548376D-02  /
   data green(39,20, 6) /   2.2605336173966286D-02  /
   data green(39,20, 7) /   2.2530565531946226D-02  /
   data green(39,20, 8) /   2.2445208998784588D-02  /
   data green(39,20, 9) /   2.2349637240156546D-02  /
   data green(39,20,10) /   2.2244256272621322D-02  /
   data green(39,20,11) /   2.2129502920943359D-02  /
   data green(39,20,12) /   2.2005840120806672D-02  /
   data green(39,20,13) /   2.1873752157961963D-02  /
   data green(39,20,14) /   2.1733739929715355D-02  /
   data green(39,20,15) /   2.1586316307380529D-02  /
   data green(39,20,16) /   2.1432001669362676D-02  /
   data green(39,20,17) /   2.1271319664431812D-02  /
   data green(39,20,18) /   2.1104793253979875D-02  /
   data green(39,20,19) /   2.0932941071100625D-02  /
   data green(39,20,20) /   2.0756274123600267D-02  /
   data green(39,21, 0) /   2.2576551928815192D-02  /
   data green(39,21, 1) /   2.2570795703488959D-02  /
   data green(39,21, 2) /   2.2553553475467048D-02  /
   data green(39,21, 3) /   2.2524904250448365D-02  /
   data green(39,21, 4) /   2.2484978585855311D-02  /
   data green(39,21, 5) /   2.2433956942102793D-02  /
   data green(39,21, 6) /   2.2372067431888056D-02  /
   data green(39,21, 7) /   2.2299583022171922D-02  /
   data green(39,21, 8) /   2.2216818255277496D-02  /
   data green(39,21, 9) /   2.2124125564715658D-02  /
   data green(39,21,10) /   2.2021891267712545D-02  /
   data green(39,21,11) /   2.1910531319871661D-02  /
   data green(39,21,12) /   2.1790486918007165D-02  /
   data green(39,21,13) /   2.1662220035125680D-02  /
   data green(39,21,14) /   2.1526208967106518D-02  /
   data green(39,21,15) /   2.1382943964217672D-02  /
   data green(39,21,16) /   2.1232923012636698D-02  /
   data green(39,21,17) /   2.1076647822076060D-02  /
   data green(39,21,18) /   2.0914620065890989D-02  /
   data green(39,21,19) /   2.0747337910092686D-02  /
   data green(39,21,20) /   2.0575292857876539D-02  /
   data green(39,21,21) /   2.0398966926918491D-02  /
   data green(39,22, 0) /   2.2333015801221451D-02  /
   data green(39,22, 1) /   2.2327444029908442D-02  /
   data green(39,22, 2) /   2.2310753764565452D-02  /
   data green(39,22, 3) /   2.2283019837637360D-02  /
   data green(39,22, 4) /   2.2244365933089445D-02  /
   data green(39,22, 5) /   2.2194963057882246D-02  /
   data green(39,22, 6) /   2.2135027454204741D-02  /
   data green(39,22, 7) /   2.2064818002111155D-02  /
   data green(39,22, 8) /   2.1984633172958599D-02  /
   data green(39,22, 9) /   2.1894807602501874D-02  /
   data green(39,22,10) /   2.1795708358445964D-02  /
   data green(39,22,11) /   2.1687730980593904D-02  /
   data green(39,22,12) /   2.1571295372498595D-02  /
   data green(39,22,13) /   2.1446841621889302D-02  /
   data green(39,22,14) /   2.1314825823351209D-02  /
   data green(39,22,15) /   2.1175715971120931D-02  /
   data green(39,22,16) /   2.1029987982800185D-02  /
   data green(39,22,17) /   2.0878121906685795D-02  /
   data green(39,22,18) /   2.0720598356664913D-02  /
   data green(39,22,19) /   2.0557895209610239D-02  /
   data green(39,22,20) /   2.0390484591266987D-02  /
   data green(39,22,21) /   2.0218830168044478D-02  /
   data green(39,22,22) /   2.0043384754141729D-02  /
   data green(39,23, 0) /   2.2086419973567681D-02  /
   data green(39,23, 1) /   2.2081030897052437D-02  /
   data green(39,23, 2) /   2.2064887360480737D-02  /
   data green(39,23, 3) /   2.2038060152953336D-02  /
   data green(39,23, 4) /   2.2000666297197904D-02  /
   data green(39,23, 5) /   2.1952867635005899D-02  /
   data green(39,23, 6) /   2.1894868894056667D-02  /
   data green(39,23, 7) /   2.1826915281108666D-02  /
   data green(39,23, 8) /   2.1749289656348123D-02  /
   data green(39,23, 9) /   2.1662309351458887D-02  /
   data green(39,23,10) /   2.1566322699511534D-02  /
   data green(39,23,11) /   2.1461705347973571D-02  /
   data green(39,23,12) /   2.1348856427045593D-02  /
   data green(39,23,13) /   2.1228194644257993D-02  /
   data green(39,23,14) /   2.1100154373039293D-02  /
   data green(39,23,15) /   2.0965181798073818D-02  /
   data green(39,23,16) /   2.0823731174036308D-02  /
   data green(39,23,17) /   2.0676261247075684D-02  /
   data green(39,23,18) /   2.0523231880574963D-02  /
   data green(39,23,19) /   2.0365100918578159D-02  /
   data green(39,23,20) /   2.0202321312153403D-02  /
   data green(39,23,21) /   2.0035338526120199D-02  /
   data green(39,23,22) /   1.9864588236223529D-02  /
   data green(39,23,23) /   1.9690494320153618D-02  /
   data green(39,24, 0) /   2.1837402630912578D-02  /
   data green(39,24, 1) /   2.1832193928903801D-02  /
   data green(39,24, 2) /   2.1816590207420807D-02  /
   data green(39,24, 3) /   2.1790658352463944D-02  /
   data green(39,24, 4) /   2.1754508954886714D-02  /
   data green(39,24, 5) /   2.1708295003427705D-02  /
   data green(39,24, 6) /   2.1652210097594025D-02  /
   data green(39,24, 7) /   2.1586486221064252D-02  /
   data green(39,24, 8) /   2.1511391125211023D-02  /
   data green(39,24, 9) /   2.1427225379470762D-02  /
   data green(39,24,10) /   2.1334319150425673D-02  /
   data green(39,24,11) /   2.1233028774524848D-02  /
   data green(39,24,12) /   2.1123733190373372D-02  /
   data green(39,24,13) /   2.1006830295567030D-02  /
   data green(39,24,14) /   2.0882733290331159D-02  /
   data green(39,24,15) /   2.0751867065980107D-02  /
   data green(39,24,16) /   2.0614664690737777D-02  /
   data green(39,24,17) /   2.0471564039059333D-02  /
   data green(39,24,18) /   2.0323004603582794D-02  /
   data green(39,24,19) /   2.0169424521519382D-02  /
   data green(39,24,20) /   2.0011257839939081D-02  /
   data green(39,24,21) /   1.9848932037263279D-02  /
   data green(39,24,22) /   1.9682865811539857D-02  /
   data green(39,24,23) /   1.9513467139902430D-02  /
   data green(39,24,24) /   1.9341131608115014D-02  /
   data green(39,25, 0) /   2.1586563871893234D-02  /
   data green(39,25, 1) /   2.1581532733244396D-02  /
   data green(39,25, 2) /   2.1566460443196711D-02  /
   data green(39,25, 3) /   2.1541410132674854D-02  /
   data green(39,25, 4) /   2.1506486202974850D-02  /
   data green(39,25, 5) /   2.1461833119995492D-02  /
   data green(39,25, 6) /   2.1407633764597471D-02  /
   data green(39,25, 7) /   2.1344107375790268D-02  /
   data green(39,25, 8) /   2.1271507131562516D-02  /
   data green(39,25, 9) /   2.1190117418693009D-02  /
   data green(39,25,10) /   2.1100250847636632D-02  /
   data green(39,25,11) /   2.1002245071491298D-02  /
   data green(39,25,12) /   2.0896459469125616D-02  /
   data green(39,25,13) /   2.0783271751868306D-02  /
   data green(39,25,14) /   2.0663074550886496D-02  /
   data green(39,25,15) /   2.0536272038720903D-02  /
   data green(39,25,16) /   2.0403276633651959D-02  /
   data green(39,25,17) /   2.0264505829912838D-02  /
   data green(39,25,18) /   2.0120379190519181D-02  /
   data green(39,25,19) /   1.9971315532920808D-02  /
   data green(39,25,20) /   1.9817730331044323D-02  /
   data green(39,25,21) /   1.9660033350805448D-02  /
   data green(39,25,22) /   1.9498626530010776D-02  /
   data green(39,25,23) /   1.9333902107883733D-02  /
   data green(39,25,24) /   1.9166241004345479D-02  /
   data green(39,25,25) /   1.8996011444727453D-02  /
   data green(39,26, 0) /   2.1334465021500413D-02  /
   data green(39,26, 1) /   2.1329608211484240D-02  /
   data green(39,26, 2) /   2.1315057700295077D-02  /
   data green(39,26, 3) /   2.1290873017275357D-02  /
   data green(39,26, 4) /   2.1257152625634369D-02  /
   data green(39,26, 5) /   2.1214032811534594D-02  /
   data green(39,26, 6) /   2.1161686163399515D-02  /
   data green(39,26, 7) /   2.1100319674504861D-02  /
   data green(39,26, 8) /   2.1030172509276239D-02  /
   data green(39,26, 9) /   2.0951513479670975D-02  /
   data green(39,26,10) /   2.0864638282415925D-02  /
   data green(39,26,11) /   2.0769866550630935D-02  /
   data green(39,26,12) /   2.0667538774487126D-02  /
   data green(39,26,13) /   2.0558013145102440D-02  /
   data green(39,26,14) /   2.0441662373992028D-02  /
   data green(39,26,15) /   2.0318870537251701D-02  /
   data green(39,26,16) /   2.0190029989471357D-02  /
   data green(39,26,17) /   2.0055538387389982D-02  /
   data green(39,26,18) /   1.9915795857757386D-02  /
   data green(39,26,19) /   1.9771202337997094D-02  /
   data green(39,26,20) /   1.9622155112292539D-02  /
   data green(39,26,21) /   1.9469046559840466D-02  /
   data green(39,26,22) /   1.9312262126400085D-02  /
   data green(39,26,23) /   1.9152178525047865D-02  /
   data green(39,26,24) /   1.8989162167326705D-02  /
   data green(39,26,25) /   1.8823567821822658D-02  /
   data green(39,26,26) /   1.8655737493650341D-02  /
   data green(39,27, 0) /   2.1081628457890193D-02  /
   data green(39,27, 1) /   2.1076942381898386D-02  /
   data green(39,27, 2) /   2.1062902918473937D-02  /
   data green(39,27, 3) /   2.1039566152274838D-02  /
   data green(39,27, 4) /   2.1007024865652767D-02  /
   data green(39,27, 5) /   2.0965407516367229D-02  /
   data green(39,27, 6) /   2.0914876837461745D-02  /
   data green(39,27, 7) /   2.0855628089039912D-02  /
   data green(39,27, 8) /   2.0787886998345745D-02  /
   data green(39,27, 9) /   2.0711907429979016D-02  /
   data green(39,27,10) /   2.0627968832126305D-02  /
   data green(39,27,11) /   2.0536373507289971D-02  /
   data green(39,27,12) /   2.0437443757141481D-02  /
   data green(39,27,13) /   2.0331518950872617D-02  /
   data green(39,27,14) /   2.0218952564872775D-02  /
   data green(39,27,15) /   2.0100109238880171D-02  /
   data green(39,27,16) /   1.9975361890122126D-02  /
   data green(39,27,17) /   1.9845088922581400D-02  /
   data green(39,27,18) /   1.9709671563614511D-02  /
   data green(39,27,19) /   1.9569491354912871D-02  /
   data green(39,27,20) /   1.9424927819437731D-02  /
   data green(39,27,21) /   1.9276356320650338D-02  /
   data green(39,27,22) /   1.9124146125253788D-02  /
   data green(39,27,23) /   1.8968658675886152D-02  /
   data green(39,27,24) /   1.8810246075850433D-02  /
   data green(39,27,25) /   1.8649249784103459D-02  /
   data green(39,27,26) /   1.8485999515390639D-02  /
   data green(39,27,27) /   1.8320812337623622D-02  /
   data green(39,28, 0) /   2.0828537884023814D-02  /
   data green(39,28, 1) /   2.0824018647237541D-02  /
   data green(39,28, 2) /   2.0810478600368928D-02  /
   data green(39,28, 3) /   2.0787970541876657D-02  /
   data green(39,28, 4) /   2.0756581833250076D-02  /
   data green(39,28, 5) /   2.0716433459295230D-02  /
   data green(39,28, 6) /   2.0667678740418784D-02  /
   data green(39,28, 7) /   2.0610501723621274D-02  /
   data green(39,28, 8) /   2.0545115284938215D-02  /
   data green(39,28, 9) /   2.0471758981005336D-02  /
   data green(39,28,10) /   2.0390696691148458D-02  /
   data green(39,28,11) /   2.0302214093843189D-02  /
   data green(39,28,12) /   2.0206616022541535D-02  /
   data green(39,28,13) /   2.0104223745768854D-02  /
   data green(39,28,14) /   1.9995372215143484D-02  /
   data green(39,28,15) /   1.9880407322694867D-02  /
   data green(39,28,16) /   1.9759683205711864D-02  /
   data green(39,28,17) /   1.9633559633520020D-02  /
   data green(39,28,18) /   1.9502399506250032D-02  /
   data green(39,28,19) /   1.9366566491004888D-02  /
   data green(39,28,20) /   1.9226422816034414D-02  /
   data green(39,28,21) /   1.9082327238743611D-02  /
   data green(39,28,22) /   1.8934633198732211D-02  /
   data green(39,28,23) /   1.8783687162702529D-02  /
   data green(39,28,24) /   1.8629827164069461D-02  /
   data green(39,28,25) /   1.8473381536524891D-02  /
   data green(39,28,26) /   1.8314667837689900D-02  /
   data green(39,28,27) /   1.8153991956348561D-02  /
   data green(39,28,28) /   1.7991647394600910D-02  /
   data green(39,29, 0) /   2.0575638976821020D-02  /
   data green(39,29, 1) /   2.0571282439552495D-02  /
   data green(39,29, 2) /   2.0558229443412564D-02  /
   data green(39,29, 3) /   2.0536529659156501D-02  /
   data green(39,29, 4) /   2.0506265287573846D-02  /
   data green(39,29, 5) /   2.0467550196503197D-02  /
   data green(39,29, 6) /   2.0420528737641226D-02  /
   data green(39,29, 7) /   2.0365374267109132D-02  /
   data green(39,29, 8) /   2.0302287399181586D-02  /
   data green(39,29, 9) /   2.0231494027069283D-02  /
   data green(39,29,10) /   2.0153243148064302D-02  /
   data green(39,29,11) /   2.0067804532646456D-02  /
   data green(39,29,12) /   1.9975466278293934D-02  /
   data green(39,29,13) /   1.9876532288777341D-02  /
   data green(39,29,14) /   1.9771319718717949D-02  /
   data green(39,29,15) /   1.9660156421267970D-02  /
   data green(39,29,16) /   1.9543378434059983D-02  /
   data green(39,29,17) /   1.9421327535226747D-02  /
   data green(39,29,18) /   1.9294348897474521D-02  /
   data green(39,29,19) /   1.9162788864063869D-02  /
   data green(39,29,20) /   1.9026992866267439D-02  /
   data green(39,29,21) /   1.8887303497576707D-02  /
   data green(39,29,22) /   1.8744058755743503D-02  /
   data green(39,29,23) /   1.8597590459772037D-02  /
   data green(39,29,24) /   1.8448222845307439D-02  /
   data green(39,29,25) /   1.8296271338556491D-02  /
   data green(39,29,26) /   1.8142041505968692D-02  /
   data green(39,29,27) /   1.7985828174419603D-02  /
   data green(39,29,28) /   1.7827914714580718D-02  /
   data green(39,29,29) /   1.7668572478519558D-02  /
   data green(39,30, 0) /   2.0323340349841842D-02  /
   data green(39,30, 1) /   2.0319142178382005D-02  /
   data green(39,30, 2) /   2.0306563284617149D-02  /
   data green(39,30, 3) /   2.0285650368759064D-02  /
   data green(39,30, 4) /   2.0256480729003681D-02  /
   data green(39,30, 5) /   2.0219161469676263D-02  /
   data green(39,30, 6) /   2.0173828415004407D-02  /
   data green(39,30, 7) /   2.0120644749993172D-02  /
   data green(39,30, 8) /   2.0059799414784415D-02  /
   data green(39,30, 9) /   1.9991505282951592D-02  /
   data green(39,30,10) /   1.9915997157313309D-02  /
   data green(39,30,11) /   1.9833529618985331D-02  /
   data green(39,30,12) /   1.9744374766516092D-02  /
   data green(39,30,13) /   1.9648819882091259D-02  /
   data green(39,30,14) /   1.9547165061009168D-02  /
   data green(39,30,15) /   1.9439720839015354D-02  /
   data green(39,30,16) /   1.9326805849755067D-02  /
   data green(39,30,17) /   1.9208744541691105D-02  /
   data green(39,30,18) /   1.9085864980481412D-02  /
   data green(39,30,19) /   1.8958496759156715D-02  /
   data green(39,30,20) /   1.8826969034622696D-02  /
   data green(39,30,21) /   1.8691608705158087D-02  /
   data green(39,30,22) /   1.8552738739803537D-02  /
   data green(39,30,23) /   1.8410676666931559D-02  /
   data green(39,30,24) /   1.8265733225930822D-02  /
   data green(39,30,25) /   1.8118211182889352D-02  /
   data green(39,30,26) /   1.7968404308458284D-02  /
   data green(39,30,27) /   1.7816596513745311D-02  /
   data green(39,30,28) /   1.7663061138132329D-02  /
   data green(39,30,29) /   1.7508060381330008D-02  /
   data green(39,30,30) /   1.7351844870758144D-02  /
   data green(39,31, 0) /   2.0072014769880190D-02  /
   data green(39,31, 1) /   2.0067970482800526D-02  /
   data green(39,31, 2) /   2.0055852299057181D-02  /
   data green(39,31, 3) /   2.0035704103010650D-02  /
   data green(39,31, 4) /   2.0007598544441899D-02  /
   data green(39,31, 5) /   1.9971636312496748D-02  /
   data green(39,31, 6) /   1.9927945139214968D-02  /
   data green(39,31, 7) /   1.9876678551869927D-02  /
   data green(39,31, 8) /   1.9818014397765218D-02  /
   data green(39,31, 9) /   1.9752153168822069D-02  /
   data green(39,31,10) /   1.9679316156155546D-02  /
   data green(39,31,11) /   1.9599743466826143D-02  /
   data green(39,31,12) /   1.9513691936049304D-02  /
   data green(39,31,13) /   1.9421432968367607D-02  /
   data green(39,31,14) /   1.9323250340689397D-02  /
   data green(39,31,15) /   1.9219437998751355D-02  /
   data green(39,31,16) /   1.9110297876568958D-02  /
   data green(39,31,17) /   1.8996137765912652D-02  /
   data green(39,31,18) /   1.8877269259910023D-02  /
   data green(39,31,19) /   1.8754005791648674D-02  /
   data green(39,31,20) /   1.8626660785263264D-02  /
   data green(39,31,21) /   1.8495545933543602D-02  /
   data green(39,31,22) /   1.8360969612702126D-02  /
   data green(39,31,23) /   1.8223235441672554D-02  /
   data green(39,31,24) /   1.8082640990250629D-02  /
   data green(39,31,25) /   1.7939476637586011D-02  /
   data green(39,31,26) /   1.7794024580030463D-02  /
   data green(39,31,27) /   1.7646557985166404D-02  /
   data green(39,31,28) /   1.7497340286991448D-02  /
   data green(39,31,29) /   1.7346624615717923D-02  /
   data green(39,31,30) /   1.7194653354450494D-02  /
   data green(39,31,31) /   1.7041657814112262D-02  /
   data green(39,32, 0) /   1.9822000572906399D-02  /
   data green(39,32, 1) /   1.9818105582856439D-02  /
   data green(39,32, 2) /   1.9806434397861651D-02  /
   data green(39,32, 3) /   1.9787028237722128D-02  /
   data green(39,32, 4) /   1.9759955352508705D-02  /
   data green(39,32, 5) /   1.9725310357257334D-02  /
   data green(39,32, 6) /   1.9683213318420646D-02  /
   data green(39,32, 7) /   1.9633808609275348D-02  /
   data green(39,32, 8) /   1.9577263555462358D-02  /
   data green(39,32, 9) /   1.9513766895173590D-02  /
   data green(39,32,10) /   1.9443527081114857D-02  /
   data green(39,32,11) /   1.9366770453220550D-02  /
   data green(39,32,12) /   1.9283739312154304D-02  /
   data green(39,32,13) /   1.9194689923915098D-02  /
   data green(39,32,14) /   1.9099890485420493D-02  /
   data green(39,32,15) /   1.8999619079823633D-02  /
   data green(39,32,16) /   1.8894161648620962D-02  /
   data green(39,32,17) /   1.8783810005422001D-02  /
   data green(39,32,18) /   1.8668859913685235D-02  /
   data green(39,32,19) /   1.8549609247883212D-02  /
   data green(39,32,20) /   1.8426356254552599D-02  /
   data green(39,32,21) /   1.8299397926608764D-02  /
   data green(39,32,22) /   1.8169028501251822D-02  /
   data green(39,32,23) /   1.8035538088839179D-02  /
   data green(39,32,24) /   1.7899211437314062D-02  /
   data green(39,32,25) /   1.7760326834211995D-02  /
   data green(39,32,26) /   1.7619155145954288D-02  /
   data green(39,32,27) /   1.7475958992104940D-02  /
   data green(39,32,28) /   1.7330992050526371D-02  /
   data green(39,32,29) /   1.7184498487922601D-02  /
   data green(39,32,30) /   1.7036712509098162D-02  /
   data green(39,32,31) /   1.6887858017374659D-02  /
   data green(39,32,32) /   1.6738148377972960D-02  /
   data green(39,33, 0) /   1.9573603230238046D-02  /
   data green(39,33, 1) /   1.9569852881355244D-02  /
   data green(39,33, 2) /   1.9558614776894648D-02  /
   data green(39,33, 3) /   1.9539927619232393D-02  /
   data green(39,33, 4) /   1.9513855500710067D-02  /
   data green(39,33, 5) /   1.9480487294319865D-02  /
   data green(39,33, 6) /   1.9439935816631230D-02  /
   data green(39,33, 7) /   1.9392336778333356D-02  /
   data green(39,33, 8) /   1.9337847541345418D-02  /
   data green(39,33, 9) /   1.9276645704462837D-02  /
   data green(39,33,10) /   1.9208927541893520D-02  /
   data green(39,33,11) /   1.9134906320747064D-02  /
   data green(39,33,12) /   1.9054810524556369D-02  /
   data green(39,33,13) /   1.8968882010242223D-02  /
   data green(39,33,14) /   1.8877374125612200D-02  /
   data green(39,33,15) /   1.8780549813568401D-02  /
   data green(39,33,16) /   1.8678679727755491D-02  /
   data green(39,33,17) /   1.8572040382494437D-02  /
   data green(39,33,18) /   1.8460912357609494D-02  /
   data green(39,33,19) /   1.8345578576258560D-02  /
   data green(39,33,20) /   1.8226322671215836D-02  /
   data green(39,33,21) /   1.8103427452315610D-02  /
   data green(39,33,22) /   1.7977173485030095D-02  /
   data green(39,33,23) /   1.7847837787491681D-02  /
   data green(39,33,24) /   1.7715692650741768D-02  /
   data green(39,33,25) /   1.7581004584640191D-02  /
   data green(39,33,26) /   1.7444033389740196D-02  /
   data green(39,33,27) /   1.7305031353544196D-02  /
   data green(39,33,28) /   1.7164242567922733D-02  /
   data green(39,33,29) /   1.7021902363105048D-02  /
   data green(39,33,30) /   1.6878236852531023D-02  /
   data green(39,33,31) /   1.6733462581982179D-02  /
   data green(39,33,32) /   1.6587786275765307D-02  /
   data green(39,33,33) /   1.6441404672290848D-02  /
   data green(39,34, 0) /   1.9327097021403547D-02  /
   data green(39,34, 1) /   1.9323486622509791D-02  /
   data green(39,34, 2) /   1.9312667572820009D-02  /
   data green(39,34, 3) /   1.9294676199676897D-02  /
   data green(39,34, 4) /   1.9269572671965482D-02  /
   data green(39,34, 5) /   1.9237440442322326D-02  /
   data green(39,34, 6) /   1.9198385480480380D-02  /
   data green(39,34, 7) /   1.9152535311485121D-02  /
   data green(39,34, 8) /   1.9100037875732456D-02  /
   data green(39,34, 9) /   1.9041060230502806D-02  /
   data green(39,34,10) /   1.8975787114839212D-02  /
   data green(39,34,11) /   1.8904419401195801D-02  /
   data green(39,34,12) /   1.8827172458252998D-02  /
   data green(39,34,13) /   1.8744274449659579D-02  /
   data green(39,34,14) /   1.8655964593248538D-02  /
   data green(39,34,15) /   1.8562491404527003D-02  /
   data green(39,34,16) /   1.8464110947019458D-02  /
   data green(39,34,17) /   1.8361085110421754D-02  /
   data green(39,34,18) /   1.8253679935575633D-02  /
   data green(39,34,19) /   1.8142164003084233D-02  /
   data green(39,34,20) /   1.8026806900037277D-02  /
   data green(39,34,21) /   1.7907877776879411D-02  /
   data green(39,34,22) /   1.7785644004007631D-02  /
   data green(39,34,23) /   1.7660369935287360D-02  /
   data green(39,34,24) /   1.7532315783386206D-02  /
   data green(39,34,25) /   1.7401736609683562D-02  /
   data green(39,34,26) /   1.7268881429557559D-02  /
   data green(39,34,27) /   1.7133992432100870D-02  /
   data green(39,34,28) /   1.6997304311790324D-02  /
   data green(39,34,29) /   1.6859043708335544D-02  /
   data green(39,34,30) /   1.6719428749860422D-02  /
   data green(39,34,31) /   1.6578668693719063D-02  /
   data green(39,34,32) /   1.6436963658602873D-02  /
   data green(39,34,33) /   1.6294504441143387D-02  /
   data green(39,34,34) /   1.6151472409936120D-02  /
   data green(39,35, 0) /   1.9082726775696487D-02  /
   data green(39,35, 1) /   1.9079251629498556D-02  /
   data green(39,35, 2) /   1.9068837588719387D-02  /
   data green(39,35, 3) /   1.9051518742865016D-02  /
   data green(39,35, 4) /   1.9027351563180105D-02  /
   data green(39,35, 5) /   1.8996414392205112D-02  /
   data green(39,35, 6) /   1.8958806741871181D-02  /
   data green(39,35, 7) /   1.8914648412399218D-02  /
   data green(39,35, 8) /   1.8864078447155135D-02  /
   data green(39,35, 9) /   1.8807253941073115D-02  /
   data green(39,35,10) /   1.8744348722234710D-02  /
   data green(39,35,11) /   1.8675551927647759D-02  /
   data green(39,35,12) /   1.8601066495189389D-02  /
   data green(39,35,13) /   1.8521107594062738D-02  /
   data green(39,35,14) /   1.8435901015990963D-02  /
   data green(39,35,15) /   1.8345681548769335D-02  /
   data green(39,35,16) /   1.8250691352769084D-02  /
   data green(39,35,17) /   1.8151178359594497D-02  /
   data green(39,35,18) /   1.8047394710404639D-02  /
   data green(39,35,19) /   1.7939595249494469D-02  /
   data green(39,35,20) /   1.7828036086656934D-02  /
   data green(39,35,21) /   1.7712973239686063D-02  /
   data green(39,35,22) /   1.7594661366196522D-02  /
   data green(39,35,23) /   1.7473352591781580D-02  /
   data green(39,35,24) /   1.7349295439461506D-02  /
   data green(39,35,25) /   1.7222733863425696D-02  /
   data green(39,35,26) /   1.7093906388279109D-02  /
   data green(39,35,27) /   1.6963045353386837D-02  /
   data green(39,35,28) /   1.6830376260488139D-02  /
   data green(39,35,29) /   1.6696117221527581D-02  /
   data green(39,35,30) /   1.6560478502628744D-02  /
   data green(39,35,31) /   1.6423662159309863D-02  /
   data green(39,35,32) /   1.6285861757401789D-02  /
   data green(39,35,33) /   1.6147262173663534D-02  /
   data green(39,35,34) /   1.6008039469785407D-02  /
   data green(39,35,35) /   1.5868360833305759D-02  /
   data green(39,36, 0) /   1.8840709649752593D-02  /
   data green(39,36, 1) /   1.8837365078293620D-02  /
   data green(39,36, 2) /   1.8827342056715317D-02  /
   data green(39,36, 3) /   1.8810672568367767D-02  /
   data green(39,36, 4) /   1.8787409603673710D-02  /
   data green(39,36, 5) /   1.8757626693142387D-02  /
   data green(39,36, 6) /   1.8721417264928739D-02  /
   data green(39,36, 7) /   1.8678893837887744D-02  /
   data green(39,36, 8) /   1.8630187063663848D-02  /
   data green(39,36, 9) /   1.8575444633573409D-02  /
   data green(39,36,10) /   1.8514830067834068D-02  /
   data green(39,36,11) /   1.8448521406035492D-02  /
   data green(39,36,12) /   1.8376709818614230D-02  /
   data green(39,36,13) /   1.8299598159494578D-02  /
   data green(39,36,14) /   1.8217399480000025D-02  /
   data green(39,36,15) /   1.8130335523661376D-02  /
   data green(39,36,16) /   1.8038635220685384D-02  /
   data green(39,36,17) /   1.7942533199657745D-02  /
   data green(39,36,18) /   1.7842268332590224D-02  /
   data green(39,36,19) /   1.7738082327746946D-02  /
   data green(39,36,20) /   1.7630218382860545D-02  /
   data green(39,36,21) /   1.7518919909433966D-02  /
   data green(39,36,22) /   1.7404429336876068D-02  /
   data green(39,36,23) /   1.7286987003289157D-02  /
   data green(39,36,24) /   1.7166830137857744D-02  /
   data green(39,36,25) /   1.7044191938019464D-02  /
   data green(39,36,26) /   1.6919300742958759D-02  /
   data green(39,36,27) /   1.6792379303475706D-02  /
   data green(39,36,28) /   1.6663644146959831D-02  /
   data green(39,36,29) /   1.6533305035051359D-02  /
   data green(39,36,30) /   1.6401564510601672D-02  /
   data green(39,36,31) /   1.6268617529749296D-02  /
   data green(39,36,32) /   1.6134651174299290D-02  /
   data green(39,36,33) /   1.5999844439124049D-02  /
   data green(39,36,34) /   1.5864368088978568D-02  /
   data green(39,36,35) /   1.5728384578929024D-02  /
   data green(39,36,36) /   1.5592048032515966D-02  /
   data green(39,37, 0) /   1.8601236913509837D-02  /
   data green(39,37, 1) /   1.8598018280137431D-02  /
   data green(39,37, 2) /   1.8588372410034582D-02  /
   data green(39,37, 3) /   1.8572329306346450D-02  /
   data green(39,37, 4) /   1.8549938686132476D-02  /
   data green(39,37, 5) /   1.8521269553222287D-02  /
   data green(39,37, 6) /   1.8486409610324557D-02  /
   data green(39,37, 7) /   1.8445464520161173D-02  /
   data green(39,37, 8) /   1.8398557027722102D-02  /
   data green(39,37, 9) /   1.8345825957735783D-02  /
   data green(39,37,10) /   1.8287425103080311D-02  /
   data green(39,37,11) /   1.8223522021092509D-02  /
   data green(39,37,12) /   1.8154296755549360D-02  /
   data green(39,37,13) /   1.8079940502498799D-02  /
   data green(39,37,14) /   1.8000654238117697D-02  /
   data green(39,37,15) /   1.7916647326398152D-02  /
   data green(39,37,16) /   1.7828136123745345D-02  /
   data green(39,37,17) /   1.7735342596554338D-02  /
   data green(39,37,18) /   1.7638492966568418D-02  /
   data green(39,37,19) /   1.7537816397361058D-02  /
   data green(39,37,20) /   1.7433543733680133D-02  /
   data green(39,37,21) /   1.7325906303699849D-02  /
   data green(39,37,22) /   1.7215134792492839D-02  /
   data green(39,37,23) /   1.7101458193306628D-02  /
   data green(39,37,24) /   1.6985102841546182D-02  /
   data green(39,37,25) /   1.6866291534761105D-02  /
   data green(39,37,26) /   1.6745242740439984D-02  /
   data green(39,37,27) /   1.6622169892046422D-02  /
   data green(39,37,28) /   1.6497280772506284D-02  /
   data green(39,37,29) /   1.6370776983283049D-02  /
   data green(39,37,30) /   1.6242853496260330D-02  /
   data green(39,37,31) /   1.6113698284889153D-02  /
   data green(39,37,32) /   1.5983492030444291D-02  /
   data green(39,37,33) /   1.5852407898765117D-02  /
   data green(39,37,34) /   1.5720611382517992D-02  /
   data green(39,37,35) /   1.5588260203800292D-02  /
   data green(39,37,36) /   1.5455504271795839D-02  /
   data green(39,37,37) /   1.5322485690176180D-02  /
   data green(39,38, 0) /   1.8364475721562051D-02  /
   data green(39,38, 1) /   1.8361378449688458D-02  /
   data green(39,38, 2) /   1.8352096041561097D-02  /
   data green(39,38, 3) /   1.8336656640220880D-02  /
   data green(39,38, 4) /   1.8315106887254424D-02  /
   data green(39,38, 5) /   1.8287511532142982D-02  /
   data green(39,38, 6) /   1.8253952894364659D-02  /
   data green(39,38, 7) /   1.8214530186968822D-02  /
   data green(39,38, 8) /   1.8169358712426361D-02  /
   data green(39,38, 9) /   1.8118568943359486D-02  /
   data green(39,38,10) /   1.8062305502234401D-02  /
   data green(39,38,11) /   1.8000726055230353D-02  /
   data green(39,38,12) /   1.7934000136264741D-02  /
   data green(39,38,13) /   1.7862307917555196D-02  /
   data green(39,38,14) /   1.7785838943144554D-02  /
   data green(39,38,15) /   1.7704790841525313D-02  /
   data green(39,38,16) /   1.7619368032904412D-02  /
   data green(39,38,17) /   1.7529780445785424D-02  /
   data green(39,38,18) /   1.7436242256454598D-02  /
   data green(39,38,19) /   1.7338970663685820D-02  /
   data green(39,38,20) /   1.7238184709572370D-02  /
   data green(39,38,21) /   1.7134104155899173D-02  /
   data green(39,38,22) /   1.7026948423928351D-02  /
   data green(39,38,23) /   1.6916935603927225D-02  /
   data green(39,38,24) /   1.6804281539254524D-02  /
   data green(39,38,25) /   1.6689198988370846D-02  /
   data green(39,38,26) /   1.6571896866777191D-02  /
   data green(39,38,27) /   1.6452579569630772D-02  /
   data green(39,38,28) /   1.6331446374656344D-02  /
   data green(39,38,29) /   1.6208690923970797D-02  /
   data green(39,38,30) /   1.6084500782575326D-02  /
   data green(39,38,31) /   1.5959057070543212D-02  /
   data green(39,38,32) /   1.5832534165338023D-02  /
   data green(39,38,33) /   1.5705099470232665D-02  /
   data green(39,38,34) /   1.5576913244454739D-02  /
   data green(39,38,35) /   1.5448128490449076D-02  /
   data green(39,38,36) /   1.5318890893512271D-02  /
   data green(39,38,37) /   1.5189338809006645D-02  /
   data green(39,38,38) /   1.5059603292388877D-02  /
   data green(39,39, 0) /   1.8130570851147016D-02  /
   data green(39,39, 1) /   1.8127590440078681D-02  /
   data green(39,39, 2) /   1.8118658030127695D-02  /
   data green(39,39, 3) /   1.8103800018439028D-02  /
   data green(39,39, 4) /   1.8083060159370067D-02  /
   data green(39,39, 5) /   1.8056499207237724D-02  /
   data green(39,39, 6) /   1.8024194424199128D-02  /
   data green(39,39, 7) /   1.7986238961045502D-02  /
   data green(39,39, 8) /   1.7942741120555970D-02  /
   data green(39,39, 9) /   1.7893823514681471D-02  /
   data green(39,39,10) /   1.7839622128168291D-02  /
   data green(39,39,11) /   1.7780285302267308D-02  /
   data green(39,39,12) /   1.7715972652890818D-02  /
   data green(39,39,13) /   1.7646853937973460D-02  /
   data green(39,39,14) /   1.7573107888873989D-02  /
   data green(39,39,15) /   1.7494921020437201D-02  /
   data green(39,39,16) /   1.7412486433844309D-02  /
   data green(39,39,17) /   1.7326002625648242D-02  /
   data green(39,39,18) /   1.7235672315451523D-02  /
   data green(39,39,19) /   1.7141701303579455D-02  /
   data green(39,39,20) /   1.7044297368868991D-02  /
   data green(39,39,21) /   1.6943669215375771D-02  /
   data green(39,39,22) /   1.6840025475435207D-02  /
   data green(39,39,23) /   1.6733573775135343D-02  /
   data green(39,39,24) /   1.6624519866901478D-02  /
   data green(39,39,25) /   1.6513066832582596D-02  /
   data green(39,39,26) /   1.6399414359192354D-02  /
   data green(39,39,27) /   1.6283758088309487D-02  /
   data green(39,39,28) /   1.6166289039100164D-02  /
   data green(39,39,29) /   1.6047193103994908D-02  /
   data green(39,39,30) /   1.5926650615243457D-02  /
   data green(39,39,31) /   1.5804835979880950D-02  /
   data green(39,39,32) /   1.5681917380068862D-02  /
   data green(39,39,33) /   1.5558056535318139D-02  /
   data green(39,39,34) /   1.5433408522755310D-02  /
   data green(39,39,35) /   1.5308121651344988D-02  /
   data green(39,39,36) /   1.5182337385826463D-02  /
   data green(39,39,37) /   1.5056190316048258D-02  /
   data green(39,39,38) /   1.4929808167380857D-02  /
   data green(39,39,39) /   1.4803311847946110D-02  /
   data green(40, 0, 0) /   2.5003912871419582D-02  /
   data green(40, 1, 0) /   2.4996088116604292D-02  /
   data green(40, 1, 1) /   2.4988270726360340D-02  /
   data green(40, 2, 0) /   2.4972658039664742D-02  /
   data green(40, 2, 1) /   2.4964862673556930D-02  /
   data green(40, 2, 2) /   2.4941520485353079D-02  /
   data green(40, 3, 0) /   2.4933754507336099D-02  /
   data green(40, 3, 1) /   2.4925995617737003D-02  /
   data green(40, 3, 2) /   2.4902762515768038D-02  /
   data green(40, 3, 3) /   2.4864185219846610D-02  /
   data green(40, 4, 0) /   2.4879594999813896D-02  /
   data green(40, 4, 1) /   2.4871886698852781D-02  /
   data green(40, 4, 2) /   2.4848804888517508D-02  /
   data green(40, 4, 3) /   2.4810478174811641D-02  /
   data green(40, 4, 4) /   2.4757118682891505D-02  /
   data green(40, 5, 0) /   2.4810479238637183D-02  /
   data green(40, 5, 1) /   2.4802835172551402D-02  /
   data green(40, 5, 2) /   2.4779945467685773D-02  /
   data green(40, 5, 3) /   2.4741936945678611D-02  /
   data green(40, 5, 4) /   2.4689018800904978D-02  /
   data green(40, 5, 5) /   2.4621479407272738D-02  /
   data green(40, 6, 0) /   2.4726784610656718D-02  /
   data green(40, 6, 1) /   2.4719217844095712D-02  /
   data green(40, 6, 2) /   2.4696559321335489D-02  /
   data green(40, 6, 3) /   2.4658933730176965D-02  /
   data green(40, 6, 4) /   2.4606546759884135D-02  /
   data green(40, 6, 5) /   2.4539681982800991D-02  /
   data green(40, 6, 6) /   2.4458696619600827D-02  /
   data green(40, 7, 0) /   2.4628960524173578D-02  /
   data green(40, 7, 1) /   2.4621483436623959D-02  /
   data green(40, 7, 2) /   2.4599093125929474D-02  /
   data green(40, 7, 3) /   2.4561911822775843D-02  /
   data green(40, 7, 4) /   2.4510141179840954D-02  /
   data green(40, 7, 5) /   2.4444059238932746D-02  /
   data green(40, 7, 6) /   2.4364016311341483D-02  /
   data green(40, 7, 7) /   2.4270429890293026D-02  /
   data green(40, 8, 0) /   2.4517521860334100D-02  /
   data green(40, 8, 1) /   2.4510146056386906D-02  /
   data green(40, 8, 2) /   2.4488058672866880D-02  /
   data green(40, 8, 3) /   2.4451379189303627D-02  /
   data green(40, 8, 4) /   2.4400304738158629D-02  /
   data green(40, 8, 5) /   2.4335107166901650D-02  /
   data green(40, 8, 6) /   2.4256129046187763D-02  /
   data green(40, 8, 7) /   2.4163778738299630D-02  /
   data green(40, 8, 8) /   2.4058524663072026D-02  /
   data green(40, 9, 0) /   2.4393041702196101D-02  /
   data green(40, 9, 1) /   2.4385777936886132D-02  /
   data green(40, 9, 2) /   2.4364025657596201D-02  /
   data green(40, 9, 3) /   2.4327901330220333D-02  /
   data green(40, 9, 4) /   2.4277597135036848D-02  /
   data green(40, 9, 5) /   2.4213378131779227D-02  /
   data green(40, 9, 6) /   2.4135578406538689D-02  /
   data green(40, 9, 7) /   2.4044596309589251D-02  /
   data green(40, 9, 8) /   2.3940888915349092D-02  /
   data green(40, 9, 9) /   2.3824965851926228D-02  /
   data green(40,10, 0) /   2.4256143535046092D-02  /
   data green(40,10, 1) /   2.4249001654987584D-02  /
   data green(40,10, 2) /   2.4227613942870499D-02  /
   data green(40,10, 3) /   2.4192093621625329D-02  /
   data green(40,10, 4) /   2.4142627541665990D-02  /
   data green(40,10, 5) /   2.4079473455604766D-02  /
   data green(40,10, 6) /   2.4002956312902602D-02  /
   data green(40,10, 7) /   2.3913463678198590D-02  /
   data green(40,10, 8) /   2.3811440398201062D-02  /
   data green(40,10, 9) /   2.3697382657616731D-02  /
   data green(40,10,10) /   2.3571831574231546D-02  /
   data green(40,11, 0) /   2.4107493114583295D-02  /
   data green(40,11, 1) /   2.4100482015129890D-02  /
   data green(40,11, 2) /   2.4079485490784522D-02  /
   data green(40,11, 3) /   2.4044613326204899D-02  /
   data green(40,11, 4) /   2.3996046719985155D-02  /
   data green(40,11, 5) /   2.3934035674252199D-02  /
   data green(40,11, 6) /   2.3858895444194780D-02  /
   data green(40,11, 7) /   2.3771002145729671D-02  /
   data green(40,11, 8) /   2.3670787639621597D-02  /
   data green(40,11, 9) /   2.3558733825293363D-02  /
   data green(40,11,10) /   2.3435366486893108D-02  /
   data green(40,11,11) /   2.3301248837875686D-02  /
   data green(40,12, 0) /   2.3947790195008131D-02  /
   data green(40,12, 1) /   2.3940917793200501D-02  /
   data green(40,12, 2) /   2.3920336153782484D-02  /
   data green(40,12, 3) /   2.3886151462027112D-02  /
   data green(40,12, 4) /   2.3838538998795646D-02  /
   data green(40,12, 5) /   2.3777740648865613D-02  /
   data green(40,12, 6) /   2.3704061510610920D-02  /
   data green(40,12, 7) /   2.3617865699583240D-02  /
   data green(40,12, 8) /   2.3519571457600853D-02  /
   data green(40,12, 9) /   2.3409645693170590D-02  /
   data green(40,12,10) /   2.3288598088072002D-02  /
   data green(40,12,11) /   2.3156974908660314D-02  /
   data green(40,12,12) /   2.3015352659117558D-02  /
   data green(40,13, 0) /   2.3777760297769950D-02  /
   data green(40,13, 1) /   2.3771033520421078D-02  /
   data green(40,13, 2) /   2.3750887503743517D-02  /
   data green(40,13, 3) /   2.3717424706256925D-02  /
   data green(40,13, 4) /   2.3670814280499990D-02  /
   data green(40,13, 5) /   2.3611289702608319D-02  /
   data green(40,13, 6) /   2.3539145545629539D-02  /
   data green(40,13, 7) /   2.3454733483426277D-02  /
   data green(40,13, 8) /   2.3358457630002234D-02  /
   data green(40,13, 9) /   2.3250769332587892D-02  /
   data green(40,13,10) /   2.3132161545480534D-02  /
   data green(40,13,11) /   2.3003162915379832D-02  /
   data green(40,13,12) /   2.2864331707985661D-02  /
   data green(40,13,13) /   2.2716249700366665D-02  /
   data green(40,14, 0) /   2.3598146684937489D-02  /
   data green(40,14, 1) /   2.3591571470855922D-02  /
   data green(40,14, 2) /   2.3571878861714416D-02  /
   data green(40,14, 3) /   2.3539167494637643D-02  /
   data green(40,14, 4) /   2.3493600236936849D-02  /
   data green(40,14, 5) /   2.3435401938201678D-02  /
   data green(40,14, 6) /   2.3364856369049069D-02  /
   data green(40,14, 7) /   2.3282302427703322D-02  /
   data green(40,14, 8) /   2.3188129712487945D-02  /
   data green(40,14, 9) /   2.3082773571083223D-02  /
   data green(40,14,10) /   2.2966709745701145D-02  /
   data green(40,14,11) /   2.2840448737075202D-02  /
   data green(40,14,12) /   2.2704530009519851D-02  /
   data green(40,14,13) /   2.2559516154671441D-02  /
   data green(40,14,14) /   2.2405987123427507D-02  /
   data green(40,15, 0) /   2.3409702680224111D-02  /
   data green(40,15, 1) /   2.3403283995297224D-02  /
   data green(40,15, 2) /   2.3384059670214355D-02  /
   data green(40,15, 3) /   2.3352124457293823D-02  /
   data green(40,15, 4) /   2.3307634832973111D-02  /
   data green(40,15, 5) /   2.3250806872521550D-02  /
   data green(40,15, 6) /   2.3181913354465768D-02  /
   data green(40,15, 7) /   2.3101280170303183D-02  /
   data green(40,15, 8) /   2.3009282130916284D-02  /
   data green(40,15, 9) /   2.2906338273141053D-02  /
   data green(40,15,10) /   2.2792906777868764D-02  /
   data green(40,15,11) /   2.2669479614785003D-02  /
   data green(40,15,12) /   2.2536577028514301D-02  /
   data green(40,15,13) /   2.2394741976882442D-02  /
   data green(40,15,14) /   2.2244534624725205D-02  /
   data green(40,15,15) /   2.2086526986767552D-02  /
   data green(40,16, 0) /   2.3213184457050465D-02  /
   data green(40,16, 1) /   2.3206926320708358D-02  /
   data green(40,16, 2) /   2.3188182326687626D-02  /
   data green(40,16, 3) /   2.3157043308430027D-02  /
   data green(40,16, 4) /   2.3113659294039551D-02  /
   data green(40,16, 5) /   2.3058237502678419D-02  /
   data green(40,16, 6) /   2.2991039613501976D-02  /
   data green(40,16, 7) /   2.2912378377239084D-02  /
   data green(40,16, 8) /   2.2822613655307956D-02  /
   data green(40,16, 9) /   2.2722147982669770D-02  /
   data green(40,16,10) /   2.2611421758166297D-02  /
   data green(40,16,11) /   2.2490908169773026D-02  /
   data green(40,16,12) /   2.2361107962141132D-02  /
   data green(40,16,13) /   2.2222544150299048D-02  /
   data green(40,16,14) /   2.2075756776875659D-02  /
   data green(40,16,15) /   2.1921297801237700D-02  /
   data green(40,16,16) /   2.1759726198106021D-02  /
   data green(40,17, 0) /   2.3009344388037856D-02  /
   data green(40,17, 1) /   2.3003249909489952D-02  /
   data green(40,17, 2) /   2.2984995571985928D-02  /
   data green(40,17, 3) /   2.2954668283169027D-02  /
   data green(40,17, 4) /   2.2912411609965128D-02  /
   data green(40,17, 5) /   2.2858423894799806D-02  /
   data green(40,17, 6) /   2.2792955686626034D-02  /
   data green(40,17, 7) /   2.2716306551572668D-02  /
   data green(40,17, 8) /   2.2628821341776420D-02  /
   data green(40,17, 9) /   2.2530886011554705D-02  /
   data green(40,17,10) /   2.2422923077232940D-02  /
   data green(40,17,11) /   2.2305386820565699D-02  /
   data green(40,17,12) /   2.2178758335881579D-02  /
   data green(40,17,13) /   2.2043540518093582D-02  /
   data green(40,17,14) /   2.1900253082940399D-02  /
   data green(40,17,15) /   2.1749427702746408D-02  /
   data green(40,17,16) /   2.1591603331150330D-02  /
   data green(40,17,17) /   2.1427321779219861D-02  /
   data green(40,18, 0) /   2.2798925025245866D-02  /
   data green(40,18, 1) /   2.2792996447823496D-02  /
   data green(40,18, 2) /   2.2775238502957883D-02  /
   data green(40,18, 3) /   2.2745734190310125D-02  /
   data green(40,18, 4) /   2.2704620643468360D-02  /
   data green(40,18, 5) /   2.2652087363323657D-02  /
   data green(40,18, 6) /   2.2588373807689955D-02  /
   data green(40,18, 7) /   2.2513766396888315D-02  /
   data green(40,18, 8) /   2.2428595007768212D-02  /
   data green(40,18, 9) /   2.2333229038530609D-02  /
   data green(40,18,10) /   2.2228073133483638D-02  /
   data green(40,18,11) /   2.2113562660410660D-02  /
   data green(40,18,12) /   2.1990159033637956D-02  /
   data green(40,18,13) /   2.1858344973375122D-02  /
   data green(40,18,14) /   2.1718619786810577D-02  /
   data green(40,18,15) /   2.1571494749208494D-02  /
   data green(40,18,16) /   2.1417488654358833D-02  /
   data green(40,18,17) /   2.1257123593687279D-02  /
   data green(40,18,18) /   2.1090921012632945D-02  /
   data green(40,19, 0) /   2.2582653756362384D-02  /
   data green(40,19, 1) /   2.2576892508303803D-02  /
   data green(40,19, 2) /   2.2559635254363181D-02  /
   data green(40,19, 3) /   2.2530961126229553D-02  /
   data green(40,19, 4) /   2.2491000888523707D-02  /
   data green(40,19, 5) /   2.2439935285995322D-02  /
   data green(40,19, 6) /   2.2377992787317884D-02  /
   data green(40,19, 7) /   2.2305446780343251D-02  /
   data green(40,19, 8) /   2.2222612285469626D-02  /
   data green(40,19, 9) /   2.2129842262983757D-02  /
   data green(40,19,10) /   2.2027523596615344D-02  /
   data green(40,19,11) /   2.1916072838998508D-02  /
   data green(40,19,12) /   2.1795931805326924D-02  /
   data green(40,19,13) /   2.1667563099408371D-02  /
   data green(40,19,14) /   2.1531445651868084D-02  /
   data green(40,19,15) /   2.1388070343803585D-02  /
   data green(40,19,16) /   2.1237935781188468D-02  /
   data green(40,19,17) /   2.1081544276215714D-02  /
   data green(40,19,18) /   2.0919398082013256D-02  /
   data green(40,19,19) /   2.0751995917175789D-02  /
   data green(40,20, 0) /   2.2361238159766868D-02  /
   data green(40,20, 1) /   2.2355644909835659D-02  /
   data green(40,20, 2) /   2.2338890373249747D-02  /
   data green(40,20, 3) /   2.2311049873324048D-02  /
   data green(40,20, 4) /   2.2272247902353059D-02  /
   data green(40,20, 5) /   2.2222656578738911D-02  /
   data green(40,20, 6) /   2.2162493539789180D-02  /
   data green(40,20, 7) /   2.2092019320440578D-02  /
   data green(40,20, 8) /   2.2011534279037036D-02  /
   data green(40,20, 9) /   2.1921375139836028D-02  /
   data green(40,20,10) /   2.1821911227909618D-02  /
   data green(40,20,11) /   2.1713540475456036D-02  /
   data green(40,20,12) /   2.1596685279285197D-02  /
   data green(40,20,13) /   2.1471788287549125D-02  /
   data green(40,20,14) /   2.1339308189917080D-02  /
   data green(40,20,15) /   2.1199715579679906D-02  /
   data green(40,20,16) /   2.1053488949097424D-02  /
   data green(40,20,17) /   2.0901110871079123D-02  /
   data green(40,20,18) /   2.0743064411423845D-02  /
   data green(40,20,19) /   2.0579829806714682D-02  /
   data green(40,20,20) /   2.0411881433920906D-02  /
   data green(40,21, 0) /   2.2135362061533142D-02  /
   data green(40,21, 1) /   2.2129936777958055D-02  /
   data green(40,21, 2) /   2.2113684889243861D-02  /
   data green(40,21, 3) /   2.2086677986915007D-02  /
   data green(40,21, 4) /   2.2049034415594856D-02  /
   data green(40,21, 5) /   2.2000917835736725D-02  /
   data green(40,21, 6) /   2.1942535259653749D-02  /
   data green(40,21, 7) /   2.1874134606762786D-02  /
   data green(40,21, 8) /   2.1796001833955459D-02  /
   data green(40,21, 9) /   2.1708457704920161D-02  /
   data green(40,21,10) /   2.1611854267852875D-02  /
   data green(40,21,11) /   2.1506571114224507D-02  /
   data green(40,21,12) /   2.1393011492146622D-02  /
   data green(40,21,13) /   2.1271598346532423D-02  /
   data green(40,21,14) /   2.1142770354910363D-02  /
   data green(40,21,15) /   2.1006978022709268D-02  /
   data green(40,21,16) /   2.0864679895436295D-02  /
   data green(40,21,17) /   2.0716338937775285D-02  /
   data green(40,21,18) /   2.0562419121607670D-02  /
   data green(40,21,19) /   2.0403382256646056D-02  /
   data green(40,21,20) /   2.0239685089086291D-02  /
   data green(40,21,21) /   2.0071776685695774D-02  /
   data green(40,22, 0) /   2.1905682280404650D-02  /
   data green(40,22, 1) /   2.1900424291756589D-02  /
   data green(40,22, 2) /   2.1884673067294001D-02  /
   data green(40,22, 3) /   2.1858496557777319D-02  /
   data green(40,22, 4) /   2.1822007108663336D-02  /
   data green(40,22, 5) /   2.1775360123782641D-02  /
   data green(40,22, 6) /   2.1718752238376906D-02  /
   data green(40,22, 7) /   2.1652419043342312D-02  /
   data green(40,22, 8) /   2.1576632411696838D-02  /
   data green(40,22, 9) /   2.1491697485593514D-02  /
   data green(40,22,10) /   2.1397949387445762D-02  /
   data green(40,22,11) /   2.1295749721832215D-02  /
   data green(40,22,12) /   2.1185482935821457D-02  /
   data green(40,22,13) /   2.1067552604318949D-02  /
   data green(40,22,14) /   2.0942377704180212D-02  /
   data green(40,22,15) /   2.0810388936414355D-02  /
   data green(40,22,16) /   2.0672025150119976D-02  /
   data green(40,22,17) /   2.0527729915171867D-02  /
   data green(40,22,18) /   2.0377948283438717D-02  /
   data green(40,22,19) /   2.0223123770769291D-02  /
   data green(40,22,20) /   2.0063695584422749D-02  /
   data green(40,22,21) /   1.9900096113284502D-02  /
   data green(40,22,22) /   1.9732748691307329D-02  /
   data green(40,23, 0) /   2.1672826032748289D-02  /
   data green(40,23, 1) /   2.1667734089519178D-02  /
   data green(40,23, 2) /   2.1652479815467546D-02  /
   data green(40,23, 3) /   2.1627127623619159D-02  /
   data green(40,23, 4) /   2.1591784028619831D-02  /
   data green(40,23, 5) /   2.1546596406472284D-02  /
   data green(40,23, 6) /   2.1491751298035149D-02  /
   data green(40,23, 7) /   2.1427472294335767D-02  /
   data green(40,23, 8) /   2.1354017550142441D-02  /
   data green(40,23, 9) /   2.1271676978969405D-02  /
   data green(40,23,10) /   2.1180769187574510D-02  /
   data green(40,23,11) /   2.1081638210972548D-02  /
   data green(40,23,12) /   2.0974650110035842D-02  /
   data green(40,23,13) /   2.0860189492982812D-02  /
   data green(40,23,14) /   2.0738656019629074D-02  /
   data green(40,23,15) /   2.0610460943418580D-02  /
   data green(40,23,16) /   2.0476023741225802D-02  /
   data green(40,23,17) /   2.0335768875009524D-02  /
   data green(40,23,18) /   2.0190122722891949D-02  /
   data green(40,23,19) /   2.0039510710415920D-02  /
   data green(40,23,20) /   1.9884354665852252D-02  /
   data green(40,23,21) /   1.9725070416717325D-02  /
   data green(40,23,22) /   1.9562065638304864D-02  /
   data green(40,23,23) /   1.9395737959184912D-02  /
   data green(40,24, 0) /   2.1437388958485154D-02  /
   data green(40,24, 1) /   2.1432461294294443D-02  /
   data green(40,24, 2) /   2.1417698709438736D-02  /
   data green(40,24, 3) /   2.1393162191939874D-02  /
   data green(40,24, 4) /   2.1358952610063868D-02  /
   data green(40,24, 5) /   2.1315209563088269D-02  /
   data green(40,24, 6) /   2.1262109808515599D-02  /
   data green(40,24, 7) /   2.1199865300261730D-02  /
   data green(40,24, 8) /   2.1128720880014100D-02  /
   data green(40,24, 9) /   2.1048951670137872D-02  /
   data green(40,24,10) /   2.0960860221048927D-02  /
   data green(40,24,11) /   2.0864773468792917D-02  /
   data green(40,24,12) /   2.0761039559670756D-02  /
   data green(40,24,13) /   2.0650024598211206D-02  /
   data green(40,24,14) /   2.0532109372750046D-02  /
   data green(40,24,15) /   2.0407686109526431D-02  /
   data green(40,24,16) /   2.0277155301779363D-02  /
   data green(40,24,17) /   2.0140922655070833D-02  /
   data green(40,24,18) /   1.9999396184233808D-02  /
   data green(40,24,19) /   1.9852983491193316D-02  /
   data green(40,24,20) /   1.9702089246668133D-02  /
   data green(40,24,21) /   1.9547112892636505D-02  /
   data green(40,24,22) /   1.9388446576613920D-02  /
   data green(40,24,23) /   1.9226473323384154D-02  /
   data green(40,24,24) /   1.9061565444949316D-02  /
   data green(40,25, 0) /   2.1199933720879912D-02  /
   data green(40,25, 1) /   2.1195168112399503D-02  /
   data green(40,25, 2) /   2.1180890587207065D-02  /
   data green(40,25, 3) /   2.1157158828614589D-02  /
   data green(40,25, 4) /   2.1124068255504575D-02  /
   data green(40,25, 5) /   2.1081750959035105D-02  /
   data green(40,25, 6) /   2.1030374246723622D-02  /
   data green(40,25, 7) /   2.0970138825182508D-02  /
   data green(40,25, 8) /   2.0901276659765728D-02  /
   data green(40,25, 9) /   2.0824048555057342D-02  /
   data green(40,25,10) /   2.0738741504341924D-02  /
   data green(40,25,11) /   2.0645665858870198D-02  /
   data green(40,25,12) /   2.0545152368868103D-02  /
   data green(40,25,13) /   2.0437549147894349D-02  /
   data green(40,25,14) /   2.0323218610451374D-02  /
   data green(40,25,15) /   2.0202534429861998D-02  /
   data green(40,25,16) /   2.0075878559539240D-02  /
   data green(40,25,17) /   1.9943638356118808D-02  /
   data green(40,25,18) /   1.9806203837719406D-02  /
   data green(40,25,19) /   1.9663965105067651D-02  /
   data green(40,25,20) /   1.9517309947582339D-02  /
   data green(40,25,21) /   1.9366621650939733D-02  /
   data green(40,25,22) /   1.9212277017301362D-02  /
   data green(40,25,23) /   1.9054644604403984D-02  /
   data green(40,25,24) /   1.8894083185187750D-02  /
   data green(40,25,25) /   1.8730940425639032D-02  /
   data green(40,26, 0) /   2.0960989127615990D-02  /
   data green(40,26, 1) /   2.0956382952467070D-02  /
   data green(40,26, 2) /   2.0942582662121012D-02  /
   data green(40,26, 3) /   2.0919642761079476D-02  /
   data green(40,26, 4) /   2.0887653425185754D-02  /
   data green(40,26, 5) /   2.0846739519177274D-02  /
   data green(40,26, 6) /   2.0797059250790405D-02  /
   data green(40,26, 7) /   2.0738802489691819D-02  /
   data green(40,26, 8) /   2.0672188785867088D-02  /
   data green(40,26, 9) /   2.0597465127288297D-02  /
   data green(40,26,10) /   2.0514903480575618D-02  /
   data green(40,26,11) /   2.0424798160892623D-02  /
   data green(40,26,12) /   2.0327463078466064D-02  /
   data green(40,26,13) /   2.0223228908943534D-02  /
   data green(40,26,14) /   2.0112440233402095D-02  /
   data green(40,26,15) /   1.9995452691335347D-02  /
   data green(40,26,16) /   1.9872630186551832D-02  /
   data green(40,26,17) /   1.9744342181802908D-02  /
   data green(40,26,18) /   1.9610961113326289D-02  /
   data green(40,26,19) /   1.9472859951536543D-02  /
   data green(40,26,20) /   1.9330409929005903D-02  /
   data green(40,26,21) /   1.9183978451823916D-02  /
   data green(40,26,22) /   1.9033927205549108D-02  /
   data green(40,26,23) /   1.8880610462390094D-02  /
   data green(40,26,24) /   1.8724373592068654D-02  /
   data green(40,26,25) /   1.8565551775089065D-02  /
   data green(40,26,26) /   1.8404468913906735D-02  /
   data green(40,27, 0) /   2.0721049717492598D-02  /
   data green(40,27, 1) /   2.0716600009523047D-02  /
   data green(40,27, 2) /   2.0703268099161072D-02  /
   data green(40,27, 3) /   2.0681105441835077D-02  /
   data green(40,27, 4) /   2.0650197183133111D-02  /
   data green(40,27, 5) /   2.0610661252172330D-02  /
   data green(40,27, 6) /   2.0562647118946526D-02  /
   data green(40,27, 7) /   2.0506334241174859D-02  /
   data green(40,27, 8) /   2.0441930231947660D-02  /
   data green(40,27, 9) /   2.0369668784208872D-02  /
   data green(40,27,10) /   2.0289807391706938D-02  /
   data green(40,27,11) /   2.0202624908423450D-02  /
   data green(40,27,12) /   2.0108418989639332D-02  /
   data green(40,27,13) /   2.0007503457760915D-02  /
   data green(40,27,14) /   1.9900205634888118D-02  /
   data green(40,27,15) /   1.9786863681983301D-02  /
   data green(40,27,16) /   1.9667823981544454D-02  /
   data green(40,27,17) /   1.9543438597063822D-02  /
   data green(40,27,18) /   1.9414062838441961D-02  /
   data green(40,27,19) /   1.9280052958099505D-02  /
   data green(40,27,20) /   1.9141763997952952D-02  /
   data green(40,27,21) /   1.8999547802847832D-02  /
   data green(40,27,22) /   1.8853751211604784D-02  /
   data green(40,27,23) /   1.8704714432641888D-02  /
   data green(40,27,24) /   1.8552769607277215D-02  /
   data green(40,27,25) /   1.8398239560350090D-02  /
   data green(40,27,26) /   1.8241436734771403D-02  /
   data green(40,27,27) /   1.8082662304039227D-02  /
   data green(40,28, 0) /   2.0480575756013308D-02  /
   data green(40,28, 1) /   2.0476279257508077D-02  /
   data green(40,28, 2) /   2.0463405998325437D-02  /
   data green(40,28, 3) /   2.0442004516819851D-02  /
   data green(40,28, 4) /   2.0412155144954185D-02  /
   data green(40,28, 5) /   2.0373969172561755D-02  /
   data green(40,28, 6) /   2.0327587701294066D-02  /
   data green(40,28, 7) /   2.0273180211260025D-02  /
   data green(40,28, 8) /   2.0210942868599407D-02  /
   data green(40,28, 9) /   2.0141096606558852D-02  /
   data green(40,28,10) /   2.0063885015947937D-02  /
   data green(40,28,11) /   1.9979572083083506D-02  /
   data green(40,28,12) /   1.9888439814468813D-02  /
   data green(40,28,13) /   1.9790785787530896D-02  /
   data green(40,28,14) /   1.9686920665825010D-02  /
   data green(40,28,15) /   1.9577165715311103D-02  /
   data green(40,28,16) /   1.9461850355744870D-02  /
   data green(40,28,17) /   1.9341309778046926D-02  /
   data green(40,28,18) /   1.9215882654875258D-02  /
   data green(40,28,19) /   1.9085908967679595D-02  /
   data green(40,28,20) /   1.8951727969412451D-02  /
   data green(40,28,21) /   1.8813676297944466D-02  /
   data green(40,28,22) /   1.8672086251202712D-02  /
   data green(40,28,23) /   1.8527284231220721D-02  /
   data green(40,28,24) /   1.8379589360738670D-02  /
   data green(40,28,25) /   1.8229312272781607D-02  /
   data green(40,28,26) /   1.8076754070814095D-02  /
   data green(40,28,27) /   1.7922205454641224D-02  /
   data green(40,28,28) /   1.7765946005207655D-02  /
   data green(40,29, 0) /   2.0239993583746296D-02  /
   data green(40,29, 1) /   2.0235846794253439D-02  /
   data green(40,29, 2) /   2.0223421729516928D-02  /
   data green(40,29, 3) /   2.0202764143694429D-02  /
   data green(40,29, 4) /   2.0173949773317622D-02  /
   data green(40,29, 5) /   2.0137083567666882D-02  /
   data green(40,29, 6) /   2.0092298632863041D-02  /
   data green(40,29, 7) /   2.0039754910387661D-02  /
   data green(40,29, 8) /   1.9979637615485157D-02  /
   data green(40,29, 9) /   1.9912155464840531D-02  /
   data green(40,29,10) /   1.9837538725967992D-02  /
   data green(40,29,11) /   1.9756037122832498D-02  /
   data green(40,29,12) /   1.9667917633342486D-02  /
   data green(40,29,13) /   1.9573462214520981D-02  /
   data green(40,29,14) /   1.9472965490441682D-02  /
   data green(40,29,15) /   1.9366732436493870D-02  /
   data green(40,29,16) /   1.9255076091326005D-02  /
   data green(40,29,17) /   1.9138315325037528D-02  /
   data green(40,29,18) /   1.9016772688976301D-02  /
   data green(40,29,19) /   1.8890772368991490D-02  /
   data green(40,29,20) /   1.8760638260319290D-02  /
   data green(40,29,21) /   1.8626692178564160D-02  /
   data green(40,29,22) /   1.8489252217589233D-02  /
   data green(40,29,23) /   1.8348631261639221D-02  /
   data green(40,29,24) /   1.8205135655762514D-02  /
   data green(40,29,25) /   1.8059064035632733D-02  /
   data green(40,29,26) /   1.7910706315234337D-02  /
   data green(40,29,27) /   1.7760342828594889D-02  /
   data green(40,29,28) /   1.7608243619826665D-02  /
   data green(40,29,29) /   1.7454667874180662D-02  /
   data green(40,30, 0) /   1.9999696263276197D-02  /
   data green(40,30, 1) /   1.9995695484846098D-02  /
   data green(40,30, 2) /   1.9983707565209896D-02  /
   data green(40,30, 3) /   1.9963775606883231D-02  /
   data green(40,30, 4) /   1.9935970968744571D-02  /
   data green(40,30, 5) /   1.9900392557915435D-02  /
   data green(40,30, 6) /   1.9857165857768182D-02  /
   data green(40,30, 7) /   1.9806441710684929D-02  /
   data green(40,30, 8) /   1.9748394878481658D-02  /
   data green(40,30, 9) /   1.9683222406994390D-02  /
   data green(40,30,10) /   1.9611141824115708D-02  /
   data green(40,30,11) /   1.9532389202516636D-02  /
   data green(40,30,12) /   1.9447217119374225D-02  /
   data green(40,30,13) /   1.9355892545667181D-02  /
   data green(40,30,14) /   1.9258694697046785D-02  /
   data green(40,30,15) /   1.9155912877013345D-02  /
   data green(40,30,16) /   1.9047844341223190D-02  /
   data green(40,30,17) /   1.8934792209326069D-02  /
   data green(40,30,18) /   1.8817063447905961D-02  /
   data green(40,30,19) /   1.8694966944986791D-02  /
   data green(40,30,20) /   1.8568811693287529D-02  /
   data green(40,30,21) /   1.8438905096074331D-02  /
   data green(40,30,22) /   1.8305551406160959D-02  /
   data green(40,30,23) /   1.8169050305435423D-02  /
   data green(40,30,24) /   1.8029695629310945D-02  /
   data green(40,30,25) /   1.7887774237766330D-02  /
   data green(40,30,26) /   1.7743565032192977D-02  /
   data green(40,30,27) /   1.7597338115127936D-02  /
   data green(40,30,28) /   1.7449354088135100D-02  /
   data green(40,30,29) /   1.7299863481599908D-02  /
   data green(40,30,30) /   1.7149106309017108D-02  /
   data green(40,31, 0) /   1.9760044473517993D-02  /
   data green(40,31, 1) /   1.9756185852252629D-02  /
   data green(40,31, 2) /   1.9744623560062805D-02  /
   data green(40,31, 3) /   1.9725398179030074D-02  /
   data green(40,31, 4) /   1.9698576906048183D-02  /
   data green(40,31, 5) /   1.9664252901797789D-02  /
   data green(40,31, 6) /   1.9622544396699783D-02  /
   data green(40,31, 7) /   1.9573593570577641D-02  /
   data green(40,31, 8) /   1.9517565226635163D-02  /
   data green(40,31, 9) /   1.9454645283610485D-02  /
   data green(40,31,10) /   1.9385039112524777D-02  /
   data green(40,31,11) /   1.9308969746255098D-02  /
   data green(40,31,12) /   1.9226675991209433D-02  /
   data green(40,31,13) /   1.9138410470679334D-02  /
   data green(40,31,14) /   1.9044437629031841D-02  /
   data green(40,31,15) /   1.8945031724838239D-02  /
   data green(40,31,16) /   1.8840474839404313D-02  /
   data green(40,31,17) /   1.8731054925058332D-02  /
   data green(40,31,18) /   1.8617063915070304D-02  /
   data green(40,31,19) /   1.8498795914324063D-02  /
   data green(40,31,20) /   1.8376545486944364D-02  /
   data green(40,31,21) /   1.8250606054091938D-02  /
   data green(40,31,22) /   1.8121268412167352D-02  /
   data green(40,31,23) /   1.7988819378786731D-02  /
   data green(40,31,24) /   1.7853540571172097D-02  /
   data green(40,31,25) /   1.7715707319087729D-02  /
   data green(40,31,26) /   1.7575587712187243D-02  /
   data green(40,31,27) /   1.7433441779639562D-02  /
   data green(40,31,28) /   1.7289520798188050D-02  /
   data green(40,31,29) /   1.7144066723368215D-02  /
   data green(40,31,30) /   1.6997311737459466D-02  /
   data green(40,31,31) /   1.6849477906862566D-02  /
   data green(40,32, 0) /   1.9521367603830067D-02  /
   data green(40,32, 1) /   1.9517647167722574D-02  /
   data green(40,32, 2) /   1.9506498630247837D-02  /
   data green(40,32, 3) /   1.9487960182053211D-02  /
   data green(40,32, 4) /   1.9462095070182209D-02  /
   data green(40,32, 5) /   1.9428990999943364D-02  /
   data green(40,32, 6) /   1.9388759313118901D-02  /
   data green(40,32, 7) /   1.9341533957529900D-02  /
   data green(40,32, 8) /   1.9287470266472359D-02  /
   data green(40,32, 9) /   1.9226743569491677D-02  /
   data green(40,32,10) /   1.9159547658302929D-02  /
   data green(40,32,11) /   1.9086093133344400D-02  /
   data green(40,32,12) /   1.9006605657459628D-02  /
   data green(40,32,13) /   1.8921324143541193D-02  /
   data green(40,32,14) /   1.8830498902674537D-02  /
   data green(40,32,15) /   1.8734389778440320D-02  /
   data green(40,32,16) /   1.8633264291640521D-02  /
   data green(40,32,17) /   1.8527395817884800D-02  /
   data green(40,32,18) /   1.8417061818299843D-02  /
   data green(40,32,19) /   1.8302542141194620D-02  /
   data green(40,32,20) /   1.8184117409920486D-02  /
   data green(40,32,21) /   1.8062067509492054D-02  /
   data green(40,32,22) /   1.7936670181861129D-02  /
   data green(40,32,23) /   1.7808199737131887D-02  /
   data green(40,32,24) /   1.7676925885528091D-02  /
   data green(40,32,25) /   1.7543112692620776D-02  /
   data green(40,32,26) /   1.7407017658231869D-02  /
   data green(40,32,27) /   1.7268890917570317D-02  /
   data green(40,32,28) /   1.7128974561546679D-02  /
   data green(40,32,29) /   1.6987502071854007D-02  /
   data green(40,32,30) /   1.6844697865293832D-02  /
   data green(40,32,31) /   1.6700776940957857D-02  /
   data green(40,32,32) /   1.6555944623231629D-02  /
   data green(40,33, 0) /   1.9283965004494096D-02  /
   data green(40,33, 1) /   1.9280378697607826D-02  /
   data green(40,33, 2) /   1.9269631789339273D-02  /
   data green(40,33, 3) /   1.9251760204983186D-02  /
   data green(40,33, 4) /   1.9226823449156025D-02  /
   data green(40,33, 5) /   1.9194904056579372D-02  /
   data green(40,33, 6) /   1.9156106837148249D-02  /
   data green(40,33, 7) /   1.9110557928751864D-02  /
   data green(40,33, 8) /   1.9058403674469478D-02  /
   data green(40,33, 9) /   1.8999809343440002D-02  /
   data green(40,33,10) /   1.8934957716841200D-02  /
   data green(40,33,11) /   1.8864047561971676D-02  /
   data green(40,33,12) /   1.8787292018389173D-02  /
   data green(40,33,13) /   1.8704916920427736D-02  /
   data green(40,33,14) /   1.8617159080218843D-02  /
   data green(40,33,15) /   1.8524264554621736D-02  /
   data green(40,33,16) /   1.8426486918283067D-02  /
   data green(40,33,17) /   1.8324085563466111D-02  /
   data green(40,33,18) /   1.8217324045389553D-02  /
   data green(40,33,19) /   1.8106468489676328D-02  /
   data green(40,33,20) /   1.7991786076211960D-02  /
   data green(40,33,21) /   1.7873543611327061D-02  /
   data green(40,33,22) /   1.7752006197818120D-02  /
   data green(40,33,23) /   1.7627436009968953D-02  /
   data green(40,33,24) /   1.7500091178484653D-02  /
   data green(40,33,25) /   1.7370224788143252D-02  /
   data green(40,33,26) /   1.7238083989043567D-02  /
   data green(40,33,27) /   1.7103909220601343D-02  /
   data green(40,33,28) /   1.6967933545937640D-02  /
   data green(40,33,29) /   1.6830382093017800D-02  /
   data green(40,33,30) /   1.6691471597835232D-02  /
   data green(40,33,31) /   1.6551410044087837D-02  /
   data green(40,33,32) /   1.6410396393149250D-02  /
   data green(40,33,33) /   1.6268620397683731D-02  /
   data green(40,34, 0) /   1.9048107354508199D-02  /
   data green(40,34, 1) /   1.9044651067599269D-02  /
   data green(40,34, 2) /   1.9034293501926678D-02  /
   data green(40,34, 3) /   1.9017068440022569D-02  /
   data green(40,34, 4) /   1.8993031845835940D-02  /
   data green(40,34, 5) /   1.8962261360677437D-02  /
   data green(40,34, 6) /   1.8924855610053130D-02  /
   data green(40,34, 7) /   1.8880933333458590D-02  /
   data green(40,34, 8) /   1.8830632352046840D-02  /
   data green(40,34, 9) /   1.8774108391507491D-02  /
   data green(40,34,10) /   1.8711533779443918D-02  /
   data green(40,34,11) /   1.8643096037975072D-02  /
   data green(40,34,12) /   1.8568996393200816D-02  /
   data green(40,34,13) /   1.8489448223557726D-02  /
   data green(40,34,14) /   1.8404675468976721D-02  /
   data green(40,34,15) /   1.8314911022169769D-02  /
   data green(40,34,16) /   1.8220395122370624D-02  /
   data green(40,34,17) /   1.8121373770492312D-02  /
   data green(40,34,18) /   1.8018097183008081D-02  /
   data green(40,34,19) /   1.7910818299981427D-02  /
   data green(40,34,20) /   1.7799791360634760D-02  /
   data green(40,34,21) /   1.7685270557721428D-02  /
   data green(40,34,22) /   1.7567508779815694D-02  /
   data green(40,34,23) /   1.7446756448515312D-02  /
   data green(40,34,24) /   1.7323260455510706D-02  /
   data green(40,34,25) /   1.7197263202552181D-02  /
   data green(40,34,26) /   1.7069001745576840D-02  /
   data green(40,34,27) /   1.6938707042658124D-02  /
   data green(40,34,28) /   1.6806603304033316D-02  /
   data green(40,34,29) /   1.6672907441251331D-02  /
   data green(40,34,30) /   1.6537828611468276D-02  /
   data green(40,34,31) /   1.6401567852096802D-02  /
   data green(40,34,32) /   1.6264317800376484D-02  /
   data green(40,34,33) /   1.6126262491967771D-02  /
   data green(40,34,34) /   1.5987577232362420D-02  /
   data green(40,35, 0) /   1.8814038112090536D-02  /
   data green(40,35, 1) /   1.8810707709820292D-02  /
   data green(40,35, 2) /   1.8800727120519341D-02  /
   data green(40,35, 3) /   1.8784128102605414D-02  /
   data green(40,35, 4) /   1.8760963274706353D-02  /
   data green(40,35, 5) /   1.8731305653240313D-02  /
   data green(40,35, 6) /   1.8695248016225602D-02  /
   data green(40,35, 7) /   1.8652902104134225D-02  /
   data green(40,35, 8) /   1.8604397671161829D-02  /
   data green(40,35, 9) /   1.8549881402480352D-02  /
   data green(40,35,10) /   1.8489515714816038D-02  /
   data green(40,35,11) /   1.8423477459023847D-02  /
   data green(40,35,12) /   1.8351956544192399D-02  /
   data green(40,35,13) /   1.8275154503212310D-02  /
   data green(40,35,14) /   1.8193283019691659D-02  /
   data green(40,35,15) /   1.8106562435634747D-02  /
   data green(40,35,16) /   1.8015220258455752D-02  /
   data green(40,35,17) /   1.7919489684728767D-02  /
   data green(40,35,18) /   1.7819608156635015D-02  /
   data green(40,35,19) /   1.7715815965418099D-02  /
   data green(40,35,20) /   1.7608354914359272D-02  /
   data green(40,35,21) /   1.7497467051895577D-02  /
   data green(40,35,22) /   1.7383393483580999D-02  /
   data green(40,35,23) /   1.7266373269683180D-02  /
   data green(40,35,24) /   1.7146642413362011D-02  /
   data green(40,35,25) /   1.7024432942625681D-02  /
   data green(40,35,26) /   1.6899972087636862D-02  /
   data green(40,35,27) /   1.6773481553466415D-02  /
   data green(40,35,28) /   1.6645176887081260D-02  /
   data green(40,35,29) /   1.6515266936212948D-02  /
   data green(40,35,30) /   1.6383953396789840D-02  /
   data green(40,35,31) /   1.6251430444823113D-02  /
   data green(40,35,32) /   1.6117884448010277D-02  /
   data green(40,35,33) /   1.5983493751850087D-02  /
   data green(40,35,34) /   1.5848428534735472D-02  /
   data green(40,35,35) /   1.5712850726294726D-02  /
   data green(40,36, 0) /   1.8581975017674678D-02  /
   data green(40,36, 1) /   1.8578766362589778D-02  /
   data green(40,36, 2) /   1.8569150375648796D-02  /
   data green(40,36, 3) /   1.8553156905520285D-02  /
   data green(40,36, 4) /   1.8530835413873743D-02  /
   data green(40,36, 5) /   1.8502254551294658D-02  /
   data green(40,36, 6) /   1.8467501573583643D-02  /
   data green(40,36, 7) /   1.8426681608122612D-02  /
   data green(40,36, 8) /   1.8379916782294271D-02  /
   data green(40,36, 9) /   1.8327345227923855D-02  /
   data green(40,36,10) /   1.8269119977329752D-02  /
   data green(40,36,11) /   1.8205407767792871D-02  /
   data green(40,36,12) /   1.8136387772068090D-02  /
   data green(40,36,13) /   1.8062250272963470D-02  /
   data green(40,36,14) /   1.7983195300017865D-02  /
   data green(40,36,15) /   1.7899431245938083D-02  /
   data green(40,36,16) /   1.7811173479749656D-02  /
   data green(40,36,17) /   1.7718642972612438D-02  /
   data green(40,36,18) /   1.7622064951001982D-02  /
   data green(40,36,19) /   1.7521667590513549D-02  /
   data green(40,36,20) /   1.7417680761958366D-02  /
   data green(40,36,21) /   1.7310334839746217D-02  /
   data green(40,36,22) /   1.7199859580830607D-02  /
   data green(40,36,23) /   1.7086483080780653D-02  /
   data green(40,36,24) /   1.6970430811875756D-02  /
   data green(40,36,25) /   1.6851924746528212D-02  /
   data green(40,36,26) /   1.6731182567854434D-02  /
   data green(40,36,27) /   1.6608416967856605D-02  /
   data green(40,36,28) /   1.6483835032460044D-02  /
   data green(40,36,29) /   1.6357637711583679D-02  /
   data green(40,36,30) /   1.6230019371508749D-02  /
   data green(40,36,31) /   1.6101167426051510D-02  /
   data green(40,36,32) /   1.5971262042434309D-02  /
   data green(40,36,33) /   1.5840475917280380D-02  /
   data green(40,36,34) /   1.5708974117819099D-02  /
   data green(40,36,35) /   1.5576913983169483D-02  /
   data green(40,36,36) /   1.5444445080457373D-02  /
   data green(40,37, 0) /   1.8352111623389820D-02  /
   data green(40,37, 1) /   1.8349020596869371D-02  /
   data green(40,37, 2) /   1.8339756893249171D-02  /
   data green(40,37, 3) /   1.8324348561285501D-02  /
   data green(40,37, 4) /   1.8302842086656501D-02  /
   data green(40,37, 5) /   1.8275302003132543D-02  /
   data green(40,37, 6) /   1.8241810357174513D-02  /
   data green(40,37, 7) /   1.8202466034627416D-02  /
   data green(40,37, 8) /   1.8157383960249004D-02  /
   data green(40,37, 9) /   1.8106694182604074D-02  /
   data green(40,37,10) /   1.8050540858326584D-02  /
   data green(40,37,11) /   1.7989081150877152D-02  /
   data green(40,37,12) /   1.7922484059686931D-02  /
   data green(40,37,13) /   1.7850929195979467D-02  /
   data green(40,37,14) /   1.7774605521609931D-02  /
   data green(40,37,15) /   1.7693710066975000D-02  /
   data green(40,37,16) /   1.7608446643457461D-02  /
   data green(40,37,17) /   1.7519024565012729D-02  /
   data green(40,37,18) /   1.7425657392422658D-02  /
   data green(40,37,19) /   1.7328561712479339D-02  /
   data green(40,37,20) /   1.7227955962964697D-02  /
   data green(40,37,21) /   1.7124059312806632D-02  /
   data green(40,37,22) /   1.7017090605261524D-02  /
   data green(40,37,23) /   1.6907267370438014D-02  /
   data green(40,37,24) /   1.6794804911972605D-02  /
   data green(40,37,25) /   1.6679915471224948D-02  /
   data green(40,37,26) /   1.6562807471005324D-02  /
   data green(40,37,27) /   1.6443684839598741D-02  /
   data green(40,37,28) /   1.6322746414722672D-02  /
   data green(40,37,29) /   1.6200185426059673D-02  /
   data green(40,37,30) /   1.6076189054145321D-02  /
   data green(40,37,31) /   1.5950938062667147D-02  /
   data green(40,37,32) /   1.5824606500638674D-02  /
   data green(40,37,33) /   1.5697361470448569D-02  /
   data green(40,37,34) /   1.5569362957439452D-02  /
   data green(40,37,35) /   1.5440763716435865D-02  /
   data green(40,37,36) /   1.5311709210503585D-02  /
   data green(40,37,37) /   1.5182337597173757D-02  /
   data green(40,38, 0) /   1.8124618826982446D-02  /
   data green(40,38, 1) /   1.8121641347364954D-02  /
   data green(40,38, 2) /   1.8112717717326209D-02  /
   data green(40,38, 3) /   1.8097874290856356D-02  /
   data green(40,38, 4) /   1.8077154750938143D-02  /
   data green(40,38, 5) /   1.8050619753107244D-02  /
   data green(40,38, 6) /   1.8018346434448079D-02  /
   data green(40,38, 7) /   1.7980427795781809D-02  /
   data green(40,38, 8) /   1.7936971966666316D-02  /
   data green(40,38, 9) /   1.7888101364445134D-02  /
   data green(40,38,10) /   1.7833951759920592D-02  /
   data green(40,38,11) /   1.7774671263258504D-02  /
   data green(40,38,12) /   1.7710419244448047D-02  /
   data green(40,38,13) /   1.7641365203033537D-02  /
   data green(40,38,14) /   1.7567687601916588D-02  /
   data green(40,38,15) /   1.7489572679811063D-02  /
   data green(40,38,16) /   1.7407213256444766D-02  /
   data green(40,38,17) /   1.7320807543873006D-02  /
   data green(40,38,18) /   1.7230557976333907D-02  /
   data green(40,38,19) /   1.7136670069974609D-02  /
   data green(40,38,20) /   1.7039351322548947D-02  /
   data green(40,38,21) /   1.6938810161873975D-02  /
   data green(40,38,22) /   1.6835254950469763D-02  /
   data green(40,38,23) /   1.6728893052433379D-02  /
   data green(40,38,24) /   1.6619929967243409D-02  /
   data green(40,38,25) /   1.6508568533885205D-02  /
   data green(40,38,26) /   1.6395008207452638D-02  /
   data green(40,38,27) /   1.6279444409237197D-02  /
   data green(40,38,28) /   1.6162067950274413D-02  /
   data green(40,38,29) /   1.6043064527390651D-02  /
   data green(40,38,30) /   1.5922614289984070D-02  /
   data green(40,38,31) /   1.5800891475085680D-02  /
   data green(40,38,32) /   1.5678064107676269D-02  /
   data green(40,38,33) /   1.5554293762779971D-02  /
   data green(40,38,34) /   1.5429735385507850D-02  /
   data green(40,38,35) /   1.5304537164977886D-02  /
   data green(40,38,36) /   1.5178840457881303D-02  /
   data green(40,38,37) /   1.5052779757390864D-02  /
   data green(40,38,38) /   1.4926482703102097D-02  /
   data green(40,39, 0) /   1.7899646391800302D-02  /
   data green(40,39, 1) /   1.7896778429910203D-02  /
   data green(40,39, 2) /   1.7888182819563093D-02  /
   data green(40,39, 3) /   1.7873884320347280D-02  /
   data green(40,39, 4) /   1.7853923978017407D-02  /
   data green(40,39, 5) /   1.7828358797785548D-02  /
   data green(40,39, 6) /   1.7797261294091152D-02  /
   data green(40,39, 7) /   1.7760718923790123D-02  /
   data green(40,39, 8) /   1.7718833411378067D-02  /
   data green(40,39, 9) /   1.7671719976322792D-02  /
   data green(40,39,10) /   1.7619506473796075D-02  /
   data green(40,39,11) /   1.7562332461042027D-02  /
   data green(40,39,12) /   1.7500348202287325D-02  /
   data green(40,39,13) /   1.7433713625482735D-02  /
   data green(40,39,14) /   1.7362597244272058D-02  /
   data green(40,39,15) /   1.7287175058427126D-02  /
   data green(40,39,16) /   1.7207629445585782D-02  /
   data green(40,39,17) /   1.7124148056511374D-02  /
   data green(40,39,18) /   1.7036922725286126D-02  /
   data green(40,39,19) /   1.6946148404892281D-02  /
   data green(40,39,20) /   1.6852022137556967D-02  /
   data green(40,39,21) /   1.6754742068076030D-02  /
   data green(40,39,22) /   1.6654506507121677D-02  /
   data green(40,39,23) /   1.6551513050310447D-02  /
   data green(40,39,24) /   1.6445957757590942D-02  /
   data green(40,39,25) /   1.6338034396330019D-02  /
   data green(40,39,26) /   1.6227933750353789D-02  /
   data green(40,39,27) /   1.6115842996151424D-02  /
   data green(40,39,28) /   1.6001945146491926D-02  /
   data green(40,39,29) /   1.5886418560842105D-02  /
   data green(40,39,30) /   1.5769436521217099D-02  /
   data green(40,39,31) /   1.5651166871443495D-02  /
   data green(40,39,32) /   1.5531771717269003D-02  /
   data green(40,39,33) /   1.5411407184309637D-02  /
   data green(40,39,34) /   1.5290223230479617D-02  /
   data green(40,39,35) /   1.5168363509294980D-02  /
   data green(40,39,36) /   1.5045965280270487D-02  /
   data green(40,39,37) /   1.4923159362534461D-02  /
   data green(40,39,38) /   1.4800070127756476D-02  /
   data green(40,39,39) /   1.4676815528511371D-02  /
   data green(40,40, 0) /   1.7677324437795450D-02  /
   data green(40,40, 1) /   1.7674562030089814D-02  /
   data green(40,40, 2) /   1.7666282580824578D-02  /
   data green(40,40, 3) /   1.7652509350737811D-02  /
   data green(40,40, 4) /   1.7633280905936917D-02  /
   data green(40,40, 5) /   1.7608650818458253D-02  /
   data green(40,40, 6) /   1.7578687253453740D-02  /
   data green(40,40, 7) /   1.7543472449213587D-02  /
   data green(40,40, 8) /   1.7503102097739441D-02  /
   data green(40,40, 9) /   1.7457684634899797D-02  /
   data green(40,40,10) /   1.7407340450302514D-02  /
   data green(40,40,11) /   1.7352201027886600D-02  /
   data green(40,40,12) /   1.7292408028858029D-02  /
   data green(40,40,13) /   1.7228112328965824D-02  /
   data green(40,40,14) /   1.7159473022240190D-02  /
   data green(40,40,15) /   1.7086656403205731D-02  /
   data green(40,40,16) /   1.7009834939254589D-02  /
   data green(40,40,17) /   1.6929186244342004D-02  /
   data green(40,40,18) /   1.6844892064473108D-02  /
   data green(40,40,19) /   1.6757137284616974D-02  /
   data green(40,40,20) /   1.6666108965739050D-02  /
   data green(40,40,21) /   1.6571995419619481D-02  /
   data green(40,40,22) /   1.6474985328050321D-02  /
   data green(40,40,23) /   1.6375266911908635D-02  /
   data green(40,40,24) /   1.6273027154509873D-02  /
   data green(40,40,25) /   1.6168451082580693D-02  /
   data green(40,40,26) /   1.6061721107171006D-02  /
   data green(40,40,27) /   1.5953016425868324D-02  /
   data green(40,40,28) /   1.5842512486795902D-02  /
   data green(40,40,29) /   1.5730380514078969D-02  /
   data green(40,40,30) /   1.5616787093755805D-02  /
   data green(40,40,31) /   1.5501893818496742D-02  /
   data green(40,40,32) /   1.5385856988973809D-02  /
   data green(40,40,33) /   1.5268827369294606D-02  /
   data green(40,40,34) /   1.5150949993573363D-02  /
   data green(40,40,35) /   1.5032364020454317D-02  /
   data green(40,40,36) /   1.4913202632220410D-02  /
   data green(40,40,37) /   1.4793592975008745D-02  /
   data green(40,40,38) /   1.4673656136603446D-02  /
   data green(40,40,39) /   1.4553507158281083D-02  /
   data green(40,40,40) /   1.4433255077233973D-02  /


end block data blkgreen
