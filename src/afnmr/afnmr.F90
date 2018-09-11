module comafnmr

!  some global variables go here

   implicit none
   integer,parameter::MAXAT=50000,MAXRES=10000,MAXIQM=600
   real(kind=8)::coord(3,MAXAT)
   real(kind=8):: qmcharge(MAXAT),rad(MAXAT)
   logical::connect(MAXRES,MAXRES)
   logical::atomsign(MAXAT)
   logical::ter(0:MAXRES)
   character(len=2):: element(MAXAT),dummyl
   character(len=3):: residue(MAXAT), residuename(MAXRES)
   character(len=4):: atomname(MAXAT)
   character(len=5):: dlabel(MAXAT)
!     here resno() goes sequentially from 1 to nres;
!          resno_user() are the residue numbers in the input pdb file
!          (note: does not recognize chainID's, so all resno_user() values
!           should be unique)
   integer ::resno(MAXAT), resno_user(MAXAT), list(MAXRES), prev_resno_user
   integer :: modnum, ier, listsize
   logical :: gaussian, orca, demon, demon5, qchem, terachem, &
              qopt, solinprot, xtb

end module

program afnmr_x

!   (Note: generally the input for afnmr.x is created by the shell script
!      "afnmr", which is in the AFNMRHOME/bin directory.)
!
!   Usage: afnmr.x program basis solinprot qopt basename list
!
!      where program  is G (Gaussian), O (Orca), D (Demon v3,4), D5 (Demon v5),
!                     Q (Qchem), T (TeraChem)
!            basis is D (double-zeta) or T (triple-zeta) 
!                  or M (primary res. T, rest D)
!            solinprot is T or F
!            qopt is T or F, to turn on or off quantum geometry optimization,
!                  or X (for internal optimization with xtb)
!            <basename>.pqr gives the input structure; put all "general"
!                  residues (water, ligands, etc.) after protein or
!                  nucleic acid residues
!            list is a list of residues to fragment; if empty, do all
!                residues up to lastprotres
!
      use comafnmr
      implicit none

      double precision :: a,b,c,d,e,x,y,z,tempdis
      integer :: ttnumber,strandno(MAXAT)
      character(len=8) :: lpchar
      character(len=1) :: program,basis,solinprotb,qoptb
      integer :: selectC(0:MAXRES+2),charge(MAXRES),cfrag
      integer :: selectCA(MAXRES+2), resmap(MAXRES)
      integer :: lastprotres
      integer :: i,j,k,ki,kk,kuser,m,iatfinal,iatstart,iitemp,iqm,iqmprot
      integer :: kfinal,ktemp,kkatom,kstart,kcount,kbas
      integer :: n1,n2,nprotc,nres,nsf,nptemp,nhighatom,nlowatom
      character(len=6) :: sn(MAXAT)
      character(len=3) :: rn
      character(len=1):: restype(MAXRES)
      character(len=30) :: pqrstart
      character(len=16) :: pqrend

      double precision total,pe,ee,ex,ec,dis,cuspfree,kinetic,nbcut
      double precision chargef(MAXRES)
      character(len=80) line,basename,pdbfile,filek,afnmrhome,version
      integer lengthb,lengthc,natom
#ifdef __INTEL_COMPILER
      integer system
#endif

      integer, parameter ::MAXNRES=31,MAXPRES=31
      character(len=3) :: nresn(MAXNRES), presn(MAXPRES)

      nresn(1) = '  G'
      nresn(2) = '  A'
      nresn(3) = '  U'
      nresn(4) = '  C'
      nresn(5) = ' DG'
      nresn(6) = ' DA'
      nresn(7) = ' DT'
      nresn(8) = ' DC'
      nresn(9) = 'GUA'
      nresn(10) = 'ADE'
      nresn(11) = 'THY'
      nresn(12) = 'CYT'
      nresn(13) = 'URA'
      nresn(14) = ' G5'
      nresn(15) = ' A5'
      nresn(16) = ' U5'
      nresn(17) = ' C5'
      nresn(18) = 'DG5'
      nresn(19) = 'DA5'
      nresn(20) = 'DT5'
      nresn(21) = 'DC5'
      nresn(22) = ' G3'
      nresn(23) = ' A3'
      nresn(24) = ' U3'
      nresn(25) = ' C3'
      nresn(26) = 'DG3'
      nresn(27) = 'DA3'
      nresn(28) = 'DT3'
      nresn(29) = 'DC3'
      nresn(30) = 'gtp'
      nresn(31) = 'MDA'

      presn(1) = 'ALA'
      presn(2) = 'ARG'
      presn(3) = 'ASH'
      presn(4) = 'ASN'
      presn(5) = 'ASP'
      presn(6) = 'CYS'
      presn(7) = 'CYX'
      presn(8) = 'GLH'
      presn(9) = 'GLN'
      presn(10) = 'GLU'
      presn(11) = 'GLY'
      presn(12) = 'HID'
      presn(13) = 'HIE'
      presn(14) = 'HIP'
      presn(15) = 'ILE'
      presn(16) = 'LEU'
      presn(17) = 'LYN'
      presn(18) = 'LYS'
      presn(19) = 'MET'
      presn(20) = 'PHE'
      presn(21) = 'PRO'
      presn(22) = 'SER'
      presn(23) = 'THR'
      presn(24) = 'TRP'
      presn(25) = 'TYR'
      presn(26) = 'VAL'
      presn(27) = 'HIS'
      presn(28) = 'ACE'
      presn(29) = 'NHE'
      presn(30) = 'NME'
      presn(31) = 'CYB'

      gaussian = .false.
      orca = .false.
      qchem = .false.
      demon = .false.
      demon5 = .false.
      terachem = .false.
      xtb = .false.
      solinprot = .false.
      qopt = .false.
      listsize = 0
      version = '5.4.beta1'

      print*
      print*,'**********************************************'
      print*,'AF-NMR'
      print*
      print*,'Authors : Xiao He, Kenneth Merz,'
      print*,'          Sishi Tang and David A. Case'
      print*,'Version: ', trim(version)
      print*,'**********************************************'

      call get_environment_variable('AFNMRHOME', afnmrhome )

      if( command_argument_count() .lt. 5 ) then
         write(6,*) &
       'Usage: af-nmr program basis solinprot qopt <basename> {list}'
         stop 1
      end if

      call get_command_argument( 1, program, lengthb )
      if( program .eq. 'G' ) then
         gaussian = .true.
      else if( program .eq. 'D' ) then
         demon = .true.
      else if( program .eq. 'E' ) then
         demon = .true.
         demon5 = .true.
      else if( program .eq. 'O' ) then
         orca = .true.
      else if( program .eq. 'Q' ) then
         qchem = .true.
      else if( program .eq. 'T' ) then
         terachem = .true.
      else
         write(0,*) 'Bad input for program: ', program
         stop 1
      endif

      call get_command_argument( 2, basis, lengthb )

      call get_command_argument( 3, solinprotb, lengthb )
      if( solinprotb .eq. 'T' .or. solinprotb .eq. 'S' ) solinprot = .true.

      call get_command_argument( 4, qoptb, lengthb )
      if( qoptb .eq. 'T' ) then
         qopt = .true.
      else if (qoptb .eq. 'X' ) then
         xtb = .true.
      else if (qoptb .ne. 'F' ) then
         write(0,*) 'Bad input for qopt: ', qoptb
         stop 1
      endif

      call get_command_argument( 5, basename, lengthb )
      pdbfile = basename(1:lengthb) // '.pqr'

      if( command_argument_count() > 5 ) then
         do i=6,command_argument_count()
            call get_command_argument( i, lpchar, lengthc )
            read( lpchar, '(i4)' ) list(i-5)
         end do
         listsize = command_argument_count() - 5
      endif

      nbcut = 3.3

      write(6,'(a,a)') 'Running afnmr, version ',trim(version)
      write(6,'(a,a1,1x,a1,1x,a1,1x,a1,1x,a)') 'Input arguments: ', &
           program, basis, solinprotb, qoptb, basename(1:lengthb)
      write(6,'(a,f6.3)') &
         'Fragments will be based on heavy atom contacts < ',nbcut

      do i=1,MAXRES
        ter(i) = .false.
        chargef(i) = 0.d0
      end do
      ter(0) = .false.
!
!     Read the input pqr file:
!
      open(10,file=pdbfile)

      nptemp = 0
      i = 0
      prev_resno_user = 0
      do k=1,MAXAT
        read(10,'(a80)',end = 101) line
        if( line(1:5) .eq. 'ATOM ' .or. line(1:5) .eq. 'HETAT') then
           i = i + 1
           read(line,100)sn(i),ttnumber,atomname(i),residue(i),  &
             resno_user(i),(coord(j,i),j=1,3),qmcharge(i),rad(i),element(i)
100        format(a6,1x,I4,1x,a4,1x,a3,2x,I4,4x,3f8.3,f8.4,f8.3,6x,a2)
           element(i) = adjustl(element(i))
           if( element(i)(1:1) .eq. 'E' ) then   ! skip "extra points"
!              .or. element(i).eq.'NA' .or. element(i).eq.'CL' &
!              .or. element(i).eq.'K ' .or. element(i).eq.'BR' ) then
              i = i - 1
              cycle
           endif
           if( resno_user(i) .ne. prev_resno_user ) then ! found new residue
              nptemp = nptemp + 1
              prev_resno_user = resno_user(i)
           endif
           resno(i) = nptemp
           resmap(resno_user(i)) = nptemp  ! maps user-resno to sequential ones
           chargef(nptemp) = chargef(nptemp) + qmcharge(i)
           if( nptemp.gt.MAXRES ) then
              write(0,*) 'too many residues: ', nptemp, MAXRES
              stop 1
           endif
           residuename(nptemp)=residue(i)

           rn = residue(i)
           restype(nptemp) = 'G'
           do m=1,MAXNRES
             if( rn.eq.nresn(m) ) restype(nptemp) = 'N'
           end do
           do m=1,MAXPRES
             if( rn.eq.presn(m) ) restype(nptemp) = 'P'
           end do

        else if( line(1:3) .eq. 'TER' ) then
!          mark residue nptemp as a terminal residue:
           ter(nptemp) = .true.
        endif
      end do
101   natom = i
      nres = max(1,nptemp)
      !  figure out the last protein residue:
      lastprotres = nres
      do i=1,nres
        if( restype(i) .eq. 'G' ) then
           lastprotres = i-1
           exit
        endif
      end do

!     lastprotres has to be a terminal residue:
      if( lastprotres > 0 ) ter(lastprotres) = .true.

      close(10)

!     real kludge, for pqr files with no element info and no 2-letter elements:
      do i=1,natom
        if( element(i).eq.'  ' ) then
           element(i)(1:1) = atomname(i)(2:2)
           element(i)(2:2) = ' '
        endif
      enddo
!
!     Populate arrays to identify C and CA atoms (for proteins); equivalent
!     atoms (O3' and C3') for nucleotides:
!
      do i=2,natom
        if( restype(resno(i)).eq.'P' ) then
           if (atomname(i).eq.' C  ')  selectC(resno(i))=i
           if (atomname(i).eq.' CA ')  selectCA(resno(i))=i
        else if( restype(resno(i)).eq.'N' ) then
           if (atomname(i).eq.' O3''')  selectC(resno(i))=i
           if (atomname(i).eq.' C3''')  selectCA(resno(i))=i
        else
           if ( resno(i-1).ne.resno(i) ) then
              selectC(resno(i))=i  ! selectC becomes first atom of a gen. res.
              selectC(resno(i-1)+1)=i  ! if residue numbers are not consecutive
           endif
        endif
      enddo
      if( residuename(lastprotres).eq.'NHE' .or. &
          residuename(lastprotres).eq.'NME' ) then
        selectC(lastprotres) = natom-2  !! will make kfinal=natom
      endif
      if( restype(1).eq.'G' ) selectC(1) = 1
      selectC(0) = 1
      selectC(nres+1) = natom+1
!
!     Set up the connectivity arrays:
!
      do i=1,nres
        do j=1,nres
          connect(i,j)=.false.
        enddo
      enddo

      do i=1,nres
        connect(i,i)=.true.
        charge(i) = nint( chargef(i) )
      enddo

!----------------------------------------------------------------------------

      do i=1,natom
        do j=i+1,natom

            if( resno(i).eq.resno(j) ) cycle

            dis=dsqrt((coord(1,i)-coord(1,j))**2   &
             +(coord(2,i)-coord(2,j))**2+(coord(3,i)-coord(3,j))**2)
!
!            nonbond distance < nbcut between heavy atom pairs
!
            if( (dis.le.nbcut .and. element(i)(1:1).ne.'H' .and.  &
                                    element(j)(1:1).ne.'H') )then

!                                   element(j)(1:1).ne.'H') .or.  &
!               (dis.le.(nbcut-0.5) .and. (element(i)(1:1).eq.'H' .or.  &
!                                         element(j)(1:1).eq.'H')) )then
#if 0
                write(6,'(a4,i4,5x,a4,i4,5x,f8.3)') &
                           atomname(i), resno_user(i), &
                           atomname(j),  resno_user(j), dis
#endif
                connect(resno(i),resno(j))=.true.
                connect(resno(j),resno(i))=.true.
!
!               need to make sure that net residue is also connected if 
!               one of the atoms is beyond selectC:
!
                if( i.ge.selectC(resno(i)) .and. restype(resno(i)).ne.'G') then
                   connect(resno(i)+1,resno(j))=.true.
                   connect(resno(j),resno(i)+1)=.true.
                endif
                if( j.ge.selectC(resno(j)) .and. restype(resno(j)).ne.'G') then
                   connect(resno(i),resno(j)+1)=.true.
                   connect(resno(j)+1,resno(i))=.true.
                endif
            endif

        enddo  !  j=i+1,natom
      enddo    !  i=1,natom

#if 0
      ! debug connectivity table:
      write(6,*) 'connectivities:'
      do i=1,nres
        do j=1,nres
          if( connect(i,j) ) write(6,*)  i,j, connect(i,j)
        end do
      end do
      stop 1
#endif

!
!     Big loop over residues to create fragments:
!
      do kcount=1,nres

        if( listsize > 0 ) then
           if( kcount > listsize ) exit
           kuser = list( kcount )
        else
           if( kcount > lastprotres ) exit
           kuser = kcount
        endif
        ! We need both "kuser" (user's file number) and "k" (a
        !    sequential residue number:
        k = resmap(kuser)

        iqm = 0   !  keep track of the atom number in the output file
        modnum = 0  ! get unique names for added hydrogens

        write(6,*)
        write(6,'(a,i4,a)') '==== Residue ',kuser, &
           ' ==========================================================='
        filek = basename(1:lengthb) // char(48+kuser/100)  &
            // char(48+(kuser-kuser/100*100)/10)  &
            // char(48+(kuser-kuser/10*10))
!
!       write header info:
!
        if( gaussian ) then
          open(30,file=filek(1:lengthb+3)//'.com')
        else if ( orca ) then
          open(30,file=filek(1:lengthb+3)//'.inp')
          open(32,file=filek(1:lengthb+3)//'.pos')
        else if ( demon ) then
          open(30,file=filek(1:lengthb+3)//'.inp')
        else if ( qchem ) then
          open(30,file=filek(1:lengthb+3)//'.in')
        else if ( terachem ) then
          open(30,file=filek(1:lengthb+3)//'.opt')
          open(32,file=filek(1:lengthb+3)//'.pos')
          open(34,file=filek(1:lengthb+3)//'.xyz1')
        end if
        if ( xtb ) then
          open(44,file=filek(1:lengthb+3)//'.xyz1')
        end if

        open(31,file=filek(1:lengthb+3)//'.pqr')
        if( solinprot ) then
           open(33,file=filek(1:lengthb+3)//'.prot.pqr')
        end if

        if( gaussian ) then
          write(30,'(a)') '%mem=800mw'
          write(30,'(a)') '%nprocshared=4'
          write(30,'(a)', advance='no')  '# OLYP/Gen charge nosymm '
          if (qopt) then
            write(30,'(a)') 'Opt ReadOptimize '
          else
            write(30,'(a)') 'nmr(printeigenvectors) integral(grid=ultrafine)'
          endif
          write(30,*)
          write(30,'(a,i4,a,a,a,f5.2)') ' AF-NMR fragment for residue ', &
               kuser, '; version = ',trim(version), '; nbcut = ', nbcut
          write(30,*)

        else if ( orca ) then
!         write(30,'(a)') '! PAL4'
          if( basis .eq. 'T' ) then
            write(30,'(a)', advance='no') '! OLYP TZVP '
          else
            write(30,'(a)', advance='no') '! OLYP VDZP '
          end if
          if( qopt ) then
            write(30,'(a)')  'TightSCF RI KDIIS Opt '
          else
            write(30,'(a)')  'TightSCF RI KDIIS '
          endif
          write(30,'(a)') ''
          write(30,'(a,a,a)') '%pointcharges "', filek(1:lengthb+3), &
               '.pos"'
          write(30,'(a)') ''
          write(30,'(a)') '%output'
          write(30,'(a)') '  PrintLevel Mini'
          write(30,'(a)') '  Print [ P_Cartesian ] 0'
          write(30,'(a)') '  Print [ P_Internal ] 0'
          write(30,'(a)') '  Print [ P_Basis ] 0'
          write(30,'(a)') '  Print [ P_Mulliken ] 0'
          write(30,'(a)') '  Print [ P_Loewdin ] 0'
          write(30,'(a)') '  Print [ P_Mayer ] 0'
          write(30,'(a)') '  Print [ P_OrbEn ] 0'
          write(30,'(a)') '  end'
          write(30,'(a)') ''
          write(30,'(a)') '%scf'
          write(30,'(a)') '  SOSCF Start 0.001 MaxIt 40 end'
          write(30,'(a)') '  DIISMaxEq 15'
          write(30,'(a)') '  directresetfreq 1'
          write(30,'(a)') '  end'
          write(30,'(a)') ''

        else if ( qchem ) then
          write(30,'(a)') '$comment'
          write(30,'(a,i4,a,a,a,f5.2)') ' AF-NMR fragment for residue ', &
               kuser, '; version = ',trim(version), '; nbcut = ', nbcut
          write(30,'(a)') '$end'
          write(30,'(a)') ''
          write(30,'(a)') '$molecule'

        else if ( demon ) then
          write(30,'(a,i4,a,a,a,f5.2)') 'TITLE AF-NMR fragment for residue ', &
               kuser, '; version = ',trim(version), '; nbcut = ', nbcut
          write(30,'(a)') 'SCFTYPE  RKS Tol=3.0e-6 MAX=200'
          write(30,'(a)') 'GUESS TB'
          write(30,'(a)') 'ORBITALS CARTESIAN'
          write(30,'(a)') 'ERIS MULTIPOLE'
          write(30,'(a)') 'VxcType Auxis OLYP'
          write(30,'(a)') 'GRID FINE'
          write(30,'(a)') 'QUADRATURE RANDOM'
          write(30,'(a)') 'MIXING -0.05'
          write(30,'(a)') 'SHIFT -0.2'
          write(30,'(a)') 'DIIS ON TOL=0.002'

        else if ( terachem ) then
          write(30,'(a)') 'basis 6-31gs'
          write(30,'(a,a)') 'coordinates ',filek(1:lengthb+3)//'.xyz'
          write(30,'(a,a)') 'pointcharges ',filek(1:lengthb+3)//'.pos'
          write(30,'(a,a)') 'scrdir ',filek(1:lengthb+3)
          write(30,'(a)') 'method pbe'
          write(30,'(a)') 'maxit 250'
          write(30,'(a)') 'levelshift yes'
          write(30,'(a)') 'levelshiftvala 0.2'
          write(30,'(a)') 'min_coordinates hdlc'
          write(30,'(a)') 'run minimize'
          write(30,'(a)') 'nstep 10'
          write(34,'(a)') 'put number of atoms here!'
          write(34,'(a)') filek(1:lengthb+3)

        end if

        if ( xtb ) then
          write(44,'(a)') 'put number of atoms here!'
          write(44,'(a,i4,a,a,a,f5.2)') 'AF-NMR fragment for residue ', &
               kuser, '; version = ',trim(version), '; nbcut = ', nbcut
        endif

        do i=1,natom
          atomsign(i)=.false.
        enddo

!       cycle through all "connected" fragments to get the charge on the
!       quantum region (cfrag), and to mark each quantum atom by atomsign(i)

        cfrag = 0
        do ktemp=1,nres
          if(connect(k,ktemp))then
            call get_atom_range( kstart, kfinal, ktemp, &
                   selectC, restype(ktemp), ter)
            atomsign(kstart:kfinal) = .true.
            cfrag = cfrag + charge(ktemp)
            write(6,'(20x,i5,a4,5i6)') kuser,residuename(ktemp), &
                 resno_user(kfinal),kstart,kfinal,charge(ktemp),cfrag
          endif
        enddo
        write(6,'(2i5)') kuser, cfrag

        if ( gaussian .or. qchem ) then
          write(30,50)cfrag,1
        else if ( orca ) then
          write(30,'(a,i3,1x,i3)')'* xyz ',cfrag,1
        else if ( demon ) then
          write(30,'(a,i3)')'CHARGE  ',cfrag
          write(30,'(a,i3)')'MULTIPLCITY  ',1
          write(30,'(a)') 'GEOMETRY CARTESIAN ANGSTROM'
        else if ( terachem ) then
          write(30,'(a,i3)')'charge  ',cfrag
          write(30,'(a)') 'end'
        end if

        nhighatom=0
        nlowatom=0
!
!       get starting, ending atoms for residue "k":
!
        call get_atom_range( iatstart, iatfinal, k, &
            selectC, restype(k), ter )
!
!       write out atoms in principal residue:
!
        do kk=iatstart,iatfinal
          iqm = iqm + 1
          nhighatom = nhighatom + 1
          call addatom( kk, iqm, basis)
        enddo

        do ktemp=1,nres
          if( ktemp.ne.k .and. connect(k,ktemp) )then

              call get_atom_range( kstart, kfinal, ktemp, &
                   selectC, restype(ktemp), ter)
!
!             write out coordinates in a "connected" residue
!
              do kk=kstart,kfinal
                nlowatom = nlowatom + 1
                iqm = iqm + 1
                call addatom( kk, iqm, ' ')
              enddo
!
!             add hydrogens to dangling residues:
!
              if(ktemp.ne.1 .and. ktemp.le.lastprotres  &
                  .and. .not. ter(ktemp-1) )then
!
!               ---look for "backwards" links to previous residue:
!
                if(.not. atomsign(kstart-1) )then
                  n1=kstart
                  n2=selectCA(ktemp-1)
                  if(atomname(n1).eq.' C  ') then
                    call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2),  &
                      coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
                  else
                    call xyzchangeO(coord(1,n2),coord(2,n2),coord(3,n2), &
                      coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
                  endif
                  nlowatom = nlowatom + 1
                  iqm = iqm + 1
                  call addH( iqm, x, y, z)
                endif
              endif

              if(ktemp.lt.lastprotres .and. .not. ter(ktemp)) then
!
!             ---look for "forwards" links to the next residue:
!
                if(.not. atomsign(kfinal+3) )then
                  n1=selectCA(ktemp)
                  n2=selectC(ktemp)
                  call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2),  &
                    coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
                  nlowatom = nlowatom + 1
                  iqm = iqm + 1
                  call addH( iqm, x, y, z)
                endif
              endif

          endif
          if( ktemp == lastprotres ) iqmprot = iqm  ! only freeze protein
        enddo  ! ktemp

        !  check for dangling S--S bonds:

        do kk=1,natom
          if( atomname(kk)(2:3) == 'SG' .and. atomsign(kk) ) then
!           write(0,*) 'found SG: ', resno(kk), atomsign(kk)
            do i=1,natom
              if( i.eq.kk) cycle
              if( atomname(i)(2:3) .ne. 'SG') cycle
              tempdis=dsqrt((coord(1,i)-coord(1,kk))**2.0d0  &
                 +(coord(2,i)-coord(2,kk))**2.0d0            &
                 +(coord(3,i)-coord(3,kk))**2.0d0)
              if(tempdis.le.2.5d0 .and. kk.ne.i   &
                   .and.  .not. atomsign(i) )then
                !we have found a "dangling" S--S bond; saturate:
                call xyzchangeS(coord(1,i),coord(2,i),coord(3,i),  &
                    coord(1,kk),coord(2,kk),coord(3,kk),x,y,z)
                iqm = iqm + 1
                call addH( iqm, x, y, z)
              endif
            enddo
          endif
        enddo
        close(31)

!       Deal with the charges representing the protein and solvent polarization:

        if( gaussian) then
          write(30,*)
        else if ( qchem ) then
          write(30,'(a)') '$end'
          write(30,'(a)') ''
          write(30,'(a)') '$external_charges'
        else if ( demon ) then
          write(30,'(a)') 'END'
          write(30,'(a)') 'EMBED READ'
        endif

!       write out the protein positions and charges:

        do kk=1,natom
          if( .not. atomsign(kk) ) then
            if( solinprot ) then
               if ( restype(resno(kk)).ne.'G' ) &
               write(33,'(a,i5,1x,a4,1x,a3,i6,4x,3f8.3,f8.4,f8.3)') 'ATOM  ', &
               kk,atomname(kk),residue(kk),resno_user(kk),(coord(j,kk),j=1,3), &
               qmcharge(kk),rad(kk)
            else if( gaussian .or. qchem .or. demon ) then
              write(30,1315)(coord(j,kk),j=1,3),qmcharge(kk)
            endif
          endif
        enddo

        if( solinprot ) then
          close(33)
#if 0
          ier = system('./runsolinprot ' // filek(1:lengthb+3))
#else  /* use this for GCC version 4.6 and above, esp. if above fails */
          call execute_command_line('./runsolinprot ' // filek(1:lengthb+3), &
               exitstat = ier)
#endif
          if( ier .ne. 0 )then
             write(0,*) "error in solinprot: check solinprot.out"
             write(0,*) "error code was ", ier
             write(0,*) "command line was:"
             write(0,*) './runsolinprot ' // filek(1:lengthb+3)
             stop 1
          end if
        endif

        if( orca ) then
          nsf=0  ! number of surface points, just count them here
          open(23,file='srfchg.pos')
          do iitemp=1,999999
            read(23,*,end=59)a,b,c,d
            nsf = nsf + 1
          enddo
59        continue
          close(23)

          nprotc = 0
          if( .not.solinprot ) then
             do kk=1,natom
               if(.not. atomsign(kk)) nprotc = nprotc + 1
             enddo
          endif
          write(32,'(i7)') nprotc+nsf
          do kk=1,natom
             if(.not. atomsign(kk) )then
                write(32,1317)qmcharge(kk), (coord(j,kk),j=1,3)
             endif
          end do
        endif

        open(23,file='srfchg.pos')
        do iitemp=1,999999
          read(23,*,end=60)a,b,c,d
          if( gaussian .or. qchem .or. demon ) then
            write(30,1316)a,b,c,d
          else if( orca .or. terachem) then
            write(32,1318)d,a,b,c
          endif
        enddo
60      continue
        close(23)

        !  More program-dependent keywords and instructions:

        if( gaussian ) then
          write(30,*)
          if( basis .eq. 'M' ) then
            !    write local basis set
            do i=1,nhighatom
              write(30,'(i3,a)') i,' 0'
              open( UNIT=11, FILE=trim(afnmrhome) // '/basis/pcSseg-1/' &
                    // trim(dlabel(i)) // '.gbs')
              rewind(11)
              do kbas=1,9999
                 read(11,'(a80)',end = 62) line
                 write(30,'(a)') trim(line)
              end do
              close(11)
62            write(30,'(a)') '****'
            end do
            if(nhighatom .ge. 9) then
               if((nhighatom+nlowatom) .ge. 100) then
                 write(30,'(i2,a1,i3,a2)') nhighatom+1,'-',nhighatom+nlowatom,' 0'
               else if((nhighatom+nlowatom) .lt. 100) then
                 write(30,'(i2,a1,i2,a2)') nhighatom+1,'-',nhighatom+nlowatom,' 0'
               end if
            else if(nhighatom .ne. 0) then
               if((nhighatom+nlowatom) .ge. 100) then
                 write(30,'(i1,a1,i3,a2)') nhighatom+1,'-',nhighatom+nlowatom,' 0'
               else if((nhighatom+nlowatom) .lt. 100) then
                 write(30,'(i1,a1,i2,a2)') nhighatom+1,'-',nhighatom+nlowatom,' 0'
               end if
            else
               if((nhighatom+nlowatom) .ge. 100) then
                  write(30,'(i1,a1,i3,a2)') 5,'-',nhighatom+nlowatom,' 0'
               else if((nhighatom+nlowatom) .lt. 100) then
                  write(30,'(i1,a1,i2,a2)') 5,'-',nhighatom+nlowatom,' 0'
               end if
            end if
            write(30,'(A)') 'SVP'
            write(30,'(A)') '****'
            write(30,*)

          else if( basis .eq. 'T' ) then
            open( UNIT=11, FILE=trim(afnmrhome) // &
                '/basis/pcSseg-1/pcSseg-1.gbs')
            rewind(11)
            do kbas=1,9999
               read(11,'(a80)',end = 63) line
               write(30,'(a)') trim(line)
            end do
            close(11)
63          continue

          else if( basis .eq. 'D' ) then
            open( UNIT=11, FILE=trim(afnmrhome) // &
                '/basis/pcSseg-0/pcSseg-0.gbs')
            rewind(11)
            do kbas=1,9999
               read(11,'(a80)',end = 64) line
               write(30,'(a)') trim(line)
            end do
            close(11)
64          continue

          endif

          if (qopt) then
            write(30,'(a,i4,a,i4,a,i4)') 'noatoms  atoms=1-', nhighatom, &
                 ', ', iqmprot+1, '-', natom
            write(30,*)
          endif

        else if ( demon ) then
          write(30,'(a)') 'END'
          if( nhighatom == 0 ) nhighatom = 4
          if( basis .eq. 'M' )then
            write(30,'(a)') 'BASIS  (DZVP)'
            do i=1,nhighatom
              write(30,'(a,a)') dlabel(i),' (pcSseg-1)'
            end do
            write(30,'(a)') 'AUXIS  (A2)'
            do i=1,nhighatom
              if( dlabel(i)(1:1) .eq. 'H' ) then
                write(30,'(a,a)') dlabel(i),' (GEN-A2)'
              else
                write(30,'(a,a)') dlabel(i),' (GEN-A2*)'
              endif
            end do
          else if (basis .eq. 'T' ) then
            write(30,'(a)') 'BASIS  (pcSseg-1)'
            write(30,'(a)') 'AUXIS  (GEN-A2*)'
          else if (basis .eq. 'D' ) then
            write(30,'(a)') 'BASIS  (pcSseg-0)'
            write(30,'(a)') 'AUXIS  (GEN-A2)'
          endif

          if( demon5 ) then
            write(30,'(a)') 'NMR READ SHIELDING'
          else
            write(30,'(a)') 'NMR READ'
          endif
          do i = 1,nhighatom
            write(30,'(a,a)', advance='no' ) dlabel(i),' '
            if( mod(i,12).eq.0 .and. i.ne.nhighatom ) write(30,'(a)') ' &'
          end do
          write(30,'(a)') ' '
          write(30,'(a)') 'PRINT NMR'

          if( qopt ) then
            write(30,'(a)') 'OPTIMIZE CARTESIAN MAX=5'
            write(30,'(a)') 'CONSTANTS'
            do i=nhighatom+1,iqmprot
               write(30,'(a,a)') dlabel(i), '  XYZ'
            end do
          endif

        else if( orca ) then
          write(30,'(a)') '*'
          if( qopt ) then
            write(30,'(a)') '%geom MaxIter=5'
            write(30,'(a)') '      Constraints'
            do i=nhighatom+1,iqmprot
               write(30,'(a,i4,a)') '        { C ', i-1, ' C }'
            end do
             write(30,'(a)') '      end'
             write(30,'(a)') '    end'
          endif
          !  next two lines are for versions of Orca up to 3.0.1
          ! write(30,'(a)') '%eprnmr  ori IGLO'
          ! write(30,'(a)') '   LocMet PM'
          !  following line is for Orca 4:
          write(30,'(a)') '%eprnmr  ori GIAO'
          write(30,'(a)', advance='no') '   nuclei = 1'
          do i=2,nhighatom
             write(30,'(a1,i3)', advance='no') ',', i
          end do
          write(30,'(a)') ' {shift}'
          write(30,'(a)') 'end'

        else if ( qchem ) then
          write(30,'(a)') '$end'
          write(30,'(a)') '   '
          write(30,'(a)') '$rem'
          write(30,'(a)') 'JOBTYPE       NMR'
          write(30,'(a)') 'EXCHANGE      OPTX'
          write(30,'(a)') 'CORRELATION   LYP'
          write(30,'(a)') 'BASIS         def2-SVP'
          write(30,'(a)') 'SYMMETRY      FALSE'
          write(30,'(a)') 'WAVEFUNCTION_ANALYSIS  FALSE'
          write(30,'(a)') '$end'

        else if ( terachem ) then  !  terachem is only for qopt calcs.
          write(30,'(a)') '$constraints'
          do i=nhighatom+1,iqmprot
             write(30,'(a,i4)') '  atom ', i
          end do
          write(30,'(a)') '$end'

        end if

        if ( xtb ) then  !  xtb is only for qopt calcs.
          write(44,'(a)') '$set'
          write(44,'(a,i4,a,i4)') 'fix ',nhighatom+1,'-',iqmprot
          write(44,'(a)') 'fixfc 0'
          write(44,'(a)') '$end'
          close(44)
        end if

        close(30)
        if( orca ) close(32)
        if( terachem ) then
          close(32)
          close(34)
          open(34,file=filek(1:lengthb+3)//'.xyz1')
          open(35,file=filek(1:lengthb+3)//'.xyz')
          read(34,*)  ! skip dummy first line
          write(35,'(i4)') iqm   ! number of atoms goes on first line
          do i=1,9999
            read(34,'(a)', end=105) line
            write(35,'(a)') trim(line)
          end do
  105     close(34)
          close(35)
          call execute_command_line( '/bin/rm -f ' // filek(1:lengthb+3) &
                // '.xyz1' )
        endif
        if( xtb ) then
          open(44,file=filek(1:lengthb+3)//'.xyz1')
          open(45,file=filek(1:lengthb+3)//'.xyz2')
          read(44,*)  ! skip dummy first line
          write(45,'(i4)') iqm   ! number of atoms goes on first line
          do i=1,9999
            read(44,'(a)', end=205) line
            write(45,'(a)') trim(line)
          end do
  205     close(44)
          close(45)
          call execute_command_line( '/bin/rm -f ' // filek(1:lengthb+3) &
                // '.xyz1' )

          !  do the xtb minimization here:
          write(6,*)
          write(6,*) 'Doing geometry optimization with xtb'
          call execute_command_line( 'xtb ' // filek(1:lengthb+3) &
                // '.xyz2 -opt > ' // filek(1:lengthb+3) // '.xtb.log' )
          call execute_command_line( '/bin/rm -f ' // filek(1:lengthb+3) &
                // '.xyz2' )
          call execute_command_line( &
            '/bin/rm -f energy charges wbo xtbrestart xtbopt.log' )

          !  extract the coordinates from the xtb output file:
          open(46,file='xtbopt.xyz')
          read(46,*) iqm
          read(46,*)   ! title line
          do i=1,iqm
             read(46,*) dummyl, coord(1,i), coord(2,i), coord(3,i)
          end do
          close(46)
          call execute_command_line( '/bin/rm -f xtbopt.xyz' )

          !  transfer the minimized coordinates to the .pqr file
          open(47,file=filek(1:lengthb+3)//'.pqr')
          open(48,file=filek(1:lengthb+3)//'.pqr1')
          do i=1,iqm
             read(47,'(a30,24x,a16)') pqrstart,pqrend
             write(48,'(a30,3f8.3,a16)') pqrstart, &
                 coord(1,i), coord(2,i), coord(3,i), pqrend
          end do
          close(47)
          close(48)
          call execute_command_line( '/bin/mv ' // filek(1:lengthb+3) &
             // '.pqr1 ' // filek(1:lengthb+3) // '.pqr' )

          if( demon ) then

             !  transfer the minimized coordinates to the .inp file
             open(47,file=filek(1:lengthb+3)//'.inp')
             open(48,file=filek(1:lengthb+3)//'.inp1')
             do i=1,9999
                read(47,'(a)', end=210) line
                write(48,'(a)') trim(line)
                if( line(1:9) == 'GEOMETRY ' ) then
                   do j=1,iqm
                      read (47,*) 
                      write(48,'(a,2x,3f12.5)') dlabel(j), &
                             coord(1,j), coord(2,j), coord(3,j)
                   end do
                end if
             end do
  210        close(47)
             close(48)
             call execute_command_line( '/bin/mv ' // filek(1:lengthb+3) &
                // '.inp1 ' // filek(1:lengthb+3) // '.inp' )

          else
             write(6,*) "xtb only works with demon for now"
             call exit(1)
          end if

        end if

1315    format(3f10.4,2x,f12.8)
1317    format(f12.8,2x,3f10.4)
50      format(1x,I3,2x,I2)
1316    format(3f15.6,2x,f12.8)
1318    format(f12.8,2x,3f15.6)

        if( solinprot ) write(0,'(a,i4)') '    done with residue ', kuser

      enddo  ! big loop over residues

end program afnmr_x

subroutine xyzchange(xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew)

        implicit none
        double precision xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew,grad

         grad=dsqrt(1.09d0**2/((xold-xzero)**2+(yold-yzero)**2   &
           +(zold-zzero)**2))
         xnew=xzero+grad*(xold-xzero)
         ynew=yzero+grad*(yold-yzero)
         znew=zzero+grad*(zold-zzero)
         return
end subroutine xyzchange

subroutine xyzchangeS(xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew)

        implicit none
        double precision xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew,grad

         grad=dsqrt(1.33d0**2/((xold-xzero)**2+(yold-yzero)**2   &
           +(zold-zzero)**2))
         xnew=xzero+grad*(xold-xzero)
         ynew=yzero+grad*(yold-yzero)
         znew=zzero+grad*(zold-zzero)
         return
end subroutine xyzchangeS

subroutine xyzchangeO(xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew)

        implicit none
        double precision xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew,grad

         grad=dsqrt(0.96d0**2/((xold-xzero)**2+(yold-yzero)**2  &
           +(zold-zzero)**2))
         xnew=xzero+grad*(xold-xzero)
         ynew=yzero+grad*(yold-yzero)
         znew=zzero+grad*(zold-zzero)
         return
end subroutine xyzchangeO

subroutine get_atom_range( kstart, kfinal, ktemp,  &
           selectC, restype, ter )
        implicit none

        integer, intent(in)   :: ktemp, selectC(0:*)
        integer, intent(out)  :: kstart, kfinal
        logical, intent(in)   :: ter(0:*)
        character(len=1), intent(in) ::  restype

        if( restype.eq.'P' ) then
          kstart=selectC(ktemp-1)
          kfinal=selectC(ktemp)-1
          if( ter(ktemp) )then
            kfinal=selectC(ktemp)+2
          endif
          if( ter(ktemp-1) ) then
            kstart=selectC(ktemp-1)+3
          endif
        else if( restype.eq.'N' ) then
          kstart=selectC(ktemp-1)
          kfinal=selectC(ktemp)-1
          if( ter(ktemp) )then
            kfinal=selectC(ktemp)+1
          endif
          if( ter(ktemp-1) ) then
            kstart=selectC(ktemp-1)+2
          endif
        else    !  general (G) residue type: should be ligand or water
          kstart=selectC(ktemp)
          kfinal=selectC(ktemp+1)-1
        endif

        return
end subroutine get_atom_range

subroutine addatom( kk, iqm, basis )

      use comafnmr
      implicit none
      integer, intent(in) ::  kk,iqm
      character*1, intent(in) ::  basis
      integer j
      character*3  i_char

      if ( demon ) then
        write( i_char, '(i3)' ) iqm
        dlabel(iqm) = trim(element(kk)) // adjustl(i_char)
        if( len_trim(element(kk)) == 1 ) dlabel(iqm)(5:5) = ' ' 
        write(30,'(a,2x,3f12.5)') dlabel(iqm),(coord(j,kk),j=1,3)
      else if ( terachem ) then
        write(34,1000)element(kk),(coord(j,kk),j=1,3)
      else if ( gaussian ) then
        write(30,1000)element(kk),(coord(j,kk),j=1,3)
        dlabel(iqm) = trim(element(kk))
      else
        write(30,1000)element(kk),(coord(j,kk),j=1,3)
        if( orca .and. basis .eq. 'M' ) then
          write(30,'(a)') 'NewGTO'
          write(30,'(a)') '"TZVP"'
          write(30,'(a)') 'end;'
        endif
      endif
      if( xtb ) then
        write(44,1000)element(kk),(coord(j,kk),j=1,3)
      endif

      write(31,'(a,i5,1x,a4,1x,a3,i6,4x,3f8.3,f8.4,f8.3)') 'ATOM  ', &
        kk,atomname(kk),residue(kk),resno_user(kk),(coord(j,kk),j=1,3),   &
        qmcharge(kk),rad(kk)

      return
1000  format(a2,4x,3f10.4)
end subroutine addatom

subroutine addH( iqm, x, y, z)

      use comafnmr
      implicit none
      integer, intent(in) :: iqm
      real(kind=8), intent(in) ::  x,y,z
      character(len=3) i_char

      if( demon ) then
        write( i_char, '(i3)' ) iqm
        dlabel(iqm) = 'H' // adjustl(i_char)
        dlabel(iqm)(5:5) = ' ' 
        write(30,'(a,2x,3f12.5)') dlabel(iqm),x,y,z
      else if ( terachem ) then
        write(34,1000)'H ',x,y,z
      else if ( gaussian ) then
        write(30,1000)'H ',x,y,z
        dlabel(iqm) = 'H'
      else
        write(30,1000)'H ',x,y,z
      endif
      if( xtb ) then
        write(44,1000)'H ',x,y,z
      end if

      modnum = modnum + 1
      if( modnum < 10 ) then
         write(31,'(a,i5,2x,a,i1,a,3f8.3,f8.4,f8.3)') 'ATOM  ',  &
           iqm,'H',modnum,'  MOD  9999    ',x,y,z,0.0,1.2
      else
         write(31,'(a,i5,2x,a,i2,a,3f8.3,f8.4,f8.3)') 'ATOM  ',  &
           iqm,'H',modnum,' MOD  9999    ',x,y,z,0.0,1.2
      endif

      return
1000    format(a2,4x,3f10.4)
end subroutine addH

