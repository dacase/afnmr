#include "../include/dprec.fh"

!The 3D-RISM-KH software found here is copyright (c) 2012 by
!Andriy Kovalenko, Tyler Luchko and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999);
!ibid. 112:10391-10417 (2000).
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010).


!> This program takes 3D solvent distribution functions from
!! 3D-RISM calculations and calculates thermodynamic properties.
module thermo_m
  use safemem
  use getopts_c
  use rism_report_c
  use rism3d_solvent_c
  use rism3d_solute_c
  use rism3d_potential_c
  use rism3d_grid_c
  use rism3d_closure_c
  implicit none

  integer, parameter :: CLEN = 1024

  !>  Filename for the input guv file.
  character(len = CLEN), pointer :: guvfile(:) => NULL()
  !> Filename for the input huv file.
  character(len = CLEN), pointer :: huvfile(:) => NULL()
  !> Filename for the input cuv file.
  character(len = CLEN), pointer :: cuvfile(:) => NULL()
  !> Filename for the input hLongRangeAsymR file.
  character(len = CLEN) :: ahrfile = ''
  !> Filename for the input dcfLongRangeAsympR file.
  character(len = CLEN) :: acrfile = ''
  !> Filename for the input xvvfile file.
  character(len = CLEN) :: xvvfile = ''

  !> Solute parameter file.
  character(len = CLEN) :: parmfile
  !> Solute coordinate file.
  character(len = CLEN) :: crdfile

  !> Indicates a volume to integrate over.
  character(len = CLEN), pointer :: volspec(:) => NULL()

  !> Name of the closure to use for calculations.
  character(len = CLEN) :: closureName = "KH"

  type boxside
     integer :: start, stop, stride
  end type boxside
  !> Box volumes to integrate over (dim, number) with each dimension.
  !! Has properties (start, stop, stride).
  type(boxside), pointer :: volbox(:,:)

  !> Origin from the DX file.  Tells us how to translate the solute and
  !! how to interpret volume specifications
  _REAL_ :: dxOrigin(3)

  !> Cut-off used to accelerate calculations.
  _REAL_ :: cut = huge(1d0)

  !> Solvent object.
  type(rism3d_solvent), save :: solvent
  !> Solute object.
  type(rism3d_solute), save :: solute

  !> Potential object.
  type(rism3d_potential), save :: potential
  !> Grid object.
  type(rism3d_grid), save :: grid
  !> Closure object.
  type(rism3d_closure), save :: closure

  !> Solute-solvent radial distribution function. (product(nr), nsites) 
  _REAL_, pointer :: guv(:,:) => NULL()
  !> Solute-solvent total correlation function. (product(nr), nsites)
  _REAL_, pointer :: huv(:,:) => NULL()
  !> Solute-solvent direct correlation function. (nr(1), nr(2), nr(3), nsite)
  _REAL_, pointer :: cuv(:,:,:,:) => NULL()

  !> Number results (xvv%ntype,# volspec).
  _REAL_, pointer :: numRes(:,:) => NULL()
  !> Number results with long range correction (xvv%ntype, # volspec).
  _REAL_, pointer :: numResLR(:,:) => NULL()

  !> PARTIALMOLARVOLUME results (# volspec).
  _REAL_, pointer :: partialMolarVolumeRes(:) => NULL()

  !> ExcessChemicalPotential results (xvv%ntype,# volspec).
  _REAL_, pointer :: excessChemicalPotentialRes(:,:) => NULL()
  !> ExcessChemicalPotential results with long range correction (xvv%ntype, # volspec).
  _REAL_, pointer :: excessChemicalPotentialResLR(:,:) => NULL()

  ! Calculation flags.
  logical :: all = .false.
  logical :: excessChemicalPotential = .false.
  logical :: number = .false.
  logical :: partialMolarVolume = .false.
  logical :: force = .false.
  !> Perform long range asymptotics corrected calculation.
  logical :: longRange = .false.

  
contains
  

  !> Get command line options.
  subroutine getOptions()
    implicit none
    integer :: err,i
    character(len=CLEN),pointer :: extra(:) => NULL()
    character(len=10*CLEN), pointer :: tempC => NULL()
    !closureOrder : order of the closure, only relevent for PSE-n
    integer :: closureOrder=1

    call getopts_add("all",all)
    call getopts_add("exchem",excessChemicalPotential)
    call getopts_add("number",number)
    call getopts_add("pmv",partialMolarVolume)
    call getopts_add("force",force)
    call getopts_add("vol",'')
    call getopts_add("xvv",'','x',min=1,max=1)
    call getopts_add("guv",'','g')
    call getopts_add("huv",'','h')
    call getopts_add("cuv",'','c')
    call getopts_add("ahr",'',min=0,max=1)
    call getopts_add("acr",'',min=0,max=1)
    call getopts_add("parm",'',min=0,max=1)
    call getopts_add("crd",'',min=0,max=1)
    err = getopts_process()
    if (err < 0) then
       call usage()
       call rism_report_error('(a,i3)', "command line read failed: ",err)
    end if
    call getopts_get("all",all)
    call getopts_get("exchem",excessChemicalPotential)
    call getopts_get("number",number)
    call getopts_get("pmv",partialMolarVolume)
    call getopts_get("force",force)
    if (all) then
       excessChemicalPotential=.true.
       number=.true.
       partialMolarVolume=.true.
       force=.true.
    end if
    call getopts_get("xvv",1,xvvfile)
    call getopts_get("ahr",1,ahrfile)
    call getopts_get("acr",1,acrfile)
    call getopts_get("parm",1,parmfile)
    call getopts_get("crd",1,crdfile)
    guvfile => getopts_getAll("guv",guvfile,CLEN)
    huvfile => getopts_getAll("huv",huvfile,CLEN)
    cuvfile => getopts_getAll("cuv",cuvfile,CLEN)
    volspec => getopts_getAll("vol",volspec,CLEN)
    extra => getopts_unread()
    if (size(extra) /= 0) then
       call usage()
       !concatenate all extra options into temp
       allocate(tempC)
       tempC=""
       do i=1, size(extra)
          tempC(len_trim(tempC):) =  trim(extra(i))
       end do
       write(0,'(a)') "ERROR: unknown options:"//tempC
    end if
    err = safemem_dealloc(extra)
    call sanityCheck1()
    call getopts_cleanup()
    !    call process_volspec()

  end subroutine getOptions


  !> Check user specified options.
  subroutine sanityCheck1()
    implicit none

    ! For each option ensure we have the necessary files.  We will
    ! check the numbers of files later.
    if (number) then
       if (size(guvfile)==0 .and. size(huvfile)==0) then
          call rism_report_error("--guv or --huv files for all solvent types required for --number calculations")
       end if
    end if
    if (excessChemicalPotential) then
       if (size(guvfile)==0 .and. size(huvfile)==0) then
          call rism_report_error("--guv or --huv files for all solvent types required for --exchem calculations")
       end if
       if (size(cuvfile)==0) then
          call rism_report_error("--cuv files for all solvent types required for --exchem calculations")
       end if
    end if
    if (partialMolarVolume) then
       if (size(cuvfile)==0) then
          call rism_report_error("--cuv files for all solvent types required for --pmv calculations")
       end if
    end if
    if (force) then
       call rism_report_error("--force calculations not supported at this time")
    end if
    ! Check solute input files.
    if (len_trim(parmfile) == 0 .and. len_trim(crdfile) /= 0) then
       call rism_report_error("--crd requires --parm")
    else if (len_trim(parmfile) /= 0 .and. len_trim(crdfile) == 0) then
       call rism_report_error("--parm requires --crd")
    end if
  end subroutine sanityCheck1


  !> Check we have the appropriate files after they have been read in.
  !! Call readFiles() first.
  subroutine sanityCheck2()
    implicit none
    integer :: i,j,k
    if (number) then
       if (.not. associated(guv)) then
          guv => safemem_realloc(guv,ubound(huv,1),ubound(huv,2),.false.)
          guv = huv + 1d0
       end if
    end if
    if (excessChemicalPotential) then
       if (.not. associated(huv)) then
          huv => safemem_realloc(huv,ubound(guv,1),ubound(guv,2),.false.)
          huv = guv - 1d0
       end if
    end if
    ! If we have integration volumes, check that they fit.  If there
    ! are none, create one.
    if (associated(volbox)) then
       do j=1,ubound(volbox,2)
          do i=1,3
             if (volbox(i,j)%start <1 .or. volbox(i,j)%stop > grid%localDimsR(i)) then
                write(0,'(a)') "ERROR: integration volume out-of-bounds"
                write(0,'(a)') "Specification         Actual Size"
                do k=1,3
                   write(0,'(a1,i3,a1,i3,a1,i3,10x,i3,a1,i3)') "x",&
                        volbox(k,j)%start,":",volbox(k,j)%stop,":",volbox(k,j)%stride, &
                        1,":",grid%localDimsR(k)
                end do
                stop
             end if
          end do
       end do
    else
       allocate(volbox(3,1))
       do i=1,3
          volbox(i,1)%start=1
          volbox(i,1)%stop=grid%localDimsR(i)
          volbox(i,1)%stride=1
       end do
    end if
  end subroutine sanityCheck2


  !> Read in input files.
  subroutine readFiles()
    use rism3d_opendx
    implicit none
    integer :: i
    _REAL_, pointer :: tmp1(:)
    _REAL_ :: delta(3)
    integer, parameter :: zeroR3(3) = 0d0
    integer :: npos(3), nkpos(3), tmppos(3)

    !Xvv is always present
    call rism3d_solvent_new(solvent,xvvfile)

    !Solute molecule
    if (len_trim(parmfile) /= 0 .and. len_trim(crdfile) /= 0) then
       call rism3d_solute_new(solute,parmfile,solvent%temperature)
       call rism3d_solute_setCoord(solute,crdfile)
    end if

    !initialize objects
    call rism3d_grid_new(grid)
    !FIXME: temporarily disabled for periodic code
    ! call rism3d_potential_new(potential,grid,solv,solu,cut,NULL())
    call rism3d_closure_new(closure,closureName,potential)

    !get grid size
    if (ubound(guvfile,1) >0) then
       call readDXHeader(guvfile(1),dxOrigin,delta,npos)
    else if (ubound(huvfile,1) >0) then
       call readDXHeader(huvfile(1),dxOrigin,delta,npos)
    else if (ubound(cuvfile,1) >0) then
       call readDXHeader(cuvfile(1),dxOrigin,delta,npos)
    else if (len_trim(ahrfile) >0) then
       call readDXHeader(ahrfile,dxOrigin,delta,npos)
    else if (len_trim(acrfile) >0) then
       call readDXHeader(acrfile,dxOrigin,delta,npos)
    end if

    !set grid size
    nkpos = npos
    nkpos(1) = nkpos(1)+2
    call rism3d_grid_resize(grid,delta,npos,nkpos,npos,nkpos,zeroR3,zeroR3)
    !Guv
    if (ubound(guvfile,1) > 0) then
       guv => safemem_realloc(guv,grid%totalLocalPointsK,solvent%numAtomTypes,.false.)
       if (solvent%numAtomTypes /= ubound(guvfile,1)) &
            call rism_report_error('(a,i4,a,i4)',"Different number of --guv files (",ubound(guvfile,1)&
            ,"than solvent sites", solvent%numAtomTypes)
       do i = 1, solvent%numAtomTypes
          !check file size
          call readDXHeader(guvfile(i),dxOrigin,delta,tmppos)
          if (sum(abs(tmppos- npos))/=0) call rism_report_error(trim(guvfile(i))//" is the wrong size")
          call readDX(guvfile(i),guv(:,i),grid%localDimsR,dxOrigin,delta)
       end do
    end if
    if (ubound(huvfile,1) > 0) then
       huv => safemem_realloc(huv,grid%totalLocalPointsK,solvent%numAtomTypes,.false.)
       if (solvent%numAtomTypes /= ubound(huvfile,1)) &
            call rism_report_error('(a,i4,a,i4)',"Different number of --huv files (",ubound(huvfile,1)&
            ,"than solvent sites", solvent%numAtomTypes)
       do i = 1, solvent%numAtomTypes
          call readDXHeader(huvfile(i),dxOrigin,delta,tmppos)
          if (sum(abs(tmppos- npos))/=0) call rism_report_error(trim(huvfile(i))//" is the wrong size")
          call readDX(huvfile(i),huv(:,i),grid%localDimsR,dxOrigin,delta)
       end do
    end if

    !Cuv
    if (ubound(cuvfile,1) > 0) then
       cuv => safemem_realloc(cuv,grid%localDimsR(1),grid%localDimsR(2),grid%localDimsR(3),&
            solvent%numAtomTypes,.false.)
       if (solvent%numAtomTypes /= ubound(cuvfile,1)) &
            call rism_report_error('(a,i4,a,i4)',"Different number of --cuv files (",ubound(guvfile,1)&
            ,"than solvent sites", solvent%numAtomTypes)
       do i = 1, solvent%numAtomTypes
          call readDXHeader(cuvfile(i),dxOrigin,delta,tmppos)
          if (sum(abs(tmppos- npos))/=0) call rism_report_error(trim(cuvfile(i))//" is the wrong size")
          call readDX(cuvfile(i),cuv(:,:,:,i),dxOrigin,delta)
       end do
    end if

    !hLongRangeAsymR and dcfLongRangeAsympR
    if (len_trim(ahrfile) /= 0) then
       tmp1 => readDX_p(ahrfile,dxOrigin,delta,tmppos)
       if (sum(abs(tmppos- npos))/=0) call rism_report_error(trim(ahrfile)//" is the wrong size")
       call rism3d_potential_setAsympTCF(potential,tmp1)
       longRange = .true.
    end if
    if (len_trim(acrfile) /= 0) then
       tmp1 => readDX_p(acrfile,dxOrigin,delta,tmppos)
       if (sum(abs(tmppos- npos))/=0) call rism_report_error(trim(acrfile)//" is the wrong size")
       call rism3d_potential_setAsympDCF(potential,tmp1)
       longRange = .true.
    end if
    call sanityCheck2()
  end subroutine readFiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculates requested thermodynamic quanities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc()
    use rism_util, only : translate
    implicit none
    _REAL_ :: cm(3)
    integer :: i
    !start with long range asymptotics if we have a solute object
    if (solute%numAtoms /= 0) then
       call translate(solute%position,solute%numAtoms,-dxOrigin)
       if (.not. associated(potential%tcfLongRangeAsympR) .or. .not. associated(potential%dcfLongRangeAsympR))&
            call rism3d_potential_dcf_tcf_long_range_asymptotics(potential)
       longRange = .true.
    end if
    if (partialMolarVolume) then
       partialMolarVolumeRes => safemem_realloc(partialMolarVolumeRes, ubound(volbox,2))
       do i=1,ubound(volbox,2)
          partialMolarVolumeRes(i)  = rism3d_closure_partialMolarVolume(closure,cuv)
       end do
    end if
    if (number) then
       numRes => safemem_realloc(numRes,solvent%numAtomTypes, ubound(volbox,2))
       numResLR => safemem_realloc(numResLR,solvent%numAtomTypes, ubound(volbox,2))
       do i=1,ubound(volbox,2)
          numRes(:,i) = rism3d_closure_excessParticles(closure,guv)
          if (longRange)&
               numResLR(:,i) = rism3d_closure_aexcessParticles(closure,guv)
       end do
    end if
    if (excessChemicalPotential) then
       excessChemicalPotentialRes => safemem_realloc(excessChemicalPotentialRes,solvent%numAtomTypes, ubound(volbox,2))
       excessChemicalPotentialResLR => safemem_realloc(excessChemicalPotentialResLR,solvent%numAtomTypes, ubound(volbox,2))
       do i=1,ubound(volbox,2)
          excessChemicalPotentialRes(:,i) = rism3d_closure_excessChemicalPotential(closure,huv,cuv)
          if (longRange)&
               excessChemicalPotentialResLR(:,i) = rism3d_closure_aexcessChemicalPotential(closure,huv,cuv)
       end do
    end if
    if (force) then
       write(0,'(a)') "ERROR: --force calculations not supported at this time"
       stop
    end if

  end subroutine calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!print results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeOutput()
    use constants, only : KB ,COULOMB_CONST_E
    implicit none
    integer :: ivol, iatom, idim
    integer :: unit
    !chCon :: convert internal charge units [sqrt(kT A)]
    _REAL_ :: chCon, vol
    character(len=16) :: name_fmt="(a15)", unit_fmt="(a10)", value_fmt="(1p,1x,e15.8)"
    chCon= sqrt((KB *solvent%temperature)/COULOMB_CONST_E)

    unit = rism_report_getMUnit()

    write(unit,*)

    if (solute%numAtoms/=0) then
       write(unit,'(1p,a,e15.8,a)') "Solute charge ",sum(solute%charge)*chCon, " e"
    end if


    if (number) then
       write(unit,*)
       write(unit,'(a)') "Number of atoms in"
       do ivol=1,ubound(numRes,2)
          call writeVolBox(ivol)
          vol = volume(ivol)

          write(unit,name_fmt,advance='no') ""
          write(unit,unit_fmt,advance='no') "Units"
          do iatom=1,solvent%numAtomTypes
             write(unit,name_fmt,advance='no') solvent%atomName(iatom)
          end do
          write(unit,name_fmt) "Total"
          write(unit,name_fmt,advance='no') "Bulk"
          write(unit,unit_fmt,advance='no') "#"
          do iatom=1,solvent%numAtomTypes
             write(unit,value_fmt,advance='no') solvent%density(iatom)*vol
          end do
          write(unit,value_fmt) sum(solvent%density(:)*vol)
          write(unit,name_fmt,advance='no') "Number"
          write(unit,unit_fmt,advance='no') "#"
          do iatom=1,solvent%numAtomTypes
             write(unit,value_fmt,advance='no') numRes(iatom,ivol)+solvent%density(iatom)*vol
          end do
          write(unit,value_fmt) sum(numRes(:,ivol)+solvent%density(:)*vol)
          write(unit,name_fmt,advance='no') "Excess"
          write(unit,unit_fmt,advance='no') "#"
          do iatom=1,solvent%numAtomTypes
             write(unit,value_fmt,advance='no') numRes(iatom,ivol)
          end do
          write(unit,value_fmt) sum(numRes(:,ivol))
          if (longRange) then
             write(unit,name_fmt,advance='no') "ExcessLR"
             write(unit,unit_fmt,advance='no') "#"
             do iatom=1,solvent%numAtomTypes
                write(unit,value_fmt,advance='no') numResLR(iatom,ivol)
             end do
             write(unit,value_fmt) sum(numResLR(:,ivol))
          end if
          write(unit,name_fmt,advance='no') "Charge"
          write(unit,unit_fmt,advance='no') "e"
          do iatom=1,solvent%numAtomTypes
             write(unit,value_fmt,advance='no') numRes(iatom,ivol)*solvent%charge(iatom)*chCon
          end do
          write(unit,value_fmt) sum(numRes(:,ivol)*solvent%charge(:)*chCon)
          if (longRange) then
             write(unit,name_fmt,advance='no') "ChargeLR"
             write(unit,unit_fmt,advance='no') "e"
             do iatom=1,solvent%numAtomTypes
                write(unit,value_fmt,advance='no') numResLR(iatom,ivol)*solvent%charge(iatom)*chCon
             end do
             write(unit,value_fmt) sum(numResLR(:,ivol)*solvent%charge(:)*chCon)
          end if
       end do
    end if
    if (excessChemicalPotential) then
       write(unit,*)
       write(unit,'(a)') "Excess chemical potential from"
       do ivol=1,ubound(numRes,2)
          call writeVolBox(ivol)
          vol = volume(ivol)

          write(unit,name_fmt,advance='no') ""
          write(unit,unit_fmt,advance='no') "Units"
          do iatom=1,solvent%numAtomTypes
             write(unit,name_fmt,advance='no') solvent%atomName(iatom)
          end do
          write(unit,name_fmt) "Total"
          write(unit,name_fmt,advance='no') "ExcessChemicalPotential"
          write(unit,unit_fmt,advance='no') "kcal/mol"
          do iatom=1,solvent%numAtomTypes
             write(unit,value_fmt,advance='no') excessChemicalPotentialRes(iatom,ivol)*KB*solvent%temperature
          end do
          write(unit,value_fmt) sum(excessChemicalPotentialRes(:,ivol))*KB*solvent%temperature
          if (longRange) then
             write(unit,name_fmt,advance='no') "ExcessChemicalPotentialLR"
             write(unit,unit_fmt,advance='no') "kcal/mol"
             do iatom=1,solvent%numAtomTypes
                write(unit,value_fmt,advance='no') excessChemicalPotentialResLR(iatom,ivol)*KB*solvent%temperature
             end do
             write(unit,value_fmt) sum(excessChemicalPotentialResLR(:,ivol))*KB*solvent%temperature
          end if
       end do
    end if
    if (partialMolarVolume) then
       write(unit,*)
       write(unit,'(1p,a,e24.16,a)') "Partial molar volume ",partialMolarVolumeRes(1)," A^3/mol"
    end if
    if (force) then
       write(0,'(a)') "ERROR: --force calculations not supported at this time"
       stop
    end if

    !combined results
    if (partialMolarVolume .and. number) then
       write(unit,*)
       write(unit,'(a)') "PARTIALMOLARVOLUME 'corrected' Number of atoms in"
       do ivol=1,ubound(numRes,2)
          call writeVolBox(ivol)
          vol = volume(ivol)

          write(unit,name_fmt,advance='no') ""
          write(unit,unit_fmt,advance='no') "Units"
          do iatom=1,solvent%numAtomTypes
             write(unit,name_fmt,advance='no') solvent%atomName(iatom)
          end do
          write(unit,name_fmt) "Total"
          write(unit,name_fmt,advance='no') "Number(PARTIALMOLARVOLUME)"
          write(unit,unit_fmt,advance='no') "#"
          do iatom=1,solvent%numAtomTypes
             write(unit,value_fmt,advance='no') numRes(iatom,ivol)+solvent%density(iatom)*(vol +partialMolarVolumeRes(1))
          end do
          write(unit,value_fmt) sum(numRes(:,ivol)+solvent%density(:)*(vol +partialMolarVolumeRes(1)))
          write(unit,name_fmt,advance='no') "Excess(PARTIALMOLARVOLUME)"
          write(unit,unit_fmt,advance='no') "#"
          do iatom=1,solvent%numAtomTypes
             write(unit,value_fmt,advance='no') numRes(iatom,ivol)&
                  +solvent%density(iatom)*partialMolarVolumeRes(1)
          end do
          write(unit,value_fmt) sum( numRes(:,ivol)&
               +solvent%density(:)*partialMolarVolumeRes(1))
          if (longRange) then
             write(unit,name_fmt,advance='no') "ExcessLR(PARTIALMOLARVOLUME)"
             write(unit,unit_fmt,advance='no') "#"
             do iatom=1,solvent%numAtomTypes
                write(unit,value_fmt,advance='no') numResLR(iatom,ivol)&
                     +solvent%density(iatom)*partialMolarVolumeRes(1)
             end do
             write(unit,value_fmt) sum(numResLR(:,ivol)&
                  +solvent%density(:)*partialMolarVolumeRes(1))
          end if
          write(unit,name_fmt,advance='no') "Charge(PARTIALMOLARVOLUME)"
          write(unit,unit_fmt,advance='no') "e"
          do iatom=1,solvent%numAtomTypes
             write(unit,value_fmt,advance='no') solvent%charge(iatom)*chCon&
                  *(numRes(iatom,ivol)+solvent%density(iatom)*partialMolarVolumeRes(1))
          end do
          write(unit,value_fmt) sum(solvent%charge(:)*chCon&
               *(numRes(:,ivol)+solvent%density(:)*partialMolarVolumeRes(1)))
          if (longRange) then
             write(unit,name_fmt,advance='no') "ChargeLR(PARTIALMOLARVOLUME)"
             write(unit,unit_fmt,advance='no') "e"
             do iatom=1,solvent%numAtomTypes
                write(unit,value_fmt,advance='no') solvent%charge(iatom)*chCon&
                     *(numResLR(iatom,ivol)+solvent%density(iatom)*partialMolarVolumeRes(1))
             end do
             write(unit,value_fmt) sum(solvent%charge(:)*chCon&
                  *(numResLR(:,ivol)+solvent%density(:)*partialMolarVolumeRes(1)))
          end if
       end do
    end if
  end subroutine writeOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Writes out description of the ith volume box
!!!IN:
!!!   ivol : volume box index
!!!SIDE-EFFECTS:
!!!   write box info to message unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeVolBox(ivol)
    implicit none
    integer, intent(in) :: ivol
    integer :: idim, unit
    _REAL_ :: vol
    unit = rism_report_getMUnit()
    write(unit,'(a,i4)') "Volume ",ivol
    do idim = 1, 3
       write(unit,'(a,a8,2(a1,a8),a)',advance='no') "(","start",":","end",":","stride",")"
       if (idim == 3) then
          write(unit,'(a)') "  "
       else
          write(unit,'(a)', advance='no') "   X "
       end if
    end do
    do idim = 1, 3
       write(unit,'(a,f8.3,2(a1,f8.3),a)', advance='no') &
            "(",&
            (volbox(idim,ivol)%start -1) * grid%spacing(idim) + dxOrigin(idim),":",&
            (volbox(idim,ivol)%stop  -1) * grid%spacing(idim) + dxOrigin(idim),":",&
            volbox(idim,ivol)%stride    * grid%spacing(idim) ,&
            ")"
       if (idim == 3) then
          write(unit,'(a)') " A"
       else
          write(unit,'(a)', advance='no') " A X "
       end if
    end do
  end subroutine writeVolBox

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the volume of the ith volume box
!!!IN:
!!!   ivol : volume box index
!!!OUT:
!!!   volume [A]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function volume(ivol)
    implicit none
    integer, intent(in) :: ivol
    integer :: idim
    _REAL_ :: volume
    volume = 1d0
    do idim=1,3
       volume=volume*(volbox(idim,ivol)%stop-volbox(idim,ivol)%start + 1)*grid%spacing(idim)
    end do
  end function volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Clean up memory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cleanup()
    implicit none
    call rism3d_solvent_destroy(solvent)
    call rism3d_solute_destroy(solute)
    call rism3d_closure_destroy(closure)
    call rism3d_potential_destroy(potential)
    call rism3d_grid_destroy(grid)
    if (safemem_dealloc(numRes)/=0) call rism_report_error("Failed to deallocate numRes")
    if (safemem_dealloc(numResLR)/=0) call rism_report_error("Failed to deallocate numResLR")
    if (safemem_dealloc(partialMolarVolumeRes)/=0) call rism_report_error("Failed to deallocate partialMolarVolumeRes")
    if (safemem_dealloc(guvfile)/=0) call rism_report_error("Failed to deallocate guvfile")
    if (safemem_dealloc(huvfile)/=0) call rism_report_error("Failed to deallocate huvfile")
    if (safemem_dealloc(cuvfile)/=0) call rism_report_error("Failed to deallocate cuvfile")
    if (safemem_dealloc(volspec)/=0) call rism_report_error("Failed to deallocate volspec")
    if (associated(volbox)) then
       deallocate(volbox)
    end if
  end subroutine cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Usage description
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine usage()
    implicit none
    integer :: munit
    munit = rism_report_getMUnit()
    call rism_report_setMUnit(rism_report_getEUnit())
    call rism_report_message("USAGE:   thermo [--all | [--exchem] [--number] [--pmv] [--force]]")
    !    call rism_report_message("                [--vol vol1 [vol2 [...]]]")
    call rism_report_message("                (-x|--xvv) xvv [(-g|--guv) guvs | (-h|--huv) huvs] [(-c|--cuv) cuvs]")
    call rism_report_message("                [--parm parm --crd pdb|rst7] [--ahr ahr --acr acr]")
    call rism_report_message("")
    call rism_report_message("Select the thermodynamic quantities of interest and the vol over")
    call rism_report_message("which to calculate them. The default is all space with long range ")
    call rism_report_message("corrections.")
    call rism_report_message("")
    call rism_report_message("For each thermodynamic quantity specific input files are required.")
    call rism_report_message("Xvv is required in all cases and Guv and Huv files are interchangable")
    call rism_report_message("as long as they are correclty labelled.")
    call rism_report_message("For long range corrections, solute parameter and coordinate files are")
    call rism_report_message("required.")
    call rism_report_message("--exchem (excess chemical potential)    :           " // &
         "huv|guv, cuv [, asymhr, asymcr]")
    call rism_report_message("                                                    [, parm, crd]")
    call rism_report_message("--number (number of solvent atoms)      :           " // &
         "huv|guv [, parm, crd] [, cuv] [, asymhr, asymcr]")
    call rism_report_message("--pmv (partial molar volume)            :           cuv")
    call rism_report_message("--force (mean solvation force per atom) :           huv|guv, parm, crd")
    call rism_report_setMUnit(munit)
  end subroutine usage
end module thermo_m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program thermo
  use thermo_m
  implicit none

  call getOptions()
  call readFiles()
  call calc()
  call writeOutput()
  call cleanup()

end program thermo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Timer stubs required by 3D-RISM code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine timer_start( label )
  integer label
end subroutine timer_start

subroutine timer_stop( label )
  integer label
end subroutine timer_stop
