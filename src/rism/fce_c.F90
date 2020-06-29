!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2014 by Andriy Kovalenko,
!Tyler Luchko, Igor Omelyan, and David A. Case.
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
!
!4) I. Omelyan and A. Kovalenko, J. Chem. Phys., 139, 244106 (2013).

#include "../include/dprec.fh"

!> Object for force-coordinate-extrapolation (FCE) multiple time step (MTS).
!! 
!! N total previous frames are held in memory for both position and
!! solvation force.  When we are extrapolating the force for a given
!! time step we first express the current position as a linear
!! combination of the previous steps:
!! 
!!   R^(N+1) ~ \sum^N_i a_i R^i
!! 
!! the a_i that best satisfy this equation are then used to calculate
!! an approximation of the solvation:
!! 
!!   F^(N+1) ~ \sum^N_i a_i F^i
module fce_c
  implicit none
  type fce
     sequence
     !enumerated values for coordinate type to use as the basis for extrapolation
     integer :: CRDPOS=0, CRDDIST=1, CRDXYZ=2

     !nbasis  :: total number of basic coordinate-force pairs used for interpolation
     !nbase   :: number of leading FCE basis vectors

     ! Note that if nbasis > nbase, then selecting of the best nbase points
     ! among all nbasis coordinates is carried out

     !crd     :: 0-position, 1-distance, 2-xyz distance

     !weigh   ::  0-usual (default) or /=0-weighted coordinate minimization
     !            deviations [expensive but more precise]

     !sort    ::  0-no sorting
     !        ::/=0-permorm sorting even through nbasis=nbase or additionally
     !              to selecting if nbasis>nbase 

     ! It should be pointed out that contrary to selecting, the sorting
     ! is not so important because it should lead to the same results as
     ! for no sorting, provided the round-off errors are neglected

     ! Note that no weighting (weigh=0) are performed for trans=0 but
     ! selecting (if nbasis>nbase) and sorting (if sort/=0) are possible

     !ifreq   :: the extended to best basic mapping list will be
     !           updated every ifreq impulsive step (ifreq >= 1)

     !ntfrcor :: 0-no net force correction
     !        ::/=0-perform net force correction
     ! Note that the net force is not equal to zero due to the approximate
     ! character of the individual extrapolation of force acting on each atom
      
     !enormsw ::  0-no minimization of the norms of the solutions
     !        ::/=0-minimization with enormsw-weight

     !trans ::
     !
     ! 0-no coordinate transformation at force extrapolation (very fast
     !   but not precise, can be used only for small outer time steps
     !   which do not exceed 200 fs) [default, not recommended for larger
     !   spacing and solvent molecules greater than 10 A in size]
     !
     ! 1-transformation with respect to the first basic point, while
     !   selecting according to the position of the current coordinate
     !   [recommended for large steps up to order of several picoseconds,
     !    because fast and precise, but for small macromolecules only]
     !
     ! 2-transformation with respect to the first basic point, while
     !   selecting according to the position of the current coordinate,
     !   i.e. like trans=1, but now using the normal equations method
     !   for the linear least squares coordinate deviation minimization
     !   instead of the QR decomposition approach as in the cases
     !   trans = 0 and 1. Moreover, the normal equations method is
     !   complemeted here by the e-minimization of the norm of the
     !   solutions and the ifreq-scheme to accelerate the calculations.
     !   Note that if enormsw=0 and ifreq=1, then trans=2 is equivalent
     !   to trans=1. However, an extra precision can be reached when
     !   enormsw accepts small positive values and a significant speedup
     !   can be observed if ifreq >> 1, especially at a large number of
     !   fcenbase (several hundreds). The case trans=2 is recommended
     !   as the best choice for relatively small macromolecules, since
     !   it gives the possibility to apply huge outer steps (up to ten
     !   picoseconds).
     
     ! 3-Same as 2 above (place holder / buffer)

     ! 4-no normalization, transformation, selecting, sorting, and 
     !   balancing, but with individual force extrapolation and 
     !   possible cutting-off of the neighbours relatively to each
     !   current atom [default for the original AMBER11 version]

     ! 5-individual transformation and selecting with respect to the
     !   current coordinate of each atom using a neighbouring scheme
     !   compemented by the e-minimization and ifreq-scheme as well
     !   as all other developed techniques. It is recommended for large
     !   macromolecules of greater than 10 A in size and can be used
     !   with huge outer steps (up to order of several picoseconds).
     !
     ! 6-individual transformation and selecting with respect to the
     !   post coordinate of each atom using a neighbouring scheme
     !   compemented by the e-minimization and the full ifreq-support.
     !   It is recommended for large macromolecules and can be used
     !   with huge outer steps (up to order of several picoseconds).
     !   It appears to be better than the above case trans=5 (partial
     !   ifreq-support version) because can be exploited with larger
     !   number (up to N~100-200) of basic points providing a higher
     !   accuracy (with nearly the same computational efforts as the
     !   trans=5-version at N~30), but may require more memory. Note
     !   that at any values of ifreq, both the approaches have the
     !   same scheme for building the index mask which maps the
     !   extended set to the best subset and differ in the way of 
     !   constructing the transformation matrix. At ifreq=1, these
     !   two approaches are equivalent.
     !
     !

     !nsample :: number of samples collected
     !natom   :: number of atoms
     integer :: nbasis=-1, crd=-1, nsample=0, natom=-1

     !mpirank - MPI rank
     !mpicomm - MPI mpicomm
     !mpisize - number of MPI processes
     !atom0 - first atom in MPI decomposition
     !atomF - final atom in MPI decomposition
     integer :: mpirank=0, mpicomm=0, mpisize=1
     integer :: atom0, atomf

     !cut :: distance cutoff used for creating the basis set
     _REAL_ :: cut

     !cutlist :: for each atoms, list of atoms to include in force extrapolation (natom+1,natom).
     !           the first row element indicates how many atoms follow in this row.
     integer,pointer :: cutlist(:,:) => NULL()
     !refit :: TBA
     logical,pointer :: refit(:) => NULL()

     !force :: previous forces (dimensions:natom:nbasis)
     !coord :: previous coordinates (dimensions:natom:nbasis)
     _REAL_,pointer :: force(:,:,:) => NULL(), coord(:,:,:) => NULL()

     !sforce :: transformed forces (dimensions:natom:nbasis)
     !scoord :: transformed coordinates (dimensions:natom:nbasis)
     _REAL_,pointer :: sforce(:,:,:) => NULL(), scoord(:,:,:) => NULL()

     ! Additional pointers which are necessary for trans=2,3,5, and 6

     _REAL_,pointer :: Ad(:,:) => NULL()
     integer,pointer :: irdd(:,:) => NULL()
     integer,pointer :: ipvs(:) => NULL()

     _REAL_,pointer :: Ads(:,:,:) => NULL()
     integer,pointer :: ipvss(:,:) => NULL()

     integer,pointer :: tscuts(:,:) => NULL()
     _REAL_,pointer :: coordps(:,:,:) => NULL()
     _REAL_,pointer :: srwrwr2d(:,:) => NULL()

     _REAL_,pointer :: sfis(:,:,:,:) => NULL()

     !! coeff :: the a_i coefficients in the above equation
     _REAL_, pointer :: coeff(:,:) => NULL()

     integer :: nbase=-1, weigh=-1, trans=-1, sort=-1, ifreq=-1, ntfrcor=-1
     integer :: bwrtunit=1117, cstep=0

     _REAL_  :: enormsw=-1.d0

!    integer :: paddings

  end type fce

  interface fce_new
     module procedure fce_new_std
  end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!public subroutines and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
public fce_new, fce_destroy, fce_update, fce_forcea, fce_forceb, fce_forcebm, &
       fce_forcesa, fce_forcesan, fce_force, fce_transformi, fce_estimate, &
       fce_wrtbasis, fce_readbasis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private subroutines and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
private nlist, selects, hsort

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!public subroutines and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Constructor. If the this is MPI, then only the rank 0 process parameters will
!!!be used.
!!!IN:
!!!   this :: FCE object
!!!   natom :: number of atoms
!!!   nbasis :: number of all basis vectors
!!!   nbase  :: number of leading basis vectors
!!!   weigh  :: weighted minimization flag
!!!   trans  :: transforming coordinate flag
!!!   sort   :: sorting flag
!!!   ifreq  :: extended to basic maping list updating frequency
!!!   ntfrcor:: net force correction flag (for individual extrapolation)
!!!   enormsw:: solution norm minimization weight
!!!   crd :: coordinate orientation method. 0-position, 1-distance, 2-xyz distance
!!!   cut :: cut off distance for dependence
!!!   o_mpicomm :: (optional) MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_new_std(this,natom,nbasis,nbase,crd,weigh,trans,sort,ifreq, &
                         ntfrcor,enormsw,cut,o_mpicomm)
    use safemem
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    type(fce),intent(inout) :: this
    integer, intent(in) :: natom, nbasis, crd
    integer, intent(in) :: nbase, weigh, trans, sort, ifreq, ntfrcor
    _REAL_, intent(in) :: cut, enormsw
    integer, optional, intent(in) :: o_mpicomm
    integer :: err
#ifdef MPI
    this%mpicomm = 0
    if(present(o_mpicomm)) this%mpicomm = o_mpicomm
    if(this%mpicomm == MPI_COMM_NULL)&
       call rism_report_error("FCE: received NULL MPI communicator")
    call mpi_comm_rank(this%mpicomm,this%mpirank,err)
    if(err /=0) call rism_report_error&
         ("(a,i8)","FCE: could not get MPI rank for communicator ",this%mpicomm)
    call mpi_comm_size(this%mpicomm,this%mpisize,err)
    if(err /=0) call rism_report_error&
         ("(a,i8)","FCE: could not get MPI size for communicator ",this%mpicomm)
#endif /*MPI*/
    if(this%mpirank==0)then
       this%natom = natom
       this%nbasis = nbasis

       this%nbase = nbase
       this%weigh = weigh
       this%trans = trans
       this%sort  = sort
       this%ifreq = ifreq
       this%ntfrcor = ntfrcor
       this%enormsw = enormsw

       this%crd = crd
       this%cut = cut**2

       this%cutlist => safemem_realloc(this%cutlist,natom,natom,.false.)
       this%force => safemem_realloc(this%force,3,natom,nbasis,.false.)
       this%coord => safemem_realloc(this%coord,3,natom,nbasis,.false.)

       this%sforce => safemem_realloc(this%sforce,3,natom,nbasis,.false.)
       this%scoord => safemem_realloc(this%scoord,3,natom,nbasis,.false.)

       this%Ad => safemem_realloc(this%Ad,nbase+1,nbase+1,.false.)
       this%irdd => safemem_realloc(this%irdd,natom,nbase,.false.)
       this%ipvs => safemem_realloc(this%ipvs,nbase+1,.false.)

       this%coeff => safemem_realloc(this%coeff,3,natom,.false.)
    end if

#ifdef MPI
    !first distribute the pieces of information needed to allocate memory
    call mpi_bcast(this%natom,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast NATOM")
    call mpi_bcast(this%nbasis,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast NBASIS")

    call mpi_bcast(this%nbase,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast NBASE")

    call mpi_bcast(this%weigh,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast WEIGH")

    call mpi_bcast(this%trans,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast TRANS")

    call mpi_bcast(this%sort,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast SORT")

    call mpi_bcast(this%ifreq,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast IFREQ")

    call mpi_bcast(this%ntfrcor,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast NTFRCOR")

    call mpi_bcast(this%enormsw,1,mpi_double_precision,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast ENORMSW")

    call mpi_bcast(this%crd,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast CRD")
    call mpi_bcast(this%cut,1,mpi_double_precision,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast CUT")

    !non-master processes should now allocate memory
    if(this%mpirank /= 0) then
       this%cutlist => safemem_realloc(this%cutlist,this%natom,this%natom,.false.)
       this%force => safemem_realloc(this%force,3,this%natom,this%nbasis,.false.)
       this%coord => safemem_realloc(this%coord,3,this%natom,this%nbasis,.false.)

       this%sforce => safemem_realloc(this%sforce,3,this%natom,this%nbasis,.false.)
       this%scoord => safemem_realloc(this%scoord,3,this%natom,this%nbasis,.false.)

       this%Ad => safemem_realloc(this%Ad,this%nbase+1,this%nbase+1,.false.)
       this%irdd => safemem_realloc(this%irdd,this%natom,this%nbase,.false.)
       this%ipvs => safemem_realloc(this%ipvs,this%nbase+1,.false.)

       this%coeff => safemem_realloc(this%coeff,3,this%natom,.false.)
    end if

    !Arrays contain no data yet so there is nothing to transfer

#endif /*MPI*/

    !set local atom range for this process
!!  this%atom0 = this%natom/this%mpisize*this%mpirank+1
!!  this%atomF = min(int(this%natom*dble(this%mpirank+1)/dble(this%mpisize)),this%natom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A new scheme for distribution of all processors among atoms
! (this avoids dublications inherent in the previous approach)

    this%atom0 = idnint(dble(this%natom)/this%mpisize)*this%mpirank+1

    if(this%mpirank+1.lt.this%mpisize) then

    this%atomF = idnint(dble(this%natom)/this%mpisize)*(this%mpirank+1)

    else

    this%atomF = this%natom

    end if

! Make the memory used as small as possible

    if(trans.eq.6) then        ! If GSFE

    this%Ads => safemem_realloc(this%Ads,this%atomF-this%atom0+1,this%nbase+1,this%nbase+1,.false.)
    this%ipvss => safemem_realloc(this%ipvss,this%atomF-this%atom0+1,this%nbase+1,.false.)

    this%tscuts => safemem_realloc(this%tscuts,this%atomF-this%atom0+1,this%natom)
    this%coordps => safemem_realloc(this%coordps,this%atomF-this%atom0+1,3,this%natom)
    this%srwrwr2d => safemem_realloc(this%srwrwr2d,this%atomF-this%atom0+1,this%natom)

    this%sfis => safemem_realloc(this%sfis,3,3,this%atomF-this%atom0+1,this%nbase,.false.)

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine fce_new_std

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Destroyer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_destroy(this)
    use safemem
    implicit none
    type(fce),intent(inout) :: this
    integer :: err
    this%natom = -1
    this%nbasis = -1

    this%nbase = -1
    this%weigh = -1
    this%sort  = -1
    this%ifreq = -1
    this%ntfrcor = -1
    this%enormsw = -1.d0

    this%crd = -1

    err = safemem_dealloc(this%cutlist)
    err = safemem_dealloc(this%refit)
    err = safemem_dealloc(this%force)
    err = safemem_dealloc(this%coord)
    err = safemem_dealloc(this%coeff)

    err = safemem_dealloc(this%sforce)
    err = safemem_dealloc(this%scoord)

    err = safemem_dealloc(this%Ad)
    err = safemem_dealloc(this%irdd)
    err = safemem_dealloc(this%ipvs)

    if(this%trans.eq.6) then

    err = safemem_dealloc(this%Ads)
    err = safemem_dealloc(this%ipvss)

    err = safemem_dealloc(this%tscuts)
    err = safemem_dealloc(this%coordps)
    err = safemem_dealloc(this%srwrwr2d)

    err = safemem_dealloc(this%sfis)

    end if

    this%trans = -1

  end subroutine fce_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!update the force and coordinate basis vectors
!!!IN:
!!!   this  : FCE object
!!!   force : forces on atoms.  In the case of MPI, it is assumed that the forces
!!!           are distributed across the processes and must be reduced internally
!!!   coord : coordinates of the atoms, assume to the the same for all processes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_update(this, force, coord)
    use safemem
    implicit none
#if defined(MPI)
      include 'mpif.h'
#endif /*defined(MPI)*/
    type(fce), intent(inout) :: this
    _REAL_, intent(in) :: force(3,this%natom), coord(3,this%natom)

    !workspace for the case of MPI 1.1
#ifdef MPI
    integer :: err
#  ifndef USE_MPI_IN_PLACE
    _REAL_, pointer :: tforce(:,:)=>NULL()
#  endif /*ifdef USE_MPI_IN_PLACE*/
#endif /*MPI*/

#ifdef RISM_DEBUG
    write(6,*) "FCE_UPDATE"; call flush(6)
    write(6,*) "FCE", this%atom0, this%atomF, this%mpisize, this%mpirank, this%mpicomm, this%cut
#endif /*RISM_DEBUG*/

    !
    !For now we will shift the data through the storage arrays as more data is
    !added.  Later, this should be modified to just update a pointer to the most
    !recent entry
    !
    if(this%nsample > 0)then
       this%coord(:,:,2:min(this%nsample+1,this%nbasis)) = this%coord(:,:,1:min(this%nsample, this%nbasis-1))
       this%force(:,:,2:min(this%nsample+1,this%nbasis)) = this%force(:,:,1:min(this%nsample, this%nbasis-1))
    end if
    this%coord(:,:,1) = coord
    this%force(:,:,1) = force
    !
    !reduce the atoms locally
    !
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE,this%force(:,:,1),3*this%natom,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not reduce force")
#  else /*ifdef USE_MPI_IN_PLACE*/
    tforce => safemem_realloc(tforce,3,this%natom,.false.)
    call mpi_allreduce(this%force(:,:,1),tforce,3*this%natom,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not reduce force")
    this%force(:,:,1) = tforce
    if(safemem_dealloc(tforce) /=0) call rism_report_error("FCE_UPDATE: deallocate TFORCE failed")
#  endif /*ifdef USE_MPI_IN_PLACE*/
#endif /*MPI*/

    this%nsample = min(this%nsample+1,this%nbasis)

  end subroutine fce_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Transform the force and coordinate basis vectors
!!!IN:
!!!   this  : FCE object
!!!   force : forces on atoms
!!!   coord : coordinates of the atoms
!!!OUT:
!!!  sforce : transformed forces
!!!  scoord : transformed coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_transformi(this)

    use safemem
    implicit none

#ifdef MPI
    include 'mpif.h'
#endif

    type(fce), intent(inout) :: this

    integer :: lswork, infos
    _REAL_, pointer :: swork(:)=>NULL()

    DOUBLE PRECISION hmu(4,4),ddd,q0,q1,q2,q3,fis(3,3),wjjj(4)
    DOUBLE PRECISION xiii1,xiii2,xiii3,fiii1,fiii2,fiii3

    integer :: iatm, isa

    integer :: err
#ifdef MPI
    DOUBLE PRECISION hmus(this%nbasis,4,4)
#endif

    ! The transformation is carried out with respect to the first pair

    ! The coordinates are known for all the processors
    do iatm = 1, this%natom
       this%scoord(1,iatm,1)=this%coord(1,iatm,1)
       this%scoord(2,iatm,1)=this%coord(2,iatm,1)
       this%scoord(3,iatm,1)=this%coord(3,iatm,1)
    end do

    ! The forces are distributed across the processors
    do iatm = this%atom0, this%atomF
       this%sforce(1,iatm,1)=this%force(1,iatm,1)
       this%sforce(2,iatm,1)=this%force(2,iatm,1)
       this%sforce(3,iatm,1)=this%force(3,iatm,1)
    end do

    ! Find optimal rotation of a solute as a whole for the next pairs to
    ! provide the minimum for the distances between the basic coordinates

#ifdef MPI

    ! Summing up the quaternion matrix over atoms across the processes

    hmus=0.d0
    do iatm = this%atom0, this%atomF
       fiii1=this%coord(1,iatm,1)
       fiii2=this%coord(2,iatm,1)
       fiii3=this%coord(3,iatm,1)
       do isa=2,this%nsample
          xiii1=this%coord(1,iatm,isa)
          xiii2=this%coord(2,iatm,isa)
          xiii3=this%coord(3,iatm,isa)
          ddd=(xiii1+fiii1)**2+(xiii2+fiii2)**2+(xiii3+fiii3)**2
          hmus(isa,1,1)=hmus(isa,1,1)+((xiii1-fiii1)**2+(xiii2-fiii2)**2+ &
               (xiii3-fiii3)**2)
          hmus(isa,1,2)=hmus(isa,1,2)+2.d0*(xiii2*fiii3-xiii3*fiii2)
          hmus(isa,1,3)=hmus(isa,1,3)+2.d0*(xiii3*fiii1-xiii1*fiii3)
          hmus(isa,1,4)=hmus(isa,1,4)+2.d0*(xiii1*fiii2-xiii2*fiii1)
          hmus(isa,2,2)=hmus(isa,2,2)+(ddd-2.d0*(xiii1*fiii1+fiii1*xiii1))
          hmus(isa,2,3)=hmus(isa,2,3)     -2.d0*(xiii1*fiii2+fiii1*xiii2)
          hmus(isa,2,4)=hmus(isa,2,4)     -2.d0*(xiii1*fiii3+fiii1*xiii3)
          hmus(isa,3,3)=hmus(isa,3,3)+(ddd-2.d0*(xiii2*fiii2+fiii2*xiii2))
          hmus(isa,3,4)=hmus(isa,3,4)     -2.d0*(xiii2*fiii3+fiii2*xiii3)
          hmus(isa,4,4)=hmus(isa,4,4)+(ddd-2.d0*(xiii3*fiii3+fiii3*xiii3))
       end do
    end do

    ! Collecting the partial sums from different processes and placing the
    ! result into the same quaternion matrix known for all processors

    call mpi_allreduce(MPI_IN_PLACE,hmus,this%nbasis*16,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)

    do isa=2,this%nsample

       hmu(:,:)=hmus(isa,:,:)

#else

       ! Summing up the quaternion matrix over atoms in the single process mode

       do isa=2,this%nsample
          hmu=0.d0
          do iatm = 1, this%natom
             fiii1=this%coord(1,iatm,1)
             fiii2=this%coord(2,iatm,1)
             fiii3=this%coord(3,iatm,1)
             xiii1=this%coord(1,iatm,isa)
             xiii2=this%coord(2,iatm,isa)
             xiii3=this%coord(3,iatm,isa)
             ddd=(xiii1+fiii1)**2+(xiii2+fiii2)**2+(xiii3+fiii3)**2
             hmu(1,1)=hmu(1,1)+((xiii1-fiii1)**2+(xiii2-fiii2)**2+ &
                  (xiii3-fiii3)**2)
             hmu(1,2)=hmu(1,2)+2.d0*(xiii2*fiii3-xiii3*fiii2)
             hmu(1,3)=hmu(1,3)+2.d0*(xiii3*fiii1-xiii1*fiii3)
             hmu(1,4)=hmu(1,4)+2.d0*(xiii1*fiii2-xiii2*fiii1)
             hmu(2,2)=hmu(2,2)+(ddd-2.d0*(xiii1*fiii1+fiii1*xiii1))
             hmu(2,3)=hmu(2,3)     -2.d0*(xiii1*fiii2+fiii1*xiii2)
             hmu(2,4)=hmu(2,4)     -2.d0*(xiii1*fiii3+fiii1*xiii3)
             hmu(3,3)=hmu(3,3)+(ddd-2.d0*(xiii2*fiii2+fiii2*xiii2))
             hmu(3,4)=hmu(3,4)     -2.d0*(xiii2*fiii3+fiii2*xiii3)
             hmu(4,4)=hmu(4,4)+(ddd-2.d0*(xiii3*fiii3+fiii3*xiii3))
          end do

#endif

          ! The problem is reduced to find eigenvalues and eigenvectors
          ! of the corresponding quaternion symmetric matrices

          ! First query the optimal workspace

          lswork = -1
          swork => safemem_realloc(swork,1,.false.)
          CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
          lswork = INT(swork(1))
          swork => safemem_realloc(swork,lswork,.false.)

          ! Now solve the eigenproblem

          CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

          if(infos.ne.0) then
             print *,infos
             stop 'Message from DSYEV'
          end if

          ! Define the 4-dimensional quaternion corresponding to the smallest eigenvalue

          q0=hmu(1,1)
          q1=hmu(2,1)
          q2=hmu(3,1)
          q3=hmu(4,1)

          ! Determine the rotational matrix in 3-dimensional space

          fis(1,1)=q0**2+q1**2-q2**2-q3**2
          fis(1,2)=2.d0*(-q0*q3+q1*q2)
          fis(1,3)=2.d0*( q0*q2+q1*q3)
          fis(2,1)=2.d0*( q0*q3+q1*q2)
          fis(2,2)=q0**2+q2**2-q1**2-q3**2
          fis(2,3)=2.d0*(-q0*q1+q2*q3)
          fis(3,1)=2.d0*(-q0*q2+q1*q3)
          fis(3,2)=2.d0*( q0*q1+q2*q3)
          fis(3,3)=q0**2+q3**2-q1**2-q2**2

          ! Perform the coordinate and force transformations

          do iatm = 1, this%natom
             this%scoord(1,iatm,isa)=fis(1,1)*this%coord(1,iatm,isa)+fis(2,1)*this%coord(2,iatm,isa)+fis(3,1)*this%coord(3,iatm,isa)
             this%scoord(2,iatm,isa)=fis(1,2)*this%coord(1,iatm,isa)+fis(2,2)*this%coord(2,iatm,isa)+fis(3,2)*this%coord(3,iatm,isa)
             this%scoord(3,iatm,isa)=fis(1,3)*this%coord(1,iatm,isa)+fis(2,3)*this%coord(2,iatm,isa)+fis(3,3)*this%coord(3,iatm,isa)
          end do

          do iatm = this%atom0, this%atomF
             this%sforce(1,iatm,isa)=fis(1,1)*this%force(1,iatm,isa)+fis(2,1)*this%force(2,iatm,isa)+fis(3,1)*this%force(3,iatm,isa)
             this%sforce(2,iatm,isa)=fis(1,2)*this%force(1,iatm,isa)+fis(2,2)*this%force(2,iatm,isa)+fis(3,2)*this%force(3,iatm,isa)
             this%sforce(3,iatm,isa)=fis(1,3)*this%force(1,iatm,isa)+fis(2,3)*this%force(2,iatm,isa)+fis(3,3)*this%force(3,iatm,isa)
          end do

       end do
       err = safemem_dealloc(swork)
     end subroutine fce_transformi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Using previous forces and coordinates, predicts solvation forces
!!! for the current set of coordinates.  Like  linprojpredict, except
!!! we are finding the forces on each atom individually and rotate the
!!! solute for each to optimize the prediction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 0-no transformation [not recommended] and weighting (weigh=0) but
     ! with possible selecting (if nbasis>nbase) and sorting (if sort/=0) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_forcea(this,force,coord)

    use safemem
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    integer :: iatm,id,err

    !LAPACK variables
    integer :: M, N, NN, NRHS=1, LDA, LDB, LWORK, RANK, INFO
    integer, pointer :: JPVT(:)=>NULL()

    !uncertainties of input coordinates in A and B are assumed
    !to be negligible small
    _REAL_ :: RCOND=0.d0
    _REAL_, pointer :: A(:,:)=>NULL(), B(:,:)=>NULL(), work(:)=>NULL()
    _REAL_ :: rrr(this%nbasis),rar(this%nbase)

    integer :: irr(this%nbasis),ird(this%nbase),irs(this%nbase)

    integer :: isa,jsa,is1,iisap1

     N = this%nsample

! No selecting and sorting if nbasis=nbase and sort=0

          if(this%nbasis.eq.this%nbase.and.this%sort.eq.0) then

          do isa=1,this%nbase
          ird(isa)=isa
          end do

          end if

! Calculate the distances between the current solute position and its
! counterparts from the extended basis set if nbasis > nbase or sort/=0 

          if(this%nbasis.gt.this%nbase.or.this%sort.ne.0) then

          do isa=1,N
          rrr(isa)=0.d0
          do iatm = this%atom0, this%atomF
          rrr(isa)=rrr(isa)+(coord(1,iatm)-this%coord(1,iatm,isa))**2+ &
                            (coord(2,iatm)-this%coord(2,iatm,isa))**2+ &
                            (coord(3,iatm)-this%coord(3,iatm,isa))**2
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,rrr,N,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

! Normalizing the distances

          do isa=1,N
          rrr(isa)=dsqrt(rrr(isa)/this%natom)
          end do

          end if

! Select the first nbase minimal values among all N=nsample elements
! if N > nbase (note that nbase <= nsample <= nbasis). In other words,
! choose the best force-coordinates pairs.

          if(this%nbasis.gt.this%nbase.and.this%sort.ne.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) then

          call selects(this%nbase,N,rrr,irr)

          if(this%nbase.ne.1) then

          do isa=1,this%nbase-1
          rar(isa)=rrr(isa)
          end do

! Now perform sorting of the selected minimal values in the ascending order

          call hsort(this%nbase-1,rar,irs)

          end if

          irs(this%nbase)=this%nbase

          do isa=1,this%nbase
          ird(isa)=irr(irs(isa))
          end do

          end if

! On the beginning when N=nsample=nbase perform only sorting because
! then selecting is not necessary

          if(N.eq.this%nbase) then

          call hsort(N,rrr,ird)

          end if

          end if

! Only select the first nbase minimal values among all N=nsample elements
! if N > nbase without sorting (sort=0)

          if(this%nbasis.gt.this%nbase.and.this%sort.eq.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) call selects(this%nbase,N,rrr,irr)

          do isa=1,this%nbase
          ird(isa)=irr(isa)
          end do

          end if

! If sort/=0, perform sorting even for nbasis=nbase, i.e. when selecting
! is not necessary

          if(this%nbasis.eq.this%nbase.and.this%sort.ne.0) then

          call hsort(this%nbase,rrr,ird)

          end if

! The nearest pair

          is1=ird(1)

! Zero-order force approximation

          do iatm = this%atom0, this%atomF
          do id=1,3
          force(id,iatm) = this%force(id,iatm,is1)
          end do
          end do

          NN=this%nbase-1

          if(NN.gt.0) then

! High-order force extrapolation

          JPVT=>safemem_realloc(JPVT,NN,.false.)
          JPVT=0

! Solve the linear least-square problem

          M = 3 * this%natom
          LDA = M
          LDB = max(M,NN)

          A=>safemem_realloc(A,M,NN,.false.)
          B=>safemem_realloc(B,max(NN,M),1,.false.)

! Using multiprocessor technique to fill in the matrix elements

#ifdef MPI
          A=0.d0
          B=0.d0

          jsa=3*(this%atom0-1)
#else
          jsa=0
#endif

          do iatm = this%atom0, this%atomF
          B(jsa+1:jsa+3,1)=coord(:,iatm)-this%coord(:,iatm,is1)
          do isa=1,NN
          iisap1=ird(isa+1)
          A(jsa+1:jsa+3,isa)=this%coord(:,iatm,iisap1)-this%coord(:,iatm,is1)
          end do
          jsa=jsa+3
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,M*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    call mpi_allreduce(MPI_IN_PLACE,B,LDB,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DGELSY'
          end if

! Extrapolate the forces

          do iatm = this%atom0, this%atomF
          do isa=1,NN
          iisap1=ird(isa+1)
          do id=1,3
          force(id,iatm) = force(id,iatm)+&
          B(isa,1)*(this%force(id,iatm,iisap1)-this%force(id,iatm,is1))
          end do
          end do
          end do

          end if

    err = safemem_dealloc(A)
    err = safemem_dealloc(B)
    err = safemem_dealloc(WORK)
    err = safemem_dealloc(JPVT)

  end subroutine fce_forcea

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 1-transformation with respect to the first basic point, while
     !   selecting according to the position of the current coordinate
     !   with sorting (if sort/=0) (recommended, because fast and precise,
     !   but for relatively small solvent molecules)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_forceb(this,force,coord)

    use safemem
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    _REAL_ :: sforce(3,this%natom),scoord(3,this%natom)

    integer :: iatm,jatm,id,err

    !LAPACK variables
    integer :: M, N, NN, NRHS=1, LDA, LDB, LWORK, RANK, INFO
    integer, pointer :: JPVT(:)=>NULL()

    integer :: lswork, infos
   _REAL_, pointer :: swork(:)=>NULL()

    !uncertainties of input coordinates in A and B are assumed
    !to be negligible small
    _REAL_ :: RCOND=0.d0
    _REAL_, pointer :: A(:,:)=>NULL(), B(:,:)=>NULL(), work(:)=>NULL()
    _REAL_ :: rrr(this%nbasis),rar(this%nbase)

    integer :: irr(this%nbasis),ird(this%nbase),irs(this%nbase)

    integer :: isa,jsa,is1,iisap1

     DOUBLE PRECISION hmu(4,4),ddd,q0,q1,q2,q3,fis(3,3),wjjj(4),rweight
     DOUBLE PRECISION xiii1,xiii2,xiii3,fiii1,fiii2,fiii3,rwrwr2

#ifdef MPI
     integer inum,jsad(this%mpisize)
#endif

     N = this%nsample

! The transformation of the curent coordinate with respect to the first
! basic point is performed by minimization of the distance between them

       hmu=0.d0

! Summing up the quaternion matrix over atoms across the processes

       do iatm = this%atom0, this%atomF
       fiii1=this%coord(1,iatm,1)
       fiii2=this%coord(2,iatm,1)
       fiii3=this%coord(3,iatm,1)
       xiii1=coord(1,iatm)
       xiii2=coord(2,iatm)
       xiii3=coord(3,iatm)
       ddd=(xiii1+fiii1)**2+(xiii2+fiii2)**2+(xiii3+fiii3)**2
       hmu(1,1)=hmu(1,1)+((xiii1-fiii1)**2+(xiii2-fiii2)**2+ &
                          (xiii3-fiii3)**2)
       hmu(1,2)=hmu(1,2)+2.d0*(xiii2*fiii3-xiii3*fiii2)
       hmu(1,3)=hmu(1,3)+2.d0*(xiii3*fiii1-xiii1*fiii3)
       hmu(1,4)=hmu(1,4)+2.d0*(xiii1*fiii2-xiii2*fiii1)
       hmu(2,2)=hmu(2,2)+(ddd-2.d0*(xiii1*fiii1+fiii1*xiii1))
       hmu(2,3)=hmu(2,3)     -2.d0*(xiii1*fiii2+fiii1*xiii2)
       hmu(2,4)=hmu(2,4)     -2.d0*(xiii1*fiii3+fiii1*xiii3)
       hmu(3,3)=hmu(3,3)+(ddd-2.d0*(xiii2*fiii2+fiii2*xiii2))
       hmu(3,4)=hmu(3,4)     -2.d0*(xiii2*fiii3+fiii2*xiii3)
       hmu(4,4)=hmu(4,4)+(ddd-2.d0*(xiii3*fiii3+fiii3*xiii3))
       end do

#ifdef MPI

! Collecting the partial sums from different processes and placing the
! result into the same quaternion matrix known for all processors

    call mpi_allreduce(MPI_IN_PLACE,hmu,16,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)

#endif

! The problem is reduced to find eigenvalues and eigenvectors 
! of the corresponding quaternion symmetric matrix

! First query the optimal workspace

       lswork = -1
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! Define the 4-dimensional quaternion corresponding to the smallest eigenvalue

       q0=hmu(1,1)
       q1=hmu(2,1)
       q2=hmu(3,1)
       q3=hmu(4,1)

! Determine the rotational matrix in 3-dimensional space

       fis(1,1)=q0**2+q1**2-q2**2-q3**2
       fis(1,2)=2.d0*(-q0*q3+q1*q2)
       fis(1,3)=2.d0*( q0*q2+q1*q3)
       fis(2,1)=2.d0*( q0*q3+q1*q2)
       fis(2,2)=q0**2+q2**2-q1**2-q3**2
       fis(2,3)=2.d0*(-q0*q1+q2*q3)
       fis(3,1)=2.d0*(-q0*q2+q1*q3)
       fis(3,2)=2.d0*( q0*q1+q2*q3)
       fis(3,3)=q0**2+q3**2-q1**2-q2**2

! Perform the transformations

       do iatm = 1, this%natom
   scoord(1,iatm)=fis(1,1)*coord(1,iatm)+fis(2,1)*coord(2,iatm)+fis(3,1)*coord(3,iatm)
   scoord(2,iatm)=fis(1,2)*coord(1,iatm)+fis(2,2)*coord(2,iatm)+fis(3,2)*coord(3,iatm)
   scoord(3,iatm)=fis(1,3)*coord(1,iatm)+fis(2,3)*coord(2,iatm)+fis(3,3)*coord(3,iatm)
       end do

! No selecting and sorting if nbasis=nbase and sort=0

          if(this%nbasis.eq.this%nbase.and.this%sort.eq.0) then

          do isa=1,this%nbase
          ird(isa)=isa
          end do

          end if

! Calculate the distances between the current solute position and its
! counterparts from the extended basis set if nbasis > nbase or sort/=0

          if(this%nbasis.gt.this%nbase.or.this%sort.ne.0) then

          do isa=1,N
          rrr(isa)=0.d0
          do iatm = this%atom0, this%atomF
          rrr(isa)=rrr(isa)+(scoord(1,iatm)-this%scoord(1,iatm,isa))**2+ &
                            (scoord(2,iatm)-this%scoord(2,iatm,isa))**2+ &
                            (scoord(3,iatm)-this%scoord(3,iatm,isa))**2
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,rrr,N,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

! Normalizing the distances

          do isa=1,N
          rrr(isa)=dsqrt(rrr(isa)/this%natom)
          end do

          end if

! Select the first nbase minimal values among all N=nsample elements
! if N > nbase (note that nbase <= nsample <= nbasis). In other words,
! choose the best transformed force-coordinates pairs.

          if(this%nbasis.gt.this%nbase.and.this%sort.ne.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) then

          call selects(this%nbase,N,rrr,irr)

          if(this%nbase.ne.1) then

          do isa=1,this%nbase-1
          rar(isa)=rrr(isa)
          end do

! Now perform sorting of the selected minimal values in the ascending order

          call hsort(this%nbase-1,rar,irs)

          end if

          irs(this%nbase)=this%nbase

          do isa=1,this%nbase
          ird(isa)=irr(irs(isa))
          end do

          end if

! On the beginning when N=nsample=nbase perform only sorting because
! then selecting is not necessary

          if(N.eq.this%nbase) then

          call hsort(N,rrr,ird)

          end if

          end if

! Only select the first nbase minimal values among all N=nsample elements
! if N > nbase without sorting (sort=0)

          if(this%nbasis.gt.this%nbase.and.this%sort.eq.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) call selects(this%nbase,N,rrr,irr)

          do isa=1,this%nbase
          ird(isa)=irr(isa)
          end do

          end if

! If sort/=0, perform sorting even for nbasis=nbase, i.e. when selecting
! is not necessary

          if(this%nbasis.eq.this%nbase.and.this%sort.ne.0) then

          call hsort(this%nbase,rrr,ird)

          end if

! The nearest pair

          is1=ird(1)

! Zero-order force approximation

          do iatm = this%atom0, this%atomF
          do id=1,3
          sforce(id,iatm) = this%sforce(id,iatm,is1)
          end do
          end do

          NN=this%nbase-1

          if(NN.gt.0) then

! High-order force extrapolation

          JPVT=>safemem_realloc(JPVT,NN,.false.)
          JPVT=0

! Solve the linear least-square problem with or without weighting

          if(this%weigh.eq.0) then

          M = 3 * this%natom
          LDA = M
          LDB = max(M,NN)

          A=>safemem_realloc(A,M,NN,.false.)
          B=>safemem_realloc(B,max(NN,M),1,.false.)

! Using multiprocessor technique to fill in the matrix elements

#ifdef MPI
          A=0.d0
          B=0.d0

          jsa=3*(this%atom0-1)
#else
          jsa=0
#endif

          do iatm = this%atom0, this%atomF
          B(jsa+1:jsa+3,1)=scoord(:,iatm)-this%scoord(:,iatm,is1)
          do isa=1,NN
          iisap1=ird(isa+1)
          A(jsa+1:jsa+3,isa)=this%scoord(:,iatm,iisap1)-this%scoord(:,iatm,is1)
          end do
          jsa=jsa+3
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,M*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    call mpi_allreduce(MPI_IN_PLACE,B,LDB,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DGELSY'
          end if

          else

          jsa=0

#ifdef MPI

          jsad=0
          do iatm = this%atom0, this%atomF
          do jatm = iatm+1, this%natom
          rwrwr2=sum((scoord(1:3,iatm)-scoord(1:3,jatm))**2)
          if(rwrwr2 < this%cut) then
          jsa=jsa+3
          end if
          end do
          end do

    call MPI_Gather(jsa,1,MPI_INTEGER,jsad,1,MPI_INTEGER,0,this%mpicomm,err);

      if(this%mpirank == 0) then
      do inum=2,this%mpisize
      jsad(inum)=jsad(inum-1)+jsad(inum)
      end do
      end if

    call mpi_allreduce(MPI_IN_PLACE,jsad,this%mpisize,MPI_INTEGER,MPI_SUM,this%mpicomm,err)

      if(this%mpirank == 0) jsa=0
      do inum=1,this%mpisize-1
      if(this%mpirank == inum) jsa=jsad(inum)
      end do

          M=jsad(this%mpisize)
#else
          M = 3 * ((this%natom*(this%natom-1))/2)
#endif

          LDA = M
          LDB = max(M,NN)

          A=>safemem_realloc(A,M,NN,.false.)
          B=>safemem_realloc(B,max(NN,M),1,.false.)

          A=0.d0
          B=0.d0

          do iatm = this%atom0, this%atomF
          do jatm = iatm+1, this%natom

          rwrwr2=sum((scoord(1:3,iatm)-scoord(1:3,jatm))**2)

! Use cutting for weighting

          if(rwrwr2 < this%cut) then

          rweight=1.d0/dsqrt(rwrwr2)

          B(jsa+1:jsa+3,1) = rweight*(&
                             scoord(:,iatm)-this%scoord(:,iatm,is1) &
                            -scoord(:,jatm)+this%scoord(:,jatm,is1))
          do isa=1,NN
          iisap1=ird(isa+1)
          A(jsa+1:jsa+3,isa) = rweight*(&
                          this%scoord(:,iatm,iisap1)-this%scoord(:,iatm,is1) &
                         -this%scoord(:,jatm,iisap1)+this%scoord(:,jatm,is1))
          end do
          jsa=jsa+3
          end if
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,M*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    call mpi_allreduce(MPI_IN_PLACE,B,LDB,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#else
          M=jsa
#endif

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DGELSY'
          end if

          end if

! Extrapolate the forces in the transformed space

          do iatm = this%atom0, this%atomF
          do isa=1,NN
          iisap1=ird(isa+1)
          do id=1,3
          sforce(id,iatm) = sforce(id,iatm)+&
          B(isa,1)*(this%sforce(id,iatm,iisap1)-this%sforce(id,iatm,is1))
          end do
          end do
          end do

          end if

! Performing the inverse transformation to obtain the extrapolated forces
! in the usual coordinates

          do iatm = this%atom0, this%atomF
    force(1,iatm)=fis(1,1)*sforce(1,iatm)+fis(1,2)*sforce(2,iatm)+fis(1,3)*sforce(3,iatm)
    force(2,iatm)=fis(2,1)*sforce(1,iatm)+fis(2,2)*sforce(2,iatm)+fis(2,3)*sforce(3,iatm)
    force(3,iatm)=fis(3,1)*sforce(1,iatm)+fis(3,2)*sforce(2,iatm)+fis(3,3)*sforce(3,iatm)
          end do

    err = safemem_dealloc(A)
    err = safemem_dealloc(B)
    err = safemem_dealloc(WORK)
    err = safemem_dealloc(JPVT)
    err = safemem_dealloc(swork)

  end subroutine fce_forceb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 2-like trans=1 but within the normal equations method complemeted
     !   by the e-minimization of the norm of the solutions as well as by
     !   the ifreq-scheme to speed up the calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_forcebm(this,force,coord,iupdate,idirom)

    use safemem
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    _REAL_ :: sforce(3,this%natom),scoord(3,this%natom)

    integer :: iatm,jatm,id,err

    !LAPACK variables
    integer :: N, NN, NRHS=1, LDA, LDB, LWORK, INFO
    integer :: IPIV(this%nbase+1)
   _REAL_, pointer :: WORK(:)=>NULL()

    integer :: lswork, infos
   _REAL_, pointer :: swork(:)=>NULL()

    !uncertainties of input coordinates in A and B are assumed
    !to be negligible small
    _REAL_ :: RCOND=0.d0
    _REAL_ :: A(this%nbase+1,this%nbase+1),B(this%nbase+1,1)
    _REAL_ :: rrr(this%nbasis),rar(this%nbase)

    integer :: irr(this%nbasis),ird(this%nbase),irs(this%nbase)

    integer :: isa,jsa,iisa,jjsa

     DOUBLE PRECISION hmu(4,4),ddd,q0,q1,q2,q3,fis(3,3),wjjj(4)
     DOUBLE PRECISION xiii1,xiii2,xiii3,fiii1,fiii2,fiii3

    integer :: iupdate,idirom,ifreqs=0

    save ifreqs

     N = this%nsample

     NN=this%nbase+1

! Description of Lapack DSYSV variables

          LDA = NN
          LDB = NN

    if(iupdate.eq.0) ifreqs=0

! Calculating the current number of steps after the last updating

       if(ifreqs.eq.this%ifreq*(ifreqs/this%ifreq)) then
       idirom=0
       else
       idirom=idirom+1
       end if

! The transformation of the curent coordinate with respect to the first
! basic point is performed by minimization of the distance between them

       hmu=0.d0

! Summing up the quaternion matrix over atoms across the processes

       do iatm = this%atom0, this%atomF
       fiii1=this%coord(1,iatm,1)
       fiii2=this%coord(2,iatm,1)
       fiii3=this%coord(3,iatm,1)
       xiii1=coord(1,iatm)
       xiii2=coord(2,iatm)
       xiii3=coord(3,iatm)
       ddd=(xiii1+fiii1)**2+(xiii2+fiii2)**2+(xiii3+fiii3)**2
       hmu(1,1)=hmu(1,1)+((xiii1-fiii1)**2+(xiii2-fiii2)**2+ &
                          (xiii3-fiii3)**2)
       hmu(1,2)=hmu(1,2)+2.d0*(xiii2*fiii3-xiii3*fiii2)
       hmu(1,3)=hmu(1,3)+2.d0*(xiii3*fiii1-xiii1*fiii3)
       hmu(1,4)=hmu(1,4)+2.d0*(xiii1*fiii2-xiii2*fiii1)
       hmu(2,2)=hmu(2,2)+(ddd-2.d0*(xiii1*fiii1+fiii1*xiii1))
       hmu(2,3)=hmu(2,3)     -2.d0*(xiii1*fiii2+fiii1*xiii2)
       hmu(2,4)=hmu(2,4)     -2.d0*(xiii1*fiii3+fiii1*xiii3)
       hmu(3,3)=hmu(3,3)+(ddd-2.d0*(xiii2*fiii2+fiii2*xiii2))
       hmu(3,4)=hmu(3,4)     -2.d0*(xiii2*fiii3+fiii2*xiii3)
       hmu(4,4)=hmu(4,4)+(ddd-2.d0*(xiii3*fiii3+fiii3*xiii3))
       end do

#ifdef MPI

! Collecting the partial sums from different processes and placing the
! result into the same quaternion matrix known for all processors

    call mpi_allreduce(MPI_IN_PLACE,hmu,16,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)

#endif

! The problem is reduced to find eigenvalues and eigenvectors
! of the corresponding quaternion symmetric matrix

! First query the optimal workspace

       lswork = -1
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! Define the 4-dimensional quaternion corresponding to the smallest eigenvalue

       q0=hmu(1,1)
       q1=hmu(2,1)
       q2=hmu(3,1)
       q3=hmu(4,1)

! Determine the rotational matrix in 3-dimensional space

       fis(1,1)=q0**2+q1**2-q2**2-q3**2
       fis(1,2)=2.d0*(-q0*q3+q1*q2)
       fis(1,3)=2.d0*( q0*q2+q1*q3)
       fis(2,1)=2.d0*( q0*q3+q1*q2)
       fis(2,2)=q0**2+q2**2-q1**2-q3**2
       fis(2,3)=2.d0*(-q0*q1+q2*q3)
       fis(3,1)=2.d0*(-q0*q2+q1*q3)
       fis(3,2)=2.d0*( q0*q1+q2*q3)
       fis(3,3)=q0**2+q3**2-q1**2-q2**2

! Perform the transformations

       do iatm = 1, this%natom
   scoord(1,iatm)=fis(1,1)*coord(1,iatm)+fis(2,1)*coord(2,iatm)+fis(3,1)*coord(3,iatm)
   scoord(2,iatm)=fis(1,2)*coord(1,iatm)+fis(2,2)*coord(2,iatm)+fis(3,2)*coord(3,iatm)
   scoord(3,iatm)=fis(1,3)*coord(1,iatm)+fis(2,3)*coord(2,iatm)+fis(3,3)*coord(3,iatm)
       end do

! Forming the extended coordinate set (every ifreq steps only)

          if(ifreqs.eq.this%ifreq*(ifreqs/this%ifreq)) then

! No selecting and sorting if nbasis=nbase and sort=0

          if(this%nbasis.eq.this%nbase.and.this%sort.eq.0) then

          do isa=1,this%nbase
          ird(isa)=isa
          end do

          end if

! Calculate the distances between the current solute position and its
! counterparts from the extended basis set if nbasis > nbase or sort/=0 

          if(this%nbasis.gt.this%nbase.or.this%sort.ne.0) then

          do isa=1,N
          rrr(isa)=0.d0
          do iatm = this%atom0, this%atomF
          rrr(isa)=rrr(isa)+(scoord(1,iatm)-this%scoord(1,iatm,isa))**2+ &
                            (scoord(2,iatm)-this%scoord(2,iatm,isa))**2+ &
                            (scoord(3,iatm)-this%scoord(3,iatm,isa))**2
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,rrr,N,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

! Normalizing the distances

          do isa=1,N
          rrr(isa)=dsqrt(rrr(isa)/this%natom)
          end do

          end if

! Select the first nbase minimal values among all N=nsample elements
! if N > nbase (note that nbase <= nsample <= nbasis). In other words,
! choose the best transformed force-coordinates pairs.

          if(this%nbasis.gt.this%nbase.and.this%sort.ne.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) then

          call selects(this%nbase,N,rrr,irr)

          if(this%nbase.ne.1) then

          do isa=1,this%nbase-1
          rar(isa)=rrr(isa)
          end do

! Now perform sorting of the selected minimal values in the ascending order

          call hsort(this%nbase-1,rar,irs)

          end if

          irs(this%nbase)=this%nbase

          do isa=1,this%nbase
          ird(isa)=irr(irs(isa))
          end do

          end if

! On the beginning when N=nsample=nbase perform only sorting because
! then selecting is not necessary

          if(N.eq.this%nbase) then

          call hsort(N,rrr,ird)

          end if

          end if

! Only select the first nbase minimal values among all N=nsample elements
! if N > nbase without sorting (sort=0)

          if(this%nbasis.gt.this%nbase.and.this%sort.eq.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) call selects(this%nbase,N,rrr,irr)

          do isa=1,this%nbase
          ird(isa)=irr(isa)
          end do

          end if

! If sort/=0, perform sorting even for nbasis=nbase, i.e. when selecting
! is not necessary

          if(this%nbasis.eq.this%nbase.and.this%sort.ne.0) then

          call hsort(this%nbase,rrr,ird)

          end if

! Using multiprocessor technique to fill in the symmetric matrix elements

          do isa=1,NN-1
          iisa=ird(isa)
          do jsa=isa,NN-1
          jjsa=ird(jsa)
          A(isa,jsa)=0.d0

       do iatm = this%atom0, this%atomF
       A(isa,jsa)=A(isa,jsa)+this%scoord(1,iatm,iisa)*this%scoord(1,iatm,jjsa) &
                            +this%scoord(2,iatm,iisa)*this%scoord(2,iatm,jjsa) &
                            +this%scoord(3,iatm,iisa)*this%scoord(3,iatm,jjsa)
          end do
       A(isa,jsa)=A(isa,jsa)/this%natom
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,NN*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

! Complemented minimization of the norms of the solutions with enormsw-weight

          do isa=1,NN-1
          A(isa,isa)=A(isa,isa)+this%enormsw
          end do

! Using Lagrange method to normalize the solutions

          do isa=1,NN-1
          A(isa,NN)=-1.d0
          end do
          A(NN,NN)=0.d0

! Compute the factorization of the symmetric matrix A
! using the Bunch-Kaufman diagonal pivoting method

          LWORK = -1   ! First query the optimal workspace
          WORK => safemem_realloc(WORK,1,.false.)
          call DSYTRF( 'U', NN, A, LDA, IPIV, WORK, LWORK, INFO )
          LWORK = INT(WORK(1))
          WORK => safemem_realloc(WORK,LWORK,.false.)

! Perform the actual calculations

          call DSYTRF( 'U', NN, A, LDA, IPIV, WORK, LWORK, INFO )
 
          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DSYTRF'
          end if

! Put in memory the best subset indexes and the factorized matrix (if ifreq > 1)

          if(this%ifreq.ne.1) then
          this%irdd(1,:)=ird(:)
          this%ipvs=IPIV
          this%Ad=A
          end if
                
          else

! Take them from memory (if ifreq > 1)

          if(this%ifreq.ne.1) then
          ird(:)=this%irdd(1,:)
          IPIV=this%ipvs
          A=this%Ad
          end if

          end if

! Solve the linear least-square problem 

          do isa=1,NN-1
          iisa=ird(isa)
          B(isa,1)=0.d0
          do iatm = this%atom0, this%atomF
          B(isa,1)=B(isa,1)+scoord(1,iatm)*this%scoord(1,iatm,iisa) &
                           +scoord(2,iatm)*this%scoord(2,iatm,iisa) &
                           +scoord(3,iatm)*this%scoord(3,iatm,iisa)
          end do
          B(isa,1)=B(isa,1)/this%natom
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,B,NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

! Using Lagrange method to normalize the solutions

          B(NN,1)=-1.d0

! Solving the extended system of linear equations
! using the factorization computed above by DSYTRF

          CALL DSYTRS( 'U', NN, NRHS, A, LDA, IPIV, B, LDB, INFO )

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DSYTRS'
          end if

! Zero-order force approximation

          do iatm = this%atom0, this%atomF
          do id=1,3
          sforce(id,iatm) = 0.d0
          end do
          end do

! High-order force approximation by extrapolating 
! the forces in the transformed space

          do iatm = this%atom0, this%atomF
          do isa=1,NN-1
          iisa=ird(isa)
          do id=1,3
          sforce(id,iatm) = sforce(id,iatm)+B(isa,1)*this%sforce(id,iatm,iisa)
          end do
          end do
          end do

! Performing the inverse transformation to obtain the extrapolated forces
! in the usual coordinates

          do iatm = this%atom0, this%atomF
    force(1,iatm)=fis(1,1)*sforce(1,iatm)+fis(1,2)*sforce(2,iatm)+fis(1,3)*sforce(3,iatm)
    force(2,iatm)=fis(2,1)*sforce(1,iatm)+fis(2,2)*sforce(2,iatm)+fis(2,3)*sforce(3,iatm)
    force(3,iatm)=fis(3,1)*sforce(1,iatm)+fis(3,2)*sforce(2,iatm)+fis(3,3)*sforce(3,iatm)
          end do

! Current number of calls to this subroutine after the last updating

      if(iupdate.eq.0) iupdate=1
      ifreqs=ifreqs+1

    err = safemem_dealloc(swork)
    err = safemem_dealloc(WORK)

  end subroutine fce_forcebm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The super advanced force extrapolation method (partial ifreq-support)
!
! It includes all the techniques: normalization, neighbouring transformation
! with weighting, basic set extension, selection with possible sorting,
! as well as balancing and netforce correction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 5-individual transformation and selecting with respect to the
     !   current coordinate of each atom using a neighbouring scheme
     !   (recommended for large macromolecules of greater than 10 A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_forcesa(this,force,coord,forcnetr,iupdate,idirom)

    use safemem
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

     integer :: N, NN, NNd
     integer :: iatm,jatm,jatms,isa,jsa,iisa,jjsa,err

    !LAPACK DSYSV variables
     integer :: NRHS=1, LDA, LDB, LWORK, INFO
     integer, pointer :: IPIV(:)=>NULL()
    _REAL_, pointer :: WORK(:)=>NULL()

    !uncertainties of input coordinates in A and B are negligible
    _REAL_, pointer :: A(:,:)=>NULL(), B(:,:)=>NULL()
    _REAL_ :: RCOND=0.d0

    !LAPACK DSYEV variables
     integer :: lswork, infos
    _REAL_, pointer :: swork(:)=>NULL()

     DOUBLE PRECISION hmu(4,4),ddd,q0,q1,q2,q3,fis(3,3),wjjj(4)
     DOUBLE PRECISION xiii(3),fiii(3),fiiid(3,this%natom)
     DOUBLE PRECISION hmud(this%nbasis,4,4)

    _REAL_ :: thsscoord(3),sscoord(3),cmu(4,4)=0.d0

    !Selection and sortion arrays
     integer :: irr(this%nbasis),ird(this%nbase),irs(this%nbase)
    _REAL_ :: rrr(this%nbasis),rar(this%nbase)

     DOUBLE PRECISION rwrwr2,rwrwr2d(this%natom)

    _REAL_ :: forcnet(3),forcnetr,forcorgd,forcorgs=0.d0,forcors=0.d0

     integer :: tscut(this%natom),iupdate,idirom,ifreqs=0

     save ifreqs,forcorgs,forcors

       N = this%nsample   ! Current number of samples (NNd.le.N.le.this%nbasis)

       NNd = this%nbase   ! Number of the best basic points

       NN=NNd+1           ! Number of linear equations

! Description of Lapack DSYSV variables and memory

          LDA = NN
          LDB = NN

          A=>safemem_realloc(A,NN,NN,.false.)
          B=>safemem_realloc(B,NN,1,.false.)
          IPIV=>safemem_realloc(IPIV,NN,.false.)

   if(iupdate.eq.0) ifreqs=0

! Calculating the current number of steps after the last updating

       if(ifreqs.eq.this%ifreq*(ifreqs/this%ifreq)) then
       idirom=0
       else
       idirom=idirom+1
       end if

! Using multiprocessing to perform the individual force extrapolation
! for each separate atom of the solvent macromolecule

       do iatm = this%atom0, this%atomF

! Definition of the list of neighbours for each current atom

       tscut=0
 
       do jatm = 1,this%natom
       fiiid(:,jatm)=coord(:,jatm)-coord(:,iatm)
       rwrwr2=sum(fiiid(1:3,jatm)**2)
       if(rwrwr2<this%cut.and.jatm.ne.iatm) then
       tscut(1)=tscut(1)+1
       tscut(tscut(1)+1)=jatm
       rwrwr2d(tscut(1)+1)=1.d0/rwrwr2
       end if
       end do

! Forming the extended coordinate set (every ifreq steps only)

       if(ifreqs.eq.this%ifreq*(ifreqs/this%ifreq)) then

! Finding the optimal rotational transformation which minimizes the
! distances between the molecule in the current state and previous 
! positions from the extended coordinate set

       do isa=1,N

       hmu=0.d0

! Summing up the corresponding quaternion matrix over the neighbours atoms

       do jatms=2,tscut(1)+1
       jatm = tscut(jatms)

       fiii(:)=fiiid(:,jatm)

       xiii(:)=this%coord(:,jatm,isa)-this%coord(:,iatm,isa)

       ddd=sum((xiii(1:3)+fiii(1:3))**2)

       cmu(1,1)=sum((xiii(1:3)-fiii(1:3))**2)
       cmu(1,2)=2.d0*(xiii(2)*fiii(3)-xiii(3)*fiii(2))
       cmu(1,3)=2.d0*(xiii(3)*fiii(1)-xiii(1)*fiii(3))
       cmu(1,4)=2.d0*(xiii(1)*fiii(2)-xiii(2)*fiii(1))
       cmu(2,2)=ddd-4.d0*xiii(1)*fiii(1)
       cmu(2,3)=-2.d0*(xiii(1)*fiii(2)+xiii(2)*fiii(1))
       cmu(2,4)=-2.d0*(xiii(1)*fiii(3)+xiii(3)*fiii(1))
       cmu(3,3)=ddd-4.d0*xiii(2)*fiii(2)
       cmu(3,4)=-2.d0*(xiii(2)*fiii(3)+xiii(3)*fiii(2))
       cmu(4,4)=ddd-4.d0*xiii(3)*fiii(3)

       hmu=hmu+cmu*rwrwr2d(jatms)

       end do

       hmud(isa,:,:)=hmu(:,:) ! Put it in memory to avoid further recalculation

! The problem is reduced to find eigenvalues of the quaternion matrice

       lswork = -1      ! First query the optimal workspace
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'N', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'N', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! The minimal eigenvalue will correspond to the desired distance

       if(this%nbasis.gt.this%nbase.or.this%sort.ne.0) rrr(isa)=wjjj(1)

       end do

! Carrying out the selection and possible sorting of the extended set
! to obtain the best basic subset

! No selecting and sorting if nbasis=nbase and sort=0

          if(this%nbasis.eq.this%nbase.and.this%sort.eq.0) then

          do isa=1,this%nbase
          ird(isa)=isa
          end do

          end if

! Select the first nbase minimal values among all N=nsample elements
! if N > nbase (note that nbase <= nsample <= nbasis). In other words,
! choose the best transformed force-coordinates pairs.

          if(this%nbasis.gt.this%nbase.and.this%sort.ne.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) then

          call selects(this%nbase,N,rrr,irr)

          if(this%nbase.ne.1) then

          do isa=1,this%nbase-1
          rar(isa)=rrr(isa)
          end do

! Now perform sorting of the selected minimal values in the ascending order

          call hsort(this%nbase-1,rar,irs)

          end if

          irs(this%nbase)=this%nbase

          do isa=1,this%nbase
          ird(isa)=irr(irs(isa))
          end do

          end if

! On the beginning when N=nsample=nbase perform only sorting because
! then selecting is not necessary

          if(N.eq.this%nbase) then

          call hsort(N,rrr,ird)

          end if

          end if

! Only select the first nbase minimal values among all N=nsample elements
! if N > nbase without sorting (sort=0)

          if(this%nbasis.gt.this%nbase.and.this%sort.eq.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) call selects(this%nbase,N,rrr,irr)

          do isa=1,this%nbase
          ird(isa)=irr(isa)
          end do

          end if

! If sort/=0, perform sorting even for nbasis=nbase, i.e. when selecting
! is not necessary

          if(this%nbasis.eq.this%nbase.and.this%sort.ne.0) then

          call hsort(this%nbase,rrr,ird)

          end if

! Put in memory the best subset indexes (if ifreq > 1)

         if(this%ifreq.ne.1) then
         do isa=1,NNd
         this%irdd(iatm,isa)=ird(isa)
         end do
         end if

         else

! Take from memory the best subset indexes (if ifreq > 1)

         if(this%ifreq.ne.1) then
         do isa=1,NNd
         ird(isa)=this%irdd(iatm,isa)
         end do
         end if

         end if
 
! Having the best subset, each basic force-coordinate pairs will now be
! rotationally transformed with respect to the current position to
! minimize the distance between them and that point and thus to
! improve the quality of the force extrapolation

       do isa=1,NNd

       iisa=ird(isa)    ! Mapping from the extended set to the best subset

! Using the already calculated quaternion matrices

       if(ifreqs.eq.this%ifreq*(ifreqs/this%ifreq)) then

       hmu(:,:)=hmud(iisa,:,:)

       else

! Calculate the quaternion matrices for the new position if ifreq > 1

       hmu=0.d0

       do jatms=2,tscut(1)+1
       jatm = tscut(jatms)

       fiii(:)=fiiid(:,jatm)

       xiii(:)=this%coord(:,jatm,iisa)-this%coord(:,iatm,iisa)

       ddd=sum((xiii(1:3)+fiii(1:3))**2)

       cmu(1,1)=sum((xiii(1:3)-fiii(1:3))**2)
       cmu(1,2)=2.d0*(xiii(2)*fiii(3)-xiii(3)*fiii(2))
       cmu(1,3)=2.d0*(xiii(3)*fiii(1)-xiii(1)*fiii(3))
       cmu(1,4)=2.d0*(xiii(1)*fiii(2)-xiii(2)*fiii(1))
       cmu(2,2)=ddd-4.d0*xiii(1)*fiii(1)
       cmu(2,3)=-2.d0*(xiii(1)*fiii(2)+xiii(2)*fiii(1))
       cmu(2,4)=-2.d0*(xiii(1)*fiii(3)+xiii(3)*fiii(1))
       cmu(3,3)=ddd-4.d0*xiii(2)*fiii(2)
       cmu(3,4)=-2.d0*(xiii(2)*fiii(3)+xiii(3)*fiii(2))
       cmu(4,4)=ddd-4.d0*xiii(3)*fiii(3)

       hmu=hmu+cmu*rwrwr2d(jatms)

       end do

       end if

! The problem is reduced to find eigenvalues as well as eigenvectors
! of the obtained symmetrical quaternion matrices

       lswork = -1      ! First query the optimal workspace
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! Define the 4-dimensional quaternion corresponding to the smallest eigenvalue

       q0=hmu(1,1)
       q1=hmu(2,1)
       q2=hmu(3,1)
       q3=hmu(4,1)

! Determine the rotational matrix in 3-dimensional space

       fis(1,1)=q0**2+q1**2-q2**2-q3**2
       fis(1,2)=2.d0*(-q0*q3+q1*q2)
       fis(1,3)=2.d0*( q0*q2+q1*q3)
       fis(2,1)=2.d0*( q0*q3+q1*q2)
       fis(2,2)=q0**2+q2**2-q1**2-q3**2
       fis(2,3)=2.d0*(-q0*q1+q2*q3)
       fis(3,1)=2.d0*(-q0*q2+q1*q3)
       fis(3,2)=2.d0*( q0*q1+q2*q3)
       fis(3,3)=q0**2+q3**2-q1**2-q2**2

! Perform the coordinate transformation of each relative neigbour position

       do jatms=2,tscut(1)+1
       jatm = tscut(jatms)

       xiii(:)=(this%coord(:,jatm,iisa)-this%coord(:,iatm,iisa))* &
                dsqrt(rwrwr2d(jatms))

       this%scoord(1,jatm,iisa)=sum(fis(1:3,1)*xiii(1:3))
       this%scoord(2,jatm,iisa)=sum(fis(1:3,2)*xiii(1:3))
       this%scoord(3,jatm,iisa)=sum(fis(1:3,3)*xiii(1:3))

       end do

! Perform the force transformation for the current atom

    this%sforce(1,iatm,iisa)=sum(fis(1:3,1)*this%force(1:3,iatm,iisa))
    this%sforce(2,iatm,iisa)=sum(fis(1:3,2)*this%force(1:3,iatm,iisa))
    this%sforce(3,iatm,iisa)=sum(fis(1:3,3)*this%force(1:3,iatm,iisa))

       end do

! Solve the linear least-square problem in the transformed space

       A=0.d0
       B=0.d0

! Filling in the corresponding elements of symmetric matrix A and vector B

          do jatms=2,tscut(1)+1
          jatm = tscut(jatms)

          sscoord(:)=fiiid(:,jatm)*dsqrt(rwrwr2d(jatms))

          do isa=1,NN-1

          iisa=ird(isa)

          B(isa,1)=B(isa,1)+sum(sscoord(1:3)*this%scoord(1:3,jatm,iisa))

          thsscoord(:)=this%scoord(:,jatm,iisa)

          do jsa=isa,NN-1

          jjsa=ird(jsa)

          A(isa,jsa)=A(isa,jsa)+sum(thsscoord(1:3)*this%scoord(1:3,jatm,jjsa))

          end do
          end do

          end do

! Renormalization of the matrices on the number of neighbours

          A=A/tscut(1)
          B=B/tscut(1)

! Balancing minimization of the norms of the solutions with enormsw-weight

          do isa=1,NN-1
          A(isa,isa)=A(isa,isa)+this%enormsw
          end do

! Using Lagrange method to normalize the solutions

          do isa=1,NN-1
          A(isa,NN)=-1.d0
          end do
          A(NN,NN)=0.d0

          B(NN,1)=-1.d0

! Solving the extended system of linear equations

          LWORK = -1      ! First query the optimal workspace
          WORK => safemem_realloc(WORK,1,.false.)
          call DSYSV( 'U', NN, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
          LWORK = INT(WORK(1))
          WORK => safemem_realloc(WORK,LWORK,.false.)

! Perform the actual calculations

          call DSYSV( 'U', NN, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DSYSV'
          end if

! Extrapolate the forces in the transformed space

          force(:,iatm) = 0.d0     ! Zero-order force approximation

! High-order force extrapolation

          do isa=1,NN-1
          iisa=ird(isa)

          force(:,iatm) = force(:,iatm)+B(isa,1)*this%sforce(:,iatm,iisa)

          end do

          end do

! Calculate the net force

       forcnet=0.d0

       forcorgd=0.d0

! Summing up the partial forces and their squared values over atoms
 
       do iatm = this%atom0, this%atomF

       forcnet(:)=forcnet(:)+force(:,iatm)

       forcorgd=forcorgd+sum(force(1:3,iatm)**2)

       end do

#ifdef MPI

! Collecting the partial sums from different processes and placing
! the result into the same variables known for all processors

    call mpi_allreduce(MPI_IN_PLACE,forcnet,3,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    call mpi_allreduce(MPI_IN_PLACE,forcorgd,1,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)

#endif

! Evaluate the correcting force

       forcnet(:)=forcnet(:)/this%natom

! Accumulate the squared correcting and initial forces per atom

       forcors=forcors+sum(forcnet(1:3)**2)

       forcorgs=forcorgs+forcorgd/this%natom

! Calculate the ratio of averaged correcting to initial forces

       forcnetr=dsqrt(forcors/forcorgs)

! Perform the net force correction

       if(this%ntfrcor.ne.0) then

       do iatm = this%atom0, this%atomF
       force(:,iatm)=force(:,iatm)-forcnet(:)
       end do

       end if

! Current number of calls to this subroutine after the last updating

      if(iupdate.eq.0) iupdate=1
      ifreqs=ifreqs+1

    err = safemem_dealloc(A)
    err = safemem_dealloc(B)
    err = safemem_dealloc(IPIV)
    err = safemem_dealloc(WORK)
    err = safemem_dealloc(swork)

  end subroutine fce_forcesa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The super advanced force extrapolation method with full ifreq-support
!
! It includes all the techniques: normalization, neighbouring transformation
! with weighting, basic set extension, selection with possible sorting,
! as well as balancing, netforce correction, and full ifreq-technique
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 6-individual transformation and selecting with respect to the
     !   post-coordinate of each atom using a neighbouring scheme
     !   (recommended for large macromolecules and large basic set)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_forcesan(this,force,coord,forcnetr,iupdate,idirom)

    use safemem
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    _REAL_ :: sforce(3,this%natom),scoord(3,this%natom)
    _REAL_ :: coordp(3,this%natom)

     integer :: N, NN, NNd
     integer :: iatm,jatm,jatms,isa,jsa,iisa,err

    !LAPACK DSYSV variables
     integer :: NRHS=1, LDA, LDB, LWORK, INFO
     integer :: IPIV(this%nbase+1)
    _REAL_, pointer :: WORK(:)=>NULL()

    !uncertainties of input coordinates in A and B are negligible
    _REAL_ :: A(this%nbase+1,this%nbase+1),B(this%nbase+1,1)
    _REAL_ :: RCOND=0.d0

    !LAPACK DSYEV variables
     integer :: lswork, infos
    _REAL_, pointer :: swork(:)=>NULL()

     DOUBLE PRECISION hmu(4,4),ddd,q0,q1,q2,q3,fis(3,3),wjjj(4)
     DOUBLE PRECISION xiii(3),fiii(3),fiiid(3,this%natom)
     DOUBLE PRECISION hmud(this%nbasis,4,4)

    _REAL_ :: tscoords(3),thsscoord(3),cmu(4,4)=0.d0

    !Selection and sortion arrays
     integer :: irr(this%nbasis),ird(this%nbase),irs(this%nbase)
    _REAL_ :: rrr(this%nbasis),rar(this%nbase)

     DOUBLE PRECISION rwrwr2,rwrwr2d(this%natom)

    _REAL_ :: forcnet(3),forcnetr,forcorgd,forcorgs=0.d0,forcors=0.d0

     integer :: tscut(this%natom),iupdate,idirom,ifreqs=0

     save ifreqs,forcorgs,forcors

       N = this%nsample   ! Current number of samples (NNd.le.N.le.this%nbasis)

       NNd = this%nbase   ! Number of the best basic points

       NN=NNd+1           ! Number of linear equations

! Description of Lapack DSYSV variables

          LDA = NN
          LDB = NN

   if(iupdate.eq.0) ifreqs=0

! Calculating the current number of steps after the last updating

       if(ifreqs.eq.this%ifreq*(ifreqs/this%ifreq)) then
       idirom=0
       else
       idirom=idirom+1
       end if

! Using multiprocessing to perform the individual force extrapolation
! for each separate atom of the solvent macromolecule

       do iatm = this%atom0, this%atomF

! Forming the extended coordinate set (every ifreq steps only)

       if(ifreqs.eq.this%ifreq*(ifreqs/this%ifreq)) then

! Definition of the list of neighbours for each current atom

       tscut=0
 
       do jatm = 1,this%natom
       fiiid(:,jatm)=coord(:,jatm)-coord(:,iatm)
       rwrwr2=sum(fiiid(1:3,jatm)**2)
       if(rwrwr2<this%cut.and.jatm.ne.iatm) then
       tscut(1)=tscut(1)+1
       tscut(tscut(1)+1)=jatm
       rwrwr2d(tscut(1)+1)=1.d0/rwrwr2
       end if
       end do

! Finding the optimal rotational transformation which minimizes the
! distances between the molecule in the current post-state and
! previous positions from the extended coordinate set

       do isa=1,N

       hmu=0.d0

! Summing up the corresponding quaternion matrix over the neighbours atoms

       do jatms=2,tscut(1)+1
       jatm = tscut(jatms)

       fiii(:)=fiiid(:,jatm)

       xiii(:)=this%coord(:,jatm,isa)-this%coord(:,iatm,isa)

       ddd=sum((xiii(1:3)+fiii(1:3))**2)

       cmu(1,1)=sum((xiii(1:3)-fiii(1:3))**2)
       cmu(1,2)=2.d0*(xiii(2)*fiii(3)-xiii(3)*fiii(2))
       cmu(1,3)=2.d0*(xiii(3)*fiii(1)-xiii(1)*fiii(3))
       cmu(1,4)=2.d0*(xiii(1)*fiii(2)-xiii(2)*fiii(1))
       cmu(2,2)=ddd-4.d0*xiii(1)*fiii(1)
       cmu(2,3)=-2.d0*(xiii(1)*fiii(2)+xiii(2)*fiii(1))
       cmu(2,4)=-2.d0*(xiii(1)*fiii(3)+xiii(3)*fiii(1))
       cmu(3,3)=ddd-4.d0*xiii(2)*fiii(2)
       cmu(3,4)=-2.d0*(xiii(2)*fiii(3)+xiii(3)*fiii(2))
       cmu(4,4)=ddd-4.d0*xiii(3)*fiii(3)

       hmu=hmu+cmu*rwrwr2d(jatms)

       end do

       hmud(isa,:,:)=hmu(:,:) ! Put it in memory to avoid further recalculation

! The problem is reduced to find eigenvalues of the quaternion matrice

       lswork = -1      ! First query the optimal workspace
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'N', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'N', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! The minimal eigenvalue will correspond to the desired minimal distance

       if(this%nbasis.gt.this%nbase.or.this%sort.ne.0) rrr(isa)=wjjj(1)

       end do

! Carrying out the selection and possible sorting of the extended set
! to obtain the best basic subset

          if(this%nbasis.eq.this%nbase.and.this%sort.eq.0) then

          do isa=1,this%nbase
          ird(isa)=isa
          end do

          end if

! Select the first nbase minimal values among all N=nsample elements
! if N > nbase (note that nbase <= nsample <= nbasis). In other words,
! choose the best transformed force-coordinates pairs.

          if(this%nbasis.gt.this%nbase.and.this%sort.ne.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) then

          call selects(this%nbase,N,rrr,irr)

          if(this%nbase.ne.1) then

          do isa=1,this%nbase-1
          rar(isa)=rrr(isa)
          end do

! Now perform sorting of the selected minimal values in the ascending order

          call hsort(this%nbase-1,rar,irs)

          end if

          irs(this%nbase)=this%nbase

          do isa=1,this%nbase
          ird(isa)=irr(irs(isa))
          end do

          end if

! On the beginning when N=nsample=nbase perform only sorting because
! then selecting is not necessary

          if(N.eq.this%nbase) then

          call hsort(N,rrr,ird)

          end if

          end if

! Only select the first nbase minimal values among all N=nsample elements
! if N > nbase without sorting (sort=0)

          if(this%nbasis.gt.this%nbase.and.this%sort.eq.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) call selects(this%nbase,N,rrr,irr)

          do isa=1,this%nbase
          ird(isa)=irr(isa)
          end do

          end if

! If sort/=0, perform sorting even for nbasis=nbase, i.e. when selecting
! is not necessary

          if(this%nbasis.eq.this%nbase.and.this%sort.ne.0) then

          call hsort(this%nbase,rrr,ird)

          end if

! Having the best subset, each basic force-coordinate pairs 
! will now be rotationally transformed accordingly

       do isa=1,NNd

       iisa=ird(isa)    ! Mapping from the extended set to the best subset

! Using the already calculated quaternion matrices

       hmu(:,:)=hmud(iisa,:,:)

! The problem is reduced to find eigenvalues as well as eigenvectors
! of the obtained symmetrical quaternion matrices

       lswork = -1      ! First query the optimal workspace
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! Define the 4-dimensional quaternion corresponding to the smallest eigenvalue

       q0=hmu(1,1)
       q1=hmu(2,1)
       q2=hmu(3,1)
       q3=hmu(4,1)

! Determine the rotational matrix in 3-dimensional space

       fis(1,1)=q0**2+q1**2-q2**2-q3**2
       fis(1,2)=2.d0*(-q0*q3+q1*q2)
       fis(1,3)=2.d0*( q0*q2+q1*q3)
       fis(2,1)=2.d0*( q0*q3+q1*q2)
       fis(2,2)=q0**2+q2**2-q1**2-q3**2
       fis(2,3)=2.d0*(-q0*q1+q2*q3)
       fis(3,1)=2.d0*(-q0*q2+q1*q3)
       fis(3,2)=2.d0*( q0*q1+q2*q3)
       fis(3,3)=q0**2+q3**2-q1**2-q2**2

! Perform the coordinate transformation of each relative neigbour position

       do jatms=2,tscut(1)+1
       jatm = tscut(jatms)

       xiii(:)=(this%coord(:,jatm,iisa)-this%coord(:,iatm,iisa))* &
                dsqrt(rwrwr2d(jatms))

       this%scoord(1,jatm,isa)=sum(fis(1:3,1)*xiii(1:3))
       this%scoord(2,jatm,isa)=sum(fis(1:3,2)*xiii(1:3))
       this%scoord(3,jatm,isa)=sum(fis(1:3,3)*xiii(1:3))

       end do

! Perform the force transformation for the current atom

       this%sforce(1,iatm,isa)=sum(fis(1:3,1)*this%force(1:3,iatm,iisa))
       this%sforce(2,iatm,isa)=sum(fis(1:3,2)*this%force(1:3,iatm,iisa))
       this%sforce(3,iatm,isa)=sum(fis(1:3,3)*this%force(1:3,iatm,iisa))

! Remember the transformation matrices

       this%sfis(:,:,iatm-this%atom0+1,isa)=fis(:,:)

       end do

! Solve the linear least-square problem 

          A=0.d0

! Filling in the corresponding elements of symmetric matrix A

          do jatms=2,tscut(1)+1
          jatm = tscut(jatms)

          do isa=1,NN-1

          thsscoord(:)=this%scoord(:,jatm,isa)

          do jsa=isa,NN-1

          A(isa,jsa)=A(isa,jsa)+sum(thsscoord(1:3)*this%scoord(1:3,jatm,jsa))

          end do
          end do

          end do

! Renormalization of the matrices on the number of neighbours

          A=A/tscut(1)

! Balancing minimization of the norms of the solutions with enormsw-weight

          do isa=1,NN-1
          A(isa,isa)=A(isa,isa)+this%enormsw
          end do

! Using Lagrange method to normalize the solutions

          do isa=1,NN-1
          A(isa,NN)=-1.d0
          end do
          A(NN,NN)=0.d0

! Compute the factorization of the symmetric matrix A
! using the Bunch-Kaufman diagonal pivoting method

          LWORK = -1   ! First query the optimal workspace
          WORK => safemem_realloc(WORK,1,.false.)
          call DSYTRF( 'U', NN, A, LDA, IPIV, WORK, LWORK, INFO )
          LWORK = INT(WORK(1))
          WORK => safemem_realloc(WORK,LWORK,.false.)

! Perform the actual calculations

          call DSYTRF( 'U', NN, A, LDA, IPIV, WORK, LWORK, INFO )
 
          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DSYTRF'
          end if

! Put in memory the best subset indexes and factorized matrix (if ifreq > 1)

       if(this%ifreq.ne.1) then
       this%irdd(iatm,:)=ird(:)
       this%ipvss(iatm-this%atom0+1,:)=IPIV(:)
       this%Ads(iatm-this%atom0+1,:,:)=A(:,:)
       end if

! Put in memory other necessary quantities (if ifreq > 1)

     if(this%ifreq.ne.1) then
     this%tscuts(iatm-this%atom0+1,:)=tscut(:)
     this%coordps(iatm-this%atom0+1,:,:)=coord(:,:)
     this%srwrwr2d(iatm-this%atom0+1,:)=rwrwr2d(:)
     end if

! On the beginning of each current frequency interval
! the post-state coincides with the current one

       fis=0.d0
       fis(1,1)=1.d0
       fis(2,2)=1.d0
       fis(3,3)=1.d0

       do jatms=2,tscut(1)+1
       jatm = tscut(jatms)

       scoord(:,jatm)=(coord(:,jatm)-coord(:,iatm))*dsqrt(rwrwr2d(jatms))

       end do

       else

! Get from memory the necessary quantities (if ifreq > 1)

     if(this%ifreq.ne.1) then
     tscut(:)=this%tscuts(iatm-this%atom0+1,:)
     coordp(:,:)=this%coordps(iatm-this%atom0+1,:,:)
     rwrwr2d(:)=this%srwrwr2d(iatm-this%atom0+1,:)
     end if

! The transformation of the curent coordinate with respect to the post-
! state is performed by minimization of the distance between them

       hmu=0.d0

! Summing up the corresponding quaternion matrix over the neighbours atoms

       do jatms=2,tscut(1)+1
       jatm = tscut(jatms)

       fiii(:)=coordp(:,jatm)-coordp(:,iatm)

       xiii(:)=coord(:,jatm)-coord(:,iatm)

       ddd=sum((xiii(1:3)+fiii(1:3))**2)

       cmu(1,1)=sum((xiii(1:3)-fiii(1:3))**2)
       cmu(1,2)=2.d0*(xiii(2)*fiii(3)-xiii(3)*fiii(2))
       cmu(1,3)=2.d0*(xiii(3)*fiii(1)-xiii(1)*fiii(3))
       cmu(1,4)=2.d0*(xiii(1)*fiii(2)-xiii(2)*fiii(1))
       cmu(2,2)=ddd-4.d0*xiii(1)*fiii(1)
       cmu(2,3)=-2.d0*(xiii(1)*fiii(2)+xiii(2)*fiii(1))
       cmu(2,4)=-2.d0*(xiii(1)*fiii(3)+xiii(3)*fiii(1))
       cmu(3,3)=ddd-4.d0*xiii(2)*fiii(2)
       cmu(3,4)=-2.d0*(xiii(2)*fiii(3)+xiii(3)*fiii(2))
       cmu(4,4)=ddd-4.d0*xiii(3)*fiii(3)

       hmu=hmu+cmu*rwrwr2d(jatms)

       end do

! The problem is reduced to find eigenvalues and eigenvectors
! of the quaternion symmetric matrice

       lswork = -1      ! First query the optimal workspace
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! Define the 4-dimensional quaternion corresponding to the smallest eigenvalue

       q0=hmu(1,1)
       q1=hmu(2,1)
       q2=hmu(3,1)
       q3=hmu(4,1)

! Determine the rotational matrix in 3-dimensional space

       fis(1,1)=q0**2+q1**2-q2**2-q3**2
       fis(1,2)=2.d0*(-q0*q3+q1*q2)
       fis(1,3)=2.d0*( q0*q2+q1*q3)
       fis(2,1)=2.d0*( q0*q3+q1*q2)
       fis(2,2)=q0**2+q2**2-q1**2-q3**2
       fis(2,3)=2.d0*(-q0*q1+q2*q3)
       fis(3,1)=2.d0*(-q0*q2+q1*q3)
       fis(3,2)=2.d0*( q0*q1+q2*q3)
       fis(3,3)=q0**2+q3**2-q1**2-q2**2

! Perform the coordinate transformation of each relative neigbour position

       do jatms=2,tscut(1)+1
       jatm = tscut(jatms)

       xiii(:)=(coord(:,jatm)-coord(:,iatm))*dsqrt(rwrwr2d(jatms))

       scoord(1,jatm)=sum(fis(1:3,1)*xiii(1:3))
       scoord(2,jatm)=sum(fis(1:3,2)*xiii(1:3))
       scoord(3,jatm)=sum(fis(1:3,3)*xiii(1:3))

       end do

! Take from memory the best subset indexes and factorized matrix (if ifreq > 1)

     if(this%ifreq.ne.1) then
     ird(:)=this%irdd(iatm,:)
     IPIV(:)=this%ipvss(iatm-this%atom0+1,:)
     A(:,:)=this%Ads(iatm-this%atom0+1,:,:)
     end if

       end if

          B=0.d0

! Filling in the corresponding elements of vector B

          do jatms=2,tscut(1)+1
          jatm = tscut(jatms)

          do isa=1,NN-1

! First reproduce the transformed coordinate set via the rotational matrices

          iisa=ird(isa)    
          xiii(:)=(this%coord(:,jatm,iisa)-this%coord(:,iatm,iisa))* &
                   dsqrt(rwrwr2d(jatms))

          tscoords(1)=sum(this%sfis(1:3,1,iatm-this%atom0+1,isa)*xiii(1:3))
          tscoords(2)=sum(this%sfis(1:3,2,iatm-this%atom0+1,isa)*xiii(1:3))
          tscoords(3)=sum(this%sfis(1:3,3,iatm-this%atom0+1,isa)*xiii(1:3))

! Calculate the vector B

          B(isa,1)=B(isa,1)+sum(scoord(1:3,jatm)*tscoords(1:3))

          end do

          end do

! Renormalization of the vector on the number of neighbours

          B=B/tscut(1)

! Using Lagrange method to normalize the solutions

          B(NN,1)=-1.d0

! Solving the extended system of linear equations
! using the factorization computed above by DSYTRF

          CALL DSYTRS( 'U', NN, NRHS, A, LDA, IPIV, B, LDB, INFO )

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DSYTRS'
          end if

! Extrapolate the forces in the transformed space

          sforce(:,iatm) = 0.d0     ! Zero-order force approximation

! High-order force approximation by extrapolating 
! the forces in the transformed space

          do isa=1,NN-1

          sforce(:,iatm) = sforce(:,iatm)+B(isa,1)*this%sforce(:,iatm,isa)

          end do

! Performing the inverse transformation to obtain the extrapolated forces
! in the usual coordinates

    force(1,iatm)=fis(1,1)*sforce(1,iatm)+fis(1,2)*sforce(2,iatm)+fis(1,3)*sforce(3,iatm)
    force(2,iatm)=fis(2,1)*sforce(1,iatm)+fis(2,2)*sforce(2,iatm)+fis(2,3)*sforce(3,iatm)
    force(3,iatm)=fis(3,1)*sforce(1,iatm)+fis(3,2)*sforce(2,iatm)+fis(3,3)*sforce(3,iatm)

          end do

! Calculate the net force

       forcnet=0.d0

       forcorgd=0.d0

! Summing up the partial forces and their squared values over atoms
 
       do iatm = this%atom0, this%atomF

       forcnet(:)=forcnet(:)+force(:,iatm)

       forcorgd=forcorgd+sum(force(1:3,iatm)**2)

       end do

#ifdef MPI

! Collecting the partial sums from different processes and placing
! the result into the same variables known for all processors

    call mpi_allreduce(MPI_IN_PLACE,forcnet,3,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    call mpi_allreduce(MPI_IN_PLACE,forcorgd,1,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)

#endif

! Evaluate the correcting force

       forcnet(:)=forcnet(:)/this%natom

! Accumulate the squared correcting and initial forces per atom

       forcors=forcors+sum(forcnet(1:3)**2)

       forcorgs=forcorgs+forcorgd/this%natom

! Calculate the ratio of averaged correcting to initial forces

       forcnetr=dsqrt(forcors/forcorgs)

! Perform the net force correction

       if(this%ntfrcor.ne.0) then

       do iatm = this%atom0, this%atomF
       force(:,iatm)=force(:,iatm)-forcnet(:)
       end do

       end if

! Current number of calls to this subroutine after the last updating

      if(iupdate.eq.0) iupdate=1
      ifreqs=ifreqs+1

    err = safemem_dealloc(swork)
    err = safemem_dealloc(WORK)

  end subroutine fce_forcesan

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 4-no normalization, coordinate transformation, weighting, selecting,
     !   and sorting, but with individual force extrapolation and possible
     !   cutting-off the neighbours [original AMBER11 version]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_force(this,force,coord)
    use safemem
    implicit none
    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    integer :: iatm,jatm,id,err
    
    !LAPACK variables
    integer :: M,N,NRHS, LDA,LDB, LWORK,RANK=0, INFO=0
    integer, pointer :: JPVT(:)=>NULL()
    !not really sure what value this should be...
    _REAL_ :: RCOND=0d0
    _REAL_, pointer :: A(:,:)=>NULL(),B(:,:)=>NULL(),work(:)=>NULL()

    _REAL_,external :: ddot

!!$    write(6,*) "FCE_FORCE", this%atom0, this%atomF, this%cut

    force=0
    N = this%nbasis
    JPVT=>safemem_realloc(JPVT,N,.false.)
    JPVT=0
    !the cutoff and non-cutoff situations are handled differently
    if(this%cut>0) then
       !update cutoff list if requested
       call nlist(this)

       !iterate over each atom and calculate the extrapolated force
       do iatm = this%atom0, this%atomF
!          call orient(this)

          M = 3*(this%cutlist(1,iatm)+1)
          LDA = M
          LDB = max(M,N)
          NRHS = 1

          A=>safemem_realloc(A,M,N,.false.)
          B=>safemem_realloc(B,max(N,M),1,.false.)

          !the target atom is never in the cutlist so put it in the first index
          B(1:3,1) = coord(:,iatm)
          A(1:3,:) = this%coord(:,iatm,:)

          !now add the atoms from within the cutoff
          do jatm=2,this%cutlist(1,iatm)+1
             B(jatm*3-2: jatm*3,1) = coord(:,this%cutlist(jatm,iatm))
             A(jatm*3-2: jatm*3,:) = this%coord(:,this%cutlist(jatm,iatm),:)
          end do

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )
          
          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )

          !the force in each dimension is now a dot product
          do id=1,3
             force(id,iatm) = ddot(N,this%force(id,iatm,:),1,B(1:N,1),1)
          end do
!!$          write(6,*) "FCE FORCE",iatm, force(:,iatm)
 !         call unorient(this)
       end do
    else !cutoff
       
    end if !cutoff
    err = safemem_dealloc(A)
    err = safemem_dealloc(B)
    err = safemem_dealloc(WORK)
    err = safemem_dealloc(JPVT)
  end subroutine fce_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private subroutines and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Very simple method to generate a cutoff list.  The resulting matrix is 2D is
!!!lists the number of neighbour atoms for each atom (first row) and then the
!!!atom numbers themselves.
!!!E.g.
!!!   1 2 1
!!!   2 1 2
!!!   0 3 0
!!!   0 0 0
!!!For this three atom system, atoms 1 and 2, 3 and 2 are neighbours.  I.e.
!!!atoms 1 and 3 each have one neighbour and 2 has two.  The number of
!!!neighbours is list in the first row.  The subsequent rows list the ids of the
!!!neighbours.
!!!IN:
!!!   this : FCE object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nlist(this)
  implicit none
  type(fce),intent(inout) :: this
  integer :: id,iatm1,iatm2
#ifdef RISM_DEBUG
  write(6,*) "NLIST"; call flush(6)
#endif /*RISM_DEBUG*/

  this%cutlist=0
  do iatm1 = 1,this%natom
     do iatm2 = iatm1+1,this%natom
        if( sum((this%coord(1:3,iatm1,1) - this%coord(1:3,iatm2,1))**2) < this%cut)then
           this%cutlist(1,iatm1) = this%cutlist(1,iatm1) +1
           this%cutlist(1,iatm2) = this%cutlist(1,iatm2) +1
           this%cutlist(this%cutlist(1,iatm1)+1,iatm1) = iatm2
           this%cutlist(this%cutlist(1,iatm2)+1,iatm2) = iatm1
        end if
     end do
  end do
end subroutine nlist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Estimate the accuracy of the extrapolation by measuring the     !
!      difference between the exact and extrapolated forces at the     !
!      end of the outer time interval                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_estimate(this,force,coord,forcem,deviat,forcnetr,iupdate,idirom)
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(in) :: force(3,this%natom), coord(3,this%natom)

    _REAL_ forcem(3,this%natom),stforce(3,this%natom),sas,sasa

    _REAL_ :: forcnetr

    DOUBLE PRECISION deviat

     integer :: iaa

     integer :: iupdate,idirom

#ifdef MPI
    _REAL_ ssas,ssasa
     integer :: err
#endif

    save sas,sasa

        if(this%nsample.eq.0) then
        sas=0.d0
        sasa=0.d0
        end if

! Calculating the extrapolated 3D-RISM forces

       if(this%nsample.ge.this%nbase) then

       if(this%trans.eq.0) call fce_forcea(this,forcem,coord)

       if(this%trans.eq.1) call fce_forceb(this,forcem,coord)

       if(this%trans.eq.2) call fce_forcebm(this,forcem,coord,iupdate,idirom)

       if(this%trans.eq.3) call fce_forcebm(this,forcem,coord,iupdate,idirom)

       if(this%trans.eq.4) call fce_force(this,forcem,coord)

       if(this%trans.eq.5) call fce_forcesa(this,forcem,coord,forcnetr,iupdate,idirom)

       if(this%trans.eq.6) call fce_forcesan(this,forcem,coord,forcnetr,iupdate,idirom)

! Summing up the exact 3D-RISM forces from different processes and
! making the result to be known for all processors

#ifdef MPI
    call mpi_allreduce(force,stforce,3*this%natom,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#else
    stforce=force
#endif

! Calculating the relative force deviations using multiprocessing

       do iaa = this%atom0, this%atomF

          sas=sas+(forcem(1,iaa)-stforce(1,iaa))**2 &
                 +(forcem(2,iaa)-stforce(2,iaa))**2 &
                 +(forcem(3,iaa)-stforce(3,iaa))**2

          sasa=sasa+stforce(1,iaa)**2+stforce(2,iaa)**2+stforce(3,iaa)**2

       end do

#ifdef MPI

! Summing up the deviations and collecting them within the root processes

    call mpi_reduce(sas,ssas,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,this%mpicomm,err)
    call mpi_reduce(sasa,ssasa,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,this%mpicomm,err)

       if(this%mpirank==0) then

       deviat=dsqrt(ssas/ssasa)

       end if
#else
       deviat=dsqrt(sas/sasa)
#endif

       end if

end subroutine fce_estimate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
! Modified Quick Select Algorithm
!************************************************************************

! The input array arr will be rearranged to have the kth smallest value
! in location arr(k), with all smaller elements moved to arr(1:k-1)
! (in arbitrary order) and all larger elements in arr(k+1:) 

      subroutine selects(k,n,arr,irr)
      INTEGER k,n
      double precision arr(*)
      INTEGER i,ir,j,l,mid,irr(*)
      double precision a,temp
      INTEGER ia,itemp
      l=1
      ir=n
1     if(ir-l.le.1)then
        if(ir-l.eq.1)then
          if(arr(ir).lt.arr(l))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp

            itemp=irr(l)
            irr(l)=irr(ir)
            irr(ir)=itemp

          endif
        endif
        return
      else
        mid=(l+ir)/2
        temp=arr(mid)
        arr(mid)=arr(l+1)
        arr(l+1)=temp

        itemp=irr(mid)
        irr(mid)=irr(l+1)
        irr(l+1)=itemp

        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp

         itemp=irr(l+1)
         irr(l+1)=irr(ir)
         irr(ir)=itemp

        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp

            itemp=irr(l)
            irr(l)=irr(ir)
            irr(ir)=itemp

        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp

          itemp=irr(l+1)
          irr(l+1)=irr(l)
          irr(l)=itemp

        endif
        i=l+1
        j=ir
        a=arr(l)

        ia=irr(l)

3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp

        itemp=irr(i)
        irr(i)=irr(j)
        irr(j)=itemp

        goto 3
5       arr(l)=arr(j)

        irr(l)=irr(j)

        arr(j)=a

        irr(j)=ia

        if(j.ge.k)ir=j-1
        if(j.le.k)l=i
      endif
      goto 1
end subroutine selects
!  (C) Copr. 1986-92 Numerical Recipes Software

!************************************************************************
! Super Efficient Sorting (Adv. Eng. Software, 1984, Vol.6, No. 4, p.198
!************************************************************************
      subroutine hsort(n,data,list)
      integer n
      integer list(*)
      double precision data(*)
      integer maxstk, ncut
      parameter (maxstk=32, ncut=12)
      integer ll, lr, lm, nl, nr, ltemp, ist, i, j, k
      integer lstack(maxstk), rstack(maxstk)
      double precision guess, value
      do i=1,n
      list(i)=i
      end do
      ll=1
      lr=n
      ist=0
   10 continue
      if((lr-ll).ge.ncut) then
      nl=ll
      nr=lr
      lm=(ll+lr)/2
      guess=data(list(lm))
   20 continue
      if(data(list(nl)).lt.guess) then
      nl=nl+1
      goto 20
      end if
   30 continue
      if(guess.lt.data(list(nr))) then
      nr=nr-1
      goto 30
      end if
      if(nl.lt.(nr-1)) then
      ltemp=list(nl)
      list(nl)=list(nr)
      list(nr)=ltemp
      nl=nl+1
      nr=nr-1
      goto 20
      end if
      if(nl.le.nr) then
      if(nl.lt.nr) then
      ltemp=list(nl)
      list(nl)=list(nr)
      list(nr)=ltemp
      end if
      nl=nl+1
      nr=nr-1
      end if
      ist=ist+1
      if(nr.lt.lm) then
      lstack(ist)=nl
      rstack(ist)=lr
      lr=nr
      else
      lstack(ist)=ll
      rstack(ist)=nr
      ll=nl
      end if
      goto 10
      end if
      if(ist.ne.0) then
      ll=lstack(ist)
      lr=rstack(ist)
      ist=ist-1
      goto 10
      end if
      j=1
      k=list(1)
      value=data(k)
      do 40 i=2,min(n,ncut)
      if(data(list(i)).lt.value) then
      j=i
      value=data(list(i))
      end if
   40 continue
      list(1)=list(j)
      list(j)=k
      do 60 i=2,n
      j=i
      k=list(i)
      value=data(k)
   50 continue
      if(value.lt.data(list(j-1))) then
      list(j)=list(j-1)
      j=j-1
      goto 50
      end if
      list(j)=k
   60 continue
end subroutine hsort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fce_readbasis(this)
  implicit none
  type(fce),intent(inout) :: this
  integer :: i, j, k
  character(80) :: buf

  !#ifdef MPI
  !  include 'mpif.h'

  if(this%mpirank == 0) then
     
     open(unit=this%bwrtunit,FILE='fcecoord.rst',STATUS='OLD')
     open(unit=this%bwrtunit+1,FILE='fceforce.rst',STATUS='OLD')

     read(this%bwrtunit,*) this%nsample, this%cstep
     read(this%bwrtunit+1,*) buf

     do i=1,this%nbasis
        
        read(this%bwrtunit,*) buf
        read(this%bwrtunit+1,*) buf
        
        do j=1,this%natom
           read(this%bwrtunit,*) this%coord(1,j,i), this%coord(2,j,i), this%coord(3,j,i)
           read(this%bwrtunit+1,*) this%force(1,j,i), this%force(2,j,i), this%force(3,j,i)
        enddo

     enddo

     close(this%bwrtunit)
     close(this%bwrtunit+1)
     
  end if

  !#else

  !#endif

end subroutine fce_readbasis

subroutine fce_wrtbasis(this, cstepInInterval)
  implicit none
  type(fce),intent(inout) :: this
  integer :: cstepInInterval, i, j, k

!#ifdef MPI
!  include 'mpif.h'

  if(this%mpirank == 0) then
 
     open(unit=this%bwrtunit,FILE='fcecoord.rst')
     open(unit=this%bwrtunit+1,FILE='fceforce.rst')

     write(this%bwrtunit,'(I5,A1,I10)') this%nsample, " ", cstepInInterval
     write(this%bwrtunit+1,'(I5,A1,I10)') this%nsample, " ", cstepInInterval

     do i=1,this%nbasis
        write(this%bwrtunit,'(A2,I5)') "# ", i
        write(this%bwrtunit+1,'(A2,I5)') "# ", i
        do j=1,this%natom
           write(this%bwrtunit,'(F12.5,A1,F12.5,A1,F12.5)') this%coord(1,j,i), " ", this%coord(2,j,i), " ", this%coord(3,j,i)
           write(this%bwrtunit+1,'(F12.5,A1,F12.5,A1,F12.5)') this%force(1,j,i), " ", this%force(2,j,i), " ", this%force(3,j,i)
        enddo
     enddo

     close(this%bwrtunit)
     close(this%bwrtunit+1)

  endif

!#else

!#endif

end subroutine fce_wrtbasis
  
end module fce_c
