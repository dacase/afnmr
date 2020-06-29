!
! Author:  Leighton W. Wilson (lwwilson@umich.edu)
!          Department of Mathematics
!          University of Michigan, Ann Arbor
!
! The cluster-particle treecode software found here is copyright
! (c) 2013-2019 by the Regents of the University of Michigan.
!
! The 3D-RISM-KH software found here is copyright (c) 2011-2012 by
! Andriy Kovalenko, Tyler Luchko and David A. Case.
!
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU General Public License in the
! ../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
! Users of the 3D-RISM capability found here are requested to acknowledge
! use of the software in reports and publications.  Such acknowledgement
! should include the following citations:
!
! 1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
! ibid. 112:10391-10417 (2000).
!
! 2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by
! F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.
!
! 3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
! and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010).
!
! 4) I. Omelyan and A. Kovalenko, Mol. Simul. 39:25-48 (2013).

#include "../include/dprec.fh"

!> Module for treecode calculations in 3D-RISM. Used to calculate
!! the Coulomb potential and the real space long-range asymptotics
!! for the TCF and DCF.
!!
!! Additionally contains routines for treecode-accelerated 3D-RISM
!! force calculations, if the interest arises in the future.
      MODULE treecode_procedures

      USE safemem
      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(15)
      _REAL_,PARAMETER :: pi = 4_r8*ATAN(1.0_r8)

! global variables for taylor expansions

      INTEGER :: torder, torderlim, torderflat
      _REAL_,ALLOCATABLE,DIMENSION(:) :: cf,cf1,cf2,cf3
      _REAL_,ALLOCATABLE,DIMENSION(:,:,:) :: a, b
      
! global variables used when computing potential/force

      _REAL_,DIMENSION(3) :: tarpos
      _REAL_ :: thetasq,tarposq

! global variables for postition and charge storage

      _REAL_,DIMENSION(3) :: xyz_ddglob
      INTEGER,DIMENSION(3) :: xyz_dimglob

! node pointer and node type declarations

      TYPE tnode_pointer
           TYPE(tnode), POINTER :: p_to_tnode
      END TYPE tnode_pointer
      TYPE tnode
           INTEGER          :: numpar
           _REAL_,DIMENSION(3) :: xyz_min, xyz_max, xyz_mid

           _REAL_           :: radius, sqradius, aspect, min_len
           INTEGER          :: level, num_children, exist_ms
           _REAL_,POINTER :: ms(:) => NULL()
           TYPE(tnode_pointer) :: child(8)

           INTEGER,DIMENSION(3) :: xyz_dim, xyz_lowindex, xyz_highindex
      END TYPE tnode

      CONTAINS


!!!!!!!!!!!!!!!


      SUBROUTINE SETUP(xyzminmax,xyzdim,xyzind,order,theta)
      IMPLICIT NONE
!
! SETUP allocates and initializes arrays needed for the Taylor expansion.
! Also, global variables are set and the Cartesian coordinates of
! the smallest box containing the particles is determined.
!
      INTEGER,INTENT(IN) :: order
      _REAL_,DIMENSION(6),INTENT(IN) :: xyzminmax
      INTEGER,DIMENSION(3),INTENT(IN) :: xyzdim
      INTEGER,DIMENSION(6),INTENT(INOUT) :: xyzind
      _REAL_,INTENT(IN) :: theta

! local variables

      INTEGER :: err,i
      _REAL_ :: t1

! global integers and reals:  TORDER, TORDERLIM and THETASQ

      torder=order
      torderlim=torder+1
      torderflat=torderlim*(torderlim+1)*(torderlim+2)/6
      thetasq=theta*theta

      xyz_dimglob = xyzdim

      xyz_ddglob(1) = (xyzminmax(2)-xyzminmax(1)) / (xyz_dimglob(1)-1)
      xyz_ddglob(2) = (xyzminmax(4)-xyzminmax(3)) / (xyz_dimglob(2)-1)
      xyz_ddglob(3) = (xyzminmax(6)-xyzminmax(5)) / (xyz_dimglob(3)-1)

! allocate global Taylor expansion variables

      ALLOCATE(cf(0:torder), cf1(torderlim), cf2(torderlim), cf3(torderlim), &
               a(-2:torderlim, -2:torderlim, -2:torderlim), &
               b(-2:torderlim, -2:torderlim, -2:torderlim), STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating Taylor variables! '
         STOP
      END IF

! initialize arrays for Taylor sums and coeffs
      DO i = 0, torder
         cf(i) = REAL(i,KIND=r8) + 1.0_r8
      END DO

      DO i = 1, torderlim
         t1 = 1.0_r8 / REAL(i,KIND=r8)
         cf1(i) = t1;
         cf2(i) = 1.0_r8 - 0.5_r8 * t1
         cf3(i) = 1.0_r8 - t1
      END DO

! find bounds of Cartesian box enclosing the particles

      xyzind(1) = 0
      xyzind(2) = xyz_dimglob(1)-1
      xyzind(3) = 0
      xyzind(4) = xyz_dimglob(2)-1
      xyzind(5) = 0
      xyzind(6) = xyz_dimglob(3)-1

      RETURN
      END SUBROUTINE SETUP


!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE CREATE_TREE_N0(p,maxparnode,xyzmm,xyzdim,xyzind,level)
      IMPLICIT NONE
!
! CREATE_TREE_N0 recursively creates the tree structure. Node P is
! input, which contains particles indexed from IBEG to IEND. After
! the node parameters are set subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
! Real array XYZMM contains the min and max values of the coordinates
! of the particle in P, thus defining the box. The division of a cluster terminates
! when the number of particles in a cluster are is less or equal to maxparnode
!
      TYPE(tnode),POINTER :: p
      INTEGER,INTENT(IN) :: level,maxparnode
      _REAL_,DIMENSION(6),INTENT(IN) :: xyzmm
      INTEGER,DIMENSION(3),INTENT(IN) :: xyzdim
      INTEGER,DIMENSION(6),INTENT(IN) :: xyzind

! local variables
      _REAL_, DIMENSION(3) :: xyz_len
      _REAL_ :: lmax
      INTEGER :: i, err, loclev, numposchild

      _REAL_, DIMENSION(6,8) :: xyzmms
      INTEGER, DIMENSION(3,8) :: xyzdims
      INTEGER, DIMENSION(6,8) :: xyzinds

      _REAL_, DIMENSION(6) ::  lxyzmm
      INTEGER, DIMENSION(3) :: lxyzdim
      INTEGER, DIMENSION(6) :: lxyzind

      xyzmms = 0.0_r8
      xyzdims = 0
      xyzinds = 0
      lxyzmm = 0.0_r8
      lxyzdim = 0
      lxyzind = 0
     
! allocate pointer 

      ALLOCATE(p,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating pointer! '
         STOP
      END IF

! set node fields: number of particles, exist_ms
! and xyz bounds 

      p%numpar=product(xyzdim)
      p%exist_ms=0

      p%xyz_min=xyzmm(1:5:2)
      p%xyz_max=xyzmm(2:6:2)

      p%xyz_dim=xyzdim

      p%xyz_lowindex=xyzind(1:5:2)
      p%xyz_highindex=xyzind(2:6:2)

! compute aspect ratio

      xyz_len=p%xyz_max-p%xyz_min

      lmax=MAXVAL(xyz_len)
      p%min_len=MINVAL(xyz_len)

      IF (p%min_len .GT. 0.0_r8) THEN
         p%aspect=lmax/p%min_len
      ELSE
         p%aspect=0.0_r8
      END IF

! midpoint coordinates , RADIUS and SQRADIUS 

      p%xyz_mid=(p%xyz_max+p%xyz_min)/2.0_r8
      p%sqradius=SUM(xyz_len**2)/4.0_r8
      p%radius=SQRT(p%sqradius)

! set particle limits, tree level of node, and nullify children pointers

      p%level=level
      p%num_children=0

      DO i=1,8
         NULLIFY(p%child(i)%p_to_tnode)
      END DO

      IF (p%numpar .GT. maxparnode) THEN
!
! set IND array to 0 and then call PARTITION routine.  IND array holds indices
! of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
!
         xyzmms(:,1)=xyzmm
         xyzdims(:,1)=xyzdim
         xyzinds(:,1)=xyzind

         CALL PARTITION_8(xyzmms,xyzdims,xyzinds,xyz_len,lmax,numposchild)
!
! create children if indicated and store info in parent
!
         loclev=level+1

         DO i=1,numposchild
            IF (((xyzinds(1,i) .LE. xyzinds(2,i)) .AND. &
                 (xyzinds(3,i) .LE. xyzinds(4,i))) .AND. &
                 (xyzinds(5,i) .LE. xyzinds(6,i))) THEN

               p%num_children=p%num_children+1

               lxyzmm=xyzmms(:,i)
               lxyzdim=xyzdims(:,i)
               lxyzind=xyzinds(:,i)

               CALL CREATE_TREE_N0(p%child(p%num_children)%p_to_tnode, &
                                   maxparnode,lxyzmm,lxyzdim,lxyzind,loclev)
            END IF
         END DO

      END IF   

      END SUBROUTINE CREATE_TREE_N0


!!!!!!!!!!!!!!!


      SUBROUTINE PARTITION_8(xyzmms,xyzdims,xyzinds,xyz_len,lmax,numposchild)

      IMPLICIT NONE
!
! PARTITION_8 determines the particle indices of the eight sub boxes
! containing the particles after the box defined by particles I_BEG
! to I_END is divided by its midpoints in each coordinate direction.
! The determination of the indices is accomplished by the subroutine
! PARTITION. A box is divided in a coordinate direction as long as the
! resulting aspect ratio is not too large. This avoids the creation of
! "narrow" boxes in which Talyor expansions may become inefficient.
! On exit the INTEGER array IND (dimension 8 x 2) contains
! the indice limits of each new box (node) and NUMPOSCHILD the number 
! of possible children.  If IND(J,1) > IND(J,2) for a given J this indicates
! that box J is empty.
!
      _REAL_,DIMENSION(6,8),INTENT(INOUT) :: xyzmms
      INTEGER,DIMENSION(3,8),INTENT(INOUT) :: xyzdims
      INTEGER,DIMENSION(6,8),INTENT(INOUT) :: xyzinds

      _REAL_,DIMENSION(3),INTENT(IN) :: xyz_len
      _REAL_,INTENT(IN) :: lmax
      INTEGER,INTENT(INOUT) :: numposchild

! local variables

      INTEGER :: i
      _REAL_ :: critlen
      INTEGER :: xdim,ydim,zdim,xn,yn,zn
      INTEGER :: xlowind,xhighind,ylowind,yhighind,zlowind,zhighind
      _REAL_ :: xlowmid,xhighmid,ylowmid,yhighmid,zlowmid,zhighmid


      numposchild=1
      critlen=lmax/sqrt(2.0_r8)

      xdim=xyzdims(1,1)
      ydim=xyzdims(2,1)
      zdim=xyzdims(3,1)

      xn=xdim/2
      yn=ydim/2
      zn=zdim/2

      IF (xyz_len(1) .GE. critlen) THEN

         xlowmid=xyzmms(1,1)+(xn-1)*xyz_ddglob(1)
         xhighmid=xyzmms(2,1)-(xdim-xn-1)*xyz_ddglob(1)

         xlowind=xyzinds(1,1)+(xn-1)
         xhighind=xyzinds(2,1)-(xdim-xn-1)

         xyzmms(:,2)=xyzmms(:,1)
         xyzinds(:,2)=xyzinds(:,1)

         xyzmms(2,1)=xlowmid
         xyzmms(1,2)=xhighmid

         xyzinds(2,1)=xlowind
         xyzinds(1,2)=xhighind

         numposchild=2*numposchild

      END IF 

      IF (xyz_len(2) .GE. critlen) THEN

         ylowmid=xyzmms(3,1)+(yn-1)*xyz_ddglob(2)
         yhighmid=xyzmms(4,1)-(ydim-yn-1)*xyz_ddglob(2)

         ylowind=xyzinds(3,1)+(yn-1)
         yhighind=xyzinds(4,1)-(ydim-yn-1)

         DO i=1,numposchild
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzinds(:,numposchild+i)=xyzinds(:,i)

            xyzmms(4,i)=ylowmid
            xyzmms(3,numposchild+i)=yhighmid

            xyzinds(4,i)=ylowind
            xyzinds(3,numposchild+i)=yhighind
         END DO

         numposchild=2*numposchild

      END IF

      IF (xyz_len(3) .GE. critlen) THEN

         zlowmid=xyzmms(5,1)+(zn-1)*xyz_ddglob(3)
         zhighmid=xyzmms(6,1)-(zdim-zn-1)*xyz_ddglob(3)

         zlowind=xyzinds(5,1)+(zn-1)
         zhighind=xyzinds(6,1)-(zdim-zn-1)

         DO i=1,numposchild

            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzinds(:,numposchild+i)=xyzinds(:,i)

            xyzmms(6,i)=zlowmid
            xyzmms(5,numposchild+i)=zhighmid

            xyzinds(6,i)=zlowind
            xyzinds(5,numposchild+i)=zhighind

         END DO

         numposchild=2*numposchild

      END IF

      xyzdims(1,:)=xyzinds(2,:)-xyzinds(1,:)+1
      xyzdims(2,:)=xyzinds(4,:)-xyzinds(3,:)+1
      xyzdims(3,:)=xyzinds(6,:)-xyzinds(5,:)+1

      RETURN 
      END SUBROUTINE PARTITION_8


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE CP_TREECODE_TCF(p, solutePosition, soluteCharge, EnP, &
                                 numparsS, numparsT, kappa, eta, eps)
      IMPLICIT NONE
!
! CP_TREECODE is the driver routine which calls COMPUTE_CP1 for each
! source particle, setting the global variable TARPOS before the call. After
! the calls to COMPUTE_CP1, CP_TREECODE calls COMPUTE_CP2 to compute
! the energy using power series. 
!
 
      INTEGER,INTENT(IN) :: numparsS,numparsT
      TYPE(tnode),POINTER :: p  
      _REAL_,DIMENSION(3,numparsS),INTENT(IN) :: solutePosition
      _REAL_,DIMENSION(numparsS),INTENT(IN) :: soluteCharge
      _REAL_,INTENT(IN) :: kappa, eta, eps
      _REAL_,DIMENSION(numparsT),INTENT(INOUT) :: EnP
 
! local variables

      INTEGER :: i,j

      EnP = 0.0_r8
      DO i = 1, numparsS
         tarpos(1) = solutePosition(1, i)
         tarpos(2) = solutePosition(2, i)
         tarpos(3) = solutePosition(3, i)
         tarposq = soluteCharge(i)

         DO j = 1, p%num_children
            CALL COMPUTE_CP1_TCF(p%child(j)%p_to_tnode, EnP, &
                                 numparsT, kappa, eta)
         END DO
      END DO

      CALL COMPUTE_CP2(p, EnP, numparsT)
 
      EnP = EnP * exp((kappa*eta)**2 / 4_r8) / (2_r8*eps)

      RETURN
      END SUBROUTINE CP_TREECODE_TCF


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE CP_TREECODE_DCF(p, solutePosition, soluteCharge, EnP, &
                                 numparsS, numparsT, eta)
      IMPLICIT NONE
!
! CP_TREECODE is the driver routine which calls COMPUTE_CP1 for each
! source particle, setting the global variable TARPOS before the call. After
! the calls to COMPUTE_CP1, CP_TREECODE calls COMPUTE_CP2 to compute
! the energy using power series. 
!
 
      INTEGER,INTENT(IN) :: numparsS,numparsT
      TYPE(tnode),POINTER :: p  
      _REAL_,DIMENSION(3,numparsS),INTENT(IN) :: solutePosition
      _REAL_,DIMENSION(numparsS),INTENT(IN) :: soluteCharge
      _REAL_,INTENT(IN) :: eta
      _REAL_,DIMENSION(numparsT),INTENT(INOUT) :: EnP
 
! local variables

      INTEGER :: i,j

      EnP = 0.0_r8
      DO i = 1, numparsS
         tarpos(1) = solutePosition(1, i)
         tarpos(2) = solutePosition(2, i)
         tarpos(3) = solutePosition(3, i)
         tarposq = soluteCharge(i)

         DO j = 1, p%num_children
            CALL COMPUTE_CP1_DCF(p%child(j)%p_to_tnode, EnP, &
                                 numparsT, eta)
         END DO
      END DO

      CALL COMPUTE_CP2(p,EnP,numparsT)

      RETURN
      END SUBROUTINE CP_TREECODE_DCF


!!!!!!!!!!!!!!


      SUBROUTINE CP_TREECODE_COULOMB(p, solutePosition, soluteCharge, EnP, &
                                     numparsS, numparsT)
      IMPLICIT NONE
!
! CP_TREECODE is the driver routine which calls COMPUTE_CP1 for each
! source particle, setting the global variable TARPOS before the call. After
! the calls to COMPUTE_CP1, CP_TREECODE calls COMPUTE_CP2 to compute
! the energy using power series. 
!
 
      INTEGER,INTENT(IN) :: numparsS,numparsT
      TYPE(tnode),POINTER :: p  
      _REAL_,DIMENSION(3,numparsS),INTENT(IN) :: solutePosition
      _REAL_,DIMENSION(numparsS),INTENT(IN) :: soluteCharge
      _REAL_,DIMENSION(numparsT),INTENT(INOUT) :: EnP
 
! local variables

      INTEGER :: i,j

      EnP = 0.0_r8
      DO i = 1, numparsS
         tarpos(1) = solutePosition(1, i)
         tarpos(2) = solutePosition(2, i)
         tarpos(3) = solutePosition(3, i)
         tarposq = soluteCharge(i)

         DO j = 1, p%num_children
            CALL COMPUTE_CP1_COULOMB(p%child(j)%p_to_tnode,EnP, &
                                     numparsT)
         END DO
      END DO

      CALL COMPUTE_CP2(p,EnP,numparsT)
 
      !EnP = EnP / (2_r8*eps)

      RETURN
      END SUBROUTINE CP_TREECODE_COULOMB


!!!!!!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE COMPUTE_CP1_TCF(p, EnP, arrdim, kappa, eta)

      IMPLICIT NONE
!
! COMPUTE_CP1 is the recursive routine for computing the interaction
! between a target particle and a source cluster. If the MAC is
! satisfied the power series coefficients for the current target 
! cluster are updated. If the MAC is not satisfied then the algorithm 
! descends to the children of the current cluster, unless the
! current cluster is a leaf then the interaction is done exactly
! via a call to the routine COMP_DIRECT
!

      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p      
      _REAL_,DIMENSION(arrdim),INTENT(INOUT) :: EnP
      _REAL_,INTENT(IN) :: kappa, eta

! local variables

      _REAL_,DIMENSION(3) :: xyz_t
      _REAL_ :: distsq, dist
      INTEGER :: i, err

! determine DISTSQ for MAC test

      xyz_t=tarpos-p%xyz_mid
      distsq=SUM(xyz_t**2)
      dist=SQRT(distsq)

! intialize potential energy and force 

! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.
      IF ((    p%sqradius .LT. distsq * thetasq) .AND. &
          (    p%sqradius .GT. 0.0_r8          ) .AND. &
          (dist-p%min_len .GE. 5.0_r8 * eta    ) .AND. &
          (    torderflat .LE. 4 * p%numpar    )) THEN

         CALL COMP_TCOEFF_TCF(xyz_t(1),xyz_t(2),xyz_t(3), kappa)

         IF (p%exist_ms .EQ. 0) THEN
             p%ms => safemem_realloc(p%ms, torderflat, .false.)
             p%ms=0.0_r8
             p%exist_ms=1
         END IF

         CALL COMP_CMS(p)   
    
      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for each.
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT_TCF(p, EnP, arrdim, kappa, eta)
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_CP1_TCF(p%child(i)%p_to_tnode, EnP, arrdim, kappa, eta)
            END DO  
         END IF 
      END IF

      RETURN
      END SUBROUTINE COMPUTE_CP1_TCF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE COMPUTE_CP1_DCF(p,EnP,arrdim, eta)

      IMPLICIT NONE
!
! COMPUTE_CP1 is the recursive routine for computing the interaction
! between a target particle and a source cluster. If the MAC is
! satisfied the power series coefficients for the current target 
! cluster are updated. If the MAC is not satisfied then the algorithm 
! descends to the children of the current cluster, unless the
! current cluster is a leaf then the interaction is done exactly
! via a call to the routine COMP_DIRECT
!

      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p      
      _REAL_,DIMENSION(arrdim),INTENT(INOUT) :: EnP
      _REAL_,INTENT(IN) :: eta

! local variables

      _REAL_,DIMENSION(3) :: xyz_t
      _REAL_ :: distsq, dist
      INTEGER :: i, err

! determine DISTSQ for MAC test

      xyz_t=tarpos-p%xyz_mid
      distsq=SUM(xyz_t**2)
      dist=SQRT(distsq)

! intialize potential energy and force 

! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.
      IF ((    p%sqradius .LT. distsq * thetasq) .AND. &
          (    p%sqradius .GT. 0.0_r8          ) .AND. &
          (dist-p%min_len .GE. 5.0_r8 * eta    ) .AND. &
          (    torderflat .LE. p%numpar        )) THEN

         CALL COMP_TCOEFF_DCF(xyz_t(1),xyz_t(2),xyz_t(3), eta)

         IF (p%exist_ms .EQ. 0) THEN
             p%ms => safemem_realloc(p%ms, torderflat, .false.)
             p%ms=0.0_r8
             p%exist_ms=1
         END IF

         CALL COMP_CMS(p)   
    
      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for each.
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT_DCF(p, EnP, arrdim, eta)
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_CP1_DCF(p%child(i)%p_to_tnode, EnP, arrdim, eta)
            END DO  
         END IF 
      END IF

      RETURN
      END SUBROUTINE COMPUTE_CP1_DCF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE COMPUTE_CP1_COULOMB(p,EnP,arrdim)

      IMPLICIT NONE
!
! COMPUTE_CP1 is the recursive routine for computing the interaction
! between a target particle and a source cluster. If the MAC is
! satisfied the power series coefficients for the current target 
! cluster are updated. If the MAC is not satisfied then the algorithm 
! descends to the children of the current cluster, unless the
! current cluster is a leaf then the interaction is done exactly
! via a call to the routine COMP_DIRECT
!

      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p      
      _REAL_,DIMENSION(arrdim),INTENT(INOUT) :: EnP

! local variables

      _REAL_,DIMENSION(3) :: xyz_t
      _REAL_ :: distsq
      INTEGER :: i, err

! determine DISTSQ for MAC test

      xyz_t=tarpos-p%xyz_mid
      distsq=SUM(xyz_t**2)

! intialize potential energy and force 

! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.
      IF ((p%sqradius .LT. distsq*thetasq) .AND. &
         (p%sqradius .GT. 0.0_r8)) THEN

         CALL COMP_TCOEFF_COULOMB(xyz_t(1),xyz_t(2),xyz_t(3))

         IF (p%exist_ms .EQ. 0) THEN
             p%ms => safemem_realloc(p%ms, torderflat, .false.)
             p%ms=0.0_r8
             p%exist_ms=1
         END IF

         CALL COMP_CMS(p)   
    
      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for each.
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT_COULOMB(p, EnP, arrdim)
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_CP1_COULOMB(p%child(i)%p_to_tnode, EnP, arrdim)
            END DO  
         END IF 
      END IF

      RETURN
      END SUBROUTINE COMPUTE_CP1_COULOMB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE COMPUTE_CP2(ap,EnP,arrdim)

        IMPLICIT NONE
!
! COMPUTE_CP2 is a recursive routine that evaluates the power series  
! approximation of the potential at the targets in a cluster via 
! a 3-D Horner's rule.  
!
        INTEGER,INTENT(IN) :: arrdim
        TYPE(tnode),POINTER :: ap
        _REAL_,DIMENSION(arrdim),INTENT(INOUT) :: EnP

! local variables

        _REAL_ :: tx,ty,peng,dx,dy,dz,xl,yl,zl,xm,ym,zm
        INTEGER :: xlind, ylind, zlind, xhind, yhind, zhind
        INTEGER :: i,nn,j,k,k1,k2,k3,kk,porder,porder1,nj,nk,xyind

        porder=torder
        porder1=porder-1

        IF (ap%exist_ms==1) THEN

          xl=ap%xyz_min(1)
          yl=ap%xyz_min(2)
          zl=ap%xyz_min(3)

          xm=ap%xyz_mid(1)
          ym=ap%xyz_mid(2)
          zm=ap%xyz_mid(3)

          xlind=ap%xyz_lowindex(1)
          ylind=ap%xyz_lowindex(2)
          zlind=ap%xyz_lowindex(3)

          xhind=ap%xyz_highindex(1)
          yhind=ap%xyz_highindex(2)
          zhind=ap%xyz_highindex(3)

          xyind=xyz_dimglob(1)*xyz_dimglob(2)
          DO k=zlind,zhind
             dz=zl+(k-zlind)*xyz_ddglob(3)-zm
             nk=k*xyind
             DO j=ylind,yhind
                dy=yl+(j-ylind)*xyz_ddglob(2)-ym
                nj=nk+j*xyz_dimglob(1)
                DO i=xlind,xhind
                   dx=xl+(i-xlind)*xyz_ddglob(1)-xm
                   nn=nj+i+1

                   kk=1
                   peng=ap%ms(kk)
                   DO k3=porder1,0,-1
                      kk=kk+1
                      ty=ap%ms(kk)
                      DO k2=porder1-k3,0,-1
                         kk=kk+1
                         tx=ap%ms(kk)
                         DO k1=porder1-k3-k2,0,-1
                            kk=kk+1
                            tx  = dx * tx + ap%ms(kk)
                         END DO
                         ty  = dy * ty + tx
                      END DO
                      peng = dz * peng + ty
                   END DO
                   EnP(nn)=EnP(nn)+peng
                END DO
             END DO
          END DO

        END IF

        DO j=1,ap%num_children
           CALL COMPUTE_CP2(ap%child(j)%p_to_tnode,EnP,arrdim)
        END DO

        RETURN
      END SUBROUTINE COMPUTE_CP2


!!!!!!!!!!!!!!!


      SUBROUTINE COMP_TCOEFF_TCF(dx,dy,dz,kappa)
      IMPLICIT NONE
!
! COMP_TCOEFF computes the Taylor coefficients of the potential
! using a recurrence formula.  The center of the expansion is the
! midpoint of the node P.  TARPOS and TORDERLIM are globally defined.
!
      _REAL_,INTENT(IN) :: dx,dy,dz
      _REAL_,INTENT(IN)  :: kappa

! local variables

      _REAL_ :: ddx,ddy,ddz,dist,fac
      _REAL_ :: kappax,kappay,kappaz
      INTEGER :: i,j,k,i1,i2,j1,j2,k1,k2

! setup variables

      ddx = 2.0_r8*dx
      ddy = 2.0_r8*dy
      ddz = 2.0_r8*dz

      kappax = kappa*dx
      kappay = kappa*dy
      kappaz = kappa*dz

      dist = dx*dx + dy*dy + dz*dz
      fac = 1.0_r8 / dist
      dist = SQRT(dist)

! 0th coeff or function val

      b(0,0,0)=-2_r8*EXP(-kappa*dist)
      a(0,0,0)=b(0,0,0)/dist

! 2 indices are 0

      b(1,0,0)=kappax*a(0,0,0)
      b(0,1,0)=kappay*a(0,0,0)
      b(0,0,1)=kappaz*a(0,0,0)

      a(1,0,0)=fac*dx*(a(0,0,0)+kappa*b(0,0,0))
      a(0,1,0)=fac*dy*(a(0,0,0)+kappa*b(0,0,0))
      a(0,0,1)=fac*dz*(a(0,0,0)+kappa*b(0,0,0))


      DO i=2,torderlim
         i1=i-1; i2=i-2;

         b(i,0,0)=cf1(i)*kappa*(dx*a(i1,0,0)-a(i2,0,0))
         b(0,i,0)=cf1(i)*kappa*(dy*a(0,i1,0)-a(0,i2,0))
         b(0,0,i)=cf1(i)*kappa*(dz*a(0,0,i1)-a(0,0,i2))

         a(i,0,0)=fac*(ddx*cf2(i)*a(i1,0,0)-cf3(i)*a(i2,0,0)+&
                  cf1(i)*kappa*(dx*b(i1,0,0)-b(i2,0,0)))
         a(0,i,0)=fac*(ddy*cf2(i)*a(0,i1,0)-cf3(i)*a(0,i2,0)+&
                  cf1(i)*kappa*(dy*b(0,i1,0)-b(0,i2,0)))
         a(0,0,i)=fac*(ddz*cf2(i)*a(0,0,i1)-cf3(i)*a(0,0,i2)+&
                  cf1(i)*kappa*(dz*b(0,0,i1)-b(0,0,i2)))
      END DO

! 1 index 0, 1 index 1, other >=1

      b(1,1,0)=kappax*a(0,1,0)
      b(1,0,1)=kappax*a(0,0,1)
      b(0,1,1)=kappay*a(0,0,1)

      a(1,1,0)=fac*(dx*a(0,1,0)+ddy*a(1,0,0)+kappax*b(0,1,0))
      a(1,0,1)=fac*(dx*a(0,0,1)+ddz*a(1,0,0)+kappax*b(0,0,1))
      a(0,1,1)=fac*(dy*a(0,0,1)+ddz*a(0,1,0)+kappay*b(0,0,1))

      DO i=2,torderlim-1
         i1=i-1; i2=i-2;

         b(1,0,i)=kappax*a(0,0,i)
         b(0,1,i)=kappay*a(0,0,i)
         b(0,i,1)=kappaz*a(0,i,0)
         b(1,i,0)=kappax*a(0,i,0)
         b(i,1,0)=kappay*a(i,0,0)
         b(i,0,1)=kappaz*a(i,0,0)

         a(1,0,i)=fac*(dx*a(0,0,i)+ddz*a(1,0,i1)-a(1,0,i2)+&
                  kappax*b(0,0,i))
         a(0,1,i)=fac*(dy*a(0,0,i)+ddz*a(0,1,i1)-a(0,1,i2)+&
                  kappay*b(0,0,i))
         a(0,i,1)=fac*(dz*a(0,i,0)+ddy*a(0,i1,1)-a(0,i2,1)+&
                  kappaz*b(0,i,0))
         a(1,i,0)=fac*(dx*a(0,i,0)+ddy*a(1,i1,0)-a(1,i2,0)+&
                  kappax*b(0,i,0))
         a(i,1,0)=fac*(dy*a(i,0,0)+ddx*a(i1,1,0)-a(i2,1,0)+&
                  kappay*b(i,0,0))
         a(i,0,1)=fac*(dz*a(i,0,0)+ddx*a(i1,0,1)-a(i2,0,1)+&
                  kappaz*b(i,0,0))
      END DO

! 1 index 0, others >= 2

      DO i=2,torderlim-2
         i1=i-1; i2=i-2;

         DO j=2,torderlim-i
            j1=j-1; j2=j-2;

            b(i,j,0)=cf1(i)*kappa*(dx*a(i1,j,0)-a(i2,j,0))
            b(i,0,j)=cf1(i)*kappa*(dx*a(i1,0,j)-a(i2,0,j))
            b(0,i,j)=cf1(i)*kappa*(dy*a(0,i1,j)-a(0,i2,j))

            a(i,j,0)=fac*(ddx*cf2(i)*a(i1,j,0)+ddy*a(i,j1,0) &
                     -cf3(i)*a(i2,j,0)-a(i,j2,0)+&
                     cf1(i)*kappa*(dx*b(i1,j,0)-b(i2,j,0)))
            a(i,0,j)=fac*(ddx*cf2(i)*a(i1,0,j)+ddz*a(i,0,j1)&
                     -cf3(i)*a(i2,0,j)-a(i,0,j2)+&
                     cf1(i)*kappa*(dx*b(i1,0,j)-b(i2,0,j)))
            a(0,i,j)=fac*(ddy*cf2(i)*a(0,i1,j)+ddz*a(0,i,j1)&
                     -cf3(i)*a(0,i2,j)-a(0,i,j2)+&
                     cf1(i)*kappa*(dy*b(0,i1,j)-b(0,i2,j)))
         END DO
      END DO

! 2 indices 1, other >= 1
! b(1,1,1) is correct, but a little tricky!
!      b(1,1,1)=5.0*dz*fac*b(1,1,0)

      b(1,1,1)=kappax*a(0,1,1)
      a(1,1,1)=fac*(dx*a(0,1,1)+ddy*a(1,0,1)+ddz*a(1,1,0)+&
               kappax*b(0,1,1))

      DO i=2,torderlim-2
         i1=i-1; i2=i-2;

         b(1,1,i)=kappax*a(0,1,i)
         b(1,i,1)=kappax*a(0,i,1)
         b(i,1,1)=kappay*a(i,0,1)

         a(1,1,i)=fac*(dx*a(0,1,i)+ddy*a(1,0,i)+ddz*a(1,1,i1)&
                 -a(1,1,i2)+kappax*b(0,1,i))
         a(1,i,1)=fac*(dx*a(0,i,1)+ddy*a(1,i1,1)+ddz*a(1,i,0)&
                 -a(1,i2,1)+kappax*b(0,i,1))
         a(i,1,1)=fac*(dy*a(i,0,1)+ddx*a(i1,1,1)+ddz*a(i,1,0)&
                 -a(i2,1,1)+kappay*b(i,0,1))
      END DO

! 1 index 1, others >=2

      DO i=2,torderlim-3
         i1=i-1; i2=i-2;

         DO j=2,torderlim-1-i
            j1=j-1; j2=j-2;

            b(1,i,j)=kappax*a(0,i,j)
            b(i,1,j)=kappay*a(i,0,j)
            b(i,j,1)=kappaz*a(i,j,0)

            a(1,i,j)=fac*(dx*a(0,i,j)+ddy*a(1,i1,j)+ddz*a(1,i,j1)&
                    -a(1,i2,j)-a(1,i,j2)+kappax*b(0,i,j))
            a(i,1,j)=fac*(dy*a(i,0,j)+ddx*a(i1,1,j)+ddz*a(i,1,j1)&
                    -a(i2,1,j)-a(i,1,j2)+kappay*b(i,0,j))
            a(i,j,1)=fac*(dz*a(i,j,0)+ddx*a(i1,j,1)+ddy*a(i,j1,1)&
                    -a(i2,j,1)-a(i,j2,1)+kappaz*b(i,j,0))

         END DO
      END DO

! all indices >=2

      DO k=2,torderlim-4
         k1=k-1; k2=k-2;

         DO j=2,torderlim-2-k
            j1=j-1; j2=j-2;

            DO i=2,torderlim-k-j
               i1=i-1; i2=i-2;

               !b(i,j,k)=cf1(i)*kappa*(dx*a(i1,j,k)-a(i2,j,k))
               ! 
               !a(i,j,k)=fac*(ddx*cf2(i)*a(i1,j,k)+ddy*a(i,j1,k)&
               !        +ddz*a(i,j,k1)-cf3(i)*a(i2,j,k)&
               !        -a(i,j2,k)-a(i,j,k2)+&
               !        cf1(i)*kappa*(dx*b(i1,j,k)-b(i2,j,k)))

               b(i,j,k) = cf1(i+j+k) * kappa &
                        * (dx*a(i1,j,k) + dy*a(i,j1,k) + dz*a(i,j,k1) &
                            - a(i2,j,k)    - a(i,j2,k)    - a(i,j,k2))

               a(i,j,k) = fac &
                        * (cf2(i+j+k) * (ddx*a(i1,j,k) + ddy*a(i,j1,k) + ddz*a(i,j,k1)) &
                         - cf3(i+j+k) * (    a(i2,j,k) +     a(i,j2,k) +     a(i,j,k2)) &
                         + cf1(i+j+k) * kappa &
                                       * (dx*b(i1,j,k) +  dy*b(i,j1,k) +  dz*b(i,j,k1) &
                                           - b(i2,j,k)     - b(i,j2,k)     - b(i,j,k2)))
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE COMP_TCOEFF_TCF


!!!!!!!!!!!!!!!


      SUBROUTINE COMP_TCOEFF_DCF(dx,dy,dz,eta)
      IMPLICIT NONE
!
! COMP_TCOEFF computes the Taylor coefficients of the potential
! using a recurrence formula.  The center of the expansion is the
! midpoint of the node P.  TARPOS and TORDERLIM are globally defined.
!
      _REAL_,INTENT(IN) :: dx,dy,dz
      _REAL_,INTENT(IN)  :: eta

! local varaibles

      _REAL_ :: ddx,ddy,ddz,dist,fac
      _REAL_ :: two_etasq, two_etasqx
      _REAL_ :: two_etasqy, two_etasqz
      INTEGER :: i,j,k,i1,i2,j1,j2,k1,k2

! setup variables

      ddx = 2.0_r8*dx
      ddy = 2.0_r8*dy
      ddz = 2.0_r8*dz

      two_etasq = 2.0_r8 / (eta*eta)
      two_etasqx = two_etasq*dx
      two_etasqy = two_etasq*dy
      two_etasqz = two_etasq*dz

      dist = dx*dx + dy*dy + dz*dz
      fac = 1.0_r8 / dist
      dist = SQRT(dist)

! 0th coeff or function val

      b(0,0,0) = -eta/sqrt(pi) * exp(-dist*dist/(eta*eta))
      a(0,0,0) = -erf(dist/eta)/dist

! 2 indices are 0

      b(1,0,0) = two_etasqx * b(0,0,0)
      b(0,1,0) = two_etasqy * b(0,0,0)
      b(0,0,1) = two_etasqz * b(0,0,0)

      a(1,0,0) = fac*dx*(a(0,0,0) + b(1,0,0))
      a(0,1,0) = fac*dy*(a(0,0,0) + b(0,1,0))
      a(0,0,1) = fac*dz*(a(0,0,0) + b(0,0,1))


      DO i=2,torderlim
         i1=i-1; i2=i-2;

         b(i,0,0) = cf1(i)*two_etasq * (dx*b(i1,0,0) - b(i2,0,0))
         b(0,i,0) = cf1(i)*two_etasq * (dy*b(0,i1,0) - b(0,i2,0))
         b(0,0,i) = cf1(i)*two_etasq * (dz*b(0,0,i1) - b(0,0,i2))

         a(i,0,0) = fac*(ddx*cf2(i)*a(i1,0,0)-cf3(i)*a(i2,0,0)+&
                    b(i,0,0))
         a(0,i,0) = fac*(ddy*cf2(i)*a(0,i1,0)-cf3(i)*a(0,i2,0)+&
                    b(0,i,0))
         a(0,0,i) = fac*(ddz*cf2(i)*a(0,0,i1)-cf3(i)*a(0,0,i2)+&
                    b(0,0,i))
      END DO

! 1 index 0, 1 index 1, other >=1

      b(1,1,0) = two_etasqx * b(0,1,0)
      b(1,0,1) = two_etasqx * b(0,0,1)
      b(0,1,1) = two_etasqy * b(0,0,1)

      a(1,1,0) = fac*(dx*a(0,1,0)+ddy*a(1,0,0) + b(1,1,0))
      a(1,0,1) = fac*(dx*a(0,0,1)+ddz*a(1,0,0) + b(1,0,1))
      a(0,1,1) = fac*(dy*a(0,0,1)+ddz*a(0,1,0) + b(0,1,1))

      DO i=2,torderlim-1
         i1=i-1; i2=i-2;

         b(1,0,i) = two_etasqx * b(0,0,i)
         b(0,1,i) = two_etasqy * b(0,0,i)
         b(0,i,1) = two_etasqz * b(0,i,0)
         b(1,i,0) = two_etasqx * b(0,i,0)
         b(i,1,0) = two_etasqy * b(i,0,0)
         b(i,0,1) = two_etasqz * b(i,0,0)

         a(1,0,i) = fac*(dx*a(0,0,i)+ddz*a(1,0,i1)-a(1,0,i2)+&
                    b(1,0,i))
         a(0,1,i) = fac*(dy*a(0,0,i)+ddz*a(0,1,i1)-a(0,1,i2)+&
                    b(0,1,i))
         a(0,i,1) = fac*(dz*a(0,i,0)+ddy*a(0,i1,1)-a(0,i2,1)+&
                    b(0,i,1))
         a(1,i,0) = fac*(dx*a(0,i,0)+ddy*a(1,i1,0)-a(1,i2,0)+&
                    b(1,i,0))
         a(i,1,0) = fac*(dy*a(i,0,0)+ddx*a(i1,1,0)-a(i2,1,0)+&
                    b(i,1,0))
         a(i,0,1) = fac*(dz*a(i,0,0)+ddx*a(i1,0,1)-a(i2,0,1)+&
                    b(i,0,1))
      END DO

! 1 index 0, others >= 2

      DO i=2,torderlim-2
         i1=i-1; i2=i-2;

         DO j=2,torderlim-i
            j1=j-1; j2=j-2;

            b(i,j,0) = cf1(i)*two_etasq * (dx*b(i1,j,0) - b(i2,j,0))
            b(i,0,j) = cf1(i)*two_etasq * (dx*b(i1,0,j) - b(i2,0,j))
            b(0,i,j) = cf1(i)*two_etasq * (dx*b(0,i1,j) - b(0,i2,j))

            a(i,j,0) = fac*(ddx*cf2(i)*a(i1,j,0)+ddy*a(i,j1,0) &
                       -cf3(i)*a(i2,j,0)-a(i,j2,0)+&
                       b(i,j,0))
            a(i,0,j) = fac*(ddx*cf2(i)*a(i1,0,j)+ddz*a(i,0,j1)&
                       -cf3(i)*a(i2,0,j)-a(i,0,j2)+&
                       b(i,0,j))
            a(0,i,j) = fac*(ddy*cf2(i)*a(0,i1,j)+ddz*a(0,i,j1)&
                       -cf3(i)*a(0,i2,j)-a(0,i,j2)+&
                       b(0,i,j))
         END DO
      END DO

! 2 indices 1, other >= 1
! b(1,1,1) is correct, but a little tricky!
!      b(1,1,1)=5.0*dz*fac*b(1,1,0)

      b(1,1,1) = two_etasqx * b(0,1,1)
      a(1,1,1) = fac*(dx*a(0,1,1)+ddy*a(1,0,1)+ddz*a(1,1,0)+&
                 b(1,1,1))

      DO i=2,torderlim-2
         i1=i-1; i2=i-2;

         b(1,1,i) = two_etasqx * b(0,1,i)
         b(1,i,1) = two_etasqx * b(0,i,1)
         b(i,1,1) = two_etasqy * b(i,0,1)

         a(1,1,i)=fac*(dx*a(0,1,i)+ddy*a(1,0,i)+ddz*a(1,1,i1)&
                 -a(1,1,i2) + b(1,1,i))
         a(1,i,1)=fac*(dx*a(0,i,1)+ddy*a(1,i1,1)+ddz*a(1,i,0)&
                 -a(1,i2,1) + b(1,i,1))
         a(i,1,1)=fac*(dy*a(i,0,1)+ddx*a(i1,1,1)+ddz*a(i,1,0)&
                 -a(i2,1,1) + b(i,1,1))
      END DO

! 1 index 1, others >=2

      DO i=2,torderlim-3
         i1=i-1; i2=i-2;

         DO j=2,torderlim-1-i
            j1=j-1; j2=j-2;

            b(1,i,j) = two_etasqx * b(0,i,j)
            b(i,1,j) = two_etasqy * b(i,0,j)
            b(i,j,1) = two_etasqz * b(i,j,0)

            a(1,i,j) = fac*(dx*a(0,i,j)+ddy*a(1,i1,j)+ddz*a(1,i,j1)&
                      -a(1,i2,j)-a(1,i,j2) + b(1,i,j))
            a(i,1,j) = fac*(dy*a(i,0,j)+ddx*a(i1,1,j)+ddz*a(i,1,j1)&
                      -a(i2,1,j)-a(i,1,j2) + b(i,1,j))
            a(i,j,1) = fac*(dz*a(i,j,0)+ddx*a(i1,j,1)+ddy*a(i,j1,1)&
                      -a(i2,j,1)-a(i,j2,1) + b(i,j,1))

         END DO
      END DO

! all indices >=2

      DO k=2,torderlim-4
         k1=k-1; k2=k-2;

         DO j=2,torderlim-2-k
            j1=j-1; j2=j-2;

            DO i=2,torderlim-k-j
               i1=i-1; i2=i-2;

               !b(i,j,k) = cf1(i)*two_etasq * (dx*b(i1,j,k) - b(i2,j,k))
               !
               !a(i,j,k) = fac * (ddx*cf2(i)*a(i1,j,k) + ddy*a(i,j1,k) &
               !                    + ddz*a(i,j,k1) - cf3(i)*a(i2,j,k) &
               !                    - a(i,j2,k) - a(i,j,k2) + b(i,j,k))

               b(i,j,k) = cf1(i+j+k) * two_etasq &
                        * (dx*b(i1,j,k) + dy*b(i,j1,k) + dz*b(i,j,k1) &
                            - b(i2,j,k)    - b(i,j2,k)    - b(i,j,k2))

               a(i,j,k) = fac &
                        * (cf2(i+j+k) * (ddx*a(i1,j,k) + ddy*a(i,j1,k) + ddz*a(i,j,k1)) &
                         - cf3(i+j+k) * (    a(i2,j,k) +     a(i,j2,k) +     a(i,j,k2)) &
                           + b(i,j,k))
            END DO
         END DO
      END DO

      RETURN

      END SUBROUTINE COMP_TCOEFF_DCF


!!!!!!!!!!!!!!!


      SUBROUTINE COMP_TCOEFF_COULOMB(dx,dy,dz)
      IMPLICIT NONE
!
! COMP_TCOEFF computes the Taylor coefficients of the potential
! using a recurrence formula.  The center of the expansion is the
! midpoint of the node P.  TARPOS and TORDERLIM are globally defined.
!
      _REAL_,INTENT(IN) :: dx,dy,dz

! local variables

      _REAL_ :: ddx,ddy,ddz,fac,sqfac
      INTEGER :: i,j,k,i1,i2,j1,j2,k1,k2

! setup variables

      ddx=2.0_r8*dx
      ddy=2.0_r8*dy
      ddz=2.0_r8*dz
      fac=1.0_r8/(dx*dx+dy*dy+dz*dz)
      sqfac=SQRT(fac)

! 0th coeff or function val 

      a(0,0,0)=-sqfac

! 2 indices are 0

      a(1,0,0)=fac*dx*a(0,0,0)
      a(0,1,0)=fac*dy*a(0,0,0)
      a(0,0,1)=fac*dz*a(0,0,0)
  
      DO i=2,torderlim
         i1=i-1; i2=i-2
         a(i,0,0)=fac*(ddx*cf2(i)*a(i1,0,0)-cf3(i)*a(i2,0,0))
         a(0,i,0)=fac*(ddy*cf2(i)*a(0,i1,0)-cf3(i)*a(0,i2,0))
         a(0,0,i)=fac*(ddz*cf2(i)*a(0,0,i1)-cf3(i)*a(0,0,i2))
      END DO

! 1 index 0, 1 index 1, other >=1

      a(1,1,0)=fac*(dx*a(0,1,0)+ddy*a(1,0,0))
      a(1,0,1)=fac*(dx*a(0,0,1)+ddz*a(1,0,0))
      a(0,1,1)=fac*(dy*a(0,0,1)+ddz*a(0,1,0))

      DO i=2,torderlim-1
         i1=i-1; i2=i-2
         a(1,0,i)=fac*(dx*a(0,0,i)+ddz*a(1,0,i1)-a(1,0,i2))
         a(0,1,i)=fac*(dy*a(0,0,i)+ddz*a(0,1,i1)-a(0,1,i2))
         a(0,i,1)=fac*(dz*a(0,i,0)+ddy*a(0,i1,1)-a(0,i2,1))
         a(1,i,0)=fac*(dx*a(0,i,0)+ddy*a(1,i1,0)-a(1,i2,0))
         a(i,1,0)=fac*(dy*a(i,0,0)+ddx*a(i1,1,0)-a(i2,1,0))
         a(i,0,1)=fac*(dz*a(i,0,0)+ddx*a(i1,0,1)-a(i2,0,1))
      END DO
  
! 1 index 0, others >= 2
      
      DO i=2,torderlim-2 
         i1=i-1; i2=i-2
         DO j=2,torderlim-i
            j1=j-1; j2=j-2
            a(i,j,0)=fac*(ddx*cf2(i)*a(i1,j,0)+ddy*a(i,j1,0) &
                            -cf3(i)*a(i2,j,0)-a(i,j2,0))
            a(i,0,j)=fac*(ddx*cf2(i)*a(i1,0,j)+ddz*a(i,0,j1) &
                            -cf3(i)*a(i2,0,j)-a(i,0,j2))
            a(0,i,j)=fac*(ddy*cf2(i)*a(0,i1,j)+ddz*a(0,i,j1) &
                            -cf3(i)*a(0,i2,j)-a(0,i,j2))
         END DO
      END DO

 
! 2 indices 1, other >= 1

      a(1,1,1)=fac*(dx*a(0,1,1)+ddy*a(1,0,1)+ddz*a(1,1,0))

      DO i=2,torderlim-2
         i1=i-1; i2=i-2
         a(1,1,i)=fac*(dx*a(0,1,i)+ddy*a(1,0,i)+ddz*a(1,1,i1) &
                        -a(1,1,i2))
         a(1,i,1)=fac*(dx*a(0,i,1)+ddy*a(1,i1,1)+ddz*a(1,i,0) &
                        -a(1,i2,1))
         a(i,1,1)=fac*(dy*a(i,0,1)+ddx*a(i1,1,1)+ddz*a(i,1,0) &
                        -a(i2,1,1))
      END DO

! 1 index 1, others >=2
      DO i=2,torderlim-3
         i1=i-1; i2=i-2 
         DO j=2,torderlim-1-i
            j1=j-1; j2=j-2
            a(1,i,j)=fac*(dx*a(0,i,j)+ddy*a(1,i1,j)+ddz*a(1,i,j1) &
                           -a(1,i2,j)-a(1,i,j2))
            a(i,1,j)=fac*(dy*a(i,0,j)+ddx*a(i1,1,j)+ddz*a(i,1,j1) &
                           -a(i2,1,j)-a(i,1,j2))
            a(i,j,1)=fac*(dz*a(i,j,0)+ddx*a(i1,j,1)+ddy*a(i,j1,1) &
                           -a(i2,j,1)-a(i,j2,1))
         END DO
      END DO

! all indices >=2
      DO k=2,torderlim-4
         k1=k-1; k2=k-2
         DO j=2,torderlim-2-k
            j1=j-1; j2=j-2
            DO i=2,torderlim-k-j
               i1=i-1; i2=i-2

               !a(i,j,k)=fac*(ddx*cf2(i)*a(i-1,j,k)+ddy*a(i,j1,k) &
               !            +ddz*a(i,j,k1)-cf3(i)*a(i-2,j,k) &
               !            -a(i,j2,k)-a(i,j,k2)) 

               a(i,j,k) = fac &
                        * (cf2(i+j+k) * (ddx*a(i1,j,k) + ddy*a(i,j1,k) + ddz*a(i,j,k1)) &
                         - cf3(i+j+k) * (    a(i2,j,k) +     a(i,j2,k) +     a(i,j,k2)))
            END DO
         END DO
      END DO

      RETURN

      END SUBROUTINE COMP_TCOEFF_COULOMB


!!!!!!!!!!!!!!!


      SUBROUTINE COMP_CMS(p)
      IMPLICIT NONE
!
! COMP_CMS computes the moments for node P needed in the Taylor approximation
!
      TYPE(tnode),POINTER :: p 

! local variables

      INTEGER :: k1,k2,k3,kk

      kk=0
        
      DO k3=torder,0,-1
         DO k2=torder-k3,0,-1
            DO k1=torder-k3-k2,0,-1
               kk=kk+1
               p%ms(kk)=p%ms(kk)+tarposq*a(k1,k2,k3)
           END DO
         END DO
      END DO
         
      RETURN
      END SUBROUTINE COMP_CMS


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_DIRECT_TCF(p, EnP, arrdim, kappa, eta)

      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the targets
! in the current cluster due to the 
! current source  (determined by the global variable TARPOS). 
!
      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p
      _REAL_,DIMENSION(arrdim),INTENT(INOUT) :: EnP
      _REAL_,INTENT(IN) :: kappa, eta

! local variables

      INTEGER :: i,j,k,nn,nk,nj
      _REAL_ :: tx,ty,tz,xl,yl,zl,rad,radk,radj
      INTEGER :: xlind,ylind,zlind,xhind,yhind,zhind,xyind
      _REAL_ :: kap_eta_2, kap_rad, rad_eta

      kap_eta_2 = kappa * eta / 2_r8

      xl=p%xyz_min(1)
      yl=p%xyz_min(2)
      zl=p%xyz_min(3)

      xlind=p%xyz_lowindex(1)
      ylind=p%xyz_lowindex(2)
      zlind=p%xyz_lowindex(3)

      xhind=p%xyz_highindex(1)
      yhind=p%xyz_highindex(2)
      zhind=p%xyz_highindex(3)

      xyind=xyz_dimglob(1)*xyz_dimglob(2)
      DO k=zlind,zhind
         tz=zl+(k-zlind)*xyz_ddglob(3)-tarpos(3)
         nk=k*xyind
         radk=tz*tz

         DO j=ylind,yhind
            ty=yl+(j-ylind)*xyz_ddglob(2)-tarpos(2)
            nj=nk+j*xyz_dimglob(1)
            radj=radk+ty*ty

            DO i=xlind,xhind
               tx=xl+(i-xlind)*xyz_ddglob(1)-tarpos(1)
               nn=nj+i+1

               rad = SQRT(radj + tx*tx)
               kap_rad = kappa * rad
               rad_eta = rad / eta

               EnP(nn) = EnP(nn) - tarposq / rad &
                                 * (exp(-kap_rad) * erfc(kap_eta_2 - rad_eta) &
                                 -  exp( kap_rad) * erfc(kap_eta_2 + rad_eta))

            END DO
         END DO
      END DO   

      RETURN
      END SUBROUTINE COMP_DIRECT_TCF


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_DIRECT_DCF(p, EnP, arrdim, eta)

      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the targets
! in the current cluster due to the 
! current source  (determined by the global variable TARPOS). 
!
      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p
      _REAL_,DIMENSION(arrdim),INTENT(INOUT) :: EnP
      _REAL_,INTENT(IN) :: eta

! local variables

      INTEGER :: i,j,k,nn,nk,nj
      _REAL_ :: tx,ty,tz,xl,yl,zl,rad,radk,radj
      INTEGER :: xlind,ylind,zlind,xhind,yhind,zhind,xyind
      _REAL_ :: rad_eta,yzhind

      xl=p%xyz_min(1)
      yl=p%xyz_min(2)
      zl=p%xyz_min(3)

      xlind=p%xyz_lowindex(1)
      ylind=p%xyz_lowindex(2)
      zlind=p%xyz_lowindex(3)

      xhind=p%xyz_highindex(1)
      yhind=p%xyz_highindex(2)
      zhind=p%xyz_highindex(3)

      xyind=xyz_dimglob(1)*xyz_dimglob(2)
      DO k=zlind,zhind
         tz=zl+(k-zlind)*xyz_ddglob(3)-tarpos(3)
         nk=k*xyind
         radk=tz*tz

         DO j=ylind,yhind
            ty=yl+(j-ylind)*xyz_ddglob(2)-tarpos(2)
            nj=nk+j*xyz_dimglob(1)
            radj=radk+ty*ty

            DO i=xlind,xhind
               tx=xl+(i-xlind)*xyz_ddglob(1)-tarpos(1)
               nn=nj+i+1

               rad = SQRT(radj+tx*tx)
               rad_eta = rad / eta

               EnP(nn) = EnP(nn) - tarposq / rad * erf(rad_eta)

            END DO
         END DO
      END DO   

      RETURN
      END SUBROUTINE COMP_DIRECT_DCF


!!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_DIRECT_COULOMB(p, EnP, arrdim)

      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the targets
! in the current cluster due to the 
! current source  (determined by the global variable TARPOS). 
!
      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p
      _REAL_,DIMENSION(arrdim),INTENT(INOUT) :: EnP

! local variables

      INTEGER :: i,j,k,nn,nk,nj
      _REAL_ :: tx,ty,tz,xl,yl,zl,radk,radj
      INTEGER :: xlind,ylind,zlind,xhind,yhind,zhind,xyind

      xl=p%xyz_min(1)
      yl=p%xyz_min(2)
      zl=p%xyz_min(3)

      xlind=p%xyz_lowindex(1)
      ylind=p%xyz_lowindex(2)
      zlind=p%xyz_lowindex(3)

      xhind=p%xyz_highindex(1)
      yhind=p%xyz_highindex(2)
      zhind=p%xyz_highindex(3)

      xyind=xyz_dimglob(1)*xyz_dimglob(2)
      DO k=zlind,zhind
         tz=zl+(k-zlind)*xyz_ddglob(3)-tarpos(3)
         nk=k*xyind
         radk=tz*tz

         DO j=ylind,yhind
            ty=yl+(j-ylind)*xyz_ddglob(2)-tarpos(2)
            nj=nk+j*xyz_dimglob(1)
            radj=radk+ty*ty

            DO i=xlind,xhind
               tx=xl+(i-xlind)*xyz_ddglob(1)-tarpos(1)
               nn=nj+i+1

               EnP(nn) = EnP(nn) - tarposq / SQRT(radj+tx*tx)

            END DO
         END DO
      END DO   

      RETURN
      END SUBROUTINE COMP_DIRECT_COULOMB


!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE PC_TREECODE_FORCES(p, guv, &                       !source info
                                    solutePosition, soluteCharge, & !target info
                                    numparsS, numparsT, &  !number sources, targets
                                    forcesT)               !output
      IMPLICIT NONE
!
! PC_TREECODE_FORCES is the driver routine which calls COMPUTE_PC_FORCES for each
! source particle, setting the global variable TARPOS before the call.
! P is the root node of the source tree, solutePosition are the coordinates of
! the target particles, and guv contains source grid information.
!

      INTEGER,INTENT(IN) :: numparsS,numparsT
      TYPE(tnode),POINTER :: p
      _REAL_,DIMENSION(3,numparsS),INTENT(IN) :: solutePosition
      _REAL_,DIMENSION(numparsS),INTENT(IN) :: soluteCharge
      _REAL_,DIMENSION(numparsS),INTENT(IN) :: guv
      _REAL_,DIMENSION(3,numparsT),INTENT(INOUT) :: forcesT

! local variables

      INTEGER :: i
      _REAL_,DIMENSION(3) :: force(3)


      forcesT = 0.0_r8

      DO i = 1, numparsT
          force = 0.0_r8
          tarpos(1) = solutePosition(1, i)
          tarpos(2) = solutePosition(2, i)
          tarpos(3) = solutePosition(3, i)
          tarposq = soluteCharge(i)

          CALL COMPUTE_PC_FORCES(p, guv, numparsS, force)

! We have to do this backwards because we invert our grid directions to 
! agree with RISM's directions
          forcesT(1,i) = -tarposq * force(1)
          forcesT(2,i) = -tarposq * force(2)
          forcesT(3,i) = -tarposq * force(3)

      END DO

      RETURN
      END SUBROUTINE PC_TREECODE_FORCES


!!!!!!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE COMPUTE_PC_FORCES(p,guv,numpars,force)
      IMPLICIT NONE

! COMPUTE_PC_FORCES is the recursive routine for computing the interaction
! between a target particle and a source cluster. If the MAC is
! satisfied the power series coefficients for the current target
! cluster are updated. If the MAC is not satisfied then the algorithm
! descends to the children of the current cluster, unless the
! current cluster is a leaf then the interaction is done exactly
! via a call to the routine COMP_DIRECT
!

      INTEGER,INTENT(IN) :: numpars
      TYPE(tnode),POINTER :: p
      _REAL_,DIMENSION(numpars),INTENT(IN) :: guv
      _REAL_,DIMENSION(3),INTENT(OUT) :: force(3)

! local variables

      _REAL_,DIMENSION(3) :: xyz_t, forcetemp
      _REAL_ :: distsq
      INTEGER :: i, err, kk, k1, k2, k3

! determine DISTSQ for MAC test

      xyz_t=tarpos-p%xyz_mid
      distsq=SUM(xyz_t**2)

! intialize potential energy and force

      force = 0.0_r8

! If MAC is accepted and there is more than 1 particle in the
! box use the expansion for the approximation.
      IF ((p%sqradius .LT. distsq*thetasq) .AND. &
         (p%sqradius .NE. 0.0_r8)) THEN

         IF (p%exist_ms .EQ. 0) THEN
             p%ms => safemem_realloc(p%ms, torderflat, .false.)
             p%ms=0.0_r8
             p%exist_ms=1
             CALL COMP_MS_FORCES(p, guv, numpars)
         END IF

         CALL COMP_TCOEFF_COULOMB(xyz_t(1), xyz_t(2), xyz_t(3))

         kk = 0
         DO k3 = 0, torder
            DO k2 = 0, torder-k3
               DO k1 = 0, torder-k3-k2
                  kk = kk + 1
                  force(1) = force(1) + cf(k1) * a(k1+1,k2,k3) * p%ms(kk)
                  force(2) = force(2) + cf(k2) * a(k1,k2+1,k3) * p%ms(kk)
                  force(3) = force(3) + cf(k3) * a(k1,k2,k3+1) * p%ms(kk)
               END DO
            END DO
         END DO

      ELSE

! If MAC fails check to see if the are children. If not, perform direct
! calculation.  If there are children, call routine recursively for each.
!
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT_FORCES(p, guv, numpars, force)
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_PC_FORCES(p%child(i)%p_to_tnode, guv, &
                               numpars, forcetemp)
               force = force + forcetemp
            END DO
         END IF
      END IF

      RETURN
      END SUBROUTINE COMPUTE_PC_FORCES


!!!!!!!!!!!!!!!


      SUBROUTINE COMP_MS_FORCES(ap,guv,numpars)
      IMPLICIT NONE
!
! COMP_CMS computes the moments for node P needed in the Taylor approximation
!
      INTEGER,INTENT(IN) :: numpars
      TYPE(tnode),POINTER :: ap
      _REAL_,DIMENSION(numpars),INTENT(IN) :: guv

      _REAL_ :: tx,ty,tz,dx,dy,dz,xl,yl,zl,xm,ym,zm
      INTEGER :: xlind, ylind, zlind, xhind, yhind, zhind, xyind
      INTEGER :: i,nn,j,k,k1,k2,k3,kk


      xl=ap%xyz_min(1)
      yl=ap%xyz_min(2)
      zl=ap%xyz_min(3)

      xm=ap%xyz_mid(1)
      ym=ap%xyz_mid(2)
      zm=ap%xyz_mid(3)

      xlind=ap%xyz_lowindex(1)
      ylind=ap%xyz_lowindex(2)
      zlind=ap%xyz_lowindex(3)

      xhind=ap%xyz_highindex(1)
      yhind=ap%xyz_highindex(2)
      zhind=ap%xyz_highindex(3)

      xyind=xyz_dimglob(1)*xyz_dimglob(2)
      ap%ms=0.0_r8

      DO k=zlind,zhind
         dz=zl+(k-zlind)*xyz_ddglob(3)-zm
         DO j=ylind,yhind
            dy=yl+(j-ylind)*xyz_ddglob(2)-ym
            DO i=xlind,xhind
               dx=xl+(i-xlind)*xyz_ddglob(1)-xm

               nn=(k*xyind)+(j*xyz_dimglob(1))+i+1

               kk = 0
               tz = 1.0_r8

               DO k3 = 0, torder
                  ty = 1.0_r8
                  DO k2 = 0, torder-k3
                     tx = 1.0_r8
                     DO k1 = 0, torder-k3-k2
                        kk = kk + 1
                        ap%ms(kk) = ap%ms(kk) + guv(nn)*tx*ty*tz
                        tx = dx * tx
                     END DO
                     ty = dy * ty
                  END DO
                  tz = tz * dz
               END DO

            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE COMP_MS_FORCES


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_DIRECT_FORCES(p, guv, numpars, force)

      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the targets
! in the current cluster due to the current source 
! (determined by the global variable TARPOS).
!
      INTEGER,INTENT(IN) :: numpars
      TYPE(tnode),POINTER :: p
      _REAL_,DIMENSION(numpars),INTENT(IN) :: guv
      _REAL_,DIMENSION(3),INTENT(OUT) :: force(3)

! local variables

      INTEGER :: i,j,k,nn
      _REAL_ :: tx,ty,tz,xl,yl,zl,dist
      INTEGER :: xlind,ylind,zlind,xhind,yhind,zhind,xyind

      xl=p%xyz_min(1)
      yl=p%xyz_min(2)
      zl=p%xyz_min(3)

      xlind=p%xyz_lowindex(1)
      ylind=p%xyz_lowindex(2)
      zlind=p%xyz_lowindex(3)

      xhind=p%xyz_highindex(1)
      yhind=p%xyz_highindex(2)
      zhind=p%xyz_highindex(3)

      force = 0.0_r8
      xyind=xyz_dimglob(1)*xyz_dimglob(2)

      DO k=zlind,zhind
         tz=zl+(k-zlind)*xyz_ddglob(3)-tarpos(3)
         DO j=ylind,yhind
            ty=yl+(j-ylind)*xyz_ddglob(2)-tarpos(2)
            DO i=xlind,xhind
               tx=xl+(i-xlind)*xyz_ddglob(1)-tarpos(1)

               nn=(k*xyind)+(j*xyz_dimglob(1))+i+1

               dist = 1.0_r8 / (tx*tx + ty*ty + tz*tz)**(1.5_r8)

               force(1) = force(1) + guv(nn) * tx * dist
               force(2) = force(2) + guv(nn) * ty * dist
               force(3) = force(3) + guv(nn) * tz * dist

            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE COMP_DIRECT_FORCES


!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE CLEANUP(p)
      IMPLICIT NONE
!
! CLEANUP deallocates allocated global variables and then
! calls recursive routine REMOVE_NODE to delete the tree.
!
      TYPE(tnode),POINTER :: p      

! local variables
  
      INTEGER :: err

      DEALLOCATE(cf,cf1,cf2,cf3,a,b, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating Taylor variables! '
         STOP
      END IF

      CALL REMOVE_NODE(p)
      
      DEALLOCATE(p, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating root node! '
         STOP
      END IF 
      NULLIFY(p)         

      RETURN
      END SUBROUTINE CLEANUP


!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE REMOVE_NODE(p)
      IMPLICIT NONE
!
! REMOVE_NODE recursively removes each node from the
! tree and deallocates its memory for MS array if it
! exits.
!
      TYPE(tnode),POINTER :: p 

! local variables

      INTEGER :: i,err

      IF (p%exist_ms .EQ. 1) THEN
         IF (safemem_dealloc(p%ms) /= 0) &
            CALL rism_report_error("p%ms deallocation failed")
      END IF

      IF (p%num_children .GT. 0) THEN
          DO i=1,p%num_children
             CALL REMOVE_NODE(p%child(i)%p_to_tnode)

             DEALLOCATE(p%child(i)%p_to_tnode,STAT=err)
             IF (err .NE. 0) THEN
                WRITE(6,*) 'Error deallocating node child!'
                STOP
             END IF                           
          END DO
      END IF 

      RETURN                
      END SUBROUTINE REMOVE_NODE      

      END MODULE treecode_procedures


!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE TREECODE(solutePosition, soluteCharge, &
                          targetGridLimits, targetGridDims, &
                          numSources, numTargets, potentialGrid, &
                          treeOrder, treeTheta, &
                          maxLeafTargets, potentialType, &
                          paramKappa, paramEta, paramEpsilon)

      USE treecode_procedures
      IMPLICIT NONE

!====================================================================
!                                                                   
! solutePosition   :: x,y,z coordinates of solute atoms (sources)
! soluteCharge     :: partial charge of solute atoms (sources)
! targetGridDims   :: target grid dimensions
! targetGridLimits :: min, max target grid limits
! numSources       :: number of sources
! numTargets       :: number of targets
! potentialGrid    :: array dimension numparsT for storing potential
!                     at each target
! maxLeafTargets   :: maximum number of particles in a leaf
!                     (only employed in the cluster-particle version)
!=====================================================================

      INTEGER,INTENT(IN) :: numSources, numTargets
      INTEGER,INTENT(IN) :: treeOrder, maxLeafTargets
      _REAL_,DIMENSION(3,numSources),INTENT(IN) :: solutePosition
      _REAL_,DIMENSION(numSources),INTENT(IN) :: soluteCharge

      _REAL_,DIMENSION(6),INTENT(IN) :: targetGridLimits
      INTEGER,DIMENSION(3),INTENT(IN) :: targetGridDims

      _REAL_,DIMENSION(numTargets),INTENT(OUT) :: potentialGrid
      _REAL_,INTENT(IN) :: treeTheta

      INTEGER,INTENT(IN) :: potentialType
      _REAL_,INTENT(IN) :: paramKappa, paramEta, paramEpsilon

! local variables

      TYPE(tnode),POINTER :: treeRootNode
      INTEGER :: treeLevel
      INTEGER,DIMENSION(6) :: targetGridIndices

! Call SETUP to allocate arrays for Taylor expansions
! and setup global variables. 

      CALL SETUP(targetGridLimits, targetGridDims, targetGridIndices, &
                 treeOrder, treeTheta)

! nullify pointer to root of tree (TROOT) and create tree

      NULLIFY(treeRootNode)

! construct tree

      treeLevel = 0
      CALL CREATE_TREE_N0(treeRootNode, maxLeafTargets, &
                          targetGridLimits, targetGridDims, targetGridIndices, &
                          treeLevel)

!Call driver routine for cluster-particle
      IF (potentialType == 0) THEN
          CALL CP_TREECODE_TCF(treeRootNode, solutePosition, soluteCharge, potentialGrid, &
                               numSources, numTargets, paramKappa, paramEta, paramEpsilon)
      ELSE IF (potentialType == 1) THEN
          CALL CP_TREECODE_DCF(treeRootNode, solutePosition, soluteCharge, potentialGrid, &
                               numSources, numTargets, paramEta)
      ELSE IF (potentialType == 2) THEN
          CALL CP_TREECODE_COULOMB(treeRootNode, solutePosition, soluteCharge, potentialGrid, &
                                   numSources, numTargets)
      END IF

      CALL CLEANUP(treeRootNode)

      END SUBROUTINE TREECODE


!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE TREECODE_FORCES(solutePosition, soluteCharge, &
                                 guv, sourceGridLimits, sourceGridDims, &
                                 numSources, numTargets, &
                                 treeOrder, treeTheta, maxLeafSources, &
                                 voxelVolume, soluteForces)

      USE treecode_procedures
      IMPLICIT NONE

!====================================================================
!
! solutePosition   :: x,y,z coordinates of solute atoms (sources)
! soluteCharge     :: partial charge of solute atoms (sources)
! sourceGridDims   :: source grid dimensions
! sourceGridLimits :: min, max source grid limits
! numSources       :: number of sources
! numTargets       :: number of targets
! soluteForces     :: array of dimension (3,numTargets) for storing
!                     forces at each target
! maxLeafSources   :: maximum number of particles in a leaf
!=====================================================================

      INTEGER,INTENT(IN) :: numSources, numTargets
      INTEGER,INTENT(IN) :: treeOrder, maxLeafSources
      _REAL_,DIMENSION(3,numSources),INTENT(IN) :: solutePosition
      _REAL_,DIMENSION(numSources),INTENT(IN) :: soluteCharge

      _REAL_,DIMENSION(numSources),INTENT(IN) :: guv
      _REAL_,DIMENSION(6),INTENT(IN) :: sourceGridLimits
      INTEGER,DIMENSION(3),INTENT(IN) :: sourceGridDims

      _REAL_,INTENT(IN) :: treeTheta
      _REAL_,INTENT(IN) :: voxelVolume

      _REAL_,DIMENSION(3,numTargets),INTENT(OUT) :: soluteForces

! local variables

      TYPE(tnode),POINTER :: treeRootNode
      INTEGER :: treeLevel
      INTEGER,DIMENSION(6) :: sourceGridIndices


! Call SETUP to allocate arrays for Taylor expansions
! and setup global variables.

      CALL SETUP(sourceGridLimits, sourceGridDims, sourceGridIndices, &
                 treeOrder, treeTheta)

! nullify pointer to root of tree (TROOT) and create tree

      NULLIFY(treeRootNode)

! set global variables to track tree levels during construction

      treeLevel = 0
      CALL CREATE_TREE_N0(treeRootNode, maxLeafSources, &
                          sourceGridLimits, sourceGridDims, sourceGridIndices, &
                          treeLevel)

!Call driver routine for particle-cluster forces
      CALL PC_TREECODE_FORCES(treeRootNode, guv, &              !source info
                              solutePosition, soluteCharge,  &  !target info
                              numSources, numTargets, &         !number sources, targets
                              soluteForces)                     !output
      soluteForces = soluteForces * voxelVolume

! Call CLEANUP to deallocate global variables and tree structure.
      CALL CLEANUP(treeRootNode)

      END SUBROUTINE TREECODE_FORCES

