#include "../include/dprec.fh"
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)

module density_surface

#  include "pb_constants.h"

   integer, allocatable :: index(:,:,:)
   integer, allocatable :: index2(:,:,:)
   integer, allocatable :: tempcnt(:)
   integer, allocatable :: ndenatm(:,:,:)
   integer, allocatable :: bndatmptr(:)
   integer, allocatable :: bndatmlst(:)
   _REAL_, allocatable :: cirreg(:,:)
   _REAL_, allocatable :: x(:), y(:), z(:)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ set level set function using the revised density approach
subroutine density_lvlset(eneopt,memopt,outlvlset,outmlvlset,natom,xm,ym,zm,gox,goy,goz,h,&
              dprobh,cutfd,mzmin,mzmax,gcrd,radi,poretype,poreradi,&
              insas,inmem,lvlset,mlvlset)

   implicit none

   integer eneopt
   integer memopt
   logical outlvlset
   logical outmlvlset
   integer natom
   integer xm,ym,zm
   _REAL_ h,gox,goy,goz
   _REAL_ dprobh
   _REAL_ cutfd
   _REAL_ mzmin,mzmax
   _REAL_ gcrd(3,natom)
   _REAL_ radi(natom)
   integer poretype
   _REAL_ poreradi
   integer insas(0:xm+1,0:ym+1,0:zm+1)
   integer inmem(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ lvlset(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ mlvlset(0:xm+1,0:ym+1,0:zm+1)

   integer lvlsetfn, mlvlfn
   character*16 lvlsetdataname, mlvldataname
   character*16 lvlsetfilename
   character*16 mlvlfilename

   integer i, j, k, l, ier
   integer lowi, lowj, lowk
   integer highi, highj, highk
   _REAL_ xi, yi, zi, rh
   _REAL_ range0, range1, range2, range3
   _REAL_ dist, cutoffh
   _REAL_ poredata(3,0:zm+1)
   _REAL_ factor

   ! WMBS: output file logistics to gen_dx

   lvlsetfilename="pbsa_lvlset.dx"; mlvlfilename="pbsa_mlvlset.dx"
   lvlsetdataname="Level Set"; mlvldataname="Level Set"
   lvlsetfn=3334; mlvlfn=3335

   ! save some info for atom list of irregular grid points

   if ( allocated(ndenatm) ) then
      deallocate(ndenatm, stat = ier); REQUIRE(ier==0)
   end if
   allocate(ndenatm(1:xm,1:ym,1:zm), stat = ier ); REQUIRE(ier==0)
   ndenatm = 0

   ! WMSB: Add atomic density contribution using the inkblot method as in SES/SAS
   ! LX: the sphere is now set to be max of cutfd and radi+2dprob

   rh = ONE/h
   factor = HALF/dprobh
   cutoffh = sqrt(cutfd)*rh

   do l = 1, natom
      range0 = radi(l)*rh ! atom radius in gridpoint units
      if ( range0 == ZERO ) cycle ! don't waste time on zero radius atoms
      xi = gcrd(1,l); yi = gcrd(2,l); zi = gcrd(3,l)
      if (eneopt == 4) then
         range1 = max(cutoffh, range0 + dprobh * 2.d0) ! distance till atomic lvlset contrib -> 0
      else
         range1 = max(cutoffh, range0 + dprobh * 2.d0 + ONE) ! H + distance till atomic lvlset contrib -> 0
      end if
      if ( zi+range1<0 .or. zi-range1>zm+1 ) cycle ! this shouldn't happen...
      lowk = max(1,ceiling(zi - range1)); highk = min(zm,floor(zi + range1))
      do k = lowk, highk ! z indices (grid line)
         range2 = sqrt(range1**2-(zi-dble(k))**2)
         if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle ! this shouldn't happen...
         lowj = max(1,ceiling(yi - range2)); highj = min(ym,floor(yi + range2))
         do j = lowj, highj ! y indices (grid disc)
            range3 = sqrt(range2**2-(yi-dble(j))**2)
            if ( range3 > ZERO ) then ! sanity check on range3
               lowi = max(1,ceiling(xi-range3)); highi = min(xm,floor(xi+range3))
               do i = lowi, highi ! x indices (grid sphere)
                  dist = sqrt((xi-i)**2+(yi-j)**2+(zi-k)**2) - range0
                  dist = dist *factor
                  lvlset(i,j,k) = lvlset(i,j,k) + density_atom(dist)
                  ndenatm(i,j,k) = ndenatm(i,j,k) + 1
               end do ! loop over x indices (grid sphere)
            end if ! sanity check on range3
         end do ! loop over y indicies (grid disc)
      end do ! loop over z indicies (grid line)
   end do

   ! WMBS: set insas using the solute level set function only due to the fractional
   ! edge algorithm in epsmap(), which is a central place to map all surface
   ! definitions into esp edges.
   ! CHW: the insas assignment is revised so boundary grid is inside as well.

   insas(0:xm+1,0:ym+1,0:zm+1) = int(sign(1.001d0, lvlset(0:xm+1,0:ym+1,0:zm+1) - ONE))

   ! WMSB: compute membrane density from mzmin-2*dprob to mzmax+2*dprob only
   ! set up the membrane pore if needed

   ! RL: the inmem label will be set to 1 for nodes outside the solute but
   ! inside the membrane, otherwise it is 0. This is consistent with the ses
   ! method, see exvwslab()

   ! convert density to level set function

   if ( memopt > 0 ) then
      if ( poretype == 1 ) call gen_cylinder_data(xm,ym,zm,mzmin,mzmax,poreradi,natom,gcrd,poredata)
      call membrane_density(xm,ym,zm,dprobh,mzmin,mzmax,poretype,poredata,mlvlset)
      lowk = max(0,ceiling(mzmin)); highk = min(zm+1,floor(mzmax))
      do k = lowk, highk
      do j = 0, ym+1; do i = 0, xm+1
         if ( insas(i,j,k) < 0 ) inmem(i,j,k) = 1
      end do; end do
      end do
      mlvlset = ONE - lvlset - mlvlset
   end if
   lvlset = ONE - lvlset

   if ( outmlvlset ) then
      !call gen_dx_file(xm,ym,zm,h,gox,goy,goz,mlvlset(1:xm,1:ym,1:zm),mlvlfilename,mlvlfn,mlvldataname)
      call gen_integer_dx_file(xm,ym,zm,h,gox,goy,goz,inmem(1:xm,1:ym,1:zm),mlvlfilename,mlvlfn,mlvldataname)
   end if
   if ( outlvlset ) then
      !call gen_dx_file(xm,ym,zm,h,gox,goy,goz,lvlset(1:xm,1:ym,1:zm),lvlsetfilename,lvlsetfn,lvlsetdataname)
      call gen_integer_dx_file(xm,ym,zm,h,gox,goy,goz,insas(1:xm,1:ym,1:zm),lvlsetfilename,lvlsetfn,lvlsetdataname)
   end if

end subroutine density_lvlset
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ set up the atom neighbor list for boundary grid points
subroutine bndatm(eneopt,nbnd,natom,xm,ym,zm,h,dprobh,cutfd,gcrd,radi,iepsav,lvlset)

   implicit none

   ! passed variables

   integer eneopt
   integer nbnd
   integer natom
   integer xm,ym,zm
   _REAL_ h
   _REAL_ dprobh
   _REAL_ cutfd
   _REAL_ gcrd(3,natom)
   _REAL_ radi(natom)
   integer iepsav(4,xm*ym*zm)
   _REAL_ lvlset(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer maxlst
   integer ier
   integer i, j, k, l
   integer lowi, lowj, lowk
   integer highi, highj, highk
   integer ibnd
   _REAL_ xi, yi, zi, rh, cutoffh
   _REAL_ range0, range1, range2, range3

   ! allocation and preparation
   ! maxlst is the maximum neighbors of an atom within cutfd
   ! with double counting. Each atom is repeated 8 times due to FD

   if ( allocated(bndatmptr) ) then
      deallocate(bndatmptr, stat = ier); REQUIRE(ier==0)
   end if
   allocate(bndatmptr(0:nbnd+8*natom), stat = ier); REQUIRE(ier==0)
   maxlst = 16*nint( dble(natom)*(THIRD*sqrt(cutfd)**3+TEN) )

   bndatmptr(0) = 0
   do l = 1, nbnd
      i = iepsav(1,l); j = iepsav(2,l); k = iepsav(3,l)
      maxlst = maxlst + ndenatm(i,j,k)
      bndatmptr(l) = bndatmptr(l-1) + ndenatm(i,j,k)
   end do

   if ( allocated(bndatmlst) ) then
      deallocate(bndatmlst, stat = ier); REQUIRE(ier==0)
   end if
   allocate(bndatmlst(maxlst), stat = ier); REQUIRE(ier==0)
   if ( allocated(tempcnt)   ) then
      deallocate(tempcnt  , stat = ier); REQUIRE(ier==0)
   end if
   allocate(tempcnt(1:nbnd+8*natom)  , stat = ier); REQUIRE(ier==0)
   bndatmlst = 0
   tempcnt = 0

   ! find all atom neighbors of a given irregular point

   rh = ONE/h
   cutoffh = sqrt(cutfd)*rh
   do l = 1, natom
      range0 = radi(l)*rh ! atom radius in gridpoint units
      if ( range0 == ZERO ) cycle ! don't waste time on zero radius atoms
      xi = gcrd(1,l); yi = gcrd(2,l); zi = gcrd(3,l)
      if (eneopt == 4) then
         range1 = max(cutoffh, range0 + dprobh * 2.d0) ! distance till atomic lvlset contrib -> 0
      else
         range1 = max(cutoffh, range0 + dprobh * 2.d0 + ONE) ! H + distance till atomic lvlset contrib -> 0
      end if
      if ( zi+range1<0 .or. zi-range1>zm+1 ) cycle ! this shouldn't happen...
      lowk = max(1,ceiling(zi - range1)); highk = min(zm,floor(zi + range1))
      do k = lowk, highk ! z indices (grid line)
         range2 = sqrt(range1**2-(zi-dble(k))**2)
         if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle ! this shouldn't happen...
         lowj = max(1,ceiling(yi - range2)); highj = min(ym,floor(yi + range2))
         do j = lowj, highj ! y indices (grid disc)
            range3 = sqrt(range2**2-(yi-dble(j))**2)
            if ( range3 > ZERO ) then ! sanity check on range3
               lowi = max(1,ceiling(xi-range3)); highi = min(xm,floor(xi+range3))
               do i = lowi, highi ! x indices (grid sphere)
                  ibnd = index2(i,j,k)
                  if ( ibnd > 0 ) then
                     tempcnt(ibnd) = tempcnt(ibnd) + 1
                     bndatmlst(bndatmptr(ibnd-1)+tempcnt(ibnd)) = l ! this is the atom no. to be saved
                  end if
               end do ! loop over x indices (grid sphere)
            end if ! sanity check on range3
         end do ! loop over y indicies (grid disc)
      end do ! loop over z indices (grid line)
   end do


end subroutine bndatm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!+ membrane density function
subroutine membrane_density(xm,ym,zm,dprobh,mzmin,mzmax,poretype,poredata,u)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! calculates the levelset value contribution due to the membrane
!   mzmin    specifies the range of grid nodes whose level set function will
!   mzmax    be influenced by the membrane,
!   poretype specifies an exclusion region
!   poretype meanings:
!       value       meaning
!       -----       ------
!       0           no pore
!       1           cylindrical pore
!       2-4         reserved for future implementations
!   poredata is set up before hand and passed in.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! passed variables

   integer xm,ym,zm    ! grid array index dimensions
   _REAL_ dprobh       ! water probe radius in grid unit
   _REAL_ mzmin, mzmax ! membrane extrema in grid unit
   integer poretype    ! pore option
   _REAL_ poredata(3,0:zm+1) ! pore parameters
   _REAL_ mdist(0:xm+1,0:ym+1,0:zm+1) ! input working distance array
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1) ! output density, should be initialized to zero

   ! local variables

   integer lowk, highk, k
   _REAL_  mzctr, mthick2 ! center and half thickness in the grid unit

   ! part a: compute mdist
   ! mdist is signed distance to slab region boundaries in the grid unit
   ! it is positive if the node is outside the membrane bounds and
   ! it is negative if it is inside

   mzctr = HALF*(mzmin + mzmax); mthick2 = HALF*(mzmax - mzmin)
   lowk = max(0,ceiling(mzmin-TWO*dprobh)); highk = min(zm+1,floor(mzmax+TWO*dprobh))

   mdist(0:xm+1,0:ym+1,0:lowk-1) = 999.9d0
   do k = lowk, highk
      mdist(0:xm+1,0:ym+1,k) = abs(k-mzctr) - mthick2
   end do
   mdist(0:xm+1,0:ym+1,highk+1:zm+1) = 999.9d0

   ! part b: correct mdist due to the existence of pore
   ! the signed membrane distance is corrected using the distance
   ! to the cylinder

   if ( poretype > 0 ) then
      call distance_correct_cylinder(xm,ym,lowk,highk,mdist)
   end if

   ! part c: compute density function with mdist

   ! Note that the probe radius is always less than the membrane
   ! thicknes. Otherwise, program would stop. This ensures the membrane
   ! levelset will
   ! 1) Not have a cusp at the center (ensured by smoothing function)
   ! 2) Not to be larger than mprobe

   mdist = mdist/(TWO*dprobh)
   call distance_to_density(xm,ym,lowk,highk,mdist,u)

   contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!+ correct signed distance value when a cylinder core region presents
subroutine distance_correct_cylinder(xm,ym,lowk,highk,mdist)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine modifies the slab-like membrane signed distance when a cylinder
! exclusion presents. The poredata array contains the radius of this region at
! z plane of the grid. This is to potentially allow more complicated surfaces,
! such as slanted solids of revolution in later versions but modification of
! the cylinder surface distance calculation will be needed if they are used.
! All distances and coordinates are supposed in the grid unit.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! passed variables

   integer xm, ym, lowk, highk
   _REAL_ mdist(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer i,j,k                  ! grid location indices for calculation
   _REAL_ cdist                   ! distance to cylinder surface
   _REAL_ ddist                   ! scratch variable
   _REAL_ xi,yi,zi                ! grid cooridinates from indicies
   _REAL_ dx, dy

   do k = lowk, highk
   do j = 0, ym+1
   do i = 0, xm+1

      ! part a. find the signed distance to the cylinder

      ! calculates the distance between node i,j,k and the surface of the cylinderical
      ! exclusion region. poredata(k,1) and poredata(k,2) are x/y location of center
      ! of the cylinderical region at z=k. poredata(k,3) is the radius of the region
      ! at k. The node is within the cylincer if the distance is negative and the node
      ! is outside the cylinder if the distnace is positive.

      dx = i - poredata(1,k)
      dy = j - poredata(2,k)
      cdist = sqrt(dx**2 + dy**2) - poredata(3,k)

      ! part b. correct mdist based on the distance to the cylinder

      ddist = mdist(i,j,k)

      if ( ddist < 0 ) then

         ! membrane interior

         if ( cdist <= 0 ) then
            ! we are inside the cylinder. Use the distance to the cylinder surface
            ! Note it is positive since it is outside membrane.
            ddist = -cdist
         else
            ! we are outside the cylinder. Distance is calculated as the
            ! geometric mean of membrane and cylinder distances when the
            ! cylinder surface is closer. If the membrane bound is closer,
            ! the cylinder distance is ignored.
            ! Note it is negative since it is inside membrane.
            ddist = -sqrt(abs(ddist*min(abs(ddist),abs(cdist))))
         end if

      else

         ! membrane exterior

         if ( cdist <= 0 ) then
            ! we are within the cylinder so we will need to compute the
            ! distance to the cirucular edge of the pore, which ends at the
            ! membrane bounds. The distance is the root square sum of the distance
            ! to the cylinder surface and the distance to the membrane bound.
            ddist = sqrt(ddist**2 + cdist**2)
         end if
            ! We are outside the cylinder, always use the distance to membrane
            ! whether it is closer to the membrane bound or not.
      end if

      ! part c. returning

      mdist(i,j,k) = ddist

   end do
   end do
   end do

end subroutine distance_correct_cylinder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!+ Compute density value for the membrane interior
_REAL_ function interior_density(dist)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute density value for the membrane interior given distance. Units of
! distance here should have been converted into the unit of solvent probe diameter as
! in other density functions before passing in. The output is density.
!
! The distance will be transformed into a density with a cubic polynomial from 1
! (lvlset of 0) to a constant of 2 (lvlset of -1) over a range of 1/2 (i.e. radius).
! The slope @ d = 0 matchs that of the spline function for the external density
! @ d = 0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! passed variables

   _REAL_ dist

   ! local variables

   _REAL_, parameter:: sdata(6) = (/1.0,-4.527143,-7.909921,-6.566905,-2.512619,-0.328492/)

   if ( dist < -ONE ) then
      interior_density = TWO
   else
      interior_density =      &
         sdata(1)           + &
         sdata(2)*(dist   ) + &
         sdata(3)*(dist**2) + &
         sdata(4)*(dist**3) + &
         sdata(5)*(dist**4) + &
         sdata(6)*(dist**5)
   end if

end function interior_density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!+ Compute density value for the membrane exterior
subroutine distance_to_density(xm,ym,lowk,highk,dist,u)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Similar to the density function in pb_exmol. The input (dist) should be
! converted into the unit of solvent probe diameter before passing in. The
! output (u) density is computed with the same spline based fitting
! as in the atomic density function. Since the membrane density is done after the
! atom density is done, it is assumed that the density function array has
! been initialized, i.e. zeroed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! passed variables

   integer xm, ym, lowk, highk
   _REAL_ dist(0:xm+1,0:ym+1,0:zm+1) ! scaled signed distance
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)    ! output density

   ! local variables

   integer i,j,k,m
   _REAL_, parameter :: dash(6) = (/0.00d0,0.20d0,0.40d0,0.60d0,0.80d0,1.00d0 /)
   _REAL_, parameter :: spcoef(5,4) = &
       reshape((/1.000000  ,0.2100000  ,0.1500000  ,0.0500000  ,0.010000   ,  &
                -4.527143  ,-2.067608  ,0.0475730  ,-0.522686  ,-0.056828  ,  &
                -3.640532  ,15.938209  ,-5.362303  ,2.5110050  ,-0.181716  ,  &
                32.631235  ,-35.500854 ,13.122180  ,-4.487867  ,1.079289/), (/5,4/))

   do k = lowk, highk
   do j = 0, ym+1
   do i = 0, xm+1
      if ( dist(i,j,k) <= ONE ) then
         if ( dist(i,j,k) <= ZERO ) then
            u(i,j,k) = u(i,j,k) + interior_density(dist(i,j,k))
         else
            do m = 1, 5
               if ( dist(i,j,k) > dash(m) .and. dist(i,j,k) <= dash(m+1) ) then
                  u(i,j,k) = u(i,j,k)                     +&
                     spcoef(m,1)                          +&
                     spcoef(m,2)*(dist(i,j,k)-dash(m))    +&
                     spcoef(m,3)*(dist(i,j,k)-dash(m))**2 +&
                     spcoef(m,4)*(dist(i,j,k)-dash(m))**3
               end if
            end do
         end if
      end if
   end do
   end do
   end do

end subroutine distance_to_density

end subroutine membrane_density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!+ Generates the poredata needed to describe a right cylinder exclusion region
subroutine gen_cylinder_data(xm,ym,zm,mzmin,mzmax,poreradih,natom,gcrd,poredata)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Generates the poredata needed to describe a right cylinder exclusion
! region, centered on centeroid of the solute inside the membrane region.
! be sure that the poredata array has been allocated as a 1D array
! -> use allocate( poredata(3*(zm+2)) ).
! All input ouput geometry/distance are in the grid unit.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! passed variables

   integer xm,ym,zm
   _REAL_ mzmin, mzmax
   integer natom
   _REAL_ gcrd(3,natom)
   _REAL_ poreradih
   _REAL_ poredata(3,0:zm+1)

   ! local variables

   integer k,iatm,acount
   _REAL_ centeroid(2)

   write(6,'(a)') "======== Setting up Membrane Cylindrical Exclusion Region ========"
   write(6,'(a,f12.4)') " Cylinder Radius (in h) = ", poreradih

   ! find centeroid of membrane bound solute

   acount = 0
   centeroid = ZERO
   do iatm = 1, natom
      if ( gcrd(3,iatm) >= mzmin .and. gcrd(3,iatm) <= mzmax ) then
         centeroid(1) = centeroid(1) + gcrd(1,iatm)
         centeroid(2) = centeroid(2) + gcrd(2,iatm)
         acount = acount + 1
      end if
   end do
   centeroid = centeroid / dble(acount)
   write(6,'(a,2f12.4)') " Transmembrane center (in h) = ",centeroid; flush(6)

   ! build cylinder data

   do k = 0,zm+1
      poredata(1,k) = centeroid(1)
      poredata(2,k) = centeroid(2)
      poredata(3,k) = poreradih
   end do

end subroutine gen_cylinder_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ irregular grid point index and data structure initialization
subroutine irreg_init(xm,ym,zm,nbnd,h,gox,goy,goz,insas,iepsav,phi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! initialize irregular grid point index and data structure for later use
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none
   integer xm,ym,zm,nbnd
   _REAL_ h, gox, goy, goz
   integer insas(0:xm+1, 0:ym+1, 0:zm+1)
   integer iepsav(4,xm*ym*zm)
   _REAL_ phi(0:xm+1,0:ym+1,0:zm+1)

   integer i, j, k, ii, ier

   ! allocating working arrays for irregular grid points

   if ( allocated(index )  ) then
      deallocate(index , stat = ier); REQUIRE(ier==0)
   end if
   allocate (index (1:xm,1:ym,1:zm), stat = ier); REQUIRE(ier==0)
   if ( allocated(index2)  ) then
      deallocate(index2, stat = ier); REQUIRE(ier==0)
   end if
   allocate (index2(1:xm,1:ym,1:zm), stat = ier); REQUIRE(ier==0)
   if ( allocated(x     )  ) then
      deallocate(x     , stat = ier); REQUIRE(ier==0)
   end if
   allocate (x     (0:xm+1        ), stat = ier); REQUIRE(ier==0)
   if ( allocated(y     )  ) then
      deallocate(y     , stat = ier); REQUIRE(ier==0)
   end if
   allocate (y     (0:ym+1        ), stat = ier); REQUIRE(ier==0)
   if ( allocated(z     )  ) then
      deallocate(z     , stat = ier); REQUIRE(ier==0)
   end if
   allocate (z     (0:zm+1        ), stat = ier); REQUIRE(ier==0)
   if ( allocated(cirreg)  ) then
      deallocate(cirreg, stat = ier); REQUIRE(ier==0)
   end if
   allocate (cirreg(1:15,1:nbnd   ), stat = ier); REQUIRE(ier==0)

   ! setting coordinates of grid points

   do i = 0, xm+1
      x(i) = gox + i*h
   end do
   do j = 0, ym+1
      y(j) = goy + j*h
   end do
   do k = 0, zm+1
      z(k) = goz + k*h
   end do

   ! setting up index array:

   ! index = 1 for interior regular grid points;
   ! index = 2 for interior irregular grid points;
   ! index = 3 for boundary grid points right on the interface;
   ! index = 4 for exterior irregular grid points;
   ! index = 5 for exterior regular grid points;

   do k = 1, zm
   do j = 1, ym
   do i = 1, xm
      if ( insas(i,j,k) > 0 ) then
         index(i,j,k) = 1
      else
         index(i,j,k) = 5
      end if
   end do
   end do
   end do

   index2 = 0

   do ii = 1, nbnd
      i = iepsav(1,ii); j = iepsav(2,ii) ; k = iepsav(3,ii)
      index2(i,j,k) = ii
      if ( phi(i,j,k) > 0 ) then
         index(i,j,k) = 4
      else if ( phi(i,j,k) < 0 ) then
         index(i,j,k) = 2
      else
         index(i,j,k) = 3
      end if
   end do

end subroutine irreg_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine indexg here]
subroutine num_irreg_grd(xm,ym,zm,nbnd,h,iepsav,phi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! set up irregular grid point index and data structure for later use including:
! 1) projection point on the interface
! 2) normal vector
! 3) transformation matrix to local frame
! 4) curvature
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none
   integer xm,ym,zm
   integer nbnd
   _REAL_ h
   integer iepsav(4,xm*ym*zm)
   _REAL_ phi(0:xm+1,0:ym+1,0:zm+1)

   integer l,i,j,k,i0,j0,k0
   _REAL_ x1,y1,z1,px,py,pz,xyy,xzz,xyz
   _REAL_ t(3,3)

   do l = 1, nbnd
      i = iepsav(1,l); j = iepsav(2,l) ; k = iepsav(3,l)

         ! find the projection point of the irregular grid point

         call project(xm,ym,zm,h,h,h,h,i,j,k,x,y,z,phi,x1,y1,z1)

         ! find the local coordinate transformation matrix
         ! this is the nearest grid point to the projection point

         i0 = nint((x1 - x(0))/h)
         j0 = nint((y1 - y(0))/h)
         k0 = nint((z1 - z(0))/h)
         call transf(xm,ym,zm,h,h,h,h,x1,y1,z1,i0,j0,k0,x,y,z,phi,t)

         ! find the curvertures at the projection

         call curv(xm,ym,zm,h,h,h,h,x1,y1,z1,i0,j0,k0,x,y,z,phi,t,xyy,xzz,xyz)

         ! saving all for later

         cirreg(1 , l) = x1
         cirreg(2 , l) = y1
         cirreg(3 , l) = z1
         cirreg(4 , l) = xyy
         cirreg(5 , l) = xzz
         cirreg(6 , l) = xyz
         cirreg(7 , l) = t(1,1)
         cirreg(8 , l) = t(1,2)
         cirreg(9 , l) = t(1,3)
         cirreg(10, l) = t(2,1)
         cirreg(11, l) = t(2,2)
         cirreg(12, l) = t(2,3)
         cirreg(13, l) = t(3,1)
         cirreg(14, l) = t(3,2)
         cirreg(15, l) = t(3,3)

   end do

end subroutine num_irreg_grd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Single atom revised density function
_REAL_ function density_atom(dist)

   implicit none

   _REAL_ dist

   integer m
   _REAL_, parameter :: dash(6) = (/0.00d0,0.20d0,0.40d0,0.60d0,0.80d0,1.00d0 /)
   _REAL_, parameter :: spcoef(5,4) = &
       reshape((/1.000000  ,0.2100000  ,0.1500000  ,0.0500000  ,0.010000   ,  &
                -4.527143  ,-2.067608  ,0.0475730  ,-0.522686  ,-0.056828  ,  &
                -3.640532  ,15.938209  ,-5.362303  ,2.5110050  ,-0.181716  ,  &
                32.631235  ,-35.500854 ,13.122180  ,-4.487867  ,1.079289/), (/5,4/))
!   _REAL_ spcoef(4,4) , dash(5)
!   data spcoef /1.000000,0.1000000,6.4999998E-02,2.9999999E-02, &
!                -4.527143,-1.745714,0.2900000,-0.2542857,       &
!                0.0000000,11.12571,-2.982857,0.8057143,         &
!                14.83429,-18.81143,5.051429,-1.074286 /
!   data dash /0.d0,0.25d0,0.5d0,0.75d0,1.d0 /

   ! The distance is in the unit of solvent probe diameter

   density_atom = 0.0d0
   if ( dist > 1.d0 ) then
   else if ( dist <= 0.d0 ) then
      density_atom = 1.0d0 - 4.527143d0 * dist
   else
      do m = 1, 5
         if ( dist > dash(m) .and. dist <= dash(m+1) ) then
            density_atom = density_atom   + &
            spcoef(m,1)                   + &
            spcoef(m,2)*(dist-dash(m))    + &
            spcoef(m,3)*(dist-dash(m))**2 + &
            spcoef(m,4)*(dist-dash(m))**3
         end if
      end do
   end if

end function density_atom
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Single atom revised density function, the derivative
_REAL_ function density_deriv_atom(dist,dist0,dx)

   implicit none

   _REAL_ dist, dist0, dx

   integer m
   _REAL_, parameter :: dash(6) = (/0.00d0,0.20d0,0.40d0,0.60d0,0.80d0,1.00d0 /)
   _REAL_, parameter :: spcoef(5,4) = &
       reshape((/1.000000  ,0.2100000  ,0.1500000  ,0.0500000  ,0.010000   ,  &
                -4.527143  ,-2.067608  ,0.0475730  ,-0.522686  ,-0.056828  ,  &
                -3.640532  ,15.938209  ,-5.362303  ,2.5110050  ,-0.181716  ,  &
                32.631235  ,-35.500854 ,13.122180  ,-4.487867  ,1.079289/), (/5,4/))

   ! The distance is in the unit of solvent probe diameter

   density_deriv_atom = 0.0d0
   if ( dist > 1.d0 ) then
   else if ( dist <= 0.d0 ) then
      density_deriv_atom = 4.527143d0 * dx/dist0
   else
      do m = 1, 5
         if ( dist > dash(m) .and. dist <= dash(m+1) ) then
            density_deriv_atom = density_deriv_atom          - &
                    spcoef(m,2)                   * dx/dist0 - &
            2.0d0 * spcoef(m,3)*(dist-dash(m))    * dx/dist0 - &
            3.0d0 * spcoef(m,4)*(dist-dash(m))**2 * dx/dist0
         end if
      end do
   end if

end function density_deriv_atom

end module density_surface
