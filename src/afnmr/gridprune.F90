program gridprune
      implicit none
      integer, parameter :: MAXGRID = 5000000
      integer, parameter :: MAXATOM = 10000
      integer, parameter :: MAXSURF = 12000

      double precision x(MAXGRID), y(MAXGRID), z(MAXGRID), c(MAXGRID), &
          mx(MAXATOM), my(MAXATOM), mz(MAXATOM), b(MAXSURF), dis(MAXATOM), &
          a(MAXATOM,MAXSURF), closest, s(MAXATOM), bcopy(MAXATOM), &
          acopy(MAXATOM,MAXSURF)
!     double precision pbc(MAXSURF)

      double precision closec, mmx, mmy, mmz, mmc, molcharge, bsum,svdprec
      integer i,j,ngrid,natom,nsurf,info,rank,nclose,lengthb
      character*80 arg, line
      logical solinprot

!     defaults:
      svdprec=1.d-3

      if( command_argument_count() .gt. 1 ) then
         write(0,*) &
       'Usage: gridprune [ -solinprot ]'
         stop
      end if

      if( command_argument_count() .eq. 1 ) then
         call get_command_argument( 1, arg, lengthb )
         if( arg(1:lengthb) == '-solinprot') then
            solinprot = .true.
         else
            write(0,*) 'Unknown argument: ',arg(1:lengthb)
            stop
         end if
      else
         solinprot = .false.
      end if

!     read the set of external charges (from, e.g. the rism3d.snglpnt 
!     output grid, or the pbsa "surface charges":
!     (For solinprot, these will be interpreted as the points and the
!     corresponding ProteinField + ReactionField potentials to be fit.)

      open( unit=13, file='input.xyzv', status='OLD' )
      do i=1,MAXGRID
         read(13,*, end=99) x(i),y(i),z(i),c(i)
      end do
      write(0,*) 'Error: more than ',MAXGRID,' grid points in input.xyzv'
      stop
   99 ngrid = i-1
      write(6,*) 'ngrid = ', ngrid
      close(13)

!     read the molecular coordinates (and charges):
!         (Note that the coordinates are ignored for solinprot, the
!          the molcharge value is still used for the constraints)
      j=0
      molcharge = 0.d0
      do i=1,MAXATOM
         read(5,'(a80)',end=199) line
         if( line(1:4) .eq. 'ATOM') then
            read(line,'(30x,4f8.3)') mmx, mmy, mmz, mmc
            j = j + 1
            mx(j) = mmx
            my(j) = mmy
            mz(j) = mmz
            molcharge = molcharge + mmc
         end if
      end do
      write(0,*) 'Error: more than ',MAXATOM,' atoms in input file'
      stop
  199 natom = j
      write(6,*) 'natom = ', natom, ', molcharge = ', molcharge

      if( solinprot ) then

!     put the target potential into the "b" array, and the corresponding
!     coordinates into mx,my,mz:

         natom = ngrid  ! for solinprot, target points become the "atoms"
         b(1:natom) = c(1:natom)
         mx(1:natom) = x(1:natom)
         my(1:natom) = y(1:natom)
         mz(1:natom) = z(1:natom)

      else

!     get closest atom distance for each grid point, then fill the b array
!     with the Coulomb potential at each atom center arising from the external
!     charges in the input.xyzv file:

         closec = 0.0
         nclose = 0
         do j=1,natom
            b(j) = 0.d0
         end do
         do i=1,ngrid
            closest = 1000.
            do j=1,natom
               dis(j)=sqrt((x(i)-mx(j))**2 + (y(i)-my(j))**2 + (z(i)-mz(j))**2)
               if( dis(j) < closest ) then
                  closest = dis(j)
               end if
            end do
            if (closest .gt. 0.5  .and. abs(c(i)) .gt. 0.00000008 )  then
               do j=1,natom
                  b(j) = b(j) + c(i)/dis(j)
               end do
            else
               closec = closec + c(i)
               nclose = nclose + 1
            end if
         end do
         write(6,'(a,f6.3,a,i6,a,i6)') 'closec = ', closec, &
               ', nclose = ', nclose, '/',ngrid

      end if ! (solinprot)

      bcopy(1:natom) = b(1:natom)

!     get the surface points (from ms):
!       (N.B.: this overwrites the previous values of x,y,z!)

      open(unit=10, file='surf.pos', status='OLD')
      do i=1,MAXSURF
         read(10,*,end=299) x(i),y(i),z(i)
      end do
      write(0,*) 'Error: more than ',MAXSURF,' points in surf.pos'
      stop
  299 nsurf = i-1
      write(6,*) 'nsurf = ', nsurf
      close(10)

!     construct the matrix of inverse distances between the atomic centers
!     and the MS surface points:

      do i=1,natom
         do j=1,nsurf
            a(i,j) = 1.0/sqrt( (x(j)-mx(i))**2 + (y(j)-my(i))**2 + &
                (z(j)-mz(i))**2 )
            acopy(i,j) = a(i,j)
         end do
      end do

!     add extra row for constraint on total surface charge:
      a(natom+1,1:nsurf) = 0.1d0
      b(natom+1) = -molcharge*0.1d0

!     solve the linear equations, A*c=b, using SVD:
      call dgelss( natom+1, nsurf, 1, a, MAXATOM, b, MAXSURF, &
          s, svdprec, rank, c, MAXGRID, info)
      write(6,*) 'dgelss returns ', info, rank
      if( info .ne. 0 ) stop

!     surface charges, for afnmr:
      bsum = sum(b(1:nsurf))
      write(6,*) 'bsum = ', bsum
      open( unit=11, file='srfchg.pos', status="REPLACE" )
      do i=1,nsurf
        write(11,'(3f10.3,f15.8)') x(i),y(i),z(i),b(i)
      end do
      close(11)

!     check accuracy of fit:
      call dgemm( 'N', 'N', natom, 1, nsurf, 1.d0, acopy, MAXATOM, &
         b, MAXSURF, 0.d0, c, MAXATOM )
      open( unit=12, file='check.dat', status="REPLACE" )
      do i=1,natom
        write(12,'(2f15.8)') bcopy(i), c(i)
      end do
      close(12)

!     compare fields at atoms for PB vs RISM:
!     call dgemm( 'N', 'N', natom, 1, nsurf, 1.d0, acopy, MAXATOM, &
!        pbc, MAXSURF, 0.d0, c, MAXATOM )
!     open( unit=12, file='pb_vs_rism.dat', status="REPLACE" )
!     do i=1,natom
!       write(12,'(2f15.8)') bcopy(i), c(i)
!     end do
!     close(12)

end program gridprune

