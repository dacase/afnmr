*****************************************************************************
*
*	3D (slow) Fourier Transform
*   this 1d->3d code is brute force approach
*   the 1d code is a double precision version of fftpack from netlib
*   due to Paul N Swartztrauber at NCAR Boulder Coloraso
*
*****************************************************************************

      subroutine pubz3di(n1,n2,n3,table,ntable)
      implicit none
      integer n1,n2,n3,ntable
      double precision table(ntable,3)
c ntable should be 4*max(n1,n2,n3) +15


      call cffti(n1,table(1,1))
      call cffti(n2,table(1,2))
      call cffti(n3,table(1,3))

      return
      end
*****************************************************************************
      subroutine pubz3d(isign,n1,n2,n3,w,ld1,ld2,table,ntable,
     $    work,nwork)
      implicit none

      integer n1,n2,n3,ld1,ld2,isign,ntable,nwork
      double complex w(ld1,ld2,n3)
      double complex work( nwork)
      double precision table(ntable,3)

      integer i,j,k
c ntable should be 4*max(n1,n2,n3) +15
c nwork should be max(n1,n2,n3)
c
c   transform along X  first ...
c
      do 100 k = 1, n3
       do 90 j = 1, n2
        do 70 i = 1,n1
          work(i) = w(i,j,k)
70      continue
        if ( isign .eq. -1) call cfftf(n1,work,table(1,1))
        if ( isign .eq. 1) call cfftb(n1,work,table(1,1))
        do 80 i = 1,n1
          w(i,j,k) = work(i)
80      continue
90     continue
100   continue
c
c   transform along Y then ...
c
      do 200 k = 1,n3
       do 190 i = 1,n1
        do 170 j = 1,n2
          work(j) = w(i,j,k)
170     continue
        if ( isign .eq. -1) call cfftf(n2,work,table(1,2))
        if ( isign .eq. 1) call cfftb(n2,work,table(1,2))
        do 180 j = 1,n2
          w(i,j,k) = work(j)
180     continue
190    continue
200   continue
c
c   transform along Z finally ...
c
      do 300 i = 1, n1
       do 290 j = 1, n2
        do 270 k = 1,n3
          work(k) = w(i,j,k)
270     continue
        if ( isign .eq. -1) call cfftf(n3,work,table(1,3))
        if ( isign .eq. 1) call cfftb(n3,work,table(1,3))
        do 280 k = 1,n3
          w(i,j,k) = work(k)
280     continue
290    continue
300   continue

      return
      end
c----------------------------------------------------
