c---------------------------------------------------------------------
      subroutine get_bspline_coeffs(
     $           numatoms,fr1,fr2,fr3,order,
     $           theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
c---------------------------------------------------------------------
c INPUT:
c      numatoms: number of atoms
c      fr1,fr2,fr3 the scaled and shifted fractional coords
c      order: the order of spline interpolation
c OUTPUT
c      theta1,theta2,theta3: the spline coeff arrays
c      dtheta1,dtheta2,dtheta3: the 1st deriv of spline coeff arrays
c---------------------------------------------------------------------
      implicit none
      integer numatoms,order
      double precision fr1(numatoms),fr2(numatoms),fr3(numatoms)
      double precision theta1(order,numatoms),theta2(order,numatoms),
     $     theta3(order,numatoms),dtheta1(order,numatoms),
     $     dtheta2(order,numatoms),dtheta3(order,numatoms)

      call fill_bspline(fr1,order,theta1,dtheta1,numatoms)
      call fill_bspline(fr2,order,theta2,dtheta2,numatoms)
      call fill_bspline(fr3,order,theta3,dtheta3,numatoms)

      return
      end
c---------------------------------------------------
      subroutine fill_bspline(fr,order,array,darray,numatoms)
c---------- use standard B-spline recursions: see doc file
      implicit none
      integer order,numatoms
      double precision fr(numatoms),array(numatoms,order),
     $      darray(numatoms,order)

      integer k
c do linear case
      call init(array,fr,numatoms)
c compute standard b-spline recursion
      do 10 k = 3,order-1
       call one_pass(array,fr,k,numatoms)
10    continue
c perform standard b-spline differentiation
      call diff(array,darray,order,numatoms)
c one more recursion
      call one_pass(array,fr,order,numatoms)
      return
      end
c---------------------------------------------------
      subroutine init(c,fr,numatoms)
      implicit none
      integer numatoms
      double precision c(numatoms,*),fr(numatoms)

      double precision w
      integer n
      do 100 n = 1,numatoms
        w = fr(n)-int(fr(n))
        c(n,2) = w
        c(n,1) = 1.d0 - w
100   continue
      return
      end
c-------------------------------------
      subroutine one_pass(c,fr,k,numatoms)
      implicit none
      double precision c(numatoms,*),fr(numatoms)
      integer k,numatoms

      double precision div,w
      integer n,j

      div = 1.d0 / (k-1)
      do 100 n = 1,numatoms
        w = fr(n)-int(fr(n))
        c(n,k) = div*w*c(n,k-1)
100   continue
      do 300 j = 1,k-2
        do 200 n = 1,numatoms
          w = fr(n)-int(fr(n))
          c(n,k-j) = div*((w+j)*c(n,k-j-1) + (k-j-w)*c(n,k-j))
200     continue
300   continue
      do 400 n = 1,numatoms
        w = fr(n)-int(fr(n))
        c(n,1) = div*(1-w)*c(n,1)
400   continue
      return
      end
c-------------------------------------
      subroutine diff(c,d,order,numatoms)
      implicit none
      double precision c(numatoms,order),d(numatoms,order)
      integer order,numatoms

      integer j,n
      do 100 n = 1,numatoms
        d(n,1) = -c(n,1)
100   continue
      do 300 j = 2,order
       do 200 n = 1,numatoms
          d(n,j) = c(n,j-1) - c(n,j)
200    continue
300   continue
      return
      end
c-------------------------------------
