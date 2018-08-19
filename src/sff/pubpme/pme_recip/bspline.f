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

      double precision w
      integer n

      do 100 n = 1,numatoms
        w = fr1(n)-int(fr1(n))
        call fill_bspline(w,order,theta1(1,n),dtheta1(1,n))
        w = fr2(n)-int(fr2(n))
        call fill_bspline(w,order,theta2(1,n),dtheta2(1,n))
        w = fr3(n)-int(fr3(n))
        call fill_bspline(w,order,theta3(1,n),dtheta3(1,n))
100   continue
      return
      end
c---------------------------------------------------
      subroutine fill_bspline(w,order,array,darray)
c---------- use standard B-spline recursions: see doc file
      implicit none
      integer order
      double precision w,array(order),darray(order)

      integer k
c do linear case
      call init(array,w,order)
c compute standard b-spline recursion
      do 10 k = 3,order-1
       call one_pass(array,w,k)
10    continue
c perform standard b-spline differentiation
      call diff(array,darray,order)
c one more recursion
      call one_pass(array,w,order)
      return
      end
c---------------------------------------------------
      subroutine init(c,x,order)
      implicit none
      integer order
      double precision c(order),x
      c(order) = 0.d0
      c(2) = x
      c(1) = 1.d0 - x
      return
      end
c-------------------------------------
      subroutine one_pass(c,x,k)
      implicit none
      double precision c(*),x
      integer k

      double precision div
      integer j

      div = 1.d0 / (k-1)
      c(k) = div*x*c(k-1)
      do 100 j = 1,k-2
       c(k-j) = div*((x+j)*c(k-j-1) + (k-j-x)*c(k-j))
100   continue
      c(1) = div*(1-x)*c(1)
      return
      end
c-------------------------------------
      subroutine diff(c,d,order)
      implicit none
      double precision c(*),d(*)
      integer order

      integer j
      d(1) = -c(1)
      do 10 j = 2,order
       d(j) = c(j-1) - c(j)
10    continue
      return
      end
c-------------------------------------
