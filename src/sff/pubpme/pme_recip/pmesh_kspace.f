      subroutine pmesh_kspace_get_sizes(
     $     nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack)
      implicit none
      integer nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack

c INPUT  
c      nfft1,nfft2,nfft3,numatoms,order
c      nfft1,nfft2,nfft3 are the dimensions of the charge grid array
c      numatoms is number of atoms
c      order is the order of B-spline interpolation

c OUTPUT
c      sizfftab,sizffwrk,siztheta,siz_Q
c      sizfftab is permanent 3d fft table storage
c      sizffwrk is temporary 3d fft work storage
c      siztheta is size of arrays theta1-3 dtheta1-3
c      sizheap is total size of permanent storage
c      sizstack is total size of temporary storage


c This routine computes the above output parameters needed for 
c heap or stack allocation.

      integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork

      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,
     $       sizfftab,sizffwrk)
      siztheta = numatoms*order
      siz_Q = 2*nfftdim1*nfftdim2*nfftdim3
      sizheap = nfft1+nfft2+nfft3+sizfftab
      sizstack = siz_Q+6*siztheta+sizffwrk+3*numatoms
c     write(6,*)'total HEAP storage needed = ',sizheap
c     write(6,*)'total STACK storage needed = ',sizstack

      return
      end
c----------------------------------------------------
      subroutine pmesh_kspace_setup(
     $    bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork,
     $    nfft1,nfft2,nfft3,order,sizfftab,sizffwrk)
      implicit none

c  see DO_PMESH_KSPACE for explanation of arguments

      integer nfft1,nfft2,nfft3,order,sizfftab,sizffwrk
      double precision bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3)
      double precision fftable(sizfftab),ffwork(sizffwrk)
   
      double precision dummy
      integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw

      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
      call load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,
     $   nfft1,nfft2,nfft3,order)
      call fft_setup(dummy,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
      return
      end
c----------------------------------------------------
      subroutine do_pmesh_kspace(
     $   numatoms,x,y,z,charge,recip,volume,ewald_coeff,
     $   order,nfft1,nfft2,nfft3,
     $   eer,dx,dy,dz,virial,
     $   sizfftab,sizffwrk,siztheta,siz_Q,
     $   bsp_mod1,bsp_mod2,bsp_mod3,fftable,
     $   Q,ffwork,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,
     $   fr1,fr2,fr3,time)
      implicit none

c INPUT 
c       numatoms:  number of atoms
c       x,y,z:   atomic coords
c       charge  atomic charges
c       recip: 3x3 array of reciprocal unit cell vectors (stored as columns)
c       volume: the volume of the unit cell
c       ewald_coeff:   ewald convergence parameter
c       order: the order of Bspline interpolation. E.g. cubic is order 4
c          fifth degree is order 6 etc. The order must be an even number 
c          and at least 4.
c       nfft1,nfft2,nfft3: the dimensions of the charge grid array
      integer numatoms,order,nfft1,nfft2,nfft3
      double precision x(numatoms),y(numatoms),z(numatoms),
     $       charge(numatoms),recip(3,3),volume,ewald_coeff

c OUTPUT
c       eer:  ewald reciprocal or k-space  energy
c       dx,dy,dz: forces incremented by k-space sum
c       virial:  virial due to k-space sum (valid for atomic scaling;
c                rigid molecule virial needs a correction term not
c                computed here
c       time: used to profile the different component routines
      double precision eer,dx(numatoms),dy(numatoms),dz(numatoms),
     $        virial(6),time(5)

c SIZES of some arrays
      integer   sizfftab,sizffwrk,siztheta,siz_Q

c HEAP STORAGE:  These arrays need to be preserved throughout simulation
      double precision bsp_mod1(nfft1),bsp_mod2(nfft2),
     $                 bsp_mod3(nfft3),fftable(sizfftab)
c STACK STORAGE: These arrays can be tossed after leaving this routine
      double precision Q(siz_Q),ffwork(sizffwrk),theta1(siztheta),
     $          theta2(siztheta),theta3(siztheta),dtheta1(siztheta),
     $          dtheta2(siztheta),dtheta3(siztheta),fr1(numatoms),
     $          fr2(numatoms),fr3(numatoms)

      integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw
      real tim1,tim2
c  get some integer array dimensions
      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
      
      call get_scaled_fractionals(
     $         numatoms,x,y,z,recip,nfft1,nfft2,nfft3,fr1,fr2,fr3)
      call second(tim1)
      call get_bspline_coeffs(
     $         numatoms,fr1,fr2,fr3,order,
     $         theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
      call second(tim2)
      time(1) = time(1) + tim2-tim1
      tim1 = tim2
      call fill_charge_grid(
     $         numatoms,charge,theta1,theta2,theta3,fr1,fr2,fr3,order,
     $         nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
      call second(tim2)
      time(2) = time(2) + tim2-tim1
      tim1 = tim2
      call fft_back(
     $         Q,fftable,ffwork,nfft1,nfft2,nfft3,
     $         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)
      call second(tim2)
      time(3) = time(3) + tim2-tim1
      tim1 = tim2
      call scalar_sum(
     $         Q,ewald_coeff,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,
     $         nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,eer,virial)
      call second(tim2)
      time(4) = time(4) + tim2-tim1
      tim1 = tim2
      call fft_forward(
     $         Q,fftable,ffwork,nfft1,nfft2,nfft3,
     $         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)
      call second(tim2)
      time(3) = time(3) + tim2-tim1
      tim1 = tim2
      call grad_sum(
     $         numatoms,charge,recip,theta1,theta2,theta3,
     $         dtheta1,dtheta2,dtheta3,dx,dy,dz,fr1,fr2,fr3,
     $         order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
      call second(tim2)
      time(5) = time(5) + tim2-tim1
      return
      end
c----------------------------------------------------------------------
      subroutine get_scaled_fractionals(
     $           numatoms,x,y,z,recip,nfft1,nfft2,nfft3,
     $           fr1,fr2,fr3)
      implicit none

c INPUT:
c      numatoms: number of atoms
c      x,y,z: arrays of cartesian coords
c      recip: the 3x3 array of reciprocal vectors stored as columns
c OUTPUT:
c     fr1,fr2,fr3 the scaled and shifted fractional coords

      integer numatoms,nfft1,nfft2,nfft3
      double precision x(numatoms),y(numatoms),z(numatoms),recip(3,3)
      double precision fr1(numatoms),fr2(numatoms),fr3(numatoms)

      integer n
      double precision w
      do 100 n = 1,numatoms
        w = x(n)*recip(1,1)+y(n)*recip(2,1)+z(n)*recip(3,1)
        fr1(n) = nfft1*(w - anint(w) + 0.5d0)
        w = x(n)*recip(1,2)+y(n)*recip(2,2)+z(n)*recip(3,2)
        fr2(n) = nfft2*(w - anint(w) + 0.5d0)
        w = x(n)*recip(1,3)+y(n)*recip(2,3)+z(n)*recip(3,3)
        fr3(n) = nfft3*(w - anint(w) + 0.5d0)
100   continue
      return
      end
c---------------------------------------------------------------
      subroutine load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,
     $   nfft1,nfft2,nfft3,order)
      implicit none
      double precision bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3)
      integer nfft1,nfft2,nfft3,order

      integer MAXORDER
      parameter (MAXORDER=25)
      integer MAXNFFT
      parameter (MAXNFFT=1000)
      double precision array(MAXORDER),darray(MAXORDER),w
      double precision bsp_arr(MAXNFFT)
      integer i,maxn

c this routine loads the moduli of the inverse DFT of the B splines
c bsp_mod1-3 hold these values, nfft1-3 are the grid dimensions,
c Order is the order of the B spline approx.

      if ( order .gt. MAXORDER )then
c      write(6,*)'order too large! check on MAXORDER'
       stop
      endif
      maxn = max(nfft1,nfft2,nfft3)
      if ( maxn .gt. MAXNFFT )then 
c      write(6,*)'nfft1-3 too large! check on MAXNFFT'
       stop
      endif
      w = 0.d0
      call fill_bspline(w,order,array,darray)
      do 100 i = 1,maxn
        bsp_arr(i) = 0.d0
100   continue
      do 150 i = 2,order+1
       bsp_arr(i) = array(i-1)
150   continue
      call DFTMOD(bsp_mod1,bsp_arr,nfft1)
      call DFTMOD(bsp_mod2,bsp_arr,nfft2)
      call DFTMOD(bsp_mod3,bsp_arr,nfft3)
      return
      end
c------------------------------------------------------------------------
      subroutine DFTMOD(bsp_mod,bsp_arr,nfft)
      implicit none
      integer nfft
      double precision bsp_mod(nfft),bsp_arr(nfft)
c Computes the modulus of the discrete fourier transform of bsp_arr,
c  storing it into bsp_mod

      integer j,k
      double precision sum1,sum2,twopi,arg,tiny
      twopi = 2.d0*3.14159265358979323846
      tiny = 1.d-7
      do 300 k = 1,nfft
       sum1 = 0.d0
       sum2 = 0.d0
       do 250 j = 1,nfft
         arg = twopi*(k-1)*(j-1)/nfft
         sum1 = sum1 + bsp_arr(j)*cos(arg)
         sum2 = sum2 + bsp_arr(j)*sin(arg)
250    continue
       bsp_mod(k) = sum1**2 + sum2**2
300   continue
      do 400 k = 1,nfft
       if ( bsp_mod(k) .lt. tiny )
     $     bsp_mod(k) = 0.5d0*(bsp_mod(k-1) + bsp_mod(k+1))
400   continue
      return
      end
c------------------------------------------------------------------------
