      subroutine fill_charge_grid(
     $         numatoms,charge,theta1,theta2,theta3,fr1,fr2,fr3,
     $         order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
c---------------------------------------------------------------------
c INPUT:
c      numatoms:  number of atoms
c      charge: the array of atomic charges
c      theta1,theta2,theta3: the spline coeff arrays
c      fr1,fr2,fr3 the scaled and shifted fractional coords
c      nfft1,nfft2,nfft3: the charge grid dimensions
c      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims
c      order: the order of spline interpolation
c OUTPUT:
c      Q the charge grid
c---------------------------------------------------------------------
      implicit none
      integer numatoms,order,nfft1,nfft2,nfft3
      integer nfftdim1,nfftdim2,nfftdim3
      double precision fr1(numatoms),fr2(numatoms),fr3(numatoms)
      double precision theta1(order,numatoms),theta2(order,numatoms),
     $     theta3(order,numatoms),charge(numatoms)
      double precision Q(2,nfftdim1,nfftdim2,nfftdim3)

      integer n,ntot,ith1,ith2,ith3,i0,j0,k0,i,j,k
      double precision prod
      ntot = 2*nfftdim1*nfftdim2*nfftdim3
      call clearQ(Q,ntot)

      do 300 n = 1,numatoms
        k0 = int(fr3(n)) - order
        do 200 ith3 = 1,order
         k0 = k0 + 1
         k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
         j0 = int(fr2(n)) - order
         do 150 ith2 = 1,order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
          prod = theta2(ith2,n)*theta3(ith3,n)*charge(n)
          i0 = int(fr1(n)) - order
          do 100 ith1 = 1,order
           i0 = i0 + 1
           i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
           Q(1,i,j,k) = Q(1,i,j,k) + theta1(ith1,n) * prod
100       continue
150      continue
200     continue
300   continue
      return
      end
c-----------------------------------------------------------
      subroutine clearQ(Q,ntot)
      integer ntot
      double precision Q(ntot)
      integer i
      do 10 i = 1,ntot
        Q(i) = 0.d0
10    continue
      return
      end
c-----------------------------------------------------------
      subroutine scalar_sum(
     $         Q,ewaldcof,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,
     $         nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,eer,vir)
      implicit none

      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      double precision Q(2,nfftdim1,nfftdim2,nfftdim3)
      double precision bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3),ewaldcof,volume
      double precision eer,vir(6)
      double precision recip(3,3)

      double precision pi,fac,denom,eterm,vterm,energy
      integer k,k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
      integer nf1,nf2,nf3
      double precision mhat1,mhat2,mhat3,msq,struc2

      indtop = nfft1*nfft2*nfft3
      pi = 3.14159265358979323846
      fac = pi**2/ewaldcof**2
      nff = nfft1*nfft2
      nf1 = nfft1/2
      if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
      nf2 = nfft2/2
      if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
      nf3 = nfft3/2
      if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1
      energy = 0.d0
      do 10 k = 1,6
        vir(k) = 0.d0
10    continue
      do 100 ind = 1,indtop-1

c get k1,k2,k3 from the relationship
c           ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1

       k3 = ind/nff + 1
       jnd = ind - (k3-1)*nff
       k2 = jnd/nfft1 + 1
       k1 = jnd - (k2-1)*nfft1 +1
       m1 = k1 - 1
       if ( k1 .gt. nf1 )m1 = k1 - 1 - nfft1
       m2 = k2 - 1
       if ( k2 .gt. nf2 )m2 = k2 - 1 - nfft2
       m3 = k3 - 1
       if ( k3 .gt. nf3 )m3 = k3 - 1 - nfft3
       mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
       mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
       mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
       msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
       denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k2)*bsp_mod3(k3)*msq
       eterm = exp(-fac*msq)/denom
       vterm = 2.d0*(fac*msq + 1.d0)/msq
       struc2 = Q(1,k1,k2,k3)**2 + Q(2,k1,k2,k3)**2
       energy = energy + eterm * struc2
       vir(1) = vir(1) + eterm * struc2 * (vterm*mhat1*mhat1 - 1.d0)
       vir(2) = vir(2) + eterm * struc2 * (vterm*mhat1*mhat2)
       vir(3) = vir(3) + eterm * struc2 * (vterm*mhat1*mhat3)
       vir(4) = vir(4) + eterm * struc2 * (vterm*mhat2*mhat2 - 1.d0)
       vir(5) = vir(5) + eterm * struc2 * (vterm*mhat2*mhat3)
       vir(6) = vir(6) + eterm * struc2 * (vterm*mhat3*mhat3 - 1.d0)
       Q(1,k1,k2,k3) = eterm * Q(1,k1,k2,k3)
       Q(2,k1,k2,k3) = eterm * Q(2,k1,k2,k3)
100   continue
      eer = 0.5d0 * energy
      do 150 k = 1,6
       vir(k) = 0.5d0*vir(k)
150   continue
      return
      end
c-----------------------------------------------------------
      subroutine grad_sum(
     $         numatoms,charge,recip,theta1,theta2,theta3,
     $         dtheta1,dtheta2,dtheta3,fx,fy,fz,fr1,fr2,fr3,
     $         order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
      implicit none
      integer numatoms,order,nfft1,nfft2,nfft3
      integer nfftdim1,nfftdim2,nfftdim3
      double precision recip(3,3)
      double precision fr1(numatoms),fr2(numatoms),fr3(numatoms)
      double precision fx(numatoms),fy(numatoms),fz(numatoms)
      double precision theta1(order,numatoms),theta2(order,numatoms),
     $     theta3(order,numatoms),charge(numatoms)
      double precision dtheta1(order,numatoms),dtheta2(order,numatoms),
     $     dtheta3(order,numatoms)
      double precision Q(2,nfftdim1,nfftdim2,nfftdim3)

      integer n,ith1,ith2,ith3,i0,j0,k0,i,j,k
      double precision f1,f2,f3,term

C$DOACROSS LOCAL(f1,f2,f3,k0,k,j0,j,i0,i,term,n,ith1,ith2,ith3),
C$&  SHARE(numatoms,fr1,fr2,fr3,charge,Q,fx,fy,fz,recip,order,
C$&   nfft1,nfft2,nfft3,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
      do 400 n = 1,numatoms
        f1 = 0.d0
        f2 = 0.d0
        f3 = 0.d0
        k0 = int(fr3(n)) - order
        do 200 ith3 = 1,order
         k0 = k0 + 1
         k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
         j0 = int(fr2(n)) - order
         do 150 ith2 = 1,order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
          i0 = int(fr1(n)) - order
          do 100 ith1 = 1,order
           i0 = i0 + 1
           i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
           term = charge(n)*Q(1,i,j,k)

c AARON: here is the change needed--- note you need to define an array to
c hold the electrostatic potentials
c assume esp(n) is initted to zero like f1-f3
c note I am putting new code in comments-- no time to make it legal
c          esp(n) = esp(n) + 
c            theta1(ith1,n)*theta2(ith2,n)*theta3(ith3,n)*Q(1,i,j,k)
c
c end of new code

c force is negative of grad
           f1 = f1 - nfft1 * term * dtheta1(ith1,n) *
     $          theta2(ith2,n) * theta3(ith3,n)
           f2 = f2 - nfft2 * term * theta1(ith1,n) *
     $          dtheta2(ith2,n) * theta3(ith3,n)
           f3 = f3 - nfft3 * term * theta1(ith1,n) *
     $          theta2(ith2,n) * dtheta3(ith3,n)
100       continue
150      continue
200     continue
        fx(n) = fx(n) + recip(1,1)*f1+recip(1,2)*f2+recip(1,3)*f3
        fy(n) = fy(n) + recip(2,1)*f1+recip(2,2)*f2+recip(2,3)*f3
        fz(n) = fz(n) + recip(3,1)*f1+recip(3,2)*f2+recip(3,3)*f3
400   continue
      return
      end
c-----------------------------------------------------------
