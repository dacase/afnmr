      subroutine fill_charge_grid(
     $         numatoms,charge,theta1,theta2,theta3,fr1,fr2,fr3,
     $         order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q,
     $         is,js,ks,ind)
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
      double precision theta1(numatoms,order),theta2(numatoms,order),
     $     theta3(numatoms,order),charge(numatoms)
      double precision fr1(numatoms),fr2(numatoms),fr3(numatoms)
      double precision Q(2,nfftdim1*nfftdim2*nfftdim3)
      double precision is(numatoms),js(numatoms),ks(numatoms)
      integer ind(numatoms)

      integer n,ntot,ith1,ith2,ith3,i0,j0,k0,i,j,k
      double precision xi,xj,xk,xith1,xith2,xith3,xnfft1,
     $     xnfft2,xnfft3,xnfftdim1,xnfftdim2,xnfftdim3,xk0,xj0,xi0
      ntot = 2*nfftdim1*nfftdim2*nfftdim3
      call clearQ(Q,ntot)

      xnfft1 = float(nfft1)
      xnfft2 = float(nfft2)
      xnfft3 = float(nfft3)
      xnfftdim1 = float(nfftdim1)
      xnfftdim2 = float(nfftdim2)
      xnfftdim3 = float(nfftdim3)
      do 50 n = 1,numatoms
       ks(n) = float(int(fr3(n)) - order)
       js(n) = float(int(fr2(n)) - order)
       is(n) = float(int(fr1(n)) - order)
50    continue
      do 400 ith3 = 1,order
       xith3 = float(ith3)
       do 300 ith2 = 1,order
        xith2 = float(ith2)
        do 200 ith1 = 1,order
         xith1 = float(ith1)
         do 100 n = 1,numatoms
          xk0 = ks(n) + xith3
          xj0 = js(n) + xith2
          xi0 = is(n) + xith1
          xk = xk0 + 1.0 + 0.5 * (xnfft3 - sign(xnfft3,xk0))
          xj = xj0 + 1.0 + 0.5 * (xnfft2 - sign(xnfft2,xj0))
          xi = xi0 + 1.0 + 0.5 * (xnfft1 - sign(xnfft1,xi0))
          ind(n) = int(xi + (xj-1.0)*xnfftdim1 + 
     $          (xk-1.0)*xnfftdim1*xnfftdim2)
100      continue
         do 150 n = 1,numatoms
          Q(1,ind(n)) = Q(1,ind(n)) + theta1(n,ith1) * 
     $           theta2(n,ith2)*theta3(n,ith3)*charge(n)
150      continue
200     continue
300    continue
400   continue
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
      integer k,k1,k2,k3,nff,ind,jnd,indtop,ii
      integer nf1,nf2,nf3
      double precision mhat1,mhat2,mhat3,msq,struc2
      integer k1s(1000),k2s(1000),k3s(1000)
      double precision m1s(1000),m2s(1000),m3s(1000),m1,m2,m3
      integer indlo,indhi,loop,numloops

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
      numloops = (indtop - 1)/512 + 1
      do 300 loop = 1,numloops
       indlo = 512*(loop - 1) + 1
       indhi = 512*loop
       if ( loop .eq. numloops )indhi = indtop - 1
       ii = 0
       do 100 ind = indlo,indhi
        ii = ii + 1

c get k1,k2,k3 from the relationship
c           ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1

        k3 = ind/nff + 1
        jnd = ind - (k3-1)*nff
        k2 = jnd/nfft1 + 1
        k1 = jnd - (k2-1)*nfft1 +1
        m1 = float(k1 - 1)
        if ( k1 .gt. nf1 )m1 = float(k1 - 1 - nfft1)
        m2 = float(k2 - 1)
        if ( k2 .gt. nf2 )m2 = float(k2 - 1 - nfft2)
        m3 = float(k3 - 1)
        if ( k3 .gt. nf3 )m3 = float(k3 - 1 - nfft3)
        k1s(ii) = k1
        k2s(ii) = k2
        k3s(ii) = k3
        m1s(ii) = m1
        m2s(ii) = m2
        m3s(ii) = m3
100    continue
       ii = 0
CDIR$ IVDEP
       do 200 ind = indlo,indhi
        ii = ii + 1
        k1 = k1s(ii)
        k2 = k2s(ii)
        k3 = k3s(ii)
        m1 = m1s(ii)
        m2 = m2s(ii)
        m3 = m3s(ii)
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
200    continue
300   continue
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
     $         order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q,
     $         is,js,ks,ind)
      implicit none
      integer numatoms,order,nfft1,nfft2,nfft3
      integer nfftdim1,nfftdim2,nfftdim3
      double precision recip(3,3)
      double precision fr1(numatoms),fr2(numatoms),fr3(numatoms)
      double precision fx(numatoms),fy(numatoms),fz(numatoms)
      double precision theta1(numatoms,order),theta2(numatoms,order),
     $     theta3(numatoms,order),charge(numatoms)
      double precision dtheta1(numatoms,order),dtheta2(numatoms,order),
     $     dtheta3(numatoms,order)
      double precision Q(2,nfftdim1*nfftdim2*nfftdim3)
      double precision is(numatoms),js(numatoms),ks(numatoms)
      integer ind(numatoms)

      integer n,ith1,ith2,ith3,i0,j0,k0,i,j,k
      double precision f1,f2,f3,term
      double precision xi,xj,xk,xith1,xith2,xith3,xnfft1,
     $     xnfft2,xnfft3,xnfftdim1,xnfftdim2,xnfftdim3,xk0,xj0,xi0

      xnfft1 = float(nfft1)
      xnfft2 = float(nfft2)
      xnfft3 = float(nfft3)
      xnfftdim1 = float(nfftdim1)
      xnfftdim2 = float(nfftdim2)
      xnfftdim3 = float(nfftdim3)
      do 50 n = 1,numatoms
       ks(n) = float(int(fr3(n)) - order)
       js(n) = float(int(fr2(n)) - order)
       is(n) = float(int(fr1(n)) - order)
50    continue
      do 400 ith3 = 1,order
       xith3 = float(ith3)
       do 300 ith2 = 1,order
        xith2 = float(ith2)
        do 200 ith1 = 1,order
         xith1 = float(ith1)
         do 100 n = 1,numatoms
          xk0 = ks(n) + xith3
          xj0 = js(n) + xith2
          xi0 = is(n) + xith1
          xk = xk0 + 1.0 + 0.5 * (xnfft3 - sign(xnfft3,xk0))
          xj = xj0 + 1.0 + 0.5 * (xnfft2 - sign(xnfft2,xj0))
          xi = xi0 + 1.0 + 0.5 * (xnfft1 - sign(xnfft1,xi0))
          ind(n) = int(xi + (xj-1.0)*xnfftdim1 + 
     $          (xk-1.0)*xnfftdim1*xnfftdim2)
100      continue
         do 150 n = 1,numatoms
          term = charge(n)*Q(1,ind(n))
          f1 = nfft1 * term * dtheta1(n,ith1) *
     $         theta2(n,ith2) * theta3(n,ith3)
          f2 = nfft2 * term * theta1(n,ith1) *
     $         dtheta2(n,ith2) * theta3(n,ith3)
          f3 = nfft3 * term * theta1(n,ith1) *
     $         theta2(n,ith2) * dtheta3(n,ith3)
c force is negative of grad
          fx(n) = fx(n) - (recip(1,1)*f1+recip(1,2)*f2+recip(1,3)*f3)
          fy(n) = fy(n) - (recip(2,1)*f1+recip(2,2)*f2+recip(2,3)*f3)
          fz(n) = fz(n) - (recip(3,1)*f1+recip(3,2)*f2+recip(3,3)*f3)
150      continue
200     continue
300    continue
400   continue
      return
      end
c-----------------------------------------------------------
