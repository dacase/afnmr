      implicit none
      integer MAXS,MAXN,MAXORD,MAXT
      parameter(MAXS=10000,MAXN=1000,MAXORD=14,MAXT=300000)
      double precision x(MAXS),y(MAXS),z(MAXS),cg(MAXS)
      double precision dir_ene,adj_ene,rec_ene,self_ene
      double precision dir_vir(6),adj_vir(6),rec_vir(6)
      double precision dfx(MAXS),dfy(MAXS),dfz(MAXS)
      double precision afx(MAXS),afy(MAXS),afz(MAXS)
      double precision rfx(MAXS),rfy(MAXS),rfz(MAXS)
      double precision fr1(MAXS),fr2(MAXS),fr3(MAXS)
      double precision fx(MAXS),fy(MAXS),fz(MAXS)
      double precision fx1(MAXS),fy1(MAXS),fz1(MAXS)
      double precision fx2(MAXS),fy2(MAXS),fz2(MAXS)
      double precision appdfx(MAXS),appdfy(MAXS),appdfz(MAXS)
      double precision bsp_mod1(MAXN),bsp_mod2(MAXN),bsp_mod3(MAXN)
      double precision box,cgh,cgo,factor
      double precision cutoff,tol,ewald_coeff
      double precision eigmin,reclng(3),recip(3,3),volume
      double precision theta1(MAXORD*MAXS),dtheta1(MAXORD*MAXS)
      double precision theta2(MAXORD*MAXS),dtheta2(MAXORD*MAXS)
      double precision theta3(MAXORD*MAXS),dtheta3(MAXORD*MAXS)
      double precision fftable(MAXN),ffwork(MAXT)
      double precision virial(6),eer,time(5),err
      integer numatoms,numwats,i,j,n,numreps
      integer nfft,nfft1,nfft2,nfft3,order
      integer sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack
      double precision Q(MAXT)
      real tim1,tim2


      factor = sqrt(332.17752)
      cgh = 0.417d0*factor
      cgo = -2.d0*cgh
      open(unit=8,file='small.pdb',status='old')
      read(8,25)box,numatoms
25    format(21x,f9.3,i6)
      write(6,*)'box = ',box,' numatoms = ',numatoms
      if ( numatoms .gt. maxs )then
        write(6,*)'Too many atoms'
        stop
      endif
      do 50 i = 1,numatoms
        fx2(i) = 0.d0
        fy2(i) = 0.d0
        fz2(i) = 0.d0
50    continue
      numwats = numatoms/3
      if ( numatoms .gt. MAXS )then
       write(6,*)'MAX too small'
       stop
      endif
      n = 0
      do 100 i = 1,numwats
       read(8,20)x(n+1),y(n+1),z(n+1)
       read(8,20)x(n+2),y(n+2),z(n+2)
       read(8,20)x(n+3),y(n+3),z(n+3)
       cg(n+1) = cgo
       cg(n+2) = cgh
       cg(n+3) = cgh
       n = n+3
100   continue
20    format(30x,3f8.3)
      open(unit=10,file='exact.bin',status='old',form='unformatted')
      read(10)numatoms
      read(10)cutoff,tol,ewald_coeff,box
      read(10)self_ene,dir_ene,adj_ene,rec_ene
      read(10)dir_vir,adj_vir,rec_vir
      read(10)(dfx(i),i=1,numatoms)
      read(10)(dfy(i),i=1,numatoms)
      read(10)(dfz(i),i=1,numatoms)
      read(10)(afx(i),i=1,numatoms)
      read(10)(afy(i),i=1,numatoms)
      read(10)(afz(i),i=1,numatoms)
      read(10)(rfx(i),i=1,numatoms)
      read(10)(rfy(i),i=1,numatoms)
      read(10)(rfz(i),i=1,numatoms)
      read(10)(appdfx(i),i=1,numatoms)
      read(10)(appdfy(i),i=1,numatoms)
      read(10)(appdfz(i),i=1,numatoms)
      close(10)
      do 35 i = 1,numatoms
       fx1(i) = dfx(i)+afx(i)+rfx(i)
       fy1(i) = dfy(i)+afy(i)+rfy(i)
       fz1(i) = dfz(i)+afz(i)+rfz(i)
35    continue
      write(6,*)'numatoms = ',numatoms
      write(6,*)'cutoff,tol,ewaldcof,box = ',
     $     cutoff,tol,ewald_coeff,box
      volume = box*box*box
      reclng(1) = box
      reclng(2) = box
      reclng(3) = box
      do 133 i = 1,3
       do 132 j = 1,3
        recip(i,j) = 0.d0
132    continue
       recip(i,i) = 1.d0/box
133   continue
      eigmin = 1.d0
c init our stuff
      write(6,*)'enter nfft and order'
      read(5,*)nfft,order
      write(6,*)'YOU ENTERED ',nfft,order
      nfft1 = nfft
      nfft2 = nfft
      nfft3 = nfft
      call pmesh_kspace_get_sizes(
     $     nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack)
      call pmesh_kspace_setup(
     $    bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork,
     $    nfft1,nfft2,nfft3,order,sizfftab,sizffwrk)
      if ( siz_Q .gt. MAXT )then
       write(6,*)'fft needs more room'
       stop
      endif
c every step do these
      write(6,*)'done setting up'
      call do_pmesh_kspace(
     $   numatoms,x,y,z,cg,recip,volume,ewald_coeff,
     $   order,nfft1,nfft2,nfft3,
     $   eer,fx,fy,fz,virial,
     $   sizfftab,sizffwrk,siztheta,siz_Q,
     $   bsp_mod1,bsp_mod2,bsp_mod3,fftable,
     $   Q,ffwork,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,
     $   fr1,fr2,fr3,time)
      write(6,*)'exact rec ene = ',rec_ene
      write(6,*)'pme recip ene = ',eer
      err = abs(rec_ene-eer)/abs(rec_ene)
      write(6,*)'relative error = ',err
      write(6,*)'log relative error = ',dlog(err)
    
      write(6,*)'checking exact total energy versus virial'
      call check_virial(self_ene,adj_ene,dir_ene,rec_ene,
     $       adj_vir,rec_vir,dir_vir)
      write(6,*)'checking approx total energy versus virial'
      call check_virial(self_ene,adj_ene,dir_ene,eer,
     $       adj_vir,virial,dir_vir)
      write(6,*)'checking the approx recip force versus exact recip',
     $   ' in terms of total force'
      call check_force(numatoms,fx,fy,fz,rfx,rfy,rfz,fx1,fy1,fz1)
      do 235 i = 1,numatoms
       fx2(i) = appdfx(i)+afx(i)+fx(i)
       fy2(i) = appdfy(i)+afy(i)+fy(i)
       fz2(i) = appdfz(i)+afz(i)+fz(i)
235   continue
      write(6,*)'checking the approx total force versus exact total',
     $   ' in terms of total force'
      call check_force(numatoms,fx1,fy1,fz1,fx2,fy2,fz2,fx1,fy1,fz1)
      write(6,*)'enter number of repeats for timing '
      read(5,*)numreps
      write(6,*)'YOU ENTERED ',numreps
      do 290 i = 1,5
        time(i) = 0.d0
290   continue
      call second(tim1)
      do 300 i = 1,numreps
       call do_pmesh_kspace(
     $   numatoms,x,y,z,cg,recip,volume,ewald_coeff,
     $   order,nfft1,nfft2,nfft3,
     $   eer,fx2,fy2,fz2,virial,
     $   sizfftab,sizffwrk,siztheta,siz_Q,
     $   bsp_mod1,bsp_mod2,bsp_mod3,fftable,
     $   Q,ffwork,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,
     $   fr1,fr2,fr3,time)
300   continue
      call second(tim2)
      write(6,*)'time/per repeat = ',(tim2-tim1)/numreps
      write(6,*)'percentage times: '
      write(6,501)100*time(1)/(tim2-tim1)
      write(6,502)100*time(2)/(tim2-tim1)
      write(6,503)100*time(3)/(tim2-tim1)
      write(6,504)100*time(4)/(tim2-tim1)
      write(6,505)100*time(5)/(tim2-tim1)
501   format(1x,'bspline    fraction: ',f6.1)
502   format(1x,'fillgrid   fraction: ',f6.1)
503   format(1x,'FFT        fraction: ',f6.1)
504   format(1x,'scalar sum fraction: ',f6.1)
505   format(1x,'grad   sum fraction: ',f6.1)
 
      stop
      end
c-------------------------------------------------------
      subroutine fix_dim(nfft,nfftdim)
      integer nfft,nfftdim,n
      nfftdim = nfft
      n = nfft/2
      if ( nfft .eq. 2*n )nfftdim = nfft+1
      return
      end
