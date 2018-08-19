      implicit none
      integer MAX,NAVAIL
      parameter(MAX=55000)
      parameter(NAVAIL=200*MAX)
      double precision dir_ene,adj_ene,rec_ene,self_ene,vir(6)
      double precision dir_vir(6),adj_vir(6),rec_vir(6)
      double precision vtot1(6)
      double precision dfx(MAX),dfy(MAX),dfz(MAX)
      double precision afx(MAX),afy(MAX),afz(MAX)
      double precision rfx(MAX),rfy(MAX),rfz(MAX)
      double precision x(MAX),y(MAX),z(MAX),cg(MAX)
      double precision fx(MAX),fy(MAX),fz(MAX),ene
      double precision fx1(MAX),fy1(MAX),fz1(MAX)
      double precision fx2(MAX),fy2(MAX),fz2(MAX)
      double precision appdfx(MAX),appdfy(MAX),appdfz(MAX)
      double precision box,cgh,cgo,factor
      double precision cutoff,tol,ewaldcof,rtol,maxexp
      double precision eigmin,reclng(3),recip(3,3),volume
      double precision force(3,MAX),frac(3,MAX),XX(NAVAIL)
      integer numatoms,numwats,i,j,n,mlimit(3),numreps
      real tim1,tim2


      factor = sqrt(332.17752)
      cgh = 0.417d0*factor
      cgo = -2.d0*cgh
      open(unit=8,file='small.pdb',status='old')
      read(8,25)box,numatoms
25    format(21x,f9.3,i6)
      write(6,*)'box = ',box,' numatoms = ',numatoms
      numwats = numatoms/3
      if ( numatoms .gt. MAX )then
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
      read(10)cutoff,tol,ewaldcof,box
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
      do 10 i = 1,6
       vtot1(i) = dir_vir(i)+adj_vir(i)+rec_vir(i)
10    continue
      do 35 i = 1,numatoms
       fx1(i) = dfx(i)+afx(i)+rfx(i)
       fy1(i) = dfy(i)+afy(i)+rfy(i)
       fz1(i) = dfz(i)+afz(i)+rfz(i)
35    continue
      write(6,*)'numatoms = ',numatoms
      write(6,*)'cutoff,tol,ewaldcof,box = ',
     $     cutoff,tol,ewaldcof,box
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
      write(6,*)'enter recip tolerance'
      read(5,*)rtol
      write(6,*)'YOU ENTERED ',rtol
      call find_maxexp(ewaldcof,rtol,maxexp)
      write(6,*)'ewaldcof,maxexp = ',ewaldcof,maxexp
      call get_mlim(maxexp,mlimit,eigmin,reclng,recip)
      write(6,*)'doing recip'
      call reg_kspace(numwats,x,y,z,cg,ewaldcof,
     $       ene,force,frac,recip,
     $       maxexp,mlimit,volume,box,XX,NAVAIL,vir)
      do 150 i = 1,numatoms
       fx(i) = force(1,i)
       fy(i) = force(2,i)
       fz(i) = force(3,i)
150   continue
      write(6,*)'exact rec  ene = ',rec_ene
      write(6,*)'approx rec ene = ',ene
      write(6,*)'checking exact total energy versus virial'
      call check_virial(self_ene,adj_ene,dir_ene,rec_ene,
     $       adj_vir,rec_vir,dir_vir)
      write(6,*)'checking approx total energy versus virial'
      call check_virial(self_ene,adj_ene,dir_ene,ene,
     $       adj_vir,vir,dir_vir)
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
      call second(tim1)
      do 300 i = 1,numreps
       call reg_kspace(numwats,x,y,z,cg,ewaldcof,
     $       ene,force,frac,recip,
     $       maxexp,mlimit,volume,box,XX,NAVAIL,vir)
300   continue
      call second(tim2)
      write(6,*)'time/per repeat = ',(tim2-tim1)/numreps
 
      stop
      end
