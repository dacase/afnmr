      implicit none
      integer MAX,NAVAIL
      parameter(MAX=5000)
      parameter(NAVAIL=300*MAX)
      double precision x(MAX),y(MAX),z(MAX),cg(MAX)
      double precision box,cgh,cgo,factor
      double precision cutoff,cut2,ewaldcof,maxexp
      double precision dtol,rtol,tol
      double precision eigmin,reclng(3),recip(3,3),volume
      double precision force(3,MAX),frac(3,MAX),XX(NAVAIL)
      double precision dir_ene,adj_ene,rec_ene,self_ene,appdene
      double precision dir_vir(6),adj_vir(6),rec_vir(6)
      double precision dfx(MAX),dfy(MAX),dfz(MAX)
      double precision appdfx(MAX),appdfy(MAX),appdfz(MAX),appdvir(6)
      double precision afx(MAX),afy(MAX),afz(MAX)
      double precision rfx(MAX),rfy(MAX),rfz(MAX)
      double precision fx(MAX),fy(MAX),fz(MAX)
      integer numatoms,numwats,i,j,n,mlimit(3)


      open(unit=10,file='exact.bin',status='new',form='unformatted')
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
      do 50 i = 1,numatoms
        dfx(i) = 0.d0
        appdfx(i) = 0.d0
        afx(i) = 0.d0
        rfx(i) = 0.d0
        dfy(i) = 0.d0
        appdfy(i) = 0.d0
        afy(i) = 0.d0
        rfy(i) = 0.d0
        dfz(i) = 0.d0
        appdfz(i) = 0.d0
        afz(i) = 0.d0
        rfz(i) = 0.d0
50    continue
      write(6,*)'enter direct space cutoff '
      read(5,*)cutoff
      write(6,*)'YOU ENTERED ',cutoff
      write(6,*)'enter tolerance '
      read(5,*)dtol
      write(6,*)'YOU ENTERED ',dtol
      tol = dtol

c For this value of cutoff and direct space tolerance dtol, first find
c the ewald coefficient ewaldcof.

      call find_ewaldcof(cutoff,dtol,ewaldcof)
      write(6,*)'ewald coefficient = ',ewaldcof

c Next find cutoff in recip space for this value of ewald coefficient
c and reciprocal space tolerance rtol. This sets up "exact" reciprocal for
c the given value of ewaldcof.

      rtol = 1.d-18
      call find_maxexp(ewaldcof,rtol,maxexp)

c Next find direct space cutoff for "exact" direct sum with given value
c of ewaldcof. 

      dtol = 1.d-18
      call find_cutoff(cut2,dtol,ewaldcof)

      write(6,*)'direct space cutoff for "exact" = ',cut2
      write(6,*)'maxexp(reciprocal space cutoff) = ',maxexp

c check if you can get all this larger sphere

      if ( cut2 .gt. 1.5d0*box )then
        write(6,*)'program will fail to get all interactions ',
     $    'within the enlarged sphere '
        stop
      endif

c Next find the number of reciprocal space vectors needed.
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
      call get_mlim(maxexp,mlimit,eigmin,reclng,recip)
c use this to do recip
      write(6,*)'doing recip'
      call reg_kspace(numwats,x,y,z,cg,ewaldcof,
     $       rec_ene,force,frac,recip,
     $       maxexp,mlimit,volume,box,XX,NAVAIL,rec_vir)
      write(6,*)'rec_ene = ',rec_ene
      do 200 i = 1,numatoms
       rfx(i) = force(1,i)
       rfy(i) = force(2,i)
       rfz(i) = force(3,i)
200   continue

c get self and adjustment values

      call self(cg,numatoms,self_ene,ewaldcof)
      write(6,*)'self ene = ',self_ene
      call adjust(numwats,x,y,z,cg,ewaldcof, 
     $    adj_ene,afx,afy,afz,adj_vir)
      write(6,*)'adjust ene = ',adj_ene

c get direct sum
      write(6,*)'doing direct'
      call direct(numwats,x,y,z,cg,cut2,ewaldcof,box,
     $      dir_ene,dfx,dfy,dfz,dir_vir)
      write(6,*)'dir ene = ',dir_ene
      call check_virial(self_ene,adj_ene,dir_ene,rec_ene,
     $       adj_vir,rec_vir,dir_vir)
      call direct(numwats,x,y,z,cg,cutoff,ewaldcof,box,
     $      appdene,appdfx,appdfy,appdfz,appdvir)
      write(6,*)'approx dir ene = ',appdene
      do 300 i = 1,numatoms
       fx(i) = dfx(i) + afx(i) + rfx(i)
       fy(i) = dfy(i) + afy(i) + rfy(i)
       fz(i) = dfz(i) + afz(i) + rfz(i)
300   continue
      call check_force(numatoms,
     $      appdfx,appdfy,appdfz,dfx,dfy,dfz,fx,fy,fz)

      write(10)numatoms
      write(10)cutoff,tol,ewaldcof,box
      write(10)self_ene,dir_ene,adj_ene,rec_ene
      write(10)dir_vir,adj_vir,rec_vir
      write(10)(dfx(i),i=1,numatoms)
      write(10)(dfy(i),i=1,numatoms)
      write(10)(dfz(i),i=1,numatoms)
      write(10)(afx(i),i=1,numatoms)
      write(10)(afy(i),i=1,numatoms)
      write(10)(afz(i),i=1,numatoms)
      write(10)(rfx(i),i=1,numatoms)
      write(10)(rfy(i),i=1,numatoms)
      write(10)(rfz(i),i=1,numatoms)
      write(10)(appdfx(i),i=1,numatoms)
      write(10)(appdfy(i),i=1,numatoms)
      write(10)(appdfz(i),i=1,numatoms)
      close(10)
      stop
      end
c------------------------------------------------------------------
