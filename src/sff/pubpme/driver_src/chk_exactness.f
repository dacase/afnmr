      implicit none
      integer MAX,NAVAIL
      parameter(MAX=25000)
      parameter(NAVAIL=75*MAX)
      double precision fx1(MAX),fy1(MAX),fz1(MAX)
      double precision fx2(MAX),fy2(MAX),fz2(MAX)
      double precision dir_ene,adj_ene,rec_ene,self_ene
      double precision dir_vir(6),adj_vir(6),rec_vir(6),svir
      double precision dfx(MAX),dfy(MAX),dfz(MAX)
      double precision afx(MAX),afy(MAX),afz(MAX)
      double precision rfx(MAX),rfy(MAX),rfz(MAX)
      integer numatoms,i
      double precision err_e
      double precision cutoff,tol,ewaldcof,box,etot1,etot2
      double precision vtot1(6),vtot2(6),err_vir


      open(unit=10,file='exact.bin',status='old',form='unformatted')
      open(unit=11,file='exact1.bin',status='old',form='unformatted')
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
      write(6,*)'numatoms=',numatoms
      write(6,*)'cutoff,tol,ewaldcof,box = ',
     $      cutoff,tol,ewaldcof,box
      etot1 = self_ene+dir_ene+adj_ene+rec_ene
      write(6,*)'etot1 = ',etot1
      do 10 i = 1,6
       vtot1(i) = dir_vir(i)+adj_vir(i)+rec_vir(i)
10    continue
      write(6,*)'virial = ',vtot1(1),vtot1(2),vtot1(3)
      write(6,*)'virial = ',vtot1(4),vtot1(5),vtot1(6)
      svir = vtot1(1) + vtot1(4) + vtot1(6)
      write(6,*)'scalar virial = ',svir
      err_vir = abs(svir+etot1)/abs(etot1)
      write(6,*)'error in virial identity = ',err_vir
      do 20 i = 1,numatoms
       fx1(i) = dfx(i)+afx(i)+rfx(i)
       fy1(i) = dfy(i)+afy(i)+rfy(i)
       fz1(i) = dfz(i)+afz(i)+rfz(i)
20    continue
      read(11)numatoms
      read(11)cutoff,ewaldcof,box
      read(11)self_ene,dir_ene,adj_ene,rec_ene
      read(11)dir_vir,adj_vir,rec_vir
      read(11)(dfx(i),i=1,numatoms)
      read(11)(dfy(i),i=1,numatoms)
      read(11)(dfz(i),i=1,numatoms)
      read(11)(afx(i),i=1,numatoms)
      read(11)(afy(i),i=1,numatoms)
      read(11)(afz(i),i=1,numatoms)
      read(11)(rfx(i),i=1,numatoms)
      read(11)(rfy(i),i=1,numatoms)
      read(11)(rfz(i),i=1,numatoms)
      write(6,*)'numatoms=',numatoms
      write(6,*)'cutoff,ewaldcof,box = ',
     $      cutoff,ewaldcof,box
      etot2 = self_ene+dir_ene+adj_ene+rec_ene
      write(6,*)'etot2 = ',etot2
      do 40 i = 1,6
       vtot2(i) = dir_vir(i)+adj_vir(i)+rec_vir(i)
40    continue
      do 50 i = 1,numatoms
       fx2(i) = dfx(i)+afx(i)+rfx(i)
       fy2(i) = dfy(i)+afy(i)+rfy(i)
       fz2(i) = dfz(i)+afz(i)+rfz(i)
50    continue
      write(6,*)'virial = ',vtot2(1),vtot2(2),vtot2(3)
      write(6,*)'virial = ',vtot2(4),vtot2(5),vtot2(6)
      svir = vtot2(1) + vtot2(4) + vtot2(6)
      write(6,*)'scalar virial = ',svir
      err_vir = abs(svir+etot2)/abs(etot2)
      write(6,*)'error in virial identity = ',err_vir
      err_e = 2.d0*abs(etot1-etot2)/(abs(etot1)+abs(etot2)) 
      write(6,*)'rel err ene = ',err_e
      call check_force(numatoms,fx1,fy1,fz1,fx2,fy2,fz2,fx1,fy1,fz1)
      end
