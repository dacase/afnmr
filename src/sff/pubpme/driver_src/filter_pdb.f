c  filter large pdb to smaller pdb
      implicit none
      double precision xx,yy,zz,cg,ww
      double precision xx1,xx2,yy1,yy2,zz1,zz2,box
      integer i,n,j,bin(100),cum_bin(100),numwats,numatoms

      open(unit=8,file='small.pdb',status='new')
      n = 0
      do 10 i = 1,100
        bin(i) = 0
        cum_bin(i) = 0
10    continue
      open(unit=7,file='big.pdb',status='old')
c pass one get center of mass
      do 25 i = 1,100000
       read(7,20,end=30)xx,yy,zz,cg
20     format(30x,3f8.3,f7.3)
       n = n + 1
       ww = dmax1(xx,yy,zz)
       j = int(ww) + 1
       bin(j) = bin(j) + 1
       read(7,20)xx,yy,zz,cg
       read(7,20)xx,yy,zz,cg
25    continue
30    continue
      do 200 i = 1,100
       do 150 j = 1,i
        cum_bin(i) = cum_bin(i) + bin(j)
150    continue
200   continue
      close(7)
      write(6,*)'num atoms = ',n
      do 300 i = 1,100
       write(6,290)i,bin(i),cum_bin(i)
300   continue
290   format(1x,'bin# and #waters in this bin,cum_bin = ',i3,2i6)
c pass 2 get clip region
      write(6,*)'how many waters in box?'
      read(5,*)numwats
      do 400 i = 1,100
       if ( cum_bin(i) .ge. numwats)then
        box = i
        numatoms = 3*cum_bin(i)
        goto 401
       endif
400   continue
401   continue
      write(6,*)'box = ',box,' numatoms = ',numatoms
      open(unit=7,file='big.pdb',status='old')
      write(8,402)box,numatoms
402   format(1x,'box size,numatoms = ',f9.3,i6)
      do 425 i = 1,100000
       read(7,20,end=430)xx,yy,zz,cg
       ww = dmax1(xx,yy,zz)
       read(7,20)xx1,yy1,zz1,cg
       read(7,20)xx2,yy2,zz2,cg
       if ( ww .lt. box )then
        write(8,420)xx,yy,zz
        write(8,420)xx1,yy1,zz1
        write(8,420)xx2,yy2,zz2
       endif
420    format(30x,3f8.3)
425   continue
430   continue
      
      stop
      end
c--------------------------------------------------------
