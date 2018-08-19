c Here is how splined erfc is used inside nonbond routine at this time
c the erftbdns is the number of grid points per unit interval
c currently set to be 10000
C NOTE the commented out explicit call forms of erfc and derivative ,used for 
c checking purposes, as well as simple coulomb for timing comparison
c Currently on sgi workstations, use of this spline is about 10-20%
c overhead wrt to simple coulomb, whereas explicit call is ~70% overhead

        del = 1.d0 / erftbdns
        delx = crds(1,n) - xk
        dely = crds(2,n) - yk
        delz = crds(3,n) - zk
        delr2 = delx*delx + dely*dely+delz*delz
        delr = sqrt(delr2)
        delr2inv = 1.d0 / delr2
        x = ewaldcof*delr
c cubic spline on erfc,derf
        ind = erftbdns*x + 1
        dx = x - (ind-1)*del
        erfcc = erf_arr(1,ind)+dx*(erf_arr(2,ind)+
     $          dx*(erf_arr(3,ind)+dx*erf_arr(4,ind)/3.d0)/2.d0)
        derfc = -(erf_arr(2,ind)+dx*(erf_arr(3,ind)+
     $          dx*erf_arr(4,ind)/3.d0))

c       call erfcfun(x,erfcc)
c       derfc = (2.d0/sqrt(pi))*exp(-x**2)

        eelt = eelt + cgi*charge(j)*erfcc/delr
        dfee = cgi*charge(j)*(ewaldcof*derfc
     $     + erfcc/delr)*delr2inv

c       eelt = eelt + cgi*charge(j)/delr
c       dfee = (cgi*charge(j)/delr)*delr2inv
