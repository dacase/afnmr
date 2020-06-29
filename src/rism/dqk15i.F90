      subroutine dqk15i(f,boun,inf,a,b,result,abserr,resabs,resasc)
!***begin prologue  dqk15i
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a3a2,h2a4a2
!***keywords  15-point transformed gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the original (infinite integration range is mapped
!            onto the interval (0,1) and (a,b) is a part of (0,1).
!            it is the purpose to compute
!            i = integral of transformed integrand over (a,b),
!            j = integral of abs(transformed integrand) over (a,b).
!***description
!
!           integration rule
!           standard fortran subroutine
!           double precision version
!
!           parameters
!            on entry
!              f      - double precision
!                       fuction subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              boun   - double precision
!                       finite bound of original integration
!                       range (set to zero if inf = +2)
!
!              inf    - integer
!                       if inf = -1, the original interval is
!                                   (-infinity,bound),
!                       if inf = +1, the original interval is
!                                   (bound,+infinity),
!                       if inf = +2, the original interval is
!                                   (-infinity,+infinity) and
!                       the integral is computed as the sum of two
!                       integrals, one over (-infinity,0) and one over
!                       (0,+infinity).
!
!              a      - double precision
!                       lower limit for integration over subrange
!                       of (0,1)
!
!              b      - double precision
!                       upper limit for integration over subrange
!                       of (0,1)
!
!            on return
!              result - double precision
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule(resk) obtained by optimal addition
!                       of abscissae to the 7-point gauss rule(resg).
!
!              abserr - double precision
!                       estimate of the modulus of the absolute error,
!                       which should equal or exceed abs(i-result)
!
!              resabs - double precision
!                       approximation to the integral j
!
!              resasc - double precision
!                       approximation to the integral of
!                       abs((transformed integrand)-i/(b-a)) over (a,b)
!
!***references  (none)
!***routines called  d1mach
!***end prologue  dqk15i
!
      double precision a,absc,absc1,absc2,abserr,b,boun,centr,dabs,dinf,&
     &  dmax1,dmin1,d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,&
     &  resabs,resasc,resg,resk,reskh,result,tabsc1,tabsc2,uflow,wg,wgk,&
     &  xgk
      integer inf,j
      external f
!
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(8)
!
!           the abscissae and weights are supplied for the interval
!           (-1,1).  because of symmetry only the positive abscissae and
!           their corresponding weights are given.
!
!           xgk    - abscissae of the 15-point kronrod rule
!                    xgk(2), xgk(4), ... abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point gauss rule
!
!           wgk    - weights of the 15-point kronrod rule
!
!           wg     - weights of the 7-point gauss rule, corresponding
!                    to the abscissae xgk(2), xgk(4), ...
!                    wg(1), wg(3), ... are set to zero.
!
      data wg(1) / 0.0d0 /
      data wg(2) / 0.129484966168869693270611432679082d0 /
      data wg(3) / 0.0d0 /
      data wg(4) / 0.279705391489276667901467771423780d0 /
      data wg(5) / 0.0d0 /
      data wg(6) / 0.381830050505118944950369775488975d0 /
      data wg(7) / 0.0d0 /
      data wg(8) / 0.417959183673469387755102040816327d0 /
!
      data xgk(1) / 0.991455371120812639206854697526329d0 /
      data xgk(2) / 0.949107912342758524526189684047851d0 /
      data xgk(3) / 0.864864423359769072789712788640926d0 /
      data xgk(4) / 0.741531185599394439863864773280788d0 /
      data xgk(5) / 0.586087235467691130294144838258730d0 /
      data xgk(6) / 0.405845151377397166906606412076961d0 /
      data xgk(7) / 0.207784955007898467600689403773245d0 /
      data xgk(8) / 0.000000000000000000000000000000000d0 /
!
      data wgk(1) / 0.022935322010529224963732008058970d0 /
      data wgk(2) / 0.063092092629978553290700663189204d0 /
      data wgk(3) / 0.104790010322250183839876322541518d0 /
      data wgk(4) / 0.140653259715525918745189590510238d0 /
      data wgk(5) / 0.169004726639267902826583426598550d0 /
      data wgk(6) / 0.190350578064785409913256402421014d0 /
      data wgk(7) / 0.204432940075298892414161999234649d0 /
      data wgk(8) / 0.209482141084727828012999174891714d0 /
!
!
!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc*  - abscissa
!           tabsc* - transformed abscissa
!           fval*  - function value
!           resg   - result of the 7-point gauss formula
!           resk   - result of the 15-point kronrod formula
!           reskh  - approximation to the mean value of the transformed
!                    integrand over (a,b), i.e. to i/(b-a)
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!
!***first executable statement  dqk15i
      epmach = d1mach(4)
      uflow = d1mach(1)
      dinf = min0(1,inf)
!
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      tabsc1 = boun+dinf*(0.1d+01-centr)/centr
      fval1 = f(tabsc1)
      if(inf.eq.2) fval1 = fval1+f(-tabsc1)
      fc = (fval1/centr)/centr
!
!           compute the 15-point kronrod approximation to
!           the integral, and estimate the error.
!
      resg = wg(8)*fc
      resk = wgk(8)*fc
      resabs = dabs(resk)
      do 10 j=1,7
        absc = hlgth*xgk(j)
        absc1 = centr-absc
        absc2 = centr+absc
        tabsc1 = boun+dinf*(0.1d+01-absc1)/absc1
        tabsc2 = boun+dinf*(0.1d+01-absc2)/absc2
        fval1 = f(tabsc1)
        fval2 = f(tabsc2)
        if(inf.eq.2) fval1 = fval1+f(-tabsc1)
        if(inf.eq.2) fval2 = fval2+f(-tabsc2)
        fval1 = (fval1/absc1)/absc1
        fval2 = (fval2/absc2)/absc2
        fv1(j) = fval1
        fv2(j) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(j)*fsum
        resabs = resabs+wgk(j)*(dabs(fval1)+dabs(fval2))
   10 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resasc = resasc*hlgth
      resabs = resabs*hlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.d0) abserr = resasc*&
     & dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1&
     & ((epmach*0.5d+02)*resabs,abserr)
      return
      end
