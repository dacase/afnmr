// -*- C++ -*-
#ifndef MomentAnalysis_tmplts_h
#define MomentAnalysis_tmplts_h

template<class RetType, class FunType>
RetType spherint (const FunType& funob,
		 unsigned ninter=1000)
     // integrate funob over the unit sphere.
     // FunType funob must be a function or callable object that
     // takes the polar coords (theta and phi), and returns RetType
     // FIXME! This is a really dumb integration algorithm
     //        (maybe use adaptation of Romberg method from libdevel?)
{
  int theta_inters = ninter;
  int phi_inters = ninter;
  double dtheta = M_PI/theta_inters;
  double dphi = 2.0*M_PI/phi_inters;
  // Trapezoid formula over theta, but using fact that integrand
  // fun*sin(theta) is formally zero at poles.
  RetType prenorm = 0.0;
  for (int itheta=1; itheta < theta_inters; ++itheta) {
    double theta = itheta*dtheta;
    // for the circular phi integration the usual factors of 1/2
    // at the ends in the trapezoid rule don't apply.
    RetType phisum = 0.0;
    for (int iphi=0; iphi < phi_inters; ++iphi) {
      double phi = iphi*dphi;
      //      cout << "theta = " << theta << ", phi = " << phi << ", f(theta,phi) = "
      //	   << funob(theta,phi) << endl;
      phisum += funob(theta, phi) * dphi;
    }
    prenorm += phisum*sin(theta)*dtheta;
  }
  return prenorm;
}

template<class FunType>
class ProdWithYstar {
public:
  ProdWithYstar(const SphericalHarmonic& y, const FunType& f) 
    : _y(y)
  {
    fp = &f;
  }
  complex<double> operator()(double theta, double phi) const
  {
    complex<double> yval = _y(theta, phi);
    return  conj(yval) * (*fp)(theta, phi);
  }
private:
  const FunType *fp;
  const SphericalHarmonic& _y;
};

template<class FunType> Moments
sphericalComponents(FunType fun, int ellmax)
  // fun must take spherical coords
  // theta and phi (i.e. it must be function on the unit sphere),
  // The coefficients for s spherical harmonic series up to ellmax
  // are returned in a Moments structure (though they aren't really moments)
{
  std::vector< std::vector< SphericalHarmonic > > vvs;
  SphericalHarmonic_series shs(ellmax);
  Moments moment(ellmax);
  for (int ell=0; ell <= ellmax; ++ell) {
    for (int m=-ell; m <= ell; ++m) {
      ProdWithYstar<FunType> integrand(shs(ell,m), fun);
      moment(ell,m) = spherint< Moments::momtype >(integrand, 100);
    }
  }
  return moment;
}
  
template<class FunType >
class USphereFunFromCartesian {
public:
  USphereFunFromCartesian(FunType cartfun, double r,
			  Coord center=Coord(0,0,0))
    : _cartfun(cartfun), _r(r), _cen(center) {}
  double operator()(double theta, double phi) const
    {
      double rsinth = _r*sin(theta);
      Coord c(float(rsinth*cos(phi)), float(rsinth*sin(phi)), float(_r*cos(theta)));
      return _cartfun(_cen + c);
    }
  FunType _cartfun;
  double _r;
  Coord _cen;
};


template<class FunType> Moments
momentsOfPotfun(FunType fun, int ellmax,
		 double radius, Coord center=Coord(0,0,0))
  // Moment analysis of function defined in cartesian 3-space
  // (taking Coord arg) by examining it's values on a sphere of
  // given radius about given center (origin by default).
{
  double fourpi = M_PI*4.01;
  USphereFunFromCartesian<FunType> usf(fun, radius, center);
  //  SphericalHarmonic Y21(2,1);
  Moments moments = sphericalComponents(usf, ellmax);
  for (int ell=0; ell <= ellmax; ++ell) {
    double rpow = pow(radius, ell+1);
    double fac = rpow*(2.0*ell + 1.0)/fourpi;
    for (int m = -ell; m <= ell; ++m) {
      moments(ell,m) *= fac;
    }
  }
  return moments;
}


#endif // MomentAnalysis_tmplts_h
