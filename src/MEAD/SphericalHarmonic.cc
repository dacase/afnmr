#include "MEAD/SphericalHarmonic.h"
#include "MEAD/globals.h"

#include <iostream>

#include <math.h>

typedef std::vector<double>::size_type v_size_type;
typedef std::vector<double>::iterator itr;
typedef std::vector<double>::const_iterator citr;

Poly::Poly()
  : coefs(), _varstring("x")
{
}

Poly::Poly(double zero_coef)
  : coefs(v_size_type(1), zero_coef), _varstring("x")
{
}

Poly::Poly(double zero_coef, double one_coef)
  : coefs(v_size_type(2), zero_coef), _varstring("x")
{
  coefs[1] = one_coef;
}

Poly::Poly(const std::vector<double> & cs)
  : coefs(cs), _varstring("x")
{
}

Poly::Poly(const Poly& o)
  :coefs(o.coefs), _varstring(o._varstring)
{
}

Poly& Poly::operator=(const Poly& other)
{
  if (&other != this) {
    coefs = other.coefs;
    _varstring = other._varstring;
  }
  return *this;
}

double Poly::operator()(double x) const
  //evaluation
{
  size_t sz = coefs.size();
  if (sz == 0) return 0.0;
  double sum = coefs[0];
  double xpow = 1.0;
  for(size_t i=1; i < sz; ++i) {
    xpow *= x;
    sum += coefs[i]*xpow;
  }
  return sum;
}

Poly Poly::derivative() const
{
  size_t sz = coefs.size();
  std::vector<double> newcoefs(v_size_type(sz-1), 0.0);
  for (size_t i=1; i<sz; ++i)
    newcoefs[i-1] = i*coefs[i];
  return Poly(newcoefs);
}

int Poly::degree() const
{
  return size() - 1;
}

int Poly::size() const
{
  return static_cast<int>(coefs.size());
}


Poly Poly::operator-() const // unary minus
{
  v_size_type sz = coefs.size();
  std::vector<double> newcoefs(sz, 0.0);
  itr j = newcoefs.begin();
  for (citr i = coefs.begin(); i != coefs.end(); ++i)
    *j++ = - *i;
  return Poly(newcoefs);
}

Poly Poly::operator+(const Poly& other) const
{
  size_t szthis = coefs.size();
  size_t szother = other.coefs.size();
  citr blarger;
  citr bsmaller;
  citr elarger;
  citr esmaller;
  if (szthis > szother) {
    blarger = coefs.begin();
    elarger = coefs.end();
    bsmaller = other.coefs.begin();
    esmaller = other.coefs.end();
  }
  else {
    blarger = other.coefs.begin();
    elarger = other.coefs.end();
    bsmaller = coefs.begin();
    esmaller = coefs.end();
  }    
  std::vector<double> newcoefs(blarger, elarger);
  itr j = newcoefs.begin();
  for (citr i = bsmaller; i != esmaller; ++i)
    *j++ += *i;
  return Poly(newcoefs);
}

Poly Poly::operator-(const Poly& other) const
{
  return ((*this) + (-other));
}

Poly Poly::operator*(const Poly& other) const
{
  v_size_type szthis = coefs.size(), szother = other.size();
  std::vector<double> newcoefs(szthis+szother-1, 0.0);
  for (v_size_type i=0; i<szthis; ++i) {
    double ct = coefs[i];
    for (v_size_type j=0; j<szother; ++j) {
      newcoefs[i+j] += ct * other.coefs[j];
    }
  }
  return Poly(newcoefs);
}

Poly Poly::operator*(double f) const
{
  v_size_type sz = coefs.size();
  std::vector<double> newcoefs(sz, 0.0);
  itr j = newcoefs.begin();
  for (citr i = coefs.begin(); i != coefs.end(); ++i)
    *j++ = (*i)*f;
  return Poly(newcoefs);
}

Poly Poly::operator/(double d) const
{
  if (d == 0.0)
    ::error("Poly division operator: Attempted division by zero");
  v_size_type sz = coefs.size();
  std::vector<double> newcoefs(sz, 0.0);
  itr j = newcoefs.begin();
  for (citr i = coefs.begin(); i != coefs.end(); ++i)
    *j++ = (*i)/d;
  return Poly(newcoefs);
}

bool Poly::operator==(const Poly& other) const
{
  return coefs == other.coefs;
}

ostream& Poly::output(ostream& ost) const
{
  if (coefs.size() == 0) return ost;
  if (coefs.size() == 1) {
    ost << coefs[0];
    return ost;
  }
  // degree 1 or more, so supress a leading zero
  //  if (coefs[0] != 0.0)
  //    ost << coefs[0];
  bool firstout = false;
  // sign indicators between terms
  const char* pos_signstring = " + ";
  const char* neg_signstring = " - ";
  // sign indicator for first printed term
  const char* init_pos_signstring = "";
  const char* init_neg_signstring = "-";

  for (v_size_type i=0; i<coefs.size(); ++i) {
    const double c = coefs[i];
    if (c==0.0) continue;
    const char* signstring;
    if (firstout) {
      signstring = c>0.0 ? pos_signstring : neg_signstring;
    }
    else {
      signstring = c>0.0 ? init_pos_signstring : init_neg_signstring;
    }
    ost << signstring;
    const double cab = fabs(c);
    if (i==0) {
      ost << cab;
    }
    else {
      if (cab != 1.0) ost << cab;
      if (i==1) ost << _varstring;
      else ost << _varstring << "**" << i;
    }
    firstout = true;
  }
  return ost;
}

Legendre::Legendre()
  : Poly(1.0)
{
}

Legendre::Legendre(int ell)
{
  if (ell < 0)
    ::error("constructor, Legendre(int ell): negative ell not allowed");
  if (ell == 0)
    operator=(Legendre());
  else if (ell == 1)
    operator=(order_one_Legendre());
  else {
    Legendre prev; // initially the ell=0
    Legendre cur(1);
    for (int i=2; i<=ell; ++i) {
      Legendre next = next_Legendre_using_prev2(cur, prev);
      prev = cur;
      cur = next;
    }
    operator=(cur);
    // operator=(next_Legendre_using_prev2(Legendre(ell-1), Legendre(ell-2)));
  }
}

Legendre::Legendre(const Legendre& other)
  : Poly(other)
{
}

Legendre& Legendre::operator=(const Legendre& other)
{
  Poly::operator=(other);
  return *this;
}


int Legendre::ell() const
{
  return degree();
}

Legendre order_one_Legendre()
{
  return Poly(0.0, 1.0);
}

Legendre
next_Legendre_using_prev2(const Legendre& prev, const Legendre& preprev)
{
  // Use the recursion relation, eq. 3.29a from Jackson 2nd ed.
  int ell = prev.ell();
  if (preprev.ell() != ell-1)
    ::error("next_Legendre_using_prev2:\n"
	    "   ell of preprev must be one less than that of prev.");
  Poly x = Poly(0.0, 1.0);
  Legendre result ( (double(2*ell + 1)/(ell+1)) * (x * prev)
		    - (double(ell)/(ell+1))*preprev );
  return result;
}

void make_legendre_series(unsigned lmax, std::vector<Legendre>* vlegp)
{
  if (lmax < 0)
    ::error("make_legendre_series: lmax must be non-negative");
  std::vector<Legendre> & vleg = *vlegp;
  vleg.clear();
  if (vleg.capacity() > lmax+1) vleg.reserve(lmax+1);
  vleg.push_back(Legendre()); // put in the ell=zero term.
  if (lmax == 0) return;
  vleg.push_back(order_one_Legendre()); // put in ell=1 term.
  if (lmax == 1) return;
  Poly x = Poly(0.0, 1.0);
  for (unsigned ell=1; ell<lmax; ++ell) {
    const Legendre& prev = vleg[ell];
    const Legendre& preprev = vleg[ell-1];
    vleg.push_back( (double(2*ell + 1)/(ell+1)) * (x * prev)
		    - (double(ell)/(ell+1))*preprev );
  }
}

// ------- associated Legendre functions -----------------



AssocLegendre::AssocLegendre(int ell, int m)
  : _polypart(Legendre(ell)), _m(m), _ell(ell)
{
  if (m > ell)
    ::error("constructor, AssogLegendre: m > ell not allowed");
  if ( m < -ell)
    ::error("constructor, AssogLegendre: m < -ell not allowed");
  if ( m < 0)
    ::error("constructor, AssogLegendre: negative m not supported!");
  for( ; m>0; --m) _polypart = _polypart.derivative();
}

AssocLegendre::AssocLegendre(const Legendre& leg, int m)
  : _polypart(leg), _m(m), _ell(leg.ell())
{
  if (m > _ell)
    ::error("constructor, AssogLegendre: m > ell not allowed");
  if ( m < -_ell)
    ::error("constructor, AssogLegendre: m < -ell not allowed");
  if ( m < 0)
    ::error("constructor, AssogLegendre: negative m not supported!");
  for( ; m>0; --m) _polypart = _polypart.derivative();
}

AssocLegendre AssocLegendre::make_next() const
  // That is, next m (not next ell)
{
  if (_m == _ell) {
    cerr << "AssocLegendre::Make_next called with m == ell" << endl;
    ::error("AssocLegendre::Make_next called with m == ell");
  }
  AssocLegendre next;
  next._polypart = _polypart.derivative();
  next._ell = _ell;
  next._m = _m + 1;
  return next;
}

double AssocLegendre::operator()(double x) const
{
  int sign = _m % 2 ? -1 : 1;
  return sign * pow((1.0 - x*x), double(_m/2.0)) * _polypart(x);
}

ostream & AssocLegendre::output(ostream &ost) const
{
  if (_m % 2 == 1)
    ost << "-(1 - " << _polypart.varstring() << "**2)**(" << _m << "/2) * ";
  else if (_m != 0)
    ost << "(1 - x**2)**" << _m/2 << " * ";
  if (_polypart.degree() > 0) 
    ost << "(" <<  _polypart << ")";
  else
    ost << _polypart;
  return ost;
}

    
ostream& operator<<(ostream& ost, const AssocLegendre& aleg)
{
  return aleg.output(ost);
}

// ===============================

// the pre-spherical harmonic class (essentially, normalized AssocLegendres

PreSpheHarm::PreSpheHarm(int ell, int m)
: AssocLegendre(ell, m)
{
  double sqrt_ell_minus_m_fac_over_ell_plus_m_fac = 1.0;
  for (int im=1; im <= m; ++im) {
    sqrt_ell_minus_m_fac_over_ell_plus_m_fac
      /= sqrt(double((ell-im+1)*(ell+im)));
  }
  _polypart = _polypart
    * sqrt_ell_minus_m_fac_over_ell_plus_m_fac;
  _polypart = _polypart * sqrt((2*ell + 1.)/4.0/M_PI);
}

PreSpheHarm::PreSpheHarm(const AssocLegendre& asl)
: AssocLegendre(asl)
{
  double sqrt_ell_minus_m_fac_over_ell_plus_m_fac = 1.0;
  for (int im=1; im <= _m; ++im) {
    sqrt_ell_minus_m_fac_over_ell_plus_m_fac
      /= sqrt(double((_ell-_m+1)*(_ell+_m)));
  }
  _polypart = _polypart
    * sqrt_ell_minus_m_fac_over_ell_plus_m_fac;
  _polypart = _polypart * sqrt((2*_ell + 1.)/4.0/M_PI);
}

void
make_PreSpheHarm_series(unsigned ellmax,
			std::vector< std::vector< PreSpheHarm > >& vvs)
{
  vvs.clear();
  if (vvs.capacity() < ellmax+1) vvs.reserve(ellmax+1);
  std::vector<Legendre> leg;
  leg.reserve(ellmax+1);
  make_legendre_series(ellmax, &leg);
  for (unsigned ell=0; ell<=ellmax; ++ell) {
    vvs.push_back(std::vector<PreSpheHarm>(ell+1));
    std::vector<PreSpheHarm> & cur_ell_ser = vvs.back();
    cur_ell_ser[0] = PreSpheHarm(leg[ell]);
    for (unsigned m=0; m < ell; ++m) {
      if ( ! cur_ell_ser[m].make_next(&cur_ell_ser[m+1]) )
	cout << "Trouble in make_PreSpheHarm_series" << endl;
    }
  }
}

// =========================

// ------- At last!  the SphericalHarmonic class ------------
// (most memfuns are inline in the .h file)

ostream & SphericalHarmonic::output(ostream &ost) const
{
  PreSpheHarm::output(ost);
  if (negm)
    ost << "*exp(-i*" << _m << "*phi)";
  else
    ost << "*exp(i*" << _m << "*phi)";
  return ost;
}

ostream& operator<<(ostream& ost, const SphericalHarmonic& s)
{return s.output(ost);}

// ----- series of them -------

SphericalHarmonic_series::SphericalHarmonic_series(unsigned maxell)
{
  std::vector<Legendre> leg;
  leg.reserve(maxell+1);
  make_legendre_series(maxell, &leg);
  vvs.reserve(maxell+1);
  for (unsigned ell = 0; ell <= maxell; ++ell) {
    vvs.push_back(std::vector<SphericalHarmonic>(2*ell+1));
    refto(ell,0) = SphericalHarmonic(leg[ell]); 
    for (unsigned m=0; m < ell; ++m) {
      const SphericalHarmonic& cur = (*this)(ell,m);
      if ( ! cur.make_next(&refto(ell,m+1)) )
	cout << "Trouble in SphericalHarmonic_series ctor" << endl;
    }
    for(int m=1; m<=static_cast<int>(ell); ++m) {
      if ( ! (*this)(ell,m).make_minus_m(&refto(ell,-m)) )
	cout << "Trouble making the minus m's" << endl;
    }
  }
}

