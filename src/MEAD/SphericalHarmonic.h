// -*- C++ -*-
#ifndef SphericalHarmonic_h
#define SphericalHarmonic_h

#include <vector>
#include <iostream>
using std::ostream;
#include <string>
using std::string;

// On msvc the define below pulls in macro defs line M_PI
#define _USE_MATH_DEFINES 1
#include <math.h>

typedef std::vector<double> vector_double;
class Legendre;
typedef std::vector<Legendre> vector_Legendre;

//!wrap!
class Poly {
public:
  Poly();
  Poly(double zero_coef);
  Poly(double zero_coef, double one_coef);
  Poly(const vector_double& coefs);
  Poly(const Poly&);
  //!nowrap!
  Poly& operator=(const Poly&);
  double operator()(double x) const;
  Poly derivative() const;
  int degree() const;
  int size() const;
  const vector_double& coefficients() const {return coefs;}
  void set_varstring (string s) {_varstring = s;}
  string varstring () const {return _varstring;}
  Poly operator-() const;
  Poly operator+(const Poly& o) const;
  Poly operator-(const Poly& o) const;
  Poly operator*(const Poly& o) const;
  Poly operator*(double f) const;
  Poly operator/(double d) const;
  bool operator==(const Poly& o) const;
  ostream& output(ostream&) const;
private:
  string _varstring;
  vector_double coefs;
};

inline Poly operator*(double f, const Poly& p)
{return p*f;}

inline ostream& operator<<(ostream& ost, const Poly& p)
{return p.output(ost);}

//!wrap!
class Legendre : public Poly {
public:
  Legendre();  // Generates P_O by default
  Legendre(int ell);
  Legendre(const Legendre&);
  //!nowrap!+
  Legendre& operator=(const Legendre&);
  // quasi-constructors
  //!nowrap!
  friend Legendre order_one_Legendre();
  friend Legendre next_Legendre_using_prev2(const Legendre& prev,
					     const Legendre& preprev);
  friend void make_legendre_series(unsigned lmax, vector_Legendre* vlegp);
  //!nowrap!-
  int ell() const;
private:
  Poly& operator=(const Poly&); // Don't allow putting in arbitrary Polys
  Legendre(const Poly& p)        // or constructing with them (except by friends)
    : Poly(p) {}
};


Legendre order_one_Legendre();
Legendre next_Legendre_using_prev2(const Legendre& prev,
				   const Legendre& preprev);
void make_legendre_series(unsigned lmax, vector_Legendre*);

// ------- associated Legendre functions -----------------

// declaration
class AssocLegendre {
public:
  AssocLegendre();
  AssocLegendre(const AssocLegendre& o);
  AssocLegendre& operator=(const AssocLegendre& o);
  AssocLegendre(int ell, int m);
  AssocLegendre(const Legendre& leg, int m);
  AssocLegendre(const Legendre& leg); // makes an m=0
  AssocLegendre make_next() const; // That is, next m (not next ell)
  // Like above, but overwrite *next.  Return true if success.
  bool make_next(AssocLegendre* next) const;
  int ell() const {return _ell;}
  int m() const {return _m;}
  double operator()(double x) const;
  ostream & output(ostream &ost) const;
protected:
  Poly _polypart;
  int _ell, _m;
};
    
ostream& operator<<(ostream& ost, const AssocLegendre& aleg);


// some inline implementation

inline AssocLegendre::AssocLegendre()
  : _polypart(Legendre()), _m(0), _ell(0)
{}

inline AssocLegendre::AssocLegendre(const AssocLegendre& o)
  : _polypart(o._polypart), _m(o._m), _ell(o._ell)
{}

inline AssocLegendre& AssocLegendre::operator=(const AssocLegendre& o)
{
  if (&o == this) return *this;
  _polypart = o._polypart;
  _ell = o._ell;
  _m = o._m;
  return *this;
}

inline AssocLegendre::AssocLegendre(const Legendre& leg)
  : _polypart(leg), _m(0), _ell(leg.ell())
{}

inline bool AssocLegendre::make_next(AssocLegendre* next) const
  // Like no-arg version but overwrite *next.  Return true if success.
  // No exception thowing, no local variables, so should be fast inlined.
{
  if (next == this) return false; // don't allow self overwriting!
  if (_m == _ell) {
    return false; 
  }
  next->_polypart = _polypart.derivative();
  next->_ell = _ell;
  next->_m = _m + 1;
  return true;
}


// ===============================

// the pre-spherical harmonic class (essentially, normalized AssocLegendres

class PreSpheHarm : public AssocLegendre {
public:
  PreSpheHarm();
  PreSpheHarm(const PreSpheHarm& o);
  PreSpheHarm& operator=(const PreSpheHarm& o);
  PreSpheHarm(int ell, int m);
  PreSpheHarm(const AssocLegendre& asl);
  // builds the ell = leg.ell(), m=0 obj.
  PreSpheHarm(const Legendre& leg);
  bool make_next(PreSpheHarm* next) const;
  friend void make_PreSpheHarm_series(int ellmax,
			      std::vector< std::vector< PreSpheHarm > >&);
};

void
make_PreSpheHarm_series(unsigned ellmax,
			std::vector< std::vector< PreSpheHarm > >& vvs);

// some inline implementations for PreSpheHarm

inline PreSpheHarm::PreSpheHarm()
: AssocLegendre()
{
  _polypart = _polypart * (1.0/sqrt(4.0*M_PI));
}

inline PreSpheHarm::PreSpheHarm(const PreSpheHarm& o)
: AssocLegendre(o)
{}

inline PreSpheHarm& PreSpheHarm::operator=(const PreSpheHarm& o)
{
  AssocLegendre::operator=(o);
  return *this;
}

inline PreSpheHarm::PreSpheHarm(const Legendre& leg)
  // builds the ell = leg.ell(), m=0 obj.;
  : AssocLegendre(leg)
{
  _polypart = _polypart * sqrt((2*_ell + 1.)/4.0/M_PI);
}

inline bool PreSpheHarm::make_next(PreSpheHarm* next) const
{
  if ( ! AssocLegendre::make_next(next) ) return false;
  next->_polypart = next->_polypart / sqrt(double((_ell+_m+1)*(_ell-_m)));
  return true;
}



// =========================

// ------- At last!  the SphericalHarmonic class ------------

#include <complex>
using std::complex;
using std::conj;
using std::norm;
using std::real;
using std::imag;

class SphericalHarmonic : public PreSpheHarm { 
public:
  SphericalHarmonic();
  SphericalHarmonic(const SphericalHarmonic&);
  SphericalHarmonic(int ell, int m);
  SphericalHarmonic(const PreSpheHarm&);
  SphericalHarmonic(const AssocLegendre&);
  SphericalHarmonic(const Legendre&); // makes m=0
  
  SphericalHarmonic& operator=(const SphericalHarmonic&);
  int m() const;
  bool make_next(SphericalHarmonic * next) const;
  bool make_minus_m(SphericalHarmonic* shpt) const;
  complex<double> operator()(double theta, double phi) const;
  ostream & output(ostream &ost) const;
private:
  bool negm;
};


ostream& operator<<(ostream& ost, const SphericalHarmonic&);

// --- implementation of some SphericalHarmonic members as inlines ---

inline SphericalHarmonic::SphericalHarmonic()
  : PreSpheHarm(), negm(false) {_polypart.set_varstring("cos(theta)");}

inline SphericalHarmonic::SphericalHarmonic(const SphericalHarmonic& o)
  : PreSpheHarm(o), negm(o.negm) {_polypart.set_varstring("cos(theta)");}

inline SphericalHarmonic::SphericalHarmonic(int ell, int m)
  : PreSpheHarm(ell, m >= 0 ? m : -m), negm(false)
{
  if(m<0) {
    negm = true;
    if (m % 2 != 0) // if odd m, sign change
      _polypart = - _polypart;
  }
  _polypart.set_varstring("cos(theta)");
}

inline SphericalHarmonic::SphericalHarmonic(const PreSpheHarm& o)
  : PreSpheHarm(o), negm(false) {_polypart.set_varstring("cos(theta)");}

inline SphericalHarmonic::SphericalHarmonic(const AssocLegendre& aleg)
  : PreSpheHarm(aleg), negm(false) {_polypart.set_varstring("cos(theta)");}

inline SphericalHarmonic::SphericalHarmonic(const Legendre& leg)
  : PreSpheHarm(leg), negm(false) {_polypart.set_varstring("cos(theta)");}

inline SphericalHarmonic&
SphericalHarmonic::operator=(const SphericalHarmonic& o)
{
  PreSpheHarm::operator=(o);
  negm = o.negm;
  return *this;
}

inline int SphericalHarmonic::m() const
{
  int am = AssocLegendre::m();
  return negm ? -am : am;
}

inline bool SphericalHarmonic::make_next(SphericalHarmonic * next) const
{
  if (negm) return false; // this operation only applies to positive m
  return PreSpheHarm::make_next(next);
}

// But this one is different!
inline complex<double>
SphericalHarmonic::operator()(double theta, double phi) const
{
  complex<double> posm_result = PreSpheHarm::operator()(cos(theta))
    * complex<double>(cos(_m*phi), sin(_m*phi));
  if (negm) 
    return conj(posm_result);
  else
    return posm_result;
}

inline bool
SphericalHarmonic::make_minus_m(SphericalHarmonic* shpt) const
{
  if (shpt == this) return false;
  if (shpt->ell() == ell() && shpt->m() == -m()) //nothing to do
    return true;
  (*shpt) = (*this);
  shpt->negm = ! negm;
  if(m() % 2 != 0)
    shpt->_polypart = - _polypart;
  return true;
}


// ------ class for a series of spherical harmonics up to some ell ---

class SphericalHarmonic_series {
public:
  SphericalHarmonic_series(unsigned ellmax);
  size_t ellmax() const {return vvs.size();}
  const SphericalHarmonic& operator()(int ell, int m) const // getting only
  {// FIXME! no range checking!  g++ library's vector doesn't provide at()
    return vvs[ell][m+ell];
  }
private:
  SphericalHarmonic& refto(int ell, int m) // for setting
  {// FIXME! no range checking!  g++ library's vector doesn't provide at()
    return vvs[ell][m+ell];
  }
  std::vector< std::vector<SphericalHarmonic> > vvs;
};

  

//======================== SWIG STUFF ===========
#if SWIGPP_LITERAL_INCLUDE
%addmethods Poly {
  Poly operator*(double f, const Poly& p);
  void write() {cout << *self << endl;}
 };

Legendre order_one_Legendre();
Legendre next_Legendre_using_prev2(const Legendre& prev,
				   const Legendre& preprev);

%{
int create_legendre_series(int lmax, vector_Legendre* vL)
{ make_legendre_series(lmax, vL); return 1; }
 %}

int create_legendre_series(int lmax, vector_Legendre* vL);


#endif // SWIGPP_LITERAL_INCLUDE


#endif
