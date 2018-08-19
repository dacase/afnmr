// -*- C++ -*-
#ifndef MomentAnalysis_h
#define MomentAnalysis_h

#include "MEAD/SphericalHarmonic.h"
#include "MEAD/Coord.h"


  
//!wrap!
class Moments {
  // Complex moments, indexed by ell and m (like spherical harmonics)
public:
  typedef complex< double > momtype; // The type of the individual moments

  Moments(unsigned maxell); // Alloc and initialize to zero,
                           // moments up to ell=maxell

  unsigned ellmax() const {return unsigned(_momvec.size() - 1);}

  //!nowrap!
  momtype& operator()(int ell, int m) // for getting and setting
  {// FIXME! no range checking!  g++ library's vector doesn't provide at()
    return _momvec[ell][m+ell];
  }

  //!nowrap!
  const momtype& operator()(int ell, int m) const // safer access
  {// FIXME! no range checking!  g++ library's vector doesn't provide at()
    return _momvec[ell][m+ell];
  }

  Moments operator-() const;
  Moments operator+(const Moments o) const;
  Moments operator-(const Moments o) const;
  Moments operator*(double f) const;
  //!nowrap!
  Moments operator*(momtype f) const;
  
private:
  std::vector< std::vector<momtype> > _momvec;
};  

inline
Moments operator*(double f, const Moments& m)
{return m*f;}
  

#include "MEAD/ElstatPot.h"
#include "MEAD/ChargeDist.h"

Moments
momentsOfElstatPot(ElstatPot esp, unsigned maxell,
		   float radius, Coord center=Coord(0,0,0));
  // Sample the values of esp on a sphere of given radius and center,
  // and return the spherical moments, up to ell=maxell, that would
  // produce that same potential.


Moments
momentsOfChargeDist(ChargeDist::const_iterator chbegin,
		    ChargeDist::const_iterator chend,
		    unsigned maxell,
		    Coord center=Coord(0,0,0));

inline Moments
momentsOfChargeDist(ChargeDist rho, unsigned maxell,
		    Coord center=Coord(0,0,0))
{
  return momentsOfChargeDist(rho.begin(), rho.end(), maxell, center);
}

inline Moments
momentsOfChargeDist(const ChargeDist_lett& rho, unsigned maxell,
		    Coord center=Coord(0,0,0))
{
  return momentsOfChargeDist(rho.pc_begin(), rho.pc_end(), maxell, center);
}

  // Return the spherical moments, up to ell=maxell, for rho.


void print_moments(Moments debmoms);
// print a little table of moments to standard out.

void compare_moments(const Moments& moma, const Moments& momb,
		     string aname, string bname);
  // See what is the ration of the moments at each ell and m
  // Displays to standard output.


//======================== SWIG STUFF ===========
#if SWIGPP_LITERAL_INCLUDE
%addmethods Moments {
  Moments operator*(double f, const Moments& m)
  double get_real(int ell, int m) const
    { return real((*self)(ell,m)); }
  double get_imag(int ell, int m) const
    { return imag((*self)(ell,m)); }
  void set(int ell, int m, double real, double imag)
    { (*self)(ell, m) = Moments::momtype(real, imag); }
};

Moments
momentsOfElstatPot(ElstatPot esp, unsigned maxell,
		   float radius, Coord center);



void print_moments(Moments debmoms);

void compare_moments(const Moments& moma, const Moments& momb,
		     string aname, string bname);

%name(momentsOfChargeDist)
Moments alt_momentsOfChargeDist(Numeric_Nby4_DOUBLE* numarr,
				unsigned maxell, Coord center);

%{

#include "Numeric/arrayobject.h"

typedef PyArrayObject Numeric_Nby4_DOUBLE;

Moments
alt_momentsOfChargeDist(Numeric_Nby4_DOUBLE* numarr, unsigned maxell,
			Coord center)
{
  ManyPointCharge * mpcp = new ManyPointCharge;
  for (int i=0; i < numarr->dimensions[0]; ++i) {
    double x = * (double*) (numarr->data
			    + i*numarr->strides[0]);
    double y = * (double*) (numarr->data
			    + i*numarr->strides[0]
			    + 1*numarr->strides[1]);
    double z = * (double*) (numarr->data
			    + i*numarr->strides[0]
			    + 2*numarr->strides[1]);
    double q = * (double*) (numarr->data
			    + i*numarr->strides[0]
			    + 3*numarr->strides[1]);
    mpcp->push_back(PointCharge(Coord(x,y,z), q));
  }
  ChargeDist rho(mpcp);
  return momentsOfChargeDist(rho, maxell, center);
}

%}

#endif // SWIGPP_LITERAL_INCLUDE

#endif // MomentAnalysis_h
