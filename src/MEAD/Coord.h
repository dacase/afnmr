// This is -*- C++ -*-
#ifndef _Coord_h
#define _Coord_h 1

/* A class for 3-space coordinates including vector operations
    Copyright (c) 1993--1995 by Donald Bashford.

    This source code file is part of the MEAD (Macroscopic
    Electrostatics with Atomic Detail) package of objects and
    programs.  MEAD is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 1, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; see the file COPYING.  If not, write to
    the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
    02139, USA.

    Donald Bashford can be contacted by electronic mail by the address,
    bashford@scripps.edu, or by paper mail at Department of Molecular
    Biology, The Scripps Research Institute, 10666 North Torrey Pines
    Road, La Jolla, California 92037.

$Id: Coord.h,v 1.29 2013/10/30 21:30:41 bashford Exp $
*/


#include <iostream>
using std::ostream;
using std::istream;

//!wrap!
class Coord {
public:
  Coord ();
  Coord (float xi, float yi, float zi);
  Coord (const Coord& c);
  ~Coord () {}
  //!nowrap!
  Coord& operator= (const Coord& c);
  Coord& operator+= (const Coord& a);
  Coord& operator-= (const Coord& a);
  Coord& operator*= (float a);
  Coord& operator/= (float a);
  Coord operator- () const;
  int operator> (const Coord& a) const;
  int operator< (const Coord& a) const;

  // Dave added these, but do they really make sense here? FIXME?
  //!nowrap!+
  double sin_th();
  double cos_th();
  double sin_phi();
  double cos_phi();
  //!nowrap!-

  //!nowrap!+
  ostream& print(ostream& ost) {
    return ost << "(" << x << ", " << y << ", " << z << ")";
  }
  istream& read(istream& ist);
  //!nowrap!-

  float x, y, z;
};

#if SWIGPP_LITERAL_INCLUDE

// Wrap these as member functions for Python
%addmethods Coord {
  void write() { self->print(cout); }
  Coord cross (const Coord& b) {return cross(*self, b);}
  double dot (const Coord& b) {return dot(*self, b);}
  Coord operator+ (const Coord& a, const Coord& b);
  Coord operator- (const Coord& a, const Coord& b);
  Coord operator* (const Coord& a, float b);
  Coord operator* (float b, const Coord& a);
  float operator* (const Coord& a, const Coord& b);
  Coord operator/ (const Coord& a, float b);
  bool operator == (const Coord& a, const Coord& b);
  bool operator != (const Coord& a, const Coord& b);
};

#endif // SWIGPP_LITERAL_INCLUDE

inline Coord::Coord ()
{x = y = z = 0;}


inline Coord::Coord (float xi, float yi, float zi)
{x=xi; y=yi; z=zi;}

inline Coord::Coord (const Coord& c)
{  x = c.x;  y = c.y;  z = c.z;}

inline Coord& Coord::operator= (const Coord& c)
{  x = c.x;  y = c.y;  z = c.z; return *this;}

inline Coord& Coord::operator+= (const Coord& a)
{x+=a.x; y+=a.y; z+=a.z; return *this;}

inline Coord& Coord::operator-= (const Coord& a)
{x-=a.x; y-=a.y; z-=a.z; return *this;}

inline Coord& Coord::operator*= (float a)
{x*=a; y*=a; z*=a; return *this;}

inline Coord& Coord::operator/= (float a)
{x/=a; y/=a; z/=a; return *this;}

inline Coord Coord::operator- () const
{Coord a; a.x = -x; a.y = -y; a.z = -z; return a;}

inline int Coord::operator> (const Coord& a) const
{return x>a.x && y>a.y && z>a.z;}

inline int Coord::operator< (const Coord& a) const
{return x<a.x && y<a.y && z<a.z;}

inline istream& Coord::read(istream& ist)
{
  // Read input of form, "f f f" where "f" is a float.
  // The expression may be surrounded by parentheses
  // and/or the "f" may be separated by commas.

  float rx, ry, rz; // hold coord values here till read is known OK.
  rx = ry = rz = 0;
  char c = 0;
  int expect_paren = 0;
  ist >> c;
  if (ist.eof()) return ist;
  if (c == '(') expect_paren = 1;
  else ist.putback(c);
  ist >> rx >> c;
  if (c != ',') ist.putback(c);
  ist >> ry >> c;
  if (c != ',') ist.putback(c);
  ist >> rz;
  if (expect_paren) {
    ist >> c;
    if (c != ')') {
      ist.setstate(std::ios::failbit);
      return ist;
    }
  }
  x = rx;
  y = ry;
  z = rz;
  return ist;
}

inline Coord operator+ (const Coord& a, const Coord& b)
{Coord sum; sum.x=a.x+b.x; sum.y=a.y+b.y; sum.z=a.z+b.z; return sum;}

inline Coord operator- (const Coord& a, const Coord& b)
{Coord sum; sum.x=a.x-b.x; sum.y=a.y-b.y; sum.z=a.z-b.z; return sum;}

inline Coord operator* (const Coord& a, float b)
{Coord prod = a; prod*=b; return prod;}

inline Coord operator* (float b, const Coord& a)
{Coord prod = a; prod*=b; return prod;}

inline float operator* (const Coord& a, const Coord& b)
{return a.x*b.x + a.y*b.y + a.z*b.z;}

inline Coord cross (const Coord& a, const Coord& b)
{Coord prod(a.y*b.z - a.z*b.y, b.x*a.z - a.x*b.z, a.x*b.y - b.x*a.y); return prod;}

// Bergsma added this 5/2/01
inline double dot (const Coord& a, const Coord& b)
{return (double) (a * b);}

inline Coord operator/ (const Coord& a, float b)
{Coord prod = a; prod/=b; return prod;}

inline bool operator == (const Coord& a, const Coord& b)
{return b.x == a.x && b.y == a.y && b.z == a.z;}

inline bool operator != (const Coord& a, const Coord& b)
{return b.x != a.x || b.y != a.y || b.z != a.z;}

inline ostream& operator<<(ostream& s, Coord c)
{return c.print(s);}

inline istream& operator>>(istream& s, Coord& c)
{return c.read(s);}

// On msvc the define below pulls in macro defs line M_PI
#define _USE_MATH_DEFINES 1
#include <math.h>

// Dave added these, but do they really make sense here? FIXME?
inline double Coord::cos_th ()
{
  double norm = sqrt(x*x + y*y + z*z);
  if (norm == 0)
    return 1.0;
  else
    return z/norm;
}
inline double Coord::sin_th ()
{
const float epsilon = 1.e-7F;
double cos = cos_th();
if (fabs(fabs(cos) - 1) < epsilon)
  return 0;
else
  return sqrt (1 - cos*cos);
}

inline double Coord::cos_phi ()
{
  double norm = sqrt(x*x + y*y);
  if (norm == 0)
    return 1.0;
  else
    return x/norm;
}
inline double Coord::sin_phi ()
{
  double norm = sqrt(x*x + y*y);
  if (norm == 0)
    return 0.0;
  else
    return y/norm;
}

#endif
