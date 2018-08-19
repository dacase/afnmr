/* Specifies dimension spacing and centering of a cubic lattice -*- C++ -*-

    Copyright (c) 1993--1995 by Donald Bashford

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

$Id: CubeLatSpec.h,v 2.6 2004/12/06 17:58:09 bashford Exp $
*/
#ifndef _CubeLatSpec_h
#define _CubeLatSpec_h 1

#include "MEAD/Coord.h"
#include "MEAD/globals.h"
#include "MEAD/CenteringStyle.h"

//!wrap!
class CubeLatSpec {
public:
// FIXME? Eliminate need (see .cc file).
  //!nowrap!
  CubeLatSpec();
  CubeLatSpec(int ngrd, float spc, Coord cntr);
  //!nowrap!
  CubeLatSpec(int ngrd, float spc, CenteringStyle sty);
  CubeLatSpec(const CubeLatSpec& k);
  ~CubeLatSpec() {}
  //!nowrap!
  CubeLatSpec& operator = (const CubeLatSpec& k);
  int operator == (const CubeLatSpec& k);
  int operator < (const CubeLatSpec& k);
  int operator > (const CubeLatSpec& k);
  void resolve (Coord geom_cent, Coord center_of_intr);
  int get_grid_dim() const;
  float get_spacing() const;
  //!nowrap!+
  const Coord& get_center() const;
  ostream& print(ostream&) const;
  CenteringStyle get_centering_style() const { return style; }
  //!nowrap!-
  bool is_resolved() const { return resolved != 0; }
private:
  int ngrid;
  float spacing;
  Coord center;
  CenteringStyle style;
  short resolved;
};

#if SWIGPP_LITERAL_INCLUDE

%addmethods CubeLatSpec {
//
// Rewrap this constructor to take an int since Python doesn't know about enums
  CubeLatSpec(int ngrd, float spc, int sty);
//
// Rewrap this method to return an int
  int get_centering_style() { return (int) self->get_centering_style(); }
//
// Python doesn't want you to redefine the print method, so call it write
  void write() { self->print(cout); }
//
// Don't wrap methods that return references to wrapped classes because
// the __del__ operator will be called on them twice!
//
// Rewrap this method to return a Coord

  Coord get_center() { Coord c = self->get_center(); return c; }
};

%wrapper %{
//
// The wrapped constructor with CenteringStyle as an int
// Note: There is no range checking for the enum assignment, so the
//       result could be undefined if an out-of-range int is used.
//
CubeLatSpec * new_CubeLatSpec(int ngrd, float spc, int sty){
  CenteringStyle style = CenteringStyle(sty);
  return new CubeLatSpec(ngrd, spc, style);
}

%}

// This wraps up the CenteringStyle enum, although it creates an entirely
// new enum here. Take care to update this CenteringStyle_enum class as well as
// the python wrapper code (CenteringStyle.py) if the C++ CenteringStyle enum
// is changed.

%inline %{

class CenteringStyle_enum {
public:
  enum { ON_ORIGIN, ON_CENT_OF_INTR, ON_GEOM_CENT, SPECIFIC };
};

%}

// Python code to create a CenteringStyle instance
%pragma(python) include="CenteringStyle.py"

#endif // SWIGPP_LITERAL_INCLUDE

ostream & operator<< (ostream& ost, const CubeLatSpec& cls);

inline CubeLatSpec& CubeLatSpec::operator = (const CubeLatSpec& k)
{
  ngrid = k.ngrid; spacing = k.spacing; center = k.center;
  style = k.style;
  resolved = k.resolved;
  return *this;
}

inline int CubeLatSpec::operator == (const CubeLatSpec& k)
{
  if (!(resolved && k.resolved))
    error ("ERROR: CubeLatSpec::operator == :",
	   " comparison not possible on unresolved objects");
  return ngrid == k.ngrid && spacing == k.spacing && center == k.center;
}

// Added by bergsma 5/30/01
// To support the Python __cmp__ operation
inline int CubeLatSpec::operator < (const CubeLatSpec& k)
{
  if (!(resolved && k.resolved))
    error ("ERROR: CubeLatSpec::operator < :",
	   " comparison not possible on unresolved objects");
  if (ngrid == k.ngrid) {
    if (spacing == k.spacing)
      return center < k.center;
    else
      return spacing < k.spacing;
  }
  else
    return ngrid < k.ngrid;
}

// Added by bergsma 5/30/01
// To support the Python __cmp__ operation
inline int CubeLatSpec::operator > (const CubeLatSpec& k)
{
  if (!(resolved && k.resolved))
    error ("ERROR: CubeLatSpec::operator < :",
	   " comparison not possible on unresolved objects");
  if (ngrid == k.ngrid) {
    if (spacing == k.spacing)
      return center > k.center;
    else
      return spacing > k.spacing;
  }
  else
    return ngrid > k.ngrid;
}

inline int CubeLatSpec::get_grid_dim() const
{
  return ngrid;
}

inline float CubeLatSpec::get_spacing() const
{
  return spacing;
}

inline const Coord& CubeLatSpec::get_center() const
{
  if (!resolved)
    error ("ERROR: CubeLatSpec::get_center called for unresolved object");
  return center;
}


#endif

// CubeLatSpec.h ends here
