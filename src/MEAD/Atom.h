/* An atom with charge and radius -*- C++ -*-
    Copyright (c) 1993--1995 by Donald Bashford and Tony You.

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

$Id: Atom.h,v 2.14 2004/12/06 17:57:55 bashford Exp $
*/
#ifndef _Atom_h
#define _Atom_h 1

#include <string>
using std::string;
#include "MEAD/Coord.h"

//!wrap!
struct Atom {
  Coord coord;
  float rad;
  float charge;
  string atname;
  string resname;
  string chainid;
  int resnum;
  Atom ();
  Atom (const Atom &);
  ~Atom() {}
  //!nowrap!
  Atom& operator = (const Atom& a);
  bool operator== (const Atom& b) const;
  bool operator!= (const Atom& a) const;
  //!nowrap!+
  void print() const;
  void print(ostream& os) const;
  //!nowrap!-
};

inline bool Atom::operator== (const Atom& b) const
{
  return atname == b.atname && resname == b.resname
    && resnum == b.resnum && chainid == b.chainid;
}

inline bool Atom::operator!= (const Atom& a) const
{
  return !(*this == a);
}

ostream& operator <<(ostream& os, const Atom& a);

extern Atom nil_atom;

// Added by Bergsma 5/14/01
// To provide the operators required by Python __cmp__
inline bool operator<(const Atom& a, const Atom& b)
{
  if (a.chainid == b.chainid) {
    if (a.resnum == b.resnum)
      if (a.resname == b.resname)
        return a.atname < b.atname;
      else
        return a.resname < b.resname;
    else
        return a.resnum < b.resnum;
    }
  else
    return a.chainid < b.chainid;
}

inline bool operator>(const Atom& a, const Atom& b)
{
  if (a.chainid == b.chainid) {
    if (a.resnum == b.resnum)
      if (a.resname == b.resname)
        return a.atname > b.atname;
      else
        return a.resname > b.resname;
    else
      return a.resnum > b.resnum;
    }
  else
    return a.chainid > b.chainid;
}

#if SWIGPP_LITERAL_INCLUDE

%addmethods Atom {
   void write() { self->print(cout); }
   bool operator<(const Atom& a, const Atom& b);
   bool operator>(const Atom& a, const Atom& b);
};

#endif // SWIGPP_LITERAL_INCLUDE

#endif

// Atom.h ends here
