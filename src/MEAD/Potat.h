/* Potentials at each atom with storage in files. -*- C++ -*-

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

$Id: Potat.h,v 2.12 2004/12/06 17:58:34 bashford Exp $
*/
#ifndef _Potat_h
#define _Potat_h 1


/*

The Potat class family associates a set of atoms with potentials centered
on their nuclii.  Its purpose is to save the "important" part of
an electrostatic calculation in a farily small data space so that
it can be used later and to read and write this info to and from binary
files.  Multimead uses it.  All Potats have an AtomSet
which dictates where the potential value will be defined, so to be
useful, Potats should be constructed with AtomSets.  OutPotat is
for when you know the potential and want to write it out.  InPotat
is for when you want to try to read the potential from a file.
The read and write functions of these subclasses return 1 on success
or zero on failure.  The multiplication by an AtomChargeSet is defined
to give the interaction energy between the potential and the charges.
Multiplication by more general ChargeDists is not possible.

TO DO: At the moment the binary file format is the old style which
wastes a lot of space writing all atomic coordinates and so on.
Implement a format with just a hash identifier, size, and the
potentials.

*/

#include "MEAD/AtomSet.h"
#include <string>
using std::string;
#include <map>
using std::map;

class AtomSet;
class AtomChargeSet;
typedef float Potval;

class Potat {
public:
  Potat();
  Potat(const AtomSet&);
  void zero();
  float operator*(const AtomChargeSet&);
  Potval operator[] (const AtomID&) const;
protected:
  int defined;
  AtomSet _atset;
  map<AtomID, Potval> _map;
};

class InPotat : public Potat {
public:
  InPotat(const AtomSet&);
  int read(const string&);
  int read_oldstyle(const string&);
};

class OutPotat : public Potat  {
public:
  OutPotat(const AtomSet&, const class ElstatPot&);
  int write(const string&) const ;
  int write_oldstyle(const string&) const ;
};

inline Potval Potat::operator[] (const AtomID& k) const
{
  std::map<AtomID,Potval>::const_iterator i = _map.find(k);
  if (i == _map.end())
    ::error("Potat::operator[]: key not found");
  return i->second;
}

#endif

// Potat.h ends here
