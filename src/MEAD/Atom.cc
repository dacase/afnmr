/* Implementation of a simple Atom class
   Copyright (C) 1990--1995 by Donald Bashford.

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

$Id: Atom.cc,v 1.10 2004/12/06 17:57:52 bashford Exp $
*/



#include "MEAD/Atom.h"
#include <iostream>

Atom nil_atom;

Atom::Atom ()
: coord(0,0,0),
  rad(0),
  charge(0),
  atname(""),
  resname(""),
  chainid(""),
  resnum(0)
{}

Atom::Atom (const Atom & a)
: resname(a.resname), atname(a.atname), chainid(a.chainid)
{
  resnum = a.resnum;
  coord = a.coord;
  rad = a.rad;
  charge = a.charge;
}

Atom& Atom::operator = (const Atom& a)
{
  chainid = a.chainid;
  resname = a.resname;
  atname = a.atname;
  resnum = a.resnum;
  coord = a.coord;
  rad = a.rad;
  charge = a.charge;
  return *this;
}



void Atom::print() const
{
  std::cout << chainid << " " << resname << " " << resnum << " " << atname
    << " rad = " << rad << " charge = " << charge << std::endl;
}

void Atom::print(ostream& os) const
{
  os << chainid << " " << resname << " " << resnum << " " << atname
    << " rad = " << rad << " charge = " << charge << "\n";
}

ostream& operator <<(ostream& os, const Atom& a)
{
  a.print(os);
  return os;
}

// Atom.cc ends here
