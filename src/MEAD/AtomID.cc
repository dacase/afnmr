/*  class AtomID is a key type for an associative array of Atoms.
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

$Id: AtomID.cc,v 1.12 2004/12/06 17:57:56 bashford Exp $
*/
#include "MEAD/AtomID.h"

AtomID::AtomID()
: resnum(0), atname(), chainid()
{
}

AtomID::AtomID(int rn, const string& an, const string& cid)
: resnum(rn), atname(an), chainid(cid)
{
}


ostream& AtomID::print(ostream &ost)
{
  if (chainid.empty())
    ost << resnum << ":" << atname;
  else
    ost << chainid << ":" << resnum << ":" << atname;
  return ost;
}

// AtomID.cc ends here
