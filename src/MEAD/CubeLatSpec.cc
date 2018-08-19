/* class carrying the specification for constructing a cubic lattice.
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

$Id: CubeLatSpec.cc,v 2.6 2004/12/06 17:58:08 bashford Exp $
*/

#include "MEAD/globals.h"
#include "MEAD/CubeLatSpec.h"

// FIXME? Eliminate the need for this in DielCubeRepMap::operator[]?
CubeLatSpec::CubeLatSpec()
: ngrid(0), spacing(0), center(0.0, 0.0, 0.0)
{
}

CubeLatSpec::CubeLatSpec(int ngrd, float spc, Coord cntr)
:ngrid(ngrd), spacing(spc), center(cntr)
{
  style = SPECIFIC;
  resolved = 1;
}

CubeLatSpec::CubeLatSpec(int ngrd, float spc, CenteringStyle sty)
:ngrid(ngrd), spacing(spc), center(0.0, 0.0, 0.0), style(sty)
{
  resolved = 0;
  if (style == SPECIFIC) {
    error ("ERROR: CubeLatSpec cannot be explicitly constructed with",
	   "style = SPECIFIC");
  }
}


CubeLatSpec::CubeLatSpec(const CubeLatSpec& k)
:  ngrid(k.ngrid), spacing(k.spacing), center(k.center),
   style(k.style), resolved(k.resolved)
{}

void CubeLatSpec::resolve (Coord geom_cent, Coord center_of_intr)
{
  Coord origin(0.0, 0.0, 0.0);
  switch (style) {
  case ON_ORIGIN:
    center = origin;
    break;
  case ON_CENT_OF_INTR:
    center = center_of_intr;
    break;
  case ON_GEOM_CENT:
    center = geom_cent;
    break;
  case SPECIFIC:
// case SPECIFIC requires no action
    break;
  }
  resolved = 1;
}

ostream & CubeLatSpec::print(ostream& ost) const
{
  if (style == SPECIFIC)
    ost << center;
  else {
    switch (style) {
    case ON_ORIGIN:
      ost << "ON_ORIGIN";
      break;
    case ON_CENT_OF_INTR:
      ost << "ON_CENT_OF_INTR";
      break;
    case ON_GEOM_CENT:
      ost << "ON_GEOM_CENT";
      break;
    }
    if (resolved)
      ost << " resolved to " << center << ")";
  }
  ost << " " << ngrid << " " << spacing;
  return ost;
}

ostream & operator<< (ostream& ost, const CubeLatSpec& cls)
{return cls.print(ost);}

// CubeLatSpec.cc ends here
