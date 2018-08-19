/* Specification of the sequence of grids to use for finite diff. proc.

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

$Id: FinDiffMethod.cc,v 2.9 2004/12/06 17:58:25 bashford Exp $
*/

#include "MEAD/globals.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/CubeLatSpec.h"

#include <string>

struct FDMLink {
  CubeLatSpec *fdm_ptr;
  FDMLink *finer;
};



FinDiffMethodRep::FinDiffMethodRep()
{
  finest = coarsest = ptr = 0;
  referenceCount = 1;
  resolved = 1;
}

// Read a file that tells the FinDiffMethod parameters for all levels.
// Each level is described by:
//    cen grid_dim spacing
// where cen can either be a Coord that specifies the center or one
// of the three CenteringStyles: ON_ORIGIN, ON_CENT_OF_INTR or ON_GEOM_CENT;
// grid_dim is the number of points along each edge of the cubic lattice;
// and spacing is the spacing between lattice points.

void FinDiffMethodRep::read (const string& filename_string)
{
  const char *filename = filename_string.c_str();
  ifstream gmin(filename);
  if (!gmin)
    ::error ("FinDiffMethodRep: cannot open file ", filename);
  string censtring;
  CenteringStyle censtl=ON_ORIGIN;
  int d;
  float sp;
  int levels_read = 0;
  while (gmin.good()) {
    gmin >> std::ws; // skip whitespace
    char c;
    gmin.get(c); // KCC needs this to trigger eof even when good is false!
    if (gmin.eof()) break;
    if (!gmin.good())
      ::error("INPUT FAILURE in FinDiffMethodRep::read from file",
	      (const char *) filename);
    gmin.putback(c);
    if (c=='O') {  // The letter, O, signals a CenteringStyle specifier
      gmin >> censtring >> d >> sp;
      if (!gmin.good())
	::error("INPUT FAILURE in FinDiffMethodRep::read from file",
		(const char *) filename, "bad CenteringStyle-type entry");
      if (censtring == "ON_ORIGIN")
	censtl = ON_ORIGIN;
      else if (censtring == "ON_CENT_OF_INTR")
	censtl = ON_CENT_OF_INTR;
      else if (censtring == "ON_GEOM_CENT")
	censtl = ON_GEOM_CENT;
      else // censtring gets a bad value so bail out..
	::error("INPUT FAILURE in FinDiffMethodRep::read from file",
		(const char *) filename, "bad centstring");
      add_level(d, sp, censtl);
    }
    else { // Otherwise we have a literal Coord for the center.
      Coord cntr;
      gmin >> cntr >> d >> sp;
      if (!gmin.good())
	::error("INPUT FAILURE in FinDiffMethodRep::read from file",
		(const char *) filename, "bad Coord-type entry");
      add_level(d, sp, cntr);
    }
    ++levels_read;
  }
  if (!levels_read) {
    cerr << "WARNING: FinDiffMethodRep::read did not find anything in the file"
      << filename << "\nProceeding with empty FinDiffMethod" << endl;
  }
}


void FinDiffMethodRep::add_level(int ngrd, float spc, Coord cntr) {
  if (ngrd%2 ==0)
    ::error("ERROR: FinDiffMethodRep::add_level: grid dimension must be odd.");
  FDMLink *lnk = new FDMLink;
  lnk->fdm_ptr = new CubeLatSpec(ngrd, spc, cntr);
  lnk->finer = 0;
  if (coarsest) {
    finest->finer = lnk;
    finest = lnk;
  }
  else
    coarsest = finest = lnk;
}

void FinDiffMethodRep::add_level(int ngrd, float spc, CenteringStyle sty)
{
  if (ngrd%2 ==0)
    ::error("ERROR: FinDiffMethodRep::add_level: grid dimension must be odd.");
  resolved = 0;
  FDMLink *lnk = new FDMLink;
  lnk->fdm_ptr = new CubeLatSpec(ngrd, spc, sty);
  lnk->finer = 0;
  if (coarsest) {
    finest->finer = lnk;
    finest = lnk;
  }
  else
    coarsest = finest = lnk;
}

void FinDiffMethodRep::resolve (Coord geom_cent, Coord center_of_intr)
{
  for (FDMLink *pl = coarsest; pl; pl=pl->finer)
    pl->fdm_ptr->resolve(geom_cent, center_of_intr);
  resolved = 1;
}

CubeLatSpec * FinDiffMethodRep::get_coarsest()
{
  if (!resolved) {
    error("ERROR: FinDiffMethodRep::get_coarsest:",
	  "attempt to access unresolved object");
  }
  if (!coarsest) {
    return 0;
  }
  ptr = coarsest; return ptr->fdm_ptr;
}

CubeLatSpec * FinDiffMethodRep::get_finer()
{
  if (!resolved) {
    error("ERROR: FinDiffMethodRep::get_coarsest:",
	  "attempt to access unresolved object");
  }
  if (ptr->finer) {ptr = ptr->finer; return ptr->fdm_ptr;}
  else return 0;
}

ostream & FinDiffMethodRep::print(ostream& ost) const
{
  if (finest == 0) {
    ost << "FindDiffMethodRep::print: No grid method currently defined\n";
    cerr << "FindDiffMethodRep::print: No grid method currently defined\n";
    return ost;
  }
  FDMLink * p = coarsest;
  ost << "                   Centering Style   grid dimension    spacing\n";
  if (p==finest) {
    ost << "Only one level: " << *p->fdm_ptr << "\n";
    return ost;
  }
  else
    ost << "Coarsest level: " << *p->fdm_ptr << "\n";
  p = p->finer;
  for (;;) {
    if (p==finest) {
      ost << "Finest level:   " << *p->fdm_ptr << "\n";
      break;
    }
    ost << "Finer level:    " << *p->fdm_ptr  << "\n";
    p = p->finer;
  }
  return ost;
}


FinDiffMethod::FinDiffMethod()
{
  rep = new FinDiffMethodRep;
}

FinDiffMethod::FinDiffMethod(const FinDiffMethod& f)
{
  f.rep->referenceCount++;
  rep = f.rep;
}

FinDiffMethod&
FinDiffMethod::operator=  (const FinDiffMethod& f)
{
  if (this != &f) {
    f.rep->referenceCount++;
    if (rep->referenceCount == 1)
      delete rep;
    else
      --rep->referenceCount;
    rep = f.rep;
  }
  return *this;
}

FinDiffMethod::~FinDiffMethod()
{
  blab3 << "FinDiffMethod destructor called with rep->referenceCount = " << rep->referenceCount << endl;
  if (rep->referenceCount == 1)
    delete rep;
  else
    --rep->referenceCount;
}

void FinDiffMethod::read (const string& filename)
{
  rep->read(filename);
}

void FinDiffMethod::add_level (int ngrd, float spc, Coord cntr)
{
  rep->add_level(ngrd, spc, cntr);
}

void FinDiffMethod::add_level(int ngrd, float spc, CenteringStyle sty)
{
  rep->add_level(ngrd, spc, sty);
}

void FinDiffMethod::resolve (Coord geom_cent, Coord center_of_intr)
{
  rep->resolve(geom_cent, center_of_intr);
}

CubeLatSpec * FinDiffMethod::get_coarsest()
{
  return rep->get_coarsest();
}

CubeLatSpec * FinDiffMethod::get_finer()
{
  return rep->get_finer();
}

ostream & FinDiffMethod::print(ostream& ost) const
{
  return rep->print(ost);
}

ostream & operator<< (ostream& ost, const FinDiffMethod& fdm)
{return fdm.print(ost);}

// FinDiffMethod.cc ends here
