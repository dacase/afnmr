/* Base class for electolyte envir, envelope for subclasses.

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

$Id: ElectrolyteEnvironment.cc,v 2.7 2012/04/17 21:35:51 bashford Exp $
*/
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/UniformElectrolyte.h"
#include "MEAD/CubeLatSpec.h"
#include "MEAD/globals.h"


int ElectrolyteEnvironment::instanceCount = 0;

// Constructs an envelope with a default letter (no electrolyte)
ElectrolyteEnvironment::ElectrolyteEnvironment()
{
  rep = new UniformElectrolyte(0.0);
}

ElectrolyteEnvironment::ElectrolyteEnvironment
(const ElectrolyteEnvironment & c)
{
  if(c.rep) { // Arg is an envelope.
    c.rep->referenceCount++;
    rep = c.rep;
  }
  else { // Arg is a letter.
    error("ERROR: ElecrolyteEnvironment copy constructor:\n",
	   "Arugments, `rep' pointer is zero\n");
  }
  ++instanceCount;
}


ElectrolyteEnvironment& ElectrolyteEnvironment::operator=
(const ElectrolyteEnvironment& c)
{
  if (this != &c) {
    if (rep) {
      if (c.rep) {
        c.rep->referenceCount++;
        if (rep->referenceCount == 1)
          delete rep;
        else
          --rep->referenceCount;
        rep = c.rep;
      }
      else {
        error("ERROR: ElectrolyteEnvironment::operator=:\n",
	      "argument has zero rep");
      }
    }
    else {
      error("ERROR ElectrolyteEnvironment::operator=: this->rep is zero");
    }
  }
  return *this;
}



// This will be an envelope containing *cdp.
ElectrolyteEnvironment::ElectrolyteEnvironment
(ElectrolyteEnvironment_lett * cdp)
{
  cdp->referenceCount++;
  rep = cdp;
  instanceCount++;
}

ElectrolyteEnvironment::~ElectrolyteEnvironment()
{
  if (rep) { // Case where this is the envelope
    blab2 << "ElectrolyteEnvironment destructor called with rep->referenceCount = " << rep->referenceCount << endl;
    if (rep->referenceCount == 1)
      delete rep;
    else
      --rep->referenceCount;
  }
  instanceCount--;
}

int ElectrolyteEnvironment_lett::instanceCount = 0;


ElectrolyteEnvironment_lett::ElectrolyteEnvironment_lett()
{
  blab2 << "ElectrolyteEnvironment_lett default constructor called" << endl;
  referenceCount=1;
  ++instanceCount;
}


ElectrolyteEnvironment_lett::ElectrolyteEnvironment_lett
(const ElectrolyteEnvironment_lett& a)
{
  blab2 << "ElectrolyteEnvironment_lett copy constructor called" << endl;
  referenceCount = 1;
  ++instanceCount;
}

// Don't delete this letter until all references have been released
void ElectrolyteEnvironment_lett::operator delete(void *p, size_t s)
{
  ElectrolyteEnvironment_lett *del = reinterpret_cast<ElectrolyteEnvironment_lett*> (p);
  blab2 << "ElectrolyteEnvironment_lett operator delete called with referenceCount = " << del->referenceCount << endl;
  if (del->referenceCount == 0) {
    blab2 << "ElectrolyteEnvironment_lett operator delete frees object" << endl;
    ::operator delete(p);
  }
}

ElectrolyteEnvironment_lett::~ElectrolyteEnvironment_lett()
{
  blab2 << "ElectrolyteEnvironment_lett destructor called with referenceCount = " << referenceCount << endl;
  if (--referenceCount)
    // It's now OK for this to happen
    blab2 << "WARNING: ElectrolyteEnvironment_lett destructor called\n"
      << "on object with referenceCount = " << referenceCount+1
	<< "\nas if the object is still in use somewhere" << endl;
  --instanceCount;
}

ElectrolyteEnvironment_lett&
ElectrolyteEnvironment_lett::operator=(const ElectrolyteEnvironment_lett& acs)
{return *this;}

ElyCubeRep::ElyCubeRep(const CubeLatSpec& cls)
{
  int n = cls.get_grid_dim();
  isgarr = new int[n*n*n];
}
ElyCubeRep::~ElyCubeRep()
{
  delete [] isgarr;
}

// ElectrolyteEnvironment.cc ends here
