/* Base class for dielectric environment and subclasses.

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

$Id: DielectricEnvironment.cc,v 2.7 2004/12/06 17:58:12 bashford Exp $
*/
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/UniformDielectric.h"
#include <iostream>


int DielectricEnvironment::instanceCount = 0;

// Constructs an envelope with a ChargDist letter, that is, a default letter.
DielectricEnvironment::DielectricEnvironment()
{
  rep = new UniformDielectric(1.0);
}

// Ctor must have a valid envelope (non-zero rep) as its argument.
DielectricEnvironment::DielectricEnvironment(const DielectricEnvironment & c)
{
  if(c.rep) { // Arg is an envelope.
    c.rep->referenceCount++;
    rep = c.rep;
  }
  else
    error("ERROR: DielectricEnvironment copy constructor:\n",
	  "argument has zero rep pointer");
  ++instanceCount;
}

// Both both sides of assignment must be envelopes (non-zero rep)

DielectricEnvironment& DielectricEnvironment::operator=
(const DielectricEnvironment& c)
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
        error("ERROR: DielectricEnvironment::operator=:\n",
	      "argument has zero rep");
      }
    }
    else {
      error("ERROR DielectricEnvironment::operator=: this->rep is zero");
    }
  }
  return *this;
}


// This will create an envelope containing *cdp.
DielectricEnvironment::DielectricEnvironment(DielectricEnvironment_lett * cdp)
{
  cdp->referenceCount++;
  rep = cdp;
  instanceCount++;
}

DielectricEnvironment::~DielectricEnvironment()
{
  if (rep) {
    blab2 << "DielectricEnvironment destructor called with rep->referenceCount = " << rep->referenceCount << endl;
    if (rep->referenceCount == 1)
      delete rep;
    else
      --rep->referenceCount;
  }
  instanceCount--;
}

int DielectricEnvironment_lett::instanceCount = 0;


DielectricEnvironment_lett::DielectricEnvironment_lett()
{
  blab2 << "DielectricEnvironment_lett default constructor called" << endl;
  referenceCount=1;
  ++instanceCount;
}


DielectricEnvironment_lett::DielectricEnvironment_lett
(const DielectricEnvironment_lett& a)
{
  blab2 << "DielectricEnvironment_lett copy constructor called" << endl;
  referenceCount = 1;
  ++instanceCount;
}

// Don't delete this letter until all references have been released
void DielectricEnvironment_lett::operator delete(void *p, size_t s)
{
  DielectricEnvironment_lett *del = reinterpret_cast<DielectricEnvironment_lett*> (p);
  blab2 << "DielectricEnvironment_lett operator delete called with referenceCount = " << del->referenceCount << endl;
  if (del->referenceCount == 0) {
    blab2 << "DielectricEnvironment_lett operator delete frees object" << endl;
    ::operator delete(p);
  }
}

DielectricEnvironment_lett::~DielectricEnvironment_lett()
{
  blab2 << "DielectricEnvironment_lett destructor called with referenceCount = " << referenceCount << endl;
  if (--referenceCount)
    // It's now OK for this to happen
    blab2 << "WARNING: DielectricEnvironment_lett destructor called\n"
      << "on object with referenceCount = " << referenceCount+1
	<< "\nas if the object is still in use somewhere" << endl;
  --instanceCount;
}

DielectricEnvironment_lett&
DielectricEnvironment_lett::operator=(const DielectricEnvironment_lett& acs)
{return *this;}

// DielectricEnvironment.cc ends here
