/* Base class for electrostatic potential, envelope for subclasses.

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

$Id: ElstatPot.cc,v 2.11 2004/12/06 17:58:16 bashford Exp $
*/
#include "MEAD/ElstatPot.h"
#include "MEAD/ElstatMaker.h"
#include "MEAD/FDElstatMaker.h"
#include "MEAD/DielectricEnvironment.h"
#include "MEAD/ElectrolyteEnvironment.h"
#include "MEAD/ChargeDist.h"
#include "MEAD/FinDiffMethod.h"
#include "MEAD/CubeLatSpec.h"
#include <iostream>


int ElstatPot::instanceCount = 0;

// This subclass (known only to this module) is a default letter
// for the ElstatPot envelope.  I think it is fairly harmless
class ZeroPotential : public AnalyticEP {
public:
  ZeroPotential() : AnalyticEP() {}
  virtual float value(Coord c) const
    {cerr << "WARNING: ZeroPotental::value will always return zero\n"
       << "This ZeroPotential object was probably created as a default\n"
	 << "For an ElstatPot that was declared with no arguments" << endl;
     return 0;}
  virtual Coord field(Coord x) const
    {cerr << "WARNING: ZeroPotental::field will always return zero\n"
       << "This ZeroPotential object was probably created as a default\n"
	 << "For an ElstatPot that was declared with no arguments" << endl;
     return Coord(0.0, 0.0, 0.0);}
  virtual Coord displacement(Coord x) const
    {cerr << "WARNING: ZeroPotental::displacement will always return zero\n"
       << "This ZeroPotential object was probably created as a default\n"
	 << "For an ElstatPot that was declared with no arguments" << endl;
     return Coord(0.0, 0.0, 0.0);}

  virtual void solve() {}
//  virtual float operator * (const ChargeDist& c) const
//    {cerr << "WARNING: ZeroPotental::operator* will always return zero\n"
//       << "This ZeroPotential object was probably created as a default\n"
//	 << "For an ElstatPot that was declared with no arguments" << endl;
//     return 0.0;}
private:
};


// The real member functions for ElstatPot start here.
ElstatPot::ElstatPot (DielectricEnvironment de, ChargeDist cd,
		      ElectrolyteEnvironment ely)
{
  rep = ElstatMaker::maker(de.rep, cd.rep, ely.rep);
  ++instanceCount;
}


ElstatPot::ElstatPot (FinDiffMethod fdm,
		      DielectricEnvironment e, ChargeDist r,
		      ElectrolyteEnvironment ely)
{
  rep = FDElstatMaker::maker(fdm, e.rep, r.rep, ely.rep);
  ++instanceCount;
}

ElstatPot::ElstatPot (ElstatPot_lett* epp)
{
  rep = epp;
  ++rep->referenceCount;
  ++instanceCount;
}

// Now the "canonical" stuff
// Default constructor puts a new ZeroPotential letter in the new envelope
ElstatPot::ElstatPot()
{
  rep = new ZeroPotential();
}

// Copy ctor argument must be a valid envelope (rep non-zero).
ElstatPot::ElstatPot(const ElstatPot & c)
{
  if(c.rep) { // Arg is an envelope.
    c.rep->referenceCount++;
    rep = c.rep;
  }
  else
    error ("ERROR: ElstatPot copy constructor:\n"
	   "argument has zero rep pointer");
  ++instanceCount;
}

// Both sides of assignment must be a valid envelope (rep non-zero).
ElstatPot& ElstatPot::operator = (const ElstatPot& c)
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
         cerr << "ERROR: ElstatPot::operator=:\n"
	   << "argument has zero rep" << endl;
       }
     }
     else {
       cerr << "ERROR ElstatPot::operator=: this->rep is zero" << endl;
     }
  }
  return *this;
}


// FIXME?  Should there be an EstatPot(ElstatPot*) constructor?


ElstatPot::~ElstatPot()
{
  if (rep) {
    if (rep->referenceCount == 1)
      delete rep;
    else
      --rep->referenceCount;
  }
  instanceCount--;
}

int ElstatPot_lett::instanceCount = 0;


ElstatPot_lett::ElstatPot_lett()
{
  rho_obj = 0;
  eps_obj = 0;
  electrolyte_obj = 0;
  referenceCount=1;
  ++instanceCount;
}


ElstatPot_lett::ElstatPot_lett(const ElstatPot_lett& a)
{
  referenceCount = 1;
  ++instanceCount;
}

ElstatPot_lett::ElstatPot_lett (DielectricEnvironment_lett* e,
				ChargeDist_lett* r,
				ElectrolyteEnvironment_lett* ely)
: rho_obj(new ChargeDist(r)),
  eps_obj(new DielectricEnvironment(e)),
  electrolyte_obj(new ElectrolyteEnvironment(ely))
{
  referenceCount=1;
  ++instanceCount;
}

// Don't delete this letter until all references have been released
void ElstatPot_lett::operator delete(void *p, size_t s)
{
  ElstatPot_lett *del = reinterpret_cast<ElstatPot_lett*> (p);
  blab2 << "ElstatPot_lett operator delete called with referenceCount = " << del->referenceCount << endl;
  if (del->referenceCount == 0) {
    blab2 << "ElstatPot_lett operator delete frees object" << endl;
    ::operator delete(p);
  }
}

ElstatPot_lett::~ElstatPot_lett()
{
  blab2 << "ElstatPot_lett destructor called with refereceCount = " << referenceCount << endl;
  if (--referenceCount)
    // It's now OK for this to happen
    blab2 << "WARNING: ElstatPot_lett destructor called\n"
      << "on object with referenceCount = " << referenceCount+1
	<< "\nas if the object is still in use somewhere" << endl;
  // It's possible to have constructed an empty letter
  if (rho_obj) delete rho_obj;
  if (eps_obj) delete eps_obj;
  if (electrolyte_obj) delete electrolyte_obj;
  --instanceCount;
}

ElstatPot_lett& ElstatPot_lett::operator=(const ElstatPot_lett& acs)
{return *this;}

// ElstatPot.cc ends here
